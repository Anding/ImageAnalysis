#define _CRT_SECURE_NO_WARNINGS

#include "image_analysis.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define XISF_HEADER_SIZE 8192u
#define MAX_DETECTED_STARS 256u
#define MAX_LOCKED_STARS 7u
#define MAX_FOCUS_FRAMES 64u

typedef struct XISFImage {
    uint16_t *pixels;
    uint32_t width;
    uint32_t height;
} XISFImage;

static void calibrate_models(const IAFocusSample *samples, uint32_t count)
{
    IAFocusModel model;
    IAFocusModel best_piecewise;
    IAFocusModel best_hyperbola;
    IAFocusFit fit;
    IAFocusFit best_piecewise_fit;
    IAFocusFit best_hyperbola_fit;
    int32_t left_slope;
    int32_t right_slope;
    int32_t radius;
    int32_t slope;

    memset(&best_piecewise, 0, sizeof(best_piecewise));
    memset(&best_hyperbola, 0, sizeof(best_hyperbola));
    memset(&best_piecewise_fit, 0, sizeof(best_piecewise_fit));
    memset(&best_hyperbola_fit, 0, sizeof(best_hyperbola_fit));
    best_piecewise_fit.rms_residual_milli = INT32_MAX;
    best_hyperbola_fit.rms_residual_milli = INT32_MAX;

    memset(&model, 0, sizeof(model));
    model.struct_size = sizeof(model);
    model.kind = IA_FOCUS_PIECEWISE_V;
    for (left_slope = 40000; left_slope <= 130000; left_slope += 2000) {
        for (right_slope = 40000; right_slope <= 130000; right_slope += 2000) {
            model.left_slope_micro_per_step = left_slope;
            model.right_slope_micro_per_step = right_slope;
            fit.struct_size = sizeof(fit);
            if (IA_FitFocusModel(
                    &model, samples, count, UINT32_MAX, 250u, &fit) == IA_OK &&
                fit.rms_residual_milli < best_piecewise_fit.rms_residual_milli) {
                best_piecewise = model;
                best_piecewise_fit = fit;
            }
        }
    }

    memset(&model, 0, sizeof(model));
    model.struct_size = sizeof(model);
    model.kind = IA_FOCUS_HYPERBOLA;
    for (radius = 1000; radius <= 6000; radius += 100) {
        for (slope = 40000; slope <= 130000; slope += 2000) {
            model.hyperbola_radius_milli = radius;
            model.hyperbola_slope_micro_per_step = slope;
            fit.struct_size = sizeof(fit);
            if (IA_FitFocusModel(
                    &model, samples, count, UINT32_MAX, 250u, &fit) == IA_OK &&
                fit.rms_residual_milli < best_hyperbola_fit.rms_residual_milli) {
                best_hyperbola = model;
                best_hyperbola_fit = fit;
            }
        }
    }

    fprintf(stderr,
        "calibrated piecewise-V slopes %d/%d, focus %.3f, residual %.3f\n",
        best_piecewise.left_slope_micro_per_step,
        best_piecewise.right_slope_micro_per_step,
        best_piecewise_fit.focus_milli_steps / 1000.0,
        best_piecewise_fit.rms_residual_milli / 1000.0);
    fprintf(stderr,
        "calibrated hyperbola radius %d, slope %d, focus %.3f, residual %.3f\n",
        best_hyperbola.hyperbola_radius_milli,
        best_hyperbola.hyperbola_slope_micro_per_step,
        best_hyperbola_fit.focus_milli_steps / 1000.0,
        best_hyperbola_fit.rms_residual_milli / 1000.0);
}

static int load_xisf(const char *path, XISFImage *image)
{
    FILE *file;
    char header[XISF_HEADER_SIZE + 1u];
    char *geometry;
    unsigned int depth;
    size_t pixels;

    memset(image, 0, sizeof(*image));
    file = fopen(path, "rb");
    if (file == NULL) {
        fprintf(stderr, "Cannot open %s\n", path);
        return 0;
    }
    if (fread(header, 1u, XISF_HEADER_SIZE, file) != XISF_HEADER_SIZE) {
        fprintf(stderr, "Cannot read XISF header from %s\n", path);
        fclose(file);
        return 0;
    }
    header[XISF_HEADER_SIZE] = '\0';
    geometry = strstr(header, "geometry=\"");
    if (geometry == NULL ||
        sscanf(geometry, "geometry=\"%u:%u:%u\"",
            &image->width, &image->height, &depth) != 3 ||
        depth != 1u) {
        fprintf(stderr, "Unsupported XISF geometry in %s\n", path);
        fclose(file);
        return 0;
    }
    pixels = (size_t)image->width * image->height;
    image->pixels = (uint16_t *)malloc(pixels * sizeof(uint16_t));
    if (image->pixels == NULL ||
        fread(image->pixels, sizeof(uint16_t), pixels, file) != pixels) {
        fprintf(stderr, "Cannot read XISF pixels from %s\n", path);
        free(image->pixels);
        image->pixels = NULL;
        fclose(file);
        return 0;
    }
    fclose(file);
    return 1;
}

static int focus_position(const char *path)
{
    const char *name = strrchr(path, '\\');
    int focus = 0;
    if (name == NULL) name = strrchr(path, '/');
    if (name == NULL) name = path; else ++name;
    if (sscanf(name, "F%d-", &focus) != 1) {
        return 0;
    }
    return focus;
}

int main(int argc, char **argv)
{
    IAConfig config;
    IAStar detected[MAX_DETECTED_STARS];
    IAStar locked[MAX_LOCKED_STARS];
    IAStar measurements[MAX_LOCKED_STARS];
    IAFrameSummary summary;
    void *workspace;
    uint32_t workspace_bytes;
    uint32_t locked_count = 0u;
    IAFocusSample samples[MAX_LOCKED_STARS][MAX_FOCUS_FRAMES];
    IAFocusSample frame_samples[MAX_FOCUS_FRAMES];
    uint32_t sample_count = 0u;
    XISFImage discovery;
    int status;
    int argument;
    uint32_t i;

    if (argc < 2 || argc - 1 > MAX_FOCUS_FRAMES) {
        fprintf(stderr, "Usage: focus_corpus central-frame.xisf [frame.xisf ...]\n");
        return 2;
    }

    IA_DefaultConfig(&config);
    config.aperture_radius_pixels = 20u;
    config.annulus_inner_radius_pixels = 22u;
    config.annulus_outer_radius_pixels = 30u;
    config.minimum_separation_pixels = 12u;
    workspace_bytes = IA_WorkspaceBytes(MAX_DETECTED_STARS);
    workspace = malloc(workspace_bytes);
    if (workspace == NULL || !load_xisf(argv[1], &discovery)) {
        free(workspace);
        return 1;
    }

    memset(&summary, 0, sizeof(summary));
    summary.struct_size = sizeof(summary);
    status = IA_AnalyzeFrame(
        discovery.pixels, discovery.width, discovery.height, discovery.width,
        &config, workspace, workspace_bytes, detected, MAX_DETECTED_STARS, &summary);
    if (status != IA_OK) {
        fprintf(stderr, "Discovery analysis failed: %d\n", status);
        free(discovery.pixels);
        free(workspace);
        return 1;
    }
    for (i = 0u; i < summary.measured_stars && locked_count < MAX_LOCKED_STARS; ++i) {
        uint32_t reject = IA_STAR_SATURATED | IA_STAR_EDGE |
            IA_STAR_TOO_SMALL | IA_STAR_LOW_SNR;
        if ((detected[i].flags & reject) == 0u) {
            locked[locked_count++] = detected[i];
        }
    }
    fprintf(stderr, "Locked %u of %u measured stars from %s\n",
        locked_count, summary.measured_stars, argv[1]);
    free(discovery.pixels);
    if (locked_count == 0u) {
        free(workspace);
        return 1;
    }

    printf("focus,milliseconds,included,background,noise,median_hfd,mad_hfd,median_fwhm");
    for (i = 0u; i < locked_count; ++i) {
        printf(",star%u_hfd", i + 1u);
    }
    putchar('\n');

    for (argument = 1; argument < argc; ++argument) {
        XISFImage image;
        clock_t start;
        clock_t finish;
        double milliseconds;

        if (!load_xisf(argv[argument], &image)) {
            free(workspace);
            return 1;
        }
        if (image.width != discovery.width && discovery.width != 0u) {
            fprintf(stderr, "Frame dimensions changed in %s\n", argv[argument]);
        }
        start = clock();
        status = IA_MeasureLockedStars(
            image.pixels, image.width, image.height, image.width, &config,
            locked, locked_count, 16u, measurements);
        finish = clock();
        if (status != IA_OK) {
            fprintf(stderr, "Locked-star measurement failed for %s: %d\n",
                argv[argument], status);
            free(image.pixels);
            free(workspace);
            return 1;
        }
        summary.struct_size = sizeof(summary);
        status = IA_SummarizeStars(
            measurements, locked_count,
            IA_STAR_SATURATED | IA_STAR_EDGE | IA_STAR_TOO_SMALL | IA_STAR_LOW_SNR,
            workspace, workspace_bytes, &summary);
        frame_samples[sample_count].focus_position = focus_position(argv[argument]);
        frame_samples[sample_count].metric_milli =
            summary.median_hfd_milli_pixels;
        frame_samples[sample_count].weight_milli =
            summary.included_stars * IA_MILLI;
        frame_samples[sample_count].flags = status == IA_OK ? 0u : 1u;
        milliseconds = 1000.0 * (finish - start) / CLOCKS_PER_SEC;
        printf("%d,%.3f,%u,%.3f,%.3f,%.3f,%.3f,%.3f",
            focus_position(argv[argument]), milliseconds, summary.included_stars,
            measurements[0].background_milli_adu / 1000.0,
            measurements[0].noise_milli_adu / 1000.0,
            summary.median_hfd_milli_pixels / 1000.0,
            summary.mad_hfd_milli_pixels / 1000.0,
            summary.median_fwhm_milli_pixels / 1000.0);
        for (i = 0u; i < locked_count; ++i) {
            printf(",%.3f", measurements[i].hfd_milli_pixels / 1000.0);
            samples[i][sample_count].focus_position = focus_position(argv[argument]);
            samples[i][sample_count].metric_milli = measurements[i].hfd_milli_pixels;
            samples[i][sample_count].weight_milli = IA_MILLI;
            samples[i][sample_count].flags = measurements[i].flags;
        }
        putchar('\n');
        ++sample_count;
        free(image.pixels);
    }

    {
        IAFocusModel models[2];
        const char *names[2] = { "piecewise-V", "hyperbola" };
        uint32_t model_index;

        memset(models, 0, sizeof(models));
        models[0].struct_size = sizeof(models[0]);
        models[0].kind = IA_FOCUS_PIECEWISE_V;
        models[0].left_slope_micro_per_step = 85000;
        models[0].right_slope_micro_per_step = 85000;
        models[1].struct_size = sizeof(models[1]);
        models[1].kind = IA_FOCUS_HYPERBOLA;
        models[1].hyperbola_radius_milli = 3850;
        models[1].hyperbola_slope_micro_per_step = 108000;

        for (model_index = 0u; model_index < 2u; ++model_index) {
            IAFocusFit fits[MAX_LOCKED_STARS];
            IAFocusSummary focus_summary;
            int32_t combine_workspace[MAX_LOCKED_STARS];
            uint32_t fitted = 0u;

            for (i = 0u; i < locked_count; ++i) {
                IAFocusFit star_fit;
                star_fit.struct_size = sizeof(star_fit);
                if (IA_FitFocusModel(
                    &models[model_index], samples[i], sample_count,
                    IA_STAR_SATURATED | IA_STAR_EDGE | IA_STAR_TOO_SMALL |
                        IA_STAR_LOW_SNR,
                    250u, &star_fit) == IA_OK) {
                    fits[fitted++] = star_fit;
                }
            }
            focus_summary.struct_size = sizeof(focus_summary);
            status = IA_CombineFocusFits(
                fits, fitted, IA_FOCUS_ONE_SIDED | IA_FOCUS_AT_BOUNDARY |
                    IA_FOCUS_INVALID_MODEL | IA_FOCUS_INSUFFICIENT_POINTS,
                combine_workspace, sizeof(combine_workspace), &focus_summary);
            if (status == IA_OK) {
                fprintf(stderr,
                    "%s focus %.3f, star MAD %.3f, median residual %.3f (%u stars)\n",
                    names[model_index],
                    focus_summary.median_focus_milli_steps / 1000.0,
                    focus_summary.mad_focus_milli_steps / 1000.0,
                    focus_summary.median_rms_residual_milli / 1000.0,
                    focus_summary.fits_included);
            } else {
                fprintf(stderr, "%s fit unavailable: %d\n",
                    names[model_index], status);
            }
        }
    }
    calibrate_models(frame_samples, sample_count);

    free(workspace);
    return 0;
}
