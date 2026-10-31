#include "focus_curve.h"
#include "analysis_parameters.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct FocusPoint {
    int32_t position;
    double hfd[MAX_STARS];
} FocusPoint;

typedef struct Identity {
    double x, y, quality;
} Identity;

typedef struct Trial {
    int valid;
    int32_t position;
    int low_risk, samples;
    double risk_cost;
    double hfd[MAX_STARS];
    double x[MAX_STARS], y[MAX_STARS];
    unsigned char seen[MAX_STARS];
} Trial;

typedef struct FitResult {
    double focus, uncertainty;
} FitResult;

static struct {
    int active, selected, limit, tracked;
    Identity identities[MAX_STARS];
    FocusPoint *points;
    size_t count, capacity;
    Trial trial;
} run;

static void reset(void)
{
    free(run.points);
    memset(&run, 0, sizeof(run));
}

int CALL start_computation(uint32_t max_stars)
{
    reset();
    if (max_stars == 0 || max_stars > MAX_STARS)
        return -1;
    run.active = 1;
    run.limit = (int)max_stars;
    return 0;
}

static int finite_centroid(const StarRecord *star)
{
    return isfinite(star->x) && isfinite(star->y);
}

static int valid_risk(const StarRecord *star)
{
    return isfinite(star->over_exposed_risk) &&
        isfinite(star->under_exposed_risk) &&
        star->over_exposed_risk >= 0.0 && star->over_exposed_risk <= 1.0 &&
        star->under_exposed_risk >= 0.0 && star->under_exposed_risk <= 1.0;
}

static int valid_scale(const StarRecord *star)
{
    return isfinite(star->exposure_scale) && star->exposure_scale > 0.0;
}

static int exposure_direction(const StarRecord *star)
{
    int over, under;
    if (!valid_risk(star) || !valid_scale(star))
        return 0;
    over = star->over_exposed_risk > MAX_EXPOSURE_RISK &&
        star->exposure_scale < 1.0;
    under = star->under_exposed_risk > MAX_EXPOSURE_RISK &&
        star->exposure_scale > 1.0;
    if (over && under)
        return star->over_exposed_risk >= star->under_exposed_risk ? -1 : 1;
    return over ? -1 : (under ? 1 : 0);
}

static int valid_sample(const StarRecord *star)
{
    return isfinite(star->hfd) && star->hfd > 0.0;
}

static double clamp_double(double value, double low, double high)
{
    return value < low ? low : (value > high ? high : value);
}

static int compare_double(const void *left, const void *right)
{
    double a = *(const double *)left, b = *(const double *)right;
    return (a > b) - (a < b);
}

static double median(double *values, size_t count)
{
    qsort(values, count, sizeof(*values), compare_double);
    return count & 1 ? values[count / 2] :
        values[count / 2 - 1] / 2.0 + values[count / 2] / 2.0;
}

static int record_count(const StarRecord *stars)
{
    int count = 0;
    while (count < MAX_STARS && stars[count].positive_pixel_count != 0)
        ++count;
    return count;
}

static int duplicate_centroid(const StarRecord *stars, int index, int count)
{
    int i;
    double radius2 = DUPLICATE_RADIUS_PIXELS * DUPLICATE_RADIUS_PIXELS;
    for (i = 0; i < index && i < count; ++i) {
        double dx, dy;
        if (!finite_centroid(stars + i))
            continue;
        dx = stars[index].x - stars[i].x;
        dy = stars[index].y - stars[i].y;
        if (dx * dx + dy * dy <= radius2)
            return 1;
    }
    return 0;
}

static int collect_candidates(const StarRecord *stars, Identity *identities,
    int *record_indices, int count)
{
    int candidate_indices[MAX_STARS];
    double fluxes[MAX_STARS];
    double flux_median;
    int i, candidates = 0;

    for (i = 0; i < count; ++i) {
        if (finite_centroid(stars + i) && isfinite(stars[i].flux) &&
            stars[i].flux > 0.0 && isfinite(stars[i].roundness) &&
            !duplicate_centroid(stars, i, count))
            candidate_indices[candidates++] = i;
    }
    if (candidates == 0)
        return 0;
    for (i = 0; i < candidates; ++i)
        fluxes[i] = stars[candidate_indices[i]].flux;
    flux_median = median(fluxes, (size_t)candidates);
    for (i = 0; i < candidates; ++i) {
        int index = candidate_indices[i];
        double roundness, brightness;
        roundness = clamp_double(stars[index].roundness, 0.0, 1.0);
        brightness = clamp_double(stars[index].flux / (2.0 * flux_median), 0.0, 1.0);
        identities[i].x = stars[index].x;
        identities[i].y = stars[index].y;
        identities[i].quality =
            ROUNDNESS_WEIGHT * roundness + BRIGHTNESS_WEIGHT * brightness;
        if (record_indices)
            record_indices[i] = index;
    }
    return candidates;
}

static void sort_indices_by_quality(Identity *identities, int *indices,
    int count)
{
    int i, j;
    for (i = 1; i < count; ++i) {
        int index = indices[i];
        Identity key = identities[i];
        j = i;
        while (j > 0 && identities[j - 1].quality < key.quality) {
            identities[j] = identities[j - 1];
            indices[j] = indices[j - 1];
            --j;
        }
        identities[j] = key;
        indices[j] = index;
    }
}

static uint32_t bounded_exposure(uint32_t exposure)
{
    if (exposure < MIN_EXPOSURE_US)
        return MIN_EXPOSURE_US;
    if (exposure > MAX_EXPOSURE_US)
        return MAX_EXPOSURE_US;
    return exposure;
}

static int recommend_exposure(const StarRecord *stars, const int *indices,
    int count, int top_count, int baseline, int low_risk, uint32_t *exposure)
{
    double scales[MAX_STARS], target;
    double over_weight = 0.0, under_weight = 0.0;
    uint32_t current = bounded_exposure(*exposure), recommended;
    int i, used = 0, needed = run.limit - low_risk;
    int over_count = 0, under_count = 0, choose_under;

    if (!baseline && needed <= 0)
        return 0;
    if (baseline)
        needed = run.limit;
    for (i = 0; i < count && i < top_count; ++i) {
        const StarRecord *star = stars + indices[i];
        if (!valid_scale(star))
            continue;
        if (baseline) {
            scales[used++] = star->exposure_scale;
            if (used == needed)
                break;
            continue;
        }
        if (!valid_risk(star))
            continue;
        if (exposure_direction(star) < 0) {
            ++over_count;
            over_weight += star->over_exposed_risk;
        }
        if (exposure_direction(star) > 0) {
            ++under_count;
            under_weight += star->under_exposed_risk;
        }
    }
    if (baseline) {
        if (used == 0)
            return 0;
        target = median(scales, (size_t)used);
    } else {
        if (over_count == 0 && under_count == 0)
            return 0;
        choose_under = under_count > over_count ||
            (under_count == over_count && under_weight > over_weight);
        used = 0;
        for (i = 0; i < count && i < top_count &&
            used < needed; ++i) {
            const StarRecord *star = stars + indices[i];
            int direction = exposure_direction(star);
            if ((choose_under && direction > 0) ||
                (!choose_under && direction < 0))
                scales[used++] = star->exposure_scale;
        }
        if (used == 0)
            return 0;
        target = median(scales, (size_t)used);
    }
    target = clamp_double((double)current * target,
        (double)MIN_EXPOSURE_US, (double)MAX_EXPOSURE_US);
    recommended = (uint32_t)floor(target + 0.5);
    if (recommended == current ||
        fabs((double)recommended - current) / current <
            MIN_EXPOSURE_CHANGE_FRACTION) {
        *exposure = current;
        return 0;
    }
    *exposure = recommended;
    return 1;
}

int CALL select_focus_stars(const StarRecord *stars, uint32_t *exposure)
{
    Identity candidates[MAX_STARS];
    int indices[MAX_STARS], count, i, selected, top_count;
    uint32_t bounded;

    if (!run.active || run.selected || !stars || !exposure)
        return -1;
    run.selected = 1;
    count = record_count(stars);
    selected = collect_candidates(stars, candidates, indices, count);
    if (selected == 0)
        return -1;
    sort_indices_by_quality(candidates, indices, selected);
    run.tracked = selected;
    for (i = 0; i < selected; ++i)
        run.identities[i] = candidates[i];

    bounded = bounded_exposure(*exposure);
    *exposure = bounded;
    top_count = AF_EXPOSURE_TOP_STARS > run.limit ?
        AF_EXPOSURE_TOP_STARS : run.limit;
    if (top_count > selected)
        top_count = selected;
    (void)recommend_exposure(stars, indices, top_count, top_count, 1, 0,
        exposure);
    return 0;
}

static int nearest_detection(const StarRecord *stars, int count, double x,
    double y, int *tie)
{
    double best = MATCH_RADIUS_PIXELS * MATCH_RADIUS_PIXELS;
    int selected = -1, ambiguous = 0, i;
    for (i = 0; i < count; ++i) {
        double dx, dy, distance;
        if (!finite_centroid(stars + i))
            continue;
        dx = stars[i].x - x;
        dy = stars[i].y - y;
        distance = dx * dx + dy * dy;
        if (distance < best || (selected < 0 && distance == best)) {
            best = distance;
            selected = i;
            ambiguous = 0;
        } else if (distance == best && selected >= 0) {
            ambiguous = 1;
        }
    }
    *tie = ambiguous;
    return ambiguous ? -1 : selected;
}

static int match_identity(const StarRecord *stars, int count, int identity)
{
    int tie, selected, nearest_identity = -1, identity_tie = 0, i;
    selected = nearest_detection(stars, count, run.identities[identity].x,
        run.identities[identity].y, &tie);
    if (selected < 0 || tie)
        return -1;
    {
        double best = MATCH_RADIUS_PIXELS * MATCH_RADIUS_PIXELS;
        for (i = 0; i < run.tracked; ++i) {
            double dx = stars[selected].x - run.identities[i].x;
            double dy = stars[selected].y - run.identities[i].y;
            double distance = dx * dx + dy * dy;
            if (distance < best || (nearest_identity < 0 && distance == best)) {
                best = distance;
                nearest_identity = i;
                identity_tie = 0;
            } else if (distance == best) {
                identity_tie = 1;
            }
        }
    }
    return !identity_tie && nearest_identity == identity ? selected : -1;
}

static int better_trial(const Trial *candidate, const Trial *best)
{
    if (!best->valid)
        return 1;
    if (candidate->low_risk != best->low_risk)
        return candidate->low_risk > best->low_risk;
    if (candidate->risk_cost != best->risk_cost)
        return candidate->risk_cost < best->risk_cost;
    return candidate->samples > best->samples;
}

int CALL add_focus_point(const StarRecord *stars, int32_t position,
    uint32_t *exposure)
{
    int indices[MAX_STARS], low_risk = 0, count, i;
    Trial candidate;
    int exposure_changed;
    uint32_t original_exposure;

    if (!run.active || !run.selected || !stars || !exposure)
        return -1;
    original_exposure = *exposure;
    *exposure = bounded_exposure(*exposure);
    count = record_count(stars);
    memset(&candidate, 0, sizeof(candidate));
    candidate.position = position;
    for (i = 0; i < run.tracked; ++i) {
        int index = match_identity(stars, count, i);
        indices[i] = index;
        if (index < 0)
            continue;
        candidate.seen[i] = 1;
        candidate.x[i] = stars[index].x;
        candidate.y[i] = stars[index].y;
        if (valid_sample(stars + index)) {
            candidate.hfd[i] = stars[index].hfd;
            ++candidate.samples;
            if (valid_risk(stars + index)) {
                candidate.risk_cost += stars[index].over_exposed_risk +
                    stars[index].under_exposed_risk;
                if (stars[index].over_exposed_risk <= MAX_EXPOSURE_RISK &&
                    stars[index].under_exposed_risk <= MAX_EXPOSURE_RISK)
                    ++low_risk;
            } else {
                candidate.risk_cost += 2.0;
            }
        }
    }
    if (candidate.samples == 0)
        return -1;
    candidate.risk_cost += 2.0 * (run.tracked - candidate.samples);
    candidate.low_risk = low_risk;
    candidate.valid = 1;
    if (better_trial(&candidate, &run.trial))
        run.trial = candidate;

    count = 0;
    for (i = 0; i < run.tracked && count < MAX_STARS; ++i) {
        if (indices[i] >= 0)
            indices[count++] = indices[i];
    }
    {
        int top_count = AF_EXPOSURE_TOP_STARS > run.limit ?
            AF_EXPOSURE_TOP_STARS : run.limit;
        if (top_count > count)
            top_count = count;
        exposure_changed = recommend_exposure(stars, indices, count,
            top_count, 0, low_risk, exposure);
    }
    return exposure_changed || original_exposure != *exposure ? 1 : 0;
}

static int commit_trial(void)
{
    FocusPoint *grown;
    size_t capacity;
    int i;
    if (!run.trial.valid)
        return -1;
    if (run.count == run.capacity) {
        capacity = run.capacity ? run.capacity * 2 : 16;
        if (capacity < run.capacity || capacity > SIZE_MAX / sizeof(*grown))
            return -1;
        grown = (FocusPoint *)realloc(run.points, capacity * sizeof(*grown));
        if (!grown)
            return -1;
        run.points = grown;
        run.capacity = capacity;
    }
    memset(&run.points[run.count], 0, sizeof(run.points[run.count]));
    run.points[run.count].position = run.trial.position;
    for (i = 0; i < run.tracked; ++i) {
        run.points[run.count].hfd[i] = run.trial.hfd[i];
        if (run.trial.seen[i]) {
            run.identities[i].x = run.trial.x[i];
            run.identities[i].y = run.trial.y[i];
        }
    }
    ++run.count;
    memset(&run.trial, 0, sizeof(run.trial));
    return 0;
}

int CALL finish_focus_point(void)
{
    if (!run.active || !run.selected)
        return -1;
    if (!run.trial.valid)
        return -1;
    return commit_trial();
}

static double model_value(double h, double slope, double focus, double x)
{
    double offset = slope * (x - focus);
    return hypot(h, offset);
}

static int solve3(double matrix[3][3], double rhs[3], double result[3])
{
    double augmented[3][4], scale = 0.0;
    int row, column, pivot;
    for (row = 0; row < 3; ++row) {
        for (column = 0; column < 3; ++column) {
            augmented[row][column] = matrix[row][column];
            if (fabs(matrix[row][column]) > scale)
                scale = fabs(matrix[row][column]);
        }
        augmented[row][3] = rhs[row];
    }
    if (!(scale > 0.0) || !isfinite(scale))
        return 0;
    for (column = 0; column < 3; ++column) {
        double largest = fabs(augmented[column][column]);
        pivot = column;
        for (row = column + 1; row < 3; ++row) {
            if (fabs(augmented[row][column]) > largest) {
                largest = fabs(augmented[row][column]);
                pivot = row;
            }
        }
        if (largest <= scale * 1e-14)
            return 0;
        if (pivot != column) {
            int j;
            for (j = column; j < 4; ++j) {
                double swap = augmented[column][j];
                augmented[column][j] = augmented[pivot][j];
                augmented[pivot][j] = swap;
            }
        }
        {
            double divisor = augmented[column][column];
            int j;
            for (j = column; j < 4; ++j)
                augmented[column][j] /= divisor;
        }
        for (row = 0; row < 3; ++row) {
            double factor;
            int j;
            if (row == column)
                continue;
            factor = augmented[row][column];
            for (j = column; j < 4; ++j)
                augmented[row][j] -= factor * augmented[column][j];
        }
    }
    for (row = 0; row < 3; ++row)
        result[row] = augmented[row][3];
    return 1;
}

static int normal_equations(const double *x, const double *y, size_t count,
    const double parameters[3], double normal[3][3], double gradient[3],
    double *sse)
{
    size_t i;
    int row, column;
    memset(normal, 0, 9 * sizeof(double));
    memset(gradient, 0, 3 * sizeof(double));
    *sse = 0.0;
    for (i = 0; i < count; ++i) {
        double h = parameters[0], slope = parameters[1], focus = parameters[2];
        double offset = x[i] - focus, model = model_value(h, slope, focus, x[i]);
        double residual = y[i] - model;
        double jacobian[3];
        if (!(model > 0.0) || !isfinite(model) || !isfinite(residual))
            return 0;
        jacobian[0] = h / model;
        jacobian[1] = slope * offset * offset / model;
        jacobian[2] = -slope * slope * offset / model;
        *sse += residual * residual;
        for (row = 0; row < 3; ++row) {
            gradient[row] += jacobian[row] * residual;
            for (column = 0; column < 3; ++column)
                normal[row][column] += jacobian[row] * jacobian[column];
        }
    }
    return isfinite(*sse);
}

static int nonlinear_fit(const double *x, const double *y, size_t count,
    double minimum_x, double maximum_x, double initial_focus,
    double parameters[3], double *sse)
{
    double normal[3][3], gradient[3], damping = FIT_INITIAL_DAMPING;
    int iteration;
    parameters[2] = clamp_double(initial_focus, minimum_x, maximum_x);
    if (!normal_equations(x, y, count, parameters, normal, gradient, sse))
        return 0;
    for (iteration = 0; iteration < FIT_MAX_ITERATIONS; ++iteration) {
        double damped[3][3], step[3], candidate[3], candidate_sse;
        double candidate_normal[3][3], candidate_gradient[3];
        double relative_step = 0.0;
        int i, j;
        memcpy(damped, normal, sizeof(damped));
        for (i = 0; i < 3; ++i)
            damped[i][i] += damping * fmax(normal[i][i], 1e-12);
        if (!solve3(damped, gradient, step))
            return 0;
        for (i = 0; i < 3; ++i) {
            candidate[i] = parameters[i] + step[i];
            relative_step = fmax(relative_step,
                fabs(step[i]) / fmax(1.0, fabs(parameters[i])));
        }
        if (!(candidate[0] > 0.0) || !(candidate[1] > 0.0) ||
            !isfinite(candidate[0]) || !isfinite(candidate[1]))
            candidate_sse = HUGE_VAL;
        else {
            candidate[2] = clamp_double(candidate[2], minimum_x, maximum_x);
            if (!normal_equations(x, y, count, candidate, candidate_normal,
                candidate_gradient, &candidate_sse))
                candidate_sse = HUGE_VAL;
        }
        if (candidate_sse < *sse) {
            double improvement = *sse - candidate_sse;
            memcpy(parameters, candidate, 3 * sizeof(double));
            *sse = candidate_sse;
            memcpy(normal, candidate_normal, sizeof(normal));
            memcpy(gradient, candidate_gradient, sizeof(gradient));
            damping = fmax(damping / 3.0, 1e-12);
            if (relative_step <= FIT_RELATIVE_TOLERANCE ||
                improvement <= FIT_RELATIVE_TOLERANCE * fmax(1.0, *sse))
                break;
        } else {
            damping *= 10.0;
            if (damping > 1e16)
                break;
        }
        for (j = 0; j < 3; ++j)
            if (!isfinite(gradient[j]))
                return 0;
    }
    return isfinite(*sse);
}

static int inverse3(const double matrix[3][3], double inverse[3][3])
{
    int column, row;
    for (column = 0; column < 3; ++column) {
        double copy[3][3], rhs[3] = {0.0, 0.0, 0.0}, result[3];
        memcpy(copy, matrix, sizeof(copy));
        rhs[column] = 1.0;
        if (!solve3(copy, rhs, result))
            return 0;
        for (row = 0; row < 3; ++row)
            inverse[row][column] = result[row];
    }
    return 1;
}

typedef struct CurveSample {
    double x, y;
} CurveSample;

static int sample_order(const void *left, const void *right)
{
    double a = ((const CurveSample *)left)->x;
    double b = ((const CurveSample *)right)->x;
    return (a > b) - (a < b);
}

static int fit_one_curve(const double *observed_x, const double *observed_y,
    size_t count, FitResult *result)
{
    double *storage, *x, *y, parameters[3], best_parameters[3];
    double minx, maxx, span, origin, miny = HUGE_VAL, maxy = 0.0;
    double best_sse = HUGE_VAL, variance, total = 0.0, mean = 0.0;
    double normal[3][3], gradient[3], inverse[3][3], residual_sse;
    size_t i, minimum_index = 0;
    int guess, fitted = 0, success = 0;

    if (count < MIN_FOCUS_POSITIONS)
        return 0;
    if (count > SIZE_MAX / (2 * sizeof(*storage)))
        return 0;
    storage = (double *)malloc(2 * count * sizeof(*storage));
    if (!storage)
        return 0;
    x = storage;
    y = storage + count;
    minx = observed_x[0];
    maxx = observed_x[count - 1];
    if (!(maxx > minx))
        goto done;
    origin = minx / 2.0 + maxx / 2.0;
    span = (maxx - minx) / 2.0;
    if (!(span > 0.0) || !isfinite(span))
        goto done;
    for (i = 0; i < count; ++i) {
        x[i] = (observed_x[i] - origin) / span;
        y[i] = observed_y[i];
        mean += y[i];
        if (y[i] < miny) {
            miny = y[i];
            minimum_index = i;
        }
        if (y[i] > maxy)
            maxy = y[i];
    }
    mean /= count;
    for (i = 0; i < count; ++i) {
        double delta = y[i] - mean;
        total += delta * delta;
    }
    if (!(total > 0.0) || !isfinite(total))
        goto done;
    for (guess = 0; guess < 2; ++guess) {
        double focus0 = guess == 0 ? x[minimum_index] : 0.0;
        double distance = fmax(fabs(x[0] - focus0), fabs(x[count - 1] - focus0));
        double h = fmax(miny * 0.9, DBL_MIN);
        double slope = sqrt(fmax(maxy * maxy - h * h, 0.0)) /
            fmax(distance, 0.1);
        double sse;
        parameters[0] = h;
        parameters[1] = fmax(slope, h * 1e-3);
        if (nonlinear_fit(x, y, count, x[0], x[count - 1], focus0,
            parameters, &sse) && sse < best_sse) {
            memcpy(best_parameters, parameters, sizeof(best_parameters));
            best_sse = sse;
            fitted = 1;
        }
    }
    if (!fitted || !(best_parameters[0] > 0.0) ||
        !(best_parameters[1] > 0.0) ||
        best_parameters[2] <= x[0] || best_parameters[2] >= x[count - 1] ||
        !isfinite(best_sse) || 1.0 - best_sse / total < MIN_FIT_R_SQUARED)
        goto done;
    if (!normal_equations(x, y, count, best_parameters, normal, gradient,
        &residual_sse) || !inverse3(normal, inverse))
        goto done;
    variance = residual_sse / (count - 3) * inverse[2][2];
    if (!(variance >= 0.0) || !isfinite(variance))
        goto done;
    result->focus = origin + best_parameters[2] * span;
    result->uncertainty = sqrt(variance) * span;
    success = isfinite(result->focus) && isfinite(result->uncertainty);
done:
    free(storage);
    return success;
}

static int fit_identity(int identity, FitResult *result)
{
    CurveSample *samples;
    double *storage, *positions, *values;
    size_t count = 0, combined = 0, i, j;
    int fitted;

    for (i = 0; i < run.count; ++i) {
        if (run.points[i].hfd[identity] > 0.0)
            ++count;
    }
    if (count < MIN_FOCUS_POSITIONS ||
        count > SIZE_MAX / sizeof(*samples) ||
        count > SIZE_MAX / (2 * sizeof(*storage)))
        return 0;
    samples = (CurveSample *)malloc(count * sizeof(*samples));
    storage = (double *)malloc(2 * count * sizeof(*storage));
    if (!samples || !storage) {
        free(samples);
        free(storage);
        return 0;
    }
    positions = storage;
    values = storage + count;
    for (i = 0, j = 0; i < run.count; ++i) {
        if (run.points[i].hfd[identity] > 0.0) {
            samples[j].x = run.points[i].position;
            samples[j].y = run.points[i].hfd[identity];
            ++j;
        }
    }
    qsort(samples, count, sizeof(*samples), sample_order);
    for (i = 0; i < count;) {
        double sum = 0.0;
        size_t repeats = 0;
        j = i;
        while (j < count && samples[j].x == samples[i].x) {
            sum += samples[j++].y;
            ++repeats;
        }
        positions[combined] = samples[i].x;
        values[combined] = sum / repeats;
        ++combined;
        i = j;
    }
    free(samples);
    fitted = fit_one_curve(positions, values, combined, result);
    free(storage);
    return fitted;
}

int CALL finish_computation(int32_t *best_focus, uint32_t *error)
{
    FitResult results[MAX_STARS];
    double focuses[MAX_STARS], deviations[MAX_STARS], uncertainties[MAX_STARS];
    double central, spread, uncertainty, rounded;
    size_t count = 0, i;
    int identity, status = -1;

    if (error)
        *error = UINT32_MAX;
    if (!run.active)
        return -1;
    if (run.trial.valid && commit_trial() != 0)
        goto done;
    if (!best_focus || !error || !run.selected)
        goto done;
    for (identity = 0; identity < run.tracked && count < (size_t)run.limit;
        ++identity) {
        if (fit_identity(identity, &results[count])) {
            focuses[count] = results[count].focus;
            uncertainties[count] = results[count].uncertainty;
            ++count;
        }
    }
    if (count == 0)
        goto done;
    central = median(focuses, count);
    for (i = 0; i < count; ++i)
        deviations[i] = fabs(focuses[i] - central);
    spread = 1.4826 * median(deviations, count);
    uncertainty = ceil(fmax(1.0, fmax(median(uncertainties, count), spread)));
    rounded = round(central);
    if (!isfinite(rounded) || rounded < INT32_MIN || rounded > INT32_MAX ||
        !isfinite(uncertainty) || uncertainty > UINT32_MAX)
        goto done;
    *best_focus = (int32_t)rounded;
    *error = (uint32_t)uncertainty;
    status = 0;
done:
    reset();
    return status;
}
