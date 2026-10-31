#include "focus_curve.h"
#include "analysis_parameters.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

typedef struct FocusPoint {
    int32_t position;
    double hfd[MAX_STARS]; /* zero means this identity was missing */
} FocusPoint;

typedef struct CurvePoint {
    double x, y;
} CurvePoint;

typedef struct Line {
    double x, y, slope, sxx, sse;
} Line;

static struct {
    int active, limit, tracked;
    double x[MAX_STARS], y[MAX_STARS];
    FocusPoint *points;
    size_t count, capacity;
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

static int usable(const StarRecord *star)
{
    return isfinite(star->x) && isfinite(star->y) &&
        isfinite(star->hfd) && star->hfd >= MIN_HFD_PIXELS &&
        star->hfd <= MAX_HFD_PIXELS;
}

/* Mutual nearest neighbours within the configured radius; never share a detection. */
static void match(const StarRecord *stars, int count, int *indices, int *clipped)
{
    int i, j;
    for (i = 0; i < run.tracked; ++i) {
        double best = MATCH_RADIUS_PIXELS * MATCH_RADIUS_PIXELS;
        int selected = -1, tied = 0;
        for (j = 0; j < count; ++j) {
            double dx, dy, distance;
            if (!usable(stars + j))
                continue;
            dx = stars[j].x - run.x[i];
            dy = stars[j].y - run.y[i];
            distance = dx * dx + dy * dy;
            if (stars[j].saturated_pixel_count > MAX_SATURATED_PIXELS) {
                if (distance <= MATCH_RADIUS_PIXELS * MATCH_RADIUS_PIXELS)
                    *clipped = 1;
                continue;
            }
            if (distance < best || (distance == best && selected < 0)) {
                best = distance;
                selected = j;
                tied = 0;
            } else if (distance == best) {
                tied = 1;
            }
        }
        if (selected >= 0 && !tied) {
            for (j = 0; j < run.tracked; ++j) {
                double dx = stars[selected].x - run.x[j];
                double dy = stars[selected].y - run.y[j];
                if (j != i && dx * dx + dy * dy <= best) {
                    selected = -1;
                    break;
                }
            }
        }
        indices[i] = tied ? -1 : selected;
    }
}

/* Only widespread clipping warrants a shorter exposure; dim is not unusable. */
static int exposure_status(uint32_t *exposure)
{
    uint32_t recommended = *exposure / EXPOSURE_REDUCTION_DIVISOR;
    if (recommended < MIN_EXPOSURE_US)
        recommended = MIN_EXPOSURE_US;
    if (recommended >= *exposure)
        return -1;
    *exposure = recommended;
    return 1;
}

int CALL add_focus_point(const StarRecord *stars, int32_t focus_position,
    uint32_t *exposure)
{
    int indices[MAX_STARS], count = 0, tracked, i, matched = 0, clipped = 0;
    FocusPoint point = {0};

    if (!run.active || !stars || !exposure || *exposure == 0)
        return -1;
    while (count < MAX_STARS && stars[count].positive_pixel_count != 0)
        ++count;
    tracked = run.tracked;
    if (tracked == 0) {
        for (i = 0; i < count; ++i) {
            if (!usable(stars + i))
                continue;
            if (stars[i].saturated_pixel_count > MAX_SATURATED_PIXELS) {
                clipped = 1;
                continue;
            }
            if (tracked < run.limit)
                indices[tracked++] = i;
        }
    } else {
        match(stars, count, indices, &clipped);
        /* Fill unused slots when new stars appear; never replace a lost identity. */
        for (i = 0; i < count && tracked < run.limit; ++i) {
            int j;
            if (!usable(stars + i) || stars[i].saturated_pixel_count > MAX_SATURATED_PIXELS)
                continue;
            for (j = 0; j < run.tracked; ++j) {
                double dx = stars[i].x - run.x[j], dy = stars[i].y - run.y[j];
                if (dx * dx + dy * dy <= MATCH_RADIUS_PIXELS * MATCH_RADIUS_PIXELS)
                    break;
            }
            if (j == run.tracked)
                indices[tracked++] = i;
        }
    }
    point.position = focus_position;
    for (i = 0; i < tracked; ++i) {
        const StarRecord *star;
        if (indices[i] < 0)
            continue;
        star = stars + indices[i];
        if (star->saturated_pixel_count > MAX_SATURATED_PIXELS) {
            indices[i] = -1;
            clipped = 1;
            continue;
        }
        point.hfd[i] = star->hfd;
        ++matched;
    }
    if (matched == 0)
        return clipped ? exposure_status(exposure) : -1;
    if (run.count == run.capacity) {
        size_t capacity = run.capacity ? run.capacity * 2 : 16;
        FocusPoint *grown;
        if (capacity < run.capacity || capacity > SIZE_MAX / sizeof(*grown))
            return -1;
        grown = realloc(run.points, capacity * sizeof(*grown));
        if (!grown)
            return -1;
        run.points = grown;
        run.capacity = capacity;
    }
    run.points[run.count++] = point;
    run.tracked = tracked;
    for (i = 0; i < tracked; ++i) {
        if (indices[i] >= 0) {
            run.x[i] = stars[indices[i]].x;
            run.y[i] = stars[indices[i]].y;
        }
    }
    return 0;
}

static int point_order(const void *left, const void *right)
{
    double a = ((const CurvePoint *)left)->x, b = ((const CurvePoint *)right)->x;
    return (a > b) - (a < b);
}

static int number_order(const void *left, const void *right)
{
    double a = *(const double *)left, b = *(const double *)right;
    return (a > b) - (a < b);
}

static double median(double *values, size_t count)
{
    qsort(values, count, sizeof(*values), number_order);
    return count % 2 ? values[count / 2] :
        (values[count / 2 - 1] + values[count / 2]) / 2;
}

static int fit_line(const CurvePoint *points, size_t count, Line *line)
{
    double sxy = 0;
    size_t i;
    memset(line, 0, sizeof(*line));
    for (i = 0; i < count; ++i) {
        line->x += points[i].x;
        line->y += points[i].y;
    }
    line->x /= count;
    line->y /= count;
    for (i = 0; i < count; ++i) {
        double dx = points[i].x - line->x;
        line->sxx += dx * dx;
        sxy += dx * (points[i].y - line->y);
    }
    if (line->sxx == 0)
        return 0;
    line->slope = sxy / line->sxx;
    for (i = 0; i < count; ++i) {
        double residual = points[i].y - line->y -
            line->slope * (points[i].x - line->x);
        line->sse += residual * residual;
    }
    return 1;
}

/* Fit both linear branches; the intersection must be bracketed, not extrapolated. */
static int fit_curve(CurvePoint *points, size_t count, double *focus, double *error)
{
    size_t i, split, distinct = 1;
    double best_sse = HUGE_VAL;
    int fitted = 0;

    if (count < MIN_FOCUS_POSITIONS || count < 2 * MIN_BRANCH_POINTS)
        return 0;
    qsort(points, count, sizeof(*points), point_order);
    for (i = 1; i < count; ++i)
        distinct += points[i].x != points[i - 1].x;
    if (distinct < MIN_FOCUS_POSITIONS)
        return 0;
    for (split = MIN_BRANCH_POINTS; split <= count - MIN_BRANCH_POINTS; ++split) {
        Line left, right;
        double x, y, contrast, sse, variance, dx_left, dx_right, slope_difference;
        if (points[split - 1].x == points[split].x ||
            !fit_line(points, split, &left) ||
            !fit_line(points + split, count - split, &right) ||
            left.slope >= 0 || right.slope <= 0)
            continue;
        slope_difference = left.slope - right.slope;
        x = (right.y - right.slope * right.x -
             left.y + left.slope * left.x) / slope_difference;
        if (x < points[split - 1].x || x > points[split].x)
            continue;
        y = left.y + left.slope * (x - left.x);
        contrast = fmax(MIN_BRANCH_RISE_PIXELS, MIN_BRANCH_RISE_FRACTION * y);
        if (y <= 0 ||
            left.slope * (points[0].x - x) < contrast ||
            right.slope * (points[count - 1].x - x) < contrast)
            continue;
        sse = left.sse + right.sse;
        if (sse >= best_sse)
            continue;
        dx_left = x - left.x;
        dx_right = x - right.x;
        variance = sse / (count - 4) *
            (1.0 / split + dx_left * dx_left / left.sxx +
             1.0 / (count - split) + dx_right * dx_right / right.sxx) /
            (slope_difference * slope_difference);
        if (!isfinite(x) || !isfinite(variance))
            continue;
        *focus = x;
        *error = sqrt(variance);
        best_sse = sse;
        fitted = 1;
    }
    return fitted;
}

int CALL finish_computation(int32_t *best_focus, uint32_t *error)
{
    double focuses[MAX_STARS], errors[MAX_STARS], deviations[MAX_STARS];
    double focus, uncertainty, origin;
    CurvePoint *points;
    size_t i, count, fitted = 0;
    int star, status = -1;

    if (!run.active)
        return -1;
    if (!best_focus || !error || run.count < MIN_FOCUS_POSITIONS)
        goto done;
    points = malloc(run.count * sizeof(*points));
    if (!points)
        goto done;
    origin = run.points[0].position;
    for (star = 0; star < run.tracked; ++star) {
        count = 0;
        for (i = 0; i < run.count; ++i) {
            if (run.points[i].hfd[star] > 0) {
                points[count].x = run.points[i].position - origin;
                points[count++].y = run.points[i].hfd[star];
            }
        }
        if (fit_curve(points, count, focuses + fitted, errors + fitted)) {
            focuses[fitted] += origin;
            ++fitted;
        }
    }
    free(points);
    if (fitted == 0)
        goto done;
    focus = median(focuses, fitted);
    for (i = 0; i < fitted; ++i)
        deviations[i] = fabs(focuses[i] - focus);
    uncertainty = ceil(fmax(1, fmax(median(errors, fitted),
                                    1.4826 * median(deviations, fitted))));
    focus = round(focus);
    if (!isfinite(focus) || focus < INT32_MIN || focus > INT32_MAX ||
        !isfinite(uncertainty) || uncertainty > UINT32_MAX)
        goto done;
    *best_focus = (int32_t)focus;
    *error = (uint32_t)uncertainty;
    status = 0;
done:
    reset();
    return status;
}
