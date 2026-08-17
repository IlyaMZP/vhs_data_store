/*
 * decode.c – Decode a data file from a PAL composite video signal (8-bit raw samples).
 *
 * Streaming version: reads continuously from stdin (or a file), processes in
 * overlapping windows, and stops once the MAGIC header has been found and the
 * declared payload length has been written.
 *
 * Usage: ./decode [<input.raw>] <output_file>
 *        (omit input or use "-" to read from stdin)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <errno.h>
#include <zlib.h>

#include "rs_codec.h"
#include "vhs_params.h"

typedef enum {
    MODE_DATA,
    MODE_AUDIO
} Mode;

Mode mode = MODE_DATA;

/* ------------------------------------------------------------------ */
/* DPLL tracking gains (early-late gate)                               */
/* ------------------------------------------------------------------ */
#define ALPHA  0.12f   /* phase */
#define BETA   0.015f  /* frequency */

/* Streaming parameters */
#define WINDOW_SAMPLES     2000000u   /* ~50 ms @ 40 MHz */
#define OVERLAP_SAMPLES    200000u

/* ------------------------------------------------------------------ */
/* Levels                                                              */
/* ------------------------------------------------------------------ */
typedef struct {
    float sync;
    float white;
    float sync_threshold;
    float bit_threshold;
} Levels;

static int cmp_float(const void *a, const void *b)
{
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

/* Exact percentile (linear interpolation), O(n log n).
 * Called only once during level calibration. */
static float percentile(const float *samples, size_t n, float p)
{
    if (n == 0)
        return 0.0f;

    float *tmp = malloc(n * sizeof(float));
    if (!tmp)
        return 0.0f;
    memcpy(tmp, samples, n * sizeof(float));

    qsort(tmp, n, sizeof(float), cmp_float);

    float rank = (p / 100.0f) * (n - 1);
    size_t lo  = (size_t)floorf(rank);
    size_t hi  = (size_t)ceilf(rank);
    float  frac = rank - (float)lo;
    float  res  = tmp[lo] * (1.0f - frac) + tmp[hi] * frac;

    free(tmp);
    return res;
}

static void levels_init(Levels *L, const float *samples, size_t n)
{
    L->sync  = percentile(samples, n, 0.2f);
    L->white = percentile(samples, n, 99.8f);
    float span = L->white - L->sync;
    L->sync_threshold = L->sync + span * 0.10f;
    L->bit_threshold  = L->sync + span * 0.65f;
}

/* ------------------------------------------------------------------ */
/* Pulse detection / cleaning / classification                         */
/* ------------------------------------------------------------------ */

/* Returns number of pulses written into starts[]/ends[].
 * Caller must free the arrays. */
static size_t detect_pulses(const float *samples, size_t n,
                            float threshold,
                            uint32_t **starts_out, uint32_t **ends_out)
{
    if (n < 2) {
        *starts_out = NULL;
        *ends_out = NULL;
        return 0;
    }

    /* first pass: count transitions */
    size_t n_start = 0, n_end = 0;
    int prev_below = samples[0] < threshold;
    for (size_t i = 1; i < n; i++) {
        int below = samples[i] < threshold;
        if (below && !prev_below) n_start++;
        if (!below && prev_below) n_end++;
        prev_below = below;
    }

    uint32_t *starts = malloc(n_start * sizeof(uint32_t));
    uint32_t *ends   = malloc(n_end   * sizeof(uint32_t));
    if (!starts || !ends) {
        free(starts); free(ends);
        *starts_out = NULL; *ends_out = NULL;
        return 0;
    }

    size_t si = 0, ei = 0;
    prev_below = samples[0] < threshold;
    for (size_t i = 1; i < n; i++) {
        int below = samples[i] < threshold;
        if (below && !prev_below) starts[si++] = (uint32_t)i;
        if (!below && prev_below) ends[ei++]   = (uint32_t)i;
        prev_below = below;
    }

    /* drop leading end that precedes first start */
    size_t e_off = 0;
    if (ei > 0 && si > 0 && ends[0] <= starts[0])
        e_off = 1;

    size_t n_pulses = si < (ei - e_off) ? si : (ei - e_off);
    if (e_off) {
        memmove(ends, ends + e_off, n_pulses * sizeof(uint32_t));
    }

    *starts_out = starts;
    *ends_out   = ends;
    return n_pulses;
}

/* Merge pulses split by noise (gap < 1 us) and drop spikes < 1.2 us.
 * Returns new count; arrays are reallocated / updated in place. */
static size_t clean_pulses(uint32_t **starts, uint32_t **ends, size_t n,
                           const Geometry *g)
{
    if (n == 0) return 0;

    uint32_t gap_min = us_to_samples(g->sr, 1.0);
    uint32_t width_min = us_to_samples(g->sr, 1.2);

    /* merge */
    uint32_t *s = malloc(n * sizeof(uint32_t));
    uint32_t *e = malloc(n * sizeof(uint32_t));
    if (!s || !e) { free(s); free(e); return 0; }

    size_t m = 0;
    s[0] = (*starts)[0];
    e[0] = (*ends)[0];
    for (size_t i = 1; i < n; i++) {
        if ((*starts)[i] - e[m] < gap_min) {
            /* merge into current */
            e[m] = (*ends)[i];
        } else {
            m++;
            s[m] = (*starts)[i];
            e[m] = (*ends)[i];
        }
    }
    m++; /* number of merged segments */

    /* drop narrow spikes */
    size_t k = 0;
    for (size_t i = 0; i < m; i++) {
        if (e[i] - s[i] >= width_min) {
            s[k] = s[i];
            e[k] = e[i];
            k++;
        }
    }

    free(*starts);
    free(*ends);
    *starts = s;
    *ends   = e;
    return k;
}

/* 0 = short, 1 = hsync, 2 = broad */
static void classify(const uint32_t *starts, const uint32_t *ends, size_t n,
                     const Geometry *g, int8_t *kinds)
{
    float mid_eq_hs = (g->eq + g->hsync) * 0.5f;
    float mid_hs_br = (g->hsync + g->broad) * 0.5f;
    for (size_t i = 0; i < n; i++) {
        float w = (float)(ends[i] - starts[i]);
        if (w < mid_eq_hs)
            kinds[i] = 0;
        else if (w >= mid_hs_br)
            kinds[i] = 2;
        else
            kinds[i] = 1;
    }
}

/* ------------------------------------------------------------------ */
/* Line decoder (TBC + preamble lock + early-late DPLL)                */
/* ------------------------------------------------------------------ */

/* Returns 1 on success and writes LINE_BYTES into out[]; 0 on failure. */
static int decode_line(const float *samples, size_t n_samples,
                       uint32_t start, float ratio,
                       const Geometry *g, const Levels *levels,
                       uint8_t *out)
{
    int32_t ext_off = (int32_t)lroundf((g->active_off - g->pad) * ratio);
    int32_t ext_len = (int32_t)lroundf((g->pad + g->active) * ratio);
    if (ext_len < 10) return 0;

    int64_t a = (int64_t)start + ext_off;
    int64_t b = a + ext_len;
    if (a < 0 || b > (int64_t)n_samples) return 0;

    const float *raw = samples + a;
    size_t raw_len = (size_t)ext_len;
    size_t grid_len = g->pad + g->active;

    /* linear resample to nominal grid (per-line TBC) */
    float *tbc = malloc(grid_len * sizeof(float));
    if (!tbc) return 0;
    for (size_t i = 0; i < grid_len; i++) {
        float t = (float)i / (float)(grid_len - 1);          /* 0..1 */
        float src = t * (raw_len - 1);
        size_t lo = (size_t)floorf(src);
        size_t hi = lo + 1 < raw_len ? lo + 1 : lo;
        float frac = src - lo;
        tbc[i] = raw[lo] * (1.0f - frac) + raw[hi] * frac;
    }

    /* light 3-tap box smooth */
    float *sm = malloc(grid_len * sizeof(float));
    if (!sm) { free(tbc); return 0; }
    for (size_t i = 0; i < grid_len; i++) {
        float s = tbc[i];
        if (i > 0) s += tbc[i - 1];
        if (i + 1 < grid_len) s += tbc[i + 1];
        int cnt = 1 + (i > 0) + (i + 1 < grid_len);
        sm[i] = s / (float)cnt;
    }
    free(tbc);

    /* binarize + edge positions (mid-sample) */
    uint8_t *binar = malloc(grid_len);
    if (!binar) { free(sm); return 0; }
    for (size_t i = 0; i < grid_len; i++)
        binar[i] = sm[i] > levels->bit_threshold ? 1 : 0;

    float *edges = malloc(grid_len * sizeof(float));
    if (!edges) { free(sm); free(binar); return 0; }
    size_t n_edges = 0;
    for (size_t i = 1; i < grid_len; i++) {
        if (binar[i] != binar[i - 1])
            edges[n_edges++] = (float)i - 0.5f;
    }
    free(binar);

    float spb = (float)g->nominal_spb;

    /* preamble edges near g.pad */
    float pre_lim = g->pad + (PREAMBLE_BITS + 0.6f) * spb;
    size_t n_pre = 0;
    float *pre = malloc(n_edges * sizeof(float));
    if (!pre) { free(sm); free(edges); return 0; }
    for (size_t i = 0; i < n_edges; i++) {
        if (edges[i] < pre_lim)
            pre[n_pre++] = edges[i];
    }
    if (n_pre < 10) {
        free(sm); free(edges); free(pre);
        return 0;
    }

    /* phase estimate from median residual */
    float *resid = malloc(n_pre * sizeof(float));
    if (!resid) { free(sm); free(edges); free(pre); return 0; }
    for (size_t i = 0; i < n_pre; i++) {
        float k = roundf((pre[i] - g->pad) / spb);
        resid[i] = pre[i] - g->pad - k * spb;
    }
    /* median */
    for (size_t i = 1; i < n_pre; i++) { /* insertion sort */
        float key = resid[i];
        size_t j = i;
        while (j > 0 && resid[j - 1] > key) {
            resid[j] = resid[j - 1];
            j--;
        }
        resid[j] = key;
    }
    float med_resid = resid[n_pre / 2];
    float e0 = g->pad + med_resid;

    if (fabsf(e0 - g->pad) > us_to_samples(g->sr, 0.5)) {
        free(sm); free(edges); free(pre); free(resid);
        return 0;
    }

    /* re-check spacing against refined phase */
    size_t unique_k = 0;
    float abs_resid_sum = 0.0f;
    int seen[64] = {0}; /* enough for preamble */
    for (size_t i = 0; i < n_pre; i++) {
        float k = roundf((pre[i] - e0) / spb);
        resid[i] = pre[i] - e0 - k * spb;
        abs_resid_sum += fabsf(resid[i]);
        int ki = (int)k;
        if (ki >= 0 && ki < 64 && !seen[ki]) {
            seen[ki] = 1;
            unique_k++;
        }
    }
    float med_abs = abs_resid_sum / n_pre; /* approx; exact median not critical */
    /* better median of abs */
    for (size_t i = 0; i < n_pre; i++) resid[i] = fabsf(resid[i]);
    for (size_t i = 1; i < n_pre; i++) {
        float key = resid[i];
        size_t j = i;
        while (j > 0 && resid[j - 1] > key) {
            resid[j] = resid[j - 1];
            j--;
        }
        resid[j] = key;
    }
    med_abs = resid[n_pre / 2];

    free(pre); free(resid);
    if (med_abs > 0.30f * spb || unique_k < 10) {
        free(sm); free(edges);
        return 0;
    }

    /* early-late DPLL for data bits */
    float center = e0 + (PREAMBLE_BITS + 0.5f) * spb;
    uint8_t bits[DATA_BITS_PER_LINE];

    for (int i = 0; i < DATA_BITS_PER_LINE; i++) {
        int lo = (int)lroundf(center - 0.25f * spb);
        int hi = (int)lroundf(center + 0.25f * spb);
        if (lo < 0) lo = 0;
        if (hi >= (int)grid_len) hi = (int)grid_len - 1;

        float sum = 0.0f;
        int cnt = 0;
        for (int x = lo; x <= hi; x++) {
            sum += sm[x];
            cnt++;
        }
        bits[i] = (sum / cnt) > levels->bit_threshold ? 1 : 0;

        float expected_edge = center + 0.5f * spb;
        center += spb;

        /* find edges in [expected - 0.35*spb, expected + 0.35*spb] */
        float win_lo = expected_edge - 0.35f * spb;
        float win_hi = expected_edge + 0.35f * spb;
        float best_err = 0.0f;
        int found = 0;
        for (size_t e = 0; e < n_edges; e++) {
            if (edges[e] >= win_lo && edges[e] <= win_hi) {
                float err = edges[e] - expected_edge;
                if (!found || fabsf(err) < fabsf(best_err)) {
                    best_err = err;
                    found = 1;
                }
            }
        }
        if (found) {
            center += ALPHA * best_err;
            spb    += BETA  * best_err;
        }
    }

    free(sm); free(edges);

    if (!(0.9f * g->nominal_spb < spb && spb < 1.1f * g->nominal_spb))
        return 0;

    /* pack bits MSB first (matches np.packbits default) */
    for (int i = 0; i < LINE_BYTES; i++) {
        uint8_t b = 0;
        for (int j = 0; j < 8; j++)
            b = (b << 1) | bits[i * 8 + j];
        out[i] = b;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Field splitting / parity / decode                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t lo, hi;   /* pulse index range [lo, hi) of the field body */
} FieldRange;

/* Collect field body ranges. Returns count; caller frees the array. */
static size_t split_fields(const int8_t *kinds, size_t n, FieldRange **out)
{
    FieldRange *fields = NULL;
    size_t n_fields = 0, cap = 0;

    size_t i = 0;
    while (i + 4 < n) {
        if (kinds[i] == 2 && kinds[i+1] == 2 && kinds[i+2] == 2 &&
            kinds[i+3] == 2 && kinds[i+4] == 2) {
            size_t j = i;
            while (j < n && kinds[j] == 2) j++;

            size_t k = j;
            while (k + 4 < n) {
                if (kinds[k] == 2 && kinds[k+1] == 2 && kinds[k+2] == 2 &&
                    kinds[k+3] == 2 && kinds[k+4] == 2)
                    break;
                k++;
            }
            if (k + 4 >= n) k = n;

            if (n_fields >= cap) {
                cap = cap ? cap * 2 : 8;
                fields = realloc(fields, cap * sizeof(FieldRange));
            }
            fields[n_fields].lo = j;
            fields[n_fields].hi = k;
            n_fields++;
            i = k;
        } else {
            i++;
        }
    }
    *out = fields;
    return n_fields;
}

/* Returns VBI offset: 17 (odd) or 16 (even) */
static int detect_field_parity(const uint32_t *starts, const int8_t *kinds,
                               size_t lo, const Geometry *g)
{
    /* last eq/broad before lo */
    int prev = -1;
    size_t from = lo > 15 ? lo - 15 : 0;
    for (size_t i = from; i < lo; i++) {
        if (kinds[i] == 0 || kinds[i] == 2)
            prev = (int)i;
    }
    if (prev < 0) return VBI_LINES;

    uint32_t gap = starts[lo] - starts[prev];
    if (gap < (uint32_t)(0.75 * g->spl))
        return VBI_LINES;       /* odd */
    else
        return VBI_LINES - 1;   /* even */
}

/* De-interleave: stream[FIELD_STREAM_BYTES] -> codewords[CODEWORDS][RS_N]
 * Only the first CODED_BYTES of stream are used. */
static void deinterleave(const uint8_t *stream, uint8_t codewords[CODEWORDS][RS_N])
{
    for (size_t s = 0; s < CODED_BYTES; s++) {
        size_t cw  = s % CODEWORDS;
        size_t sym = s / CODEWORDS;
        codewords[cw][sym] = stream[s];
    }
}

typedef struct {
    size_t lines_total;
    size_t lines_ok;
    size_t symbols_corrected;
    size_t codewords_failed;
} Stats;

/* Decode one field. Returns 1 on success, 0 on failure.
 * The payload buffer must be allocated by the caller and be FIELD_PAYLOAD bytes. */
static int decode_field(const float *samples, size_t n_samples,
                        const uint32_t *starts, const int8_t *kinds,
                        size_t lo, size_t hi,
                        const Geometry *g, const Levels *levels,
                        Stats *stats,
                        uint8_t *payload)  /* pre-allocated buffer */
{
    /* collect hsync indices */
    size_t *hs_idx = malloc((hi - lo) * sizeof(size_t));
    if (!hs_idx) return 0;
    size_t n_hs = 0;
    for (size_t i = lo; i < hi; i++) {
        if (kinds[i] == 1)
            hs_idx[n_hs++] = i;
    }
    if (n_hs < DATA_LINES / 2) {
        free(hs_idx);
        return 0;
    }

    int vbi_offset = detect_field_parity(starts, kinds, hs_idx[0], g);

    /* median line period from plausible diffs */
    double P = (double)g->spl;
    if (n_hs >= 2) {
        double *diffs = malloc((n_hs - 1) * sizeof(double));
        size_t n_plaus = 0;
        for (size_t j = 1; j < n_hs; j++) {
            double d = (double)starts[hs_idx[j]] - starts[hs_idx[j-1]];
            if (d > 0.9 * g->spl && d < 1.1 * g->spl)
                diffs[n_plaus++] = d;
        }
        if (n_plaus) {
            /* median */
            for (size_t i = 1; i < n_plaus; i++) {
                double key = diffs[i];
                size_t j = i;
                while (j > 0 && diffs[j-1] > key) {
                    diffs[j] = diffs[j-1];
                    j--;
                }
                diffs[j] = key;
            }
            P = diffs[n_plaus / 2];
        }
        free(diffs);
    }

    uint8_t stream[FIELD_STREAM_BYTES];
    memset(stream, 0, sizeof(stream));
    uint8_t line_ok[DATA_LINES];
    memset(line_ok, 0, sizeof(line_ok));

    int n = 0;
    uint32_t prev = starts[hs_idx[0]];

    for (size_t j = 0; j < n_hs; j++) {
        uint32_t s = starts[hs_idx[j]];
        if (j > 0) {
            double d = (double)s - prev;
            if (d < 0.5 * P)
                continue; /* spurious */
            int k = (int)lround(d / P);
            if (k < 1) k = 1;
            n += k;
            if (fabs(d / k - P) < 0.02 * P)
                P = 0.9 * P + 0.1 * (d / k);
            prev = s;
        }

        int data_idx = n - vbi_offset;
        if (data_idx < 0 || data_idx >= DATA_LINES || line_ok[data_idx])
            continue;

        uint32_t nxt = (j + 1 < n_hs) ? starts[hs_idx[j + 1]]
                                      : (uint32_t)(s + P);
        float ratio = (float)((double)(nxt - s) / g->spl);
        if (!(0.95f < ratio && ratio < 1.05f))
            ratio = (float)(P / g->spl);

        uint8_t line_bytes[LINE_BYTES];
        if (decode_line(samples, n_samples, s, ratio, g, levels, line_bytes)) {
            memcpy(stream + data_idx * LINE_BYTES, line_bytes, LINE_BYTES);
            line_ok[data_idx] = 1;
        }
    }
    free(hs_idx);

    stats->lines_total += DATA_LINES;
    size_t ok_count = 0;
    for (int i = 0; i < DATA_LINES; i++)
        if (line_ok[i]) ok_count++;
    stats->lines_ok += ok_count;

    if (ok_count == 0)
        return 0;

    /* de-whiten */
    for (size_t i = 0; i < FIELD_STREAM_BYTES; i++)
        stream[i] ^= whitening_bytes[i];

    /* de-interleave */
    uint8_t codewords[CODEWORDS][RS_N];
    deinterleave(stream, codewords);

    /* per-symbol validity (line-level erasure) */
    uint8_t byte_ok[FIELD_STREAM_BYTES];
    for (int i = 0; i < DATA_LINES; i++) {
        for (int b = 0; b < LINE_BYTES; b++)
            byte_ok[i * LINE_BYTES + b] = line_ok[i];
    }

    size_t payload_pos = 0;

    for (size_t c = 0; c < CODEWORDS; c++) {
        size_t erase_pos[RS_N];
        size_t erase_count = 0;
        for (size_t sym = 0; sym < RS_N; sym++) {
            size_t stream_idx = c + CODEWORDS * sym;
            if (stream_idx < FIELD_STREAM_BYTES && !byte_ok[stream_idx])
                erase_pos[erase_count++] = sym;
        }

        uint8_t out_data[RS_K];
        size_t corrected = 0;
        rs_status_t st = rs_decode(codewords[c], RS_N, NSYM,
                                   erase_pos, erase_count,
                                   out_data, RS_K, &corrected);
        if (st == RS_OK) {
            stats->symbols_corrected += corrected;
            memcpy(payload + payload_pos, out_data, RS_K);
        } else {
            /* Fallback: output zeroes for uncorrectable Reed-Solomon blocks */
            memset(payload + payload_pos, 0x00, RS_K);
            stats->codewords_failed++;
        }
        payload_pos += RS_K;
    }

    return 1;
}

/* ------------------------------------------------------------------ */
/* Streaming I/O helpers                                               */
/* ------------------------------------------------------------------ */

static size_t read_chunk(FILE *f, uint8_t *buf, size_t n)
{
    return fread(buf, 1, n, f);
}

static const uint8_t SILENCE[FIELD_PAYLOAD] = {0};

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    const char *input_path = "-";
    const char *output_path = NULL;

    int data_set = 0;
    int audio_set = 0;

    // 1. Parse flags
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--data") == 0) {
            data_set = 1;
        } else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--audio") == 0) {
            audio_set = 1;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            fprintf(stderr, "Usage: %s [-d|--data | -a|--audio] [<input.raw>] <output_file>\n", argv[0]);
            return 1;
        }
        i++;
    }

    // 2. Validate flags and set mode
    if (data_set && audio_set) {
        fprintf(stderr, "Error: cannot specify both -d/--data and -a/--audio\n");
        return 1;
    }
    mode = audio_set ? MODE_AUDIO : MODE_DATA;

    // 3. Handle remaining positional arguments (the last one or two)
    int remaining = argc - i;
    if (remaining == 1) {
        output_path = argv[argc - 1];
    } else if (remaining == 2) {
        input_path  = argv[argc - 2];
        output_path = argv[argc - 1];
    } else {
        fprintf(stderr, "Usage: %s [-d|--data | -a|--audio] [<input.raw>] <output_file>\n", argv[0]);
        return 1;
    }

    Geometry g;
    geometry_init(&g, DEFAULT_SAMPLE_RATE);

    FILE *infile;
    if (strcmp(input_path, "-") == 0) {
        infile = stdin;
        fprintf(stderr, "Reading from stdin …\n");
        fflush(stderr);
    } else {
        infile = fopen(input_path, "rb");
        if (!infile) {
            perror(input_path);
            return 1;
        }
        fprintf(stderr, "Reading from %s\n", input_path);
        fflush(stderr);
    }

    FILE *outf = fopen(output_path, "wb");
    if (!outf) {
        perror(output_path);
        if (infile != stdin) fclose(infile);
        return 1;
    }

    /* rolling sample buffer (float) */
    size_t buf_cap = WINDOW_SAMPLES + OVERLAP_SAMPLES + 100000;
    float *buf = malloc(buf_cap * sizeof(float));
    if (!buf) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    size_t buf_len = 0;

    Levels levels;
    int levels_ready = 0;
    Stats stats = {0};
    uint8_t *blob = NULL;
    size_t blob_len = 0, blob_cap = 0;
    int header_found = 0;
    uint32_t payload_length = 0;
    uint32_t expected_crc = 0;
    size_t written = 0;
    size_t field_count = 0;
    int64_t last_processed_end = 0; /* absolute sample index */

    size_t LEVEL_CALIB_SAMPLES = (size_t)g.sr; /* ≈ 1 s */

    uint8_t *payload = malloc(FIELD_PAYLOAD);
    if (!payload) {
        fprintf(stderr, "OOM\n");
        return 1;
    }

    int done = 0;

    while (!done) {
        /* ---- fill rolling buffer ---- */
        size_t need;
        if (buf_len > OVERLAP_SAMPLES)
            need = WINDOW_SAMPLES - (buf_len - OVERLAP_SAMPLES);
        else
            need = WINDOW_SAMPLES;
        if (need <= 0) need = WINDOW_SAMPLES / 2;
        if (need > buf_cap - buf_len) {
            /* grow */
            buf_cap = buf_len + need + 100000;
            float *nb = realloc(buf, buf_cap * sizeof(float));
            if (!nb) { fprintf(stderr, "OOM\n"); break; }
            buf = nb;
        }
        uint8_t *raw_chunk = malloc(need);
        if (!raw_chunk) { fprintf(stderr, "OOM\n"); break; }
        size_t nread = read_chunk(infile, raw_chunk, need);
        if (nread == 0 && buf_len == 0) {
            fprintf(stderr, "EOF with no data\n");
            free(raw_chunk);
            break;
        }
        if (nread) {
            for (size_t i = 0; i < nread; i++)
                buf[buf_len++] = (float)raw_chunk[i];
        }
        free(raw_chunk);

        /* ---- calibrate levels once ---- */
        if (!levels_ready) {
            size_t need_calib = LEVEL_CALIB_SAMPLES < 100000 ? LEVEL_CALIB_SAMPLES : 100000;
            if (buf_len < need_calib) {
                if (nread == 0) {
                    fprintf(stderr, "EOF before enough samples for level calibration\n");
                    break;
                }
                continue;
            }
            size_t calib_n = buf_len < LEVEL_CALIB_SAMPLES ? buf_len : LEVEL_CALIB_SAMPLES;
            levels_init(&levels, buf, calib_n);
            levels_ready = 1;
            fprintf(stderr,
                    "Levels (first %.2fs): sync=%.1f white=%.1f "
                    "sync_th=%.1f bit_th=%.1f\n",
                    (double)calib_n / g.sr,
                    levels.sync, levels.white,
                    levels.sync_threshold, levels.bit_threshold);
            fflush(stderr);
        }
        /* ---- detect pulses ---- */
        uint32_t *starts = NULL, *ends = NULL;
        size_t n_pulses = detect_pulses(buf, buf_len, levels.sync_threshold,
                                        &starts, &ends);
        n_pulses = clean_pulses(&starts, &ends, n_pulses, &g);

        if (n_pulses == 0) {
            free(starts); free(ends);
            if (nread == 0) break;
            /* discard oldest part */
            if (buf_len > WINDOW_SAMPLES + OVERLAP_SAMPLES) {
                size_t drop = buf_len - WINDOW_SAMPLES;
                memmove(buf, buf + drop, (buf_len - drop) * sizeof(float));
                buf_len -= drop;
                last_processed_end = last_processed_end > (int64_t)drop
                                     ? last_processed_end - drop : 0;
            }
            continue;
        }

        int8_t *kinds = malloc(n_pulses * sizeof(int8_t));
        if (!kinds) {
            free(starts); free(ends);
            break;
        }
        classify(starts, ends, n_pulses, &g, kinds);

        FieldRange *fields = NULL;
        size_t n_fields = split_fields(kinds, n_pulses, &fields);

        size_t safe_end = buf_len > OVERLAP_SAMPLES
                          ? buf_len - OVERLAP_SAMPLES : buf_len;

        for (size_t f = 0; f < n_fields; f++) {
            size_t lo = fields[f].lo;
            size_t hi = fields[f].hi;
            if (lo >= n_pulses) continue;

            uint32_t field_start_sample = starts[lo];
            if ((int64_t)field_start_sample < last_processed_end)
                continue; /* already handled */

            size_t last_idx = hi > 0 ? hi - 1 : 0;
            if (last_idx >= n_pulses) last_idx = n_pulses - 1;
            if (starts[last_idx] > safe_end)
                continue; /* not fully inside safe region */

            field_count++;
            fprintf(stderr, "Field %zu\n", field_count);
            fflush(stderr);

            int success = decode_field(buf, buf_len,
                                      starts, kinds, lo, hi,
                                      &g, &levels, &stats,
                                      payload);

            if (mode == MODE_AUDIO) {
                if (success) {
                    fwrite(payload, 1, FIELD_PAYLOAD, outf);
                } else {
                    /* Zero out only on whole field failure */
                    memset(payload, 0, FIELD_PAYLOAD);
                    fwrite(payload, 1, FIELD_PAYLOAD, outf);
                }
                fflush(outf);
                free(payload);
            } else {
                if (!success) {
                    fprintf(stderr, "Field %zu failed to decode\n", field_count);
                    memset(payload, 0, FIELD_PAYLOAD);
                }

                /* append to blob */
                if (blob_len + FIELD_PAYLOAD > blob_cap) {
                    blob_cap = (blob_len + FIELD_PAYLOAD) * 2 + 4096;
                    uint8_t *nb = realloc(blob, blob_cap);
                    if (!nb) {
                        fprintf(stderr, "OOM\n");
                        done = 1;
                        break;
                    }
                    blob = nb;
                }
                memcpy(blob + blob_len, payload, FIELD_PAYLOAD);
                blob_len += FIELD_PAYLOAD;
            }

            last_processed_end = (int64_t)starts[last_idx] + g.spl;

            /* ---- header / payload streaming ---- */
            if (!header_found) {
                /* search for MAGIC */
                size_t magic_offset = (size_t)-1;
                for (size_t i = 0; i + 4 <= blob_len; i++) {
                    if (memcmp(blob + i, MAGIC, 4) == 0) {
                        magic_offset = i;
                        break;
                    }
                }
                if (magic_offset != (size_t)-1) {
                    if (magic_offset > 0) {
                        fprintf(stderr,
                                "Notice: Discarded %zu bytes of leading garbage before MAGIC.\n",
                                magic_offset);
                        fflush(stderr);
                        memmove(blob, blob + magic_offset, blob_len - magic_offset);
                        blob_len -= magic_offset;
                    }
                    if (blob_len >= HEADER_LEN) {
                        memcpy(&payload_length, blob + 4, 4);
                        memcpy(&expected_crc,   blob + 8, 4);
                        /* little-endian already if host is LE; otherwise swap */
                        /* assume LE host for simplicity (common for this use) */
                        header_found = 1;
                        fprintf(stderr, "Header found: length=%u  crc=%08x\n",
                                payload_length, expected_crc);
                        fflush(stderr);

                        size_t data_len = blob_len - HEADER_LEN;
                        size_t to_write = data_len < payload_length
                                          ? data_len : payload_length;
                        if (to_write) {
                            fwrite(blob + HEADER_LEN, 1, to_write, outf);
                            fflush(outf);
                            written = to_write;
                            fprintf(stderr, "Flushed %zu payload bytes\n", written);
                            fflush(stderr);
                        }
                        /* keep remaining unwritten payload */
                        size_t remain = data_len - to_write;
                        memmove(blob, blob + HEADER_LEN + to_write, remain);
                        blob_len = remain;
                    }
                }
            } else {
                if (written < payload_length) {
                    size_t needw = payload_length - written;
                    size_t take = blob_len < needw ? blob_len : needw;
                    if (take) {
                        fwrite(blob, 1, take, outf);
                        fflush(outf);
                        written += take;
                        memmove(blob, blob + take, blob_len - take);
                        blob_len -= take;
                        fprintf(stderr, "Flushed %zu/%u payload bytes\n",
                                written, payload_length);
                        fflush(stderr);
                    }
                }
            }

            if (header_found && written >= payload_length) {
                fprintf(stderr, "Payload complete – stopping.\n");
                fflush(stderr);
                done = 1;
                break;
            }
        }

        free(fields);
        free(kinds);
        free(starts);
        free(ends);

        if (done) break;

        /* ---- slide buffer ---- */
        if (nread == 0) {
            /* final pass already attempted via safe_end = whole buffer */
            if (last_processed_end >= (int64_t)buf_len - (int64_t)g.spl * 2)
                break;
            /* force one more iteration with full buffer */
            continue;
        }

        if (buf_len > WINDOW_SAMPLES + OVERLAP_SAMPLES) {
            size_t drop = buf_len - WINDOW_SAMPLES;
            memmove(buf, buf + drop, (buf_len - drop) * sizeof(float));
            buf_len -= drop;
            last_processed_end = last_processed_end > (int64_t)drop
                                 ? last_processed_end - drop : 0;
        }
    }

    free(buf);
    free(payload);
    free(blob);
    if (infile != stdin) fclose(infile);
    fclose(outf);

    /* ------------------------------------------------------------------ */
    /* Final report & CRC                                                  */
    /* ------------------------------------------------------------------ */
    fprintf(stderr,
            "\nLines decoded     : %zu/%zu (%zu erased)\n",
            stats.lines_ok, stats.lines_total,
            stats.lines_total - stats.lines_ok);
    fprintf(stderr, "RS symbols fixed  : %zu\n", stats.symbols_corrected);
    fprintf(stderr, "RS codeword fails : %zu\n", stats.codewords_failed);

    if (!header_found) {
        fprintf(stderr, "ERROR: no valid header found (MAGIC sequence not found in stream)\n");
        return 1;
    }
    if (written < payload_length) {
        fprintf(stderr, "ERROR: truncated payload (%zu/%u bytes written)\n",
                written, payload_length);
        return 1;
    }

    /* re-read written file for CRC */
    FILE *cf = fopen(output_path, "rb");
    if (!cf) {
        perror(output_path);
        return 1;
    }
    fseek(cf, 0, SEEK_END);
    long fsz = ftell(cf);
    fseek(cf, 0, SEEK_SET);
    uint8_t *data = malloc(fsz > 0 ? (size_t)fsz : 1);
    size_t got = fread(data, 1, (size_t)fsz, cf);
    fclose(cf);

    uint32_t actual_crc = crc32(0L, data, (uInt)got) & 0xffffffffu;
    free(data);

    if (actual_crc == expected_crc) {
        fprintf(stderr, "CRC OK — recovered %u bytes -> %s\n",
                payload_length, output_path);
        return 0;
    } else {
        fprintf(stderr,
                "CRC MISMATCH (expected %08x, got %08x) — "
                "wrote %u bytes to %s anyway\n",
                expected_crc, actual_crc, payload_length, output_path);
        return 2;
    }
}
