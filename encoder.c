/* encode.c - Encode a data file into a PAL composite video signal (8-bit raw samples).
 *
 * Data is carried in the active (visible) portion of 288 lines per field:
 * preamble bits for decoder PLL lock + DATA_BITS_PER_LINE NRZ data bits.
 * Per field, payload is protected by interleaved RS(255,223) codewords so
 * that head-switching noise wiping out whole lines costs each codeword only
 * a few symbols.
 *
 * Usage: encode <input_file> <output.raw> [--sample-rate 40000000]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <zlib.h>          /* crc32 */

#include "rs_codec.h"
#include "vhs_params.h"

typedef enum {
    MODE_DATA,
    MODE_AUDIO
} Mode;

Mode mode = MODE_DATA;

/* ---------- helpers ---------- */

static void die(const char *msg)
{
    fprintf(stderr, "error: %s\n", msg);
    exit(1);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n);
    if (!p) die("out of memory");
    return p;
}

/* CRC-32 identical to Python zlib.crc32 (unsigned 32-bit) */
static uint32_t crc32_zlib(const uint8_t *data, size_t len)
{
    return crc32(0L, data, (uInt)len) & 0xFFFFFFFFu;
}

/* ---------- field stream construction ---------- */

/*
 * RS-encode + interleave + whiten one field payload -> FIELD_STREAM_BYTES.
 * payload must be at most FIELD_PAYLOAD bytes; it is zero-padded.
 */
static void build_field_stream(const uint8_t *payload, size_t payload_len,
                               uint8_t *out_stream)
{
    uint8_t padded[FIELD_PAYLOAD];
    uint8_t codewords[CODEWORDS][RS_N];
    size_t c, s;

    memset(padded, 0, sizeof padded);
    if (payload_len > FIELD_PAYLOAD)
        payload_len = FIELD_PAYLOAD;
    memcpy(padded, payload, payload_len);

    for (c = 0; c < CODEWORDS; c++) {
        rs_status_t st = rs_encode(padded + c * RS_K, RS_K, NSYM,
                                   codewords[c], RS_N);
        if (st != RS_OK) {
            fprintf(stderr, "rs_encode failed: %s\n", rs_status_string(st));
            exit(1);
        }
    }

    /* interleave: stream[s] = codewords[s % CODEWORDS][s / CODEWORDS]
     * only the first CODED_BYTES are filled; the rest stay zero then
     * are whitened together with the coded part.
     */
    memset(out_stream, 0, FIELD_STREAM_BYTES);
    for (s = 0; s < CODED_BYTES; s++)
        out_stream[s] = codewords[s % CODEWORDS][s / CODEWORDS];

    /* whiten */
    for (s = 0; s < FIELD_STREAM_BYTES; s++)
        out_stream[s] ^= whitening_bytes[s];
}

/* ---------- line templates ---------- */

typedef struct {
    const Geometry *g;
    uint8_t *blank_line;
    uint8_t *black_line;
    uint8_t *short_half;
    uint8_t *broad_half;
    uint8_t *half_line_blank;
    uint8_t *half_line_hsync;
    uint8_t *data_tmpl;
} LineWriter;

static void lw_init(LineWriter *lw, const Geometry *g)
{
    lw->g = g;

    lw->blank_line = xmalloc(g->spl);
    memset(lw->blank_line, U8_BLANK, g->spl);

    lw->black_line = xmalloc(g->spl);
    memcpy(lw->black_line, lw->blank_line, g->spl);
    memset(lw->black_line, U8_SYNC, g->hsync);
    memset(lw->black_line + g->active_off, U8_BLACK, g->active);

    lw->short_half = xmalloc(g->half);
    memset(lw->short_half, U8_BLANK, g->half);
    memset(lw->short_half, U8_SYNC, g->eq);

    lw->broad_half = xmalloc(g->half);
    memset(lw->broad_half, U8_BLANK, g->half);
    memset(lw->broad_half, U8_SYNC, g->broad);

    lw->half_line_blank = xmalloc(g->half);
    memset(lw->half_line_blank, U8_BLANK, g->half);

    lw->half_line_hsync = xmalloc(g->half);
    memset(lw->half_line_hsync, U8_BLANK, g->half);
    memset(lw->half_line_hsync, U8_SYNC, g->hsync);
    /* active region up to -pad (Python: g.active_off:-g.pad) */
    if (g->half > g->active_off + g->pad) {
        size_t len = g->half - g->active_off - g->pad;
        memset(lw->half_line_hsync + g->active_off, U8_BLACK, len);
    }

    lw->data_tmpl = xmalloc(g->spl);
    memcpy(lw->data_tmpl, lw->blank_line, g->spl);
    memset(lw->data_tmpl, U8_SYNC, g->hsync);
}

static void lw_free(LineWriter *lw)
{
    free(lw->blank_line);
    free(lw->black_line);
    free(lw->short_half);
    free(lw->broad_half);
    free(lw->half_line_blank);
    free(lw->half_line_hsync);
    free(lw->data_tmpl);
}

/* Unpack 8 bytes into 72 bits (MSB first, matching numpy.unpackbits) */
static void unpack_bits(const uint8_t *bytes, uint8_t *bits /* 72 */)
{
    for (int i = 0; i < LINE_BYTES; i++) {
        uint8_t b = bytes[i];
        for (int k = 7; k >= 0; k--)
            bits[i * 8 + (7 - k)] = (b >> k) & 1;
    }
}

static void data_line(const LineWriter *lw, const uint8_t *line_bytes,
                      uint8_t *out /* g->spl */)
{
    const Geometry *g = lw->g;
    uint8_t bits[TOTAL_BITS_PER_LINE];
    uint32_t i;

    /* preamble + data bits */
    memcpy(bits, PREAMBLE, PREAMBLE_BITS);
    unpack_bits(line_bytes, bits + PREAMBLE_BITS);

    memcpy(out, lw->data_tmpl, g->spl);

    for (i = 0; i < g->active; i++) {
        uint32_t b = bit_map(g, i);
        out[g->active_off + i] = (bits[b] == 1) ? U8_BIT1 : U8_BIT0;
    }
}

/* ---------- write one field ---------- */

static void write_field(FILE *f, const LineWriter *lw,
                        const uint8_t *stream, int field_index)
{
    const Geometry *g = lw->g;
    uint8_t *line_buf = xmalloc(g->spl);   /* temporary for data lines */
    int i;

    if (field_index % 2 == 1) {            /* odd field */
        for (i = 0; i < 5; i++)
            fwrite(lw->broad_half, 1, g->half, f);
        for (i = 0; i < 5; i++)
            fwrite(lw->short_half, 1, g->half, f);
        for (i = 0; i < VBI_LINES; i++)
            fwrite(lw->black_line, 1, g->spl, f);
    } else {                               /* even field */
        for (i = 0; i < 5; i++)
            fwrite(lw->short_half, 1, g->half, f);
        for (i = 0; i < 5; i++)
            fwrite(lw->broad_half, 1, g->half, f);
        for (i = 0; i < 5; i++)
            fwrite(lw->short_half, 1, g->half, f);
        /* extra half-line */
        fwrite(lw->half_line_blank, 1, g->half, f);
        for (i = 0; i < VBI_LINES - 1; i++)
            fwrite(lw->black_line, 1, g->spl, f);
    }

    for (i = 0; i < DATA_LINES; i++) {
        data_line(lw, stream + i * LINE_BYTES, line_buf);
        fwrite(line_buf, 1, g->spl, f);
    }

    if (field_index % 2 == 0) {            /* even field end */
        fwrite(lw->half_line_hsync, 1, g->half, f);
        for (i = 0; i < 5; i++)
            fwrite(lw->short_half, 1, g->half, f);
    }

    free(line_buf);
}

/* ---------- main ---------- */

int main(int argc, char **argv)
{
    const char *input_path  = NULL;
    const char *output_path = NULL;
    uint32_t sample_rate    = DEFAULT_SAMPLE_RATE;
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
    int remaining = argc - i;
    if (remaining == 2) {
        input_path  = argv[argc - 2];
        output_path = argv[argc - 1];
    } else {
        fprintf(stderr, "usage: encode <input_file> <output.raw>\n", argv[0]);
        return 1;
    }

    Geometry g;
    geometry_init(&g, sample_rate);

    /* read whole input */
    FILE *inf = fopen(input_path, "rb");
    if (!inf) {
        perror(input_path);
        return 1;
    }
    if (fseek(inf, 0, SEEK_END) != 0) die("fseek failed");
    long fsz = ftell(inf);
    if (fsz < 0) die("ftell failed");
    rewind(inf);

    uint8_t *data = xmalloc((size_t)fsz);
    if (fread(data, 1, (size_t)fsz, inf) != (size_t)fsz)
        die("short read on input");
    fclose(inf);

    size_t data_len = (size_t)fsz;
    uint8_t *blob = NULL;
    size_t blob_len = 0;

    /* blob = MAGIC + length(4 LE) + crc32(4 LE) + data */
    if (mode == MODE_DATA) {
        blob_len = HEADER_LEN + data_len;
        blob = xmalloc(blob_len);
        memcpy(blob, MAGIC, 4);
        uint32_t len32 = (uint32_t)data_len;
        uint32_t crc   = crc32_zlib(data, data_len);
        blob[4] = (uint8_t)(len32);
        blob[5] = (uint8_t)(len32 >> 8);
        blob[6] = (uint8_t)(len32 >> 16);
        blob[7] = (uint8_t)(len32 >> 24);
        blob[8]  = (uint8_t)(crc);
        blob[9]  = (uint8_t)(crc >> 8);
        blob[10] = (uint8_t)(crc >> 16);
        blob[11] = (uint8_t)(crc >> 24);
        memcpy(blob + HEADER_LEN, data, data_len);
        free(data);
    } else {
        blob_len = data_len;
        blob = xmalloc(blob_len);
        memcpy(blob, data, data_len);
    }

    /* number of fields (round up, then force even) */
    size_t n_fields = (blob_len + FIELD_PAYLOAD - 1) / FIELD_PAYLOAD;
    if (n_fields == 0) n_fields = 1;
    n_fields += n_fields % 2;          /* whole frames */
    size_t n_frames = n_fields / 2;

    LineWriter lw;
    lw_init(&lw, &g);

    FILE *outf = fopen(output_path, "wb");
    if (!outf) {
        perror(output_path);
        return 1;
    }

    uint8_t stream[FIELD_STREAM_BYTES];

    for (size_t fi = 0; fi < n_fields; fi++) {
        size_t off = fi * FIELD_PAYLOAD;
        size_t plen = (off < blob_len) ? blob_len - off : 0;
        if (plen > FIELD_PAYLOAD) plen = FIELD_PAYLOAD;

        build_field_stream(blob + off, plen, stream);
        write_field(outf, &lw, stream, (int)(fi + 1));

        fprintf(stderr, "Encoded field %zu/%zu\n", fi + 1, n_fields);
    }

    fclose(outf);
    lw_free(&lw);
    free(blob);

    uint64_t total = (uint64_t)n_fields * g.field_samples;

    printf("\nInput           : %zu bytes (+%zu header)\n",
           data_len, blob_len - data_len);
    printf("Fields / frames : %zu / %zu\n", n_fields, n_frames);
    printf("Capacity used   : %zu/%zu payload bytes\n",
           blob_len, n_fields * FIELD_PAYLOAD);
    printf("Sample rate     : %u Hz, %u samples/line, %.2f samples/bit\n",
           g.sr, g.spl, g.nominal_spb);
    printf("Output          : %s (%llu samples, %.0f ms of video)\n",
           output_path, (unsigned long long)total,
           total / (double)g.sr * 1000.0);
    printf("Effective rate  : %.3f kB/s of payload\n",
           FIELD_PAYLOAD * 50.0 / 1000.0);

    return 0;
}
