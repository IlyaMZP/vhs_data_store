#include "rs_codec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die(const char *what, rs_status_t st) {
    fprintf(stderr, "%s: %s (%d)\n", what, rs_status_string(st), st);
    exit(1);
}

static size_t rand_below(size_t n) {
    return (size_t)(rand() % (int)n);
}

static void shuffle(size_t *a, size_t n) {
    for (size_t i = n; i > 1; --i) {
        size_t j = rand_below(i);
        size_t t = a[i - 1]; a[i - 1] = a[j]; a[j] = t;
    }
}

int main(void) {
    enum { N = 255, K = 223, NSYM = N - K };
    srand(7);

    uint8_t msg[K], cw[N], damaged[N], decoded[K];
    size_t positions[N];

    for (int trial = 0; trial < 200; ++trial) {
        for (size_t i = 0; i < K; ++i) msg[i] = (uint8_t)(rand() & 255);
        rs_status_t st = rs_encode(msg, K, NSYM, cw, sizeof cw);
        if (st != RS_OK) die("encode", st);
        memcpy(damaged, cw, sizeof damaged);
        for (size_t i = 0; i < N; ++i) positions[i] = i;
        shuffle(positions, N);

        size_t n_erase = rand_below(20);
        size_t n_err = rand_below((NSYM - n_erase) / 2 + 2);
        if (n_err > (NSYM - n_erase) / 2) n_err = (NSYM - n_erase) / 2;
        size_t total = n_erase + n_err;
        size_t erase_pos[NSYM];

        for (size_t i = 0; i < total; ++i) {
            size_t p = positions[i];
            uint8_t e;
            do e = (uint8_t)(rand() & 255); while (e == 0);
            damaged[p] ^= e;
            if (i < n_erase) erase_pos[i] = p;
        }

        size_t corrected = 0;
        st = rs_decode(damaged, N, NSYM, erase_pos, n_erase,
                       decoded, sizeof decoded, &corrected);
        if (st != RS_OK) die("decode", st);
        if (memcmp(decoded, msg, K) != 0) {
            fprintf(stderr, "trial %d: decoded data mismatch\n", trial);
            return 1;
        }
        if (corrected != total) {
            fprintf(stderr, "trial %d: expected %zu corrections, got %zu\n",
                    trial, total, corrected);
            return 1;
        }
    }

    /* Clean codeword path. */
    {
        rs_status_t st = rs_encode(msg, K, NSYM, cw, sizeof cw);
        if (st != RS_OK) die("encode clean", st);
        size_t corrected = 99;
        st = rs_decode(cw, N, NSYM, NULL, 0, decoded, sizeof decoded, &corrected);
        if (st != RS_OK || corrected != 0 || memcmp(decoded, msg, K) != 0)
            die("clean decode", st);
    }

    printf("rs_codec self-test: 200/200 random error+erasure trials passed\n");
    return 0;
}
