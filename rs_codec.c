#include "rs_codec.h"

#include <stdlib.h>
#include <string.h>

static uint8_t gf_exp[512];
static uint8_t gf_log[256];
static int gf_initialized;

static void gf_init(void) {
    if (gf_initialized) return;
    uint16_t x = 1;
    for (int i = 0; i < 255; ++i) {
        gf_exp[i] = (uint8_t)x;
        gf_log[x] = (uint8_t)i;
        x <<= 1;
        if (x & 0x100) x ^= 0x11d;
    }
    for (int i = 255; i < 512; ++i)
        gf_exp[i] = gf_exp[i - 255];
    gf_initialized = 1;
}

static uint8_t gf_mul(uint8_t x, uint8_t y) {
    if (!x || !y) return 0;
    return gf_exp[gf_log[x] + gf_log[y]];
}

static uint8_t gf_div(uint8_t x, uint8_t y) {
    if (!y) return 0; /* callers check this */
    if (!x) return 0;
    int d = (int)gf_log[x] - (int)gf_log[y];
    if (d < 0) d += 255;
    return gf_exp[d];
}

static uint8_t gf_pow(uint8_t x, int power) {
    if (!x) return power == 0 ? 1 : 0;
    int p = ((int)gf_log[x] * power) % 255;
    if (p < 0) p += 255;
    return gf_exp[p];
}

static uint8_t gf_inv(uint8_t x) {
    return x ? gf_exp[255 - gf_log[x]] : 0;
}

static uint8_t poly_eval(const uint8_t *p, size_t len, uint8_t x) {
    uint8_t y = p[0];
    for (size_t i = 1; i < len; ++i)
        y = gf_mul(y, x) ^ p[i];
    return y;
}

static void poly_add(const uint8_t *p, size_t plen,
                    const uint8_t *q, size_t qlen,
                    uint8_t *r, size_t rlen) {
    memset(r, 0, rlen);
    memcpy(r + rlen - plen, p, plen);
    for (size_t i = 0; i < qlen; ++i)
        r[i + rlen - qlen] ^= q[i];
}

static rs_status_t poly_mul(const uint8_t *p, size_t plen,
                            const uint8_t *q, size_t qlen,
                            uint8_t *r, size_t rlen) {
    if (rlen != plen + qlen - 1) return RS_ERR_INVALID_ARGUMENT;
    memset(r, 0, rlen);
    for (size_t j = 0; j < qlen; ++j) {
        if (!q[j]) continue;
        for (size_t i = 0; i < plen; ++i) {
            if (p[i]) r[i + j] ^= gf_mul(p[i], q[j]);
        }
    }
    return RS_OK;
}

/* Generator polynomial for nsym, built into caller-owned buffer. */
static rs_status_t generator_poly(size_t nsym, uint8_t *g) {
    size_t len = 1;
    g[0] = 1;
    for (size_t i = 0; i < nsym; ++i) {
        uint8_t next[257];
        uint8_t factor[2] = {1, gf_pow(2, (int)i)};
        rs_status_t st = poly_mul(g, len, factor, 2, next, len + 1);
        if (st != RS_OK) return st;
        memcpy(g, next, len + 1);
        ++len;
    }
    return RS_OK;
}

static void calc_syndromes(const uint8_t *msg, size_t mlen, size_t nsym,
                           uint8_t *synd) {
    synd[0] = 0;
    for (size_t i = 0; i < nsym; ++i)
        synd[i + 1] = poly_eval(msg, mlen, gf_pow(2, (int)i));
}

static rs_status_t find_errata_locator(const size_t *pos, size_t n,
                                       uint8_t *loc, size_t *loc_len,
                                       size_t msg_len) {
    (void)msg_len;
    loc[0] = 1;
    *loc_len = 1;
    for (size_t i = 0; i < n; ++i) {
        size_t coef_pos = msg_len - 1 - pos[i];
        uint8_t factor[3] = {gf_pow(2, (int)coef_pos), 0, 0};
        uint8_t one[1] = {1};
        uint8_t sum[2] = {0, 0};
        poly_add(one, 1, factor, 2, sum, 2);
        uint8_t tmp[257];
        if (poly_mul(loc, *loc_len, sum, 2, tmp, *loc_len + 1) != RS_OK)
            return RS_ERR_INVALID_ARGUMENT;
        memcpy(loc, tmp, *loc_len + 1);
        ++*loc_len;
    }
    return RS_OK;
}

static rs_status_t find_error_evaluator(const uint8_t *synd, size_t synd_len,
                                        const uint8_t *err_loc, size_t loc_len,
                                        size_t nsym, uint8_t *eval,
                                        size_t *eval_len) {
    size_t prod_len = synd_len + loc_len - 1;
    uint8_t prod[512];
    if (poly_mul(synd, synd_len, err_loc, loc_len, prod, prod_len) != RS_OK)
        return RS_ERR_INVALID_ARGUMENT;
    size_t div_len = nsym + 2;
    uint8_t divisor[257];
    memset(divisor, 0, div_len);
    divisor[0] = 1;

    uint8_t out[512];
    memcpy(out, prod, prod_len);
    for (size_t i = 0; i + div_len <= prod_len; ++i) {
        uint8_t coef = out[i];
        if (!coef) continue;
        for (size_t j = 1; j < div_len; ++j)
            if (divisor[j]) out[i + j] ^= gf_mul(divisor[j], coef);
    }
    size_t rem_len = div_len - 1;
    if (rem_len > prod_len) rem_len = prod_len;
    memcpy(eval, out + prod_len - rem_len, rem_len);
    *eval_len = rem_len;
    return RS_OK;
}

static rs_status_t correct_errata(uint8_t *msg, size_t mlen,
                                  const uint8_t *synd, size_t synd_len,
                                  const size_t *err_pos, size_t nerr) {
    if (!nerr) return RS_OK;

    size_t coef_pos[256];
    for (size_t i = 0; i < nerr; ++i)
        coef_pos[i] = mlen - 1 - err_pos[i];

    uint8_t err_loc[257], err_eval[257];
    size_t err_loc_len, err_eval_len;
    rs_status_t st = find_errata_locator(err_pos, nerr, err_loc,
                                         &err_loc_len, mlen);
    if (st != RS_OK) return st;

    /* Python uses synd[::-1] here. */
    uint8_t rev_synd[256];
    for (size_t i = 0; i < synd_len; ++i)
        rev_synd[i] = synd[synd_len - 1 - i];

    st = find_error_evaluator(rev_synd, synd_len, err_loc, err_loc_len,
                              err_loc_len - 1, err_eval, &err_eval_len);
    if (st != RS_OK) return st;

    /* Python reverses the evaluator before evaluating. */
    uint8_t X[256];
    for (size_t i = 0; i < nerr; ++i)
        X[i] = gf_pow(2, -(255 - (int)coef_pos[i]));

    for (size_t i = 0; i < nerr; ++i) {
        uint8_t Xi = X[i];
        uint8_t Xi_inv = gf_inv(Xi);
        uint8_t prime = 1;
        for (size_t j = 0; j < nerr; ++j) {
            if (j != i)
                prime = gf_mul(prime, (uint8_t)(1 ^ gf_mul(Xi_inv, X[j])));
        }
        if (!prime) return RS_ERR_CORRECT;

        uint8_t y = poly_eval(err_eval, err_eval_len, Xi_inv);
        y = gf_mul(Xi, y);
        if (!prime) return RS_ERR_CORRECT;
        msg[err_pos[i]] ^= gf_div(y, prime);
    }
    return RS_OK;
}

static rs_status_t find_error_locator(const uint8_t *synd, size_t synd_len,
                                      size_t nsym, const uint8_t *erase_loc,
                                      size_t erase_len, size_t erase_count,
                                      uint8_t *out, size_t *out_len) {
    uint8_t err_loc[257], old_loc[257];
    size_t err_len, old_len;
    if (erase_loc && erase_len) {
        memcpy(err_loc, erase_loc, erase_len);
        memcpy(old_loc, erase_loc, erase_len);
        err_len = old_len = erase_len;
    } else {
        err_loc[0] = old_loc[0] = 1;
        err_len = old_len = 1;
    }

    size_t synd_shift = synd_len - nsym;
    for (size_t i = 0; i < nsym - erase_count; ++i) {
        size_t K = erase_loc && erase_len ? erase_count + i + synd_shift
                                          : i + synd_shift;
        uint8_t delta = synd[K];
        for (size_t j = 1; j < err_len; ++j)
            delta ^= gf_mul(err_loc[err_len - (j + 1)], synd[K - j]);

        old_loc[old_len++] = 0;
        if (delta) {
            if (old_len > err_len) {
                size_t prev_err_len = err_len;
                size_t prev_old_len = old_len;
                uint8_t new_loc[257];
                for (size_t j = 0; j < prev_old_len; ++j)
                    new_loc[j] = gf_mul(old_loc[j], delta);

                uint8_t inv = gf_inv(delta);
                for (size_t j = 0; j < prev_err_len; ++j)
                    old_loc[j] = gf_mul(err_loc[j], inv);
                old_len = prev_err_len;

                memcpy(err_loc, new_loc, prev_old_len);
                err_len = prev_old_len;
            }
            size_t rlen = err_len > old_len ? err_len : old_len;
            uint8_t scaled[257], sum[257];
            for (size_t j = 0; j < old_len; ++j)
                scaled[j] = gf_mul(old_loc[j], delta);
            poly_add(err_loc, err_len, scaled, old_len, sum, rlen);
            memcpy(err_loc, sum, rlen);
            err_len = rlen;
        }
    }
    size_t first = 0;
    while (first < err_len && err_loc[first] == 0) ++first;
    err_len -= first;
    if (first) memmove(err_loc, err_loc + first, err_len);
    size_t errs = err_len ? err_len - 1 : 0;
    long capacity_used = 2L * ((long)errs - (long)erase_count) + (long)erase_count;
    if (capacity_used > (long)nsym)
        return RS_ERR_TOO_MANY_ERRORS;
    memcpy(out, err_loc, err_len);
    *out_len = err_len;
    return RS_OK;
}

static rs_status_t find_errors(const uint8_t *err_loc, size_t err_len,
                               size_t nmess, size_t *err_pos, size_t *nfound) {
    size_t errs = err_len - 1;
    size_t n = 0;
    for (size_t i = 0; i < nmess; ++i) {
        if (poly_eval(err_loc, err_len, gf_pow(2, (int)i)) == 0) {
            if (n >= errs) return RS_ERR_LOCATE_ERRORS;
            err_pos[n++] = nmess - 1 - i;
        }
    }
    if (n != errs) return RS_ERR_LOCATE_ERRORS;
    *nfound = n;
    return RS_OK;
}

static void forney_syndromes(const uint8_t *synd, size_t nsym,
                             const size_t *pos, size_t npos, size_t nmess,
                             uint8_t *fsynd) {
    size_t len = nsym;
    memcpy(fsynd, synd + 1, nsym);
    for (size_t i = 0; i < npos; ++i) {
        size_t reversed = nmess - 1 - pos[i];
        uint8_t x = gf_pow(2, (int)reversed);
        for (size_t j = 0; j + 1 < len; ++j)
            fsynd[j] = gf_mul(fsynd[j], x) ^ fsynd[j + 1];
    }
}

rs_status_t rs_encode(const uint8_t *data, size_t data_len,
                      size_t nsym, uint8_t *out_codeword, size_t out_len) {
    gf_init();
    if ((!data && data_len) || !out_codeword || nsym == 0 || nsym >= 256)
        return RS_ERR_INVALID_ARGUMENT;
    if (data_len + nsym > 255) return RS_ERR_MESSAGE_TOO_LONG;
    if (out_len < data_len + nsym) return RS_ERR_INVALID_ARGUMENT;

    uint8_t gen[257];
    rs_status_t st = generator_poly(nsym, gen);
    if (st != RS_OK) return st;

    size_t total = data_len + nsym;
    uint8_t *work = (uint8_t *)calloc(total, 1);
    if (!work) return RS_ERR_ALLOC;
    memcpy(work, data, data_len);
    for (size_t i = 0; i < data_len; ++i) {
        uint8_t coef = work[i];
        if (!coef) continue;
        for (size_t j = 1; j < nsym + 1; ++j)
            work[i + j] ^= gf_mul(gen[j], coef);
    }
    memcpy(out_codeword, data, data_len);
    memcpy(out_codeword + data_len, work + data_len, nsym);
    free(work);
    return RS_OK;
}

rs_status_t rs_decode(uint8_t *codeword, size_t codeword_len, size_t nsym,
                      const size_t *erase_pos, size_t erase_count,
                      uint8_t *out_data, size_t out_data_len,
                      size_t *corrected_count) {
    gf_init();
    if (!codeword || !out_data || !corrected_count || nsym == 0 || nsym >= 256)
        return RS_ERR_INVALID_ARGUMENT;
    if (codeword_len > 255 || codeword_len < nsym)
        return RS_ERR_MESSAGE_TOO_LONG;
    if (out_data_len < codeword_len - nsym)
        return RS_ERR_INVALID_ARGUMENT;
    if (erase_count > nsym) return RS_ERR_TOO_MANY_ERASURES;
    if (erase_count && !erase_pos) return RS_ERR_INVALID_ARGUMENT;

    for (size_t i = 0; i < erase_count; ++i) {
        if (erase_pos[i] >= codeword_len) return RS_ERR_INVALID_ARGUMENT;
        codeword[erase_pos[i]] = 0;
    }

    uint8_t synd[256];
    calc_syndromes(codeword, codeword_len, nsym, synd);
    int clean = 1;
    for (size_t i = 1; i <= nsym; ++i)
        if (synd[i]) { clean = 0; break; }
    if (clean) {
        memcpy(out_data, codeword, codeword_len - nsym);
        *corrected_count = 0;
        return RS_OK;
    }

    uint8_t fsynd[256];
    forney_syndromes(synd, nsym, erase_pos, erase_count, codeword_len, fsynd);

    uint8_t err_loc[257];
    size_t err_loc_len;
    rs_status_t st = find_error_locator(fsynd, nsym, nsym,
                                        NULL, 0, erase_count,
                                        err_loc, &err_loc_len);
    if (st != RS_OK) return st;

    /* Python calls rs_find_errors(err_loc[::-1], ...). */
    uint8_t rev_loc[257];
    for (size_t i = 0; i < err_loc_len; ++i)
        rev_loc[i] = err_loc[err_loc_len - 1 - i];

    size_t err_pos[256], err_count;
    st = find_errors(rev_loc, err_loc_len, codeword_len,
                     err_pos, &err_count);
    if (st != RS_OK) return st;

    size_t all_pos[256];
    for (size_t i = 0; i < erase_count; ++i) all_pos[i] = erase_pos[i];
    for (size_t i = 0; i < err_count; ++i) all_pos[erase_count + i] = err_pos[i];

    st = correct_errata(codeword, codeword_len, synd, nsym,
                        all_pos, erase_count + err_count);
    if (st != RS_OK) return st;

    calc_syndromes(codeword, codeword_len, nsym, synd);
    for (size_t i = 1; i <= nsym; ++i)
        if (synd[i]) return RS_ERR_CORRECT;

    memcpy(out_data, codeword, codeword_len - nsym);
    *corrected_count = erase_count + err_count;
    return RS_OK;
}

const char *rs_status_string(rs_status_t status) {
    switch (status) {
    case RS_OK: return "ok";
    case RS_ERR_INVALID_ARGUMENT: return "invalid argument";
    case RS_ERR_MESSAGE_TOO_LONG: return "message too long";
    case RS_ERR_TOO_MANY_ERASURES: return "too many erasures";
    case RS_ERR_TOO_MANY_ERRORS: return "too many errors";
    case RS_ERR_LOCATE_ERRORS: return "could not locate errors";
    case RS_ERR_CORRECT: return "could not correct message";
    case RS_ERR_DIV_ZERO: return "division by zero";
    case RS_ERR_ALLOC: return "allocation failed";
    default: return "unknown error";
    }
}
