#ifndef RS_CODEC_H
#define RS_CODEC_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RS_CODEC_MAX_SYMBOLS 255

typedef enum {
    RS_OK = 0,
    RS_ERR_INVALID_ARGUMENT = -1,
    RS_ERR_MESSAGE_TOO_LONG = -2,
    RS_ERR_TOO_MANY_ERASURES = -3,
    RS_ERR_TOO_MANY_ERRORS = -4,
    RS_ERR_LOCATE_ERRORS = -5,
    RS_ERR_CORRECT = -6,
    RS_ERR_DIV_ZERO = -7,
    RS_ERR_ALLOC = -8
} rs_status_t;

/*
 * Systematic Reed-Solomon codec over GF(2^8), primitive polynomial 0x11d,
 * first consecutive root (fcr) = 0.
 *
 * nsym is the number of parity symbols. A codeword has data_len + nsym
 * symbols and must contain at most 255 symbols total.
 */

/* Encode data into out_codeword. out_len must be >= data_len + nsym. */
rs_status_t rs_encode(const uint8_t *data, size_t data_len,
                      size_t nsym, uint8_t *out_codeword, size_t out_len);

/*
 * Decode/correct a codeword in place.
 * erase_pos contains known-bad symbol positions (0 = first symbol).
 * On success, the corrected data portion is copied to out_data.
 * corrected_count receives the number of corrected symbols.
 */
rs_status_t rs_decode(uint8_t *codeword, size_t codeword_len, size_t nsym,
                      const size_t *erase_pos, size_t erase_count,
                      uint8_t *out_data, size_t out_data_len,
                      size_t *corrected_count);

/* Human-readable status string. */
const char *rs_status_string(rs_status_t status);

#ifdef __cplusplus
}
#endif

#endif
