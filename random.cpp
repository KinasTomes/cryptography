/**
 * Random number generation using libsodium
 */

#include "ecc_generator.h"
#include <sodium.h>

// ─── Random bytes via libsodium ─────────────────────────────────────────────

/**
 * Fill `len` bytes with cryptographically secure random data.
 * Thin wrapper so the rest of the code never calls randombytes_buf directly.
 */
void csprng_bytes(void *buf, size_t len) {
    randombytes_buf(buf, len);
}

/**
 * Build a PARI integer of exactly `bits` bits whose top bit and bottom bit
 * are both set (odd, correct magnitude), using libsodium for the raw bytes.
 *
 * The integer is allocated on the current PARI stack frame.
 */
GEN random_odd_of_bits(int bits) {
    int bytes = (bits + 7) / 8;

    unsigned char *buf = (unsigned char *)pari_malloc(bytes);
    csprng_bytes(buf, bytes);

    // Force top bit → correct magnitude.
    int top_bit_pos = (bits - 1) % 8;
    buf[0] |=  (unsigned char)(1u << top_bit_pos);
    buf[0] &= (unsigned char)((1u << (top_bit_pos + 1)) - 1u);

    // Force bottom bit → odd.
    buf[bytes - 1] |= 0x01u;

    // Convert bytes -> hex string -> GEN in one PARI call.
    char *hex = (char *)pari_malloc(bytes * 2 + 3);
    hex[0] = '0';
    hex[1] = 'x';
    for (int i = 0; i < bytes; i++)
        sprintf(hex + 2 + i * 2, "%02x", buf[i]);
    hex[bytes * 2 + 2] = '\0';

    GEN result = strtoi(hex);

    pari_free(hex);
    pari_free(buf);
    return result;
}

/**
 * Generate a random field element in [1, p-1] using libsodium.
 * Rejection-samples until the value falls in range.
 */
GEN random_field_element(GEN p) {
    long bits = expi(p) + 1;   // bit-length of p
    while (true) {
        pari_sp av = avma;
        GEN candidate = random_odd_of_bits(bits);
        GEN result = gmod(candidate, p);
        if (equaliu(result, 0)) {
            avma = av;
            continue;
        }
        return result;
    }
}
