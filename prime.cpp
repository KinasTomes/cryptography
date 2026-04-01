/**
 * Prime number generation
 */

#include "ecc_generator.h"

// ─── Prime generation ────────────────────────────────────────────────────────

/**
 * Generate a random prime of exactly `bits` bits using libsodium for entropy
 * and PARI's isprime() for certification.
 *
 * Returns a GEN allocated on the PARI stack.
 */
GEN generate_random_prime(int bits) {
    while (true) {
        pari_sp av = avma;

        GEN candidate = random_odd_of_bits(bits);

        // PARI isprime: 0 = composite, 1 = probable prime, 2 = proven prime
        // For speed we use ispseudoprime (BPSW) then confirm with isprime.
        if (!ispseudoprime(candidate, 0)) {
            avma = av;   // discard
            continue;
        }
        if (!isprime(candidate)) {
            avma = av;
            continue;
        }
        return candidate;   // caller owns this stack object
    }
}
