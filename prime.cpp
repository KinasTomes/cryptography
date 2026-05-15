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

        // Filter composites quickly before running the full primality proof.
        if (!ispseudoprime(candidate, 0)) {
            avma = av;
            continue;
        }
        if (!isprime(candidate)) {
            avma = av;
            continue;
        }
        return candidate;
    }
}
