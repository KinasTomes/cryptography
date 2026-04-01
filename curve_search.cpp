/**
 * Elliptic curve searching logic
 */

#include "ecc_generator.h"
#include <iostream>
#include <iomanip>

// ─── Curve search ────────────────────────────────────────────────────────────

/**
 * Search for random (a, b) pairs over F_p such that:
 *   1. The curve y² = x³ + ax + b is non-singular (discriminant ≠ 0 mod p).
 *   2. The group order N = #E(F_p) is a prime (cofactor h = 1).
 *   3. N ≠ p  (not anomalous — prevents SSSA attack).
 *
 * Logs progress to stdout.  Returns on first success.
 */
CurveResult search_curve(GEN p, long attempt_number,
                          std::chrono::steady_clock::time_point t0) {
    long ab_attempt = 0;

    while (true) {
        ab_attempt++;
        pari_sp av = avma;  // save stack for this (a,b) trial

        // Pick random a, b in [1, p-1]
        GEN a = random_field_element(p);
        GEN b = random_field_element(p);

        // ── Non-singularity: 4a³ + 27b² ≢ 0 (mod p) ──
        // Compute mod p to keep numbers small.
        GEN a3   = Fp_powu(a, 3, p);                    // a³ mod p
        GEN b2   = Fp_powu(b, 2, p);                    // b² mod p
        GEN disc = Fp_add(
                       Fp_mulu(a3, 4, p),
                       Fp_mulu(b2, 27, p),
                       p);

        if (equaliu(disc, 0)) {
            avma = av;
            continue;   // singular curve, skip
        }

        // ── Build PARI elliptic curve object ──
        // ellinit([a, b], p) — short Weierstrass form over F_p
        GEN curve_coeffs = mkvec2(a, b);
        GEN E = ellinit(curve_coeffs, p, 0);
        if (gequal0(E)) {
            avma = av;
            continue;   // degenerate
        }

        // ── Count points: N = #E(F_p) via Schoof-Elkies-Atkin ──
        GEN N = ellcard(E, p);

        // ── Anti-anomalous filter: N ≠ p ──
        if (equalii(N, p)) {
            avma = av;
            continue;
        }

        // ── Check N is prime ──
        if (!isprime(N)) {
            avma = av;
            continue;
        }

        // ── Success! Clone onto a longer-lived stack frame before returning ──
        CurveResult res;
        res.p = gcopy(p);
        res.a = gcopy(a);
        res.b = gcopy(b);
        res.N = gcopy(N);

        double t = elapsed_since(t0);
        std::cout
            << CYAN << "[Testing]" << RESET
            << " | Attempt #" << attempt_number
            << " | (a,b) try #" << ab_attempt
            << " | t=" << std::fixed << std::setprecision(1) << t << "s"
            << " | " << GREEN << "PRIME ORDER FOUND" << RESET << "\n";

        return res;
    }
}
