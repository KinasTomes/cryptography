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
static std::string trunc(const std::string &s, int w) {
    if ((int)s.size() <= w) return std::string(w - s.size(), ' ') + s;
    return s.substr(0, w - 3) + "...";
}

CurveResult search_curve(GEN p, long attempt_number,
                          std::chrono::steady_clock::time_point t0,
                          bool verbose) {
    long ab_attempt = 0;

    if (verbose) {
        std::cout
            << CYAN << std::left
            << std::setw(7)  << "pAtt"
            << std::setw(7)  << "abTry"
            << std::setw(20) << "a"
            << std::setw(20) << "b"
            << std::setw(20) << "N"
            << std::setw(8)  << "t"
            << std::setw(20) << "Result"
            << RESET << "\n";
    }

    while (true) {
        ab_attempt++;
        pari_sp av = avma;

        GEN a = random_field_element(p);
        GEN b = random_field_element(p);

        // ── Non-singularity: 4a³ + 27b² ≢ 0 (mod p) ──
        GEN a3   = Fp_powu(a, 3, p);
        GEN b2   = Fp_powu(b, 2, p);
        GEN disc = Fp_add(
                       Fp_mulu(a3, 4, p),
                       Fp_mulu(b2, 27, p),
                       p);

        if (equaliu(disc, 0)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << trunc(gen_to_str(a), 18)
                    << std::setw(20) << trunc(gen_to_str(b), 18)
                    << std::setw(20) << "-"
                    << std::setw(8)  << (std::to_string((int)t) + "s")
                    << RED << std::setw(20) << "SINGULAR" << RESET
                    << "\n";
            }
            avma = av;
            continue;
        }

        GEN curve_coeffs = mkvec2(a, b);
        GEN E = ellinit(curve_coeffs, p, 0);
        if (gequal0(E)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << trunc(gen_to_str(a), 18)
                    << std::setw(20) << trunc(gen_to_str(b), 18)
                    << std::setw(20) << "-"
                    << std::setw(8)  << (std::to_string((int)t) + "s")
                    << RED << std::setw(20) << "DEGENERATE" << RESET
                    << "\n";
            }
            avma = av;
            continue;
        }

        GEN N = ellcard(E, p);
        std::string N_str;
        if (verbose) N_str = gen_to_str(N);

        if (equalii(N, p)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << trunc(gen_to_str(a), 18)
                    << std::setw(20) << trunc(gen_to_str(b), 18)
                    << std::setw(20) << trunc(N_str, 18)
                    << std::setw(8)  << (std::to_string((int)t) + "s")
                    << YELLOW << std::setw(20) << "ANOMALOUS (N=p)" << RESET
                    << "\n";
            }
            avma = av;
            continue;
        }

        if (!isprime(N)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << trunc(gen_to_str(a), 18)
                    << std::setw(20) << trunc(gen_to_str(b), 18)
                    << std::setw(20) << trunc(N_str, 18)
                    << std::setw(8)  << (std::to_string((int)t) + "s")
                    << RED << std::setw(20) << "N NOT PRIME" << RESET
                    << "\n";
            }
            avma = av;
            continue;
        }

        // ── Success! ──
        CurveResult res;
        res.p = gcopy(p);
        res.a = gcopy(a);
        res.b = gcopy(b);
        res.N = gcopy(N);

        double t = elapsed_since(t0);
        if (verbose) {
            if (N_str.empty()) N_str = gen_to_str(N);
            std::cout
                << std::left
                << std::setw(7)  << attempt_number
                << std::setw(7)  << ab_attempt
                << std::setw(20) << trunc(gen_to_str(a), 18)
                << std::setw(20) << trunc(gen_to_str(b), 18)
                << std::setw(20) << trunc(N_str, 18)
                << std::setw(8)  << (std::to_string((int)t) + "s")
                << GREEN << std::setw(20) << "PRIME ORDER FOUND!" << RESET
                << "\n";
        } else {
            std::cout
                << CYAN << "[Testing]" << RESET
                << " | Attempt #" << attempt_number
                << " | (a,b) try #" << ab_attempt
                << " | t=" << std::fixed << std::setprecision(1) << t << "s"
                << " | " << GREEN << "PRIME ORDER FOUND" << RESET << "\n";
        }

        return res;
    }
}
