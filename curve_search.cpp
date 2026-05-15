/**
 * Elliptic curve searching logic
 */

#include "ecc_generator.h"
#include <iostream>
#include <iomanip>

// ─── Curve search ────────────────────────────────────────────────────────────

/**
 * Search curves over F_p in the fixed-a family y^2 = x^3 - 3x + b such that:
 *   1. The curve is non-singular.
 *   2. The group order N = #E(F_p) is prime (cofactor h = 1).
 *   3. N != p (not anomalous).
 *
 * Logs progress to stdout. Returns on first success.
 */
static std::string trunc(const std::string &s, int w) {
    if ((int)s.size() <= w) return std::string(w - s.size(), ' ') + s;
    return s.substr(0, w - 3) + "...";
}

CurveResult search_curve(GEN p, long attempt_number,
                          std::chrono::steady_clock::time_point t0,
                          bool verbose) {
    long ab_attempt = 0;
    GEN a = subii(p, stoi(3));   // a = -3 mod p

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

        GEN b = random_field_element(p);

        // For a = -3, discriminant is singular iff b^2 == 4 (mod p).
        GEN b2   = Fp_powu(b, 2, p);
        if (equaliu(b2, 4)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << "-3"
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
                    << std::setw(20) << "-3"
                    << std::setw(20) << trunc(gen_to_str(b), 18)
                    << std::setw(20) << "-"
                    << std::setw(8)  << (std::to_string((int)t) + "s")
                    << RED << std::setw(20) << "DEGENERATE" << RESET
                    << "\n";
            }
            avma = av;
            continue;
        }

        // Early-abort SEA: returns 0 as soon as a small bad factor is detected.
        GEN N = ellsea(E, 1);
        std::string N_str;
        if (verbose && !gequal0(N)) N_str = gen_to_str(N);

        if (gequal0(N)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << "-3"
                    << std::setw(20) << trunc(gen_to_str(b), 18)
                    << std::setw(20) << "-"
                    << std::setw(8)  << (std::to_string((int)t) + "s")
                    << RED << std::setw(20) << "SEA EARLY ABORT" << RESET
                    << "\n";
            }
            avma = av;
            continue;
        }

        if (equalii(N, p)) {
            if (verbose) {
                double t = elapsed_since(t0);
                std::cout
                    << std::left
                    << std::setw(7)  << attempt_number
                    << std::setw(7)  << ab_attempt
                    << std::setw(20) << "-3"
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
                    << std::setw(20) << "-3"
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
                << std::setw(20) << "-3"
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
