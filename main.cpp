/**
 * High-Performance Elliptic Curve Parameter Generator
 * Uses libpari for number theory and libsodium for CSPRNG.
 *
 * Build (Linux/MSYS2/WSL):
 *   g++ -O2 -o ecc_pari_generator main.cpp helpers.cpp random.cpp prime.cpp curve_search.cpp \
 *       $(pkg-config --cflags --libs pari) -lsodium -lstdc++
 *
 * Windows (MSYS2 MinGW64):
 *   g++ -O2 -o ecc_pari_generator main.cpp helpers.cpp random.cpp prime.cpp curve_search.cpp \
 *       -I/mingw64/include -L/mingw64/lib -lpari -lsodium -lstdc++
 */

#include "ecc_generator.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sodium.h>

// ─── Main ────────────────────────────────────────────────────────────────────

int main() {
    // ── Configuration ──────────────────────────────────────────────────────
    const int TARGET_BITS = 128;   // <-- change this: e.g. 64, 128, 192, 256
    // ───────────────────────────────────────────────────────────────────────

    // Initialise libpari: 500 MB stack, enough for 256-bit+ SEA point counts.
    pari_init(500000000UL, 2);

    // Initialise libsodium
    if (sodium_init() < 0) {
        std::cerr << RED << "[FATAL] libsodium initialisation failed.\n" << RESET;
        return 1;
    }

    sep('=');
    std::cout
        << BOLD << "  ECC Parameter Generator  (libpari + libsodium)\n" << RESET
        << "  Target: " << TARGET_BITS << "-bit prime p, cofactor h=1, N prime, N≠p\n"
        << "  Curve:  y² ≡ x³ + ax + b  (mod p)\n";
    sep('=');
    std::cout << "\n";

    auto t0 = std::chrono::steady_clock::now();

    long attempt = 0;
    CurveResult found;

    while (true) {
        attempt++;
        pari_sp av = avma;   // stack frame for this prime candidate

        // ── 1. Generate prime p ──
        std::cout
            << YELLOW << "[Searching]" << RESET
            << " | Attempt #" << attempt
            << " | Generating " << TARGET_BITS << "-bit prime p...\r"
            << std::flush;

        GEN p = generate_random_prime(TARGET_BITS);

        std::string p_str = gen_to_str(p);
        std::cout
            << YELLOW << "[Searching]" << RESET
            << " | Attempt #" << attempt
            << " | p = " << p_str.substr(0, 20)
            << (p_str.size() > 20 ? "..." : "")
            << " | Searching (a,b)...    \n";

        // ── 2. Search for a curve with prime order over this p ──
        found = search_curve(p, attempt, t0);

        // search_curve returns only on success; we're done.
        (void)av;   // stack state is managed inside search_curve
        break;
    }

    double elapsed = elapsed_since(t0);

    // ── Print results ──
    std::cout << "\n";
    sep('=');
    std::cout << BOLD << GREEN << "  *** SECURE CURVE FOUND ***\n" << RESET;
    sep('=');
    std::cout << "  Bits  : " << TARGET_BITS << "\n";
    std::cout << "  p     = " << gen_to_str(found.p) << "\n";
    std::cout << "  a     = " << gen_to_str(found.a) << "\n";
    std::cout << "  b     = " << gen_to_str(found.b) << "\n";
    std::cout << "  N     = " << gen_to_str(found.N) << "\n";
    std::cout << "  h     = 1  (N is prime, cofactor = 1)\n";
    std::cout << "  N ≠ p : " << (equalii(found.N, found.p) ? "NO (anomalous!)" : "YES (safe)") << "\n";
    sep('-');
    std::cout << "  isprime(N) : " << (isprime(found.N) ? GREEN + "TRUE" + RESET : RED + "FALSE" + RESET) << "\n";
    std::cout << "  isprime(p) : " << (isprime(found.p) ? GREEN + "TRUE" + RESET : RED + "FALSE" + RESET) << "\n";
    sep('-');
    std::cout
        << "  Time elapsed: "
        << BOLD << std::fixed << std::setprecision(3) << elapsed << "s" << RESET
        << "\n";
    sep('=');

    pari_close();
    return 0;
}
