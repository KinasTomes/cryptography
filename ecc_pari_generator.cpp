/**
 * High-Performance Elliptic Curve Parameter Generator
 * Uses libpari for number theory and libsodium for CSPRNG.
 *
 * Build (Linux/MSYS2/WSL):
 *   g++ -O2 -o ecc_pari_generator ecc_pari_generator.cpp \
 *       $(pkg-config --cflags --libs pari) -lsodium -lstdc++
 *
 * Windows (MSYS2 MinGW64):
 *   g++ -O2 -o ecc_pari_generator ecc_pari_generator.cpp \
 *       -I/mingw64/include -L/mingw64/lib -lpari -lsodium -lstdc++
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>
#include <cstring>

// libpari
#include <pari/pari.h>

// libsodium
#include <sodium.h>

// ─── Helpers ────────────────────────────────────────────────────────────────

static const std::string RESET  = "\033[0m";
static const std::string BOLD   = "\033[1m";
static const std::string GREEN  = "\033[32m";
static const std::string YELLOW = "\033[33m";
static const std::string CYAN   = "\033[36m";
static const std::string RED    = "\033[31m";

// Return elapsed seconds since a reference point.
static double elapsed_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Print a GEN value as a decimal string (caller must own a PARI stack frame).
static std::string gen_to_str(GEN g) {
    // GENtostr allocates on the PARI stack; wrap in std::string immediately.
    char *s = GENtostr(g);
    std::string result(s);
    return result;
}

// Separator lines for the UI.
static void sep(char c = '-', int w = 72) {
    std::cout << std::string(w, c) << "\n";
}

// ─── Random bytes via libsodium ─────────────────────────────────────────────

/**
 * Fill `len` bytes with cryptographically secure random data.
 * Thin wrapper so the rest of the code never calls randombytes_buf directly.
 */
static void csprng_bytes(void *buf, size_t len) {
    randombytes_buf(buf, len);
}

/**
 * Build a PARI integer of exactly `bits` bits whose top bit and bottom bit
 * are both set (odd, correct magnitude), using libsodium for the raw bytes.
 *
 * The integer is allocated on the current PARI stack frame.
 */
static GEN random_odd_of_bits(int bits) {
    int bytes = (bits + 7) / 8;

    // Allocate a temporary byte buffer.
    unsigned char *buf = (unsigned char *)pari_malloc(bytes);
    csprng_bytes(buf, bytes);

    // Force top bit → correct magnitude.
    int top_bit_pos = (bits - 1) % 8;
    buf[0] |=  (unsigned char)(1u << top_bit_pos);   // set bit (bits-1)
    buf[0] &= (unsigned char)((1u << (top_bit_pos + 1)) - 1u); // clear higher bits

    // Force bottom bit → odd.
    buf[bytes - 1] |= 0x01u;

    // Convert raw bytes (big-endian) to a PARI integer.
    GEN result = strtoi("0");  // start at 0
    for (int i = 0; i < bytes; i++) {
        result = addii(mulis(result, 256), stoi(buf[i]));
    }

    pari_free(buf);
    return result;
}

/**
 * Generate a random field element in [1, p-1] using libsodium.
 * Rejection-samples until the value falls in range.
 */
static GEN random_field_element(GEN p) {
    long bits = expi(p) + 1;   // bit-length of p
    GEN result;
    do {
        pari_sp av = avma;      // save stack
        GEN candidate = random_odd_of_bits(bits);
        // candidate may be >= p; reduce mod p
        result = gmod(candidate, p);
        if (equaliu(result, 0)) {
            avma = av;          // discard and retry
            continue;
        }
        // Keep result (do NOT restore avma — we want to return it)
        break;
    } while (true);
    return result;
}

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

// ─── Curve search ────────────────────────────────────────────────────────────

struct CurveResult {
    GEN p, a, b, N;   // all on PARI stack — clone before freeing outer frame
};

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
