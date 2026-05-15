/**
 * High-Performance Elliptic Curve Parameter Generator
 * Header file with shared declarations
 */

#ifndef ECC_GENERATOR_H
#define ECC_GENERATOR_H

#include <string>
#include <chrono>
#include <pari/pari.h>

// ─── Color constants ────────────────────────────────────────────────────────
extern const std::string RESET;
extern const std::string BOLD;
extern const std::string GREEN;
extern const std::string YELLOW;
extern const std::string CYAN;
extern const std::string RED;

// ─── Helper functions ───────────────────────────────────────────────────────
double elapsed_since(std::chrono::steady_clock::time_point t0);
std::string gen_to_str(GEN g);
void sep(char c = '-', int w = 72);

// ─── Random number generation ───────────────────────────────────────────────
void csprng_bytes(void *buf, size_t len);
GEN random_odd_of_bits(int bits);
GEN random_field_element(GEN p);

// ─── Prime generation ───────────────────────────────────────────────────────
GEN generate_random_prime(int bits);

// ─── Curve search ───────────────────────────────────────────────────────────
struct CurveResult {
    GEN p, a, b, N;   // all on PARI stack — clone before freeing outer frame
};

CurveResult search_curve(GEN p, long attempt_number,
                         std::chrono::steady_clock::time_point t0,
                         bool verbose = false);

#endif // ECC_GENERATOR_H
