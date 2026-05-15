/**
 * Helper functions for UI and utilities
 */

#include "ecc_generator.h"
#include <iostream>
#include <chrono>

// ─── Color constants ────────────────────────────────────────────────────────
const std::string RESET  = "\033[0m";
const std::string BOLD   = "\033[1m";
const std::string GREEN  = "\033[32m";
const std::string YELLOW = "\033[33m";
const std::string CYAN   = "\033[36m";
const std::string RED    = "\033[31m";

// Return elapsed seconds since a reference point.
double elapsed_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Print a GEN value as a decimal string (caller must own a PARI stack frame).
std::string gen_to_str(GEN g) {
    char *s = GENtostr(g);
    std::string result(s);
    pari_free(s);
    return result;
}

// Separator lines for the UI.
void sep(char c, int w) {
    std::cout << std::string(w, c) << "\n";
}
