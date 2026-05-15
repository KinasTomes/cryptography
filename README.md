# ECC Parameter Generator

Generates cryptographically secure elliptic curve parameters using **libpari** (number theory / SEA point counting) and **libsodium** (CSPRNG).

Given a target bit size, it finds a prime `p`, then searches the fixed-`a` family:
- `y² ≡ x³ - 3x + b (mod p)`
- Group order `N` is prime (cofactor h = 1)
- `N ≠ p` (not anomalous)

The search uses `ellsea(E, 1)` and benefits substantially from PARI `seadata`.

---

## Prerequisites

Install the required libraries (Ubuntu/Debian/WSL):

```bash
sudo apt install libpari-dev libsodium-dev pari-seadata cmake build-essential
```

If you prefer a local copy instead of the system package, extract the PARI `seadata` archive into `pari-data/seadata` at the repo root. The program auto-detects that directory at startup.

---

## Configuration

Open `main.cpp` and set the desired key size on this line:

```cpp
const int TARGET_BITS = 128;   // e.g. 64, 128, 192, 256
```

---

## Build

```bash
mkdir build && cd build
cmake ..
make
```

The binary will be at `build/app`.

---

## Run

```bash
./app
./app --verbose
```

Example output:

```
================================================
  ECC Parameter Generator  (libpari + libsodium)
  Target: 128-bit prime p, cofactor h=1, N prime, N≠p
  Curve:  y² ≡ x³ - 3x + b  (mod p)
  seadata: LOADED
================================================

[Searching] | Attempt #1 | p = 28134...  | Searching b (with a = -3)...

================================================
  *** SECURE CURVE FOUND ***
================================================
  Bits  : 128
  p     = 281341...
  a     = p - 3
  b     = 991234...
  N     = 281341...
  h     = 1  (N is prime, cofactor = 1)
  N ≠ p : YES (safe)
--------------------------------------------------
  isprime(N) : TRUE
  isprime(p) : TRUE
--------------------------------------------------
  Time elapsed: 2.341s
================================================
```
