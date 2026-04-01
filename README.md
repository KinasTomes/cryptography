# ECC Parameter Generator

Generates cryptographically secure elliptic curve parameters using **libpari** (number theory / SEA point counting) and **libsodium** (CSPRNG).

Given a target bit size, it finds a prime `p`, then searches for curve coefficients `a`, `b` such that:
- `y² ≡ x³ + ax + b (mod p)`
- Group order `N` is prime (cofactor h = 1)
- `N ≠ p` (not anomalous)

---

## Prerequisites

Install the required libraries (Ubuntu/Debian/WSL):

```bash
sudo apt install libpari-dev libsodium-dev cmake build-essential
```

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
```

Example output:

```
================================================
  ECC Parameter Generator  (libpari + libsodium)
  Target: 128-bit prime p, cofactor h=1, N prime, N≠p
  Curve:  y² ≡ x³ + ax + b  (mod p)
================================================

[Searching] | Attempt #1 | p = 28134...  | Searching (a,b)...

================================================
  *** SECURE CURVE FOUND ***
================================================
  Bits  : 128
  p     = 281341...
  a     = 174823...
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
