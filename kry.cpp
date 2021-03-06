/* kry.cpp
 * KRY 2021/22 projekt 2
 * Samuel Olekšák (xoleks00)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>

#include <gmp.h>

#include "kry.hpp"

#define MILLER_RABIN_ITERATIONS 40

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Missing mode argument, use  -g, -e, -d or -b.\n");
        return EXIT_FAILURE;
    }

    char *mode = argv[1];

    // Key generation mode: ./kry -g B
    // B (dec) - key size in bits
    // --------------------
    // Outputs: P Q N E D
    // P, Q - random integers
    // N - public modulo
    // E - public exponent
    // D - private exponent
    if (strcmp(mode, "-g") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Incorrect number of arguments after -g, expected 1 (./kry -g B)\n");
            return EXIT_FAILURE;
        }

        char* endptr;

        int keySize = strtol(argv[2], &endptr, 10);

        if (*endptr != '\0' || keySize < 5) {
            fprintf(stderr, "Failed to parse key size parameter, expected positive integer larger than 4\n");
            return EXIT_FAILURE;
        }

        generateKeys(keySize);
    }
    // Message encryption mode ./kry -e E N M
    // E (hex) - public exponent
    // N (hex) - public modulo
    // M (hex) - plaintext message to be encrypted
    // --------------------
    // Outputs: C
    // C - cyphertext of input message
    else if (strcmp(mode, "-e") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Incorrect number of arguments after -e, expected 3 (./kry -e E N M)\n");
            return EXIT_FAILURE;
        }

        mpz_t publicExponent, modulus, message;

        mpz_inits(publicExponent, modulus, message, 0);
        
        if (mpz_set_str(publicExponent, argv[2], 0) || mpz_cmp_si(publicExponent, 1) < 0) {
            fprintf(stderr, "Failed to parse public exponent parameter, expected positive integer\n");
            
            mpz_clears(publicExponent, modulus, message, 0);
            return EXIT_FAILURE;
        }

        if (mpz_set_str(modulus, argv[3], 0) || mpz_cmp_si(modulus, 1) < 0) {
            fprintf(stderr, "Failed to parse modulus parameter, expected positive integer\n");

            mpz_clears(publicExponent, modulus, message, 0);
            return EXIT_FAILURE;
        }

        if (mpz_set_str(message, argv[4], 0) || mpz_cmp_si(message, 1) < 0) {
            fprintf(stderr, "Failed to parse message parameter, expected positive integer\n");
            
            mpz_clears(publicExponent, modulus, message, 0);
            return EXIT_FAILURE;
        }

        encrypt(publicExponent, modulus, message);

        mpz_clears(publicExponent, modulus, message, 0);
    }
    // Message decryption mode ./kry -d D N C
    // D (hex) - private exponent
    // N (hex) - public modulo
    // C (hex) - cyphertext to be decrypted
    // --------------------
    // Outputs: M
    // M - plaintext of input cypher
    else if (strcmp(mode, "-d") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Incorrect number of arguments after -d, expected 3 (./kry -d D N C)\n");
            return EXIT_FAILURE;
        }

        mpz_t privateExponent, modulus, cypher;

        mpz_inits(privateExponent, modulus, cypher, 0);

        if (mpz_set_str(privateExponent, argv[2], 0) || mpz_cmp_si(privateExponent, 1) < 0) {
            fprintf(stderr, "Failed to parse private exponent parameter, expected positive integer\n");

            mpz_clears(privateExponent, modulus, cypher, 0);
            return EXIT_FAILURE;
        }

        if (mpz_set_str(modulus, argv[3], 0) || mpz_cmp_si(modulus, 1) < 0) {
            fprintf(stderr, "Failed to parse modulus parameter, expected positive integer\n");
            
            mpz_clears(privateExponent, modulus, cypher, 0);
            return EXIT_FAILURE;
        }

        if (mpz_set_str(cypher, argv[4], 0) || mpz_cmp_si(cypher, 1) < 0) {
            fprintf(stderr, "Failed to parse cypher parameter, expected positive integer\n");
            
            mpz_clears(privateExponent, modulus, cypher, 0);
            return EXIT_FAILURE;
        }

        decrypt(privateExponent, modulus, cypher);

        mpz_clears(privateExponent, modulus, cypher, 0);
    }
    // Key factorization mode ./kry -b N
    // N (hex) - public modulo
    // --------------------
    // Outputs: P
    // P - one of the prime factors of input N
    else if (strcmp(mode, "-b") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Incorrect number of arguments after -b, expected 1 (./kry -b N)\n");
            return EXIT_FAILURE;
        }

        mpz_t modulus;
        mpz_init(modulus);

        if (mpz_set_str(modulus, argv[2], 0) || mpz_cmp_si(modulus, 1) < 0) {
            fprintf(stderr, "Failed to parse modulus parameter, expected positive integer\n");
            
            mpz_clear(modulus);
            return EXIT_FAILURE;
        }

        breakCypher(modulus);

        mpz_clear(modulus);
    }
    else {
        fprintf(stderr, "Unknown mode, expected -g, -e, -d or -b.\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// prints both prime factors, modulus and, public exponent and private exponent
void generateKeys(int keyBits) {
    mpz_t p, q, n, phi, e, d, gcd, a, b;

    mpz_inits(p, q, n, phi, e, d, gcd, a, b, 0);

    int primeBits = keyBits % 2 == 0 ? keyBits / 2 : (keyBits + 1) / 2;

    // p and q primes generation
    do {
        do {
            generateRandom(primeBits, &p);
        }
        while (!isPrime(p));

        do {
            generateRandom(primeBits, &q);
        }
        while (!isPrime(q));

        // public modulo n = p * q
        mpz_mul(n, p, q);
    }
    // check whether product of two primes has correct bit size
    while (static_cast<int>(mpz_sizeinbase(n, 2)) != keyBits);

    // phi(n) = (p - 1) * (q - 1)
    mpz_sub_ui(p, p, 1u);
    mpz_sub_ui(q, q, 1u);
    mpz_mul(phi, p, q);
    mpz_add_ui(p, p, 1u);
    mpz_add_ui(q, q, 1u);

    gmp_randstate_t randomizer;
    gmp_randinit_default(randomizer);
    gmp_randseed_ui(randomizer, getRandomSeed());

    do {
        // generate random e such that 1 < e < phi(n)
        mpz_sub_ui(phi, phi, 2u);
        mpz_urandomm(e, randomizer, phi);
        mpz_add_ui(e, e, 2u);
        mpz_add_ui(phi, phi, 2u);

        // a and b are Bézout coefficients
        extendedEuclid(e, phi, &a, &b, &gcd);
    }
    // generate e until gcd(e, phi) = 1
    while (mpz_cmp_si(gcd, 1) != 0);

    // calculate modular multiplicative inverse
    // d = (a % phi + phi) % phi
    mpz_mod(d, a, phi);
    mpz_add(d, d, phi);
    mpz_mod(d, d, phi);

    gmp_printf("0x%Zx 0x%Zx 0x%Zx 0x%Zx 0x%Zx\n", p, q, n, e, d);

    gmp_randclear(randomizer);
    mpz_clears(p, q, n, phi, e, d, gcd, a, b, 0);
}

// prints encrypted cyphertext
void encrypt(mpz_t publicExponent, mpz_t modulus, mpz_t message) {
    mpz_t cypher;
    mpz_init(cypher);

    // c = m^e % n
    mpz_powm(cypher, message, publicExponent, modulus);

    gmp_printf("0x%Zx\n", cypher);

    mpz_clear(cypher);
}

// prints decrypted message
void decrypt(mpz_t privateExponent, mpz_t modulus, mpz_t cypher) {
    mpz_t message;
    mpz_init(message);

    // m = c^d % n
    mpz_powm(message, cypher, privateExponent, modulus);

    gmp_printf("0x%Zx\n", message);

    mpz_clear(message);
}

// prints public modulus prime factor
void breakCypher(mpz_t modulus) {
    mpz_t factors[2], remainder;

    mpz_inits(factors[0], factors[1], remainder, 0);

    // check if 2 is factor, optimizes following trivial division
    if (mpz_even_p(modulus)) {
        printf("0x2\n");

        mpz_clears(factors[0], factors[1], remainder, 0);
        return;
    }

    // trivial division for first 1 milion factors
    for (unsigned int divisor = 3u; divisor < 1'000'000u; divisor += 2) {
        mpz_cdiv_r_ui(remainder, modulus, divisor);

        if (mpz_cmp_si(remainder, 0) == 0) {
            printf("0x%x\n", divisor);

            mpz_clears(factors[0], factors[1], remainder, 0);
            return;
        }
    }

    // Pollard Rho factorization if first 1 milion numbers are not factors
    PollardRho(modulus, factors);

    gmp_printf("0x%Zx\n", factors[0]);

    mpz_clears(factors[0], factors[1], remainder, 0);
}

// generates random number with specific bit size
// MSB must be one - result is in the range [2^(n-1),2^n)
void generateRandom(int bits, mpz_t *out) {
    mpz_t randomPrime;
    mpz_t lowerBound;

    mpz_inits(randomPrime, lowerBound, 0);

    gmp_randstate_t randomizer;
    gmp_randinit_default(randomizer);
    gmp_randseed_ui(randomizer, getRandomSeed());

    // get random number in range [0, 2^(n-1)) and add 2^(n-1) to it
    mpz_urandomb(randomPrime, randomizer, bits - 1);
    mpz_ui_pow_ui(lowerBound, 2, bits - 1);
    mpz_add(randomPrime, randomPrime, lowerBound);

    mpz_set(*out, randomPrime);

    gmp_randclear(randomizer);

    mpz_clears(randomPrime, lowerBound, 0);
}

// Extended Euclid algorithm for calculating GCD and Bézout coefficients
// https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
// ax + by = gcd(a, b)
void extendedEuclid(mpz_t a, mpz_t b, mpz_t *x, mpz_t *y, mpz_t *gcd) {
    // from basic euclidean algorithm we know that if a reaches 0 then b = gcd(a, b)
    // from Bézout's identity ax + by = gcd(a, b) we know that if b = gcd(a, b) then y must be 1
    // as a = 0, x can be basically anything
    if (mpz_cmp_si(a, 0) == 0) {
        mpz_set_si(*x, 0);
        mpz_set_si(*y, 1);
        mpz_set(*gcd, b);
    }
    else {
        mpz_t xNested, yNested, gcdNested, temp;
        mpz_inits(xNested, yNested, gcdNested, temp, 0);

        // find the remainder to be used as a base in the next round
        // b' = b % a
        mpz_mod(temp, b, a);

        extendedEuclid(temp, a, &xNested, &yNested, &gcdNested);

        // x = y' - b / a * x'
        mpz_tdiv_q(temp, b, a);
        mpz_mul(temp, temp, xNested);
        mpz_sub(*x, yNested, temp);

        // y = x'
        mpz_set(*y, xNested);

        mpz_set(*gcd, gcdNested);

        mpz_clears(xNested, yNested, gcdNested, temp, 0);
    }
}


// expects n > 0, returns both factors as array of two mpz_t
// https://en.wikipedia.org/wiki/Fermat%27s_factorization_method
void FermatFactorization(mpz_t n, mpz_t *factors) {
    mpz_t a, b, bSquared;
    mpf_t tempFloat;

    mpz_inits(a, b, bSquared, 0);
    mpf_init(tempFloat);

    // if n is even
    if (mpz_even_p(n)) {
        // set factors to 2 and n / 2
        mpz_set_si(factors[0], 2);
        mpz_cdiv_q(factors[1], n, factors[0]);
    }
    // else n is odd
    else {
        // a = √n rounded up
        mpf_set_z(tempFloat, n);
        mpf_sqrt(tempFloat, tempFloat);
        mpf_ceil(tempFloat, tempFloat);
        mpz_set_f(a, tempFloat);

        // bSquared = a^2 - n
        mpz_mul(bSquared, a, a);
        mpz_sub(bSquared, bSquared, n);

        while (mpz_perfect_square_p(bSquared) == 0) {
            // a += 1
            mpz_add_ui(a, a, 1);

            // bSquared = a^2 - n
            mpz_mul(bSquared, a, a);
            mpz_sub(bSquared, bSquared, n);
        }

        //b = sqrt(a * a - n);
        mpz_mul(bSquared, a, a);
        mpz_sub(bSquared, bSquared, n);
        mpf_set_z(tempFloat, bSquared);
        mpf_sqrt(tempFloat, tempFloat);
        mpz_set_f(b, tempFloat);

        //factors[0] = a - b;
        mpz_sub(factors[0], a, b);

        //factors[1] = a + b;
        mpz_add(factors[1], a, b);
    }

    mpz_clears(a, b, bSquared, 0);
    mpf_clear(tempFloat);
}

// binary GCD algorithm
// https://en.wikipedia.org/wiki/Binary_GCD_algorithm
void gcd(mpz_t in_a, mpz_t in_b, mpz_t *out) {
    mpz_t a, b;

    mpz_inits(a, b, 0);

    // copy input params, so they are not mutated
    mpz_set(a, in_a);
    mpz_set(b, in_b);

    int doubles = 0;
 
    // if one of the factors is 0, gcd is the other one
    if (mpz_cmp_si(a, 0) == 0) {
        mpz_set(*out, b);

        mpz_clears(a, b, 0);
        return;
    }

    if (mpz_cmp_si(b, 0) == 0) {
        mpz_set(*out, a);

        mpz_clears(a, b, 0);
        return;
    }

    // a and b is even divide both by 2 (right bit-shift)
    // remember how many times we divided in doubles variable
    // so that we can revert this
    for (; mpz_even_p(a) && mpz_even_p(b); ++doubles) {
        mpz_fdiv_q_2exp(a, a, 1);
        mpz_fdiv_q_2exp(b, b, 1);
    }

    // dividing a by 2 (right bit-shift) until a becomes odd
    while (mpz_even_p(a)) {
        mpz_fdiv_q_2exp(a, a, 1);
    }

    do {
        // if b is even, remove all factor of 2 in b
        // so that both a an b are odd
        while (mpz_even_p(b)) {
            mpz_fdiv_q_2exp(b, b, 1);
        }
 
        // make sure larger operand is saved in b
        // swap a and b if necessary
        if (mpz_cmp(a, b) > 0) {
            mpz_swap(b, a);
        }

        // set b = b - a, odd - odd = even, so b is now even
        mpz_sub(b, b, a);
    // repeat until b = 0, GCD is saved in a
    } while (mpz_cmp_si(b, 0) != 0);

    // bitshift left by how many times we divided
    // both operands by two
    mpz_mul_2exp(*out, a, doubles);

    mpz_clears(a, b, 0);
}

// Pollard's Rho algorithm for integer factorization
// https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
void PollardRho(mpz_t n, mpz_t *out) {
    mpz_t x, y, divisor, candidate, temp;

    mpz_inits(x, y, divisor, candidate, temp, 0);

    // one has no prime factors
    if (mpz_cmp_si(n, 1) == 0) {
        mpz_set(*out, n);

        mpz_clears(x, y, divisor, candidate, temp, 0);
        return;
    }
 
    // if number is even, one of the factors is definitely 2
    if (mpz_even_p(n)) {
        mpz_set_si(*out, 2);

        mpz_clears(x, y, divisor, candidate, temp, 0);
        return;
    }

    gmp_randstate_t randomizer;
    gmp_randinit_default(randomizer);
    gmp_randseed_ui(randomizer, getRandomSeed());
 
    // x = y = random [2, n)
    mpz_sub_ui(n, n, 2u);
    mpz_urandomm(x, randomizer, n);
    mpz_add_ui(x, x, 2u);
    mpz_set(y, x);
    mpz_add_ui(n, n, 2u);
 
    // candidate = random [1, n)
    mpz_sub_ui(n, n, 1u);
    mpz_urandomm(candidate, randomizer, n);
    mpz_add_ui(candidate, candidate, 1u);
    mpz_add_ui(n, n, 1u);
 
    mpz_set_si(divisor, 1);

    while (mpz_cmp_si(divisor, 1) == 0) {
        // baby step
        // we apply polynomial x' = (x^2 + 1) % n once
        // x = (x^2 % n + candidate + n) % n
        mpz_powm_ui(temp, x, 2, n);
        mpz_add(temp, temp, candidate);
        mpz_add(temp, temp, n);
        mpz_mod(x, temp, n);
 
        // giant step
        // we apply polynomial x' = (x^2 + 1) % n once
        // y = (y^2 % n + candidate + n) % n
        mpz_powm_ui(temp, y, 2, n);
        mpz_add(temp, temp, candidate);
        mpz_add(temp, temp, n);
        mpz_mod(y, temp, n);

        mpz_powm_ui(temp, y, 2, n);
        mpz_add(temp, temp, candidate);
        mpz_add(temp, temp, n);
        mpz_mod(y, temp, n);

        // divisor = gcd(|x - y|, n)
        mpz_sub(temp, x, y);
        mpz_abs(temp, temp);
        gcd(temp, n, &divisor);
 
        // try again if fail
        if (mpz_cmp(divisor, n) == 0) {
            PollardRho(n, out);
        }
    }
    
    // sucessfully found factor
    mpz_set(*out, divisor);

    gmp_randclear(randomizer);
    mpz_clears(x, y, divisor, candidate, temp, 0);
}

// checks for primality using multiple rounds of Miller Rabin test
bool isPrime(mpz_t n) {
    for (int i = 0; i < MILLER_RABIN_ITERATIONS; ++i) {
        if (!MillerRabinTest(n)) {
            return false;
        }
    }

    return true;
}

// Miller Rabin primality test
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
// definitely returns false if number is compount
// returns true when number is probably a prime
bool MillerRabinTest(mpz_t n) {
    mpz_t nMinus1, m, a, b;

    mpz_inits(nMinus1, m, a, b, 0);

    mpz_sub_ui(nMinus1, n, 1u);
    mpz_set(m, nMinus1);

    // n - 1 = 2^k * m
    while (mpz_even_p(m)) {
        mpz_fdiv_q_2exp(m, m, 1);
    }

    gmp_randstate_t randomizer;
    gmp_randinit_default(randomizer);
    gmp_randseed_ui(randomizer, getRandomSeed());

    // random a between [2, n-2]
    mpz_sub_ui(n, n, 3u);
    mpz_urandomm(a, randomizer, n);
    mpz_add_ui(a, a, 2u);
    mpz_add_ui(n, n, 3u);

    // b0 = a^m mod n
    mpz_powm(b, a, m, n);

    // if b0 = +- 1 then n is probably prime
    if (mpz_cmp(b, nMinus1) == 0 || mpz_cmp_si(b, 1) == 0) {
        gmp_randclear(randomizer);
        mpz_clears(nMinus1, m, a, b, 0);

        return true; 
    }

    // b' = b^2 mod n
    // if b' = 1 then n is definitely not a prime
    // else if b' = -1 the n is probably a prime
    // else generate new b' and try again
    while (mpz_cmp(m, nMinus1) != 0) {
        mpz_powm_ui(b, b, 2u, n);

        if (mpz_cmp_si(b, 1) == 0) {
            gmp_randclear(randomizer);
            mpz_clears(nMinus1, m, a, b, 0);

            return false;
        }

        if (mpz_cmp(b, nMinus1) == 0) {
            gmp_randclear(randomizer);
            mpz_clears(nMinus1, m, a, b, 0);

            return true;
        }

        mpz_mul_si(m, m, 2);
    }

    gmp_randclear(randomizer);
    mpz_clears(nMinus1, m, a, b, 0);

    return false;
}

// returns a random number in the range of unsigned long int using /dev/urandom
unsigned long int getRandomSeed() {
    unsigned long int randomSeed;

    if (syscall(SYS_getrandom, &randomSeed, sizeof(unsigned long int), 0) == sizeof(unsigned long int)) {
        return randomSeed;
    }
    else {
        fprintf(stderr, "Failed to generate random number using getrandom() system call\n");
        exit(EXIT_FAILURE);
    }
}
