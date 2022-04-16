/* kry.cpp
 * KRY 2021/22 projekt 2
 * Samuel Olekšák (xoleks00)
 */

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include <iostream>

#include "kry.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* ARGUMENT PARSING */
    if (argc < 2) {
        fprintf(stderr, "Missing mode argument, use  -g, -e, -d or -b.\n");
        return EXIT_FAILURE;
    }

    string mode = argv[1];

    if (mode == "-g") {
        if (argc != 3) {
            fprintf(stderr, "Incorrect number of arguments after -g, expected 1 (./kry -g B)\n");
            return EXIT_FAILURE;
        }

        char* endptr;

        int keySize = strtol(argv[2], &endptr, 10);

        if (*endptr != '\0' || keySize < 1) {
            fprintf(stderr, "Failed to parse key size parameter, expected positive integer\n");
            return EXIT_FAILURE;
        }

        generateKeys(keySize);
    }
    else if (mode == "-e") {
        if (argc != 5) {
            fprintf(stderr, "Incorrect number of arguments after -e, expected 3 (./kry -e E N M)\n");
            return EXIT_FAILURE;
        }

        char* endptr;

        int publicExponent = strtol(argv[2], &endptr, 16);

        if (*endptr != '\0' || publicExponent < 1) {
            fprintf(stderr, "Failed to parse public exponent parameter, expected positive integer\n");
            return EXIT_FAILURE;
        }

        int modulus = strtol(argv[3], &endptr, 16);

        if (*endptr != '\0' || modulus < 1) {
            fprintf(stderr, "Failed to parse modulus parameter, expected positive integer\n");
            return EXIT_FAILURE;
        }

        string message = argv[4];

        encrypt(publicExponent, modulus, message);
    }
    else if (mode == "-d") {
        if (argc != 5) {
            fprintf(stderr, "Incorrect number of arguments after -d, expected 3 (./kry -d D N C)\n");
            return EXIT_FAILURE;
        }

        char* endptr;

        int privateExponent = strtol(argv[2], &endptr, 16);

        if (*endptr != '\0' || privateExponent < 1) {
            fprintf(stderr, "Failed to parse private exponent parameter, expected positive integer\n");
            return EXIT_FAILURE;
        }

        int modulus = strtol(argv[3], &endptr, 16);

        if (*endptr != '\0' || modulus < 1) {
            fprintf(stderr, "Failed to parse modulus parameter, expected positive integer\n");
            return EXIT_FAILURE;
        }

        string cypher = argv[4];

        decrypt(privateExponent, modulus, cypher);
    }
    else if (mode == "-b") {
        if (argc != 3) {
            fprintf(stderr, "Incorrect number of arguments after -b, expected 1 (./kry -b N)\n");
            return EXIT_FAILURE;
        }

        mpz_t modulus;
        mpz_init(modulus);

        mpz_set_str(modulus, argv[2], 10);

        if (mpz_set_str(modulus, argv[2], 10) || mpz_cmp_si(modulus, 1) < 0) {
            fprintf(stderr, "Failed to parse modulus parameter, expected positive integer\n");
            return EXIT_FAILURE;
        }

        breakCypher(modulus);
    }
    else {
        fprintf(stderr, "Unknown mode, expected -g, -e, -d or -b.\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// prints both prime factors, modulus and, public exponent and private exponent
void generateKeys(int keySize) {

}

// prints encrypted cyphertext
void encrypt(int publicExponent, int modulus, string message) {

}

// prints decrypted message
void decrypt(int privateExponent, int modulus, string cypher) {

}

// prints prime factor
void breakCypher(mpz_t modulus) {
    mpz_t factors[2];

    mpz_init(factors[0]);
    mpz_init(factors[1]);

    PollardRho(modulus, factors);
    gmp_printf("%Zd\n", factors[0]);
}

// expects n > 0, returns array of two mpz_t
void fermatFactorization(mpz_t n, mpz_t *factors) {
    mpz_t a;
    mpz_t b;
    mpz_t square;
    mpf_t tempFloat;

    mpz_init(a);
    mpz_init(b);
    mpz_init(square);
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

        // square = a^2 - n
        mpz_mul(square, a, a);
        mpz_sub(square, square, n);

        while (mpz_perfect_square_p(square) == 0) {
            // a += 1
            mpz_add_ui(a, a, 1);

            // square = a^2 - n
            mpz_mul(square, a, a);
            mpz_sub(square, square, n);
        }

        //b = sqrt(a * a - n);
        mpz_mul(square, a, a);
        mpz_sub(square, square, n);
        mpf_set_z(tempFloat, square);
        mpf_sqrt(tempFloat, tempFloat);
        mpz_set_f(b, tempFloat);

        //factors[0] = a - b;
        mpz_sub(factors[0], a, b);

        //factors[1] = a + b;
        mpz_add(factors[1], a, b);
    }
}

// binary GCD algorithm
// https://en.wikipedia.org/wiki/Binary_GCD_algorithm
void gcd(mpz_t in_a, mpz_t in_b, mpz_t *out) {
    mpz_t a, b;

    mpz_init(a);
    mpz_init(b);

    // copy input params, so they are not mutated
    mpz_set(a, in_a);
    mpz_set(b, in_b);

    int doubles = 0;
 
    // if one of the factors is 0, gcd is the other one
    if (mpz_cmp_si(a, 0) == 0) {
        mpz_set(*out, b);
        return;
    }

    if (mpz_cmp_si(b, 0) == 0) {
        mpz_set(*out, a);
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
}

// out = base ^ 2 % modulus
void modulatedSquare(mpz_t base, mpz_t modulus, mpz_t *out) {
    mpz_t exponent;

    mpz_init(exponent);

    mpz_set_si(exponent, 2);
    mpz_set_si(*out, 1);
 
    while (mpz_cmp_si(exponent, 0) > 0) {
        if (!mpz_even_p(exponent)) {
            mpz_mul(*out, *out, base);
            mpz_mod(*out, *out, modulus);
        }

        // bit shift exponent left (divide by 2)
        mpz_fdiv_q_2exp(exponent, exponent, 1);

        mpz_mul(base, base, base);
        mpz_mod(base, base, modulus);
    }
}

// Pollard's Rho algorithm for integer factorization
// https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
void PollardRho(mpz_t n, mpz_t *out) {
    mpz_t x;
    mpz_t y;
    mpz_t divisor;
    mpz_t candidate;
    mpz_t temp;

    mpz_init(x);
    mpz_init(y);
    mpz_init(divisor);
    mpz_init(candidate);
    mpz_init(temp);

    gmp_randstate_t randomizer;
    gmp_randinit_default(randomizer);
    gmp_randseed_ui(randomizer, time(NULL));
 
    // one has no prime factors
    if (mpz_cmp_si(n, 1) == 0) {
        mpz_set(*out, n);
        return;
    }
 
    // if number is even, one of the factors is definitely 2
    if (mpz_even_p(n)) {
        mpz_set_si(*out, 2);
        return;
    }
 
    // random [2, N)
    mpz_sub_ui(n, n, 2u);
    mpz_urandomm(x, randomizer, n);
    mpz_add_ui(x, x, 2u);
    mpz_set(y, x);
    mpz_add_ui(n, n, 2u);
 
    // random [1, N)
    mpz_sub_ui(n, n, 1u);
    mpz_urandomm(candidate, randomizer, n);
    mpz_add_ui(candidate, candidate, 1u);
    mpz_add_ui(n, n, 1u);
 
    mpz_set_si(divisor, 1);

    while (mpz_cmp_si(divisor, 1) == 0) {
        // tortoise step
        modulatedSquare(x, n, &temp);
        mpz_add(temp, temp, candidate);
        mpz_add(temp, temp, n);
        mpz_mod(x, temp, n);
 
        // hare step
        modulatedSquare(y, n, &temp);
        mpz_add(temp, temp, candidate);
        mpz_add(temp, temp, n);
        mpz_mod(y, temp, n);

        modulatedSquare(y, n, &temp);
        mpz_add(temp, temp, candidate);
        mpz_add(temp, temp, n);
        mpz_mod(y, temp, n);
 
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
}
