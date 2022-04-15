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

        mpz_t a,b,c;
        mpz_init(a);
        mpz_init(b);
        mpz_init(c);

        mpz_set_si(a, 7);
        mpz_set_si(b, 14);

        gcd(a, b, &c);
        gmp_printf("%Zd\n", c);

        //breakCypher(modulus);
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

    fermatFactorization(modulus, factors);
}

// expects n > 0
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
    // else n is off
    else {
        // a = ceil(sqrt(n))
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

void gcd(mpz_t a, mpz_t b, mpz_t *out) {
    int doubles = 0;
 
    // if one of the factors is 0, gcd is the second one
    if (mpz_cmp_si(a, 0) == 0) {
        mpz_set(*out, b);
        return;
    }
    if (mpz_cmp_si(b, 0) == 0) {
        mpz_set(*out, a);
        return;
    }

    // a and b is even divide both by 2 (right bit-shift)
    for (; mpz_even_p(a) && mpz_even_p(b); ++doubles) {
        mpz_fdiv_q_2exp(a, a, 1);
        mpz_fdiv_q_2exp(b, b, 1);
    }

    // dividing a by 2 (right bit-shift) until a becomes odd
    while (mpz_even_p(a)) {
        mpz_fdiv_q_2exp(a, a, 1);
    }

    // From here on, 'a' is always odd.
    do {
        // If b is even, remove all factor of 2 in b
        while (mpz_even_p(b)) {
            mpz_fdiv_q_2exp(b, b, 1);
        }
 
        // Now a and b are both odd.
        //  Swap if necessary so a <= b,
        //   then set b = b - a (which is even).
        if (mpz_cmp(a, b) > 0) {
            mpz_swap(b, a);
        }

        mpz_sub(b, b, a);
    } while (mpz_cmp_si(b, 0) != 0);

    mpz_mul_2exp(*out, a, doubles);
}
