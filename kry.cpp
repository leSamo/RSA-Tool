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

        gmp_printf("%Zd\n", modulus);

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

    fermatFactorization(modulus, factors);

    gmp_printf("%Zd %Zd\n", factors[0], factors[1]);
}

// expects n > 0
void fermatFactorization(mpz_t n, mpz_t *factors) {
    mpz_t a;
    mpz_t b;
    mpz_t tempInt;
    mpf_t tempFloat;
    mpf_t square;

    mpz_init(a);
    mpz_init(b);
    mpz_init(tempInt);
    mpf_init(tempFloat);
    mpf_init(square);

    if (mpz_even_p(n) != 0) {
        mpz_set_si(factors[0], 2);
        mpz_cdiv_q(factors[1], n, factors[0]);
    }

    mpf_set_z(tempFloat, n);
    mpf_sqrt(tempFloat, tempFloat);
    mpf_ceil(tempFloat, tempFloat);
    mpz_set_f(a, tempFloat);

    mpz_mul(square, a, a);
    mpz_sub(square, square, n);

    gmp_printf("%Ff\n", tempFloat);
}

    /*
    ...

    while (int(square) != square) {
        ++a;

        square = (a * a - n);
    }

    b = sqrt(a * a - n);

    factors[0] = a - b;
    factors[1] = a + b;
    */
