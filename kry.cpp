#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <iostream>

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
    }
    else if (mode == "-e") {
        if (argc != 5) {
            fprintf(stderr, "Incorrect number of arguments after -e, expected 3 (./kry -e E N M)\n");
            return EXIT_FAILURE;
        }
    }
    else if (mode == "-d") {
        if (argc != 5) {
            fprintf(stderr, "Incorrect number of arguments after -d, expected 3 (./kry -d D N C)\n");
            return EXIT_FAILURE;
        }
    }
    else if (mode == "-b") {
        if (argc != 3) {
            fprintf(stderr, "Incorrect number of arguments after -b, expected 1 (./kry -b N)\n");
            return EXIT_FAILURE;
        }
    }
    else {
        fprintf(stderr, "Unknown mode, expected -g, -e, -d or -b.\n");
        return EXIT_FAILURE;
    }
    /*
    mpz_t a,b,c;
    mpz_inits(a,b,c,NULL);

    mpz_set_str(a, "1234", 10);
    mpz_set_str(b,"-5678", 10); //Decimal base

    mpz_add(c,a,b);

    cout<<"\nThe exact result is:";
    mpz_out_str(stdout, 10, c); //Stream, numerical base, var
    cout<<endl;

    mpz_abs(c, c);
    cout<<"The absolute value result is:";
    mpz_out_str(stdout, 10, c);
    cout<<endl;

    cin.get();
    */

    return EXIT_SUCCESS;
}