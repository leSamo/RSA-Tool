/* kry.h
 * KRY 2021/22 projekt 2
 * Samuel Olekšák (xoleks00)
 */

void generateKeys(int keySize);

void encrypt(int publicExponent, int modulus, std::string message);

void decrypt(int privateExponent, int modulus, std::string cypher);

void breakCypher(mpz_t modulus);

void fermatFactorization(mpz_t n, mpz_t *factors);

void gcd(mpz_t a, mpz_t b, mpz_t *out);

void modulatedSquare(mpz_t base, mpz_t modulus, mpz_t *out);

void PollardRho(mpz_t n, mpz_t *out);
