/* kry.hpp
 * KRY 2021/22 projekt 2
 * Samuel Olekšák (xoleks00)
 */

void generateKeys(int keyBits);

void encrypt(mpz_t publicExponent, mpz_t modulus, mpz_t message);

void decrypt(mpz_t privateExponent, mpz_t modulus, mpz_t cypher);

void breakCypher(mpz_t modulus);

void FermatFactorization(mpz_t n, mpz_t *factors);

void gcd(mpz_t a, mpz_t b, mpz_t *out);

void PollardRho(mpz_t n, mpz_t *out);

void generateRandom(int bytes, mpz_t *out);

bool MillerRabinTest(mpz_t n);

bool isPrime(mpz_t n);

unsigned long int getRandomSeed();

void extendedEuclid(mpz_t a, mpz_t b, mpz_t *x, mpz_t *y, mpz_t *gcd);
