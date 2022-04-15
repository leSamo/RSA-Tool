/* kry.h
 * KRY 2021/22 projekt 2
 * Samuel Olekšák (xoleks00)
 */

void generateKeys(int keySize);

void encrypt(int publicExponent, int modulus, string message);

void decrypt(int privateExponent, int modulus, string cypher);

void breakCypher(int modulus);
