Tool for RSA key generation, encryption, decryption and weak key cracking.


Modes:
  -g - generates public and private RSA keys
    usage: ./kry -g B
    output: P Q N E D

  -e - encrypts input number using provided public exponent and modulus
    usage: ./kry -e E N M
    output: C

  -d - decripts input number using provided private exponent and modulus
    usage: ./kry -d D N C
    output: M

  -b - cracks provided weak modulo and outputs its prime factors
    usage: ./kry -b N
    output: P

  B - key size in bits
  P, Q - prime factors
  N - modulus
  E - public exponent
  D - private exponent
  M - plaintext number
  C - cyphertext number


Example usage:
  Generating 96-bit key

    $ ./kry -g 96
    > 0xf13ad36fd55f 0xf53ca2fe8ec9 0xe7166fbaae2083b7dd6b3997 0x2becfc96fa007153e0b88a33 0x36df69f8529fa8816bae6f9b

    First two numbers represent randomly generated primes used in key creation.
    Third and fourth number are public and private exponents. Last number represents modulus.
    Public key: [0x2becfc96fa007153e0b88a33, 0xe7166fbaae2083b7dd6b3997]
    Private key: [0x36df69f8529fa8816bae6f9b, 0xe7166fbaae2083b7dd6b3997]

  Encrypting "message" number 0xabcd with public key generated in previous step

    $ ./kry -e 0x2becfc96fa007153e0b88a33 0xe7166fbaae2083b7dd6b3997 0xabcd
    > 0xa67dd43054bc1351556cd4b0

  Decrypting cypher from previous step using private key from the first step

    $ ./kry -d 0x36df69f8529fa8816bae6f9b 0xe7166fbaae2083b7dd6b3997 0xa67dd43054bc1351556cd4b0
    > 0xabcd

  Try to crack the key by finding one of the prime factors of modulus
    $ ./kry -b 0xe7166fbaae2083b7dd6b3997
    > 0xf13ad36fd55f

    We found prime factor of modulus in just a few moments which is why key with only 96-bit lenght is nowadays insecure.


Technical details:
  Tool is implemented in C++ using library for manipulating with large numbers GMP. Primality of generated primes is ensured
  using 40 iterations of Miller-Rabin algorithm. Finding of modular inverse used in private exponent generation is implemented
  using extended Euclidean algorithm. Modulus factorization is done using Pollard's Rho algorithm, which means most personal
  computers should crack modulus of 96-bits or less in just a few minutes. 
