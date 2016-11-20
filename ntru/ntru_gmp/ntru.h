#include "poly.h"

#define NTRU_TEST
#define NUM_TESTS 100

// #define PRIME_P "263"
// #define PRIME_Q "17179869143"
#define PRIME_P "13"
#define PRIME_Q "10007"

#define D 2*(N/3)

void ntru_keygen(poly f, poly h);
void ntru_encrypt(poly c, poly m, poly h);
void ntru_decrypt(poly m, poly c, poly f);
