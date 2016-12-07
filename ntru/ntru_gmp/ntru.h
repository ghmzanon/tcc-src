#include "params.h"
#include "poly.h"


void ntru_keygen(poly f, poly h, mpz_t P, mpz_t Q);
void ntru_encrypt(poly c, poly m, poly h, mpz_t Q);
void ntru_decrypt(poly m, poly c, poly f, mpz_t P, mpz_t Q);
