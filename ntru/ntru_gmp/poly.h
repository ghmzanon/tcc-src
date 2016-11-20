#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

#define N 512

typedef mpz_t* poly;

poly poly_new();
void poly_zero(poly p);
void poly_clear(poly p);
void poly_print(poly p);
void poly_rand(poly p,  mpz_t m);
uint16_t poly_degree(poly p);
void poly_swap(poly *p, poly *q);
void poly_mod(poly p, mpz_t m);
void poly_conv(poly r, poly p, poly q, mpz_t m);
void poly_add(poly r, poly p, poly q, mpz_t m);
void poly_set(poly p, poly q);
bool poly_cmp_eq(poly p, poly q);
void poly_center(poly p, mpz_t m);
bool poly_invert(poly r, poly p, mpz_t m);
