#include <stdio.h>
#include <immintrin.h>
#include "ntru.h"

static void poly_rand_ternary(poly p, mpz_t val) {
    uint64_t j;
    mpz_t neg_val;

    mpz_init(neg_val);
    mpz_neg(neg_val, val);
    poly_zero(p);
    for (uint64_t i = 0; i <= D; i++) {
        do {
            _rdrand64_step((unsigned long long int *)&j);
            j %= N;
        } while (mpz_sgn(p[j]));
        if (i <= D / 2)
            mpz_set(p[j], val);
        else
            mpz_set(p[j], neg_val);
    }

    mpz_clear(neg_val);
}

void ntru_keygen(poly f, poly h, mpz_t P, mpz_t Q) {
    poly f_q = poly_new();
    do {
        poly_rand_ternary(f, P);
        mpz_add_ui(f[0], f[0], 1);
    } while (!poly_invert(f_q, f, Q));

    poly g = poly_new();
    do {
        poly_rand_ternary(g, P);
    } while (!poly_invert(NULL, g, Q));
    poly_conv(h, f_q, g, Q);

    poly_clear(f_q); poly_clear(g);
}

void ntru_encrypt(poly c, poly m, poly h, mpz_t Q) {
    mpz_t ONE;

    mpz_init_set_ui(ONE, 1);
    poly_zero(c);
    poly_rand_ternary(c, ONE);
    poly_conv(c, c, h, Q); // c = r*h (mod Q)
    poly_add(c, c, m, Q); // c = h*r + m (mod Q)

    mpz_clear(ONE);
}

void ntru_decrypt(poly m, poly c, poly f, mpz_t P, mpz_t Q) {
    poly_conv(m, f, c, Q); //m = f*c (mod Q)
    poly_center(m, Q);
    poly_mod(m, P); // m = m (mod P)
}
