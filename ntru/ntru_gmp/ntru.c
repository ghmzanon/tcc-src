#include <stdio.h>
#include <immintrin.h>
#include "ntru.h"

#ifdef NTRU_TEST
#include <string.h>
#include <time.h>
#include <stdlib.h>
#endif

static void poly_rand_ternary(poly p, mpz_t val) {
    FILE *file = fopen("/dev/urandom", "r");
    uint64_t j;
    mpz_t neg_val;
    mpz_init(neg_val);
    mpz_neg(neg_val, val);

    poly_zero(p);

    for (uint64_t i = 0; i <= D; i++) {
        do {
            /* _rdrand64_step((unsigned long long int *)&j); */
            fread(&j, sizeof(uint64_t), 1, file);
            j %= N;
        } while (mpz_sgn(p[j]));
        if (i <= D / 2)
            mpz_set(p[j], val);
        else
            mpz_set(p[j], neg_val);
    }
    fclose(file);
    mpz_clear(neg_val);
}

void ntru_keygen(poly f, poly h) {
    mpz_t P, Q;
    mpz_init_set_str(P, PRIME_P, 10);
    mpz_init_set_str(Q, PRIME_Q, 10);

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

    mpz_clear(P);
    mpz_clear(Q);
    poly_clear(f_q);
    poly_clear(g);
}

void ntru_encrypt(poly c, poly m, poly h) {
    mpz_t Q, ONE;
    mpz_init_set_str(Q, PRIME_Q, 10);
    mpz_init_set_ui(ONE, 1);

    poly_zero(c);
    poly_rand_ternary(c, ONE);

    poly_conv(c, c, h, Q); // c = r*h (mod Q)
    poly_add(c, c, m, Q); // c = h*r + m (mod Q)

    mpz_clear(Q);
    mpz_clear(ONE);
}

void ntru_decrypt(poly m, poly c, poly f) {
    mpz_t P, Q;
    mpz_init_set_str(P, PRIME_P, 10);
    mpz_init_set_str(Q, PRIME_Q, 10);

    poly_conv(m, f, c, Q); //m = f*c (mod Q)
    poly_center(m, Q);
    poly_mod(m, P); // m = m (mod P)

    mpz_clear(P); mpz_clear(Q);
}


#ifdef NTRU_TEST
void poly_rand_binary(poly p, uint64_t n) {
    FILE *file = fopen("/dev/urandom", "r");
    uint16_t bit;
    poly_zero(p);

    for (uint16_t i = 0; i < n; i++) {
        /* _rdrand16_step((unsigned short *)&bit); */
        fread(&bit, sizeof(uint16_t), 1, file);
        mpz_set_ui(p[i], bit & 0x1);
    }
    fclose(file);
}

void ntru_homomorphic_add_test() {
    uint32_t err = 0;
    poly f = poly_new();  // private key
    poly h = poly_new();  // public key
    poly m = poly_new();
    poly m_sum = poly_new();
    poly c = poly_new();
    poly c_sum = poly_new();

    mpz_t Q, P;
    mpz_init_set_str(Q, PRIME_Q, 10);
    mpz_init_set_str(P, PRIME_P, 10);

    mpz_t i;
    for (int n = 0 ; n < NUM_TESTS; n++) {
        /* printf("Generating keys...\n"); */
        ntru_keygen(f, h);
        /* printf("Generated!\n"); */

        /* printf("Generating message 1...\n"); */
        poly_rand_binary(m_sum, N);
        /* printf("Generated!\n"); */
        /* printf("Encrypting message 1...\n"); */
        ntru_encrypt(c_sum, m_sum, h);
        /* printf("Encrypted!\n"); */

        mpz_init_set_ui(i, 2);
        while (mpz_cmp(i, P) < 0) {
            /* gmp_printf("Generating message %Zd...\n", i); */
            poly_rand_binary(m, N);
            /* printf("Generated!\n"); */
            /* gmp_printf("Encrypting message %Zd...\n", i); */
            ntru_encrypt(c, m, h);
            /* printf("Encrypted!\n"); */
            /* printf("Performing homomorphic addition...\n"); */
            poly_add(c_sum, c_sum, c, Q);
            /* printf("Added!\n"); */
            /* printf("Adding plain message...\n"); */
            poly_add(m_sum, m_sum, m, Q);
            /* printf("Added!\n"); */
            mpz_add_ui(i, i, 1);
        }

        /* printf("Decrypting homomorphic added message...\n"); */
        ntru_decrypt(c, c_sum, f);
        /* printf("Decrypted!\n"); */

        /* poly_print(m_sum); */
        /* poly_print(c); */
        /* if (poly_cmp_eq(m_sum, c)) */
        /*     printf("OK\n"); */
        if (!poly_cmp_eq(m_sum, c))
            err++;
    }
    if (err == 0)
        printf("ALL TESTS HAVE PASSED\n");
    else
        printf("FAILED: %d\n", err);

    mpz_clear(i);
    mpz_clear(Q); mpz_clear(P);
    poly_clear(f); poly_clear(h);
    poly_clear(c); poly_clear(c_sum);
    poly_clear(m); poly_clear(m_sum);
}

void ntru_homomorphic_mul_test() {
    uint32_t err = 0;
    poly f = poly_new();  // private key
    poly h = poly_new();  // public key
    poly m = poly_new();
    poly m_mul = poly_new();
    poly c = poly_new();
    poly c_mul = poly_new();

    mpz_t Q, P;
    mpz_init_set_str(Q, PRIME_Q, 10);
    mpz_init_set_str(P, PRIME_P, 10);

    for (int n = 0 ; n < NUM_TESTS; n++) {
        ntru_keygen(f, h);

        poly_rand_binary(m_mul, N/2);
        ntru_encrypt(c_mul, m_mul, h);

        poly_rand_binary(m, N/2);
        ntru_encrypt(c, m, h);
        poly_conv(c_mul, c_mul, c, Q);
        poly_conv(m_mul, m_mul, m, Q);

        poly_conv(f, f, f, Q);
        ntru_decrypt(c, c_mul, f);

        if (!poly_cmp_eq(m_mul, c)) {
            err++;
        }
    }
    if (err == 0)
        printf("ALL TESTS HAVE PASSED\n");
    else
        printf("FAILED: %d\n", err);

    mpz_clear(Q); mpz_clear(P);
    poly_clear(f); poly_clear(h);
    poly_clear(c); poly_clear(c_mul);
    poly_clear(m); poly_clear(m_mul);
}

int main() {
    ntru_homomorphic_add_test();
    /* ntru_homomorphic_mul_test(); */
    return 0;
}
#endif
