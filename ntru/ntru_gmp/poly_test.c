#include <stdio.h>
#include <immintrin.h>
#include "poly.h"

void poly_rand(poly p, mpz_t m) {
    gmp_randstate_t state;
    uint64_t seed;

    gmp_randinit_mt(state);
    _rdrand64_step((unsigned long long int *)&seed);
    gmp_randseed_ui(state, seed);
    poly_zero(p);
    for (uint64_t i = 0; i < N; i++)
        mpz_urandomm(p[i], state, m);

    gmp_randclear(state);
}

void invert_test() {
    uint32_t err = 0;
    poly p = poly_new();
    poly p_inv = poly_new();
    poly r = poly_new();
    mpz_t Q;

    mpz_init_set_str(Q, PRIME_Q, 10);
    printf("Performing %d tests...\n", N_TESTS);
    for (uint32_t i = 0; i < N_TESTS; i++) {
        do {
            poly_rand(p, Q);
        } while (!poly_invert(p_inv, p, Q));
        poly_conv(r, p, p_inv, Q);
        if (poly_degree(r) != 0 || mpz_cmp_ui(r[0], 1) != 0) {
            err++;
            printf("p := "); poly_print(p);
            printf("p_inv := "); poly_print(p_inv);
            printf("r := "); poly_print(r);
        }
    }
    if (err) {
        printf("Test Failed. Errors: %d\n", err);
    } else {
        printf("All tests have passed.\n");
    }

    poly_clear(p); poly_clear(p_inv);
    poly_clear(r); mpz_clear(Q);
}

int main() {
    invert_test();
}