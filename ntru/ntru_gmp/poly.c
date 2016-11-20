#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>
#include "poly.h"

/* #define POLY_TEST */

#ifdef POLY_TEST
#include <time.h>
#define N_TESTS 1000
#endif

poly poly_new() {
    poly p = malloc((N + 1) * sizeof(mpz_t));
    if (p == NULL) perror(NULL);
    for (uint16_t i = 0; i <= N; i++)
        mpz_init(p[i]);
    return p;
}

void poly_zero(poly p) {
    for (uint16_t i = 0; i <= N; i++)
        mpz_set_ui(p[i], 0);
}

void poly_clear(poly p) {
    for (uint16_t i = 0; i <= N; i++)
        mpz_clear(p[i]);
    free(p);
}

void poly_print(poly p) {
    gmp_printf("%Zd", p[0]);
    for (uint16_t i = 1; i <= N; i++)
        if (mpz_sgn(p[i]))
            gmp_printf(" + %Zd*x^%d", p[i], i);
    printf("\n");
}

#ifdef POLY_TEST
void poly_rand(poly p, mpz_t m) {
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    /* uint64_t seed; */
    /* _rdrand64_step((unsigned long long int *)&seed); */
    /* gmp_randseed_ui(state, seed); */
    srand(time(NULL));
    gmp_randseed_ui(state, rand());

    poly_zero(p);
    for (uint64_t i = 0; i < N; i++)
        mpz_urandomm(p[i], state, m);

    gmp_randclear(state);
}
#endif

uint16_t poly_degree(poly p) {
    for (uint16_t i = N; i > 0; i--)
        if (mpz_sgn(p[i]))
            return i;
    return 0;
}

void poly_swap(poly *p, poly *q) {
    poly ptr;

    ptr = *p;
    *p = *q;
    *q = ptr;
}

static void swap(uint16_t *a, uint16_t *b) {
    uint16_t tmp = *a;
    *a = *b;
    *b = tmp;
}

void poly_mod(poly p, mpz_t m) {
    for (uint16_t i = 0; i < N; i++)
        mpz_mod(p[i], p[i], m);
}

void poly_conv(poly r, poly p, poly q, mpz_t m) {
    poly t = poly_new();
    mpz_t tmp;
    mpz_init(tmp);

    for (uint16_t i = 0; i < N; i++)
        for (uint16_t j = 0; j < N; j++) {
            mpz_mul(tmp, p[j], q[(N + i - j) % N]);
            mpz_add(t[i], t[i], tmp);
            /* if (mpz_cmp(t[i], m) >= 0 || mpz_sgn(t[i]) < 0) */
            mpz_mod(t[i], t[i], m);
        }

    if (r != NULL)
        poly_set(r, t);

    mpz_clear(tmp);
    poly_clear(t);
}

void poly_add(poly r, poly p, poly q, mpz_t m) {
    for (uint16_t i = 0; i < N; i++) {
        mpz_add(r[i], p[i], q[i]);
        mpz_mod(r[i], r[i], m);
    }
}

void poly_set(poly p, poly q) {
    for (uint16_t i = 0; i <= N; i++)
        mpz_set(p[i], q[i]);
}

bool poly_cmp_eq(poly p, poly q) {
    for (uint16_t i = 0; i <= N; i++)
        if (mpz_cmp(p[i], q[i]))
            return 0;
    return 1;
}

void poly_center(poly p, mpz_t m) {
    mpz_t u, l;
    mpz_init(u);
    mpz_init(l);
    mpz_fdiv_q_2exp(u, m, 1); // u = m/2^1
    mpz_set(l, u);
    mpz_neg(l, u); // l = -u

    for (int64_t i = 0; i <= N; i++) {
        if (mpz_cmp(p[i], u) > 0)
            mpz_sub(p[i], p[i], m);
        else if (mpz_cmp(p[i], l) < 0)
            mpz_add(p[i], p[i], m);
    }

    mpz_clear(u);
    mpz_clear(l);
}

bool poly_invert(poly r, poly p, mpz_t m) {
    poly u = poly_new();
    poly_set(u, p);
    uint16_t du = poly_degree(u);

    poly v = poly_new();
    mpz_set_ui(v[N], 1);
    mpz_sub_ui(v[0], m, 1);
    uint16_t dv = N;

    poly s = poly_new();
    uint16_t ds = 0;

    poly t = poly_new();
    mpz_set_ui(t[0], 1);
    uint16_t dt = 0;

    mpz_t tmp, c;
    mpz_init(tmp);
    mpz_init(c);

    while (du > 0) {
        if (du < dv) {
            poly_swap(&u, &v);
            swap(&du, &dv);
            poly_swap(&s, &t);
            swap(&ds, &dt);
        }

        if (mpz_invert(tmp, v[dv], m)) {
            mpz_mul(c, u[du], tmp);
            mpz_mod(c, c, m);
        } else {
            mpz_clear(tmp);
            mpz_clear(c);
            poly_clear(u);
            poly_clear(v);
            poly_clear(s);
            poly_clear(t);
            return false;
        }

        uint16_t j = du - dv;

        for (uint16_t i = 0; i <= dv; i++) {
            mpz_mul(tmp , c, v[i]);
            mpz_sub(u[i + j], u[i + j], tmp);
            mpz_mod(u[i + j], u[i + j], m);
        }
        du = poly_degree(u);

        for (uint16_t i = 0; i <= ds; i++) {
            mpz_mul(tmp , c, s[i]);
            mpz_sub(t[i + j], t[i + j], tmp);
            mpz_mod(t[i + j], t[i + j], m);
        }
        dt = poly_degree(t);
    }

    if (du == 0 && mpz_invert(tmp, u[0], m)) {
        for (uint16_t i = 0; i <= dt; i++) {
            mpz_mul(t[i], t[i], tmp);
            mpz_mod(t[i], t[i], m);
        }
        if (r != NULL)
            poly_set(r, t);
        mpz_clear(tmp); mpz_clear(c);
        poly_clear(u); poly_clear(v);
        poly_clear(s); poly_clear(t);
        return true;
    }

    mpz_clear(tmp); mpz_clear(c);
    poly_clear(u); poly_clear(v);
    poly_clear(s); poly_clear(t);
    return false;
}

#ifdef POLY_TEST
int main() {
    uint32_t err = 0;
    poly p = poly_new();
    poly p_inv = poly_new();
    poly r = poly_new();

    mpz_t Q;
    mpz_init_set_str(Q, "11", 10);

    for (uint32_t i = 0; i < N_TESTS; i++) {
        do {
            poly_rand(p, Q);
        } while (!poly_invert(p_inv, p, Q));
        poly_conv(r, p, p_inv, Q);
        if (poly_degree(r) != 0 || mpz_cmp_ui(r[0], 1) != 0) {
            err++;
            printf("p := ");
            poly_print(p);
            printf("p_inv := ");
            poly_print(p_inv);
            printf("r := ");
            poly_print(r);
        }
    }
    poly_clear(p);
    poly_clear(p_inv);
    poly_clear(r);
    mpz_clear(Q);
    printf("ERRORS: %d\n", err);
}
#endif
