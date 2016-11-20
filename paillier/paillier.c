#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include <immintrin.h>

#define VERBOSE
#define PAILLIER_TEST

#ifdef PAILLIER_TEST
#define N_TESTS 1000000
#ifdef VERBOSE
#include <stdio.h>
#endif
#endif

#define ceil16(x) { (x >> 4) + ((x & 0x7) == 0 ? 0 : 1) }
#define PRIME_BITS 512
#define REPS 50

/* p represents the prime, k it's bit size and t is the number of
 * repetitions that miller-rabin primality test makes */
void mpz_random_prime(mpz_t p, uint32_t k, size_t t) {
    size_t size = ceil16(k);
    uint16_t *short_stream = malloc(size * sizeof(uint16_t));
    do {
        for (int i = 0; i < size; i++)
            _rdrand16_step((uint16_t *)&short_stream[i]);
        /* Makes sure that the integer is odd */
        short_stream[size - 1] |= 0x0100;
        mpz_import(p, size, 1, sizeof(uint16_t), 1, 0, short_stream);
    } while (mpz_probab_prime_p(p, t) == 0);
    free(short_stream);
}

void func_L(mpz_t *u, const mpz_t n) {
    mpz_sub_ui(*u, *u, 1); // u = u-1
    mpz_fdiv_q(*u, *u, n); // u = floor(u/n)
}

void paillier_keygen(mpz_t *n, mpz_t *g, mpz_t *lamb, mpz_t *me) {
    mpz_t p, q, n2;
    mpz_init(p); mpz_init(q); mpz_init(n2);

    gmp_randstate_t state;
    gmp_randinit_mt(state);

    uint64_t seed;
    _rdrand64_step((unsigned long long int *)&seed);
    gmp_randseed_ui(state, seed);

    do {
        mpz_random_prime(p, PRIME_BITS, REPS);
        mpz_random_prime(q, PRIME_BITS, REPS);
        mpz_mul(*n, p ,q);          // n = pq
        mpz_sub_ui(*g, p, 1);       // g = p-1
        mpz_sub_ui(*me, q, 1);      // me = q-1
        mpz_mul(*lamb, *g, *me);    // lamb = (p-1)(q-1)
        mpz_gcd(*lamb, *n, *lamb);  // lamb = gcd(p*q, (p-1)(q-1))
    } while(mpz_cmp_ui(*lamb, 1));

    mpz_lcm(*lamb, *g, *me);
    mpz_mul(n2, *n, *n);

    do {
        do {
            mpz_urandomm(*g, state, n2); // g \in Z_n^2
        } while(!mpz_sgn(*g));           // g \in Z*_n^2
        mpz_powm(*me, *g, *lamb, n2);
        func_L(me, *n);
    } while(!mpz_invert(*me, *me, *n));

    mpz_clear(p); mpz_clear(q); mpz_clear(n2);
    gmp_randclear(state);
}


void paillier_encrypt(mpz_t *c, const mpz_t m, const mpz_t n, const mpz_t g) {
    mpz_t r, n2;
    mpz_init(r); mpz_init(n2);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    uint64_t seed;
    _rdrand64_step((unsigned long long int *)&seed);
    gmp_randseed_ui(state, seed);

    mpz_urandomm(r, state, n); // r \in Z_n
    mpz_mul(n2, n, n);
    mpz_powm(*c, g, m, n2);
    mpz_powm(r, r, n, n2);
    mpz_mul(*c, *c, r);
    mpz_mod(*c, *c, n2);

    mpz_clear(r); mpz_clear(n2);
    gmp_randclear(state);
}

void paillier_decrypt(mpz_t *m, const mpz_t c, const mpz_t n, const mpz_t lamb,
        const mpz_t me) {
    mpz_t n2;
    mpz_init(n2);

    mpz_mul(n2, n, n);
    mpz_powm(*m, c, lamb, n2);
    func_L(m, n);
    mpz_mul(*m, *m, me);
    mpz_mod(*m, *m, n);

    mpz_clear(n2);
}

/**********************************************************************/
#ifdef PAILLIER_TEST
void simple_test() {
    mpz_t n, g, lamb, me, c, m, m_d;
    mpz_init(n); mpz_init(g);     // public key
    mpz_init(lamb); mpz_init(me); // private key
    mpz_init(c); mpz_init(m); mpz_init(m_d);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    uint64_t seed;
    _rdrand64_step((unsigned long long int *)&seed);
    gmp_randseed_ui(state, seed);

    uint16_t err = 0;

    for (uint16_t i = 0; i < N_TESTS; i++) {
        paillier_keygen(&n, &g, &lamb, &me);

        mpz_urandomm(m, state, n); // m \in Z_n

        paillier_encrypt(&c, m, n, g);
        paillier_decrypt(&m_d, c, n, lamb, me);

        if (mpz_cmp(m, m_d)) {
#ifdef VERBOSE
            gmp_printf("n: %Zd\n", n);
            gmp_printf("g: %Zd\n", g);
            gmp_printf("lamb: %Zd\n", lamb);
            gmp_printf("me: %Zd\n", me);
            printf("\n");
            gmp_printf("m: %Zd\n", m);
            gmp_printf("m_d: %Zd\n", m_d);
            printf("----\n");
#endif
            err++;
        }
    }

    printf("ERRORS: %d of %d\n", err, N_TESTS);

    mpz_clear(n); mpz_clear(g);
    mpz_clear(lamb); mpz_clear(me);
    mpz_clear(c); mpz_clear(m); mpz_clear(m_d);
    gmp_randclear(state);
}

void add_homomorphic_test() {
    mpz_t n, g, lamb, me, c1, c2, m1, m2, n2;
    mpz_init(n); mpz_init(g);     // public key
    mpz_init(lamb); mpz_init(me); // private key
    mpz_init(m1); mpz_init(m2);
    mpz_init(c1); mpz_init(c2);
    mpz_init(n2);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    uint64_t seed;
    _rdrand64_step((unsigned long long int *)&seed);
    gmp_randseed_ui(state, seed);

    uint16_t err = 0;

    for (uint16_t i = 0; i < N_TESTS; i++) {
        paillier_keygen(&n, &g, &lamb, &me);

        mpz_urandomm(m1, state, n); // m \in Z_n
        mpz_urandomm(m2, state, n); // m \in Z_n

        mpz_mul(n2, n, n);

        paillier_encrypt(&c1, m1, n, g);
        paillier_encrypt(&c2, m2, n, g);

        mpz_mul(c1, c1, c2); // c = c1c2
        mpz_mod(c1, c1, n2); // c = c1c2 mod n^2 = d(m1) + d(m2) mod n

        paillier_decrypt(&c1, c1, n, lamb, me);
        mpz_add(m1, m1, m2); // m = m1 + m2
        mpz_mod(m1, m1, n);

        if (mpz_cmp(m1, c1)) {
#ifdef VERBOSE
            gmp_printf("n: %Zd\n", n);
            gmp_printf("g: %Zd\n", g);
            gmp_printf("lamb: %Zd\n", lamb);
            gmp_printf("me: %Zd\n", me);
            printf("\n");
            gmp_printf("m: %Zd\n", m1);
            gmp_printf("c: %Zd\n", c1);
            printf("----\n");
#endif
            err++;
        }
    }

    printf("ERRORS: %d of %d\n", err, N_TESTS);

    mpz_clear(n); mpz_clear(g);
    mpz_clear(lamb); mpz_clear(me);
    mpz_clear(m1); mpz_clear(m2);
    mpz_clear(c1); mpz_clear(c2);
    mpz_clear(n2);
    gmp_randclear(state);
}

int main() {
    /* simple_test(); */
    add_homomorphic_test();
    return 0;
}
#endif
