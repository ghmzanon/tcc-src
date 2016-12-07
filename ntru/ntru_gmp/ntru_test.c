#include <stdio.h>
#include <string.h>
#include <immintrin.h>
#include "ntru.h"

void poly_rand_binary(poly p, uint64_t n) {
    uint16_t bit;

    poly_zero(p);
    for (uint16_t i = 0; i < n; i++) {
        _rdrand16_step((unsigned short *)&bit);
        mpz_set_ui(p[i], bit & 0x1);
    }
}

void homo_mul_P(mpz_t P, int m) {
   int n = (N - 1)/m;
   poly v = poly_new();
   poly w = poly_new();

   for (int i = 0; i <= N; i++)
      mpz_set_ui(v[i], 1);
   for (int k = 2; k <= m; k++) {
      poly_zero(w);
      for (int i = 0; i <= (k*n + 1)/2; i++) {
         for (int j = i; j >= 0; j--) {
            if (i <= n || j >= i - n) {
               mpz_add(w[i], w[i], v[j]);
            } else {
               break;
            }
         }
      }
      for (int i = (k*n + 1)/2 + 1; i <= k*n; i++) {
         mpz_set(w[i], w[k*n - i]);
      }
      poly_swap(&v, &w);
   }
   mpz_nextprime(P, v[(m*n + 1)/2]);
   poly_clear(v); poly_clear(w);
}

void ntru_homomorphic_add_test() {
    uint32_t err = 0;
    poly f = poly_new();  // private key
    poly h = poly_new();  // public key
    poly m = poly_new();
    poly m_sum = poly_new();
    poly c = poly_new();
    poly c_sum = poly_new();
    mpz_t Q, P, i;

    mpz_init_set_str(Q, PRIME_Q, 10);
    mpz_init_set_str(P, PRIME_P, 10);
    for (int n = 0 ; n < N_TESTS; n++) {
        ntru_keygen(f, h, P, Q);
        poly_rand_binary(m_sum, N);
        ntru_encrypt(c_sum, m_sum, h, Q);
        mpz_init_set_ui(i, 2);
        while (mpz_cmp(i, P) < 0) {
            poly_rand_binary(m, N);
            ntru_encrypt(c, m, h, Q);
            poly_add(c_sum, c_sum, c, Q);
            poly_add(m_sum, m_sum, m, Q);
            mpz_add_ui(i, i, 1);
        }
        ntru_decrypt(c, c_sum, f, P, Q);
        if (!poly_cmp_eq(m_sum, c)) {
            err++;
            printf("Fail.\n");
        } else {
            printf("Pass.\n");
        }
    }
    if (err) {
        printf("Test Failed. Errors: %d\n", err);
    } else {
        printf("All tests have passed.\n");
    }

    mpz_clear(i);
    mpz_clear(Q); mpz_clear(P);
    poly_clear(f); poly_clear(h);
    poly_clear(c); poly_clear(c_sum);
    poly_clear(m); poly_clear(m_sum);
}

void ntru_homomorphic_mul_test() {
    poly f = poly_new();  // private key
    poly f_mul = poly_new();  // private key
    poly h = poly_new();  // public key
    poly m = poly_new();
    poly m_mul = poly_new();
    poly c = poly_new();
    poly c_mul = poly_new();
    mpz_t Q, P;
    uint8_t flag = 1;

    mpz_init(P); mpz_init(Q);
    int M = 2;
    homo_mul_P(P, M);
    mpz_nextprime(Q, P);
    while ((N-1)/M >= 16) {
        for (int n = 0 ; n < N_TESTS; n++) {
            ntru_keygen(f, h, P, Q);
            poly_set(f_mul, f);

            poly_rand_binary(m_mul, (N-1)/M);
            ntru_encrypt(c_mul, m_mul, h, Q);
            for (int i = 1; i < M; i++) {
                poly_rand_binary(m, (N-1)/M);
                ntru_encrypt(c, m, h, Q);
                poly_conv(c_mul, c, c_mul, Q);
                poly_conv(m_mul, m, m_mul, Q);
            }
            for (int j = 1; j < M; j++) {
                poly_conv(f_mul, f, f_mul, Q);
            }
            ntru_decrypt(c_mul, c_mul, f_mul, P, Q);
            if (!poly_cmp_eq(m_mul, c_mul)) {
                flag = 0;
                mpz_mul_ui(Q, Q, 2);
                mpz_nextprime(Q, Q);
                break;
            }
        }
        if (flag) {
            gmp_printf("%Zd %Zd\n", P, Q);
            M++;
            homo_mul_P(P, M);
        }
        flag = 1;
    }

    mpz_clear(Q); mpz_clear(P);
    poly_clear(f); poly_clear(h);
    poly_clear(c); poly_clear(m);
    poly_clear(c_mul); poly_clear(m_mul); poly_clear(f_mul);
}

int main() {
    // ntru_homomorphic_add_test();
    ntru_homomorphic_mul_test();
    return 0;
}