#!/usr/bin/python

from random import randrange
import consts

def new_polynomial():
    return consts.N*[0]

def clear_polynomial(p):
    for i in range(0, consts.N):
        p[i] = 0

def print_polynomial(p):
    print(p[0], end="")
    for i in range(1, consts.N):
        if (p[i] != 0):
            print(" + " + str(p[i]) + "*x^" + str(i), end="")
    print("")

def find_degree(p):
    for i in range (consts.N-1, -1, -1):
        if (p[i] != 0):
            return i
    return 0

def modularize_polynomial(p, m):
    for i in range(consts.N):
        p[i] %= m

def modular_inverse(a, m):
    t = 0; newt = 1
    r = m; newr = a
    while (newr != 0):
        q = r//newr
        (t, newt) = (newt, t - q*newt)
        (r, newr) = (newr, r - q*newr)
    if (r > 1): return -1
    if (t < 0): t += m
    return t

def polynomial_inverse(p, m):
    u = p[:]
    modularize_polynomial(u, m)
    du = find_degree(u)
    v = (consts.N+1)*[0]; v[consts.N] = 1; v[0] = m - 1; dv = consts.N
    t = new_polynomial(); t[0] = 1; dt = 0
    s = new_polynomial(); ds = 0

    while (du > 0):
        if (du < dv):
            (u, v) = (v, u);  (du, dv) = (dv, du)
            (s, t) = (t, s);  (ds, dt) = (dt, ds)

        j = du - dv

        mod_inv = modular_inverse(v[dv], m)
        if (mod_inv != -1):
            c = (u[du]*mod_inv) % m
        else: # if m is a prime number, v[dv] should always have modular inverse
            return [-1]

        for i in range(dv+1):
            u[i+j] -= c*v[i]
            u[i+j] %= m
        du = find_degree(u)

        for i in range(ds+1):
            t[i+j] -= c*s[i]
            t[i+j] %= m
        dt = find_degree(t)

    mod_inv = modular_inverse(u[0], m)
    if (du == 0 and mod_inv != -1):
        for i in range(dt+1):
            t[i] = (t[i]*mod_inv) % m
        return t

    return [-1]

def convolution(p, q, m):
    t = new_polynomial()

    for i in range(consts.N):
        for j in range(consts.N):
            t[i] += p[j]*q[(consts.N+i-j)%consts.N]
            t[i] %= m
    return t

def add_polynomials(p, q, m):
    t = new_polynomial()
    for i in range(consts.N):
        t[i] += p[i] + q[i]
        t[i] = t[i] % m
    return t

def rand_polynomial(m):
    p = new_polynomial()
    for i in range(consts.N):
        p[i] = random.randrange(-m+1, m)
    return p;

###############################################################################
# Tests

def inverse_modulo_test(): #OK
    err = 0
    for i in range(1000):
        a = random.randrange(1, consts.Q)
        a_inv = modular_inverse(a, consts.Q)
        if ((a*a_inv) % consts.Q != 1):
            print("a:     " + str(a))
            print("a_inv: " + str(a_inv))
            err += 1
    print(err)

def polynomial_inverse_test(): #OK
    err = 0
    for i in range(10):
        q = [-1]
        while (len(q) == 1 and q[0] == -1):
            p = rand_polynomial(consts.Q)
            q = polynomial_inverse(p, consts.Q)
        r = convolution(p, q, consts.Q)
        modularize_polynomial(r, consts.Q)
        if (find_degree(r) != 0 or r[0] != 1):
            err += 1
    print(err)

# def main():
    # inverse_modulo_test()
    # polynomial_inverse_test()

# if __name__ == "__main__":
    # main()

