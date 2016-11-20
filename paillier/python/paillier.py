#!/usr/bin/python
from random import randrange

def gcd(a, b):
    while a != 0:
        (a, b) = (b % a, a)
    return b

def lcm(a, b):
    return a * (b//gcd(a, b))

def mod_inv(a, m):
    t = 0; newt = 1
    r = m; newr = a
    while (newr != 0):
        q = r//newr
        (t, newt) = (newt, t - q*newt)
        (r, newr) = (newr, r - q*newr)
    if (r > 1): return None
    if (t < 0): t += m
    return t

def L(u, n):
    return (u - 1)//n

def keygen():
    while (True):
        p = 115792089237316195423570985008687907853269984665640564039457584007913129640233
        q = 115792089237316195423570985008687907853269984665640564039457584007913129640237
        n = p*q
        if (gcd(n, (p - 1)*(q - 1)) == 1):
            break

    lamb = lcm(p - 1, q - 1)

    while (True):
        g = randrange(1, n*n)
        me = mod_inv(L(pow(g, lamb, n*n), n), n)
        if (me != None):
            break

    return [[n, g], [lamb, me]] # public and private keys, respectively

def encrypt(m, pub):
    n = pub[0]
    g = pub[1]
    if (m < 0 or m >= n):
        print("Message" +  m + " is not in range(0, "+  n + ")")
        return
    r = randrange(1, n)
    c = (pow(g, m, n*n)*pow(r, n, n*n)) % (n*n)
    return c

def decrypt(c, pub, prv):
    n = pub[0]
    lamb = prv[0]
    me = prv[1]
    m = (L(pow(c, lamb, n*n), n)*me) % n
    return m

###############################################################################
# TESTS

def simple_test():
    err = 0
    for i in range(1000):
        keys = keygen()
        m = randrange(0, keys[0][0])
        c = encrypt(m, keys[0])
        if (c != None):
            m_d = decrypt(c, keys[0], keys[1])
        if (m != m_d):
            err += 1
    print(err)

def homomorphic_test():
    err = 0
    for i in range(100):
        keys = keygen()
        n = keys[0][0]

        m1 = randrange(0, n)
        m2 = randrange(0, n)
        m = (m1 + m2) % n

        c1 = encrypt(m1, keys[0])
        c2 = encrypt(m2, keys[0])
        c = (c1*c2) % (n*n) # c1 * c2 = (enc(m1) + enc(m2)) mod n

        if (c == None):
            print("ERROR")
            return

        sum_m = decrypt(c, keys[0], keys[1])
        if (m != sum_m):
            err += 1
    print(err)

def main():
    simple_test()
    # homomorphic_test()

if __name__ == "__main__":
    main()

