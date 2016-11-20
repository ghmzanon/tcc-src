#!/usr/bin/python

from polynomial import (new_polynomial, polynomial_inverse, convolution,
        add_polynomials, modularize_polynomial, print_polynomial)
from random import randrange
import consts

def rand_ternary_polynomial(val):
    p = new_polynomial()
    for i in range(0, 2*consts.D + 1):
        while (True):
            j = randrange(0, consts.N)
            if (p[j] == 0):
                break
        if (i <= consts.D):
            p[j] = 1*val
        else:
            p[j] = -1*val
    return p

def rand_binary_polynomial():
    p = new_polynomial()
    for i in range(consts.N):
        j = randrange(0, consts.N)
        if (p[j] == 0):
            p[j] = 1
    return p

def center(p):
    for i in range(consts.N):
        if (p[i] > consts.Q//2):
            p[i] -= consts.Q
        elif (p[i] < -consts.Q//2):
            p[i] += consts.Q

def keygen():
    while (True):
        f = rand_ternary_polynomial(consts.P) # private key
        f[0] += 1
        f_q = polynomial_inverse(f, consts.Q)
        if (len(f_q) != 1):
            break
    while (True):
        g = rand_ternary_polynomial(consts.P)
        if (len(polynomial_inverse(g, consts.Q)) != 1):
            break
    h = convolution(f_q, g, consts.Q) # public key
    return [h, f]

def encrypt(m, h):
    while (True):
        r = rand_ternary_polynomial(1)
        if (len(polynomial_inverse(r, consts.Q)) != 1):
            break
    # c = h*r + m (mod Q)
    c = add_polynomials(convolution(h, r, consts.Q), m, consts.Q)
    return c

def decrypt(c, f):
    m = convolution(f, c, consts.Q) # m = f*c (mod Q)
    center(m)
    modularize_polynomial(m, consts.P) # m = m (mod P)
    return m

###############################################################################

def main():
    m1 = rand_binary_polynomial()
    m2 = rand_binary_polynomial()
    # m3 = rand_binary_polynomial()
    # m4 = rand_binary_polynomial()

    # m1 = new_polynomial()
    # m2 = new_polynomial()

    # m1[4] = m1[3] = 1 # m1 = (11000) = (24)
    # m2[6] = m2[3] = m2[2] = 1 # m2 = (1001100) = (76)

    print("Generating keys..")
    keys = keygen()
    print("Generated!")

    print("Encrypting message 1...")
    c1 = encrypt(m1, keys[0])
    print("Encrypted!")
    print("Encrypting message 2...")
    c2 = encrypt(m2, keys[0])
    print("Encrypted!")
    # print("Encrypting message 3...")
    # c3 = encrypt(m3, keys[0])
    # print("Encrypted!")
    # print("Encrypting message 4...")
    # c4 = encrypt(m4, keys[0])
    # print("Encrypted!")
    
    print("Performing homomorphic addition...")
    c = add_polynomials(c1, c2, consts.Q)
    # c = add_polynomials(c, c3, consts.Q)
    # c = add_polynomials(c, c4, consts.Q)
    print("Added!")
    
    print("Decrypting homomorphic added message...")
    m = decrypt(c, keys[1])
    print("Decrypted!")
    
    print("Adding plain messages...")
    sum_m = add_polynomials(m1, m2, consts.Q)
    # sum_m = add_polynomials(sum_m, m3, consts.Q)
    # sum_m = add_polynomials(sum_m, m4, consts.Q)
    print("Added!")

    if (m == sum_m):
        print("OK")
    else:
        print("ERROR")
        print_polynomial(sum_m)
        print_polynomial(m)
        for i in range(consts.N):
            if (sum_m[i] != m[i]):
                print(i, sum_m[i], m[i])

if __name__ == "__main__":
    main()