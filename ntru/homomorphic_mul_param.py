import math

N = 64 # lattice dimension
M = 4 # circuit depth

# -------------------------------------------------------------------- #
n = (N-1)//M

v = (n+1)*[1]
for k in range(2, M + 1):
    w = (k*n+1)*[0]
    for i in range(0, (k*n + 1)//2 + 1):
        for j in range(i, -1, -1):
            if (i <= n):
                w[i] += v[j]
            elif (j >= i - n):
                w[i] += v[j]
            else:
                break
    for i in range((k*n + 1)//2 + 1, k*n + 1):
        w[i] = w[k*n-i]
    v = w
print("p > NextPrime(" + str(v[(M*n + 1)//2]) + ")")
