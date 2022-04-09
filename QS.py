from math import sqrt


def sieveOfEratosthenes(n):
    isPrime = [True for _ in range(n+1)]
    isPrime[0] = isPrime[1] = False
    
    ret = []
    for i, isprime in enumerate(isPrime):
        if isprime:
            ret.append(i)
            for j in range(i*i, n+1, i):
                isPrime[j] = False
    return ret

def legendre(a,p):
    return pow(a, (p-1) // 2, p)

#generates B-smooth factor base
def find_base(N, B): 
    factor_base = []
    primes = sieveOfEratosthenes(B)
    print("Primes under " + str(B) + ": " + str(primes))

    for p in primes:
        if legendre(N, p) == 1:
            factor_base.append(p)
    return factor_base


def tonelli(n, p): #tonelli-shanks to solve modular square root, x^2 = N (mod p)
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r,p-r
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return (r,p-r)

"""
Find B-smooth numbers, using sieve and Tonelli-Shanks
"""
def find_smooth(factor_base, N, I):
    #generates a sequence from y(x) = x^2 - N starting from sqrt(N)
    def generate_sieve(N, I):
        sieve_seq = [x**2 - N for x in range(root-I, root + I)]
        return sieve_seq
    
    sieve_seq = generate_sieve(N,I)
    sieve_list = sieve_seq.copy()
    
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        
        #now we are at a multiple of 2, so every other term will be divisible by 2
        for j in range(i, len(sieve_list), 2):
            #get rid of all powers of 2
            while sieve_list[j] % 2 == 0:
                sieve_list[j] //= 2
    
    for p in factor_base[1:]:
        residues = tonelli(N, p)

        for r in residues:
            for i in range((r-root+I) % p, len(sieve_list), p):
                #get rid of all powers of p
                while sieve_list[i] % p == 0:
                    sieve_list[i] //= p
            #negative direciton
            for i in range(((r-root+ I) % p) + I, 0, -p):
                while sieve_list[i] % p == 0:
                    sieve_list[i] //= p

    B_smooth_nums = []

    for i in range(len(sieve_list)):
        #we have enough rows to achieve a linear dependence
        if len(B_smooth_nums) >= len(factor_base) + 1:
            break
        #found a b-smooth number    
        if sieve_list[i] == 1 or sieve_list == -1:
            B_smooth_nums.append(sieve_seq[i])
    
    return B_smooth_nums

# """
# Build exponent vectors mod 2 from B-smooth numbers, then combines into a matrix
# """
# def build_matrix(smooth_nums, factor_base):
#     M = []
#     factor_base.insert(0,-1)

#     for n in smooth_nums:

#main function
def QS(n, B, I):

    global N
    global root

    N, root= n, int(sqrt(n))
    print(N, root)

    print("Attempting to factor {}...".format(N))


    print("Generating {}-smooth factor base...".format(B))
    factor_base = find_base(N,B)

    global F
    F = len(factor_base)
    print("Factor Base: " + str(factor_base))

    print("Looking for {} {}-smooth numbers...".format(F+1, B))
    #find B-smooth numbers, using sieve and Tonelli-Shanks
    smooth_nums = find_smooth(factor_base, N, I)

    print("Found {} smooth numbers.".format(len(smooth_nums)))

    print(smooth_nums)

if __name__ == "__main__":
    QS(8051, 300, 2500)
