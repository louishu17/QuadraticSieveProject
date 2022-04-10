import numpy as np
import math
from math import sqrt, exp, log
from gaussian_eliminate import gaussian_elimiate


def gcd(a,b): # Euclid's algorithm
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)
    
def isqrt(n): # Newton's method, returns exact int for large squares
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

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

def mprint(M): #prints a matrix in readable form
    for row in M:
        print(row)

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
    xlist = []

    for i in range(len(sieve_list)):
        #we have enough rows to achieve a linear dependence
        if len(B_smooth_nums) >= len(factor_base) + 1:
            break
        #found a b-smooth number    
        if sieve_list[i] == 1 or sieve_list == -1:
            B_smooth_nums.append(sieve_seq[i])
            xlist.append(i+root-I)
    
    return B_smooth_nums, xlist

"""
Finds how many times each factor in factor_base goes into n, returns list of all factors
"""
def factor(n, factor_base):
    factors = {}

    if n < 0:
        factors.append(-1)
    for p in factor_base:
        if p == -1:
            continue
        else:
            while n % p == 0:
                if p not in factors:
                    factors[p] = 1
                else:
                    factors[p] += 1
                n //= p
    
    return factors
"""
Build exponent vectors mod 2 from B-smooth numbers, then combines into a matrix
"""
def build_matrix(smooth_nums, factor_base):
    M = []
    factor_base.insert(0,-1)

    for index,n in enumerate(smooth_nums):
        exp_vector = [0] * len(factor_base)
        n_factors = factor(n, factor_base)
        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = n_factors[factor_base[i]] % 2
        
        #we found a square number already
        if 1 not in exp_vector:
            return True, n, index
        
        M.append(exp_vector)
    
    print("Matrix built:")
    mprint(M)
    # return False, transpose(M), -1
    return False, M, -1


def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)

def gaussian_elimiate(matrix):
    m = len(matrix)
    n = len(matrix[0])
    print(m)
    print(n)
    # if there are less rows than columns, then there wouldn't be a linear dependence
    if m < n:
        print("error")
        raise Exception("not enough data") 
    # set a m length array to locate the rows with pivots
    pivot = [0]*m
   
    pivot_dict = {}
    
    for j in range(n):

        # looks for the pivot in the column
        for i in range(m):
            # if a 1 is found at the i,j value
            # print(i, j, matrix[i][j])
            if(matrix[i][j] == 1):
                pivot[i] = 1
                # records the corresponding value for j where the pivot occurs
                pivot_dict[j]=i
                # adds the two rows using mod 2 addition
                for k in range(0,j):
                    if (matrix[i][k] == 1):
                        for row in range(m):
                            matrix[row][k] = (matrix[row][j] + matrix[row][k])%2

                for k in range(j+1,n):
                    if (matrix[i][k] == 1):
                        for row in range(m):
                            matrix[row][k] = (matrix[row][j] + matrix[row][k])%2
                break
    print(matrix)
    # stores which rows are dependent
    
    ret_all = []
    for i in range(m):
        ret = []
        # if there is a row that isnt a pivot, then it finds all the 1's in the row
        #  and the corresponding rows to those ones
        if pivot[i] == 0:
            ret.append(i)
            for key, col_value in enumerate(matrix[i]):
                if col_value == 1:
                    ret.append(pivot_dict[key])
            ret_all.append(ret)
                        
    return ret_all

def find_solution(dependent_rows, xlist, smooth_nums, N):
    A = 1
    b = 1
    for i in dependent_rows:
        A *= smooth_nums[i]
        b *= xlist[i]

    a = isqrt(A)
    
    factor = gcd((b-a)%N,N)
    return factor

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
    smooth_nums, xlist = find_smooth(factor_base, N, I)

    print("Found {} smooth numbers.".format(len(smooth_nums)))

    print(smooth_nums)

    if len(smooth_nums) < len(factor_base):
        return("Not enough smooth numbers. Increase the sieve interval or size of the factor base.")
    
    print("Building exponent matrix...")
    #builds exponent matrix mod 2 from B-smooth numbers
    #M_transpose_matrix is either the B-smooth matrix or just one B-smooth number that is a square
    is_square, M_transpose_matrix, index = build_matrix(smooth_nums,factor_base)

    #case when we found a B-smooth number that is a square
    if is_square:
        factor = gcd(xlist[index] - sqrt(M_transpose_matrix), N)
        print("Found a square!")
        return factor, N/factor

    #Need to find row dependency
    row_dependencies = gaussian_elimiate(M_transpose_matrix)
    print(row_dependencies)

    # iterate and check all dependent rows
    
    # for dependency in row_dependencies:
    #     a = 1
    #     b = 1
    #     for idex in dependency:
    #         a = (a * smooth_nums[idex])
    #         b = (b * xlist[idex]) % N
    #     a_rt = sqrt(a) % N
    #     print(a_rt)
    #     print((a_rt-b)%N)
    #     gcd_ab = gcd((a_rt-b)%N, N)
    #     print(gcd_ab)
    #     if not gcd_ab == 1 or gcd_ab == N:
    #         return gcd_ab, N/gcd_ab
    for dependency in row_dependencies:
        factor = find_solution(dependency, xlist, smooth_nums, N)
        if factor == 1 or factor == N:
            print('try again')
        else:
            print('factor found')
            return factor, N/factor

    # first_val = tonelli(a, N)
    
    
    # first_val = tonelli(b, N)

def calc_B_X(N,C_b,C_x):
    ln = log(N)
    B = int(exp((1/2 + C_b)*math.pow((ln*log(ln)),1/2)))
    X = int(math.pow(N, 1/2 + C_x) - isqrt(N))
    return B,X


if __name__ == "__main__":
    N = 16921456439215439701
    N = 46839566299936919234246726809
    # N = 1811706971
    C_b = .125
    C_x = .00000003
    B, I = calc_B_X(N, C_b, C_x)
    # print(B,I)
    # B = 8000
    # I = 300000
    B = 4000
    I = 25000000
    print("The two factors of " + str(N) + " are " + str(QS(N,B,I)))
    print(calc_B_X(N, C_b, C_x))
