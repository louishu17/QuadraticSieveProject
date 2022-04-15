import time
from math import sqrt, exp, log, log2
import os
from subprocess import Popen
import pickle
from itertools import count


def is_prime(n):
  if n == 2 or n == 3: return True
  if n < 2 or n%2 == 0: return False
  if n < 9: return True
  if n%3 == 0: return False
  r = int(n**0.5)
  # since all primes > 3 are of the form 6n ± 1
  # start with f=5 (which is prime)
  # and test f, f+2 for being prime
  # then loop by 6. 
  f = 5
  while f <= r:
    print('\t',f)
    if n % f == 0: return False
    if n % (f+2) == 0: return False
    f += 6
  return True   
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
def modular_sqrt(a, p):
    """ Find a quadratic residue (mod p) of 'a'. p
        must be an odd prime.

        Solve the congruence of the form:
            x^2 = a (mod p)
        And returns x. Note that p - x is also a root.

        0 is returned is no square root exists for
        these a and p.

        The Tonelli-Shanks algorithm is used (except
        for some simple cases in which the solution
        is known from an identity). This algorithm
        runs in polynomial time (unless the
        generalized Riemann hypothesis is false).
    """
    # Simple cases
    #
    if legendre(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return 0
    elif p % 4 == 3:
        return pow(a, (p + 1) / 4, p)

    # Partition p-1 to s * 2^e for an odd s (i.e.
    # reduce all the powers of 2 from p-1)
    #
    s = p - 1
    e = 0
    while s % 2 == 0:
        s /= 2
        e += 1

    # Find some 'n' with a legendre symbol n|p = -1.
    # Shouldn't take long.
    #
    n = 2
    while legendre(n, p) != -1:
        n += 1

    # Here be dragons!
    # Read the paper "Square roots from 1; 24, 51,
    # 10 to Dan Shanks" by Ezra Brown for more
    # information
    #

    # x is a guess of the square root that gets better
    # with each iteration.
    # b is the "fudge factor" - by how much we're off
    # with the guess. The invariant x^2 = ab (mod p)
    # is maintained throughout the loop.
    # g is used for successive powers of n to update
    # both a and b
    # r is the exponent - decreases with each update
    #
    x = pow(a, (s + 1) / 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

    while True:
        t = b
        m = 0
        for m in xrange(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m

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

def compute_modular_sqrts(factor_base):
    tsqrt = []
    for p in factor_base:
        tsqrt.append(int(modular_sqrt(N,p)))
    return tsqrt

def compute_log_2(factor_base):
    tlog = []
    for p in factor_base:
        tlog.append(log(p,10))
    return tlog

def mod_inv(a, m):
  a = int(a%m)
  x, u = 0, 1
  while a:
    x, u = u, x - (m/a)*u
    m, a = a, m%a
  return x
def generate_polynomials(factor_base, N, I):


    root_2N = isqrt(2*N)
    root_A = isqrt(root_2N)/I

    smooths = []
    partials = {}
    polynomials = 0
    while True:
        while True:
            root_A = nextprime(root_A)
            leg = legendre(N, root_A)
            if  leg == 1:
                break
            elif leg == 0:
                return root_A
        polynomials += 1
        a = int(root_A ** 2)
        b = modular_sqrt(N,root_A)
        b = int((b-(b*b-N) * mod_inv(2*b, root_A)) % a)
        c = int((b*b-N)/a)
        tsqrt = compute_modular_sqrts(factor_base)
        tlog = compute_log_2(factor_base)
        for i,p in enumerate(factor_base):
            ainv = pow(a,p-2,p)
            sol1 = ainv*(tsqrt[i] - b) % p
            sol2 = ainv*(-tsqrt[i] - b) % p
    return

def nextprime(n):
    if n < 2: return 2
    if n == 2: return 3
    n = (n + 1) | 1    # first odd larger than n
    m = n % 6
    if m == 3:
        if is_prime(n+2): return n+2
        n += 4
    elif m == 5:
        if is_prime(n): return n
        n += 2
    for m in count(n, 6):
        if is_prime(m  ): return m
        if is_prime(m+4): return m+4
"""
Find B-smooth numbers, using sieve and Tonelli-Shanks
"""
def find_smooth(factor_base, N, I):
    #generates a sequence from y(x) = x^2 - N starting from sqrt(N)
    def generate_sieve(N, I):
        # sieve_seq = [math.pow(x,2) - N for x in range(root-I, root + I)]
        sieve_seq = [x*x - N for x in range(root-I, root + I)]
        return sieve_seq
    t = time.process_time()
    sieve_seq = generate_sieve(N,I)
    print("time to generate seive: " + str(time.process_time() - t))
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
    
    t = time.process_time()
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
    print("time to divide by all primes in seive: " + str(time.process_time() - t))

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
        factors[-1] = 1
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
    # mprint(M)
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

def gaussian_eliminate_and_solve(matrix, xlist, smooth_nums, N):
    t = time.process_time()
    m = len(matrix)
    n = len(matrix[0])
    # print(m)
    # print(n)
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
    # print(matrix)
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

    print("gaussian elimation took: " + str(time.process_time()-t))

    for dependency in ret_all:
        factor = find_solution(dependency, xlist, smooth_nums, N)
        if factor == 1 or factor == N:
            print('try again')
        else:
            print('factor found')
            return factor, N/factor
            # return sieve_time, matrix_time, 
                        
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
    sieve_time = 0
    matrix_time = 0

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
    t = time.process_time()

    
    smooth_nums, xlist = find_smooth(factor_base, N, I)
    sieve_time = time.process_time() - t
    print("Found {} smooth numbers.".format(len(smooth_nums)))

    # print(smooth_nums)

    if len(smooth_nums) <= len(factor_base):
        return None
    
    print("Building exponent matrix...")
    #builds exponent matrix mod 2 from B-smooth numbers
    #M_transpose_matrix is either the B-smooth matrix or just one B-smooth number that is a square
    t = time.process_time()
    is_square, M_matrix, index = build_matrix(smooth_nums,factor_base)
    print("matrix loading takes " + str(time.process_time() - t))

    #case when we found a B-smooth number that is a square
    if is_square:
        factor = gcd(xlist[index] - sqrt(M_matrix), N)
        print("Found a square!")
        return sieve_time, matrix_time, factor, N/factor

    t = time.process_time()
    p_file = open('matrix.pkl', 'wb')
    pickle.dump(M_matrix, p_file)
    p_file.close()
    print("pickle dump takes " + str(time.process_time() - t))

    p_file = open('smooth.pkl', 'wb')
    pickle.dump(smooth_nums, p_file)
    p_file.close()

    p_file = open('xlist.pkl', 'wb')
    pickle.dump(xlist, p_file)
    p_file.close()

    p_file = open('N.pkl', 'wb')
    pickle.dump(N, p_file)
    p_file.close()

    # t = time.process_time()
    # p_file = open('matrix.pkl', 'rb')
    # pickle.load(p_file)
    # p_file.close()
    # print("pickle load takes " + str(time.process_time() - t))

    

    #Need to find row dependency
    t = time.process_time()
    # print("The two factors of " + str(N) + " are " + str())
    
    # row_dependencies = gaussian_eliminate(M_matrix)

    

    # os.popen('''pypy3 -c "'/Users/herbertwang/Duke 2024/Sophomore Year/Math 404/QuadraticSieveProject/QS.py' import *; row_dependencies = gaussian_eliminate(M_matrix)"''', 'w')
    cmd = '''pypy3 -c "from QS import *; import pickle; import sys; p_file = open('matrix.pkl', 'rb');\
        M_matrix = pickle.load(p_file); p_file = open('smooth.pkl', 'rb');\
        smooth_nums = pickle.load(p_file); p_file = open('xlist.pkl', 'rb');\
        xlist = pickle.load(p_file); p_file = open('N.pkl', 'rb');\
        N = pickle.load(p_file); print(gaussian_eliminate_and_solve(M_matrix, xlist, smooth_nums, N)); exit(); print(1)"'''
    p = Popen(cmd, shell=True)
    
    matrix_time = time.process_time() - t

    print("row dependencies acquired")
    
    # iterate and check all dependent rows
    
    # for dependency in row_dependencies:
    #     factor = find_solution(dependency, xlist, smooth_nums, N)
    #     if factor == 1 or factor == N:
    #         print('try again')
    #     else:
    #         print('factor found')

    #         return sieve_time, matrix_time, factor, N/factor

    # first_val = tonelli(a, N)
    
    
    # first_val = tonelli(b, N)

def calc_B_X(N,C_b,C_x):
    ln = log(N)
    # print(ln)
    B = int(exp((1/2 + C_b)*math.pow((ln*log(ln)),1/2)))
    X = int(math.pow(N, 1/2 + C_x) - isqrt(N))
    return B,X

def calc_B_X_2(N):
    ln = log(N)
    # print(ln)
    B = int(math.pow(exp(math.pow(ln*log(ln),1/2)),sqrt(2)/3))
    X = int(math.pow(exp(math.pow(ln*log(ln),1/2)),3*sqrt(2)/4))
    return B,X


if __name__ == "__main__":
    N = 16921456439215439701
    N = 46839566299936919234246726809
    # N = 6172835808641975203638304919691358469663
    # N = 1811706971
    C_b = .1
    C_x = .0000000003
    B, I = calc_B_X_2(N)
    print(B,I)
    # t = time.process_time()
    # #do some stuff
    # elapsed_time = time.process_time() - t
    # # B = 8000
    # # I = 300000
    # B = 500000
    # I = 60000000
    # B = 5000
    # I = 2500000
    # B = 400000
    # I = 1100000
    B = 5000
    I = 2500000
    B = 50000
    I = 1100000
    t = time.process_time()
    # print("The two factors of " + str(N) + " are " + str(QS(N,B,I)))
    print(str(QS(N,B,I)))
    print("sieve time, matrix time, factor, other factor: (for B = " + str(B) + " and I = " + str(I) + ")")
    elapsed_time = time.process_time() - t
    print("the elapsed_time is: " + str(elapsed_time))
    # print(calc_B_X(N, C_b, C_x))
    # time_dict = {}
    # with open('times.txt', 'w') as f:
    #     f.write("For N = " + str(N))
    #     for B in [2000, 3000, 5000, 8000, 9000]:
    #     # for B in range(5,6):
    #         for I in [100000, 300000, 5000000, 1500000, 25000000, 75000000]:
    #         # for I in [25000000]:
                
    #             total_time = 0
    #             sieve_total_time = 0
    #             gaussian_total_time = 0
    #             trials = 2
    #             time_array = []
    #             sieve_time_array = []
    #             gaussian_time_array = []
    #             total_trials = 0
    #             for _ in range(trials):
    #                 t = time.process_time()
    #                 
                    # composed = QS(N,B,I)
                    # elapsed_time = time.process_time() - t
    #                 if composed:
    #                     sieve_time, gaussian_time, factor1, factor2 = composed
    #                     sieve_total_time += sieve_time
    #                     gaussian_total_time += gaussian_time
    #                     print()
                        
    #                     total_time += elapsed_time
    #                     time_array.append(elapsed_time)
    #                     sieve_time_array.append(sieve_time)
    #                     gaussian_time_array.append(gaussian_time)
    #                     total_trials += 1
    #             if not B in time_dict:
    #                 time_dict[B] = {}

    #             if not I in time_dict[B]:
    #                 time_dict[B][I] = {}
    #             if not total_trials == 0:
    #                 time_dict[B][I]['avg'] = total_time/total_trials
    #                 time_dict[B][I]['sieve_avg'] = sieve_total_time/total_trials
    #                 time_dict[B][I]['gaussian_avg'] = gaussian_total_time/total_trials
    #             else:
    #                 time_dict[B][I]['avg'] = 0
    #                 time_dict[B][I]['sieve_avg'] = 0
    #                 time_dict[B][I]['gaussian_avg'] = 0

    #             time_dict[B][I]['times'] = time_array
                
    #             time_dict[B][I]['sieve_times'] = sieve_time_array
                
    #             time_dict[B][I]['gaussian_times'] = gaussian_time_array

    #             f.write('The average time for B = ' + str(B) + ' and I = ' + str(I) \
    #                 + " is " + str(time_dict[B][I]['avg']) + ". (values = " + str(time_array) + ")")
    #             f.write('The average time for B = ' + str(B) + ' and I = ' + str(I) \
    #                 + " is " + str(time_dict[B][I]['sieve_avg']) + ". (values = " + str(sieve_time_array) + ")")
    #             f.write('The average time for B = ' + str(B) + ' and I = ' + str(I) \
    #                 + " is " + str(time_dict[B][I]['gaussian_avg']) + ". (values = " + str(gaussian_time_array) + ")")
    
    # with open('times.txt', 'w') as f:
    #     f.write("For N = " + str(N))
    #     for B, B_values in time_dict.items():
    #         for I, I_values in B_values.items():
                
    #             f.write('The average time for B = ' + str(B) + ' and I = ' + str(I) \
    #                 + " is " + str(I_values['avg']) + ". (values = " + str(I_values['times']) + ") ")
    #             f.write('The average sieve time for B = ' + str(B) + ' and I = ' + str(I) \
    #             + " is " + str(I_values['sieve_avg']) + ". (values = " + str(I_values['sieve_times']) + ") ")
    #             f.write('The average gaussian time for B = ' + str(B) + ' and I = ' + str(I) \
    #             + " is " + str(I_values['gaussian_avg']) + ". (values = " + str(I_values['gaussian_times']) + ") ")

            
