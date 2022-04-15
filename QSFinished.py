from math import fabs, ceil, sqrt, exp, log
from itertools import chain
import time

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

def build_matrix(smooth_nums,factor_base):
# generates exponent vectors mod 2 from previously obtained smooth numbers, then builds matrix

    def factor(n,factor_base):#trial division from factor base
        factors = []
        if n < 0:
            factors.append(-1)
        for p in factor_base:
            if p == -1:
                pass
            else:
                while n % p == 0:
                    factors.append(p)
                    n //= p
        return factors

    M = []
    factor_base.insert(0,-1)
    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        n_factors = factor(n,factor_base)
        #print(n,n_factors)
        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(factor_base[i])) % 2

        #print(n_factors, exp_vector)
        if 1 not in exp_vector: #search for squares
            return True, n
        else:
            pass
        
        M.append(exp_vector)  
    #print("Matrix built:")
    #mprint(M)
    return(False, transpose(M))

    
def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)

'''def optimize(M):
    for row in M: #matrix optimization; delete factors that only occur once
        if row.count(1) == 1:
            for r in M:
                del r[row.index(1)]
            del row

    return(M)'''
        
def gauss_elim(M):
#reduced form of gaussian elimination, finds rref and reads off the nullspace
#https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    
    #M = optimize(M)
    marks = [False]*len(M[0])
    
    for i in range(len(M)): #do for all rows
        row = M[i]
        #print(row)
        
        for num in row: #search for pivot
            if num == 1:
                #print("found pivot at column " + str(row.index(num)+1))
                j = row.index(num) # column index
                marks[j] = True
                
                for k in chain(range(0,i),range(i+1,len(M))): #search for other 1s in the same column
                    if M[k][j] == 1:
                        for i in range(len(M[k])):
                            M[k][i] = (M[k][i] + row[i])%2
                break
    
    M = transpose(M)
        
    sol_rows = []
    for i in range(len(marks)): #find free columns (which have now become rows)
        if marks[i]== False:
            free_row = [M[i],i]
            sol_rows.append(free_row)
    
    if not sol_rows:
        return("No solution found. Need more smooth numbers.")
    print("Found {} potential solutions".format(len(sol_rows)))
    return sol_rows,marks,M

def solve_row(sol_rows,M,marks,K=0):
    solution_vec, indices = [],[]
    free_row = sol_rows[K][0] # may be multiple K
    for i in range(len(free_row)):
        if free_row[i] == 1: 
            indices.append(i)
    for r in range(len(M)): #rows with 1 in the same column will be dependent
        for i in indices:
            if M[r][i] == 1 and marks[r]:
                solution_vec.append(r)
                break
            
    solution_vec.append(sol_rows[K][1])       
    return(solution_vec)
    
def solve(solution_vec,smooth_nums,xlist,N):
    
    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]
    
    Asquare = 1
    for n in solution_nums:
        Asquare *= n
        
    b = 1
    for n in x_nums:
        b *= n

    a = isqrt(Asquare)
    
    factor = gcd(b-a,N)
    return factor


def QS(n,B,I):
#single polynomial version of quadratic sieve, given smoothness bound B and sieve interval I
    
    global N
    global root
    global T #tolerance factor
    N,root,K,T = n,int(sqrt(n)),0,1
    
    if isinstance(sqrt(N),int):
        return isqrt(N)
    
    #print(root)
    print("Attempting to factor {}...".format(N))
    #F,I = size_bound(N)
    
    print("Generating {}-smooth factor base...".format(B))
    factor_base = find_base(N,B) #generates a B-smooth factor base
    #print(factor_base)

    global F
    F = len(factor_base)
    
    print("Looking for {} {}-smooth relations...".format(F+T,B))
    smooth_nums,xlist = find_smooth(factor_base, N,I)
    #finds B-smooth relations, using sieving and Tonelli-Shanks
    
    print("Found {} B-smooth numbers.".format(len(smooth_nums)))
   
    print(smooth_nums)
    
    if len(smooth_nums) < len(factor_base):
        return("Not enough smooth numbers. Increase the sieve interval or size of the factor base.")
    
    print("Building exponent matrix...")
    is_square, t_matrix = build_matrix(smooth_nums,factor_base)
    #builds exponent matrix mod 2 from relations
    
    if is_square == True:
        x = smooth_nums.index(t_matrix)
        factor = gcd(xlist[x]+sqrt(t_matrix),N)
        print("Found a square!")
        return factor, N/factor
    
    print("Performing Gaussian Elimination...")
    sol_rows,marks,M = gauss_elim(t_matrix) #solves the matrix for the null space, finds perfect square
    solution_vec = solve_row(sol_rows,M,marks,0)
    
    '''vec = [0]*len(smooth_nums) # prints solution vector
    for i in solution_vec:
        vec[i] = 1
    print("Solution vector found: " + str(vec))'''
    
    print("Solving congruence of squares...")
    #print(solution_vec)
    factor = solve(solution_vec,smooth_nums,xlist,N) #solves the congruence of squares to obtain factors

    for K in range(1,len(sol_rows)):
        if (factor == 1 or factor == N):
            print("Didn't work. Trying different solution vector...")
            solution_vec = solve_row(sol_rows,M,marks,K)
            factor = solve(solution_vec,smooth_nums,xlist,N)
        else:
            print("Found factors!")
            return factor, int(N/factor)
            
            
    return("Didn't find any nontrivial factors!")
                   

print(QS(16921456439215439701,5000,2500000))    