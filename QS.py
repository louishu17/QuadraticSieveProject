
def sieveOfEratosthenes(n):
    isPrime = [True for i in range(n+1)]
    p = 2
    while p * p <= n:
        if isPrime[p]:
            for i in range(p*p, n+1, p):
                isPrime[p] = True
        p += 1
    primes = []

    for i in range(2, n+1):
        if isPrime[i]:
            primes.append(i)
    return primes



