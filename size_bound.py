from math import exp, log, sqrt

def size_bound(N): # finds optimal factor base size and interval

    F = pow(exp(sqrt(log(N)*log(log(N)))),sqrt(2)/4)
    I = F**3
    #print(F,I)
    return int(F),int(I)

print(size_bound(16921456439215439701))