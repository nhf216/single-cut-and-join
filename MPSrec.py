import MPScount
import math
import itertools

lookup = dict()
lookup2 = dict()

#A single m and a single w
def SMW(m, w):
    return Sbar(m, w) + 2**(m + w - 2) * math.comb(m + w, m)

#Number of ways to sort g (no cycles) without interactions
def Sstar(g, d = -1):
    pwr = 0
    if d < 0:
        d = 0
        for m in g.M:
            d += m * g.M[m]
        for w in g.W:
            d += w * g.W[w]
        for n in g.N:
            d += n * g.N[n]
    running_d = d
    ret = 1
    for m in g.M:
        pwr += g.M[m] * (m - 1)
        for i in range(g.M[m]):
            ret *= math.comb(running_d, m)
            running_d -= m
    for w in g.W:
        pwr += g.W[w] * (w - 1)
        for j in range(g.W[w]):
            ret *= math.comb(running_d, w)
            running_d -= w
    for n in g.N:
        for k in range(g.N[n]):
            ret *= math.comb(running_d, n)
            running_d -= n
    return 2**pwr * ret
    
#Defunct recurrence for S (no cycles)
def S(g):
    if str(g) in lookup:
        return lookup[str(g)]
    
    #Find the total distance and counts
    d = 0
    alpha = 0
    beta = 0
    gamma = 0
    for m in g.M:
        d += m * g.M[m]
        alpha += g.M[m]
    for w in g.W:
        d += w * g.W[w]
        beta += g.W[w]
    for n in g.N:
        d += n * g.N[n]
        gamma += g.N[n]
    
    if len(g.M) > 0:
        #Find the smallest M
        extra_m = min(g.M)

        #Do the recurrence
        ret = 0
        new_g = MPScount.AdjacencyGraphInfo(g)
        new_g.removeM(extra_m)
        for w in g.W:
            new_new_g = MPScount.AdjacencyGraphInfo(new_g)
            new_new_g.removeW(w)
            ret += g.W[w] * math.comb(d, extra_m + w) * SMW(extra_m, w) * S(new_new_g)
        ret -= (beta - 1) * Sstar(g, d)
    else:
        ret = Sstar(g, d)
    lookup[str(g)] = ret
    return ret

#Number of ways to sort an M and W that interact
def Sbar(m, w):
    ret = 0
    for k in range(w):
        total = 0
        for c in range(m):
            total += math.comb(m + w - 1 - k, c)
        ret += total * 2**(k + 2)
    return ret

#Exponentially large formula for counting MPS without constructing all of them
#See Section 5 (p. 7) in the Approximating Number of Double Cut And Join Scenarios
def S2(g):
    if str(g) in lookup2:
        return lookup2[str(g)]
    
    #Find the total distance and number of M's
    #Also, convert to lists from "multisets"
    d = 0
    alpha = 0
    M = []
    W = []
    for m in g.M:
        d += m * g.M[m]
        alpha += g.M[m]
        M += [m] * g.M[m]
    for w in g.W:
        d += w * g.W[w]
        W += [w] * g.W[w]
    for n in g.N:
        d += n * g.N[n]
    
    #Construct the answer
    ret = 0
    #First, choose the size of a subset of the M's
    for a in range(alpha + 1):
        #Now, iterate over all such subsets
        for A in itertools.combinations(M, a):
            #Next, iterate over all permutations of the W's of that size to match with the M's
            for B in itertools.permutations(W, a):
                #Build the term corresponding to this matching
                term = 1
                #Use this for incrementally building a multinomial coefficient
                running_d = d
                new_g = MPScount.AdjacencyGraphInfo(g)
                #Iterate over the matching
                for i in range(a):
                    #Match the current pair
                    #and deal with chunks where we may or may not be choosing to operate on that pair.
                    term *= Sbar(A[i], B[i]) * math.comb(running_d, A[i] + B[i])
                    new_g.removeM(A[i])
                    new_g.removeW(B[i])
                    running_d -= A[i]
                    running_d -= B[i]
                #Include independently sorting the rest
                term *= Sstar(new_g)
                #Accumulate the term
                ret += term
    lookup2[str(g)] = ret
    return ret

#Clear the memory for the lookup table
def free_lookup():
    global lookup, lookup2
    lookup = dict()
    lookup2 = dict()