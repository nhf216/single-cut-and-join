#Wrapper class for counting N, M, W, C for a given adjacency graph
class AdjacencyGraphInfo:
    #Constructor, with a copy option
    #If not copying, construct a blank object
    def __init__(self, other = None):
        if other is not None:
            self.N = dict(other.N)
            self.W = dict(other.W)
            self.M = dict(other.M)
            self.C = dict(other.C)
        else:
            self.N = dict()
            self.W = dict()
            self.M = dict()
            self.C = dict()
    
    #Add an N to the adjacency graph with l V's in it
    def addN(self, l):
        if l >= 1:
            if l not in self.N:
                self.N[l] = 0
            self.N[l] += 1
    
    #Add a W to the adjacency graph with l V's in it
    def addW(self, l):
        if l > 0:
            if l not in self.W:
                self.W[l] = 0
            self.W[l] += 1
    
    #Add an M to the adjacency graph with l ^'s in it
    def addM(self, l):
        if l > 0:
            if l not in self.M:
                self.M[l] = 0
            self.M[l] += 1
    
    #Add a crown to the adjacency graph with l V's in it
    def addC(self, l):
        if l >= 2:
            if l not in self.C:
                self.C[l] = 0
            self.C[l] += 1
    
    #Remove an N with l V's in it from the adjacency graph
    def removeN(self, l):
        self.N[l] -= 1
        if self.N[l] == 0:
            del self.N[l]
    
    #Remove a W with l V's in it from the adjacency graph
    def removeW(self, l):
        self.W[l] -= 1
        if self.W[l] == 0:
            del self.W[l]
    
    #Remove an M with l ^'s in it from the adjacency graph
    def removeM(self, l):
        self.M[l] -= 1
        if self.M[l] == 0:
            del self.M[l]
    
    #Remove a crown with l V's in it from the adjacency graph
    def removeC(self, l):
        self.C[l] -= 1
        if self.C[l] == 0:
            del self.C[l]

#Given an AdjacencyGraphInfo object, count the number of MPS's for
#the corresponding genomes
def countMPS(g, combine_first = False):
    if len(g.N) + len(g.W) + len(g.M) + len(g.C) == 0:
        return 1
    
    answer = 0
    if not combine_first:
        #Reduce an N
        for l in g.N:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeN(l)
            new_g.addN(l - 1)
            answer += g.N[l] * countMPS(new_g)

        #Reduce a W
        for l in g.W:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeW(l)
            if l == 1:
                answer += g.W[l] * countMPS(new_g)
            else:
                new_g.addW(l - 1)
                answer += 2 * g.W[l] * countMPS(new_g)

        #Reduce an M
        for l in g.M:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeM(l)
            for n in range(l):
                new_new_g = AdjacencyGraphInfo(new_g)
                new_new_g.addN(n)
                new_new_g.addN(l - n - 1)
                answer += g.M[l] * countMPS(new_new_g)

        #Reduce a C
        for l in g.C:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeC(l)
            new_g.addW(l)
            answer += l * g.C[l] * countMPS(new_g)
    
    #Combine a C and W
    for lc in g.C:
        for lw in g.W:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeC(lc)
            new_g.removeW(lw)
            new_g.addW(lc + lw)
            answer += 2 * lc * g.C[lc] * g.W[lw] * countMPS(new_g)
    
    #Combine a C and N
    for lc in g.C:
        for ln in g.N:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeC(lc)
            new_g.removeN(ln)
            new_g.addN(lc + ln)
            answer += lc * g.C[lc] * g.N[ln] * countMPS(new_g)
    
    #Combine an M and W
    for lm in g.M:
        for lw in g.W:
            new_g = AdjacencyGraphInfo(g)
            new_g.removeM(lm)
            new_g.removeW(lw)
            for n in range(lm):
                new_new_g = AdjacencyGraphInfo(new_g)
                new_new_g.addN(n)
                new_new_g.addN(lm + lw - n - 1)
                answer += 4 * g.M[lm] * g.W[lw] * countMPS(new_new_g)
    
    return answer