import MPScount2
import random

#T means "trivial cycle"
from MPScount2 import N, M, W, C
T = 5

#Given a list of adjacency pairs (or singletons), reverse the list and its elements
def reverseAdjacencies(adj):
    ret = list(adj)
    #Reverse the list
    ret.reverse()
    #Reverse the pairs in the list
    for i in range(len(ret)):
        if len(ret[i]) == 2:
            ret[i] = (ret[i][1], ret[i][0])
    return ret

#A class representing a genome via its gene count and adjacencies
#Genes are numbered from 1 up to numGenes
#Negative means tail; positive means head
#So, -3 would mean the tail of gene 3, and 4 would mean the head of gene 4
class Genome:
    #Constructor
    #Input: Number of genes, or a genome to copy
    def __init__(self, numGenes):
        if isinstance(numGenes, Genome):
            #Make a copy
            other = numGenes
            self.numGenes = other.numGenes
            self.adjacencies = dict(other.adjacencies)
        else:
            self.numGenes = numGenes
            self.adjacencies = dict()
    
    #Perform a cut between extremities g1 and g2
    def cut(self, g1, g2):
        if g1 not in self.adjacencies or g2 not in self.adjacencies or\
                self.adjacencies[g1] != g2 or self.adjacencies[g2] != g1:
            raise ValueError("Invalid cut: %d and %d not adjacent" % (g1, g2))
        del self.adjacencies[g1]
        del self.adjacencies[g2]
    
    #Perform a join between extremities g1 and g2
    def join(self, g1, g2):
        if g1 in self.adjacencies:
            raise ValueError("Invalid join: %d already adjacent to %d" %\
                             (g1, self.adjacencies[g1]))
        elif g2 in self.adjacencies:
            raise ValueError("Invalid join: %d already adjacent to %d" %\
                             (g2, self.adjacencies[g2]))
        self.adjacencies[g1] = g2
        self.adjacencies[g2] = g1
    
    #Perform a cut-join between extremities g1, g2, and g3. g2 goes from g1 to g3.
    def cutjoin(self, g1, g2, g3):
        try:
            self.cut(g1, g2)
            self.join(g2, g3)
        except ValueError as e:
            raise ValueError("Invalid cut-join because of %s" % str(e))
    
    #Get the component containing the extremity g
    #If the component is a cycle, the gene (g, -g) will appear twice in the output
    def getComponent(self, g):
        #Start with the gene including the extremity g
        L = [(g, -g)]
        #Traverse the component from -g
        current = -g
        while True:
            if current in self.adjacencies:
                current = self.adjacencies[current]
                L.append((current, -current))
                #It's a cycle. We've cycled around.
                if abs(current) == g:
                    break
                current = -current
            else:
                #We've reached a telomere
                break
        if abs(L[-1][0]) != abs(g) or len(L) == 1:
            #Need to go the other way too
            current = g
            while True:
                if current in self.adjacencies:
                    current = self.adjacencies[current]
                    L.insert(0, (-current, current))
                    current = -current
                else:
                    #We've reached the other telomere
                    break
        return tuple(L)
    
    #Retrieve all components of the genome, each one exactly once.
    def getComponents(self):
        components = set()
        used_genes = set()
        for g in range(1, self.numGenes + 1):
            if g not in used_genes:
                component = self.getComponent(g)
                components.add(component)
                for gene in component:
                    used_genes.add(abs(gene[0]))
        return components
    
    def __str__(self):
        return str(self.getComponents())
    
    def __repr__(self):
        return str(self)

#Class representing a component of an adjacency graph
#A Component has a type and two adjacency lists of extremities.
#The first adjacacency list is for the top genome.
#The second is for the bottom genome.
class Component:
    #Constructor
    def __init__(self, c_type, adj1 = None, adj2 = None):
        if adj1 is None:
            #Make a copy
            other = c_type
            self.type = other.type
            self.adj1 = other.adj1
            self.adj2 = other.adj2
        else:
            #Build a new Component object
            self.type = c_type
            self.adj1 = adj1
            self.adj2 = adj2
    
    #The index of a Component is the lowest numbered extremity (ignoring head/tail)
    #If both head and tail of lowest gene number occur, the tail is the index.
    def getIndex(self):
        #Track both minima, as positive numbers
        #0 is a placeholder
        min_head = 0
        min_tail = 0
        for adj in self.adj1:
            for g in adj:
                if g > 0:
                    if min_head == 0:
                        min_head = g
                    else:
                        min_head = min(g, min_head)
                else:
                    if min_tail == 0:
                        min_tail = -g
                    else:
                        min_tail = min(-g, min_tail)
        #Different cases for whether to return a head or a tail
        if min_tail == 0:
            return min_head
        elif min_head == 0:
            return -min_tail
        elif min_head < -min_tail:
            return min_head
        else:
            return -min_tail
    
    #Determine which genes are in this component
    def getGenes(self):
        ret = set()
        for adj in self.adj1:
            for g in adj:
                ret.add(g)
        return ret
    
    #Invert the component. In other words, flip it upside down.
    #This is done via a 180 degree rotation, not a reflection.
    #This way, N's always have their upper extremity on the right.
    def getInverse(self):
        #Interchange M and W; keep others the same
        if self.type == M:
            c_type = W
        elif self.type == W:
            c_type = M
        else:
            c_type = self.type
        #Reverse and interchange the adjacencies
        return Component(c_type, reverseAdjacencies(self.adj2),\
                         reverseAdjacencies(self.adj1))
    
    #The size of a component is the floor of half the number of genes
    def getSize(self):
        return len(self.getGenes()) // 2
    
    def __str__(self):
        if self.type == T:
            label = "T"
        else:
            label = MPScount2.label_dict[self.type]
        return label + " " + str((self.adj1, self.adj2))
    
    def __repr__(self):
        return str(self)

#An abstraction of an operation on genomes
#Members:
#-op: "C" for cut, "J" for join, "CJ" for cut-join
#-type1 and type2: The type (M, W, C, N) of the first component
#(M or C always first; type2 only used if a second component is involved)
#-index1 and index2: The indices (see above) of the components
#-cut_loc: For an M or C, location where to cut (indexed from 0 from the left)
#-left1 and left2: For M and C, is the cut part joined to the other component on the left or right?
#For W, which side is being joined to or cut-joined?
#-gene1, gene2, and gene3: The extremities involved in the operation.
#(gene3 only used for cut-join)
class Operation:
    #Constructor
    #Make an empty Operation
    def __init__(self):
        self.type1 = None
        self.index1 = None
        self.cut_loc = None
        self.left1 = None
        self.type2 = None
        self.index2 = None
        self.left2 = None
        self.gene1 = None
        self.gene2 = None
        self.gene3 = None
        self.op = None
    
    #Get a string representation of this operation
    #The string depends on all the different pieces and whether they exist
    def __str__(self):
        genes = []
        if self.gene1 is not None:
            genes.append(self.gene1)
        if self.gene2 is not None:
            genes.append(self.gene2)
        if self.gene3 is not None:
            genes.append(self.gene3)
        if self.type1 is None:
            comp1 = None
        else:
            comp1 = [MPScount2.label_dict[self.type1]]
            if self.index1 is not None:
                comp1.append(self.index1)
            if self.cut_loc is not None:
                comp1.append(self.cut_loc)
            if self.left1 == False:
                comp1.append("R")
            elif self.left1 == True:
                comp1.append("L")
        if self.type2 is None:
            comp2 = None
        else:
            comp2 = [MPScount2.label_dict[self.type2]]
            if self.index2 is not None:
                comp2.append(self.index2)
            if self.left2 == False:
                comp2.append("R")
            elif self.left2 == True:
                comp2.append("L")
        if comp1 is None:
            return "%s %s" % (self.op, str(tuple(genes)))
        elif comp2 is None:
            return "%s %s %s" % (self.op, str(tuple(genes)), str(tuple(comp1)))
        else:
            return "%s %s %s %s" % (self.op, str(tuple(genes)), str(tuple(comp1)),\
                                    str(tuple(comp2)))
    
    def __repr__(self):
        return str(self)

#An adjacency graph
#Tracks both the top and bottom genomes along with the components of the graph
class AdjacencyGraph:
    #Constructor
    #aginfo is either another AdjacencyGraph object to copy, or it's an AdjacencyGraphInfo
    #object from MPScount2. aginfo is stored for population purposes, but it is never updated later.
    #If shuffle is True, randomize the order of the extremities
    #If populate is True, construct the genomes right away
    def __init__(self, aginfo, shuffle = False, populate = True):
        if isinstance(aginfo, AdjacencyGraph):
            #Make a copy
            other = aginfo
            self.g = other.g
            self.genome1 = Genome(other.genome1)
            self.genome2 = Genome(other.genome2)
            self.components = dict()
            self.assoc = dict()
            for key in {N, W, M, C, T}:
                self.components[key] = dict()
                for o_comp in other.components[key].values():
                    comp = Component(o_comp)
                    self.components[key][comp.getIndex()] = comp
                    for index in comp.getGenes():
                        self.assoc[index] = comp
        else:
            if len(aginfo.components[MPScount2.N]) % 2 == 1:
                aginfo = MPScount2.AdjacencyGraphInfo(aginfo)
                aginfo.addN(0)
            self.g = aginfo
            self.components = dict()
            for key in {N, W, M, C, T}:
                self.components[key] = dict()
            self.assoc = dict()
            if populate:
                self.populateGenomes(shuffle)
    
    #Invert this graph by inverting all components and swapping the genomes
    def getInverse(self):
        g = MPScount2.AdjacencyGraphInfo(self.g)
        g.swapMW()
        ret = AdjacencyGraph(g, populate = False)
        ret.genome1 = Genome(self.genome2)
        ret.genome2 = Genome(self.genome1)
        ret.components = dict()
        ret.assoc = dict()
        new_key = {N:N, M:W, W:M, C:C, T:T}
        for key in new_key:
            ret.components[new_key[key]] = dict()
            for idx in self.components[key]:
                ret.components[new_key[key]][idx] = self.components[key][idx].getInverse()
            for comp in ret.components[new_key[key]].values():
                for index in comp.getGenes():
                    ret.assoc[index] = comp
        return ret
    
    #Use the AdjacencyGraphInfo to construct the actual genomes
    #so that we get a graph that agrees with the component counts
    def populateGenomes(self, shuffle = False):
        #How many genes do we need?
        numGenes = len(self.g.components[MPScount2.N]) // 2 +\
            sum(self.g.components[N].values()) +\
            sum(self.g.components[W].values()) +\
            sum(self.g.components[M].values()) +\
            sum(self.g.components[C].values())
        #Create two genomes with no adjacencies initially
        genome1 = Genome(numGenes)
        genome2 = Genome(numGenes)
        #Create a list of all extremities needed
        genes = []
        for g in range(1, numGenes + 1):
            genes.append(-g)
            genes.append(g)
        #Shuffle the extremities if asked to
        if shuffle:
            random.shuffle(genes)
        #Track where in the list of extremities we are
        gene_index = 0
        #Process the N's
        #Always make the upper telomere be on the right side
        for n in self.g.components[N].values():
            start_index = gene_index
            adj1 = []
            adj2 = []
            adj2.append((genes[gene_index],))
            for i in range(n):
                genome1.join(genes[gene_index], genes[gene_index + 1])
                adj1.append((genes[gene_index], genes[gene_index + 1]))
                genome2.join(genes[gene_index + 1], genes[gene_index + 2])
                adj2.append((genes[gene_index + 1], genes[gene_index + 2]))
                gene_index += 2
            adj1.append((genes[gene_index],))
            comp = Component(N, adj1, adj2)
            self.components[N][comp.getIndex()] = comp
            for i in comp.getGenes():
                self.assoc[i] = comp
            gene_index += 1
        #Process the W's
        for w in self.g.components[W].values():
            start_index = gene_index
            adj1 = []
            adj2 = []
            adj1.append((genes[gene_index],))
            for i in range(w):
                if i < w - 1:
                    genome1.join(genes[gene_index + 1], genes[gene_index + 2])
                    adj1.append((genes[gene_index + 1], genes[gene_index + 2]))
                else:
                    adj1.append((genes[gene_index + 1],))
                genome2.join(genes[gene_index], genes[gene_index + 1])
                adj2.append((genes[gene_index], genes[gene_index + 1]))
                gene_index += 2
            comp = Component(W, adj1, adj2)
            self.components[W][comp.getIndex()] = comp
            for i in comp.getGenes():
                self.assoc[i] = comp
        #Process the M's
        for m in self.g.components[M].values():
            start_index = gene_index
            adj1 = []
            adj2 = []
            adj2.append((genes[gene_index],))
            for i in range(m):
                genome1.join(genes[gene_index], genes[gene_index + 1])
                adj1.append((genes[gene_index], genes[gene_index + 1]))
                if i < m - 1:
                    genome2.join(genes[gene_index + 1], genes[gene_index + 2])
                    adj2.append((genes[gene_index + 1], genes[gene_index + 2]))
                else:
                    adj2.append((genes[gene_index + 1],))
                gene_index += 2
            comp = Component(M, adj1, adj2)
            self.components[M][comp.getIndex()] = comp
            for i in comp.getGenes():
                self.assoc[i] = comp
        #Process the C's
        for c in self.g.components[C].values():
            start_index = gene_index
            adj1 = []
            adj2 = []
            for i in range(c):
                genome1.join(genes[gene_index], genes[gene_index + 1])
                adj1.append((genes[gene_index], genes[gene_index + 1]))
                if i < c - 1:
                    genome2.join(genes[gene_index + 1], genes[gene_index + 2])
                    adj2.append((genes[gene_index + 1], genes[gene_index + 2]))
                gene_index += 2
            genome2.join(genes[start_index], genes[gene_index - 1])
            adj2.append((genes[gene_index - 1], genes[start_index]))
            comp = Component(C, adj1, adj2)
            self.components[C][comp.getIndex()] = comp
            for i in comp.getGenes():
                self.assoc[i] = comp
        self.genome1 = genome1
        self.genome2 = genome2
    
    #Get a list of all of the components in this graph
    def getComponents(self):
        return list(self.components[N].values()) + list(self.components[W].values()) +\
            list(self.components[M].values()) + list(self.components[C].values()) +\
            list(self.components[T].values())
    
    #Is the graph trivial?
    #That is, does it consist only of trivial N's and trivial cycles (T's)?
    def isTrivial(self):
        if len(self.components[W]) + len(self.components[M]) + len(self.components[C]) > 0:
            return False
        for comp in self.components[N].values():
            if comp.getSize() > 0:
                return False
        return True
    
    #Do a cut, and update the adjacency graph
    #Only works if the cut is part of an MPS
    #Returns an Operation object documenting the cut
    def cut(self, g1, g2):
        #Cut the genome
        self.genome1.cut(g1, g2)
        #What component of the adjacency graph are we cutting?
        comp = self.assoc[g1]
        #Build the operation to return
        ret = Operation()
        #Register that it's a cut on extremities g1 and g2
        ret.op = "C"
        ret.gene1 = g1
        ret.gene2 = g2
        if comp.type == M:
            #Option 1: We're cutting an M
            #The result is two N's
            adj1l = []
            adj2l = [comp.adj2[0]]
            #Find where the cut is happening
            for i in range(len(comp.adj1)):
                h1, h2 = comp.adj1[i]
                if {h1, h2} == {g1, g2}:
                    adj1l.append((h1,))
                    adj1r = [(h2,)]
                    adj1r += comp.adj1[i+1:]
                    adj2r = comp.adj2[i+1:]
                    #The right N is "backwards," so we reverse it
                    adj1r = reverseAdjacencies(adj1r)
                    adj2r = reverseAdjacencies(adj2r)
                    cut_loc = i
                    break
                else:
                    adj1l.append(comp.adj1[i])
                    adj2l.append(comp.adj2[i+1])
            #Make the N's
            comp1 = Component(N, adj1l, adj2l)
            comp2 = Component(N, adj1r, adj2r)
            #Clean up references to the old/new components
            del self.components[M][comp.getIndex()]
            self.components[N][comp1.getIndex()] = comp1
            self.components[N][comp2.getIndex()] = comp2
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            for g in comp2.getGenes():
                self.assoc[g] = comp2
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = M
            ret.index1 = comp.getIndex()
            ret.cut_loc = cut_loc
            return ret
        else:
            #Option 2: We're cutting a C
            #The result is a W
            #Find where the cut is happening
            for i in range(len(comp.adj1)):
                h1, h2 = comp.adj1[i]
                if {h1, h2} == {g1, g2}:
                    cut_loc = i
                    break
            #Build the W
            adj1 = [(h2,)]
            adj2 = [comp.adj2[cut_loc]]
            for i in range(1, len(comp.adj1)):
                adj1.append(comp.adj1[(cut_loc + i) % len(comp.adj1)])
                adj2.append(comp.adj2[(cut_loc + i) % len(comp.adj2)])
            adj1.append((h1,))
            comp1 = Component(W, adj1, adj2)
            #Clean up references to the old/new components
            del self.components[C][comp.getIndex()]
            self.components[W][comp1.getIndex()] = comp1
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = C
            ret.index1 = comp.getIndex()
            ret.cut_loc = cut_loc
            return ret
    
    #Do a join, and update the adjacency graph
    #Only works if the join is part of an MPS
    #Returns an Operation object documenting the join
    def join(self, g1, g2):
        #Join the genome
        self.genome1.join(g1, g2)
        #What component of the adjacency graph are we joining within?
        comp = self.assoc[g1]
        #Build the operation to return
        ret = Operation()
        #Register that it's a join on extremities g1 and g2
        ret.op = "J"
        ret.gene1 = g1
        ret.gene2 = g2
        #Only option is a W of size 1
        #The result is a trivial cycle (T)
        comp1 = Component(T, list(comp.adj2), list(comp.adj2))
        #Clean up references to the old/new components
        del self.components[W][comp.getIndex()]
        self.components[T][comp1.getIndex()] = comp1
        for g in comp1.getGenes():
            self.assoc[g] = comp1
        #Fill in the rest of the info in the output object, and return it
        ret.type1 = W
        ret.index1 = comp.getIndex()
        return ret
    
    #Do a cut-join, and update the adjacency graph
    #Only works if the cut-join is part of an MPS
    #Returns an Operation object documenting the cut-join
    def cutjoin(self, g1, g2, g3):
        #Cut-join the genome
        self.genome1.cutjoin(g1, g2, g3)
        #What component(s) of the adjacency graph are we joining within?
        compA = self.assoc[g1]
        compB = self.assoc[g3]
        #Build the operation to return
        ret = Operation()
        #Register that it's a cut-join on extremities g1, g2, and g3
        ret.op = "CJ"
        ret.gene1 = g1
        ret.gene2 = g2
        ret.gene3 = g3
        if compA == compB and compA.type == N:
            #We're cut-joining within an N
            #Both components are the same N
            #The result is a smaller N and a T
            comp = compA
            adj1 = comp.adj1[:-1]
            adj1[-1] = adj1[-1][:1]
            comp1 = Component(N, adj1, comp.adj2[:-1])
            comp2 = Component(T, [comp.adj2[-1]], [comp.adj2[-1]])
            #Clean up references to the old/new components
            del self.components[N][comp.getIndex()]
            self.components[N][comp1.getIndex()] = comp1
            self.components[T][comp2.getIndex()] = comp2
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            for g in comp2.getGenes():
                self.assoc[g] = comp2
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = N
            ret.index1 = comp.getIndex()
            return ret
        elif compA == compB and compA.type == W:
            #We're cut-joining within a W
            #Both components are the same W
            #The result is a smaller W
            comp = compA
            #Which side of the W are we operating on?
            left = g3 in comp.adj1[0]
            adj1 = list(comp.adj1)
            adj2 = list(comp.adj2)
            if left:
                comp2 = Component(T, [adj2[0]], [adj2[0]])
                adj1 = adj1[1:]
                adj1[0] = adj1[0][1:]
                comp1 = Component(W, adj1, adj2[1:])
            else:
                comp2 = Component(T, [adj2[-1]], [adj2[-1]])
                adj1 = adj1[:-1]
                adj1[-1] = adj1[-1][:1]
                comp1 = Component(W, adj1, adj2[:-1])
            #Clean up references to the old/new components
            del self.components[W][comp.getIndex()]
            self.components[W][comp1.getIndex()] = comp1
            self.components[T][comp2.getIndex()] = comp2
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            for g in comp2.getGenes():
                self.assoc[g] = comp2
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = W
            ret.index1 = comp.getIndex()
            ret.left1 = left
            return ret
        elif compA.type == M:
            #We're combining an M and a W
            #compB has to be a W
            adj1l = []
            adj2l = [compA.adj2[0]]
            #Figure out where we're cutting the M
            for i in range(len(compA.adj1)):
                h1, h2 = compA.adj1[i]
                if {h1, h2} == {g1, g2}:
                    adj1r = compA.adj1[i+1:]
                    adj2r = compA.adj2[i+1:]
                    adj1r = reverseAdjacencies(adj1r)
                    adj2r = reverseAdjacencies(adj2r)
                    cut_loc = i
                    break
                else:
                    adj1l.append(compA.adj1[i])
                    adj2l.append(compA.adj2[i+1])
            #Figure out which half of the M got joined
            leftM = g2 == h1
            if leftM:
                adj1r.append((h2,))
                #The result of this operation is two N's
                #Construct one of the N's now
                comp2 = Component(N, adj1r, adj2r)
                adj1 = adj1l
                adj2 = adj2l
                g = h1
            else:
                adj1l.append((h1,))
                #The result of this operation is two N's
                #Construct one of the N's now
                comp2 = Component(N, adj1l, adj2l)
                adj1 = adj1r
                adj2 = adj2r
                g = h2
            #Figure out which end of the W got joined
            leftW = g3 in compB.adj1[0]
            if leftW:
                adj1.append((g, compB.adj1[0][0]))
                adj1 += compB.adj1[1:]
                adj2 += compB.adj2
            else:
                adj1w = reverseAdjacencies(compB.adj1)
                adj2w = reverseAdjacencies(compB.adj2)
                adj1.append((g, adj1w[0][0]))
                adj1 += adj1w[1:]
                adj2 += adj2w
            #Construct the other N
            comp1 = Component(N, adj1, adj2)
            #Clean up references to the old/new components
            del self.components[M][compA.getIndex()]
            del self.components[W][compB.getIndex()]
            self.components[N][comp1.getIndex()] = comp1
            self.components[N][comp2.getIndex()] = comp2
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            for g in comp2.getGenes():
                self.assoc[g] = comp2
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = M
            ret.index1 = compA.getIndex()
            ret.cut_loc = cut_loc
            ret.left1 = leftM
            ret.type2 = W
            ret.index2 = compB.getIndex()
            ret.left2 = leftW
            return ret
        elif compB.type == W:
            #We're combining a C and a W
            #compA has to be a C
            #Figure out where we're cutting the C
            for i in range(len(compA.adj1)):
                h1, h2 = compA.adj1[i]
                if {h1, h2} == {g1, g2}:
                    cut_loc = i
                    break
            #Actually cut it
            adj1 = [(h2,)]
            adj2 = [compA.adj2[cut_loc]]
            for i in range(1, len(compA.adj1)):
                adj1.append(compA.adj1[(cut_loc + i) % len(compA.adj1)])
                adj2.append(compA.adj2[(cut_loc + i) % len(compA.adj2)])
            adj1.append((h1,))
            #Figure out which side of the C got joined
            leftC = g2 == h1
            if leftC:
                g = h1
            else:
                adj1 = reverseAdjacencies(adj1)
                adj2 = reverseAdjacencies(adj2)
                g = h2
            #Figure out which end of the W got joined
            leftW = g3 in compB.adj1[0]
            if leftW:
                adj1[-1] = (g, compB.adj1[0][0])
                adj1 += compB.adj1[1:]
                adj2 += compB.adj2
            else:
                adj1w = reverseAdjacencies(compB.adj1)
                adj2w = reverseAdjacencies(compB.adj2)
                adj1[-1] = (g, adj1w[0][0])
                adj1 += adj1w[1:]
                adj2 += adj2w
            #The result is a W. Construct it.
            comp1 = Component(W, adj1, adj2)
            #Clean up references to the old/new components
            del self.components[C][compA.getIndex()]
            del self.components[W][compB.getIndex()]
            self.components[W][comp1.getIndex()] = comp1
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = C
            ret.index1 = compA.getIndex()
            ret.cut_loc = cut_loc
            ret.left1 = leftC
            ret.type2 = W
            ret.index2 = compB.getIndex()
            ret.left2 = leftW
            return ret
        else:
            #We're combining a C and an N
            #compA has to be a C; compB has to be an N
            #Figure out where we're cutting the C
            for i in range(len(compA.adj1)):
                h1, h2 = compA.adj1[i]
                if {h1, h2} == {g1, g2}:
                    cut_loc = i
                    break
            #Actually cut it
            adj1 = [(h2,)]
            adj2 = [compA.adj2[cut_loc]]
            for i in range(1, len(compA.adj1)):
                adj1.append(compA.adj1[(cut_loc + i) % len(compA.adj1)])
                adj2.append(compA.adj2[(cut_loc + i) % len(compA.adj2)])
            adj1.append((h1,))
            #Figure out which side of the C got joined
            leftC = g2 == h1
            if leftC:
                g = h1
            else:
                adj1 = reverseAdjacencies(adj1)
                adj2 = reverseAdjacencies(adj2)
                g = h2
            #Build the resulting N
            adj1n = reverseAdjacencies(compB.adj1)
            adj2n = reverseAdjacencies(compB.adj2)
            adj1[-1] = (g, adj1n[0][0])
            adj1 += adj1n[1:]
            adj2 += adj2n
            adj1 = reverseAdjacencies(adj1)
            adj2 = reverseAdjacencies(adj2)
            comp1 = Component(N, adj1, adj2)
            #Clean up references to the old/new components
            del self.components[C][compA.getIndex()]
            del self.components[N][compB.getIndex()]
            self.components[N][comp1.getIndex()] = comp1
            for g in comp1.getGenes():
                self.assoc[g] = comp1
            #Fill in the rest of the info in the output object, and return it
            ret.type1 = C
            ret.index1 = compA.getIndex()
            ret.cut_loc = cut_loc
            ret.left1 = leftC
            ret.type2 = N
            ret.index2 = compB.getIndex()
            return ret
    
    #Reduce the N containing the given extremity (gene).
    def reduceN(self, gene):
        comp = self.assoc[gene]
        if len(comp.adj1) > 1:
            #Only one option
            g1 = comp.adj1[-2][0]
            g2 = comp.adj1[-2][1]
            g3 = comp.adj1[-1][0]
            return self.cutjoin(g1, g2, g3)
    
    #Reduce the M containing the given extremity (gene).
    #Cut it at either the given location (cut_loc) or extremity (cut_gene).
    def reduceM(self, gene, cut_loc = 0, cut_gene = None):
        comp = self.assoc[gene]
        if cut_gene is not None:
            for i in range(len(comp.adj1)):
                if cut_gene in comp.adj1[i]:
                    cut_loc = i
                    break
        g1 = comp.adj1[cut_loc][0]
        g2 = comp.adj1[cut_loc][1]
        return self.cut(g1, g2)
    
    #Reduce the W containing the given extremity (gene).
    #Cut-join it on the specified side.
    def reduceW(self, gene, leftW = True):
        comp = self.assoc[gene]
        if len(comp.adj2) == 1:
            #Only one option
            g1 = comp.adj2[0][0]
            g2 = comp.adj2[0][1]
            return self.join(g1, g2)
        else:
            if leftW:
                g1 = comp.adj1[1][1]
                g2 = comp.adj1[1][0]
                g3 = comp.adj1[0][0]
            else:
                g1 = comp.adj1[-2][0]
                g2 = comp.adj1[-2][1]
                g3 = comp.adj1[-1][0]
            return self.cutjoin(g1, g2, g3)
    
    #Reduce the C containing the given extremity (gene).
    #Cut it at either the given location (cut_loc) or extremity (cut_gene).
    def reduceC(self, gene, cut_loc = 0, cut_gene = None):
        #Same exact action as cutting an M
        return self.reduceM(gene, cut_loc, cut_gene)
    
    #Reduce the M and W together containing respectively the given extremities (geneM and geneW)
    #Cut the M at either the given location (cut_loc) or extremity (cut_gene).
    #Join the specified piece of the cut M to the specified side of the W.
    def reduceMW(self, geneM, geneW, cut_loc = 0, cut_gene = None,\
                 leftM = True, leftW = True):
        compM = self.assoc[geneM]
        compW = self.assoc[geneW]
        if cut_gene is not None:
            for i in range(len(compM.adj1)):
                if cut_gene in compM.adj1[i]:
                    cut_loc = i
                    break
        if leftM:
            g1 = compM.adj1[cut_loc][1]
            g2 = compM.adj1[cut_loc][0]
            if leftW:
                g3 = compW.adj1[0][0]
            else:
                g3 = compW.adj1[-1][0]
        else:
            g1 = compM.adj1[cut_loc][0]
            g2 = compM.adj1[cut_loc][1]
            if leftW:
                g3 = compW.adj1[0][0]
            else:
                g3 = compW.adj1[-1][0]
        return self.cutjoin(g1, g2, g3)
    
    #Reduce the C and W together containing respectively the given extremities (geneC and geneW)
    #Cut the C at either the given location (cut_loc) or extremity (cut_gene).
    #Join the specified piece of the cut C to the specified side of the W.
    def reduceCW(self, geneC, geneW, cut_loc = 0, cut_gene = None,\
                 leftC = True, leftW = True):
        #Same exact action as cut-joining an M and a W
        return self.reduceMW(geneC, geneW, cut_loc, cut_gene, leftC, leftW)
    
    #Reduce the C and N together containing respectively the given extremities (geneC and geneN)
    #Cut the C at either the given location (cut_loc) or extremity (cut_gene).
    #Join the specified piece of the cut C to the N.
    def reduceCN(self, geneC, geneN, cut_loc = 0, cut_gene = None, leftC = True):
        compC = self.assoc[geneC]
        compN = self.assoc[geneN]
        if cut_gene is not None:
            for i in range(len(compC.adj1)):
                if cut_gene in compC.adj1[i]:
                    cut_loc = i
                    break
        if leftC:
            g1 = compC.adj1[cut_loc][1]
            g2 = compC.adj1[cut_loc][0]
            g3 = compN.adj1[-1][0]
        else:
            g1 = compC.adj1[cut_loc][0]
            g2 = compC.adj1[cut_loc][1]
            g3 = compN.adj1[-1][0]
        return self.cutjoin(g1, g2, g3)
    
    def __str__(self):
        if self.genome1 is None:
            return str(self.g)
        else:
            ret = [str(self.genome1), str(self.genome2)]
            for comp in self.getComponents():
                ret.append(str(comp))
            return '\n'.join(ret)
    
    def __repr__(self):
        return str(self)

#Obtain all most parsimonious scenarios from the given adjacency graph
#This operation does NOT mutate adjGraph
def getAllMPS(adjGraph):
    #Base case. If it's trivial, return the empty scenario.
    if adjGraph.isTrivial():
        return [[]]
    #Aggregate all scenarios
    ret = []
    #Reduce all the nontrivial N's
    for idx in adjGraph.components[N]:
        if adjGraph.components[N][idx].getSize() > 0:
            g = AdjacencyGraph(adjGraph)
            red = g.reduceN(idx)
            for mps in getAllMPS(g):
                ret.append([red] + mps)
    #Reduce all the W's
    for idx in adjGraph.components[W]:
        if adjGraph.components[W][idx].getSize() == 1:
            #Its size is 1, so join it
            lst = {None}
        else:
            #Iterate through both possible reductions
            lst = {True, False}
        for left in lst:
            g = AdjacencyGraph(adjGraph)
            red = g.reduceW(idx, leftW = left)
            for mps in getAllMPS(g):
                ret.append([red] + mps)
    #Reduce all the M's
    for idx in adjGraph.components[M]:
        #Iterate through all possible cut locations
        for c in range(adjGraph.components[M][idx].getSize()):
            g = AdjacencyGraph(adjGraph)
            red = g.reduceM(idx, c)
            for mps in getAllMPS(g):
                ret.append([red] + mps)
    #Reduce all the C's
    for idx in adjGraph.components[C]:
        #Iterate through all possible cut locations
        for c in range(adjGraph.components[C][idx].getSize()):
            g = AdjacencyGraph(adjGraph)
            red = g.reduceC(idx, c)
            for mps in getAllMPS(g):
                ret.append([red] + mps)
    #Reduce all the MW pairs
    #Iterate through the M's
    for idx1 in adjGraph.components[M]:
        #Iterate through all possible cut locations
        for c in range(adjGraph.components[M][idx].getSize()):
            #Iterate through both choices of M piece
            for leftM in {True, False}:
                #Iterate through the W's
                for idx2 in adjGraph.components[W]:
                    #Iterate through both choices of W side
                    for leftW in {True, False}:
                        g = AdjacencyGraph(adjGraph)
                        red = g.reduceMW(idx1, idx2, cut_loc = c, leftM = leftM, leftW = leftW)
                        for mps in getAllMPS(g):
                            ret.append([red] + mps)
    #Reduce all the CW pairs
    #Iterate through the C's
    for idx1 in adjGraph.components[C]:
        #Iterate through all possible cut locations
        for c in range(adjGraph.components[C][idx].getSize()):
            #Iterate through both choices of C end
            for leftC in {True, False}:
                #Iterate through the W's
                for idx2 in adjGraph.components[W]:
                    #Iterate through both choices of W side
                    for leftW in {True, False}:
                        g = AdjacencyGraph(adjGraph)
                        red = g.reduceCW(idx1, idx2, cut_loc = c, leftC = leftC, leftW = leftW)
                        for mps in getAllMPS(g):
                            ret.append([red] + mps)
    #Reduce all the CN pairs
    #Iterate through the C's
    for idx1 in adjGraph.components[C]:
        #Iterate through all possible cut locations
        for c in range(adjGraph.components[C][idx].getSize()):
            #Iterate through both choices of C end
            for leftC in {True, False}:
                #Iterate through the N's
                for idx2 in adjGraph.components[N]:
                    g = AdjacencyGraph(adjGraph)
                    red = g.reduceCN(idx1, idx2, cut_loc = c, leftC = leftC)
                    for mps in getAllMPS(g):
                        ret.append([red] + mps)
    return ret

#Given an adjacency graph and list of scenarios, get the dual scenario for each one
def getReversals(adjGraph, mpsList):
    #Invert the graph
    invGraph = adjGraph.getInverse()
    #Start with an empty list
    ret = []
    #Iterate through the scenarios
    for mps in mpsList:
        #Copy the graph to avoid mutating it
        g = AdjacencyGraph(invGraph)
        rev_mps = []
        #Iterate through the operations in verse
        for i in range(len(mps) - 1, -1, -1):
            op = mps[i]
            #Do the inverse operation on the inverse graph and record the result
            if op.op == "CJ":
                rev_mps.append(g.cutjoin(op.gene3, op.gene2, op.gene1))
            elif op.op == "C":
                rev_mps.append(g.join(op.gene2, op.gene1))
            else:
                rev_mps.append(g.cut(op.gene2, op.gene1))
        ret.append(rev_mps)
    return ret