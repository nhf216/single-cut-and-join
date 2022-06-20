#Lookup table
lookup = dict()

#Global variables for the different component types
N = 1
W = 2
M = 3
C = 4

label_dict = {N:'N', W:'W', M:'M', C:'C'}

#Wrapper class for info about a single component
class Component:
    def __init__(self, c_type, c_size):
        self.type = c_type
        self.size = c_size

#Wrapper class for adjacency graph enumeration
class AdjacencyGraphInfo:
    #Constructor, with a copy option
    #If not copying, construct a blank object
    def __init__(self, other = None):
        if other is not None:
            #Make a copy
            self.components = dict()
            for key in other.components:
                self.components[key] = dict(other.components[key])
            self.mindex = other.mindex
            self.gindex = other.gindex
            self.isolates = set(other.isolates)
            self.groups = dict()
            for group in other.groups:
                self.groups[group] = set(other.groups[group])
            self.groupdict = dict(other.groupdict)
            self.mergers = set(other.mergers)
        else:
            #Create a blank object
            #The components themselves
            self.components = dict()
            for key in {N, W, M, C}:
                self.components[key] = dict()
            #Indices for tracking components and groups
            self.mindex = 0
            self.gindex = 0
            #Info bits about interactions
            self.isolates = set()
            self.groups = dict()
            self.groupdict = dict()
            self.mergers = set()
    
    #Retrieve a component based on its index
    #Returns None if there's no component with that index
    def getComponent(self, index):
        for key in self.components:
            if index in self.components[key]:
                return Component(key, self.components[key][index])
    
    #Retrieve the type of a component, based on its index
    #Returns None if there's no component with that index
    def getType(self, index):
        component = self.getComponent(index)
        if component is None:
            return None
        else:
            return component.type
    
    #Count the number of components of the given type
    def countComponent(self, c_type):
        return len(self.components[c_type])
    
    #Add a component of the given type and size
    def addComponent(self, c_type, c_size):
        self.components[c_type][self.mindex] = c_size
        self.mindex += 1
        return self.mindex - 1
    
    #Add an N to the adjacency graph with l V's in it
    def addN(self, l):
        if l >= 0:
            return self.addComponent(N, l)
    
    #Add a W to the adjacency graph with l V's in it
    def addW(self, l):
        if l >= 0:
            return self.addComponent(W, l)
    
    #Add an M to the adjacency graph with l ^'s in it
    def addM(self, l):
        if l >= 0:
            return self.addComponent(M, l)
    
    #Add a crown to the adjacency graph with l V's in it
    def addC(self, l):
        if l >= 2:
            return self.addComponent(C, l)
    
    #Remove an index from group references, possibly replacing it with some new_indices
    #Do not call this method externally.
    def groupRemove(self, index, *new_indices):
        #Check if there's any work to do
        if index in self.groupdict:
            #Get the group
            group = self.groups[self.groupdict[index]]
            #Remove the index from it
            group.remove(index)
            #Do we need to add anything in place of index?
            if len(new_indices) == 0:
                #We now have a singleton group
                if len(group) == 1:
                    #Remove the remaining element from it
                    self.groupRemove(list(group)[0])
                #We now have an empty group
                elif len(group) == 0:
                    #Delete the group
                    del self.groups[self.groupdict[index]]
            else:
                #Put the new indices in
                for new_index in new_indices:
                    group.add(new_index)
                    self.groupdict[new_index] = self.groupdict[index]
            #Delete the group reference for index
            del self.groupdict[index]
    
    #Remove all references to an index, possibly replacing it with some new_indices
    #Do not call this method externally.
    #The argument c_type is there as an optimization
    def purgeIndex(self, c_type, index, *new_indices):
        #Delete this index from the components
        del self.components[c_type][index]
        #Remove/replace group references
        self.groupRemove(index, *new_indices)
        #Remove/replace isolation requirements
        if index in self.isolates:
            self.isolates.remove(index)
            for new_index in new_indices:
                self.isolates.add(new_index)
        #Remove/replace merge requirement
        if index in self.mergers:
            self.mergers.remove(index)
            for new_index in new_indices:
                self.mergers.add(new_index)
    
    #Operate on one or two components in any legal way
    #Specify the component(s) by their indices
    #Any operation involving an M requires cut_loc, an integer from 0 to 1 less than the M's size
    #An M-W cut-join requires cj_left. If true, the left piece of the M is joined to the W
    #If there are multiple ways to obtain isomorphic objects, only one is obtained
    def operateComponents(self, index, index2 = None, cut_loc = None, cj_left = None):
        component = self.getComponent(index)
        if index2 is None:
            #Single N
            if component.type == N:
                if component.size > 0:
                    #Make it smaller
                    self.components[N][index] -= 1
            #Single W
            elif component.type == W:
                if component.size > 1:
                    #Make it smaller
                    self.components[W][index] -= 1
                else:
                    #Delete it (join to trivial cycle)
                    self.purgeIndex(W, index)
            #Single M
            elif component.type == M:
                #Cut it
                idx1 = self.addN(cut_loc)
                idx2 = self.addN(component.size - cut_loc - 1)
                self.purgeIndex(M, index, idx1, idx2)
            #Single C
            else:
                #Cut it
                idx1 = self.addW(component.size)
                self.purgeIndex(C, index, idx1)
        else:
            component2 = self.getComponent(index2)
            #M-W cut-join
            if {component.type, component2.type} == {M, W}:
                #Figure out which component is which
                if component.type == W:
                    m_index = index2
                    m_size = component2.size
                    w_index = index
                    w_size = component.size
                else:
                    m_index = index
                    m_size = component.size
                    w_index = index2
                    w_size = component2.size
                #Cut-join into two N's
                self.purgeIndex(M, m_index)
                self.purgeIndex(W, w_index)
                if cj_left:
                    self.addN(m_size - cut_loc - 1)
                    self.addN(w_size + cut_loc)
                else:
                    self.addN(cut_loc)
                    self.addN(w_size + m_size - cut_loc - 1)
            #C-W cut-join
            elif {component.type, component2.type} == {C, W}:
                #Figure out which component is which
                if component.type == W:
                    c_index = index2
                    c_size = component2.size
                    w_index = index
                    w_size = component.size
                else:
                    c_index = index
                    c_size = component.size
                    w_index = index2
                    w_size = component2.size
                #Flag these as merged
                if w_index in self.mergers:
                    self.mergers.remove(w_index)
                #Overwrite the old W with the new, bigger W
                self.components[W][w_index] += c_size
                self.purgeIndex(C, c_index)
            #C-N cut-join
            elif {component.type, component2.type} == {C, N}:
                #Figure out which component is which
                if component.type == N:
                    c_index = index2
                    c_size = component2.size
                    n_index = index
                    n_size = component.size
                else:
                    c_index = index
                    c_size = component.size
                    n_index = index2
                    n_size = component2.size
                #Flag these as merged
                if n_index in self.mergers:
                    self.mergers.remove(n_index)
                #Overwrite the old N with the new, bigger N
                self.components[N][n_index] += c_size
                self.purgeIndex(C, c_index)
    
    #Force the component with the given index not to cut-join with another component
    def isolateComponent(self, index):
        self.isolates.add(index)
    
    #For the component with the given index to cut-join with another component
    #This does not allow an M to cut and later have one of the N's cut-join to a C.
    def mergeComponent(self, index):
        self.mergers.add(index)
    
    #Force the components at the given indices to sort together
    def groupComponents(self, *indices):
        #Check if valid group
        if len(indices) <= 1:
            return
        else:
            counts = dict()
            for key in {N, M, W, C}:
                counts[key] = 0
            for index in indices:
                counts[self.getType(index)] += 1
            #Options: at most 1 M, at most 1 W, any C's; one C and one N
            if not ((counts[N] == 0 and counts[M] <= 1 and counts[W] <= 1) or\
                    (counts[M] == 0 and counts[W] == 0 and counts[C] == 1 and counts[N] == 1)):
                return
        #Build the group if allowed to, including all necessary references
        gindex = self.gindex
        self.gindex += 1
        self.groups[gindex] = set(indices)
        for index in indices:
            self.groupdict[index] = gindex
        return gindex
    
    #Can we use M-W symmetry to our advantage? Dependency information may block this.
    def canSwapMW(self):
        #Make sure have none of the following:
        #-M or W in mergers
        #-M or W in isolates
        #-W with C in group
        for index in self.isolates | self.mergers:
            if self.getType(index) in {M, W}:
                return False
        for group in self.groups:
            hasC = False
            hasW = False
            for index in group:
                c_type = self.getType(index)
                if c_type == C:
                    hasC = True
                elif c_type == W:
                    hasW = True
            if hasC and hasW:
                return False
        return True
    
    #Swap the M's and W's in this adjacency graph
    def swapMW(self):
        temp = self.components[M]
        self.components[M] = self.components[W]
        self.components[W] = temp
    
    #Can we reduce the given component(s)?
    #If two indices are given, check only if they can be cut-joined together immediately.
    def canReduce(self, index, index2 = None):
        component = self.getComponent(index)
        #Single component
        if index2 is None:
            #Can't cut an M that's forced to interact
            if component.type == M and index in (set(self.groupdict) | self.mergers):
                return False
            #Can't join a V that's forced to interact
            elif component.type == W and index in (set(self.groupdict) | self.mergers) and\
                    component.size == 1:
                return False
            #Can't operate on a trivial N
            elif component.type == N and component.size == 0:
                return False
            #Can't cut a crown that already has a W it's forced to interact with
            elif component.type == C:
                #Check for the presence of other crowns
                if index in self.mergers and self.countComponent(C) == 1:
                    return False
                elif index in self.groupdict:
                    group = self.groups[self.groupdict[index]]
                    for g_index in group:
                        if self.getType(g_index) == W:
                            return False
            return True
        else:
            #Can't interact things that are required to isolate
            #FLATTEN THE CURVE!!!
            if index in self.isolates or index2 in self.isolates:
                return False
            component2 = self.getComponent(index2)
            #An M and a W
            if {component.type, component2.type} == {M, W}:
                #Only problem is if at least one is in a group without the other
                #or if both are in a group together with crowns
                if index in self.groupdict:
                    group = self.groups[self.groupdict[index]]
                    if index2 not in group or len(group) >= 3:
                        return False
                elif index2 in self.groupdict:
                    return False
                return True
            elif {component.type, component2.type} == {W, C}:
                if component.type == W:
                    c_index = index2
                    c_size = component2.size
                    w_index = index
                    w_size = component.size
                else:
                    c_index = index
                    c_size = component.size
                    w_index = index2
                    w_size = component2.size
                #Only problem is if c_index is in a group without w_index or
                #w_index is in a group
                if c_index in self.groupdict:
                    group = self.groups[self.groupdict[c_index]]
                    if w_index not in group:
                        return False
                elif w_index in self.groupdict:
                    return False
                return True
            elif {component.type, component2.type} == {N, C}:
                if component.type == N:
                    c_index = index2
                    c_size = component2.size
                    n_index = index
                    n_size = component.size
                else:
                    c_index = index
                    c_size = component.size
                    n_index = index2
                    n_size = component2.size
                #Only problem is if c_index is in a group without n_index or
                #n_index is in a group
                if c_index in self.groupdict:
                    group = self.groups[self.groupdict[c_index]]
                    if n_index not in group:
                        return False
                elif n_index in self.groupdict:
                    return False
                return True
            else:
                #Can't merge anything else
                return False
    
    #Do we only have trivial N's and trivial cycles?
    #(Note: Trivial cycles are not represented.)
    def isTrivial(self):
        if self.countComponent(M) + self.countComponent(W) + self.countComponent(C) > 0:
            return False
        for index in self.components[N]:
            if self.components[N][index] > 0:
                return False
        return True
    
    def __eq__(self, other):
        return self.components == other.components and self.isolates == other.isolates and\
            self.groups == other.groups and self.mergers == other.mergers
    
    def __str__(self):
        pieces1 = []
        for key in {N, W, M, C}:
            if self.countComponent(key) > 0:
                pieces1.append("%s: %s" % (label_dict[key], str(self.components[key])))
        pieces2 = []
        if len(self.isolates) > 0:
            pieces2.append("%s: %s" % ("I", str(self.isolates)))
        if len(self.mergers) > 0:
            pieces2.append("%s: %s" % ("M", str(self.mergers)))
        if len(self.groups) > 0:
            pieces2.append("%s: %s" % ("G", str(self.groups)))
        if len(pieces2) == 0:
            return ', '.join(pieces1)
        else:
            return ', '.join(pieces1) + ' (' + ', '.join(pieces2) + ')'
    
    def __repr__(self):
        return str(self)

#Given an AdjacencyGraphInfo object, count the number of MPS's for
#the corresponding genomes
def countMPS(g, combine_first = False, use_lookup = True, optimize_MW = True,\
             d = 0, debug = False):
    if debug:
        print("  "*d + str(g))
    #Are we at a base case?
    if g.isTrivial():
        return 1
    
    #If optimize_MW is true, if there are more M's than W's, swap them if able
    if optimize_MW and g.countComponent(M) > g.countComponent(W) and g.canSwapMW():
        new_g = AdjacencyGraphInfo(g)
        new_g.swapMW()
        return countMPS(new_g, combine_first = combine_first, use_lookup = use_lookup,\
                        optimize_MW = optimize_MW, d = d+1, debug = debug)
    
    if not combine_first and use_lookup and str(g) in lookup:
        return lookup[str(g)]
    
    answer = 0
    if not combine_first:
        #Go through the components, reducing each if allowed
        for index in g.components[N]:
            if g.canReduce(index):
                new_g = AdjacencyGraphInfo(g)
                new_g.operateComponents(index)
                answer += countMPS(new_g, use_lookup = use_lookup,\
                                   optimize_MW = optimize_MW, d = d+1, debug = debug)
        for index in g.components[W]:
            if g.canReduce(index):
                new_g = AdjacencyGraphInfo(g)
                new_g.operateComponents(index)
                multiplier = 2
                if g.components[W][index] == 1:
                    multiplier = 1
                answer += multiplier * countMPS(new_g, use_lookup = use_lookup,\
                                   optimize_MW = optimize_MW, d = d+1, debug = debug)
        for index in g.components[M]:
            if g.canReduce(index):
                size = g.components[M][index]
                for cut_loc in range(size):
                    new_g = AdjacencyGraphInfo(g)
                    new_g.operateComponents(index, cut_loc = cut_loc)
                    answer += countMPS(new_g, use_lookup = use_lookup,\
                                       optimize_MW = optimize_MW, d = d+1, debug = debug)
        for index in g.components[C]:
            if g.canReduce(index):
                size = g.components[C][index]
                new_g = AdjacencyGraphInfo(g)
                new_g.operateComponents(index)
                answer += size * countMPS(new_g, use_lookup = use_lookup,\
                                   optimize_MW = optimize_MW, d = d+1, debug = debug)
    #Do combos
    #M and W
    for m_index in g.components[M]:
        for w_index in g.components[W]:
            if g.canReduce(m_index, w_index):
                m_size = g.components[M][m_index]
                for cut_loc in range(m_size):
                    new_g = AdjacencyGraphInfo(g)
                    new_g.operateComponents(m_index, w_index, cut_loc = cut_loc, cj_left = False)
                    answer += 2 * countMPS(new_g, use_lookup = use_lookup,\
                                       optimize_MW = optimize_MW, d = d+1, debug = debug)
                    new_g = AdjacencyGraphInfo(g)
                    new_g.operateComponents(m_index, w_index, cut_loc = cut_loc, cj_left = True)
                    answer += 2 * countMPS(new_g, use_lookup = use_lookup,\
                                       optimize_MW = optimize_MW, d = d+1, debug = debug)
    for c_index in g.components[C]:
        c_size = g.components[C][c_index]
        #C and W
        for w_index in g.components[W]:
            if g.canReduce(c_index, w_index):
                new_g = AdjacencyGraphInfo(g)
                new_g.operateComponents(c_index, w_index)
                answer += 4 * c_size * countMPS(new_g, use_lookup = use_lookup,\
                                                 optimize_MW = optimize_MW, d = d+1, debug = debug)
        #C and N
        for n_index in g.components[N]:
            if g.canReduce(c_index, n_index):
                new_g = AdjacencyGraphInfo(g)
                new_g.operateComponents(c_index, n_index)
                answer += 2 * c_size * countMPS(new_g, use_lookup = use_lookup,\
                                                 optimize_MW = optimize_MW, d = d+1, debug = debug)
    
    if use_lookup and not combine_first:
        lookup[str(g)] = answer
    return answer

#Clear the memory for the lookup table
def free_lookup():
    global lookup
    lookup = dict()