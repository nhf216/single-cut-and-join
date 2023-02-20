import itertools
import math

#Lookup table
lookup_table = dict()

#Standardized lookup table keys
def lookup_key(c = tuple(), m = tuple(), w = tuple(), n = tuple()):
    c = list(c)
    c.sort()
    m = list(m)
    m.sort()
    w = list(w)
    w.sort()
    n = list(n)
    n.sort()
    #Swap m and w?
    if len(w) < len(m) or w < m:
        #Yes
        t = w
        w = m
        m = t
    return (tuple(c), tuple(m), tuple(w), tuple(n))

#Distance of an MPS with the given component sizes
def distance(c = tuple(), m = tuple(), w = tuple(), n = tuple()):
    return sum(c) + len(c) + sum(m) + sum(w) + sum(n)

#Generate all partitions {0,1,2,...,n-1}
def partitions(n):
    if n == 0:
        return frozenset([frozenset()])
    #Start with partitions without n-1
    prev_parts = partitions(n - 1)
    #There's one partition with a single part
    parts = set([frozenset([frozenset(range(n))])])
    #Iterate through the previous partitions
    for part in prev_parts:
        #Insert n-1 into every possible spot
        parts.add(part|set([frozenset([n - 1])]))
        for i in range(len(part)):
            part2 = list(part)
            part2[i] |= set([n-1])
            parts.add(frozenset(part2))
    return frozenset(parts)

#Given a matching, identify all unmatched vertices
#There are n vertices on the side being checked
#If left is true, that's the left side, otherwise right
def unmatched(matching, n, left = True):
    unmatched_vertices = set(range(n))
    if left:
        index = 0
    else:
        index = 1
    for pair in matching:
        unmatched_vertices.remove(pair[index])
    return frozenset(unmatched_vertices)

#Generate all matchings in K_{m,n}, with vertices on
#both sides indexed from 0 to m-1 (resp. n-1)
def matchings(m, n):
    if m == 0 or n == 0:
        return frozenset([frozenset()])
    #Start with matchings without n-1
    prev_matchings = matchings(m, n - 1)
    #Insert n-1 wherever possible
    matches = set()
    for matching in prev_matchings:
        #This includes not matching n-1
        matches.add(matching)
        #Or matching it with any of the unmatched left vertices
        unmatched_left = unmatched(matching, m)
        for vtx in unmatched_left:
            matching2 = set(matching)
            matching2.add((vtx, n - 1))
            matches.add(frozenset(matching2))
    return matches

#Multinomial coefficient calculator
def multinomial(n, *k):
    top = n
    ret = 1
    for i in range(len(k)):
        ret *= math.comb(top, k[i])
        top -= k[i]
    return ret

#Generate all bipartitions of {0,1,2,...,n-1}
#Output is a set of tuples
def bipartitions(n):
    if n == 0:
        return frozenset([(frozenset(), frozenset())])
    #Generate all bipartitions without n-1
    prev_bips = bipartitions(n - 1)
    ret = set()
    for prev_bip in prev_bips:
        #Insert n-1 into each side
        ret.add((prev_bip[0] | {n - 1}, prev_bip[1]))
        ret.add((prev_bip[0], prev_bip[1] | {n - 1}))
    return ret

#Helper function for extracting chunks of a list R-style
def index_into(L):
    return lambda i: L[i]

#Count ways a bunch of crowns can sort together
#Could use crowns formula, but this doesn't
def count_C(c):
    key = lookup_key(c = c)
    if key in lookup_table:
        return lookup_table[key]
    ret = 0
    #Loop over all possible ways to cut a crown into a W
    for i in range(len(c)):
        ret += c[i] * count_W(c[:i] + c[i+1:], c[i])
    lookup_table[key] = ret
    return ret

#Count ways a bunch of crowns and an N can sort together
def count_N(c, n):
    key = lookup_key(c = c, n = (n,))
    if key in lookup_table:
        return lookup_table[key]
    ret = 0
    #Loop over all possible ways to cut-join a crown to the N
    for i in range(len(c)):
        ret += 2 * c[i] * count_N(c[:i] + c[i+1:], c[i] + n)
    #Add in the possibility of reducing the N
    if n > 0:
        ret += count_N(c, n - 1)
    #Is this the base case?
    if len(c) == 0 and n == 0:
        ret = 1
    lookup_table[key] = ret
    return ret

#Count ways a bunch of crowns and two N's can sort "together"
#This is only used as a helper for count_MW
def count_NN(c, n1, n2):
    key = lookup_key(c = c, n = (n1, n2))
    if key in lookup_table:
        return lookup_table[key]
    dist = distance(c = c, n = (n1, n2))
    #Bipartition the crowns to decide who goes with which crown
    bips = bipartitions(len(c))
    ret = 0
    for bip in bips:
        c1 = tuple(map(index_into(c), bip[0]))
        c2 = tuple(map(index_into(c), bip[1]))
        #Do both N counts, multiply by an appropriate binomial coefficient
        ret += math.comb(dist, distance(c = c1, n = (n1,))) * count_N(c1, n1) * count_N(c2, n2)
    lookup_table[key] = ret
    return ret

#Count ways a bunch of crowns and a W can sort together
def count_W(c, w):
    key = lookup_key(c = c, w = (w,))
    if key in lookup_table:
        return lookup_table[key]
    ret = 0
    #Loop over all possible ways to cut-join a crown to the W
    for i in range(len(c)):
        ret += 4 * c[i] * count_W(c[:i] + c[i+1:], c[i] + w)
    #Add in the possibility of reducing the W
    if w > 1:
        ret += 2 * count_W(c, w - 1)
    #Is this the base case?
    if len(c) == 0 and w == 1:
        ret = 1
    lookup_table[key] = ret
    return ret

#Count ways a bunch of crowns, an M, and a W can sort together
def count_MW(c, m, w):
    key = lookup_key(c = c, m = (m,), w = (w,))
    if key in lookup_table:
        return lookup_table[key]
    ret = 0
    #Loop over all possible ways to cut-join a crown to the W
    for i in range(len(c)):
        ret += 4 * c[i] * count_MW(c[:i] + c[i+1:], m, c[i] + w)
    #Loop over all possible ways to cut-join the M with the W
    for i in range(m):
        ret += 4 * count_NN(c, i, m + w - i - 1)
    #Add in the possibility of reducing the W
    if w > 1:
        ret += 2 * count_MW(c, m, w - 1)
    lookup_table[key] = ret
    return ret

#Count the number of scenarios for the given component sizes
#The size of a component is the floor of the number of edges over 2
def count_MPS(c = tuple(), m = tuple(), w  = tuple(), n = tuple()):
    dist = distance(c, m, w, n)
    #Pre-processing
    #Compute the lookup table
    if len(m) > 0:
        max_m = max(m)
    else:
        max_m = 0
    if len(w) > 0:
        max_max_w = max(w) + sum(c)
    else:
        max_max_w = sum(c)
    if len(n) > 0:
        max_max_n = max(n) + max_m + max_max_w
        biggest_n = max(n)
    else:
        max_max_n = max_m + max_max_w
        biggest_n = 0
    #Generate subsets of the crows in increasing order of cardinality
    for num_c in range(len(c)):
        for c_indices in itertools.combinations(range(len(c)), num_c):
            c_set = tuple(map(index_into(c), c_indices))
            #How big do we need to generate things?
            max_w = max_max_w - sum(c_set)
            max_n = max_max_n - sum(c_set)
            #Do the N's
            for i in range(max_n + 1):
                count_N(c_set, i)
            #Do the W's
            for i in range(1, max_w + 1):
                count_W(c_set, i)
            #Do the NN's
            for n1 in range(max_m):
                for n2 in range(n1, max_n + biggest_n - n1 + 1):
                    count_NN(c_set, n1, n2)
            #Do the MW's
            for i in range(1, max_m + 1):
                for j in range(1, max_w + 1):
                    count_MW(c_set, i, j)
            #Do the C's only
            if num_c > 0:
                count_C(c_set)
    
    #Do the partitioning and add everything up
    #In what follows, a "component" is a group of components of the
    #adjacency graph that all sort together
    ret = 0
    #Iterate over all possible matchings of M's and W's
    for mw_matching in matchings(len(m), len(w)):
        #Components consisting of an M and a W (and maybe some crowns)
        mw_compts = []
        for edge in mw_matching:
            mw_compts.append((m[edge[0]], w[edge[1]]))
        #How many components are there with an M and a W?
        mw_len = len(mw_compts)
        #Components consisting of an M (and maybe some crowns)
        m_compts = tuple(map(index_into(m), unmatched(mw_matching, len(m))))
        #How many components are there with an M?
        m_len = mw_len + len(m_compts)
        #Components consisting of a W (and maybe some crowns)
        w_compts = tuple(map(index_into(w), unmatched(mw_matching, len(w), False)))
        #How many components are there with an M or a W?
        w_len = m_len + len(w_compts)
        #Components consisting of an N (and maybe some crowns)
        n_compts = tuple(map(index_into(n), range(len(n))))
        #How many components are there?
        n_len = w_len + len(n_compts)
        #Iterate over all possible partitions of the crowns
        for crown_partition in partitions(len(c)):
            crown_partition = tuple(crown_partition)
            #Iterate over all possible matchings of crown partition parts with
            #non-crown components
            for cc_matching in matchings(len(crown_partition), n_len):
                #Keep track of all component types in the sorting
                compts = {"c":[], "n":[], "w":[], "mw":[]}
                #Unmatched clumps of crowns sort amongst themselves
                unmatched_c = unmatched(cc_matching, len(crown_partition))
                for c_compt in unmatched_c:
                    c_indices = crown_partition[c_compt]
                    compts["c"].append((tuple(map(index_into(c), c_indices)),))
                #Matched clumps of crowns sort with other stuff
                for edge in cc_matching:
                    c_indices = crown_partition[edge[0]]
                    if edge[1] < mw_len:
                        #It's an MW
                        compts["mw"].append((tuple(map(index_into(c), c_indices)),\
                                             *mw_compts[edge[1]]))
                    elif edge[1] < m_len:
                        #It's an M
                        #M's and W's are symmetric, so treat it as a W
                        compts["w"].append((tuple(map(index_into(c), c_indices)),\
                                            m_compts[edge[1] - mw_len]))
                    elif edge[1] < w_len:
                        #It's a W
                        compts["w"].append((tuple(map(index_into(c), c_indices)),\
                                            w_compts[edge[1] - m_len]))
                    else:
                        #It's an N
                        compts["n"].append((tuple(map(index_into(c), c_indices)),\
                                            n_compts[edge[1] - w_len]))
                #Unmatched other stuff sorts without crowns
                unmatched_other = unmatched(cc_matching, n_len, False)
                #print(mw_len, m_len, w_len, n_len)
                #print(mw_compts, m_compts, w_compts, n_compts)
                #print(unmatched_other)
                for o_compt in unmatched_other:
                    if o_compt < mw_len:
                        #It's an MW
                        compts["mw"].append((tuple(), *mw_compts[o_compt]))
                    elif o_compt < m_len:
                        #It's an M
                        #M's and W's are symmetric, so treat it as a W
                        compts["w"].append((tuple(), m_compts[o_compt - mw_len]))
                    elif o_compt < w_len:
                        #It's a W
                        compts["w"].append((tuple(), w_compts[o_compt - m_len]))
                    else:
                        #It's an N
                        compts["n"].append((tuple(), n_compts[o_compt - w_len]))
                
                #print(compts)
                
                #Finally, we can actually count things!
                #We need to track distances of each sorting component
                #along with each one's count
                dists = []
                counts = []
                #Process the crown-only clumps
                for compt in compts["c"]:
                    dists.append(distance(c = compt[0]))
                    counts.append(count_C(*compt))
                #Process the MW clumps
                for compt in compts["mw"]:
                    dists.append(distance(c = compt[0], m = (compt[1],), w = (compt[2],)))
                    counts.append(count_MW(*compt))
                #Process the W clumps
                for compt in compts["w"]:
                    dists.append(distance(c = compt[0], w = (compt[1],)))
                    counts.append(count_W(*compt))
                #Process the N clumps
                for compt in compts["n"]:
                    dists.append(distance(c = compt[0], n = (compt[1],)))
                    counts.append(count_N(*compt))
                #Accumulate the result
                ret += multinomial(dist, *dists) * math.prod(counts)
    return ret