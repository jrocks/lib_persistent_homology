import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import itertools as it
# import numpy.linalg as la
import phat
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import matplotlib as mpl
# import skimage.morphology as morph
# import networkx as nx
import queue
import collections as co

class Complex:
    
    def __init__(self, dim, ordered=True):
        self.dim = dim
        self.ordered = ordered
        if self.ordered:
            self.dims = []
            self.facets = []
        else:
            self.dims = {}
            self.facets = {}
        
    def add_cell(self, d, f, label=None):
        if self.ordered:
            if label is None:
                self.dims.append(d)
                self.facets.append(f)
            else:
                self.dims.insert(label, d)
                self.facets.insert(label, f)
        else:
            if label is None:
                self.dims[len(self.dims)] = d
                self.facets[len(self.facets)] = f
            else:
                self.dims[label] = d
                self.facets[label] = f
            
    def get_cells(self):
        
        if self.ordered:
            for i in range(len(self.facets)):
                yield i
                
        else:
            for i in self.facets:
                yield i
        
    def construct_cofacets(self):
        
        if self.ordered:
            self.cofacets = [set() for i in range(len(self.facets))]

            for i in range(len(self.facets)):
                for j in self.facets[i]:
                    self.cofacets[j].add(i)
                    
        else:
            self.cofacets = {i:{} for i in self.facets}
            
            for i in self.facets:
                for j in self.facets[i]:
                    self.cofacets[j].add(i) 
            

def construct_cubical_complex(shape):
        
    dim = len(shape)
    
    comp = Complex(dim)
    
    if dim == 2:
        nrows = shape[0]
        ncols = shape[1]
                
        for i in range(nrows*ncols):
            comp.add_cell(0, {})
        
        for i in range(nrows):
            for j in range(ncols-1):
                comp.add_cell(1, set([ncols*i + j, ncols*i + j+1]))

        for i in range(nrows-1):
            for j in range(ncols):
                comp.add_cell(1, set([ncols*i + j, ncols*(i+1) + j]))
        
        for i in range(nrows-1):
            for j in range(ncols-1):
                comp.add_cell(2, set([nrows*ncols + (ncols-1)*i + j, 
                             nrows*ncols + (ncols-1)*nrows + ncols*i + j+1, 
                             nrows*ncols + (ncols-1)*(i+1) + j, 
                             nrows*ncols + (ncols-1)*nrows + ncols*i + j]))         
        
    
        
    return comp
   
def get_lex_val(c, I, dims, F):
    
    lex_val = set()
        
    Q = co.deque()
    Q.append(c)

    while len(Q) > 0:
        j = Q.popleft()

        if dims[j] == 0:
            lex_val.add(F[j])

        Q.extend(I[j])
        
    return sorted(lex_val, reverse=True)
        
    
def construct_discrete_gradient(comp, F):
    
    def get_lower_star(x):
        
        fx = lex_order[x][0]
        
        lstar = set()
        
        Q = co.deque()
        Q.append(x)
        
        while len(Q) > 0:
            i = Q.popleft()
            
            if all(k <= fx for k in lex_order[i]):
                lstar.add(i)
                Q.extend(comp.cofacets[i])
        
        return lstar
        
    # get faces from a set of cells (not including self)
    def get_faces(alpha, cells):

        faces = set()
        
        Q = co.deque()
        Q.append(alpha)
        
        while len(Q) > 0:
            i = Q.popleft()
            
            if i in cells or i == alpha:
                Q.extend(comp.facets[i])
                faces.add(i)
               
        faces.discard(alpha)
            
        return faces
        
       
    # get cofaces from a set of cells (not including self)
    def get_cofaces(alpha, cells):
        
        cofaces = set()
        
        Q = co.deque()
        Q.append(alpha)
        
        while len(Q) > 0:
            i = Q.popleft()
                        
            if i in cells or i == alpha:
                Q.extend(comp.cofacets[i])
                cofaces.add(i)
               
        cofaces.discard(alpha)
                
        return cofaces
        
    
    if comp.ordered:
        lex_order = [get_lex_val(i, comp.facets, comp.dims, F) for i in comp.get_cells()]
    else:
        lex_order = {i:get_lex_val(i, comp.facets, comp.dims, F) for i in comp.get_cells()}
        
    V = {}
        
    # iterate through each vertex
    for x in range(len(comp.facets)):
        if comp.dims[x] > 0:
            break
                       
        # get lower star of vertex
        lstar = get_lower_star(x)
                
        # if the vertex is the only cell in the lower star then add new critical cell
        if len(lstar) == 1:
            V[x] = x
         
        # otherwise process lower star
        else:
            
            # minimum one cell
            delta = None
            
            # all other one cells
            one_cells = set()
            
            for i in lstar:
                if comp.dims[i] == 1:
                    one_cells.add(i)
                    
                    if delta is None or lex_order[i] < lex_order[delta]:
                        delta = i
                    
            one_cells.discard(delta)
                        
            # add (x, delta) to gradient list and remove both from lstar
            V[x] = delta
            lstar.discard(x)
            lstar.discard(delta)
            
            # priority queue for one cells
            PQone = queue.PriorityQueue()
            
            for i in one_cells:
                PQone.put((lex_order[i], i))
            
            # priority queue for higher dimension cells
            PQother = queue.PriorityQueue()
                        
            for alpha in get_cofaces(delta, lstar):
                
                # if cell is a coface of delta and only has one unpaired face that is in lstar
                if len(get_faces(alpha, lstar)) == 1:
                    PQother.put((lex_order[alpha], alpha))
                
            while not PQone.empty() or not PQother.empty():
                while not PQother.empty():
                    (lex, alpha) = PQother.get()
                    
                    faces = get_faces(alpha, lstar)
                    faces.discard(alpha)
                    # if no unpaired faces
                    if len(faces) == 0:
                        PQone.put((lex_order[alpha], alpha))
                        continue
                        
                    # add pair
                    pair_alpha = faces.pop()
                    V[pair_alpha] = alpha

                    lstar.discard(alpha)
                    lstar.discard(pair_alpha)
                    
                    for beta in get_cofaces(alpha, lstar):
                        if len(get_faces(beta, lstar)) == 1:
                            PQother.put((lex_order[beta], beta))
                            
                    for beta in get_cofaces(pair_alpha, lstar):
                        if len(get_faces(beta, lstar)) == 1:
                            PQother.put((lex_order[beta], beta))
                        
                        
                while not PQone.empty():
                    (lex, gamma) = PQone.get()
                    if gamma not in lstar:
                        continue
                       
                    V[gamma] = gamma
                    lstar.discard(gamma)
                    
                    for alpha in get_cofaces(gamma, lstar):
                        if len(get_faces(alpha, lstar)) == 1:
                            PQother.put((lex_order[alpha], alpha))
    
    return V
    
def reverse_discrete_gradient(V):
    
    coV = {}
    for x in V:
        coV[V[x]] = x
        
    return coV
    

def traverse_flow(s, V, I, coordinated=False):
    
    if coordinated:
        nrIn = {}
        k = {}
        for (a, b, c) in traverse_flow(s, V, I):
            if c not in nrIn:
                nrIn[c] = 1
            else:
                nrIn[c] += 1
            
    Q = co.deque([s])
    seen = {s}
    
    while len(Q) > 0:
        a = Q.popleft()
        for b in I[a]:
            if b in V and V[b] != a:
                c = V[b]
                yield (a, b, c)
                if c not in seen and c != b:
                    if coordinated:
                        if c not in k:
                            k[c] = 1
                        else:
                            k[c] += 1
                            
                        if k[c] != nrIn[c]:
                            continue
                            
                    Q.append(c)
                    seen.add(c)
                
                
            
def find_basins(coV, coI, dims, d):
        
    basins = {}
    
    for s in coV:
        if coV[s] == s and dims[s] == d:
            basins[s] = set([s])
            for (a, b, c) in traverse_flow(s, coV, coI):
                if dims[c] == d:
                    basins[s].add(c)
                
                
    return basins
    
# returns all critical cells that are boundaries of s in the morse complex
# counts = 1: there is a single path and (s, c) are a cancellable close pair
# counts = 2: there are an even number of paths and c is not a boundary
# counts = 3: there are an odd number of paths and (s, c) are not a cancellable close pair
def calc_morse_boundary(s, V, I):
    
    counts = {s:1}
    boundary = set()
    
    for (a, b, c) in traverse_flow(s, V, I, True):
        n = counts[a]
        if c in counts:
            n += counts[c]
            
        if n <= 3:
            counts[c] = n
        else:
            counts[c] = n % 2 + 2
        
        if b == c:
            boundary.add(c)
            
    for c in boundary:
        yield (c, counts[c])

        
# find path from s to t
def find_connections(s, t, V, coV, I, coI):
    
    active = set([t])
    for (a, b, c) in traverse_flow(t, coV, coI):
        active.add(c)
        
    for (a, b, c) in traverse_flow(s, V, I):
        if b in active:
            yield (a, b, c)

            
def construct_morse_complex(V, I, comp):
    
    mcomp = Complex(comp.dim, ordered=False)
    
    for s in V:
        if V[s] == s:
            facets = {c for (c, k) in calc_morse_boundary(s, V, I) if k % 2 == 1}
            mcomp.add_cell(comp.dims[s], facets, s)
            
    return mcomp
            
def construct_filtration(I, dims, comp, F):
    
    Q = queue.PriorityQueue()
    for i in comp.get_cells():
        Q.put((get_lex_val(i, I, dims, F), i))
        
    while not Q.empty():
        (lex, c) = Q.get()
        yield (max(lex), c)
    
def find_morse_skeleton(V, coV, I, coI, dims, d):
    
    skeleton = set()
    for s in V:
        if V[s] == s and dims[s] == d:
            for (t, count) in calc_morse_boundary(s, V, I):
                for (a, b, c) in find_connections(s, t, V, coV, I, coI):
                    skeleton.add(b)
                    
                    
    return skeleton
                
                
def compute_persistence_pairs(comp, filtration, show_zero=False):
    
    columns = []
    
    cell_to_col = {}
    col_to_cell = []
    
    cell_index = {}
    
    pi = 0
    hbins = []
       
    icol = 0
    
    for h, ci in filtration:
        
        
        
        if len(hbins) == 0: 
            hbins.append(h)
        elif h != hbins[-1]:
            hbins.append(h)
            pi += 1    
        
        d = comp.dims[ci]
        
        # print(pi, h, ci, d)
        
        col = []
        for cj in comp.facets[ci]:
            if cj not in cell_to_col:
                print(ci, cj, d, h, comp.facets[ci])
            col.append(cell_to_col[cj])
            
        columns.append((d, sorted(col)))
        
        
        cell_to_col[ci] = icol
        icol += 1
        col_to_cell.append(ci)

        cell_index[ci] = (pi, h)             
                        
    boundary_matrix = phat.boundary_matrix(columns=columns)
    
    alive = set(range(len(columns)))
        
    pairs = []
    for (i, j) in boundary_matrix.compute_persistence_pairs():
        pairs.append((col_to_cell[i], col_to_cell[j]))
        alive.discard(i)
        alive.discard(j)
        
    for i in alive:
        pairs.append((col_to_cell[i], None))
            
    
    return (pairs, cell_index)

    


# def compute_graph_segmentation(G, heights, euclidean=False, positions=None):
    
# #     print("Finding unique heights")
        
#     # basin 0 is the watersheds
#     basins = np.full(G.order(), -1, int)
#     revisit = set()
    
#     current_label = 0
    
#     min_dist = 1.0
         
#     print("Sorting heights")
        
#     # iterate through each level set
#     # for ih in range(len(level_sets)):
    
#     arg_heights = np.argsort(heights)
    
#     ih = 0
#     h = heights[arg_heights[ih]]
    
#     print("Iterating through levels")
    
#     while ih < len(arg_heights):
        
#         unvisited = revisit.copy()
                
#         while True:
            
#             if ih % 1000 == 0:
#                 print("Level:", ih, "/", len(arg_heights), h, len(revisit))
            
#             vi = arg_heights[ih]
            
#             h = heights[vi]
#             unvisited.add(vi)
            
#             ih += 1
                        
#             if ih >= len(arg_heights) or h != heights[arg_heights[ih]]:
#                 break
                    
#         dist = {}
        
#         if euclidean:
#             Q = queue.PriorityQueue()
#             closest = {}
#         else:
#             Q = queue.Queue()            
        
            
#         # find neighbors in lower level sets
#         counter = 0
#         for vi in unvisited:
            
#             for nbr in G[vi]:
#                 if basins[nbr] > 0:
                                        
#                     if euclidean:
#                         new_dist = la.norm(positions[vi] - positions[nbr])                    
#                     else:
#                         new_dist = 1.0
                    
#                     if vi not in dist or new_dist < dist[vi]:
#                         dist[vi] = new_dist
#                         Q.put([new_dist, counter, vi])
#                         counter += 1
                        
#                         if euclidean:
#                             closest[vi] = nbr
                            
     
                                    
#         # breadth first search through connected part of level set
        
#         while not Q.empty():
                        
#             [current_dist, counter, vi] = Q.get()    
            
#             if vi not in unvisited:
#                 continue
            
#             for nbr in G[vi]:
                
#                 if euclidean:
#                     new_dist = la.norm(positions[nbr] - positions[closest[vi]])
#                 else:
#                     new_dist = current_dist + 1.0   
                        
#                 # neighbor has already been assigned to a basin
#                 # 1. has been assigned a basin in a previous level sweep (not in dist)
#                 # 2. assigned a basin because it has a smaller distance (dist[nbr] < current_dist)
#                 if basins[nbr] > 0 and (nbr not in dist or dist[nbr] < current_dist):
                    
#                     # haven't visited this node yet
#                     if vi in unvisited:
#                         # set to neighbor's basin and mark as visited
#                         basins[vi] = basins[nbr]
#                         unvisited.discard(vi)
                        
#                         # if was a watershed then remove from watershed
#                         revisit.discard(vi)
                        
#                     # have already visited but the neighbor is in a different basin
#                     # if already in watershed, this does nothing
#                     elif basins[vi] != basins[nbr]:
#                         basins[vi] = 0
                        
#                         if dist[vi] > min_dist:
#                             revisit.add(vi)
                
#                 # neighbor has not been visited and not been given a distance
#                 elif nbr in unvisited and (nbr not in dist or new_dist < dist[nbr]):
#                     # append to queue
                    
#                     if vi == 5535 and nbr == 5611:
#                         print("hello")
                    
                    
#                     dist[nbr] = new_dist
#                     Q.put([new_dist, counter, nbr])
#                     counter += 1
                    
#                     if euclidean:
#                         closest[nbr] = closest[vi]
                        
#         # look for components only connected via watershed and assign them to watershed
#         counter = 0
#         for vi in unvisited:
#             for nbr in G[vi]:
#                 if basins[nbr] == 0:
#                     Q.put([counter,vi])
#                     counter += 1
                    
#         while not Q.empty():
#             (d, vi) = Q.get()
#             basins[vi] = 0
#             revisit.add(vi)
#             unvisited.discard(vi)
            
#             for nbr in G[vi]:
#                 if nbr in unvisited:
#                     Q.put([counter, nbr])
#                     counter += 1
                    
                    
        
        
#         # breadth first search through separate components corresponding to new minima
#         counter = 0
#         while len(unvisited) > 0:
                        
#             vi = unvisited.pop()
                        
#             for nbr in G[vi]:
#                 if heights[nbr] < heights[vi]:
#                     print(vi, nbr)
#                     print(heights[vi], heights[nbr])
#                     print(basins[vi], basins[nbr])
            
            
#             current_label += 1
#             basins[vi] = current_label
            
#             Q.put([counter,vi])
#             counter += 1
            
            
                    
#             while not Q.empty():
#                 (d, vj) = Q.get()
                
#                 for nbr in G[vj]:
#                     if nbr in unvisited:
#                         Q.put([counter, nbr])
#                         counter += 1
#                         basins[nbr] = current_label
#                         unvisited.discard(nbr)

    
#     segments = {i:set() for i in range(current_label+1)}
#     for i in range(G.order()):
#         segments[basins[i]].add(i)
                
#     return segments
    
    
# def construct_mesh_graph(nrows, ncols, triangulate=True, get_pos=False):
#     G = nx.Graph()
    
#     pos = {}
    
#     for i in range(nrows):
#         for j in range(ncols):
#             G.add_node(ncols*i + j)
#             if get_pos:
#                 pos[ncols*i + j] = np.array([j, i])
            
            
#     for i in range(nrows):
#         for j in range(ncols-1):
#             G.add_edge(ncols*i + j, ncols*i + j+1)

#     for i in range(nrows-1):
#         for j in range(ncols):
#             G.add_edge(ncols*i + j, ncols*(i+1) + j)

            
#     if triangulate:
        
#         for i in range(nrows-1):
#             for j in range(ncols-1):
#                 if i % 2 == j % 2:
#                     G.add_edge(ncols*i + j, ncols*(i+1) + j+1)
#                 else:
#                     G.add_edge(ncols*(i+1) + j, ncols*i + j+1)
    
#     if get_pos:
#         return (G, pos)
#     else:
#         return G
    
# def construct_mesh_complex(nrows, ncols, triangulate=True, compactify=False):
    
#     nverts = nrows*ncols
#     if compactify:
#         nverts += 1
#     comp = CellComplex(2)

#     for i in range(nverts):
#         comp.add(0, i, [])
    
#     iedge = 0
#     for i in range(nrows):
#         for j in range(ncols-1):
#             comp.add(1, iedge, [ncols*i + j, ncols*i + j+1])
#             iedge += 1

#     for i in range(nrows-1):
#         for j in range(ncols):
#             comp.add(1, iedge, [ncols*i + j, ncols*(i+1) + j])
#             iedge += 1
            
#     if triangulate:
#         for i in range(nrows-1):
#             for j in range(ncols-1):
#                 if i % 2 == j % 2:
#                     comp.add(1, iedge, [ncols*i + j, ncols*(i+1) + j+1])
#                 else:
#                     comp.add(1, iedge, [ncols*(i+1) + j, ncols*i + j+1])
#                 iedge += 1
                                
#         iface = 0
#         for i in range(nrows-1):
#             for j in range(ncols-1):
#                 if i % 2 == j % 2:
#                     comp.add(2, iface, [(ncols-1)*i + j, 
#                                  (ncols-1)*nrows + ncols*i + j+1,
#                                  (ncols-1)*nrows + ncols*(nrows-1) + (ncols-1)*i + j])
                    
#                     iface += 1
                    
#                     comp.add(2, iface, [(ncols-1)*(i+1) + j, 
#                                  (ncols-1)*nrows + ncols*i + j,
#                                  (ncols-1)*nrows + ncols*(nrows-1) + (ncols-1)*i + j])
                    
#                     iface += 1
                    
#                 else:
#                     comp.add(2, iface, [(ncols-1)*i + j, 
#                                  (ncols-1)*nrows + ncols*(nrows-1) + (ncols-1)*i + j,
#                                  (ncols-1)*nrows + ncols*i + j])
                
#                     iface += 1
                    
#                     comp.add(2, iface, [(ncols-1)*(i+1) + j, 
#                                  (ncols-1)*nrows + ncols*(nrows-1) + (ncols-1)*i + j,
#                                  (ncols-1)*nrows + ncols*i + j+1])
                
#                     iface += 1
      

#     else:
#         iface = 0
#         for i in range(nrows-1):
#             for j in range(ncols-1):
#                 comp.add(2, iface, [(ncols-1)*i + j, 
#                              (ncols-1)*nrows + ncols*i + j+1, 
#                              (ncols-1)*(i+1) + j, 
#                              (ncols-1)*nrows + ncols*i + j]) 
#                 iface += 1
                             
                                    
                                    
                                    
#                              # (ncols-1)*(i+1) + j, 
#                              # (ncols-1)*nrows + ncols*i + j])
#                 iface += 1
            
# #     if compactify:
# #         for j in range(ncols-1):
# #             comp.add(1, iedge, [j, nrows*ncols])
# #             iedge += 1
        
# #         for i in range(nrows-1):
# #             comp.add(1, iedge, [ncols*i + ncols-1, ncols*nrows])
# #             iedge += 1
            
# #         for j in range(ncols-1, 0, -1):
# #             comp.add(1, iedge, [ncols*(nrows-1) + j, ncols*nrows])
# #             iedge += 1
            
# #         for i in range(nrows-1, 0, -1):
# #             comp.add(1, iedge, [ncols*i, ncols*nrows])
# #             iedge += 1
        
# #         for j in range(ncols-1):
# #             comp.add(2, iface, [j, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + j, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + j+1])
# #             iface += 1
            
# #         for i in range(nrows-1):
# #             comp.add(2, iface, [(ncols-1)*nrows + ncols-1 + ncols*i, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + i, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + i+1])
# #             iface += 1
        
# #         for j in range(ncols-1):
# #             comp.add(2, iface, [(ncols-1)*nrows - 1 - j , 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + j, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + j+1])
        
# #         iface += 1
# #         for i in range(nrows-2):
# #             comp.add(2, iface, [(ncols-1)*nrows + ncols*(nrows-1) - ncols - ncols*i, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + i, 
# #                          (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + i+1])
# #             iface += 1
        
# #         comp.add(2, iface, [(ncols-1)*nrows + ncols*(nrows-1) - ncols - ncols*(nrows-2), 
# #                      (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + nrows-2, 
# #                      (ncols-1)*nrows + ncols*(nrows-1)])
            
#     return comp
    
# def compute_graph_dilation(G, sources, euclidean=False, positions=None):
#     dist = np.full(G.order(), -1, float)
    
#     if euclidean:
#         Q = queue.PriorityQueue()
#         closest = np.full(G.order(), -1, int)
#     else:
#         Q = queue.Queue()      
    
#     unvisited = set(np.arange(G.order()))
    
#     for i in sources:
#         dist[i] = 0
#         Q.put((0, i))
        
#         if euclidean:
#             closest[i] = i
        
        
#     while not Q.empty():
        
#         (current_dist, vi) = Q.get()
#         if vi not in unvisited:
#             continue
            
#         unvisited.discard(vi)
                
#         for nbr in G[vi]:
            
#             if nbr in unvisited:
                
#                 if euclidean:
#                     new_dist = la.norm(positions[nbr] - positions[closest[vi]])
#                 else:
#                     new_dist = current_dist + 1
                                                                
#                 if nbr in unvisited and (dist[nbr] == -1 or new_dist < dist[nbr]):
#                     dist[nbr] = new_dist 
#                     Q.put((new_dist, nbr))   
                    
#                     if euclidean:
#                         closest[nbr] = closest[vi]
                                    
#     return dist

# def laplace_dilation_filtration(mat, show=False):
    
#     image = np.zeros_like(mat, float)

#     i = 0
#     while True:

#         new_image = np.copy(mat.astype(float))

#         # bulk
#         new_image[1:mat.shape[0]-1, 1:mat.shape[1]-1] += (image[0:mat.shape[0]-2, 1:mat.shape[1]-1] 
#                                                           + image[1:mat.shape[0]-1, 0:mat.shape[1]-2] 
#                                                           + image[2:mat.shape[0], 1:mat.shape[1]-1] 
#                                                           + image[1:mat.shape[0]-1, 2:mat.shape[1]])

#         # top row
#         new_image[0, 1:mat.shape[1]-1] += (image[0, 0:mat.shape[1]-2] 
#                                           + image[1, 1:mat.shape[1]-1] 
#                                           + image[0, 2:mat.shape[1]])

#         # bottom row
#         new_image[mat.shape[0]-1, 1:mat.shape[1]-1] += (image[mat.shape[0]-2, 1:mat.shape[1]-1] 
#                                                           + image[mat.shape[0]-1, 0:mat.shape[1]-2] 
#                                                           + image[mat.shape[0]-1, 2:mat.shape[1]])

#         # left column
#         new_image[1:mat.shape[0]-1, 0] += (image[0:mat.shape[0]-2, 0] 
#                                                           + image[2:mat.shape[0], 0] 
#                                                           + image[1:mat.shape[0]-1, 1])

#         # right column
#         new_image[1:mat.shape[0]-1, mat.shape[1]-1] += (image[0:mat.shape[0]-2, mat.shape[1]-1] 
#                                                           + image[1:mat.shape[0]-1, mat.shape[1]-2] 
#                                                           + image[2:mat.shape[0], mat.shape[1]-1])

#         # upper left corner
#         new_image[0, 0] += (image[1, 0] 
#                               + image[0, 1])

#         # upper right corner
#         new_image[0, mat.shape[1]-1] += (image[0, mat.shape[1]-2] 
#                                           + image[1, mat.shape[1]-1])

#         # bottorm left corner
#         new_image[mat.shape[0]-1, 0] += (image[mat.shape[0]-2, 0] 
#                                             + image[mat.shape[0]-1, 1])

#         # bottom right corner
#         new_image[mat.shape[0]-1, mat.shape[1]-1] += (image[mat.shape[0]-2, mat.shape[1]-1] 
#                                                           + image[mat.shape[0]-1, mat.shape[1]-2])

#         new_image /= 6.0

        

#         max_delta = np.max(np.abs((image-new_image) / np.where(image==0, 1e-16, np.abs(image))))


#         if show and i %10 == 0:
#             print(i, max_delta)
            
#             # norm = mcolors.LogNorm(vmin=1e-16, vmax=1.0)
#             # cmap = mpl.cm.GnBu
            
#             norm = mcolors.SymLogNorm(vmin=-1.0, vmax=1.0, linthresh=1e-32)
#             # cmap = mpl.cm.RdBu_r
#             cmap = mpl.cm.RdYlGn
            
#             fig, ax = plt.subplots(figsize=(8,8))
#             im = ax.imshow(new_image, cmap=cmap, norm=norm)
#             # plt.colorbar(im)
#             # ax.axis('off')
#             plt.show()

        
#         if max_delta <= 1e-4:
#             image = new_image
#             break

#         image = new_image
#         i += 1
        
#     heights = list(np.sort(image.flatten()))
#     pixels = list(np.argsort(image.flatten()))
    
#     pixels.append(mat.shape[0]*mat.shape[1])
#     heights.append(heights[-1]+1)
    
#     return (pixels, heights)

        