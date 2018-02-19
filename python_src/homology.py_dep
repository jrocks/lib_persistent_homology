import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import numpy.linalg as la
import scipy.sparse as sparse
import scipy.optimize as opt
import itertools as it
import phat
import queue
import collections as co
import persist




# class CellComplex:
    
    
#     def __init__(self, dim, reserve=64, regular=True, oriented=False, ordered=True):
        
#         self.dim = dim
#         self.ncells = 0
#         self.nfacets = 0
        
#         self.ordered = ordered
#         self.regular = regular
#         self.oriented = oriented
        
        
#         if not self.ordered:
#             self.cell_to_index = {}
            
#         self.dims = np.empty(reserve, int)
            
#         self.facet_ind = np.empty(reserve+1, int)
#         self.facet_ind[0] = 0
#         self.facets = np.empty(reserve, int)
   
#         if not self.regular or self.oriented:
#             self.coeffs = np.empty(reserve, int)
            
        
#     # label only relevant if not ordered
#     # coeffs only relevant if not regular or is oriented
#     def add_cell(self, dim, facets, label=None, coeffs=None):
        
#         # check if arrays need to be resized
#         if self.ncells + 1 > len(self.dims):
#             self.dims.resize(2*self.ncells)
#             self.facet_ind.resize(2*self.ncells+1)
            
#         if self.nfacets + len(facets) > len(self.facets):
#             self.facets.resize(2*self.nfacets)
#             if not self.regular or self.oriented:
#                 self.coeffs.resize(2*self.nfacets)
        
#         if not self.ordered:
#             self.cell_to_index[label] = self.ncells
                
#         self.dims[self.ncells] = dim
#         self.facet_ind[self.ncells+1] = self.facet_ind[self.ncells] + len(facets)
        
#         self.facets[self.facet_ind[self.ncells]:self.facet_ind[self.ncells+1]] = facets
#         if not self.regular or self.oriented:
#             self.coeffs[self.facet_ind[self.ncells]:self.facet_ind[self.ncells+1]] = coeffs
        
#         self.ncells += 1
#         self.nfacets += len(facets)
        
        
#     def make_compressed(self):
#         self.dims.resize(self.ncells)
#         self.facet_ind.resize(self.ncells+1)
#         self.facets.resize(self.nfacets)
#         if not self.regular or self.oriented:
#             self.coeffs.resize(self.nfacets)
        
    
#     def get_dim(self, alpha):
#         if self.ordered:
#             return self.dims[alpha]
#         else:
#             return self.dims[self.cell_to_index[alpha]]
    
#     def get_facets(self, alpha):
#         if self.ordered:
#             index = alpha
#         else:
#             index = self.cell_to_index[alpha]
            
#         return self.facets[self.facet_ind[index]:self.facet_ind[index+1]]
        
#     def get_coeffs(self, alpha):
#         if self.ordered:
#             index = alpha
#         else:
#             index = self.cell_to_index[alpha]
            
#         facets = self.facets[self.facet_ind[index]:self.facet_ind[index+1]]
#         coeffs = self.coeffs[self.facet_ind[index]:self.facet_ind[index+1]]
#         return dict(zip(facets, coeffs))        

    
#     def get_cofacets(self, alpha):
#         if self.ordered:
#             index = alpha
#         else:
#             index = self.cell_to_index[alpha]
            
#         return self.cofacets[self.cofacet_ind[index]:self.cofacet_ind[index+1]]
    
            
#     def get_cells(self):
        
#         if self.ordered:
#             for i in range(self.ncells):
#                 yield i     
#         else:
#             for i in self.cell_to_index:
#                 yield i
        
#     def construct_cofacets(self):
        
#         self.cofacet_ind = np.empty(self.ncells+1, int)
#         self.cofacet_ind[0] = 0
#         self.cofacets = np.empty(self.nfacets, int)
        
#         cell_list = [[] for i in range(self.ncells)]
#         if self.ordered:
#             for i in range(self.ncells):
#                 for j in self.facets[self.facet_ind[i]:self.facet_ind[i+1]]:
#                         cell_list[j].append(i)
#         else:
#             for i in self.cell_to_index:
#                 for j in self.facets[self.facet_ind[self.cell_to_index[i]]:self.facet_ind[self.cell_to_index[i]+1]]:
#                     cell_list[self.cell_to_index[j]].append(i)

#         for i in range(self.ncells):
#             self.cofacet_ind[i+1] = self.cofacet_ind[i] + len(cell_list[i])
#             self.cofacets[self.cofacet_ind[i]:self.cofacet_ind[i+1]] = cell_list[i]
    
            
def construct_cubical_complex(shape, oriented=False, dual=False):
        
    return persist.construct_cubical_complex(shape, oriented, dual)
        
#     dim = len(shape)
    
#     # comp = CellComplex(dim, reserve=64, regular=True, oriented=oriented, ordered=True)
    
#     comp = persist.CellComplex(dim, True, oriented, True)
    
#     if dim == 2 and not dual:
#         nrows = shape[0]
#         ncols = shape[1]
          
#         # vertices
#         for i in range(nrows*ncols):
            
#             if oriented:
#                 coeffs = []
#             else:
#                 coeffs = []
#             comp.add_cell(0, [], -1, coeffs)
        
#         # horizontal edges
#         for i in range(nrows):
#             for j in range(ncols-1):
#                 if oriented:
#                     coeffs = [-1, 1]
#                 else:
#                     coeffs = []
#                 comp.add_cell(1, [ncols*i + j, ncols*i + j+1], -1, coeffs)

#         # vertical edges
#         for i in range(nrows-1):
#             for j in range(ncols):
#                 if oriented:
#                     coeffs = [-1, 1]
#                 else:
#                     coeffs = []
#                 comp.add_cell(1, [ncols*i + j, ncols*(i+1) + j], -1, coeffs)
        
#         # faces
#         for i in range(nrows-1):
#             for j in range(ncols-1):
#                 if oriented:
#                     coeffs = [1, 1, -1, -1]
#                 else:
#                     coeffs = []
#                 comp.add_cell(2, [nrows*ncols + (ncols-1)*i + j, 
#                              nrows*ncols + (ncols-1)*nrows + ncols*i + j+1, 
#                              nrows*ncols + (ncols-1)*(i+1) + j, 
#                              nrows*ncols + (ncols-1)*nrows + ncols*i + j], -1, coeffs)
                
#     elif dim == 2 and dual:
#         nrows = shape[0]+1
#         ncols = shape[1]+1
        
#         nverts = nrows*ncols
#         nhedges = nrows*(ncols-1)
#         nvedges = (nrows-1)*ncols
#         nfaces = (nrows-1)*(ncols-1)
                
#         # faces
#         for i in range(nrows-1):
#             for j in range(ncols-1):
#                 if oriented:
#                     coeffs = [1, 1, -1, -1]
#                 else:
#                     coeffs = []
#                 comp.add_cell(2, [nfaces + (ncols-1)*i + j, 
#                              nfaces + nhedges + ncols*i + j+1, 
#                              nfaces + (ncols-1)*(i+1) + j, 
#                              nfaces + nhedges + ncols*i + j], -1, coeffs)
                
                
#         # horizontal edges
#         for i in range(nrows):
#             for j in range(ncols-1):
#                 if oriented:
#                     coeffs = [-1, 1]
#                 else:
#                     coeffs = []
#                 comp.add_cell(1, [nfaces + nhedges + nvedges + ncols*i + j, 
#                                   nfaces + nhedges + nvedges + ncols*i + j+1], -1, coeffs)

#         # vertical edges
#         for i in range(nrows-1):
#             for j in range(ncols):
#                 if oriented:
#                     coeffs = [-1, 1]
#                 else:
#                     coeffs = []
#                 comp.add_cell(1, [nfaces + nhedges + nvedges + ncols*i + j, 
#                                   nfaces + nhedges + nvedges + ncols*(i+1) + j], -1, coeffs)
            
            
#         # vertices
#         for i in range(nrows*ncols):
            
#             if oriented:
#                 coeffs = []
#             else:
#                 coeffs = []
#             comp.add_cell(0, [], -1, coeffs)
        
        
        
#     comp.make_compressed() 
    
        
#     return comp
   
def check_boundary_op(comp):
    
    return persist.check_boundary_op(comp);
    
#     valid = True
#     for i in comp.get_cells():
#         if comp.get_dim(i) > 0:
            
#             if not comp.regular or comp.oriented:
#                 sub_faces = {}
#                 for j in comp.get_facets(i):
#                     for k in comp.get_facets(j):
#                         sub_faces[k] = sub_faces.get(k, 0) + comp.get_coeffs(i)[j] * comp.get_coeffs(j)[k]
                        
#                 for j in sub_faces:
#                     if (comp.oriented and sub_faces[j] != 0) or (not comp.oriented and sub_faces[j] % 2 != 0):
#                         print("Error:", i, j, sub_faces[j])
#                         valid = False
#             else:
#                 sub_faces = set()
#                 for j in comp.get_facets(i):
#                     sub_faces ^= set(comp.get_facets(j))
                   
#                 if len(sub_faces) != 0:
#                     print("Error:", i, sub_faces)
#                     valid = False
                    
#     return valid
            
    
    
def construct_vertex_filtration_order(comp, vertex_time, euclidean=False, positions=None, dual=False):
    
    return persist.construct_vertex_filtration_order(comp, vertex_time, dual)
    
#     vertex_order = np.zeros(len(vertex_time), int)
        
#     submerged = np.zeros(len(vertex_time), bool)
        
#     vertex_argsort = list(np.argsort(vertex_time)[::-1])
        
#     if dual:
#         coI = comp.get_facets
#         I = comp.get_cofacets
#     else:
#         coI = comp.get_cofacets
#         I = comp.get_facets
        
#     ti = 0
#     while len(vertex_argsort) > 0:
#         unvisited = set()
            
#         # find all vertices in level set
#         while True:
 
#             vi = vertex_argsort.pop()
                    
#             t = vertex_time[vi]
#             unvisited.add(vi)
                        
#             if len(vertex_argsort) == 0 or t != vertex_time[vertex_argsort[-1]]:
#                 break
          
        
#         # if only a single vertex, then take a shortcut
#         if len(unvisited) == 1:
#             vi = unvisited.pop()
#             submerged[vi] = True
#             vertex_order[vi] = ti
#             ti += 1
                        
#             continue
            
        
#         # set up queue
#         dist = {}
#         if euclidean:
#             Q = queue.PriorityQueue()
#             closest = {}
#         else:
#             Q = queue.Queue() 
            
#         # put vertices in queue if they have a previously submerged neighbor
#         for a in unvisited:
#             for b in coI(a):
#                 for c in I(b):
#                     if submerged[c]:
                        
#                         if euclidean:
#                             new_dist = la.norm(positions[a] - positions[c])                    
#                         else:
#                             new_dist = 1.0
                            
#                         if a not in dist or new_dist < dist[a]:
#                             dist[a] = new_dist
#                             Q.put([new_dist, a])
                            
#                             if euclidean:
#                                 closest[a] = c
             
                            
#         # perform BFS on level set
#         while not Q.empty():
                        
#             [current_dist, a] = Q.get()    
                                        
#             if a not in unvisited:
#                 continue
                    
#             unvisited.remove(a)
#             submerged[a] = True
#             vertex_order[a] = ti
#             ti += 1
              
#             for b in coI(a):
#                 for c in I(b):
#                     if c in unvisited:
                    
#                         if euclidean:
#                             new_dist = la.norm(positions[c] - positions[closest[a]])
#                         else:
#                             new_dist = current_dist + 1.0  
                          
#                         if c not in dist or new_dist < dist[c]:
#                             dist[c] = new_dist
#                             Q.put([new_dist, c])
                            
#                             if euclidean:
#                                 closest[c] = closest[a]
      
    
#         # search through new segments
#         while len(unvisited) > 0:

#             s = unvisited.pop()
            
#             Q.put([0, s])
            
#             if euclidean:
#                 closest[s] = s
            
#             while not Q.empty():
#                 [current_dist, a] = Q.get()
                                
#                 unvisited.discard(a)
#                 submerged[a] = True
#                 vertex_order[a] = ti
#                 ti += 1
                
#                 for b in coI(a):
#                     for c in I(b):
#                         if c in unvisited:

#                             if euclidean:
#                                 new_dist = la.norm(positions[c] - positions[closest[a]])
#                             else:
#                                 new_dist = current_dist + 1.0  

#                             if c not in dist or new_dist < dist[c]:
#                                 dist[c] = new_dist
#                                 Q.put([new_dist, c])

#                                 if euclidean:
#                                     closest[c] = closest[a]    
#     return vertex_order  
    
    
# def get_star(alpha, I, cell_dim=None, dims=None, insert_order=None):
      
#     if insert_order is not None:
#         STAR = 1
#         star_index = insert_order[alpha][STAR]

#     star = set()

#     Q = co.deque()
#     Q.append(alpha)

#     while len(Q) > 0:
#         a = Q.popleft()

#         if insert_order is None or star_index == insert_order[a][STAR]:
#             Q.extend(I(a))
#             if cell_dim is None or dims(a) == cell_dim:
#                 star.add(a)
                
#     return star
    

# # def get_lex_val(alpha, I, cell_dim, dims, F):
# def get_lex_val(alpha, co, comp, target_dim, F):
    
#     lex_val = set()
        
#     # cells = get_star(alpha, I, cell_dim=cell_dim, dims=dims)
#     cells = persist.get_star(alpha, co, comp, target_dim)
#     for c in cells:
#         lex_val.add((F[c], c))
        
#     return sorted(lex_val, reverse=True)
      
    
# construct map of cells to insertion times
# insertion times are a triple of 
# [insertion time, 
# lower star (upper costar) cell, 
# lexicographic insertion index (index ordering on all cells)]
# where the insertion index is found by computing the lexicographic order of the cells
def construct_time_of_insertion_map(comp, vertex_time, vertex_order, dual=False):
    
    
    return persist.construct_time_of_insertion_map(comp, vertex_order, vertex_time, dual)
    
#     order_to_time = vertex_time[np.argsort(vertex_order)]
                
#     lex = []
    
#     if dual:
#         # Order by:
#         # 1. function value of cell whose upper costar this cell belongs to
#         # 2. cell dimension
#         # 3. lexicographic ordering of function values of highest dimension cofaces from high to low
#         for c in comp.get_cells():
#             # lex_val = get_lex_val(c, comp.get_cofacets, comp.dim, comp.get_dim, vertex_order) 
#             lex_val = get_lex_val(c, False, comp, comp.dim, vertex_order) 
#             lex.append((lex_val[-1][0], comp.get_dim(c), lex_val, c, lex_val[-1][1]))
#     else:
#         # Order by:
#         # 1. function value of vertex whose lower star this cell belongs to
#         # 2. cell dimension
#         # 3. lexicographic ordering of function values of vertices from high to low
#         for c in comp.get_cells():
#             # lex_val = get_lex_val(c, comp.get_facets, 0, comp.get_dim, vertex_order)
#             lex_val = get_lex_val(c, True, comp, 0, vertex_order) 
#             lex.append((lex_val[0][0], comp.get_dim(c), lex_val, c, lex_val[0][1]))
      
#     lex = sorted(lex)
    
#     if comp.ordered:
#         # insert_order = [None for c in comp.get_cells()]
#         insert_order = np.zeros([comp.ncells, 3], dtype=np.float64)
#     else:
#         insert_order = {}
    
#     for i, (star_val, d, lex_val, c, star) in enumerate(lex):
#         insert_order[c] = (order_to_time[star_val], star, i)  

#     return insert_order
 
    
    
def construct_discrete_gradient(comp, insert_order, dual=False):
    
    
    return persist.construct_discrete_gradient(comp, insert_order, dual)
        
#     STAR = 1
#     LEX = 2     
        
#     V = -np.ones(comp.ncells, int)
#     # V = {}
    
#     # iterate through each vertex (or dual cell)
#     for x in comp.get_cells():
#         if (dual and comp.get_dim(x) != comp.dim) or (not dual and comp.get_dim(x) != 0):
#             continue
    
#         if dual:
#             # get upper costar
#             # star = get_star(x, comp.get_facets, insert_order=insert_order)
#             star = persist.get_lower_star(x, True, comp, -1, insert_order)
#         else:
#             # get lower star
#             # star = get_star(x, comp.get_cofacets, insert_order=insert_order)
#             star = persist.get_lower_star(x, False, comp, -1, insert_order)
            
#         unpaired = {}
            
#         if dual:
#             for alpha in star:
#                 unpaired[alpha] = star & set(comp.get_cofacets(alpha))
#         else:
#             for alpha in star:
#                 unpaired[alpha] = star & set(comp.get_facets(alpha))

#         # print("cell", x)
#         # print(star)
#         # print(unpaired)
        
#         PQone = queue.PriorityQueue()
#         PQzero = queue.PriorityQueue()

#         for alpha in star:
#             if len(unpaired[alpha]) == 0:
#                 if dual:
#                     PQzero.put((-insert_order[alpha][LEX], alpha))
#                 else:
#                     PQzero.put((insert_order[alpha][LEX], alpha))
#                     # print("zero", alpha)
#             elif len(unpaired[alpha]) == 1:
#                 if dual:
#                     PQone.put((-insert_order[alpha][LEX], alpha))
#                 else:
#                     PQone.put((insert_order[alpha][LEX], alpha))
#                     # print("one", alpha)
                
#         while not PQone.empty() or not PQzero.empty():
            
#             while not PQone.empty():
#                 (order, alpha) = PQone.get()
                
#                 if len(unpaired[alpha]) == 0:
#                     PQzero.put((order, alpha))
#                     continue
                
#                 beta = unpaired[alpha].pop()
#                 if dual:
#                     V[alpha] = beta
#                 else:
#                     V[beta] = alpha
                
#                 # print(beta, alpha)
                
#                 star.discard(alpha)
#                 star.discard(beta)
                
#                 del unpaired[alpha]
#                 for gamma in star:
#                     if alpha in unpaired[gamma] or beta in unpaired[gamma]:
#                         unpaired[gamma].discard(alpha)
#                         unpaired[gamma].discard(beta)
                        
#                         if len(unpaired[gamma]) == 1:
#                             if dual:
#                                 PQone.put((-insert_order[gamma][LEX], gamma))
#                             else:
#                                 PQone.put((insert_order[gamma][LEX], gamma))
                        
                    
#             if not PQzero.empty():
#                 (order, alpha) = PQzero.get()
#                 if alpha in star:
#                     V[alpha] = alpha
                    
#                     # print(alpha, alpha)
                    
#                     star.discard(alpha)
                    
#                     del unpaired[alpha]
#                     for gamma in star:
#                         if alpha in unpaired[gamma]:
#                             unpaired[gamma].discard(alpha)

#                             if len(unpaired[gamma]) == 1:
#                                 if dual:
#                                     PQone.put((-insert_order[gamma][LEX], gamma))
#                                 else:
#                                     PQone.put((insert_order[gamma][LEX], gamma))
                                    
          
     # coV = -np.ones(len(V), int)
#     for x in range(len(V)):
#         if V[x] != -1:
#             coV[V[x]] = x
        
        
#     return (V, coV)

    
def reverse_discrete_gradient(V):
    
    return persist.reverse_gradient(V)
    
# #     coV = {}
# #     for x in V:
# #         coV[V[x]] = x
      
#     coV = -np.ones(len(V), int)
#     for x in range(len(V)):
#         if V[x] != -1:
#             coV[V[x]] = x
    
#     return coV
    

def traverse_flow(s, V, I, coordinated=False):
    
    if coordinated:
        nrIn = {}
        k = {}
        for (a, b, c) in traverse_flow(s, V, I):
            nrIn[c] = nrIn.get(c, 0) + 1
            
    Q = co.deque([s])
    
    seen = {s}

    while len(Q) > 0:
        a = Q.popleft()
        
        for b in I(a):    
            
            # if b in V and V[b] != a:
            if V[b] != -1 and V[b] != a:
                c = V[b]     
                    
                yield (a, b, c)    
                
                if c not in seen and c != b:
                    if coordinated:
                        k[c] = k.get(c, 0) + 1
                        if k[c] != nrIn[c]:
                            continue
                            
                    Q.append(c)
                    seen.add(c)
                

# returns all critical cells that are boundaries of s in the morse complex
# return tuples (c, count, mult)
# c is the critical boundary cell
# count is the total number of V-paths to from s to c
# mult is the coefficient of the boundary operator 
def calc_morse_boundary(s, V, I, coeff, oriented=False):
    
    counts = {s:1}
    boundary = set()
    
    if oriented:
        mult = {s:1}
    
    for (a, b, c) in traverse_flow(s, V, I, True):
        # n = counts[a] + counts.get(c, 0)
        # counts[c] = n if n <= 3 else n % 2 + 2
        
        counts[c] = counts[a] + counts.get(c, 0)
        
        if b == c:
            boundary.add(c)
            
            if oriented:
                mult[c] = mult.get(c, 0) + mult[a] * ( -coeff(a)[b] )
            
        elif oriented:
            mult[c] = mult.get(c, 0) + mult[a] * ( -coeff(a)[b] * coeff(c)[b] )
            
    for c in boundary:
        if oriented:
            yield (c, counts[c], mult[c])  
        else:
            yield (c, counts[c], counts[c] % 2) 
    
def construct_morse_complex(V, comp, oriented=False):
    
    return persist.construct_morse_complex(V, comp, oriented)
    
#     mcomp = CellComplex(comp.dim, ordered=False, oriented=oriented, regular=False)
        
#     # for s in V:
#     for s in range(len(V)):
#         if V[s] == s:
#             facets = []
#             coeffs = []
#             for (c, k, m) in calc_morse_boundary(s, V, comp.get_facets, comp.get_coeffs, oriented=oriented):
#                 facets.append(c)
#                 coeffs.append(m)
                
#                 if abs(m) > 1:
#                     print("Large Coefficient", k, m, comp.get_dim(c), comp.get_dim(s))
            
         
#             # facets = [c for (c, k) in calc_morse_boundary(s, V, I) if k % 2 == 1]
                
#             mcomp.add_cell(comp.get_dim(s), facets, s, coeffs)
       
#     mcomp.make_compressed()
#     return mcomp
    
        
# find path from s to t
def find_connections(s, t, V, coV, I, coI):
    
    active = set([t])
    for (a, b, c) in traverse_flow(t, coV, coI):
        print("backwards", a, b, c)
        active.add(c)
        
    for (a, b, c) in traverse_flow(s, V, I):
        print("forwards", a, b, c)
        if b in active:
            yield (a, b, c)
    


def simplify_morse_complex(threshold, V, coV, comp, insert_order):
    
    
    persist.simplify_morse_complex(threshold, V, coV, comp, insert_order, verbose=True);
    return (V, coV)
     
    TIME = 0
    LEX = 2
        
    crit_cells = {s for s in range(len(V)) if V[s] == s}
    crit_order = sorted([(insert_order[s][LEX], s) for s in crit_cells])
        
    n = 0
    while True:
        
        n += 1
                 
        print("Pass:", n)
            
        cancel_pairs = []
        
        # print(sorted(crit_cells))
        print("Critical Cells", len(crit_cells))
              
            
        n_cancel = 0
        for (lex, s) in crit_order:
            
            if comp.get_dim(s) == 0:
                continue
                
            # print("s", s, insert_order[s], comp.get_dim(s))
            
            close_alpha = None
            close_alpha_time = None
            for (c, k, m) in calc_morse_boundary(s, V, comp.get_facets, comp.get_coeffs, oriented=comp.oriented):

                if k == 1:
                    # if close_alpha is None or ((insert_order[c][TIME], insert_order[c][LEX]) > close_alpha_time):
                    #     close_alpha = c
                    #     close_alpha_time = (insert_order[c][TIME], insert_order[c][LEX])
                        
                    if close_alpha is None or insert_order[c][LEX] > close_alpha_time:
                        close_alpha = c
                        close_alpha_time = insert_order[c][LEX]
                    
            if close_alpha is None:
                continue
                
                                
            close_beta = None
            close_beta_time = None
            for (c, k, m) in calc_morse_boundary(close_alpha, coV, comp.get_cofacets, comp.get_coeffs, oriented=comp.oriented):
                if k == 1:
                    # if close_beta is None or ((insert_order[c][TIME], insert_order[c][LEX]) < close_beta_time):
                    #     close_beta = c
                    #     close_beta_time = (insert_order[c][TIME], insert_order[c][LEX])
                        
                    if close_beta is None or insert_order[c][LEX] < close_beta_time:
                        close_beta = c
                        close_beta_time = insert_order[c][LEX]
               
           
            if s == close_beta:
                n_cancel += 1
                if insert_order[s][TIME] - insert_order[close_alpha][TIME] <= threshold:
                    cancel_pairs.append((insert_order[s][TIME] - insert_order[close_alpha][TIME], (close_alpha, close_beta)))
        
        
        cancel_pairs = sorted(cancel_pairs)
        
        print("Cancellable Pairs:", n_cancel)
        
        if len(cancel_pairs) == 0 or cancel_pairs[0][0] > threshold:
            print("Cancelled Pairs:", 0)
            break
          
        for (time, (t, s)) in cancel_pairs:
            print(s, t)
        
            # print(n_cancelled, time, t, s, comp.get_dim(t), comp.get_dim(s))
            reverse_pairs = []
            for (a, b, c) in find_connections(s, t, V, coV, comp.get_facets, comp.get_cofacets):
                reverse_pairs.append((a, b))
                
                
            print(reverse_pairs) 
            return (V, coV)
                
            for (a, b) in reverse_pairs:                
                V[b] = a
                
                coV[a] = b
                
                V[a] = -1
                coV[b] = -1
                
                
                # if a in V:
                #     del V[a]
                # if b in coV:
                #     del coV[b]
                    
            crit_cells.remove(t)
            crit_cells.remove(s)
            
        crit_order = sorted([(insert_order[s][LEX], s) for s in crit_cells])
               
        print("Cancelled Pairs:", len(cancel_pairs))
        
        print("Remaining Critical Cells", len(crit_cells))
            
        
                
    return (V, coV)
       
        
def construct_filtration(mcomp, filtration):
    
    LEX = 2
    
    Q = queue.PriorityQueue()
    for i in range(mcomp.ncells):
        # Q.put((insert_order[i][LEX], i))
        Q.put((filtration.get_filtration_order(mcomp.get_label(i)), i))
        
    while not Q.empty():
        (order, c) = Q.get()
        yield c
    
        
def get_morse_weights(mcomp, V, coV, I, coI):
    
    weights = {}
    for s in range(mcomp.ncells):
        if mcomp.get_dim(s) == 0:
            weights[s] = 1
            continue
            
        cell = set()
        for (t, m) in mcomp.get_facets(s).items():
            for (a, b, c) in find_connections(s, t, V, coV, I, coI):
                cell.add(a)
         
        weights[s] = len(cell)
    
    
    return weights
    
        

def convert_morse_to_real_complex(mfeature, V, I):

    feature = set()
    for s in mfeature:
        feature.add(s)
        for (a, b, c) in traverse_flow(s, V, I):
            if b != c:
                feature.add(c)
    
    return feature


def convert_to_pixels(feature, comp, insert_order, dual=False):
    
    return persist.convert_to_pixels(feature, comp, insert_order, dual)
    
    
#     pixels = set()
    
#     STAR = 1
#     LEX = 2
    
#     for s in feature:
        
#         if comp.get_dim(s) == 0:
#             if not dual:
                
#                 pixels.add(s)
                
#             else:
                
#                 #start replacing snippets here
                
#                 # cofaces = get_star(s, comp.get_cofacets, cell_dim=comp.dim, dims=comp.get_dim)
#                 cofaces = persist.get_star(s, False, comp, comp.dim)
                        
#                 for c in cofaces:
                    
#                     # verts = get_star(c, comp.get_facets, cell_dim=0, dims=comp.get_dim)
#                     verts = persist.get_star(c, True, comp, 0)
                    
#                     ucostar_verts = set()
#                     for v in verts:
#                         if insert_order[c][STAR] == insert_order[v][STAR]:
#                             ucostar_verts.add(v)
                
#                     # no verts in ucostar:
#                     if len(ucostar_verts) == 0:
#                         min_vert = verts.pop()
                        
#                         for v in verts:
#                             if insert_order[v][LEX] < insert_order[min_vert][LEX]:
#                                 min_vert = v
                                
#                         if min_vert in feature:
#                             pixels.add(c)
                     
#                     # all verts in ucostar are in feature
#                     elif ucostar_verts.issubset(feature):
                        
#                         pixels.add(c)
                      
#                     # only some ucostar verts are in feature
#                     else:
                        
#                         min_vert = ucostar_verts.pop()
#                         for v in ucostar_verts:
#                             if insert_order[v][LEX] < insert_order[min_vert][LEX]:
#                                 min_vert = v
                                
#                         if min_vert in feature:
#                             pixels.add(c)
                            
#         elif comp.get_dim(s) == 1:
            
#             if not dual:
                
#                 pixels.update(comp.get_facets(s))
                
#             else:
                        
#                 ucostar = int(insert_order[s][STAR])
                    
#                 pixels.add(ucostar)
                
#         elif comp.get_dim(s) == 2:
            
#             if not dual:
#                 for c in comp.get_facets(s):
#                     pixels.update(comp.get_facets(c))
                    
#             else:
                
#                 pixels.add(s)
                            
                
                        
#     return pixels                
                            
     
def find_basins(mcomp, coV, comp, insert_order, dual=False):
    
    return persist.find_basins(mcomp, coV, comp, insert_order, dual)
    
#     basins = {}
    
#     STAR = 1
    
#     for c in mcomp.get_cells():
#         if mcomp.get_dim(c) == 0:
            
#             if dual:
#                 index = int(insert_order[c][STAR])
#             else:
#                 index = c
                
#             basins[index] = convert_to_pixels(convert_morse_to_real_complex({c}, coV, comp.get_cofacets), 
#                                                 comp, insert_order, dual=dual)
            
          
#     return basins
     
    
def find_morse_skeleton(mcomp, V, comp, d, insert_order, dual=False):
    
    return persist.find_morse_skeleton(mcomp, V, comp, d, insert_order, dual)
    
    
#     skeleton = set()
    
#     # for s in V:
#     for s in range(len(V)):
#         if V[s] == s and mcomp.get_dim(s) == d:
#             feature = convert_to_pixels(convert_morse_to_real_complex({s}, V, comp.get_facets), comp, insert_order, dual=dual)
 
#             skeleton.update(feature)
            
            
#     return skeleton 

        
    
# if (i, j) represents a 1-cycle (connected component), 
# then expand cell i in morse complex up to (but not including) death cell j
# if (i, j) represents a d-cycle (cycle for 2-manifold or void for 3-manifold), 
# then expand cell j in morse complex up to (but not including) birth cell i    


def extract_persistence_feature(i, j, mcomp, comp, V, coV, insert_order):
    
    return persist.extract_persistence_feature(i, j, mcomp, comp, V, coV, insert_order)
    
#     TIME = 0
#     STAR = 1
#     LEX = 2
        
#     if comp.get_dim(i) == 0:
#         I = mcomp.get_facets
#         coI = mcomp.get_cofacets
        
#         barrier_test = lambda b: insert_order[b][TIME] >= insert_order[j][TIME]
        
#         seen = {i}
#         Q = co.deque([i])
    
#     else:
#         I = mcomp.get_cofacets
#         coI = mcomp.get_facets
        
#         barrier_test = lambda b: insert_order[b][TIME] <= insert_order[i][TIME]
        
#         seen = {j}
#         Q = co.deque([j])
    
#     while len(Q) > 0:
#         a = Q.popleft()
#         for b in coI(a):
#             if barrier_test(b):
#                 continue
                
#             for c in I(b):
                
#                 if c != a and c not in seen:
#                     Q.append(c)
#                     seen.add(c)
            
#     if comp.get_dim(i) == 0:
#         feature = convert_morse_to_real_complex(seen, coV, comp.get_cofacets)
#     else:
#         feature = convert_morse_to_real_complex(seen, V, comp.get_facets)
    
#     return feature
        


def get_boundary(cells, comp):
    
    return persist.get_boundary(cells, comp)
    
#     cycle = set()
    
#     if comp.regular:
#         for c in cells:
#             cycle ^= set(comp.get_facets(c))
        
#     return cycle

        
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

        