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
        
        
#     def compress():
#         self.dims.resize(self.ncells)
#         self.facet_ind.resize(self.ncells+1)
#         self.facets.resize(self.nfacets)
#         if not self.regular or self.oriented:
#             self.coeffs.resize(self.nfacets)
        
    
#     def get_dim(self, alpha):
#         if self.ordered:
#             return self.dims[alpha]
#         else:
#             return self.dims[cell_to_index[alpha]]
    
#     def get_facets(self, alpha):
#         if self.ordered:
#             index = alpha
#         else:
#             index = cell_to_index[alpha]
            
#         facets = self.facets[self.facet_ind[index]:self.facet_ind[index+1]]
#         if self.oriented or not self.regular:
#             coeffs = self.coeffs[self.facet_ind[index]:self.facet_ind[index+1]]
#             return dict(zip(facets, coeffs))        
#         else:
#             return facets
    
#     def get_cofacets(self, alpha):
#         if self.ordered:
#             index = alpha
#         else:
#             index = cell_to_index[alpha]
            
#         return self.cofacets[self.cofacet_ind[index]:self.cofacet_ind[index+1]]
    
            
#     def get_cells(self):
        
#         if self.ordered:
#             for i in range(self.ncells):
#                 yield i     
#         else:
#             for i in cell_to_index:
#                 yield i
        
#     def construct_cofacets(self):
        
#         self.cofacet_ind = np.empty(self.ncells+1)
#         self.cofacet_ind[0] = 0
#         self.cofacets = np.empty(self.nfacets)
        
#         if self.ordered:
#             cell_list = [[] for i in range(self.ncells)]
#             for i in range(self.ncells):
#                 for j in self.facets[self.facet_ind[i]:self.facet_ind[i+1]]:
#                     cell_list[j].append(i)
                    
#             for i in range(self.ncells):
#                 self.cofacet_ind[i+1] = self.cofacet_ind[i] + len(cell_list[i])
#                 self.cofacets[self.cofacet_ind[i]:self.cofacet_ind[i+1]] = cell_list[i]
        
        
        
        
        
#         if self.ordered:
#             self.cofacets = [set() for i in range(len(self.facets))]

#             for i in range(len(self.facets)):
#                 for j in self.facets[i]:
#                     self.cofacets[j].add(i)
                    
#         else:
#             self.cofacets = {i:set() for i in self.facets}
            
#             for i in self.facets:
#                 for j in self.facets[i]:
#                     self.cofacets[j].add(i) 

class CellComplex:
    
    def __init__(self, dim, regular=True, oriented=False, ordered=True):
        self.dim = dim
        self.ordered = ordered
        self.regular = regular
        self.oriented = oriented
        if self.ordered:
            self.dims = []
            self.facets = []
        else:
            self.dims = {}
            self.facets = {}
        
    def add_cell(self, dim, facets, label=None, coeffs=None):
        if self.oriented or not self.regular:
            f = {facet:coeff for (facet, coeff) in zip(facets, coeffs)}
        else:
            f = set(facets)            
             
        if self.ordered:
            if label is None:
                self.dims.append(dim)
                self.facets.append(f)
            else:
                self.dims.insert(label, dim)
                self.facets.insert(label, f)
        else:
            if label is None:
                self.dims[len(self.dims)] = dim
                self.facets[len(self.facets)] = f
            else:
                self.dims[label] = dim
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
            self.cofacets = {i:set() for i in self.facets}
            
            for i in self.facets:
                for j in self.facets[i]:
                    self.cofacets[j].add(i) 
        
    
            
def construct_cubical_complex(shape, oriented=False, dual=False):
        
    dim = len(shape)
    
    comp = CellComplex(dim, regular=True, oriented=oriented, ordered=True)
    
    if dim == 2 and not dual:
        nrows = shape[0]
        ncols = shape[1]
          
        # vertices
        for i in range(nrows*ncols):
            
            if oriented:
                coeffs = []
            else:
                coeffs = None
            comp.add_cell(0, [], coeffs=coeffs)
        
        # horizontal edges
        for i in range(nrows):
            for j in range(ncols-1):
                if oriented:
                    coeffs = [-1, 1]
                else:
                    coeffs = None
                comp.add_cell(1, [ncols*i + j, ncols*i + j+1], coeffs=coeffs)

        # vertical edges
        for i in range(nrows-1):
            for j in range(ncols):
                if oriented:
                    coeffs = [-1, 1]
                else:
                    coeffs = None
                comp.add_cell(1, [ncols*i + j, ncols*(i+1) + j], coeffs=coeffs)
        
        # faces
        for i in range(nrows-1):
            for j in range(ncols-1):
                if oriented:
                    coeffs = [1, 1, -1, -1]
                else:
                    coeffs = None
                comp.add_cell(2, [nrows*ncols + (ncols-1)*i + j, 
                             nrows*ncols + (ncols-1)*nrows + ncols*i + j+1, 
                             nrows*ncols + (ncols-1)*(i+1) + j, 
                             nrows*ncols + (ncols-1)*nrows + ncols*i + j], coeffs=coeffs)
                
    elif dim == 2 and dual:
        nrows = shape[0]+1
        ncols = shape[1]+1
        
        nverts = nrows*ncols
        nhedges = nrows*(ncols-1)
        nvedges = (nrows-1)*ncols
        nfaces = (nrows-1)*(ncols-1)
                
        # faces
        for i in range(nrows-1):
            for j in range(ncols-1):
                if oriented:
                    coeffs = [1, 1, -1, -1]
                else:
                    coeffs = None
                comp.add_cell(2, [nfaces + (ncols-1)*i + j, 
                             nfaces + nhedges + ncols*i + j+1, 
                             nfaces + (ncols-1)*(i+1) + j, 
                             nfaces + nhedges + ncols*i + j], coeffs=coeffs)
                
                
        # horizontal edges
        for i in range(nrows):
            for j in range(ncols-1):
                if oriented:
                    coeffs = [-1, 1]
                else:
                    coeffs = None
                comp.add_cell(1, [nfaces + nhedges + nvedges + ncols*i + j, 
                                  nfaces + nhedges + nvedges + ncols*i + j+1], coeffs=coeffs)

        # vertical edges
        for i in range(nrows-1):
            for j in range(ncols):
                if oriented:
                    coeffs = [-1, 1]
                else:
                    coeffs = None
                comp.add_cell(1, [nfaces + nhedges + nvedges + ncols*i + j, 
                                  nfaces + nhedges + nvedges + ncols*(i+1) + j], coeffs=coeffs)
            
            
        # vertices
        for i in range(nrows*ncols):
            
            if oriented:
                coeffs = []
            else:
                coeffs = None
            comp.add_cell(0, [], coeffs=coeffs)
        
        
        
        
        
    
        
    return comp
   
def check_boundary_op(comp):
    
    valid = True
    for i in comp.get_cells():
        if comp.dims[i] > 0:
            
            if not comp.regular or comp.oriented:
                sub_faces = {}
                for j in comp.facets[i]:
                    for k in comp.facets[j]:
                        sub_faces[k] = sub_faces.get(k, 0) + comp.facets[i][j] * comp.facets[j][k]
                        
                for j in sub_faces:
                    if (comp.oriented and sub_faces[j] != 0) or (not comp.oriented and sub_faces[j] % 2 != 0):
                        print("Error:", i, j, sub_faces[j])
                        valid = False
            else:
                sub_faces = set()
                for j in comp.facets[i]:
                    sub_faces ^= comp.facets[j]
                   
                if len(sub_faces) != 0:
                    print("Error:", i, sub_faces)
                    valid = False
                    
    return valid
            
    
    
def construct_vertex_filtration_order(comp, vertex_time, euclidean=False, positions=None):
    
    
    vertex_order = np.zeros(len(vertex_time), int)
        
    submerged = np.zeros(len(vertex_time), bool)
        
    vertex_argsort = list(np.argsort(vertex_time)[::-1])
        
    
    visited = set()
    count = 0
    
    ti = 0
    while len(vertex_argsort) > 0:
        unvisited = set()
            
        # find all vertices in level set
        while True:
 
            vi = vertex_argsort.pop()
                    
            t = vertex_time[vi]
            unvisited.add(vi)
                        
            if len(vertex_argsort) == 0 or t != vertex_time[vertex_argsort[-1]]:
                break
          

        count += len(unvisited)
        
        # if only a single vertex, then take a shortcut
        if len(unvisited) == 1:
            vi = unvisited.pop()
            submerged[vi] = True
            vertex_order[vi] = ti
            ti += 1
            
            visited.add(vi)
            
            continue
            
        
        # set up queue
        dist = {}
        if euclidean:
            Q = queue.PriorityQueue()
            closest = {}
        else:
            Q = queue.Queue() 
            
        # put vertices in queue if they have a previously submerged neighbor
        for a in unvisited:
            for b in comp.cofacets[a]:
                for c in comp.facets[b]:
                    if submerged[c]:
                        
                        if euclidean:
                            new_dist = la.norm(positions[a] - positions[c])                    
                        else:
                            new_dist = 1.0
                            
                        if a not in dist or new_dist < dist[a]:
                            dist[a] = new_dist
                            Q.put([new_dist, a])
                            
                            if euclidean:
                                closest[a] = c
                            

        # perform BFS on level set
        while not Q.empty():
                        
            [current_dist, a] = Q.get()    
                        
            visited.add(a)
                
            if a not in unvisited:
                continue
                    
            unvisited.remove(a)
            submerged[a] = True
            vertex_order[a] = ti
            ti += 1
              
            for b in comp.cofacets[a]:
                for c in comp.facets[b]:
                    if c in unvisited:
                    
                        if euclidean:
                            new_dist = la.norm(positions[c] - positions[closest[a]])
                        else:
                            new_dist = current_dist + 1.0  
                          
                        if c not in dist or new_dist < dist[c]:
                            dist[c] = new_dist
                            Q.put([new_dist, c])
                            
                            if euclidean:
                                closest[c] = closest[a]
      
    
        # search through new segments
        while len(unvisited) > 0:

            a = unvisited.pop()
            
            Q.put([0, a])
            
            if euclidean:
                closest[a] = a
            
            while not Q.empty():
                [current_dist, a] = Q.get()
                
                visited.add(a)
                
                unvisited.discard(a)
                submerged[a] = True
                vertex_order[a] = ti
                ti += 1
                
                for b in comp.cofacets[a]:
                    for c in comp.facets[b]:
                        if c in unvisited:

                            if euclidean:
                                new_dist = la.norm(positions[c] - positions[closest[a]])
                            else:
                                new_dist = current_dist + 1.0  

                            if c not in dist or new_dist < dist[c]:
                                dist[c] = new_dist
                                Q.put([new_dist, c])

                                if euclidean:
                                    closest[c] = closest[a]    
    return vertex_order  
    
    
def get_star(alpha, I, cell_dim=None, dims=None, insert_order=None):
      
    if insert_order is not None:
        STAR = 1
        star_index = insert_order[alpha][STAR]

    star = set()

    Q = co.deque()
    Q.append(alpha)

    while len(Q) > 0:
        a = Q.popleft()

        if insert_order is None or star_index == insert_order[a][STAR]:
            Q.extend(I[a])
            if cell_dim is None or dims[a] == cell_dim:
                star.add(a)
                
    return star
    

def get_lex_val(alpha, I, cell_dim, dims, F):
    
    lex_val = set()
        
    cells = get_star(alpha, I, cell_dim=cell_dim, dims=dims)
    for c in cells:
        lex_val.add((F[c], c))
        
    return sorted(lex_val, reverse=True)
      
    
# construct map of cells to insertion times
# insertion times are a triple of 
# [insertion time, 
# lower star (upper costar) cell, 
# lexicographic insertion index (index ordering on all cells)]
# where the insertion index is found by computing the lexicographic order of the cells
def construct_time_of_insertion_map(comp, vertex_time, vertex_order, dual=False):
    
    order_to_time = vertex_time[np.argsort(vertex_order)]
                
    lex = []
    
    if dual:
        # Order by:
        # 1. function value of cell whose upper costar this cell belongs to
        # 2. cell dimension
        # 3. lexicographic ordering of function values of highest dimension cofaces from high to low
        for c in comp.get_cells():
            lex_val = get_lex_val(c, comp.cofacets, comp.dim, comp.dims, vertex_order)
            lex.append((lex_val[-1][0], comp.dims[c], lex_val, c, lex_val[-1][1]))
    else:
        # Order by:
        # 1. function value of vertex whose lower star this cell belongs to
        # 2. cell dimension
        # 3. lexicographic ordering of function values of vertices from high to low
        for c in comp.get_cells():
            lex_val = get_lex_val(c, comp.facets, 0, comp.dims, vertex_order)
            lex.append((lex_val[0][0], comp.dims[c], lex_val, c, lex_val[0][1]))
      
    lex = sorted(lex)
    
    if comp.ordered:
        insert_order = [None for c in comp.get_cells()]
    else:
        insert_order = {}
    
    for i, (star_val, d, lex_val, c, star) in enumerate(lex):
        insert_order[c] = (order_to_time[star_val], star, i)  

    return insert_order
 
    
    
def construct_discrete_gradient(comp, insert_order, dual=False):
    
    STAR = 1
    LEX = 2     
        
    V = {}
        
    # iterate through each vertex (or dual cell)
    for x in comp.get_cells():
        if (dual and comp.dims[x] != comp.dim) or (not dual and comp.dims[x] != 0):
            continue
    
        if dual:
            # get upper costar
            star = get_star(x, comp.facets, insert_order=insert_order)
            # star = get_ulstar(x, comp.facets)
        else:
            # get lower star
            star = get_star(x, comp.cofacets, insert_order=insert_order)
            # star = get_ulstar(x, comp.cofacets)
            
        unpaired = {}
            
        if dual:
            for alpha in star:
                unpaired[alpha] = star & set(comp.cofacets[alpha])
        else:
            for alpha in star:
                unpaired[alpha] = star & set(comp.facets[alpha])

        # print("cell", x)
        # print(star)
        # print(unpaired)
        
        PQone = queue.PriorityQueue()
        PQzero = queue.PriorityQueue()

        for alpha in star:
            if len(unpaired[alpha]) == 0:
                if dual:
                    PQzero.put((-insert_order[alpha][LEX], alpha))
                else:
                    PQzero.put((insert_order[alpha][LEX], alpha))
                    # print("zero", alpha)
            elif len(unpaired[alpha]) == 1:
                if dual:
                    PQone.put((-insert_order[alpha][LEX], alpha))
                else:
                    PQone.put((insert_order[alpha][LEX], alpha))
                    # print("one", alpha)
                
        while not PQone.empty() or not PQzero.empty():
            
            while not PQone.empty():
                (order, alpha) = PQone.get()
                
                if len(unpaired[alpha]) == 0:
                    PQzero.put((order, alpha))
                    continue
                
                beta = unpaired[alpha].pop()
                if dual:
                    V[alpha] = beta
                else:
                    V[beta] = alpha
                
                # print(beta, alpha)
                
                star.discard(alpha)
                star.discard(beta)
                
                del unpaired[alpha]
                for gamma in star:
                    if alpha in unpaired[gamma] or beta in unpaired[gamma]:
                        unpaired[gamma].discard(alpha)
                        unpaired[gamma].discard(beta)
                        
                        if len(unpaired[gamma]) == 1:
                            if dual:
                                PQone.put((-insert_order[gamma][LEX], gamma))
                            else:
                                PQone.put((insert_order[gamma][LEX], gamma))
                        
                    
            if not PQzero.empty():
                (order, alpha) = PQzero.get()
                if alpha in star:
                    V[alpha] = alpha
                    
                    # print(alpha, alpha)
                    
                    star.discard(alpha)
                    
                    del unpaired[alpha]
                    for gamma in star:
                        if alpha in unpaired[gamma]:
                            unpaired[gamma].discard(alpha)

                            if len(unpaired[gamma]) == 1:
                                if dual:
                                    PQone.put((-insert_order[gamma][LEX], gamma))
                                else:
                                    PQone.put((insert_order[gamma][LEX], gamma))
                                    
            
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
            nrIn[c] = nrIn.get(c, 0) + 1
            
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
def calc_morse_boundary(s, V, I, oriented=False):
    
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
                mult[c] = mult.get(c, 0) + mult[a] * ( -I[a][b] )
            
        elif oriented:
            mult[c] = mult.get(c, 0) + mult[a] * ( -I[a][b] * I[c][b] )
            
    for c in boundary:
        if oriented:
            yield (c, counts[c], mult[c])  
        else:
            yield (c, counts[c], counts[c] % 2) 
    
def construct_morse_complex(V, I, comp, oriented=False):
    
    mcomp = CellComplex(comp.dim, ordered=False, oriented=oriented, regular=False)
    
    for s in V:
        if V[s] == s:
            facets = []
            coeffs = []
            for (c, k, m) in calc_morse_boundary(s, V, I, oriented=oriented):
                facets.append(c)
                coeffs.append(m)
                
                if abs(m) > 1:
                    print("Large Coefficient", k, m, comp.dims[c], comp.dims[s])
            
         
            # facets = [c for (c, k) in calc_morse_boundary(s, V, I) if k % 2 == 1]
                
            mcomp.add_cell(comp.dims[s], facets, label=s, coeffs=coeffs)
            
    return mcomp
    
        
# find path from s to t
def find_connections(s, t, V, coV, I, coI):
    
    active = set([t])
    for (a, b, c) in traverse_flow(t, coV, coI):
        active.add(c)
        
    for (a, b, c) in traverse_flow(s, V, I):
        if b in active:
            yield (a, b, c)
    


def simplify_morse_complex(threshold, V, coV, comp, insert_order):
     
    TIME = 0
    LEX = 2
        
    crit_cells = {s for s in V if V[s] == s}
        
    n = 0
    while True:
        
        n += 1
                 
        print("Pass:", n)
            
        close_pairs = []
        
        # print(sorted(crit_cells))
        print("Critical Cells", len(crit_cells))
                
        for s in crit_cells:
            
            if comp.dims[s] == 0:
                continue
                
            # print("s", s, insert_order[s], comp.dims[s])
            
            close_alpha = None
            close_alpha_time = None
            for (c, k, m) in calc_morse_boundary(s, V, comp.facets, oriented=comp.oriented):

                if k == 1:
                    if close_alpha is None or ((insert_order[c][TIME], insert_order[c][LEX]) > close_alpha_time):
                        close_alpha = c
                        close_alpha_time = (insert_order[c][TIME], insert_order[c][LEX])
                    
            if close_alpha is None:
                continue
                
                                
            close_beta = None
            close_beta_time = None
            for (c, k, m) in calc_morse_boundary(close_alpha, coV, comp.cofacets, oriented=comp.oriented):
                if k == 1:
                    if close_beta is None or ((insert_order[c][TIME], insert_order[c][LEX]) < close_beta_time):
                        close_beta = c
                        close_beta_time = (insert_order[c][TIME], insert_order[c][LEX])
               
           
            if s == close_beta:
                close_pairs.append((insert_order[s][TIME] - insert_order[close_alpha][TIME], (close_alpha, s)))
        
        
        close_pairs = sorted(close_pairs)
        
        print("Close Pairs:", len(close_pairs))
        
        if len(close_pairs) == 0 or close_pairs[0][0] > threshold:
            print("Cancelled Pairs:", 0)
            break
          
        i = 0
        for (time, (t, s)) in close_pairs:
            if time > threshold:
                break
        
            # print(i, time, t, s, comp.dims[t], comp.dims[s])
            i += 1
            reverse_pairs = []
            for (a, b, c) in find_connections(s, t, V, coV, comp.facets, comp.cofacets):
                reverse_pairs.append((a, b))
                
            for (a, b) in reverse_pairs:                
                V[b] = a
                
                coV[a] = b
                
                if a in V:
                    del V[a]
                if b in coV:
                    del coV[b]
                    
            crit_cells.remove(t)
            crit_cells.remove(s)
               
        print("Cancelled Pairs:", i)
        
        print("Remaining Critical Cells", len(crit_cells))
            
        
                
    return (V, coV)
       
        
def construct_filtration(comp, insert_order):
    
    LEX = 2
    
    Q = queue.PriorityQueue()
    for i in comp.get_cells():
        Q.put((insert_order[i][LEX], i))
        
    while not Q.empty():
        (order, c) = Q.get()
        yield c
    
    
           
        

    
        
def get_morse_weights(mcomp, V, coV, I, coI):
    
    weights = {}
    for s in mcomp.get_cells():
        if mcomp.dims[s] == 0:
            weights[s] = 1
            continue
            
        cell = set()
        for (t, m) in mcomp.facets[s].items():
            for (a, b, c) in find_connections(s, t, V, coV, I, coI):
                cell.add(a)
         
        weights[s] = len(cell)
    
    
    return weights
    
    
                
def compute_persistence(comp, filtration, extended=False, birth_cycles=False, optimal_cycles=False,
                       weights=None, relative_cycles=False):
    
    columns = []
    
    cell_to_col = {}
    col_to_cell = []
       
    icol = 0
        
        
    for ci in filtration:
        
        # Z2 coefficients
        if not optimal_cycles:
            col = set()
            for cj in comp.facets[ci]:
                if comp.regular:
                    col.add(cell_to_col[cj])
                else:
                    if comp.facets[ci][cj] != 0:
                        col.add(cell_to_col[cj])
                
            if birth_cycles:
                columns.append(col)
            else:
                columns.append((comp.dims[ci], sorted(col)))
                
        # general coefficients
        else:
            col = {}
            for cj in comp.facets[ci]:
                if comp.facets[ci][cj] != 0:
                    col[cell_to_col[cj]] = comp.facets[ci][cj]

            columns.append(col)
        
        cell_to_col[ci] = icol
        icol += 1
        col_to_cell.append(ci)

    
    if not optimal_cycles and not birth_cycles:
    
        boundary_matrix = phat.boundary_matrix(columns=columns)

        alive = set(range(len(columns)))

        pairs = []
        for (i, j) in boundary_matrix.compute_persistence_pairs():
            alive.discard(i)
            alive.discard(j)

            ci = col_to_cell[i]
            cj = col_to_cell[j]

            pairs.append((ci, cj))


        for ci in alive:
            pairs.append((col_to_cell[ci], None))
        
        return pairs
    
    elif not optimal_cycles and birth_cycles:
    
        # pivot row of each column if it has one
        pivot_row = {}
        # row to reduced column with pivot in that row
        pivot_col = {}
        
        g = []
            
        for j in range(len(columns)):
            
            g.append({j})
            
            if len(columns[j]) > 0:
                pivot_row[j] = max(columns[j])
            
            while len(columns[j]) > 0 and pivot_row[j] in pivot_col:
                
                l = pivot_col[pivot_row[j]]
                columns[j] ^= columns[l]
                
                g[j] ^= g[l]
                
                if len(columns[j]) > 0:
                    pivot_row[j] = max(columns[j])
                else:
                    del pivot_row[j]
                
                
                    
            if len(columns[j]) > 0:
                pivot_col[pivot_row[j]] = j
 
        
        alive = set(range(len(columns)))

        pairs = []
        for (j, i) in pivot_row.items():
            alive.discard(i)
            alive.discard(j)

            ci = col_to_cell[i]
            cj = col_to_cell[j]
            
            pairs.append((ci, cj))


        for ci in alive:
            pairs.append((col_to_cell[ci], None))

        bcycles = {}
        for i in range(len(columns)):
            if len(columns[i]) > 0 or comp.dims[col_to_cell[i]] == 0:
                continue

            ci = col_to_cell[i]

            bcycles[ci] = set()
            for j in g[i]:
                bcycles[ci].add(col_to_cell[j])
                                            
        return (pairs, bcycles)
    
        
    elif optimal_cycles:
                
        # pivot row of each column if it has one
        pivot_row = {}
        # row to reduced column with pivot in that row
        pivot_col = {}
        
        g = []
        
        ocycles = {}
        
        cell_counts = {i:0 for i in range(1, comp.dim+1)}
        
        x_to_cell = {i:{} for i in range(1, comp.dim+1)}
        cell_to_x = {i:{} for i in range(1, comp.dim+1)}
        
        B = {i+1:{} for i in range(1, comp.dim+1)}
        
        Z = {i:{} for i in range(1, comp.dim+1)}
            
        for j in range(len(columns)):
            
            g.append({j:1})
            
            if len(columns[j]) > 0:
                pivot_row[j] = max(columns[j])
               
            d = comp.dims[col_to_cell[j]]
            if d > 0:
                x_to_cell[d][cell_counts[d]] = col_to_cell[j]
                cell_to_x[d][col_to_cell[j]] = cell_counts[d]
                cell_counts[d] += 1
                
            
            while len(columns[j]) > 0 and pivot_row[j] in pivot_col:
                
                p = pivot_row[j]
                l = pivot_col[p]
                
                r = 1.0 * columns[j][p] / columns[l][p]
                
                for k in columns[l]:
                    columns[j][k] = columns[j].get(k, 0) - r * columns[l][k]
                    if columns[j][k] == 0.0:
                        del columns[j][k]
                    
                for k in g[l]:
                    g[j][k] = g[j].get(k, 0) - r * g[l][k]
                    if g[j][k] == 0.0:
                        del g[j][k]
                                    
                if len(columns[j]) > 0:
                    pivot_row[j] = max(columns[j])
                else:
                    del pivot_row[j]
                
            if len(columns[j]) == 0 and d > 0:
                
                c = np.ones(2*(cell_counts[d]+len(B[d+1])+len(Z[d])), float)
                if weights is not None:
                    for k in range(cell_counts[d]):
                        c[k] = weights[x_to_cell[d][k]]
                        c[k+cell_counts[d]] = weights[x_to_cell[d][k]]
                else:
                    c[0:2*cell_counts[d]] = 1.0
                
                A_i = []
                A_j = []
                A_val = []
                
                for k in range(cell_counts[d]):
                    A_i.append(k)
                    A_j.append(k)
                    A_val.append(1.0)
                    
                    A_i.append(k)
                    A_j.append(k+cell_counts[d])
                    A_val.append(-1.0)
                    
                                    
                for bi, kj in enumerate(B[d+1]):
                    for ki in B[d+1][kj]:
                        A_i.append(cell_to_x[d][ki])
                        A_j.append(2*cell_counts[d] + bi)
                        A_val.append(-B[d+1][kj][ki])
                        
                        A_i.append(cell_to_x[d][ki])
                        A_j.append(2*cell_counts[d] + len(B[d+1]) + bi)
                        A_val.append(B[d+1][kj][ki])
                 
                for zi, kj in enumerate(Z[d]):
                    for ki in Z[d][kj]:
                        A_i.append(cell_to_x[d][ki])
                        A_j.append(2*cell_counts[d] + 2*len(B[d+1]) + zi)
                        A_val.append(Z[d][kj][ki])

                        A_i.append(cell_to_x[d][ki])
                        A_j.append(2*cell_counts[d] + 2*len(B[d+1]) + len(Z[d]) + zi)
                        A_val.append(Z[d][kj][ki])

                
                b_eq = np.zeros(cell_counts[d], float)
                # print(g[j])
                for k in g[j]:
                    b_eq[cell_to_x[d][col_to_cell[k]]] = g[j][k]
                
                res = opt.linprog(c, A_eq=sparse.coo_matrix((A_val, (A_i, A_j))), 
                                  b_eq=b_eq, method='interior-point', 
                                  options={'disp':False, 'maxiter': 100000, 'sparse': True, 
                                           'ip': False, 'permc_spec':'COLAMD'})
                
                if res.status != 0:
                    print(res)
                                
                
                print(j, "/", len(columns), "size", len(c), "nit", res.nit, "Change", np.sum([weights[x_to_cell[d][k]] for k in np.where(b_eq != 0)[0]]), "->", int(round(res.fun)))
                
                z = res.x[0:cell_counts[d]] - res.x[cell_counts[d]:2*cell_counts[d]]
                
                # print(z)
                
                col = {}
                for k in range(cell_counts[d]):
                    h = int(round(z[k]))
                    if h != 0:
                        col[x_to_cell[d][k]] = h
                        
                if relative_cycles:
                    Z[d][j] = col
                    
                ocycles[col_to_cell[j]] = set(col.keys())
                    
            elif len(columns[j]) > 0:
                p = pivot_row[j]
                pivot_col[p] = j
                
                if d > 1:
                    col = {}
                    for k in columns[j]:
                        if columns[j][k] != 0.0:
                            col[col_to_cell[k]] = columns[j][k]
                    B[d][j] = col
                
                    if relative_cycles:
                        del Z[d-1][p]
                            
           
        alive = set(range(len(columns)))

        pairs = []
        for (j, i) in pivot_row.items():
            alive.discard(i)
            alive.discard(j)

            ci = col_to_cell[i]
            cj = col_to_cell[j]

            pairs.append((ci, cj))


        for ci in alive:
            pairs.append((col_to_cell[ci], None))

        if birth_cycles:
            
            bcycles = {}
            for i in range(len(columns)):
                if len(columns[i]) > 0 or comp.dims[col_to_cell[i]] == 0:
                    continue

                ci = col_to_cell[i]

                bcycles[ci] = set()
                for j in g[i]:
                    bcycles[ci].add(col_to_cell[j])
                
            return (pairs, bcycles, ocycles)
        else:
            return (pairs, ocycles)
        

def convert_morse_to_real_complex(mfeature, V, I):

    feature = set()
    for s in mfeature:
        feature.add(s)
        for (a, b, c) in traverse_flow(s, V, I):
            if b != c:
                feature.add(c)
    
    return feature


def convert_to_pixels(feature, comp, insert_order, dual=False):
    
    pixels = set()
    
    STAR = 1
    LEX = 2
    
    for s in feature:
        
        if comp.dims[s] == 0:
            if not dual:
                
                pixels.add(s)
                
            else:
                
                #start replacing snippets here
                
                cofaces = get_star(s, comp.cofacets, cell_dim=comp.dim, dims=comp.dims)
                        
                for c in cofaces:
                    
                    verts = get_star(c, comp.facets, cell_dim=0, dims=comp.dims)
                    
                    ucostar_verts = set()
                    for v in verts:
                        if insert_order[c][STAR] == insert_order[v][STAR]:
                            ucostar_verts.add(v)
                
                    # no verts in ucostar:
                    if len(ucostar_verts) == 0:
                        min_vert = verts.pop()
                        
                        for v in verts:
                            if insert_order[v][LEX] < insert_order[min_vert][LEX]:
                                min_vert = v
                                
                        if min_vert in feature:
                            pixels.add(c)
                     
                    # all verts in ucostar are in feature
                    elif ucostar_verts.issubset(feature):
                        
                        pixels.add(c)
                      
                    # only some ucostar verts are in feature
                    else:
                        
                        min_vert = ucostar_verts.pop()
                        for v in ucostar_verts:
                            if insert_order[v][LEX] < insert_order[min_vert][LEX]:
                                min_vert = v
                                
                        if min_vert in feature:
                            pixels.add(c)
                            
        elif comp.dims[s] == 1:
            
            if not dual:
                
                pixels.update(comp.facets[s])
                
            else:
                        
                ucostar = insert_order[s][STAR]
                    
                pixels.add(ucostar)
                
        elif comp.dims[s] == 2:
            
            if not dual:
                for c in comp.facets[s]:
                    pixels.update(comp.facets[c])
                    
            else:
                
                pixels.add(s)
                            
                
                        
    return pixels                
                            
     
def find_basins(mcomp, coV, comp, insert_order, dual=False):
    
    basins = {}
    
    STAR = 1
    
    for c in mcomp.get_cells():
        if mcomp.dims[c] == 0:
            
            if dual:
                index = insert_order[c][STAR]
            else:
                index = c
                
            basins[index] = convert_to_pixels(convert_morse_to_real_complex({c}, coV, comp.cofacets), 
                                                comp, insert_order, dual=dual)
     
   
    return basins
        
    
def find_morse_skeleton(mcomp, V, comp, d, insert_order, dual=False):
    skeleton = set()
    
    for s in V:
        if V[s] == s and mcomp.dims[s] == d:
            feature = convert_to_pixels(convert_morse_to_real_complex({s}, V, comp.facets), comp, insert_order, dual=dual)
            skeleton.update(feature)
            
            
    return skeleton 

        

# if (i, j) represents a 1-cycle (connected component), 
# then expand cell i in morse complex up to (but not including) death cell j
# if (i, j) represents a d-cycle (cycle for 2-manifold or void for 3-manifold), 
# then expand cell j in morse complex up to (but not including) birth cell i    


def extract_persistence_feature(i, j, mcomp, comp, V, coV, insert_order):
    
    TIME = 0
    STAR = 1
    LEX = 2
        
    if comp.dims[i] == 0:
        I = mcomp.facets
        coI = mcomp.cofacets
        
        barrier_test = lambda b: insert_order[b][TIME] >= insert_order[j][TIME]
        
        seen = {i}
        Q = co.deque([i])
    
    else:
        I = mcomp.cofacets
        coI = coI = mcomp.facets
        
        barrier_test = lambda b: insert_order[b][TIME] <= insert_order[i][TIME]
        
        seen = {j}
        Q = co.deque([j])
    
    while len(Q) > 0:
        a = Q.popleft()
        for b in coI[a]:
            if barrier_test(b):
                continue
                
            for c in I[b]:
                if c != a and c not in seen:
                    Q.append(c)
                    seen.add(c)
            
    if comp.dims[i] == 0:
        feature = convert_morse_to_real_complex(seen, coV, comp.cofacets)
    else:
        feature = convert_morse_to_real_complex(seen, V, comp.facets)
    
    return feature
        


def get_boundary(cells, comp):
    
    cycle = set()
    
    if comp.regular:
        for c in cells:
            cycle ^= comp.facets[c]
        
    return cycle

        
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

        