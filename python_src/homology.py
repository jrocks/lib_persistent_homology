import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.optimize as opt
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
    
    def __init__(self, dim, regular=True, oriented=False, ordered=True):
        self.dim = dim
        self.ordered = ordered
        self.regular = regular
        self.oriented=oriented
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
            

def construct_cubical_complex(shape, oriented=False):
        
    dim = len(shape)
    
    comp = Complex(dim, regular=True, oriented=oriented, ordered=True)
    
    if dim == 2:
        nrows = shape[0]
        ncols = shape[1]
                
        for i in range(nrows*ncols):
            
            if oriented:
                coeffs = []
            else:
                coeffs = None
            comp.add_cell(0, [], coeffs=coeffs)
        
        for i in range(nrows):
            for j in range(ncols-1):
                if oriented:
                    coeffs = [-1, 1]
                else:
                    coeffs = None
                comp.add_cell(1, [ncols*i + j, ncols*i + j+1], coeffs=coeffs)

        for i in range(nrows-1):
            for j in range(ncols):
                if oriented:
                    coeffs = [-1, 1]
                else:
                    coeffs = None
                comp.add_cell(1, [ncols*i + j, ncols*(i+1) + j], coeffs=coeffs)
        
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
# counts = 1: there is a single path and (s, c) are a cancellable close pair
# counts = 2: there are an even number of paths and c is not a boundary
# counts = 3: there are an odd number of paths and (s, c) are not a cancellable close pair
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
    
    mcomp = Complex(comp.dim, ordered=False, oriented=oriented, regular=False)
    
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

            

            
def construct_filtration(I, dims, comp, F):
    
    Q = queue.PriorityQueue()
    for i in comp.get_cells():
        Q.put((get_lex_val(i, I, dims, F), i))
        
    while not Q.empty():
        (lex, c) = Q.get()
        yield (max(lex), c)

                
def compute_persistence(comp, filtration, show_zero=False, extended=False, birth_cycles=False, optimal_cycles=False):
    
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
        
        # print(pi, h, ci, d)
        
        
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

        cell_index[ci] = (pi, h)             

    
    if not optimal_cycles and not birth_cycles:
    
        boundary_matrix = phat.boundary_matrix(columns=columns)

        alive = set(range(len(columns)))

        pairs = []
        for (i, j) in boundary_matrix.compute_persistence_pairs():
            alive.discard(i)
            alive.discard(j)

            ci = col_to_cell[i]
            cj = col_to_cell[j]
            if not show_zero and cell_index[ci] == cell_index[cj]:
                continue

            pairs.append((ci, cj))


        for ci in alive:
            pairs.append((col_to_cell[ci], None))
        
        return (pairs, cell_index)
    
    elif not optimal_cycles and birth_cycles:
    
        # pivot row of each column if it has one
        pivot_row = {}
        # row to reduced column with pivot in that row
        pivot_col = {}
        
        g = []
            
        for j in range(len(columns)):
            
            g.append(set([j]))
            
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
            if not show_zero and cell_index[ci] == cell_index[cj]:
                continue

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
                                            
        return (pairs, cell_index, bcycles)
    
        
    elif optimal_cycles:
                
        # pivot row of each column if it has one
        pivot_row = {}
        # row to reduced column with pivot in that row
        pivot_col = {}
        
        g = []
        
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
                    
                    
                # change B and Z to list of column dictionaries
                
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
                
#                 for (Bi, Bj, Bval) in B[d+1]:
#                     A_i.append(Bi)
#                     A_i.append(2*cell_counts[d] + Bj)
#                     A_val.append(-Bval)
                    
#                     A_i.append(Bi)
#                     A_i.append(2*cell_counts[d] + bound_counts[d+1] + Bj)
#                     A_val.append(Bval)
                    
#                 for (Zi, Zj, Zval) in Z[d+1]:
#                     A_i.append(Zi)
#                     A_i.append(2*cell_counts[d] + 2*bound_counts[d+1] + Zj)
#                     A_val.append(Zval)
                    
#                     A_i.append(Zi)
#                     A_i.append(2*cell_counts[d] + 2*bound_counts[d+1] + rel_counts[d] + Zj)
#                     A_val.append(-Zval)
                               
                
                b_eq = np.zeros(cell_counts[d], float)
                # print(g[j])
                for k in g[j]:
                    b_eq[cell_to_x[d][col_to_cell[k]]] = g[j][k]
                
                bounds = [(0, None) for k in range(2*(cell_counts[d] + len(B[d+1]) + len(Z[d])))]
                
                # print(cell_counts[d], len(B[d+1]), len(Z[d]))
                # print(sparse.coo_matrix((A_val, (A_i, A_j))).todense())
                # print(b_eq)
                
                res = opt.linprog(c, A_eq=sparse.coo_matrix((A_val, (A_i, A_j))).todense(), 
                                  b_eq=b_eq, bounds=bounds, method='simplex', options={'disp':True})
                
                print(res)
                
                print(np.sum(np.abs(b_eq)), "->", res.fun)
                
                z = res.x[0:cell_counts[d]] - res.x[cell_counts[d]:2*cell_counts[d]]
                
                # print(z)
                
                col = {}
                for k in range(cell_counts[d]):
                    if z[k] != 0.0:
                        col[x_to_cell[d][k]] = z[k]
                Z[d][j] = col
                    
                    
            elif len(columns[j]) > 0:
                p = pivot_row[j]
                pivot_col[p] = j
                
                if d > 1:
                    col = {}
                    for k in columns[j]:
                        if columns[j][k] != 0.0:
                            col[col_to_cell[k]] = columns[j][k]
                    B[d][j] = col
                
                    del Z[d-1][p]
                            
           
        alive = set(range(len(columns)))

        pairs = []
        for (j, i) in pivot_row.items():
            alive.discard(i)
            alive.discard(j)

            ci = col_to_cell[i]
            cj = col_to_cell[j]
            if not show_zero and cell_index[ci] == cell_index[cj]:
                continue

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
            
            ocycles = {}
                
            return (pairs, cell_index, bcycles, ocycles)
        else:
            return (pairs, cell_index, ocycles)
        
    
            
def find_basins(coV, coI, dims, d):
        
    basins = {}
    
    for s in coV:
        if coV[s] == s and dims[s] == d:
            basins[s] = set([s])
            for (a, b, c) in traverse_flow(s, coV, coI):
                if dims[c] == d:
                    basins[s].add(c)
                
                
    return basins
    
    
def level_bfs(s, h, I, coI, F):
    
    Q = co.deque([s])
    seen = {s}
    
    while len(Q) > 0:
        a = Q.popleft()
        for b in coI[a]:
            if F[b][1] >= h:
                continue
            for c in I[b]:
                if c != a and c not in seen:
                    Q.append(c)
                    seen.add(c)
            
            
            
    return seen
    
# might want to use calc_morse_boundary instead of the Z_2 morse complex
def find_segment(i, j, basins, cell_index, mcomp):
    
    seen = level_bfs(i, cell_index[j][1], mcomp.facets, mcomp.cofacets, cell_index)

    segment = set()
    for i in seen:
        segment.update(basins[i])
        
    return segment

def find_cycle(mcycle, mcomp, V, coV, I, coI):
    
    cycle = set()
    for s in mcycle:
        for (t, m) in mcomp.facets[s].items():
        # for (t, count) in calc_morse_boundary(s, V, I):
            for (a, b, c) in find_connections(s, t, V, coV, I, coI):
                cycle.add(b)
    
    return cycle

    
def find_morse_skeleton(mcomp, V, coV, I, coI, dims, d):
    
    skeleton = set()
    for s in V:
        if V[s] == s and dims[s] == d:
            for (t, m) in mcomp.facets[s].items():
                for (a, b, c) in find_connections(s, t, V, coV, I, coI):
                    skeleton.add(b)
                                 
    return skeleton
                
        
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

        