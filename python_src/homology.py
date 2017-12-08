import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import phat

# # simplices - main simplices in filtration. If the boundary of a simplex
# # is not specified, it is automatically added at the same level of the
# # filtration.
# # 
# # dimensions - the corresponding list of lists of simplex dimensions
# # 
# # heights - heights used for filtration
# # 
# # boundaries - boundaries for each possible simplex
# # 
# # returns lists of persistence pairs

# def compute_persistence_pairs(simplices, dimensions, heights, boundaries):
        
#     icol = 0
#     columns = []
    
#     sim_to_col = [{} for d in range(len(boundaries)+1)]
#     col_to_sim = []
    
#     sim_to_pindex = [{} for d in range(len(boundaries)+1)]
        
#     # go through groups of simplices in order from lowest to highest
#     hsort = np.unique(heights)
#     for pi, h in enumerate(hsort):
                
#         # this can be optimized
#         sindex = np.where(heights==h)[0]
        
#         # iterate through each simplex with same height
#         for si in sindex:
            
#             sim = simplices[si]
#             dim = dimensions[si]
            
#             if sim in sim_to_col[dim]:
#                 return "Redundancy Error"
                        
#             # tabulate all simplices to check at each dimension
#             all_sim = [set() for d in range(dim+1)]
#             all_sim[dim].add(sim)    
#             for d in range(dim, 0, -1):
#                 for sj in all_sim[d]:
#                     all_sim[d-1].update(boundaries[d][sj])
                               
#             # add vertices (simplices without boundary)
#             for vi in all_sim[0]:
#                 if vi not in sim_to_col[0]:
#                     columns.append((0, []))
                    
#                     sim_to_col[0][vi] = icol
#                     icol += 1
                    
#                     sim_to_pindex[0][vi] = pi
#                     col_to_sim.append(vi)
                    
#             # add simplices with boundary in order from lowest to highest dimension
#             # boundaries are guarenteed to always exist
#             for d in range(1, dim+1):
#                 for sj in all_sim[d]:
#                     if sj not in sim_to_col[d]:

#                         col = []
#                         for bi in boundaries[d][sj]:
#                             col.append(sim_to_col[d-1][bi])
                            
#                         columns.append((d, sorted(col)))
                        
#                         sim_to_col[d][sj] = icol
#                         icol += 1
                        
#                         sim_to_pindex[d][sj] = pi
#                         col_to_sim.append(sj)
                    

#     boundary_matrix = phat.boundary_matrix(columns=columns)
    
#     alive = [set() for i in range(len(boundaries)+1)]
#     for i, col in enumerate(columns):
#         alive[col[0]].add(i)
        
#     raw_pairs = boundary_matrix.compute_persistence_pairs()
    
#     ipairs = [{} for i in range(len(boundaries)+1)]
    
#     for (i, j) in raw_pairs:

#         d = columns[i][0]
#         pi = sim_to_pindex[d][col_to_sim[i]]
#         pj = sim_to_pindex[d+1][col_to_sim[j]]
        
#         if pi != pj:
            
#             if (pi, pj) not in ipairs[d]:
#                 ipairs[d][(pi, pj)] = 1
#             else:
#                 ipairs[d][(pi, pj)] += 1
            
#         alive[d].discard(i)
#         alive[d+1].discard(j)
        
#     persist = [{} for i in range(len(boundaries)+1)]
#     for d in range(len(boundaries)+1):
#         for i in alive[d]:
#             pi = sim_to_pindex[d][col_to_sim[i]]
#             if pi not in persist[d]:
#                 persist[d][pi] = 1
#             else:
#                 persist[d][pi] += 1
    
    
#     return (ipairs, hsort, persist, sim_to_pindex)


class SimplicialComplex:
    
    def __init__(self, dim, Nverts):
        self.dim = dim
        self.Nverts = Nverts
        self.faces = {i+1:[] for i in range(dim)}
        
    def add_simplex(self, dim, simplex):
        self.faces[dim].append(simplex)
        
        
def pixel_triangulation(Nx, Ny, compactify=False):
    
    Nverts = Nx*Ny
    if compactify:
        Nverts += 1
    comp = SimplicialComplex(2, Nverts)

    for i in range(Ny):
        for j in range(Nx-1):
            comp.add_simplex(1, [Nx*i + j, Nx*i + j+1])

    for i in range(Ny-1):
        for j in range(Nx):
            comp.add_simplex(1, [Nx*i + j, Nx*(i+1) + j])

#     for i in range(Ny-1):
#         for j in range(Nx-1):
#             comp.add_simplex(1, [Nx*i + j, Nx*(i+1) + j+1])
    
#     for i in range(Ny-1):
#         for j in range(Nx-1):
#             comp.add_simplex(2, [(Nx-1)*i + j, (Nx-1)*Ny + Nx*i + j+1, (Nx-1)*Ny + Nx*(Ny-1) + (Nx-1)*i + j])
#             comp.add_simplex(2, [(Nx-1)*(i+1) + j, (Nx-1)*Ny + Nx*i + j, (Nx-1)*Ny + Nx*(Ny-1) + (Nx-1)*i + j])
            
    for i in range(Ny-1):
        for j in range(Nx-1):
            comp.add_simplex(2, [(Nx-1)*i + j, (Nx-1)*Ny + Nx*i + j+1, (Nx-1)*(i+1) + j, (Nx-1)*Ny + Nx*i + j])
            
    if compactify:
        for j in range(Nx-1):
            comp.add_simplex(1, [j, Nx*Ny])
           
        for i in range(Ny-1):
            comp.add_simplex(1, [Nx*i + Nx-1, Nx*Ny])
            
        for j in range(Nx-1, 0, -1):
            comp.add_simplex(1, [Nx*(Ny-1) + j, Nx*Ny])
            
        for i in range(Ny-1, 0, -1):
            comp.add_simplex(1, [Nx*i, Nx*Ny])
        
        for j in range(Nx-1):
            comp.add_simplex(2, [j, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + j, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + j+1])
            
        for i in range(Ny-1):
            comp.add_simplex(2, [(Nx-1)*Ny + Nx*i + Nx-1, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + i, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + i+1])
        
        for j in range(Nx-1):
            comp.add_simplex(2, [(Nx-1)*Ny - 1 - j , 
                                 (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + Ny-1 + j, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + Ny-1 + j+1])
        
        for i in range(Ny-2):
            comp.add_simplex(2, [(Nx-1)*Ny + Nx*(Ny-1) - Nx - Nx*i, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + Ny-1 + Nx-1 + i, 
                                 (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + Ny-1 + Nx-1 + i+1])
        
        comp.add_simplex(2, [(Nx-1)*Ny + Nx*(Ny-1) - Nx - Nx*(Ny-2), 
                            (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 + Ny-1 + Nx-1 + Ny-2, 
                            (Nx-1)*Ny + Nx*(Ny-1)])
        
        # for j in range(Nx-1):    
        #     comp.add_simplex(1, [Nx*(Ny-1) + j, Nx*Ny])
            
#         for i in range(Ny):
#             comp.add_simplex(1, [Nx*i, Nx*Ny])
#             comp.add_simplex(1, [Nx*i + Nx-1, Nx*Ny])
            
#         for j in range(Nx-1):
#             comp.add_simplex(2, [j, (Nx-1)*Ny + Nx*(Ny-1) + j, (Nx-1)*Ny + Nx*(Ny-1) + j+1])
#             comp.add_simplex(2, [(Nx-1)*(Ny-1) + j, (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 j, (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 j+1])
            
        # for i in range(Ny-1):
        #     comp.add_simplex(2, [j, (Nx-1)*Ny + Nx*(Ny-1) + j, (Nx-1)*Ny + Nx*(Ny-1) + j+1])
        #     comp.add_simplex(2, [(Nx-1)*(Ny-1) + j, (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 j, (Nx-1)*Ny + Nx*(Ny-1) + Nx-1 j+1])
            
    return comp

def convert_to_vertices(comp, simps, dims, heights):
    
    marked = set()
    verts = []
    vheights = []
    
    
    for (si, dim, h) in zip(simps, dims, heights):
        all_sim = [set() for d in range(dim+1)]
        all_sim[dim].add(si)   
        for d in range(dim, 0, -1):
            for sj in all_sim[d]:
                all_sim[d-1].update(comp.faces[d][sj])
                
        for vi in all_sim[0]:
            if vi not in marked:
                marked.add(vi)
                verts.append(vi)
                vheights.append(h)
                
    return (verts, vheights)

# comp - simplicial complex
# verts - ordered list of vertices used in filtration. 
# Does not need to include all possible vertices.
def construct_lower_star_filtration(comp, verts):
    
    cofaces = {}
    cofaces[0] = [[] for i in range(comp.Nverts)]
    for d in range(1, comp.dim):
        cofaces[d] = [[] for i in range(len(comp.faces[d]))]
    
    for d in range(1, comp.dim+1):
        for si, simplex in enumerate(comp.faces[d]):
            for sj in simplex:
                cofaces[d-1][sj].append(si)
    
    simp_filt = []
    dims = []
    
    marked = {i: set() for i in range(comp.dim+1)}
    
    # go through vertices from lowest to highest
    for vi in verts:
        
        curr_simps = set([vi])
        new_simps = set()
        
        simp_filt.append(vi)
        dims.append(0)
        marked[0].add(vi)
                
        for d in range(comp.dim):
            for cfi in curr_simps:
                for cfj in cofaces[d][cfi]:
                    faces = set(comp.faces[d+1][cfj])
                    if cfj not in marked[d+1] and faces <= marked[d]:
                        simp_filt.append(cfj)
                        dims.append(d+1)
                        marked[d+1].add(cfj)
                        new_simps.add(cfj)
                                      
            curr_simps = new_simps
            new_simps = set()
     
    return (simp_filt, dims)


# simp_filt - ordered list of simplices used in filtration
# 
# dimensions - the corresponding list of lists of simplex dimensions
# 
# heights - ordered list of heights used for filtration (just defined for vertices)
# 
# comp - full simplicial complex
# 
# returns lists of persistence pairs
def compute_persistence_pairs(simp_filt, dims, heights, comp):
        
   
    columns = []
    
    sim_to_col = [{} for d in range(comp.dim+1)]
    col_to_sim = [] 
    
    sim_to_pindex = [{} for d in range(comp.dim+1)]
        
    # go through groups of simplices in order from lowest to highest
    hbins = np.unique(heights)
    isim = 0
    ih = 0
    for pi, h in enumerate(hbins):
                
        while isim < len(simp_filt):
            si = simp_filt[isim]
            d = dims[isim]
            
            
            
            if d == 0:
                if heights[ih] <= h:
                    ih += 1
                else:
                    break            
            col = []

            if d > 0:
                for fi in comp.faces[d][si]:
                    col.append(sim_to_col[d-1][fi])

            columns.append((d, sorted(col)))

            sim_to_col[d][si] = isim
            isim += 1

            sim_to_pindex[d][si] = pi
            col_to_sim.append(si)
                        
    boundary_matrix = phat.boundary_matrix(columns=columns)
    
    alive = [set() for i in range(comp.dim+1)]
    for i, col in enumerate(columns):
        alive[col[0]].add(i)
        
    raw_pairs = boundary_matrix.compute_persistence_pairs()
    
    ipairs = [{} for i in range(comp.dim+1)]
    
    for (i, j) in raw_pairs:

        d = columns[i][0]
        pi = sim_to_pindex[d][col_to_sim[i]]
        pj = sim_to_pindex[d+1][col_to_sim[j]]
        
        if pi != pj:
            
            if (pi, pj) not in ipairs[d]:
                ipairs[d][(pi, pj)] = 1
            else:
                ipairs[d][(pi, pj)] += 1
            
        alive[d].discard(i)
        alive[d+1].discard(j)
        
    persist = [{} for i in range(comp.dim+1)]
    for d in range(comp.dim+1):
        for i in alive[d]:
            pi = sim_to_pindex[d][col_to_sim[i]]
            if pi not in persist[d]:
                persist[d][pi] = 1
            else:
                persist[d][pi] += 1
    
    
    return (ipairs, hbins, persist, sim_to_pindex)       
        