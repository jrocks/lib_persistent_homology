import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import phat

# simplices - main simplices in filtration. If the boundary of a simplex
# is not specified, it is automatically added at the same level of the
# filtration.
# 
# dimensions - the corresponding list of lists of simplex dimensions
# 
# heights - heights used for filtration
# 
# boundaries - boundaries for each possible simplex
# 
# returns lists of persistence pairs
def compute_persistence_pairs(simplices, dimensions, heights, boundaries):
    
    icol = 0
    columns = []
    
    sim_to_col = [{} for d in range(len(boundaries)+1)]
    col_to_sim = []
    
    sim_to_pindex = [{} for d in range(len(boundaries)+1)]
        
    # go through groups of simplices in order from lowest to highest
    hsort = np.unique(heights)
    for pi, h in enumerate(hsort):
                
        # this can be optimized
        sindex = np.where(heights==h)[0]
        
        # iterate through each simplex with same height
        for si in sindex:
            
            sim = simplices[si]
            dim = dimensions[si]
            
            if sim in sim_to_col[dim]:
                return "Error"
                        
            # tabulate all simplices to check at each dimension
            all_sim = [set() for d in range(dim+1)]
            all_sim[dim].add(sim)    
            for d in range(dim, 0, -1):
                for sj in all_sim[d]:
                    all_sim[d-1].update(boundaries[d][sj])
                               
            # add vertices (simplices without boundary)
            for vi in all_sim[0]:
                if vi not in sim_to_col[0]:
                    columns.append((0, []))
                    
                    sim_to_col[0][vi] = icol
                    icol += 1
                    
                    sim_to_pindex[0][vi] = pi
                    col_to_sim.append(vi)
                    
            # add simplices with boundary in order from lowest to highest dimension
            # boundaries are guarenteed to always exist
            for d in range(1, dim+1):
                for sj in all_sim[d]:
                    if sj not in sim_to_col[d]:

                        col = []
                        for bi in boundaries[d][sj]:
                            col.append(sim_to_col[d-1][bi])
                            
                        columns.append((d, sorted(col)))
                        
                        sim_to_col[d][sj] = icol
                        icol += 1
                        
                        sim_to_pindex[d][sj] = pi
                        col_to_sim.append(sj)
                    
        
    boundary_matrix = phat.boundary_matrix(columns=columns)

    alive = [set() for i in range(len(boundaries)+1)]
    for col in boundary_matrix.columns:
        alive[col.dimension].add(col.index)
    
    raw_pairs = boundary_matrix.compute_persistence_pairs()

    ipairs = [{} for i in range(len(boundaries)+1)]
    
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
        
    persist = [{} for i in range(len(boundaries)+1)]
    for d in range(len(boundaries)+1):
        for i in alive[d]:
            pi = sim_to_pindex[d][col_to_sim[i]]
            if pi not in persist[d]:
                persist[d][pi] = 1
            else:
                persist[d][pi] += 1
    
    
    return (ipairs, hsort, persist, sim_to_pindex)
