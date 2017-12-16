import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import numpy.linalg as la
import phat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import skimage.morphology as morph
import networkx as nx
import queue


class CellComplex:
    
    def __init__(self, dim):
        self.dim = dim
        self.faces = [{} for i in range(dim+1)]
        
    def add(self, dim, label, cell):
        self.faces[dim][label] = cell
        
        
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
  

# comp - simplicial complex
# verts - ordered list of vertices used in filtration. 
# heights - list of heights for all verts
# Does not need to include all possible vertices.
# def construct_lower_star_filtration(comp, verts, heights):
    
#     cofaces = {}
#     cofaces[0] = [[] for i in range(comp.nverts)]
#     for d in range(1, comp.dim):
#         cofaces[d] = [[] for i in range(len(comp.faces[d]))]
    
#     for d in range(1, comp.dim+1):
#         for si, simplex in enumerate(comp.faces[d]):
#             for sj in simplex:
#                 cofaces[d-1][sj].append(si)
    
#     simp_filt = []
#     dims = []
#     filt_heights = []
    
#     marked = {i: set() for i in range(comp.dim+1)}
    
#     # go through vertices from lowest to highest
#     for vi in verts:
        
#         curr_simps = set([vi])
#         new_simps = set()
        
#         simp_filt.append(vi)
#         dims.append(0)
#         filt_heights.append(heights[vi])
#         marked[0].add(vi)
                
#         for d in range(comp.dim):
#             for cfi in curr_simps:
#                 for cfj in cofaces[d][cfi]:
#                     faces = set(comp.faces[d+1][cfj])
#                     if cfj not in marked[d+1] and faces <= marked[d]:
#                         simp_filt.append(cfj)
#                         dims.append(d+1)
#                         filt_heights.append(heights[vi])
#                         marked[d+1].add(cfj)
#                         new_simps.add(cfj)
                                      
#             curr_simps = new_simps
#             new_simps = set()
     
#     return (simp_filt, dims, filt_heights)


def compute_graph_segmentation(G, heights, euclidean=False, positions=None):
    
#     print("Finding unique heights")
    
    unique = np.unique(heights)
    
#     print("Finding Level Sets")
#     # find level sets
#     level_sets = []
#     for i, h in enumerate(unique):
#         if i % 10000 == 0:
#             print("Level:", i, "/", len(unique))
#         level_sets.append(set(np.where(heights==h)[0]))
    
    # basin 0 is the watersheds
    basins = np.full(G.order(), -1, int)
    revisit = set()
    
    current_label = 0
    
    min_dist = 1.0
     
    print("Sorting heights")
        
    # iterate through each level set
    # for ih in range(len(level_sets)):
    
    arg_heights = np.argsort(heights)
    
    ih = 0
    h = heights[arg_heights[ih]]
    
    print("Iterating through levels")
    
    while ih < len(arg_heights):
        
        unvisited = revisit.copy()
        
        while True:
            
            if ih % 100000 == 0:
                print("Level:", ih, "/", len(arg_heights), h)
            
            h = heights[arg_heights[ih]]
            unvisited.add(arg_heights[ih])
            
            ih += 1
            
            
            if ih >= len(arg_heights) or h != heights[arg_heights[ih]]:
                break
            
                
        
        # level = level_sets[ih]
        
        dist = {}
        
        if euclidean:
            Q = queue.PriorityQueue()
            closest = {}
        else:
            Q = queue.Queue()            
        
        
        
        
        # if ih % 10000 == 0:
        #     print("Level:", ih, "/", len(level_sets), unique[ih])
        #     print("Visiting:", len(unvisited))   

            
        # find neighbors in lower level sets
        for vi in unvisited:
            
            for nbr in G[vi]:
                if basins[nbr] > 0:
                    
                    if euclidean:
                        new_dist = la.norm(positions[vi] - positions[nbr])                    
                    else:
                        new_dist = 1.0
                    
                    if vi not in dist or new_dist < dist[vi]:
                        dist[vi] = new_dist
                        Q.put((new_dist, vi))
                        
                        if euclidean:
                            closest[vi] = nbr
                                    
        # breadth first search through connected part of level set
        while not Q.empty():
                        
            (current_dist, vi) = Q.get()    
            
            if vi not in unvisited:
                continue
            
            for nbr in G[vi]:
                
                if euclidean:
                    new_dist = la.norm(positions[nbr] - positions[closest[vi]])
                else:
                    new_dist = current_dist + 1                        
                        
                # neighbor has already been assigned to a basin
                # 1. has been assigned a basin in a previous level sweep (not in dist)
                # 2. assigned a basin because it has a smaller distance (dist[nbr] < current_dist)
                if basins[nbr] > 0 and (nbr not in dist or dist[nbr] < current_dist):
                    
                    # haven't visited this node yet
                    if vi in unvisited:
                        # set to neighbor's basin and mark as visited
                        basins[vi] = basins[nbr]
                        unvisited.discard(vi)
                        
                        # if was a watershed then remove from watershed
                        revisit.discard(vi)
                        
                    # have already visited but the neighbor is in a different basin
                    # if already in watershed, this does nothing
                    elif basins[vi] != basins[nbr]:
                        basins[vi] = 0
                        
                        if dist[vi] > min_dist:
                            revisit.add(vi)
                
                # neighbor has not been visited and not been given a distance
                elif nbr in unvisited and (nbr not in dist or new_dist < dist[nbr]):
                    # append to queue
                    
                    dist[nbr] = new_dist
                    Q.put((new_dist, nbr))
                    
                    if euclidean:
                        closest[nbr] = closest[vi]
            
            
        # breadth first search through separate components corresponding to new minima
        while len(unvisited) > 0:
                        
            vi = unvisited.pop()
            current_label += 1
            basins[vi] = current_label
            
            Q.put((0,vi))
            
            while not Q.empty():
                (d, vj) = Q.get()
                
                for nbr in G[vj]:
                    if nbr in unvisited:
                        Q.put((0,nbr))
                        basins[nbr] = current_label
                        unvisited.discard(nbr)

    
    segments = {i:set() for i in range(current_label+1)}
    for i in range(G.order()):
        segments[basins[i]].add(i)
        
    return segments
    

def construct_mesh_graph(nrows, ncols, diagonals=False, pos=False):
    G = nx.Graph()
    
    positions = {}
    
    for i in range(nrows):
        for j in range(ncols):
            G.add_node(ncols*i + j)
            if pos:
                positions[ncols*i + j] = np.array([j, i])
            
            
    for i in range(nrows):
        for j in range(ncols-1):
            G.add_edge(ncols*i + j, ncols*i + j+1)

    for i in range(nrows-1):
        for j in range(ncols):
            G.add_edge(ncols*i + j, ncols*(i+1) + j)
         
    if diagonals:
        for i in range(nrows-1):
            for j in range(ncols-1):
                G.add_edge(ncols*i + j, ncols*(i+1) + j+1)

        for i in range(nrows-1):
            for j in range(ncols-1):
                G.add_edge(ncols*(i+1) + j, ncols*i + j+1)
    
    if pos:
        return (G, positions)
    else:
        return G
    
def construct_mesh_complex(nrows, ncols, compactify=False):
    
    nverts = nrows*ncols
    if compactify:
        nverts += 1
    comp = CellComplex(2)

    for i in range(nverts):
        comp.add(0, i, [])
    
    iedge = 0
    for i in range(nrows):
        for j in range(ncols-1):
            comp.add(1, iedge, [ncols*i + j, ncols*i + j+1])
            iedge += 1

    for i in range(nrows-1):
        for j in range(ncols):
            comp.add(1, iedge, [ncols*i + j, ncols*(i+1) + j])
            iedge += 1
            
    iface = 0
    for i in range(nrows-1):
        for j in range(ncols-1):
            comp.add(2, iface, [(ncols-1)*i + j, 
                         (ncols-1)*nrows + ncols*i + j+1, 
                         (ncols-1)*(i+1) + j, 
                         (ncols-1)*nrows + ncols*i + j])
            iface += 1
            
    if compactify:
        for j in range(ncols-1):
            comp.add(1, iedge, [j, nrows*ncols])
            iedge += 1
        
        for i in range(nrows-1):
            comp.add(1, iedge, [ncols*i + ncols-1, ncols*nrows])
            iedge += 1
            
        for j in range(ncols-1, 0, -1):
            comp.add(1, iedge, [ncols*(nrows-1) + j, ncols*nrows])
            iedge += 1
            
        for i in range(nrows-1, 0, -1):
            comp.add(1, iedge, [ncols*i, ncols*nrows])
            iedge += 1
        
        for j in range(ncols-1):
            comp.add(2, iface, [j, 
                         (ncols-1)*nrows + ncols*(nrows-1) + j, 
                         (ncols-1)*nrows + ncols*(nrows-1) + j+1])
            iface += 1
            
        for i in range(nrows-1):
            comp.add(2, iface, [(ncols-1)*nrows + ncols-1 + ncols*i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + i+1])
            iface += 1
        
        for j in range(ncols-1):
            comp.add(2, iface, [(ncols-1)*nrows - 1 - j , 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + j, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + j+1])
        
        iface += 1
        for i in range(nrows-2):
            comp.add(2, iface, [(ncols-1)*nrows + ncols*(nrows-1) - ncols - ncols*i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + i+1])
            iface += 1
        
        comp.add(2, iface, [(ncols-1)*nrows + ncols*(nrows-1) - ncols - ncols*(nrows-2), 
                     (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + nrows-2, 
                     (ncols-1)*nrows + ncols*(nrows-1)])
            
    return comp



def compute_graph_dilation(G, sources, euclidean=False, positions=None):
    dist = np.full(G.order(), -1, float)
    
    if euclidean:
        Q = queue.PriorityQueue()
        closest = np.full(G.order(), -1, int)
    else:
        Q = queue.Queue()      
    
    unvisited = set(np.arange(G.order()))
    
    for i in sources:
        dist[i] = 0
        Q.put((0, i))
        
        if euclidean:
            closest[i] = i
        
        
    while not Q.empty():
        
        (current_dist, vi) = Q.get()
        if vi not in unvisited:
            continue
            
        unvisited.discard(vi)
                
        for nbr in G[vi]:
            
            if nbr in unvisited:
                
                if euclidean:
                    new_dist = la.norm(positions[nbr] - positions[closest[vi]])
                else:
                    new_dist = current_dist + 1
                                                                
                if nbr in unvisited and (dist[nbr] == -1 or new_dist < dist[nbr]):
                    dist[nbr] = new_dist 
                    Q.put((new_dist, nbr))   
                    
                    if euclidean:
                        closest[nbr] = closest[vi]
                                    
    return dist

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

        