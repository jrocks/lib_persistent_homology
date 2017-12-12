import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
import numpy as np
import scipy as sp
import phat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import skimage.morphology as morph
import networkx as nx

class CellComplex:
    
    def __init__(self, dim, nverts):
        self.dim = dim
        self.nverts = nverts
        self.faces = {i+1:[] for i in range(dim)}
        
    def add(self, dim, cell):
        self.faces[dim].append(cell)
        
        
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
# Does not need to include all possible vertices.
def construct_lower_star_filtration(comp, verts):
    
    cofaces = {}
    cofaces[0] = [[] for i in range(comp.nverts)]
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
        
    
    
def mesh_complex(nrows, ncols, compactify=False):
    
    nverts = nrows*ncols
    if compactify:
        nverts += 1
    comp = CellComplex(2, nverts)

    for i in range(nrows):
        for j in range(ncols-1):
            comp.add(1, [ncols*i + j, ncols*i + j+1])

    for i in range(nrows-1):
        for j in range(ncols):
            comp.add(1, [ncols*i + j, ncols*(i+1) + j])
            
    for i in range(nrows-1):
        for j in range(ncols-1):
            comp.add(2, [(ncols-1)*i + j, 
                         (ncols-1)*nrows + ncols*i + j+1, 
                         (ncols-1)*(i+1) + j, 
                         (ncols-1)*nrows + ncols*i + j])
            
    if compactify:
        for j in range(ncols-1):
            comp.add(1, [j, nrows*ncols])
        
        for i in range(nrows-1):
            comp.add(1, [ncols*i + ncols-1, ncols*nrows])
            
        for j in range(ncols-1, 0, -1):
            comp.add(1, [ncols*(nrows-1) + j, ncols*nrows])
            
        for i in range(nrows-1, 0, -1):
            comp.add(1, [ncols*i, ncols*nrows])
        
        for j in range(ncols-1):
            comp.add(2, [j, 
                         (ncols-1)*nrows + ncols*(nrows-1) + j, 
                         (ncols-1)*nrows + ncols*(nrows-1) + j+1])
            
        for i in range(nrows-1):
            comp.add(2, [(ncols-1)*nrows + ncols-1 + ncols*i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + i+1])
        
        for j in range(ncols-1):
            comp.add(2, [(ncols-1)*nrows - 1 - j , 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + j, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + j+1])
        
        for i in range(nrows-2):
            comp.add(2, [(ncols-1)*nrows + ncols*(nrows-1) - ncols - ncols*i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + i, 
                         (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + i+1])
        
        comp.add(2, [(ncols-1)*nrows + ncols*(nrows-1) - ncols - ncols*(nrows-2), 
                     (ncols-1)*nrows + ncols*(nrows-1) + ncols-1 + nrows-1 + ncols-1 + nrows-2, 
                     (ncols-1)*nrows + ncols*(nrows-1)])
            
    return comp




def binary_dilation_filtration(mat, show=False):
    
    pixels = []
    heights = []

    curr_height = 0

    num_zeros = np.count_nonzero(mat==0)

    print("Height:", curr_height)
    print("Zeros:", num_zeros)

    pixels.extend(np.nonzero(mat.flatten())[0])
    heights.extend(np.full_like(pixels, curr_height))

    print(len(pixels))

    if show:
        fig, ax = plt.subplots(figsize=(6,6))
        ax.imshow(mat, cmap=plt.cm.gray)
        # ax.axis('off')
        plt.show()


    max_radius = 4

    base = mat
    prev_dilated = mat
    radius = 1

    while num_zeros > 0:

        selem = morph.disk(radius)

    #     print(selem)

        curr_height += 1

        dilated = morph.dilation(base, selem)

        new_nonzeros = np.nonzero((dilated - prev_dilated).flatten())[0]

        num_zeros -= len(new_nonzeros)

        print("Height:", curr_height, "Zeros:", num_zeros, "Radius", radius)

        if show:
            fig, ax = plt.subplots(figsize=(6,6))
            ax.imshow(dilated, cmap=plt.cm.gray)
            # ax.axis('off')
            plt.show()

        pixels.extend(new_nonzeros)
        heights.extend(np.full_like(new_nonzeros, curr_height))

        prev_dilated = dilated

        if radius == max_radius:
            base = dilated
            radius = 1
        else:
            radius += 1

    #     if curr_height >= 40:
    #         break

    pixels.append(mat.shape[0]*mat.shape[1])
    heights.append(curr_height+1)
    
    return (pixels, heights)

def laplace_dilation_filtration(mat, show=False):
    
    image = np.zeros_like(mat, float)

    i = 0
    while True:

        new_image = np.copy(mat.astype(float))

        # bulk
        new_image[1:mat.shape[0]-1, 1:mat.shape[1]-1] += (image[0:mat.shape[0]-2, 1:mat.shape[1]-1] 
                                                          + image[1:mat.shape[0]-1, 0:mat.shape[1]-2] 
                                                          + image[2:mat.shape[0], 1:mat.shape[1]-1] 
                                                          + image[1:mat.shape[0]-1, 2:mat.shape[1]])

        # top row
        new_image[0, 1:mat.shape[1]-1] += (image[0, 0:mat.shape[1]-2] 
                                          + image[1, 1:mat.shape[1]-1] 
                                          + image[0, 2:mat.shape[1]])

        # bottom row
        new_image[mat.shape[0]-1, 1:mat.shape[1]-1] += (image[mat.shape[0]-2, 1:mat.shape[1]-1] 
                                                          + image[mat.shape[0]-1, 0:mat.shape[1]-2] 
                                                          + image[mat.shape[0]-1, 2:mat.shape[1]])

        # left column
        new_image[1:mat.shape[0]-1, 0] += (image[0:mat.shape[0]-2, 0] 
                                                          + image[2:mat.shape[0], 0] 
                                                          + image[1:mat.shape[0]-1, 1])

        # right column
        new_image[1:mat.shape[0]-1, mat.shape[1]-1] += (image[0:mat.shape[0]-2, mat.shape[1]-1] 
                                                          + image[1:mat.shape[0]-1, mat.shape[1]-2] 
                                                          + image[2:mat.shape[0], mat.shape[1]-1])

        # upper left corner
        new_image[0, 0] += (image[1, 0] 
                              + image[0, 1])

        # upper right corner
        new_image[0, mat.shape[1]-1] += (image[0, mat.shape[1]-2] 
                                          + image[1, mat.shape[1]-1])

        # bottorm left corner
        new_image[mat.shape[0]-1, 0] += (image[mat.shape[0]-2, 0] 
                                            + image[mat.shape[0]-1, 1])

        # bottom right corner
        new_image[mat.shape[0]-1, mat.shape[1]-1] += (image[mat.shape[0]-2, mat.shape[1]-1] 
                                                          + image[mat.shape[0]-1, mat.shape[1]-2])

        new_image /= 6.0

        

        max_delta = np.max(np.abs((image-new_image) / np.where(image==0, 1e-16, np.abs(image))))


        if show and i %10 == 0:
            print(i, max_delta)
            
            # norm = mcolors.LogNorm(vmin=1e-16, vmax=1.0)
            # cmap = mpl.cm.GnBu
            
            norm = mcolors.SymLogNorm(vmin=-1.0, vmax=1.0, linthresh=1e-32)
            # cmap = mpl.cm.RdBu_r
            cmap = mpl.cm.RdYlGn
            
            fig, ax = plt.subplots(figsize=(8,8))
            im = ax.imshow(new_image, cmap=cmap, norm=norm)
            # plt.colorbar(im)
            # ax.axis('off')
            plt.show()

        
        if max_delta <= 1e-4:
            image = new_image
            break

        image = new_image
        i += 1
        
    heights = list(np.sort(image.flatten()))
    pixels = list(np.argsort(image.flatten()))
    
    pixels.append(mat.shape[0]*mat.shape[1])
    heights.append(heights[-1]+1)
    
    return (pixels, heights)


# def triangulate_complex(comp, vert_order):
    
#     cofaces = {}
#     cofaces[0] = [[] for i in range(comp.nverts)]
#     for d in range(1, comp.dim):
#         cofaces[d] = [[] for i in range(len(comp.faces[d]))]
    
#     for d in range(1, comp.dim+1):
#         for si, simplex in enumerate(comp.faces[d]):
#             for sj in simplex:
#                 cofaces[d-1][sj].append(si)
                
#     for vi in vert_order:
#         for ei in cofaces[0][vi]:
        

# def compute_morse_complex(comp, vert_order):
    
#     cofaces = {}
#     cofaces[0] = [[] for i in range(comp.nverts)]
#     for d in range(1, comp.dim):
#         cofaces[d] = [[] for i in range(len(comp.faces[d]))]
    
#     for d in range(1, comp.dim+1):
#         for si, simplex in enumerate(comp.faces[d]):
#             for sj in simplex:
#                 cofaces[d-1][sj].append(si)
                
#     heights = {}
#     for i, vi in enumerate(vert_order):
#         heights[vi] = i
        
        
#     paths = []
#     end_points = {}
            
#     for vi in vert_order:
        
#         lstar = {d:{} for d in range(1, comp.dim+1)}
#         ustar = {d:{} for d in range(1, comp.dim+1)}
        
#         for ei in cofaces[0][vi]:
#             vj = (set(comp.faces[1][ei]) - set([vi])).pop()
            
#             if vj in heights:
#                 if heights[vj] < heights[vi]:
#                     lstar[1][ei] =  [vj]
#                 else:
#                     ustar[1][ei] =  [vj]
                
#         for d in range(1, comp.dim):
#             for si in lstar[d]:
#                 for sj in cofaces[d][si]:
#                     if sj not in lstar[d+1]:
#                         lstar[d+1][sj] = [si]
#                     else:
#                         lstar[d+1][sj].append(si)
                        
#             for si in ustar[d]:
#                 for sj in cofaces[d][si]:
#                     if sj not in ustar[d+1]:
#                         ustar[d+1][sj] = [si]
#                     else:
#                         ustar[d+1][sj].append(si)
           
                
#         for d in range(comp.dim, 1, -1):
            
#             while len(lstar[d]) > 0:
#                 c, c_list = lstar[d].popitem()
                
#                 si = c_list.pop()
#                 while len(c_list) > 0:
#                     sj = c_list.pop()
#                     if sj in lstar[d-1]:
#                         sj_list = lstar[d-1].pop(sj)
#                         lstar[d-1][si].extend(sj_list)
                 
#             lstar.pop(d)
            
#             while len(ustar[d]) > 0:
#                 c, c_list = ustar[d].popitem()
                
#                 si = c_list.pop()
#                 while len(c_list) > 0:
#                     sj = c_list.pop()
#                     if sj in ustar[d-1]:
#                         sj_list = ustar[d-1].pop(sj)
#                         ustar[d-1][si].extend(sj_list)
                 
#             ustar.pop(d)
            
            
            
            
#         lower = list(lstar[1].values())
        
#         upper = list(ustar[1].values())
                
#         print(lower, upper)
        
        
        
#         # minimum
#         if len(lower) == 0:
#             print(vi, "min")
#             pass
        
#         # regular
#         elif len(lower) == 1 and len(upper) == 1:
#             print(vi, "regular")
#             pass
        
#         # maximum
#         elif len(upper) == 0:
#             print(vi, "max")
#             pass
        
#         #saddle
#         else:
#             print(vi, "saddle")
#             pass
        
        
            
        