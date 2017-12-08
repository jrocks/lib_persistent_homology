import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')

import numpy as np
import scipy as sp
import scipy.io as spio

import skimage.morphology as morph

import homology

# from memory_profiler import profile
@profile
def test():
    mat = spio.loadmat("../sample_data/Creases17.mat")

    print(mat['ridges'].shape)
    data = mat['ridges'][2000:2250, 2000:2250]

    Nx = data.shape[0]
    Ny = data.shape[1]

    print(Nx, Ny)

    boundaries = {i+1:[] for i in range(2)}

    for i in range(Ny+1):
        for j in range(Nx):
            boundaries[1].append([(Nx+1)*i + j,(Nx+1)*i + j+1])

    for i in range(Ny):
        for j in range(Nx+1):
            boundaries[1].append([(Nx+1)*i + j,(Nx+1)*(i+1) + j])


    for i in range(Ny):
        for j in range(Nx):
            boundaries[2].append([Ny*i + j, Ny*(Nx+1) + (Ny+1)*i + j+1, Ny*(i+1) + j, Ny*(Nx+1) + (Ny+1)*i + j])


    simplices = []
    dims = []
    heights = []

    curr_height = 0

    num_zeros = np.count_nonzero(data==0)

    print("Height:", curr_height)
    print("Zeros:", num_zeros)

    simplices.extend(np.nonzero(data.flatten())[0])
    dims.extend(np.full_like(simplices, 2))
    heights.extend(np.full_like(simplices, curr_height))

    print(len(simplices))


    selem = morph.disk(4)
    print(selem)

    while num_zeros > 0:

        curr_height += 1


        dilated = morph.dilation(data, selem)

        new_nonzeros = np.nonzero((dilated - data).flatten())[0]

        num_zeros -= len(new_nonzeros)

        print("Height:", curr_height)
        print("Zeros:", num_zeros)

        simplices.extend(new_nonzeros)
        dims.extend(np.full_like(new_nonzeros, 2))
        heights.extend(np.full_like(new_nonzeros, curr_height))

        data = dilated

        if curr_height >= 32:
            break


    (ipairs, pheight, persist, sim_to_pindex) = homology.compute_persistence_pairs(simplices, dims, heights, boundaries)

    print(ipairs)
    print(persist)
    
test()