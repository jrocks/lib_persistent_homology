import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')


import numpy as np
import numpy.linalg as la
import scipy as sp
import scipy.io as spio

import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib as mpl
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.image as mpimg

import itertools as it
import collections as co
import queue
import networkx as nx
import time
import pickle

from skimage import filters as skifilters
from skimage import color as skicolor

from sklearn import linear_model

import homology


def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    dict_handler = lambda d: it.chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    co.deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = sys.getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = sys.getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=sys.stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
                
                
        if not hasattr(o.__class__, '__slots__'):
            if hasattr(o, '__dict__'):
                s+=sizeof(o.__dict__) # no __slots__ *usually* means a __dict__, but some special builtin classes (such as `type(None)`) have neither
            # else, `o` has no attributes at all, so sys.getsizeof() actually returned the correct value
        else:
            s+=sum(sizeof(getattr(o, x)) for x in o.__class__.__slots__ if hasattr(o, x))
                
        return s

    return sizeof(o)

# from memory_profiler import profile
def test():

    mat = spio.loadmat("../notebooks/Everest.mat")
    data = mat['Expression1'][:100, :100]
    print(total_size(mat))

    print(data.shape)
    print(total_size(data))

    dual = True

    print("Constructing complex")
    comp = homology.construct_cubical_complex(data.shape, oriented=False, dual=dual)

    # for c in comp.get_cells():
    #     print(c, comp.get_dim(c), comp.get_facets(c))

    print("Checking boundary operator")
    print(homology.check_boundary_op(comp))

    print("Constructing cofacets")
    comp.construct_cofacets()

    # for c in comp.get_cells():
    #     print(c, comp.get_dim(c), comp.get_cofacets(c))

    print(total_size(comp))


    print("Finding vertex order")

    vertex_time = data.flatten()
    vertex_order = homology.construct_vertex_filtration_order(comp, vertex_time, euclidean=False, positions=None, dual=dual)

    print(total_size(vertex_time))
    print(total_size(vertex_order))

    # print(vertex_order)

    # # fig, ax = plt.subplots(figsize=(16,16))
    # # im = ax.imshow(vertex_order.reshape(data.shape), cmap=plt.cm.Greys_r)

    # # ax.axis('off')
    # # # plt.colorbar(im)
    # # plt.show()

    start = time.time()

    print("Finding Insertion Times")

    insert_order = homology.construct_time_of_insertion_map(comp, vertex_time, vertex_order, dual=dual)

    print(total_size(insert_order))

    # print(insert_order)

    end = time.time()

    print("Elapsed Time:", end - start)

    start = time.time()

    print("Constructing discrete gradient")

    V = homology.construct_discrete_gradient(comp, insert_order, dual=dual)

    print(total_size(V))

    # print(V)

    end = time.time()

    print("Elapsed Time:", end - start)

    n = 0
    for v in range(len(V)):
        if V[v] == v:
            n += 1

    # n = 0
    # for v in V:
    #     if v == V[v]:
    #         n += 1

    print("Number Critical Cells:", n)

    print("Reversing gradient")
    coV = homology.reverse_discrete_gradient(V)

    print(total_size(coV))

    # # print(V)
    # # print(coV)

    print("Calculating Morse complex")
    mcomp = homology.construct_morse_complex(V, comp, oriented=False)

    print("Checking boundary operator")
    print(homology.check_boundary_op(mcomp))

    print("Constructing cofacets")
    mcomp.construct_cofacets()

    # for c in mcomp.get_cells():
    #     print(c, mcomp.get_dim(c), mcomp.get_facets(c))
    # for c in mcomp.get_cells():
    #     print(c, mcomp.get_dim(c), mcomp.get_cofacets(c))

    print(total_size(mcomp))

    start = time.time()

    print("Finding basins")
    basins = homology.find_basins(mcomp, coV, comp, insert_order, dual=dual)
    # print(basins)

    print(total_size(basins))

    print("Calculating Morse skeleton")
    skeleton = homology.find_morse_skeleton(mcomp, V, comp, 1, insert_order, dual=dual)

    print(total_size(skeleton))

    end = time.time()
    print("Elapsed Time:", end - start)


    palette1 = it.cycle(sns.color_palette("deep"))
    palette2 = it.cycle(['Blues_r', 'Greens_r', 'Purples_r', 'Oranges_r', 'RdPu_r'])

    basin_color_map1 = {}
    basin_color_map2 = {}
    for v in basins:
    # for v in mcomp.get_cells():
    #     if mcomp.dims[v]== 0:
        basin_color_map1[v] = next(palette1)
        basin_color_map2[v] = next(palette2)

    print("Complete")
    
    
    print("Calculating persistence pairs...")

    filtration = homology.construct_filtration(mcomp, insert_order)

    # weights = homology.get_morse_weights(mcomp, V, coV, comp.facets, comp.cofacets)

    # print("Weights:", weights)

    start = time.time()

    pairs = homology.compute_persistence(mcomp, filtration)
    # (pairs, bcycles) = homology.compute_persistence(mcomp, filtration, 
    #                                                             birth_cycles=True, optimal_cycles=False)
    # (pairs, bcycles, ocycles) = homology.compute_persistence(mcomp, filtration, 
    #                                                             birth_cycles=True, optimal_cycles=True,
    #                                                                     weights=weights, relative_cycles=True)
    end = time.time()

    print("Elapsed Time:", end - start)

    # print("Pairs:", pairs)
    # print("Birth Cycles:", bcycles)
    # print("Death Cycles:", ocycles)

    print("Complete")
    
    
    persistence = []
    area = []

    for pi, (i, j) in enumerate(pairs):

        if pi % 100 == 0:
            print(pi, "/", len(pairs))

        if j is None:
            persistence.append(np.inf)
            area.append(data.shape[0]*data.shape[1])
        else:
            persistence.append(insert_order[j][0] - insert_order[i][0])
            feature = homology.extract_persistence_feature(i, j, mcomp, comp, V, coV, insert_order)
            pixels = homology.convert_to_pixels(feature, comp, insert_order, dual=dual)
            area.append(len(pixels))

    print("Complete")

    
test()