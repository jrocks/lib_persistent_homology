import sys, os
sys.path.insert(0, '../')
sys.path.insert(0, '../python_src/')
sys.path.insert(0, '../../lib_persistent_homology/')
sys.path.insert(0, '../../lib_persistent_homology/python_src/')
sys.path.insert(0, '../../lib_network_tuning/')
sys.path.insert(0, '../../lib_network_tuning/python_src/')

import numpy as np
import scipy as sp
import pandas as pd

import chomology
import homology
import network_solver as ns




def get_configuration(state, index):
    
    DIM = len(state.dimensions['dim'])

    NP = len(state.dimensions['NP'])

    pos = state.variables['pos'][index]

    rad2 = state.variables['rad'][index]**2

    box_mat = np.array(state.variables['BoxMatrix'][index].reshape((DIM, DIM)), order='F')

    return (NP, pos, rad2, DIM, box_mat)

def get_neighborhood(particle, config, max_neigh_dist):
    
    (NP, pos, rad2, DIM, box_mat) = config
        
    (neighborhood, neigh_pos) = chomology.get_neighborhood(particle, max_neigh_dist, DIM, pos, box_mat)
    
    neigh_pos = np.array(neigh_pos).flatten()
    
    return (len(neighborhood), np.array(neighborhood), np.array(neigh_pos).flatten(), rad2[neighborhood])
     

def get_comp_network(NP, pos, box_mat, comp): 

    NE = 0
    edgei = []
    edgej = []
    for c in range(comp.ncells):
        if comp.get_dim(c) == 1:

            NE += 1

            facets = comp.get_facets(c)
            edgei.append(facets[0])
            edgej.append(facets[1])
            
    
    DIM = 2
    node_pos = np.zeros_like(pos)
    for i in range(NP):
        node_pos[DIM*i:DIM*i+DIM] = box_mat.dot(pos[DIM*i:DIM*i+DIM])
        

    net = ns.Network2D(NP, node_pos, NE, edgei, edgej, box_mat.diagonal())

            
    return net


def get_gap_dist(particles, config, max_tri_dist, use_angular=False, only_neighbor=False, max_neigh_dist=None, verbose=False):
    
    (NP, pos, rad2, DIM, box_mat) = config
    
    r2norm = np.min(rad2)
    
    if only_neighbor:
        
        gap_dist = {}
        
        for p in particles:
        
            (neigh_NP, neighborhood, neigh_pos, neigh_rad2) = get_neighborhood(p, config, max_neigh_dist)

            if DIM == 2:
                neigh_comp = chomology.construct_alpha_complex_2D(neigh_NP, neigh_pos, neigh_rad2, box_mat, periodic=False)
                alpha_vals = chomology.calc_alpha_vals_2D(neigh_pos, neigh_rad2, neigh_comp, box_mat, periodic=False, alpha0=-r2norm)
            elif DIM == 3:
                neigh_comp = chomology.construct_alpha_complex_3D(neigh_NP, neigh_pos, neigh_rad2, box_mat, periodic=False)
                alpha_vals = chomology.calc_alpha_vals_3D(neigh_pos, neigh_rad2, neigh_comp, box_mat, periodic=False, alpha0=-r2norm)
                
            local_p = np.argwhere(neighborhood == p)[0][0]
                             
            if use_angular:
                part_gap_dist = chomology.calc_angular_gap_distribution([local_p], alpha_vals, neigh_comp, max_dist=max_tri_dist)
            else:
                part_gap_dist = chomology.calc_radial_gap_distribution([local_p], alpha_vals, neigh_comp, max_dist=max_tri_dist)
                
                
            gap_dist[p] = part_gap_dist[local_p]
                                    
    else:
                
        if DIM == 2:
            comp = chomology.construct_alpha_complex_2D(NP, pos, rad2, box_mat, periodic=True)
            alpha_vals = chomology.calc_alpha_vals_2D(pos, rad2, comp, box_mat, periodic=True, alpha0=-r2norm)
        elif DIM == 3:
            comp = chomology.construct_alpha_complex_3D(NP, pos, rad2, box_mat, periodic=True)
            alpha_vals = chomology.calc_alpha_vals_3D(pos, rad2, comp, box_mat, periodic=True, alpha0=-r2norm)
        
        if use_angular:
            gap_dist = chomology.calc_angular_gap_distribution(particles, alpha_vals, comp, max_dist=max_tri_dist, verbose=verbose)
        else:
            gap_dist = chomology.calc_radial_gap_distribution(particles, alpha_vals, comp, max_dist=max_tri_dist, verbose=verbose)
        
        
    col_set = set()
    
    if use_angular:
        for p in gap_dist:
            for rad_dist in gap_dist[p][0]:
                for ang_dist in gap_dist[p][0][rad_dist]:
                    col_set.add((ang_dist, rad_dist, 'gg'))

            for rad_dist in gap_dist[p][1]:
                for ang_dist in gap_dist[p][1][rad_dist]:
                    col_set.add((ang_dist, rad_dist, 'oo'))

            for rad_dist in gap_dist[p][2]:
                for ang_dist in gap_dist[p][2][rad_dist]:
                    col_set.add((ang_dist, rad_dist, 'go'))

        col_index = {x:i for i, x in enumerate(sorted(col_set))}

        s_list = []

        for p in gap_dist:

            s = [0]*len(col_index)

            for rad_dist in gap_dist[p][0]:
                for ang_dist in gap_dist[p][0][rad_dist]:
                    s[col_index[(ang_dist, rad_dist, 'gg')]] =  gap_dist[p][0][rad_dist][ang_dist]

            for rad_dist in gap_dist[p][1]:
                for ang_dist in gap_dist[p][1][rad_dist]:
                    s[col_index[(ang_dist, rad_dist, 'oo')]] =  gap_dist[p][1][rad_dist][ang_dist]

            for rad_dist in gap_dist[p][2]:
                for ang_dist in gap_dist[p][2][rad_dist]:
                    s[col_index[(ang_dist, rad_dist, 'go')]] =  gap_dist[p][2][rad_dist][ang_dist]


            s_list.append([p]+s)
            
    else:
        
        for p in gap_dist:
            for rad_dist in gap_dist[p][0]:
                col_set.add((rad_dist, 'g'))

            for rad_dist in gap_dist[p][1]:
                col_set.add((rad_dist, 'o'))

        col_index = {x:i for i, x in enumerate(sorted(col_set))}

        s_list = []

        for p in gap_dist:

            s = [0]*len(col_index)

            for rad_dist in gap_dist[p][0]:
                s[col_index[(rad_dist, 'g')]] =  gap_dist[p][0][rad_dist]

            for rad_dist in gap_dist[p][1]:
                s[col_index[(rad_dist, 'o')]] =  gap_dist[p][1][rad_dist]

            s_list.append([p]+s)
                      
    columns = ['particle'] + sorted(col_set)
    
    df = pd.DataFrame(s_list, columns=columns)   
    
    return df



    
    
    
    
    
    
    
        
    
    
def get_persistence_diag(particle, config,  max_tri_dist, only_neighbor=False, max_neigh_dist=None, weighted=False):
    
    (NP, pos, rad2, DIM, box_mat) = config
    
    r2norm = np.min(rad2)
        
    if only_neighbor:

        (neigh_NP, neighborhood, neigh_pos, neigh_rad2) = get_neighborhood(particle, config, max_neigh_dist)

        if DIM == 2:
            comp = chomology.construct_alpha_complex_2D(neigh_NP, neigh_pos, neigh_rad2, box_mat, periodic=False)
            alpha_vals = chomology.calc_alpha_vals_2D(neigh_pos, neigh_rad2, comp, box_mat, periodic=False, alpha0=-r2norm)
        elif DIM == 3:
            comp = chomology.construct_alpha_complex_3D(neigh_NP, neigh_pos, neigh_rad2, box_mat, periodic=False)
            alpha_vals = chomology.calc_alpha_vals_3D(neigh_pos, neigh_rad2, comp, box_mat, periodic=False, alpha0=-r2norm)
            
        local_p = np.argwhere(neighborhood == particle)[0][0]
        
       
    else:
        
        if DIM == 2:
            comp = chomology.construct_alpha_complex_2D(NP, pos, rad2, box_mat, periodic=True)
            alpha_vals = chomology.calc_alpha_vals_2D(pos, rad2, comp, box_mat, periodic=True, alpha0=-r2norm)
        elif DIM == 3:
            comp = chomology.construct_alpha_complex_3D(NP, pos, rad2, box_mat, periodic=True)
            alpha_vals = chomology.calc_alpha_vals_3D(pos, rad2, comp, box_mat, periodic=True, alpha0=-r2norm)
            
        local_p = particle
    
    
    filt = chomology.construct_filtration(alpha_vals, comp)

    pyfilt = homology.construct_filtration(comp, filt, label=False)

    pairs = homology.compute_persistence(comp, pyfilt)

    

#     weights = np.ones(comp.ncells, float)
#     if weighted:
#         for c in range(neigh_comp.ncells):
#             if neigh_comp.get_dim(c) == 1:
#                 verts = neigh_comp.get_facets(c)

#                 posi = neigh_pos[DIM*verts[0]:DIM*verts[0]+DIM]
#                 posj = neigh_pos[DIM*verts[1]:DIM*verts[1]+DIM]

#                 bvec = posj - posi

#                 weights[c] = la.norm(bvec)
    
    cycles = {}
    for dim in [1, 2]:
#         cycles[dim] = chomology.calc_optimal_cycles(filt, comp, weights, dim=dim, verbose=False)
        
#         print(cycles[dim])

#         cycles[dim] = chomology.calc_optimal_homologous_cycles(filt, comp, weights, dim=dim, verbose=False)
        
#         print(cycles[dim])
        
        
        cycles[dim] = chomology.calc_homologous_birth_cycles(filt, comp, dim=dim)
        
#         print(cycles[dim])
    
    dists = chomology.find_all_tri_distances(local_p, comp, max_tri_dist)

    s_list = []
    for n, (i, j) in enumerate(pairs):
        if j is None:
            s_list.append([(i, j), comp.get_dim(i), 
                           filt.get_time(i) / r2norm, np.inf, np.inf, None, None, None, None])
        else:
            
            if dists[j] > max_tri_dist:
                continue
            
            verts = list(chomology.get_star(j, True, comp, 0))
            
            
            
            type_map = {v: 't' if v == local_p else 's' if neigh_rad2[v] == r2norm else 'l' for v in verts}
                        
            dpart_types = ''.join(sorted(type_map.values(), reverse=True))
            
            dedge_types = []
            

            for e in chomology.get_star(j, True, comp, 1):

                if filt.get_time(e) / r2norm < 0.0:
                    continue

                verts = comp.get_facets(e)
                label = [type_map[v] for v in comp.get_facets(e)]

                dedge_types.append(''.join(sorted(label, reverse=True)))
        
        
            dedge_gaps = ''.join(sorted(dedge_types, reverse=True))
            
            
            cycle_type = (len(dpart_types), dpart_types, len(dedge_gaps)//2, dedge_gaps)
        
            if comp.get_dim(i) in cycles:
                
                bcycle = set(cycles[comp.get_dim(i)][i])
                dcycle = set(comp.get_facets(j))
                
                bsize = len(bcycle)
                
                bdcycle_equal = (bcycle == dcycle)
                
                
            else:
                bsize = None
                bdcycle_equal = None
                
            s_list.append([(i, j), comp.get_dim(i), 
                           filt.get_time(i) / r2norm, filt.get_time(j) / r2norm, 
                           np.abs(filt.get_time(i)-filt.get_time(j)) / r2norm, dists[j],
                           cycle_type, bsize, bdcycle_equal])
            
            
            
#             d = neigh_comp.get_dim(i)
            
#             if d == 1 or d == 2:
        
            
# #                 bcycle = []
# #                 bdist = 1e6
# #                 for c in cycles[d][i]:
# #                     if dists[c] < bdist:
# #                         bdist = dists[c]
                    
# #                 bcycle = neigh_comp.get_labels(set(cycles[d][i]))
                    
#                 dcycle = neigh_comp.get_facets(j)
#                 dcycle = neigh_comp.get_labels(set(dcycle))
                    
        
#                 s_list.append([ptype, seed, index, particle, (i, j), neigh_comp.get_dim(i), 
#                                filt.get_time(i) / r2norm, filt.get_time(j) / r2norm, 
#                                np.abs(filt.get_time(i)-filt.get_time(j)) / r2norm,
#                                None, None, None, None, None, dists[j]])
# #                                 bcycle, len(bcycle), bdist, dcycle, len(dcycle), dists[j]])
#             else:
        
        
#                 feature = chomology.extract_persistence_feature(i, j, neigh_comp, filt, target_dim=0)
                            
#                 min_dist = 1e6
#                 for c in feature:
#                     if dists[c] < min_dist:
#                         min_dist = dists[c]
        
#                 s_list.append([ptype, seed, index, particle, (i, j), neigh_comp.get_dim(i), 
#                                filt.get_time(i) / r2norm, filt.get_time(j) / r2norm, 
#                                np.abs(filt.get_time(i)-filt.get_time(j)) / r2norm,
#                                None, len(feature), min_dist, None, None, dists[j]])
            
    
#     end = time.time()
#     print("filling dataframe", end-start)
    
    
    columns = ['pair', 'dim', 'birth', 'death', 'persistence', 'ddist', 'cycle_type', 'bsize', 'bdcycle_equal']
    
    df = pd.DataFrame(s_list, columns=columns)   
    
    return df





def get_cycle_dist(particles, config, max_tri_dist, only_neighbor=False, max_neigh_dist=None, verbose=False):
    
    (NP, pos, rad2, DIM, box_mat) = config
    
    r2norm = np.min(rad2)
    
    if only_neighbor:
        
        cycle_dist = {}
        
        for p in particles:
        
            (neigh_NP, neighborhood, neigh_pos, neigh_rad2) = get_neighborhood(p, config, max_neigh_dist)

            if DIM == 2:
                neigh_comp = chomology.construct_alpha_complex_2D(neigh_NP, neigh_pos, neigh_rad2, box_mat, periodic=False)
                alpha_vals = chomology.calc_alpha_vals_2D(neigh_pos, neigh_rad2, neigh_comp, box_mat, periodic=False, alpha0=-r2norm)
            elif DIM == 3:
                neigh_comp = chomology.construct_alpha_complex_3D(neigh_NP, neigh_pos, neigh_rad2, box_mat, periodic=False)
                alpha_vals = chomology.calc_alpha_vals_3D(neigh_pos, neigh_rad2, neigh_comp, box_mat, periodic=False, alpha0=-r2norm)
                
            local_p = np.argwhere(neighborhood == p)[0][0]
                      
            vtypes = np.where(neigh_rad2 == r2norm, 's', 'l')

            part_cycle_dist = chomology.calc_radial_cycle_distribution([local_p], alpha_vals, vtypes, neigh_comp, max_dist=max_tri_dist)
                
                
            cycle_dist[p] = part_cycle_dist[local_p]
                                    
    else:
                
        if DIM == 2:
            comp = chomology.construct_alpha_complex_2D(NP, pos, rad2, box_mat, periodic=True)
            alpha_vals = chomology.calc_alpha_vals_2D(pos, rad2, comp, box_mat, periodic=True, alpha0=-r2norm)
        elif DIM == 3:
            comp = chomology.construct_alpha_complex_3D(NP, pos, rad2, box_mat, periodic=True)
            alpha_vals = chomology.calc_alpha_vals_3D(pos, rad2, comp, box_mat, periodic=True, alpha0=-r2norm)
        
   
        vtypes = np.where(rad2 == r2norm, 's', 'l')
    
        cycle_dist = chomology.calc_radial_cycle_distribution(particles, alpha_vals, vtypes, comp, max_dist=max_tri_dist, verbose=verbose)
        
        
        
    print(cycle_dist)
    col_set = set()
    

    for p in cycle_dist:
        for r in cycle_dist[p]:
            for cycle_type in cycle_dist[p][r]:
                
                if len(cycle_type[0]) != 3:
                    continue
                
                col_set.add((r, len(cycle_type[0]), cycle_type[0], len(cycle_type[1])//2, cycle_type[1]))

    col_index = {x:i for i, x in enumerate(sorted(col_set))}
    
    s_list = []

    for p in cycle_dist:

        s = [0]*len(col_index)

        for r in cycle_dist[p]:
            for cycle_type in cycle_dist[p][r]:
                
                if len(cycle_type[0]) != 3:
                    continue
                
                s[col_index[(r, len(cycle_type[0]), cycle_type[0], len(cycle_type[1])//2, cycle_type[1])]] =  cycle_dist[p][r][cycle_type]

        s_list.append([p]+s)
                      
    columns = ['particle'] + sorted(col_set)
        
    df = pd.DataFrame(s_list, columns=columns)   
    
    return df
