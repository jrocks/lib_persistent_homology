import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import collections as collections
import matplotlib.patches as patches


import sys
sys.path.insert(0, '../../lib_persistent_homology/')

import phom


gray = "#969696"


def show_network(ax, comp, embed, styles={}, alpha=1.0, boundary_cutoff=0.1, zorder=None, kwargs=dict()):
    
    box_mat = embed.box_mat
    L = np.diagonal(box_mat)
        
    image_offsets = [np.array([0.0, 0.0]),
                    np.array([1.0, 0.0]),
                    np.array([-1.0, 0.0]),
                    np.array([0.0, 1.0]),
                    np.array([0.0, -1.0]),
                    np.array([1.0, 1.0]),
                    np.array([-1.0, 1.0]),
                    np.array([1.0, -1.0]),
                    np.array([-1.0, -1.0])]
    
    edges = []
    edge_index = []
    for c in range(*comp.dcell_range[1]):
        
        (vi, vj) = comp.get_facets(c)
        
        vposi = embed.get_vpos(vi)
        vposj = embed.get_vpos(vj)
        
        vbvec = embed.get_vdiff(vposi, vposj)
       
            
        posi = box_mat.dot(vposi) / L
        bvec = box_mat.dot(vbvec) / L
        posj = posi + bvec
        
        if embed.periodic:
                        
            test_duplicates = ((posi < 0.0).any() or (posi > 1.0).any()
                              or (posj < 0.0).any() or (posj > 1.0).any())
            
            
            if test_duplicates:
                for offset in image_offsets:
                    oposi = box_mat.dot(vposi+offset) / L
                    oposj = oposi+bvec
                    
                    if ((oposi > -boundary_cutoff).all() and (oposi < 1.0+boundary_cutoff).all()
                        and (oposj > -boundary_cutoff).all() and (oposj < 1.0+boundary_cutoff).all()):
                        edges.append([tuple(oposi),tuple(oposj)])
                        edge_index.append(comp.get_label(c))
                    
            else:
                edges.append([tuple(posi),tuple(posj)])
                edge_index.append(comp.get_label(c))
                
        else:
            edges.append([tuple(posi),tuple(posj)])
            edge_index.append(comp.get_label(c))
                
        
        
        
    ls = []
    colors = []
    lw = []
    
    for i, b in enumerate(edge_index):
        
        if b in styles and 'color' in styles[b]:
            colors.append(styles[b]['color'])
        else:
            colors.append(gray)
            
            
        if b in styles and 'ls' in styles[b]:
            ls.append(styles[b]['ls'])
        else:
            ls.append('solid')
            
        if b in styles and 'lw' in styles[b]:
            lw.append(styles[b]['lw'])
        else:
            lw.append(2.0)
            
            
            
    lc = collections.LineCollection(edges, linestyle=ls, lw=lw, alpha=alpha, color=colors, zorder=zorder, **kwargs)
    ax.add_collection(lc)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
def show_discs(ax, comp, embed, rad, subset=None, styles={}, alpha=1.0, boundary_cutoff=0.01, zorder=None, edgecolor='k', kwargs=dict()):
    
    box_mat = embed.box_mat
    L = np.diagonal(box_mat)
    
    image_offsets = [np.array([0.0, 0.0]),
                    np.array([1.0, 0.0]),
                    np.array([-1.0, 0.0]),
                    np.array([0.0, 1.0]),
                    np.array([0.0, -1.0]),
                    np.array([1.0, 1.0]),
                    np.array([-1.0, 1.0]),
                    np.array([1.0, -1.0]),
                    np.array([-1.0, -1.0])]    
    
    discs = []
    disc_index = []
    
    if subset is not None:
        cell_list = subset
    else:
        cell_list = range(*comp.dcell_range[0])
    
    for c in cell_list:
        
        vi = comp.get_label(c)
        
        vposi = embed.get_vpos(vi)       

        posi = box_mat.dot(vposi) / L
                    
        r = rad[vi] / L[0]
        
        
        if embed.periodic:
            
            test_duplicates = ((posi < r+boundary_cutoff).any() or (posi > 1.0-r-boundary_cutoff).any())
            
            if test_duplicates:
                for offset in image_offsets:
                    oposi = box_mat.dot(vposi+offset) / L
                    
                    if ((oposi > -r-boundary_cutoff).all() and (oposi < 1.0+r+boundary_cutoff).all()):
                        discs.append(patches.Circle(oposi, r))
                        disc_index.append(vi)
                    
            else:
                discs.append(patches.Circle(posi, r))
                disc_index.append(vi)
                
        else:
            discs.append(patches.Circle(posi, r))
            disc_index.append(vi)
        
        
        
        
    colors = []
    
    for i, b in enumerate(disc_index):
        
        if b in styles and 'color' in styles[b]:
            colors.append(styles[b]['color'])
        else:
            colors.append('white')
            
            
            
    pc = collections.PatchCollection(discs, edgecolor=edgecolor, linewidth=0.5, alpha=alpha, facecolors=colors, zorder=zorder, **kwargs)
    ax.add_collection(pc)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    
def show_vec_field(ax, comp, embed, vec_field, zorder=None, color='k', kwargs=dict()):
    
    L = np.diagonal(embed.box_mat)
    DIM = embed.dim
    
    X = []
    Y = []
    U = []
    V = []
    
    for c in range(*comp.dcell_range[0]):
        
        vi = comp.get_label(c)
        
        pos = embed.get_pos(vi) / L
        
        u = vec_field[DIM*vi:DIM*vi+DIM] / L
        
        X.append(pos[0])
        Y.append(pos[1])
        U.append(u[0])
        V.append(u[1])
        
        
    ax.quiver(X, Y, U, V, units='xy', scale=1.0, width=0.005, zorder=None, color=color, **kwargs)
        


