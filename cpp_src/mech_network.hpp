#ifndef MECHNET_HPP
#define MECHNET_HPP

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <utility>


#include "eigen_macros.hpp"  
#include "embedding.hpp"
#include "cell_complex.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;



template <int DIM> std::tuple<CellComplex, Embedding<DIM> , std::vector<int> > make_ball(double rad, CellComplex &comp, Embedding<DIM> &embed) {
    
    
    std::vector<int> rem_verts;
    
    int NV = 0;
    XVec vert_pos = XVec::Zero(DIM*embed.NV);
    
    DVec center = 0.5 * DVec::Ones();
    
    for(int vi = 0; vi < embed.NV; vi++) {
        
        DVec posi = embed.get_vpos(vi);
        
        DVec bvec = embed.get_vdiff(center, posi);
        
        if(bvec.norm() > rad / 2.0) {
            rem_verts.push_back(vi);
        } else {
            vert_pos.segment<DIM>(DIM*NV) = embed.get_vpos(vi);
            NV++;
        }
        
    }
    
    CellComplex red_comp = prune_cell_complex(rem_verts, comp);
    
    Embedding<DIM> red_embed(NV, vert_pos.segment(0, NV*DIM), embed.box_mat, false);
    
    return std::make_tuple(red_comp, red_embed, rem_verts);
    
}


template <int DIM> std::tuple<CellComplex, Embedding<DIM> , std::vector<int> > prune_zero_modes(double rad, CellComplex &comp, Embedding<DIM> &embed) {
    
    
    std::vector<int> rem_verts;
    
    int NV = 0;
    XVec vert_pos = XVec::Zero(DIM*embed.NV);
    
    DVec center = 0.5 * DVec::Ones();
    
    for(int vi = 0; vi < embed.NV; vi++) {
        
        DVec posi = embed.get_vpos(vi);
        
        DVec bvec = embed.get_vdiff(center, posi);
        
        if(bvec.norm() > rad / 2.0) {
            rem_verts.push_back(vi);
        } else {
            vert_pos.segment<DIM>(DIM*NV) = embed.get_vpos(vi);
            NV++;
        }
        
    }
    
    CellComplex red_comp = prune_cell_complex(rem_verts, comp);
    
    Embedding<DIM> red_embed(NV, vert_pos.segment(0, NV*DIM), embed.box_mat, false);
    
    return std::make_tuple(red_comp, red_embed, rem_verts);
    
}



#endif //MECHNET_HPP