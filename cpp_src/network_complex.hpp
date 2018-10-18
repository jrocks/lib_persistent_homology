#ifndef NETWORKCOMPLEX_HPP
#define NETWORKCOMPLEX_HPP
    
#include "eigen_macros.hpp"
#include "cell_complex.hpp"
#include "embedding.hpp"

#include <pybind11/pybind11.h>    
namespace py = pybind11;

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>


    
class Network {
    
public:

    // Number of vertices
    int NV;

    // Number of edges
    int NE;
    // Nodes for each edge
    std::vector<int> edgei;
    std::vector<int> edgej;

    bool has_list;
    std::vector<std::unordered_set<int> > neighbors; 
    
    Network(int NV, int NE, std::vector<int> &edgei, std::vector<int> &edgej) :
         NV(NV), NE(NE), edgei(edgei), edgej(edgej) {
             has_list = false;
         };

    void construct_neighbor_list() {

        if(has_list) {
            return;
        } else {
            has_list = true;
        }

        neighbors.resize(NV);

        for(int ei = 0; ei < NE; ei++) {
            int vi = edgei[ei];
            int vj = edgej[ei];

            neighbors[vi].insert(vj);
            neighbors[vj].insert(vi);

        }
    }
};

    
CellComplex construct_network_complex(Network &net, bool oriented) {
        
    CellComplex comp(1, true, oriented);
       
    // Vertices are guarenteed to have matching cell number and label
    for(int i = 0; i < net.NV; i++) {
        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(i, 0, facets, coeffs);
    }


    for(int i = 0; i < net.NE; i++) {
        std::vector<int> facets;
        facets.push_back(net.edgei[i]);
        facets.push_back(net.edgej[i]);

        std::vector<int> coeffs;
        if(oriented) {
            coeffs.push_back(-1);
            coeffs.push_back(1);
        }
        comp.add_cell(i, 1, facets, coeffs);
 
    }
    
    
    comp.construct_cofacets();
    comp.make_compressed(); 

    return comp;
    
}



template <int DIM> std::vector<double> calc_edge_extensions(RXVec disp, Network &net, Embedding<DIM> &embed, bool is_strain=false) {
        
    std::vector<double> ext;
    
    for(int i = 0; i < net.NE; i++) {

        int vi = net.edgei[i];
        int vj = net.edgej[i];
        
        DVec posi = embed.get_vpos(vi);
        DVec posj = embed.get_vpos(vj);
        
        DVec bvec = posj - posi;
        
        for(int d = 0; d < DIM; d++) {
            if(std::fabs(bvec(d)) > 0.5) {
                bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
            }
        }
        
        bvec = embed.box_mat * bvec;
                
        double l0 = bvec.norm();
        
        bvec.normalize();
        
        DVec ui = disp.segment<DIM>(DIM*vi);
        DVec uj = disp.segment<DIM>(DIM*vj);
        
        DVec du = uj - ui;
        
        double eps = bvec.dot(du);
        
        if(is_strain) {
            eps /= l0;
        }
        
        ext.push_back(eps);
            
        
    }
    
    return ext;
    
}



    
#endif // NETWORKCOMPLEX_HPP