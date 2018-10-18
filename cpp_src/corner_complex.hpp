#ifndef CORNERCOMPLEX_HPP
#define CORNERCOMPLEX_HPP
    
#include "eigen_macros.hpp"
#include "cell_complex.hpp"
#include "embedding.hpp"
#include "network_complex.hpp"

#include <pybind11/pybind11.h>    
namespace py = pybind11;

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>



template <int DIM> std::vector<std::vector<int> > find_corners(Network &net, Embedding<DIM> &embed) {
    
    
    net.construct_neighbor_list();
        
    std::vector<std::vector<int> > corners;
        
    // Iterate over each vertex
    for(int vi = 0; vi < net.NV; vi++) {
                
        // Get all vertices
        std::vector<int> verts;
        // Get all vertex positions relative to center vertex
        std::unordered_map<int, DVec > positions;
        
        verts.push_back(vi);
        positions[vi] = DVec::Zero();
        
        // Center vertex positions
        DVec O = embed.get_vpos(vi);
            
        for(auto vj: net.neighbors[vi]) {
            verts.push_back(vj);
            
            DVec bvec = embed.get_vpos(vj); - O;
            
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
            
            bvec =  embed.box_mat * bvec;
                                      
            positions[vj] = bvec.normalized();
        }
                        
        // Must have at least dim vertices other than the central vertex for a corner to exist
        if((int)verts.size() < DIM + 1) {
            continue;
        }

        // A mask to pick out exactly dim verts
        std::vector<bool> mask(DIM, true);
        mask.resize(verts.size(), false);
        
        // Iterate through each simplex and test if part of convex hull
        std::vector<std::unordered_set<int> > hull_simps;
        do {  

            // Find vertices in simplex
            std::unordered_set<int> simplex;
            for(std::size_t j = 0; j < verts.size(); j++) {
                if(mask[j]) {
                    simplex.insert(verts[j]);
                }
            }
                                    
            // Calculate positions of all verts relative to one vert
            
            XMat cross_mat(DIM, DIM-1);
            
            auto sit = simplex.begin();
            DVec A = positions[*sit];
                        
            sit++;
            for(int m = 0; m < DIM-1; m++, sit++) {
                cross_mat.col(m) = positions[*sit] - A;
            }
                        
            // Use cofactors of cross product matrix to calculate normal vector
            DVec normal;
            for(int m = 0; m < DIM; m++) {
                
                XMat cofactor_mat(DIM-1, DIM-1);
                cofactor_mat.block(0, 0, m, DIM-1) = cross_mat.block(0, 0, m, DIM-1);
                cofactor_mat.block(m, 0, DIM-1-m, DIM-1) = cross_mat.block(m+1, 0, DIM-1-m, DIM-1);
                
                int sign = (m % 2 == 0) ? 1 : -1;
                normal(m) = sign * cofactor_mat.determinant();
            }
                        
            // Check if all points are colinear
            double norm = normal.norm();
            if(norm < 1e-8) {
                continue;
            }
            normal /= norm;
                        
            int above = 0;
            int below = 0;
            for(std::size_t j = 0; j < verts.size(); j++) {
                if(!mask[j]) {
                    
                    DVec v = positions[verts[j]] - A;
                    
                    if(v.dot(normal) > 0) {
                        above++;
                    } else {
                        below++;
                    }
                    
                }
            }
            
            if(above == 0 || below == 0) {
                hull_simps.push_back(simplex);
            }

            
        } while(std::prev_permutation(mask.begin(), mask.end()));
        
        for(auto simplex: hull_simps) {
            
            if(!simplex.count(vi)) {
                
                corners.emplace_back(1, vi);
                corners.back().insert(corners.back().end(), simplex.begin(), simplex.end());
                
                std::sort(corners.back().begin()+1, corners.back().end());
                
            }
            
        }
        
    }
    
    std::sort(corners.begin(), corners.end());
    
    return corners;
    
}


template <int DIM> std::vector<double> calc_corner_strains(std::vector< std::vector<int> > &corners, 
                                                           RXVec disp, Embedding<DIM> &embed) {

    
    std::vector<double> corner_strains;
        
    for(auto corner: corners) {
        
        int vi = corner[0];
        
        DVec O = embed.get_vpos(vi);
        DVec uO = disp.segment<DIM>(DIM*vi);

        DMat X = DMat::Zero();
        DMat Y = DMat::Zero();
        for(int m = 0; m < DIM; m++) {  

            int vj = corner[1+m];
            
            DVec bvec = embed.get_vpos(vj) - O;  
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
        
            bvec = embed.box_mat * bvec;
            
            DVec du = disp.segment<DIM>(DIM*vj) - uO;
            
            X += bvec * bvec.transpose();
            Y += du * bvec.transpose();

        }

        DMat eps = Y * X.inverse();

        eps = 0.5 * (eps + eps.transpose());

        corner_strains.push_back(eps.norm());
         
    }
    
    return corner_strains;
    
    
}



template <int DIM> CellComplex construct_corner_complex(std::vector<std::vector<int> > &corners, Network &net) {
     
    CellComplex comp(DIM, true, false);
          
    // Map of vertices to lists of vertices of all simplices of all corners at that vertex to index of simplex in comp
    std::vector<std::map<std::vector<int> , int> > vertices_to_index(net.NV);
    
    // First process vertices
    for(int i = 0; i < net.NV; i++) {

        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(i, 0, facets, coeffs);
        
        vertices_to_index[i].emplace(std::piecewise_construct, std::forward_as_tuple(1, i), std::forward_as_tuple(i));
        
    }
    
    // Next process edges
    for(int i = 0; i < net.NE; i++) {
        int vi = net.edgei[i];
        int vj = net.edgej[i];
        
        std::vector<int> facets;
        facets.push_back(vi);
        facets.push_back(vj);

        std::vector<int> coeffs;
        
        comp.add_cell(net.NV + i, 1, facets, coeffs);
        
        std::sort(facets.begin(), facets.end());
        vertices_to_index[vi][facets] = net.NV + i;
        vertices_to_index[vj][facets] = net.NV + i;
        
        vertices_to_index[vi].emplace(std::piecewise_construct, std::forward_as_tuple(1, vj), std::forward_as_tuple(vj));
        vertices_to_index[vj].emplace(std::piecewise_construct, std::forward_as_tuple(1, vi), std::forward_as_tuple(vi));
    }
    
    
    std::vector<int> corner_to_cell;
    
    // Iterate through each corner and add all higher-dimensional faces of corner simplices
    for(int d = 1; d <= DIM; d++) {
        
        for(std::size_t i = 0; i < corners.size(); i++) {
            
            int vi = corners[i][0];
            
            // Iterate through every size d+1 subset
            
            // A mask to pick out exactly d+1 verts
            std::vector<bool> mask(d+1, true);
            mask.resize(corners[i].size(), false);
            do {  
                
                // Find vertices in simplex
                std::vector<int> simplex;
                for(std::size_t j = 0; j < corners[i].size(); j++) {
                    if(mask[j]) {
                        simplex.push_back(corners[i][j]);
                    }
                }
                
                // Sorted list of vertices of cell
                std::sort(simplex.begin(), simplex.end());
                
                // If simplex already exists in net complex, then skip                
                if(vertices_to_index[vi].count(simplex)) {
                    continue;
                }
                                
                vertices_to_index[vi][simplex] = comp.ncells;
                
                // Find facets
                std::vector<int> facets;
                std::vector<int> coeffs;
                for(std::size_t j = 0; j < simplex.size(); j++) {
                    std::vector<int> facet(simplex);
                    facet.erase(facet.begin()+j);
                    facets.push_back(vertices_to_index[vi][facet]);
                }
                
                // Label new cells -1 to indicate they don't correspond to anything in original net complex
                // Or if cell corresponds to corner, label it with corner index
                if(d == DIM) {
                    comp.add_cell(i, d, facets, coeffs);
                } else {
                    comp.add_cell(-1, d, facets, coeffs);
                }
                    
                
            } while(std::prev_permutation(mask.begin(), mask.end()));

        }
        
    }
    
    
    comp.construct_cofacets();
    comp.make_compressed(); 

    return comp;
}

        
    

    
#endif // CORNERCOMPLEX_HPP