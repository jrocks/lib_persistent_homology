#ifndef GRAPHCOMPLEX
#define GRAPHCOMPLEX
    
#include <Eigen/Core>
#include <Eigen/Geometry>
  
typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;

typedef Eigen::Matrix<double, 2, 1> D2Vec;
typedef Eigen::Matrix<double, 3, 1> D3Vec;
    
typedef Eigen::Ref<XVec > RXVec;


#include <pybind11/pybind11.h>    
namespace py = pybind11;

#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include "cell_complex.hpp"

    
CellComplex construct_graph_complex(int NV, int NE, std::vector<int> &edgei, std::vector<int> &edgej, bool oriented) {
    
    int dim = 1;
    
    CellComplex comp(dim, true, oriented);
       
    
    // Vertices are guarenteed to have matching cell number and label
    for(int i = 0; i < NV; i++) {
        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(i, 0, facets, coeffs);
    }


    for(int i = 0; i < NE; i++) {
        std::vector<int> facets;
        facets.push_back(edgei[i]);
        facets.push_back(edgej[i]);

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

std::vector<std::vector<int> > find_corners(CellComplex &comp, int dim, RXVec vert_pos, RXVec L) {
    
    std::vector<std::vector<int> > corners;
        
    // Iterate over each vertex
    for(int i = 0; i < comp.ncells; i++) {
        
        if(comp.get_dim(i) != 0) {
            continue;
        }
                
        // Get all vertices
        std::vector<int> verts;
        // Get all vertex positions relative to center vertex
        std::unordered_map<int, XVec > positions;
        
        verts.push_back(i);
        positions[i] = XVec::Zero(dim);
        
        // Center vertex positions
        XVec O = vert_pos.segment(dim*i, dim);
            
        auto rangea = comp.get_cofacet_range(i);
        for(auto it = rangea.first; it != rangea.second; it++) {
            auto rangeb = comp.get_facet_range(*it);
            int vi = *(rangeb.first);
            int vj = *(rangeb.first+1);
            
            int j = (vi == i) ? vj : vi;
            
            verts.push_back(j);
            

            XVec bvec = vert_pos.segment(dim*j, dim) - O;
                        
            bvec -= ((bvec.array() / L.array()).round() * L.array()).matrix();
                          
            positions[j] = bvec.normalized();
          
            
        }
                        
        // Must have at least dim vertices other than the central vertex for a corner to exist
        if((int)verts.size() < dim + 1) {
            continue;
        }

        // A mask to pick out exactly dim verts
        std::vector<bool> mask(dim, true);
        mask.resize(verts.size(), false);
        
        // Iterate through each simplex and test if part of convex hull
        std::vector<std::unordered_set<int> > hull_simps;
        do {  

            // Find vertices in simplex
            std::unordered_set<int> simplex;
            for(unsigned int j = 0; j < verts.size(); j++) {
                if(mask[j]) {
                    simplex.insert(verts[j]);
                }
            }
                                    
            // Calculate positions of all verts relative to one vert
            XMat cross_mat(dim, dim-1);
            
            auto sit = simplex.begin();
            XVec A = positions[*sit];
                        
            sit++;
            for(int m = 0; m < dim-1; m++, sit++) {
                cross_mat.col(m) = positions[*sit] - A;
            }
                        
            // Use cofactors of cross product matrix to calculate normal vector
            XVec normal(dim);
            for(int m = 0; m < dim; m++) {
                
                XMat cofactor_mat(dim-1, dim-1);
                cofactor_mat.block(0, 0, m, dim-1) = cross_mat.block(0, 0, m, dim-1);
                cofactor_mat.block(m, 0, dim-1-m, dim-1) = cross_mat.block(m+1, 0, dim-1-m, dim-1);
                
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
            for(unsigned int j = 0; j < verts.size(); j++) {
                if(!mask[j]) {
                    
                    XVec v = positions[verts[j]] - A;
                    
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
            
            if(!simplex.count(i)) {
                
                corners.emplace_back(1, i);
                corners.back().insert(corners.back().end(), simplex.begin(), simplex.end());
                
            }
            
        }
        
    }
    
    return corners;
    
}

std::vector<double> calc_corner_strains(std::vector< std::vector<int> > &corners, RXVec disp, CellComplex &comp, int dim, int NV, RXVec vert_pos, RXVec L) {
    
    std::vector<double> corner_strains;
    
    for(auto corner: corners) {
        
        int vi = corner[0];
        
        XVec O = vert_pos.segment(dim*vi, dim);
        XVec uO = disp.segment(dim*vi, dim);

        XMat X = XMat::Zero(dim, dim);
        XMat Y = XMat::Zero(dim, dim);

        for(int m = 0; m < dim; m++) {  

            int vj = corner[1+m];
            XVec bvec = vert_pos.segment(dim*vj, dim) - O;
            bvec -= ((bvec.array() / L.array()).round() * L.array()).matrix();

            XVec du = disp.segment(dim*vj, dim) - uO;
            
            X += bvec * bvec.transpose();
            Y += du * bvec.transpose();

        }

        XMat eps = Y * X.inverse();

        eps = 0.5 * (eps + eps.transpose());

        corner_strains.push_back(eps.norm());
         
    }
    
    return corner_strains;
    
    
}

std::vector<double> calc_edge_extension(RXVec disp, CellComplex &comp, int dim, RXVec vert_pos, RXVec L) {
    
    
    std::vector<double> ext;
    
    for(int i = 0; i < comp.ncells; i++) {
        
        if(comp.get_dim(i) != 1) {
            continue;
        }

        auto range = comp.get_facet_range(i);
        int vi = *(range.first);
        int vj = *(range.first+1);
        
        XVec posi = vert_pos.segment(dim*vi, dim);
        XVec posj = vert_pos.segment(dim*vj, dim);
        
        XVec bvec = posj - posi;
        bvec -= ((bvec.array() / L.array()).round() * L.array()).matrix();
        bvec.normalize();
        
        XVec ui = disp.segment(dim*vi, dim);
        XVec uj = disp.segment(dim*vj, dim);
        
        XVec du = uj - ui;
        
        ext.push_back(bvec.dot(du));
            
        
    }
    
    return ext;
    
}


std::vector<double> corner_to_edge_time(std::vector<std::vector<int> > corners, std::vector<double> &corner_time, CellComplex &comp, bool low_to_high) {
    
    
    // Sort each corner in order of corner insertion time
    auto cmp = [&corner_time, low_to_high](const int &lhs, const int &rhs) {
        if(low_to_high) {                    
            return corner_time[lhs] < corner_time[rhs];
        } else {
            return corner_time[lhs] > corner_time[rhs];
        }
    };
    
    std::vector<int> order(corners.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), cmp);
        
    
    // Create map of vertices to edges
    std::map<std::pair<int, int>, int> verts_to_edge;
    int index = 0;
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) != 1) {
            continue;
        }
        
        auto range = comp.get_facet_range(i);
        int vi = *(range.first);
        int vj = *(range.first+1);
        
        verts_to_edge.emplace(std::piecewise_construct, std::forward_as_tuple(vi, vj), 
                    std::forward_as_tuple(index));
        
        verts_to_edge.emplace(std::piecewise_construct, std::forward_as_tuple(vj, vi), 
                    std::forward_as_tuple(index));
    }
    
    
    
    // Calculate insertion time of each edge
    std::vector<double> edge_time;
    for(auto corner: order) {
        
    }
   
    
    
    
}
    
#endif // GRAPHCOMPLEX