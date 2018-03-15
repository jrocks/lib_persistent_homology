#ifndef GRAPHCOMPLEX_HPP
#define GRAPHCOMPLEX_HPP
    
#include <Eigen/Core>
#include <Eigen/Geometry>
  
typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;
typedef Eigen::Ref<XVec > RXVec;

#include <pybind11/pybind11.h>    
namespace py = pybind11;

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include "cell_complex.hpp"

    
class Graph {
    
public:

    // Number of vertices
    int NV;

    // Number of edges
    int NE;
    // Nodes for each edge
    std::vector<int> edgei;
    std::vector<int> edgej;

    bool has_list;
    std::vector<std::unordered_set<int> > vertex_neighbors; 
    // std::vector<std::unordered_set<int> > edge_neighbors; 
    
    Graph(int NV, int NE, std::vector<int> &edgei, std::vector<int> &edgej) :
         NV(NV), NE(NE), edgei(edgei), edgej(edgej) {
             has_list = false;
         };

    void construct_neighbor_list() {

        if(has_list) {
            return;
        } else {
            has_list = true;
        }

        vertex_neighbors.resize(NV);

        for(int ei = 0; ei < NE; ei++) {
            int vi = edgei[ei];
            int vj = edgej[ei];

            vertex_neighbors[vi].insert(vj);
            vertex_neighbors[vj].insert(vi);

        }
    }
};

template <int DIM> class Embedding {
 
public:
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    
    // Embedding information
    const int dim;
    // vertex positions
    XVec vert_pos;
    // Box dimensions
    DVec L;
    
    Embedding(RXVec vert_pos, RXVec L) :
        dim(DIM), vert_pos(vert_pos), L(L) {}
    
    
};
    
    
CellComplex construct_graph_complex(Graph &graph, bool oriented) {
        
    CellComplex comp(1, true, oriented);
       
    // Vertices are guarenteed to have matching cell number and label
    for(int i = 0; i < graph.NV; i++) {
        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(i, 0, facets, coeffs);
    }


    for(int i = 0; i < graph.NE; i++) {
        std::vector<int> facets;
        facets.push_back(graph.edgei[i]);
        facets.push_back(graph.edgej[i]);

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



template <int DIM> std::vector<std::vector<int> > find_corners(Graph &graph, Embedding<DIM> &embed) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    
    graph.construct_neighbor_list();
        
    std::vector<std::vector<int> > corners;
        
    // Iterate over each vertex
    for(int vi = 0; vi < graph.NV; vi++) {
                
        // Get all vertices
        std::vector<int> verts;
        // Get all vertex positions relative to center vertex
        std::unordered_map<int, DVec > positions;
        
        verts.push_back(vi);
        positions[vi] = DVec::Zero();
        
        // Center vertex positions
        DVec O = embed.vert_pos.template segment<DIM>(DIM*vi);
            
        for(auto vj: graph.vertex_neighbors[vi]) {
            verts.push_back(vj);
            
            DVec bvec = embed.vert_pos.template segment<DIM>(DIM*vj) - O;
                        
            bvec -= ((bvec.array() / embed.L.array()).round() * embed.L.array()).matrix();
                          
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


template <int DIM> std::vector<double> calc_corner_strains(std::vector< std::vector<int> > &corners, RXVec disp, Embedding<DIM> &embed) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    typedef Eigen::Matrix<double, DIM, DIM> DMat;
    
    std::vector<double> corner_strains;
        
    for(auto corner: corners) {
        
        int vi = corner[0];
        
        DVec O = embed.vert_pos.template segment<DIM>(DIM*vi);
        DVec uO = disp.segment<DIM>(DIM*vi);

        DMat X = DMat::Zero();
        DMat Y = DMat::Zero();
        for(int m = 0; m < DIM; m++) {  

            int vj = corner[1+m];
            DVec bvec = embed.vert_pos.template segment<DIM>(DIM*vj) - O;
            bvec -= ((bvec.array() / embed.L.array()).round() * embed.L.array()).matrix();

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

template <int DIM> std::vector<double> calc_edge_extension(RXVec disp, Graph &graph, Embedding<DIM> &embed, bool is_strain=false) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    
    std::vector<double> ext;
    
    for(int i = 0; i < graph.NE; i++) {

        int vi = graph.edgei[i];
        int vj = graph.edgej[i];
        
        DVec posi = embed.vert_pos.template segment<DIM>(DIM*vi);
        DVec posj = embed.vert_pos.template segment<DIM>(DIM*vj);
        
        DVec bvec = posj - posi;
        bvec -= ((bvec.array() / embed.L.array()).round() * embed.L.array()).matrix();
        
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



// std::tuple<std::vector<double>, std::vector<int> > perform_corner_transform(std::vector<std::vector<int> > &corners, std::vector<double> &corner_strains,
//                                           Graph &graph, bool ascend = true) {
      
    
//     // Return vertex_time which just consists of double vales
//     // and vertex_order which is a discretized version of vertex_time which is a lex sort of lex vals
//     // low to high with each lex val sorted high to low
//     // vertex_order allows for ties so that watershed alg can do its job
//     // give vertex_order to watershed alg, but vertex_time to filtration
    
//     std::vector<double> vertex_time(graph.NV);
//     std::vector<int> vertex_order(graph.NV);
    
//     std::vector<std::vector<double> > lex_val(graph.NV);
//     for(std::size_t ci = 0; ci < corners.size(); ci++) {
//         lex_val[corners[ci][0]].push_back(corner_strains[ci]);
//     }
        
//     for(int i = 0; i < graph.NV; i++) {
//         if(ascend) {
//             // Sort from largest to smallest
//             std::sort(lex_val[i].begin(), lex_val[i].end(), std::greater<int>()); 
//             vertex_time[i] = lex_val[i].back();
//         } else {
//             // Sort from smallest to largest
//             std::sort(lex_val[i].begin(), lex_val[i].end(), std::less<int>());
//             vertex_time[i] = lex_val[i].back();
//         }
        
        
//     }
    
//     // This sorting should be equivalent to sorting on edges
//     auto lex_cmp = [&lex_val, ascend](const int &lhs, const int &rhs) {
        
//         // Ascending filtraiton
//         if(ascend) {
//             // Compare smallest corner (smallest goes first) 
//             if(lex_val[lhs].back() != lex_val[rhs].back()) {
//                 return lex_val[lhs].back() < lex_val[rhs].back();
//             // Compare lexicographic values (smallest goes first)
//             } else if(lex_val[lhs] != lex_val[rhs]) {
//                 return lex_val[lhs] < lex_val[rhs];
//             // Finally, if cells have identical lexicographic orderings, 
//             // then sort by raw cell index
//             } else {
//                 py::print("Breaking tie with indices", lhs, rhs);
//                 return lhs < rhs;
//             }
//         } else {
//             // largest corner (largest goes first) 
//             if(lex_val[lhs].back() != lex_val[rhs].back()) {
//                 return lex_val[lhs].back() > lex_val[rhs].back();
//             // Compare lexicographic values (largest goes first)
//             } else if(lex_val[lhs] != lex_val[rhs]) {
//                 return lex_val[lhs] > lex_val[rhs];
//             // Finally, if cells have identical lexicographic orderings, 
//             // then sort by raw cell index
//             } else {
//                 py::print("Breaking tie with indices", lhs, rhs);
//                 return lhs > rhs;
//             }
//         }
            
//     };
    
//     std::vector<int> cells(graph.NV);
//     std::iota(cells.begin(), cells.end(), 0);
//     std::sort(cells.begin(), cells.end(), lex_cmp);
    
//     int index = 0;
//     vertex_order[cells[0]] = index;
//     for(int i = 1; i < graph.NV; i++) {
//         if(lex_val[cells[i]] > lex_val[cells[i-1]]) {
//             index++;
//         }
        
//         vertex_order[cells[i]] = index;
//     }
      
//     return std::forward_as_tuple(vertex_time, vertex_order);
    
    
// }                                         
        
    
template <int DIM> CellComplex construct_corner_complex(std::vector<std::vector<int> > &corners, Graph &graph) {
     
    CellComplex comp(DIM, true, false);
          
    // Map of vertices to lists of vertices of all simplices of all corners at that vertex to index of simplex in comp
    std::vector<std::map<std::vector<int> , int> > vertices_to_index(graph.NV);
    
    // First process vertices
    for(int i = 0; i < graph.NV; i++) {

        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(i, 0, facets, coeffs);
        
        vertices_to_index[i].emplace(std::piecewise_construct, std::forward_as_tuple(1, i), std::forward_as_tuple(i));
        
    }
    
    // Next process edges
    for(int i = 0; i < graph.NE; i++) {
        int vi = graph.edgei[i];
        int vj = graph.edgej[i];
        
        std::vector<int> facets;
        facets.push_back(vi);
        facets.push_back(vj);

        std::vector<int> coeffs;
        
        comp.add_cell(graph.NV + i, 1, facets, coeffs);
        
        std::sort(facets.begin(), facets.end());
        vertices_to_index[vi][facets] = graph.NV + i;
        vertices_to_index[vj][facets] = graph.NV + i;
        
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
                
                // If simplex already exists in graph complex, then skip                
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
                
                // Label new cells -1 to indicate they don't correspond to anything in original graph complex
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


    
#endif // GRAPHCOMPLEX_HPP