#ifndef GRAPHCOMPLEX_HPP
#define GRAPHCOMPLEX_HPP
    
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
        
        // Embedding information
        int dim;
        // vertex positions
        XVec vert_pos;
        // Box dimensions
        XVec L;
    
        bool has_list;
        std::vector<std::unordered_set<int> > vertex_neighbors; 
        // std::vector<std::unordered_set<int> > edge_neighbors; 
    
        Graph(int NV, int NE, std::vector<int> &edgei, std::vector<int> &edgej) :
             NV(NV), NE(NE), edgei(edgei), edgej(edgej) {
                 has_list = false;
             };
    
        void set_embedding(int dim, RXVec vert_pos, RXVec L) {
            this->dim = dim;
            this->vert_pos = vert_pos;
            this->L = L;
        }
    
        void construct_neighbor_list() {
            
            if(has_list) {
                return;
            } else {
                has_list = true;
            }
            
            vertex_neighbors.resize(NV);
            // edge_neighbors.resize(NV);
            
            for(int ei = 0; ei < NE; ei++) {
                int vi = edgei[ei];
                int vj = edgej[ei];
                
                vertex_neighbors[vi].insert(vj);
                vertex_neighbors[vj].insert(vi);
                
                // edge_neighbors[vi].insert(ei);
                // edge_neighbors[vj].insert(ei);
                
            }
        }
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



std::vector<std::vector<int> > find_corners(Graph &graph) {
    
    graph.construct_neighbor_list();
    
    int dim = graph.dim;
    
    std::vector<std::vector<int> > corners;
        
    // Iterate over each vertex
    for(int vi = 0; vi < graph.NV; vi++) {
                
        // Get all vertices
        std::vector<int> verts;
        // Get all vertex positions relative to center vertex
        std::unordered_map<int, XVec > positions;
        
        verts.push_back(vi);
        positions[vi] = XVec::Zero(dim);
        
        // Center vertex positions
        XVec O = graph.vert_pos.segment(dim*vi, dim);
            
        for(auto vj: graph.vertex_neighbors[vi]) {
            verts.push_back(vj);
            
            XVec bvec = graph.vert_pos.segment(dim*vj, dim) - O;
                        
            bvec -= ((bvec.array() / graph.L.array()).round() * graph.L.array()).matrix();
                          
            positions[vj] = bvec.normalized();
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
            for(std::size_t j = 0; j < verts.size(); j++) {
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
            for(std::size_t j = 0; j < verts.size(); j++) {
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


std::vector<double> calc_corner_strains(std::vector< std::vector<int> > &corners, RXVec disp, Graph &graph) {
    
    std::vector<double> corner_strains;
    
    int dim = graph.dim;
    
    for(auto corner: corners) {
        
        int vi = corner[0];
        
        XVec O = graph.vert_pos.segment(dim*vi, dim);
        XVec uO = disp.segment(dim*vi, dim);

        XMat X = XMat::Zero(dim, dim);
        XMat Y = XMat::Zero(dim, dim);

        for(int m = 0; m < dim; m++) {  

            int vj = corner[1+m];
            XVec bvec = graph.vert_pos.segment(dim*vj, dim) - O;
            bvec -= ((bvec.array() / graph.L.array()).round() * graph.L.array()).matrix();

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

std::vector<double> calc_edge_extension(RXVec disp, Graph &graph) {
    
    int dim = graph.dim;
    
    std::vector<double> ext;
    
    for(int i = 0; i < graph.NE; i++) {

        int vi = graph.edgei[i];
        int vj = graph.edgej[i];
        
        XVec posi = graph.vert_pos.segment(dim*vi, dim);
        XVec posj = graph.vert_pos.segment(dim*vj, dim);
        
        XVec bvec = posj - posi;
        bvec -= ((bvec.array() / graph.L.array()).round() * graph.L.array()).matrix();
        bvec.normalize();
        
        XVec ui = disp.segment(dim*vi, dim);
        XVec uj = disp.segment(dim*vj, dim);
        
        XVec du = uj - ui;
        
        ext.push_back(bvec.dot(du));
            
        
    }
    
    return ext;
    
}



std::tuple<std::vector<double>, std::vector<int> > perform_corner_transform(std::vector<std::vector<int> > &corners, std::vector<double> &corner_strains,
                                          Graph &graph, bool ascend = true) {
      
    
    // Return vertex_time which just consists of double vales
    // and vertex_order which is a discretized version of vertex_time which is a lex sort of lex vals
    // low to high with each lex val sorted high to low
    // vertex_order allows for ties so that watershed alg can do its job
    // give vertex_order to watershed alg, but vertex_time to filtration
    
    std::vector<double> vertex_time(graph.NV);
    std::vector<int> vertex_order(graph.NV);
    
    std::vector<std::vector<double> > lex_val(graph.NV);
    for(std::size_t ci = 0; ci < corners.size(); ci++) {
        lex_val[corners[ci][0]].push_back(corner_strains[ci]);
    }
        
    for(int i = 0; i < graph.NV; i++) {
        if(ascend) {
            // Sort from largest to smallest
            std::sort(lex_val[i].begin(), lex_val[i].end(), std::greater<int>()); 
            vertex_time[i] = lex_val[i].back();
        } else {
            // Sort from smallest to largest
            std::sort(lex_val[i].begin(), lex_val[i].end(), std::less<int>());
            vertex_time[i] = lex_val[i].back();
        }
        
        
    }
    
    // This sorting should be equivalent to sorting on edges
    auto lex_cmp = [&lex_val, ascend](const int &lhs, const int &rhs) {
        
        // Ascending filtraiton
        if(ascend) {
            // Compare smallest corner (smallest goes first) 
            if(lex_val[lhs].back() != lex_val[rhs].back()) {
                return lex_val[lhs].back() < lex_val[rhs].back();
            // Compare lexicographic values (smallest goes first)
            } else if(lex_val[lhs] != lex_val[rhs]) {
                return lex_val[lhs] < lex_val[rhs];
            // Finally, if cells have identical lexicographic orderings, 
            // then sort by raw cell index
            } else {
                py::print("Breaking tie with indices", lhs, rhs);
                return lhs < rhs;
            }
        } else {
            // largest corner (largest goes first) 
            if(lex_val[lhs].back() != lex_val[rhs].back()) {
                return lex_val[lhs].back() > lex_val[rhs].back();
            // Compare lexicographic values (largest goes first)
            } else if(lex_val[lhs] != lex_val[rhs]) {
                return lex_val[lhs] > lex_val[rhs];
            // Finally, if cells have identical lexicographic orderings, 
            // then sort by raw cell index
            } else {
                py::print("Breaking tie with indices", lhs, rhs);
                return lhs > rhs;
            }
        }
            
    };
    
    std::vector<int> cells(graph.NV);
    std::iota(cells.begin(), cells.end(), 0);
    std::sort(cells.begin(), cells.end(), lex_cmp);
    
    int index = 0;
    vertex_order[cells[0]] = index;
    for(int i = 1; i < graph.NV; i++) {
        if(lex_val[cells[i]] > lex_val[cells[i-1]]) {
            index++;
        }
        
        vertex_order[cells[i]] = index;
    }
      
    return std::forward_as_tuple(vertex_time, vertex_order);
    
    
}                                         
        

// 1. find corners and corner strains
// 2. construct corner and graph complexes
// 3. watershed on corner complex
// 4. construct ascending and descending filtrations on corner complex
// 5. reduce both to graph complex
                                          
CellComplex construct_corner_complex(std::vector<std::vector<int> > &corners, Graph &graph) {
     
    CellComplex comp(graph.dim, true, false);
          
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
    for(int d = 1; d <= graph.dim; d++) {
        
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
                if(d == graph.dim) {
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