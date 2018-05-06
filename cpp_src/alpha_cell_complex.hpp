#ifndef ALPHACOMPLEX_HPP
#define ALPHACOMPLEX_HPP
    
    
#include <Eigen/Core>
#include <Eigen/Dense>

typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;
typedef Eigen::Ref<XVec > RXVec;
    
#include <CGAL/Epick_d.h>
// #include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Regular_triangulation.h>

#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <utility>
    
#include "cell_complex.hpp"
    
#include <pybind11/pybind11.h>
namespace py = pybind11;
    

template <int DIM> double calc_power_distance(double alpha, Eigen::Matrix<double, DIM, 1> a, Eigen::Matrix<double, DIM, 1> pos, double weight) {
    return (a - pos).squaredNorm() - weight - alpha;
}


template <int DIM> std::tuple<double, Eigen::Matrix<double, DIM, 1> > calc_circumsphere(std::vector<int> &vertices, RXVec vert_pos, std::vector<double> &weights) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    typedef Eigen::Matrix<double, DIM, DIM> DMat;
    
    if(vertices.size() == 1) {
        return std::make_tuple(0.0, DVec::Zero());
    } else if(vertices.size() == DIM + 1) {
        
        XMat A = XMat::Zero(DIM+1, DIM+1);
        
        for(std::size_t i = 0; i < DIM+1; i++) {
            A.block<1,DIM>(i, 0) = 2.0*vert_pos.segment<DIM>(DIM*vertices[i]);
        }
        
        A.block<DIM+1,1>(0, DIM) = XVec::Ones(DIM+1);        
        
        XVec b = XVec::Zero(DIM+1);
        for(std::size_t i = 0; i < DIM+1; i++) {
            int vi = vertices[i];
            b(i) = vert_pos.segment<DIM>(DIM*vi).squaredNorm() - weights[vi];
        }
                
        XVec x = A.partialPivLu().solve(b);
    
        return std::make_tuple(x(DIM) + x.segment<DIM>(0).squaredNorm(), x.segment<DIM>(0));
        
    } else {
        
        XMat A = XMat::Zero(DIM+vertices.size(), DIM+vertices.size());
        
        for(std::size_t i = 0; i < vertices.size(); i++) {
            A.block<1,DIM>(i, 0) = 2.0*vert_pos.segment<DIM>(DIM*vertices[i]);
        }
        
        A.block(0, DIM, vertices.size(), 1) = XVec::Ones(vertices.size());
        
        DVec v0 = vert_pos.segment<DIM>(DIM*vertices[0]);
        for(std::size_t i = 1; i < vertices.size(); i++) {
            A.block<DIM,1>(vertices.size(), DIM+i) = v0 - vert_pos.segment<DIM>(DIM*vertices[i]);
        }
        
        A.block<DIM, DIM>(vertices.size(), 0) = DMat::Identity();
        
        
        XVec b = XVec::Zero(DIM+vertices.size());
        for(std::size_t i = 0; i < vertices.size(); i++) {
            int vi = vertices[i];
            b(i) = vert_pos.segment<DIM>(DIM*vi).squaredNorm() - weights[vi];
        }
        
        b.segment<DIM>(vertices.size()) = v0;
        
        XVec x = A.partialPivLu().solve(b);
    
        return std::make_tuple(x(DIM) + x.segment<DIM>(0).squaredNorm(), x.segment<DIM>(0));
        
        
    }
    
    
    
    
}

template <int DIM> CellComplex construct_alpha_complex(int NV, RXVec vert_pos, std::vector<double> &weights) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    
    // d-dimensional Kernel used to define Euclidean space (R^d)
    typedef CGAL::Epick_d< CGAL::Dimension_tag<DIM> > Kernel;
    
    // Triangulation data structure 
    // Explicitly defined in order to stick an integer label into each Triangulation_vertex
    typedef CGAL::Triangulation_data_structure<typename Kernel::Dimension,
                                            CGAL::Triangulation_vertex<CGAL::Regular_triangulation_traits_adapter<Kernel>, int>,
                                            CGAL::Triangulation_full_cell<CGAL::Regular_triangulation_traits_adapter<Kernel> > >
                                                Triangulation_data_structure;
    // Regular Delaunay triangulation
    typedef CGAL::Regular_triangulation<Kernel, Triangulation_data_structure> Regular_triangulation;


    // d-dimensional point
    typedef typename Kernel::Point_d Point;
    // d-dimensional weighted point
    typedef typename Kernel::Weighted_point_d WPoint;
    
    // Triangulation object
    Regular_triangulation t(DIM);
    
    // Map of lists of vertices of all simplices to index of simplex in comp
    std::map<std::vector<int> , int> simplex_to_index;
    
    CellComplex del_comp(DIM, true, true);
    
    // Add all vertices to cell complex
    for(int i = 0; i < NV; i++) {
        
        DVec pos = vert_pos.segment<DIM>(DIM*i);
        std::vector<double> coords(&vert_pos(DIM*i), &vert_pos(DIM*i)+DIM);
        WPoint w(Point(coords.begin(), coords.end()), weights[i]);
        
        auto vertex = t.insert(w);
        vertex->data() = i;
        
        
        std::vector<int> facets;
        std::vector<int> coeffs;
        del_comp.add_cell(i, 0, facets, coeffs);
        
        simplex_to_index.emplace(std::piecewise_construct, std::forward_as_tuple(1, i), std::forward_as_tuple(i));
    }
    
    
    // py::print("Regular triangulation successfully computed: " , t.number_of_vertices(), " vertices, ",
    // t.number_of_finite_full_cells()," finite cells.");

//     for(auto it = t.finite_full_cells_begin(); it != t.finite_full_cells_end(); it++) {
//         py::print("Cell:");
                
//         for(auto vit = it->vertices_begin(); vit != it->vertices_end(); vit++) {
//             // Need to dereference first, since vit is a pointer to a vertex handle
//             py::print((*vit)->data());
                        
//         }
        
//     }
    
    // Iterate through each corner and add all higher-dimensional faces of corner simplices
    for(int d = 1; d <= DIM; d++) {
        
        int label = 0;
        
        for(auto it = t.finite_full_cells_begin(); it != t.finite_full_cells_end(); it++) {
            
            std::vector<int> vertices;
            for(auto vit = it->vertices_begin(); vit != it->vertices_end(); vit++) {
                // Need to dereference first, since vit is a pointer to a vertex handle
                vertices.push_back((*vit)->data());

            }
                        
            // Iterate through every size d+1 subset
            
            // A mask to pick out exactly d+1 verts
            std::vector<bool> mask(d+1, true);
            mask.resize(vertices.size(), false);
            do {  
                
                // Find vertices in simplex
                std::vector<int> simplex;
                for(std::size_t j = 0; j < vertices.size(); j++) {
                    if(mask[j]) {
                        simplex.push_back(vertices[j]);
                    }
                }
                
                // Sorted list of vertices of cell
                std::sort(simplex.begin(), simplex.end());
                
                // If simplex already exists in graph complex, then skip                
                if(simplex_to_index.count(simplex)) {
                    continue;
                }
                                
                simplex_to_index[simplex] = del_comp.ncells;
                
                // Find facets
                std::vector<int> facets;
                std::vector<int> coeffs;
                for(std::size_t j = 0; j < simplex.size(); j++) {
                    std::vector<int> facet(simplex);
                    facet.erase(facet.begin()+j);
                    facets.push_back(simplex_to_index[facet]);
                    coeffs.push_back(-2*(j%2)+1);
                }
                del_comp.add_cell(label, d, facets, coeffs);
                label++;
                
            } while(std::prev_permutation(mask.begin(), mask.end()));

        }
        
    }
    
    del_comp.construct_cofacets();
    del_comp.make_compressed(); 
    
    return del_comp;
    
}



template <int DIM> std::vector<double> calc_alpha_vals(RXVec vert_pos, std::vector<double> &weights, CellComplex &comp) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    
    std::vector<double> alpha_vals(comp.ncells);
    
    double min_alpha = 1e10;
    
    for(int c = comp.ncells-1; c >= 0; c--) {
        
        // Skip vertices
        if(comp.get_dim(c) == 0) {
            continue;
        }
        
        std::unordered_set<int> verts = get_star(c, true, comp, 0);
        
        double alpha;
        DVec a;
        std::vector<int> tmp(verts.begin(), verts.end());
        std::tie(alpha, a) = calc_circumsphere<DIM>(tmp, vert_pos, weights);
        
        alpha_vals[c] = alpha;
        
        if(alpha < min_alpha) {
            min_alpha = alpha;
        }
        
        // Skip highest dimension triangles
        if(comp.get_dim(c) == DIM) {
            continue;
        }
        
        // For a Delaunay triangulation, only simplices of dimension DIM must have empty circumspheres
        // Some simplices between dimension 0 and DIM don't appear in the alpha shape filtration at the value of their particular alpha
        // This is because for a particular value of alpha, each simplex follows one of two rules:
        // 1. has an empty circumsphere or
        // 2. is the face of another simplex that is part of the alpha complex (see 1)
        // We must detect simplices that don't obey #1 and change their value of alpha to that of their youngest coface
        
        // Find cofaces of dimension DIM
        auto cofaces = get_star(c, false, comp, DIM);
        
        bool conflict = false;
        // For each coface get the list of its vertices and check if any fall inside the circumsphere of c
        for(auto cf: cofaces) {
            
            std::unordered_set<int> coface_verts = get_star(cf, true, comp, 0);
            
            for(auto vi: coface_verts) {

                // If inside the circumsphere, mark as conflict and continue
                DVec x = vert_pos.segment<DIM>(DIM*vi);
                if(!verts.count(vi) && calc_power_distance(alpha, a, x, weights[vi]) < 0.0) {
                    conflict = true;
                    break;
                }
            }
            
            if(conflict) {
                break;
            }
            
        }
        
        
        if(!conflict) {
            continue;
        }
        
        // If this simplex poses a conflict, then set it's value of alpha to it's youngest cofacet
        auto cofacets = comp.get_cofacets(c);
        double conflict_alpha = 1e10;
        for(auto cf: cofacets) { 
            if(alpha_vals[cf] < conflict_alpha) {
                conflict_alpha = alpha_vals[cf];
            }
        }
        
        alpha_vals[c] = conflict_alpha;
        
        
    }
    
    
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) == 0) {
            alpha_vals[i] = min_alpha;
        }
    }
    
    return alpha_vals;
    
}


    
    
#endif //ALPHACOMPLEX_HPP