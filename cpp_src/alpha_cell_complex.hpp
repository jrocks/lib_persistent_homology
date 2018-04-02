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
    
    CellComplex comp(DIM, true, false);
    
    for(int i = 0; i < NV; i++) {
        
        DVec pos = vert_pos.segment<DIM>(DIM*i);
        std::vector<double> coords(&vert_pos(DIM*i), &vert_pos(DIM*i)+DIM);
        WPoint w(Point(coords.begin(), coords.end()), weights[i]);
        
        auto vertex = t.insert(w);
        vertex->data() = i;
        
        
        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(i, 0, facets, coeffs);
        
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
                                
                simplex_to_index[simplex] = comp.ncells;
                
                // Find facets
                std::vector<int> facets;
                std::vector<int> coeffs;
                for(std::size_t j = 0; j < simplex.size(); j++) {
                    std::vector<int> facet(simplex);
                    facet.erase(facet.begin()+j);
                    facets.push_back(simplex_to_index[facet]);
                }
                comp.add_cell(comp.ncells, d, facets, coeffs);
                    
                
            } while(std::prev_permutation(mask.begin(), mask.end()));

        }
        
    }
    
    
    comp.construct_cofacets();
    comp.make_compressed(); 

    return comp;
    
}


template <int DIM> double calc_radius_squared(std::vector<int> &vertices, RXVec vert_pos, std::vector<double> &weights) {
    
    typedef Eigen::Matrix<double, DIM, 1> DVec;
    typedef Eigen::Matrix<double, DIM, DIM> DMat;
    
    if(vertices.size() == 1) {
        return 0.0;
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
    
        return x(DIM) + x.segment<DIM>(0).squaredNorm();
        
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
    
        return x(DIM) + x.segment<DIM>(0).squaredNorm();
        
        
    }
    
    
    
    
}

template <int DIM> std::vector<double> calc_alpha_vals(RXVec vert_pos, std::vector<double> &weights, CellComplex &comp) {
    
    std::vector<double> alpha_vals(comp.ncells);
    
    double min_alpha = 1e10;
    
    for(int i = 0; i < comp.ncells; i++) {
        std::unordered_set<int> star = get_star(i, true, comp, 0);
        std::vector<int> vertices;
        vertices.reserve(star.size());
        for(auto j: star) {
            vertices.push_back(comp.get_label(j));
        }
        
        if(comp.get_dim(i) != 0) {
            alpha_vals[i] = calc_radius_squared<DIM>(vertices, vert_pos, weights);
            
            if(alpha_vals[i] < min_alpha) {
                min_alpha = alpha_vals[i];
            }
        }
    }
    
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) == 0) {
            alpha_vals[i] = min_alpha;
        }
    }
    
    return alpha_vals;
    
}


    
    
#endif //ALPHACOMPLEX_HPP