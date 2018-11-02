#ifndef ALPHACOMPLEX_HPP
#define ALPHACOMPLEX_HPP
    

// d-dimensional triangulations
#include <CGAL/Epick_d.h>
// #include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Regular_triangulation.h>
    
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>
    
#include "eigen_macros.hpp"  
#include "embedding.hpp"
#include "cell_complex.hpp"
#include "filtration.hpp"
    
#include <pybind11/pybind11.h>
namespace py = pybind11;
    


template <int DIM> CellComplex construct_alpha_complex(int NV, Embedding<DIM> &embed, 
                                                       std::vector<double> &weights, bool oriented=false) {
        
    
    // d-dimensional Kernel used to define Euclidean space (R^d)
    typedef CGAL::Epick_d< CGAL::Dimension_tag<DIM> > Kernel;

    // Triangulation data structure 
    // Template info is explicitly defined in order to stick an integer label into each Triangulation_vertex
    typedef CGAL::Triangulation_data_structure<
                typename Kernel::Dimension,
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
    Regular_triangulation tri(DIM);
    
    // Map of lists of vertices of all simplices to index of simplex in alpha complex
    std::map<std::vector<int> , int> simplex_to_index;

    CellComplex alpha_comp(DIM, true, oriented);

    // Add all vertices to cell complex
    for(int i = 0; i < NV; i++) {

        DVec pos = embed.get_pos(i);
        std::vector<double> coords(pos.data(), pos.data()+DIM);
        WPoint w(Point(coords.begin(), coords.end()), weights[i]);

        auto vertex = tri.insert(w);
        vertex->data() = i;

        std::vector<int> facets;
        std::vector<int> coeffs;
        alpha_comp.add_cell(i, 0, facets, coeffs);

        simplex_to_index.emplace(std::piecewise_construct, std::forward_as_tuple(1, i), std::forward_as_tuple(i));
    }
    
    
    if(embed.periodic) {
        
        
        // First create list of images represented by the offset from the main image
        std::vector<double> offset{-1, 0, 1};
        
        std::vector<std::vector<double> > image_list(1, std::vector<double>());
        
        // Calculate cartesion product of offset
        for(int d = 0; d < DIM; d++) {
            std::vector<std::vector<double> > new_image_list;
            
            for(auto image: image_list) {
                for(auto off: offset) {
                    
                    auto copy = image;
                    copy.push_back(off);
                    
                    new_image_list.push_back(copy);
                    
                }
            }
            
            image_list = new_image_list;
            
        }
        
        // Delete the main image as it already exists
        std::vector<double> zero_image(DIM, 0);
        image_list.erase(std::remove(image_list.begin(), image_list.end(), zero_image), image_list.end());   
        
        // Add each image to the triangulation
        int vi = NV;
        for(auto image: image_list) {

            DVec offset = XVecMap(image.data(), DIM);
            
            for(int i = 0; i < NV; i++) {

                DVec pos = embed.box_mat * (embed.get_vpos(i) + offset);
                
                std::vector<double> coords(pos.data(), pos.data()+DIM);
                WPoint w(Point(coords.begin(), coords.end()), weights[i]);

                auto vertex = tri.insert(w);
                vertex->data() = vi;

                vi++;
            }

        }
        
        
    }



    // py::print("Regular triangulation successfully computed: " , tri.number_of_vertices(), " vertices, ",
    // tri.number_of_finite_full_cells()," finite cells.");

//     for(auto it = tri.finite_full_cells_begin(); it != tri.finite_full_cells_end(); it++) {
//         py::print("Cell:");

//         for(auto vit = it->vertices_begin(); vit != it->vertices_end(); vit++) {
//             // Need to dereference first, since vit is a pointer to a vertex handle
//             py::print((*vit)->data());

//         }

//     }

    // Iterate through each corner and add all higher-dimensional faces of corner simplices
    for(int d = 1; d <= DIM; d++) {

        int label = 0;

        for(auto it = tri.finite_full_cells_begin(); it != tri.finite_full_cells_end(); it++) {

            bool has_central_image = false;
            std::vector<int> vertices;
            for(auto vit = it->vertices_begin(); vit != it->vertices_end(); vit++) {
                // Need to dereference first, since vit is a pointer to a vertex handle
                // Mode by NV to take care of periodic BCs if turned on
                int vi = (*vit)->data();
                vertices.push_back(vi % NV);
                
                if(vi < NV) {
                    has_central_image = true;
                }

            }
            
            if(!has_central_image) {
                continue;
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

                simplex_to_index[simplex] = alpha_comp.ncells;

                // Find facets
                std::vector<int> facets;
                std::vector<int> coeffs;
                for(std::size_t j = 0; j < simplex.size(); j++) {
                    std::vector<int> facet(simplex);
                    facet.erase(facet.begin()+j);
                    facets.push_back(simplex_to_index[facet]);
                    coeffs.push_back(-2*(j%2)+1);
                }
                alpha_comp.add_cell(label, d, facets, coeffs);
                label++;

            } while(std::prev_permutation(mask.begin(), mask.end()));

        }

    }
        
    // }
    
    alpha_comp.construct_cofacets();
    alpha_comp.make_compressed(); 
    
    return alpha_comp;
    
}


template <int DIM> double calc_power_distance(double alpha, DVec a, DVec pos, double weight=0.0) {
    return (a - pos).squaredNorm() - weight - alpha;
}


template <int DIM> std::tuple<double, DVec > calc_circumsphere(std::vector<int> &vertices, 
                                                               Embedding<DIM> &embed, std::vector<double> &weights) {
    
    
    
    // Find vertex positions relative to first vertex
    XMat X = XMat::Zero(DIM, vertices.size());
    for(std::size_t i = 0; i < vertices.size(); i++) {
        X.block<DIM, 1>(0, i) = embed.get_vpos(vertices[i]);
        
        if(embed.periodic && i > 0) {
            DVec bvec = X.block<DIM, 1>(0, i) - X.block<DIM, 1>(0, 0);
            
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
            
            X.block<DIM, 1>(0, i) = X.block<DIM, 1>(0, 0) + bvec;
            
        }
        
    }
    
    // Scale vertex positions to box size and shape
    for(std::size_t i = 0; i < vertices.size(); i++) {
        X.block<DIM, 1>(0, i) = embed.box_mat * X.block<DIM, 1>(0, i);
    }
    
    if(vertices.size() == 1) {
        return std::make_tuple(0.0, DVec::Zero());
    } else if(vertices.size() == DIM + 1) {
        
        XMat A = XMat::Zero(DIM+1, DIM+1);
        
        A.block<DIM+1, DIM>(0, 0) = 2.0 * X.block<DIM, DIM+1>(0, 0).transpose();
        
        A.block<DIM+1,1>(0, DIM) = XVec::Ones(DIM+1);        
        
        XVec b = XVec::Zero(DIM+1);
        for(std::size_t i = 0; i < DIM+1; i++) {
            int vi = vertices[i];
            
            b(i) = X.block<DIM, 1>(0, i).squaredNorm() - weights[vi];
        }
                
        XVec x = A.partialPivLu().solve(b);
    
        return std::make_tuple(x(DIM) + x.segment<DIM>(0).squaredNorm(), x.segment<DIM>(0));
        
    } else {
        
        XMat A = XMat::Zero(DIM+vertices.size(), DIM+vertices.size());
        
        A.block(0, 0, vertices.size(), DIM) = 2.0 * X.block(0, 0, DIM, vertices.size()).transpose();
        
        A.block(0, DIM, vertices.size(), 1) = XVec::Ones(vertices.size());
        
        DVec v0 = X.block<DIM, 1>(0, 0);
        for(std::size_t i = 1; i < vertices.size(); i++) {
            A.block<DIM,1>(vertices.size(), DIM+i) = v0 - X.block<DIM, 1>(0, i);
        }
        
        A.block<DIM, DIM>(vertices.size(), 0) = DMat::Identity();
        
        XVec b = XVec::Zero(DIM+vertices.size());
        for(std::size_t i = 0; i < vertices.size(); i++) {
            int vi = vertices[i];            
            b(i) = X.block<DIM, 1>(0, i).squaredNorm() - weights[vi];
        }
        
        b.segment<DIM>(vertices.size()) = v0;
                
        XVec x = A.partialPivLu().solve(b);
            
        return std::make_tuple(x(DIM) + x.segment<DIM>(0).squaredNorm(), x.segment<DIM>(0));
        
        
    }
    
}

template <int DIM> std::vector<double> calc_alpha_vals(CellComplex &comp, Embedding<DIM> &embed, 
                                                       std::vector<double> &weights, double alpha0 = -1.0) {
        
    
    DMat box_mat_inv = embed.box_mat.inverse();
    
    std::vector<double> alpha_vals(comp.ncells);
        
    for(int c = comp.ncells-1; c >= 0; c--) {
        
        // Skip vertices
        if(comp.get_dim(c) == 0) {
            alpha_vals[c] = alpha0;
            continue;
        }
        
        std::unordered_set<int> verts = comp.get_faces(c, 0);
        
        double alpha;
        DVec a;
        std::vector<int> tmp(verts.begin(), verts.end());
        std::tie(alpha, a) = calc_circumsphere<DIM>(tmp, embed, weights);
        
        alpha_vals[c] = alpha;
        
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
        auto cofaces = comp.get_cofaces(c, DIM);
                
        bool conflict = false;
        // For each coface get the list of its vertices and check if any fall inside the circumsphere of c
        for(auto cf: cofaces) {
            
            std::unordered_set<int> coface_verts = comp.get_faces(cf, 0);
            
            for(auto vi: coface_verts) {

                // If inside the circumsphere, mark as conflict and continue
                DVec x = embed.get_vpos(vi);
                
                if(embed.periodic) {
                    DVec bvec = x - box_mat_inv * a;
            
                    for(int d = 0; d < DIM; d++) {
                        if(std::fabs(bvec(d)) > 0.5) {
                            bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                        }
                    }
                    
                    x = a + bvec;
                }
                
                x = embed.box_mat * x;
                
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
            
            // if(comp.get_dim(c) < DIM+1) {
            //     py::print(comp.get_dim(c), alpha, a);
            // }
            
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
    
    return alpha_vals;
    
}


template <int DIM> Filtration construct_alpha_filtration(CellComplex &comp, Embedding<DIM> &embed, 
                                                         std::vector<double> &weights, double alpha0 = -1.0) {
    
    auto alpha_vals = calc_alpha_vals<DIM>(comp, embed, weights, alpha0);

    return construct_filtration(comp, alpha_vals);
    
}


//////////////////////////////////////////////////////////////////////////
//Deformations
//////////////////////////////////////////////////////////////////////////

template <int DIM> std::tuple<std::vector<double>, std::vector<double> > 
            calc_strains(RXVec disp, CellComplex &comp, Embedding<DIM> &embed) {
    
    
    std::vector<double> eps_shear(comp.ncells, 0.0);
    std::vector<double> eps_comp(comp.ncells, 0.0);
    
    for(int c = 0; c < comp.ncells; c++) {
        int d = comp.get_dim(c);
        
        if(d != DIM) {
            continue;
        }
        
        auto vset = comp.get_faces(c, 0);
        
        std::vector<int> verts(vset.begin(), vset.end());
        
        
        int vi = verts[0];
        
        DVec O = embed.get_vpos(vi);
        DVec uO = disp.segment<DIM>(DIM*vi);

        DMat X = DMat::Zero();
        DMat Y = DMat::Zero();
        
        
        for(int m = 0; m < DIM; m++) {  

            int vj = verts[1+m];
            DVec bvec = embed.get_vpos(vj); - O;
            
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
        
        eps_comp[c] = eps.trace();
                
        eps_shear[c] = (eps - DMat::Identity() * 1.0/DIM * eps_comp[c]).norm();
        
    }
    
    return std::forward_as_tuple(eps_comp, eps_shear);
    
    
    
    
}


template <int DIM> std::vector<double> calc_voronoi_D2min(RXVec disp, CellComplex &comp, Embedding<DIM> &embed) {
    
    
    std::vector<double> D2min(disp.size() / DIM); 
        
    for(int vi = 0; vi < comp.ncells; vi++) {        
        if(comp.get_dim(vi) != 0) {
            continue;
        }
                
        // get edges
        auto eset = comp.get_cofacets(vi);
        
        std::unordered_set<int> verts;
        
        for(auto e: eset) {
            auto vset = comp.get_facets(e);
            verts.insert(vset.begin(), vset.end());
        }
        
        verts.erase(vi);
                
        DVec O = embed.get_vpos(vi);
        DVec uO = disp.segment<DIM>(DIM*vi);

        DMat X = DMat::Zero();
        DMat Y = DMat::Zero();
        
        
        for(auto vj: verts) {  

            DVec bvec = embed.get_vpos(vj) - O;
            
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
            
            bvec = embed.box_mat * bvec;

            DVec du = disp.segment<DIM>(DIM*vj) - uO;
            
            X += du * bvec.transpose();
            Y += bvec * bvec.transpose();

        }
        
        DMat eps = X * Y.inverse(); 
        
        
        for(auto vj: verts) {  

            DVec bvec = embed.get_vpos(vj) - O;
            
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
            
            bvec = embed.box_mat * bvec;

            DVec du = disp.segment<DIM>(DIM*vj) - uO;
            
            D2min[vi] += (du - eps*bvec).squaredNorm();

        }
        
    }
    
    return D2min;
    
    
}


//////////////////////////////////////////////////////////////////////////
//Cell type counting
//////////////////////////////////////////////////////////////////////////


std::unordered_map<int, std::tuple<std::map<int, int>, std::map<int, int> > > 
    calc_radial_gap_distribution(std::vector<int> cell_list, std::vector<double> &alpha_vals, CellComplex &comp, int max_dist=-1, bool verbose=false) {
    
    // particle->(gaps[dist]->count, overlaps[dist]->count)
    std::unordered_map<int, std::tuple<std::map<int, int>, std::map<int, int> > > gap_distribution;
    
    int index = 0;
    for(auto c: cell_list) {
        
        if(verbose && index % 500 == 0) {
            py::print(index, "/", cell_list.size(), py::arg("flush")=true);
        }
        
        index++;
        
        std::map<int, int> gaps;
        std::map<int, int> overlaps;
        
        auto dists = calc_comp_point_dists(c, comp, max_dist);
        
        for(int i = 0; i < comp.ncells; i++) {
            if(comp.get_dim(i) != 1 || dists[i] <= 0) {
                continue;
            }
            
            double alpha = alpha_vals[i];
            
            if(alpha > 0.0) {
                gaps[dists[i]]++;
            } else {
                overlaps[dists[i]]++;
            }
        }
        
        gap_distribution[c] = std::make_tuple(gaps, overlaps);
        
    }
    
    return gap_distribution;
    
    
}

std::unordered_map<int, std::tuple<std::map<int, std::map<int, int> >, 
                                    std::map<int, std::map<int, int> >, 
                                    std::map<int, std::map<int, int> > > >
    calc_angular_gap_distribution(std::vector<int> cell_list, std::vector<double> &alpha_vals, CellComplex &comp, int max_dist=-1, bool verbose=false) {
    
    std::unordered_map<int, std::tuple<std::map<int, std::map<int, int> >, 
                                    std::map<int, std::map<int, int> >, 
                                    std::map<int, std::map<int, int> > > > gap_distribution;
    
    int index = 0;
    for(auto c: cell_list) {
        
        if(verbose && index % 500 == 0) {
            py::print(index, "/", cell_list.size(), py::arg("flush")=true);
        }
        
        index++;
                
        std::map<int, std::map<int, int> > gap_gap;
        std::map<int, std::map<int, int> > overlap_overlap;
        std::map<int, std::map<int, int> > gap_overlap;
        
        auto dists = calc_comp_point_dists(c, comp, max_dist);
        
        int max_rad_dist;
        if(max_dist == -1) {
            max_rad_dist = *std::max_element(dists.begin(), dists.end());
        } else {
            max_rad_dist = max_dist;
        }
        
        std::unordered_map<int, std::unordered_set<int> > dist_sets;
        for(int i = 0; i < comp.ncells; i++) {
            if(dists[i] <= max_rad_dist) {
                dist_sets[dists[i]].insert(i);
            }
        }
                
        for(int i = 1; i <= max_rad_dist; i++) {
                
            for(auto a: dist_sets[i]) {
                
                if(comp.get_dim(a) != 1) {
                    continue;
                }
                                
                double alphaa = alpha_vals[a];
                
                auto layer_dists = calc_comp_point_dists_search_zone(a, dist_sets[i], comp);
                            
                for(auto pair: layer_dists) {
                    
                    int b = pair.first;
                    
                    if(comp.get_dim(b) != 1) {
                        continue;
                    }
                    
                    double alphab = alpha_vals[b];
                    
                    if(alphaa > 0.0 && alphab > 0.0) {
                        gap_gap[i][pair.second/2]++;
                    } else if(alphaa < 0.0 && alphab < 0.0) {
                        overlap_overlap[i][pair.second/2]++;
                    } else {
                        gap_overlap[i][pair.second/2]++;
                    }

                }

            }
            
        }
        
        for(auto rad_pair: gap_gap) {
            for(auto ang_pair: rad_pair.second) {
                if(ang_pair.first != 0) {
                    gap_gap[rad_pair.first][ang_pair.first] /= 2;
                }
            }
        }
        
        for(auto rad_pair: overlap_overlap) {
            for(auto ang_pair: rad_pair.second) {
                if(ang_pair.first != 0) {
                    overlap_overlap[rad_pair.first][ang_pair.first] /= 2;
                }
            }
        }
        
        for(auto rad_pair: gap_overlap) {
            for(auto ang_pair: rad_pair.second) {
                gap_overlap[rad_pair.first][ang_pair.first] /= 2;
            }
        }
        
        gap_distribution[c] = std::make_tuple(gap_gap, overlap_overlap, gap_overlap);
        
        
    }
    
    return gap_distribution;
    
    
}




std::unordered_map<int,  std::unordered_map<int, std::map<std::tuple<std::string, std::string >, int> > >
    calc_radial_cycle_distribution(std::vector<int> cell_list, std::vector<double> &alpha_vals, std::vector<std::string> &vertex_types, CellComplex &comp, int max_dist=-1, bool verbose=false) {
    
    
    // particle->dist->(vertex_types, edge_types)->count
    std::unordered_map<int,  std::unordered_map<int, std::map<std::tuple<std::string, std::string >, int> > > cycle_distribution;
    
    int index = 0;
    for(auto c: cell_list) {
        
        if(verbose && index % 500 == 0) {
            py::print(index, "/", cell_list.size(), py::arg("flush")=true);
        }
        
        index++;
        
        
        auto dists = calc_comp_point_dists(c, comp, max_dist);
        
        for(int i = 0; i < comp.ncells; i++) {
            if(dists[i] <= 0) {
                continue;
            }
            
            
            auto verts = comp.get_faces(i, 0);
            auto edges = comp.get_faces(i, 1);
            
            std::vector<std::string> vtypes_list;
            for(auto v: verts) {
                if(v == c) {
                    vtypes_list.push_back("t");
                } else {
                    vtypes_list.push_back(vertex_types[v]);
                }
            }
            
            std::sort(vtypes_list.begin(), vtypes_list.end(), std::greater<std::string>());
            std::string vtypes = std::accumulate(vtypes_list.begin(), vtypes_list.end(), std::string(""));
            
            std::vector<std::string> etypes_list;
            for(auto e: edges) {
                
                double alpha = alpha_vals[e];
                if(alpha > 0.0) {
                    continue;
                }
                 
                std::vector<std::string> elabel_list;
                auto everts = comp.get_facets(e);
                for(auto ev: everts) {
                    
                    if(ev == c) {
                        elabel_list.push_back("t");
                    } else {
                        elabel_list.push_back(vertex_types[ev]);
                    }
                }
                
                std::sort(elabel_list.begin(), elabel_list.end(), std::greater<std::string>());
                std::string elabel = std::accumulate(elabel_list.begin(), elabel_list.end(), std::string(""));

                etypes_list.push_back(elabel);
                
            }
            
            std::sort(etypes_list.begin(), etypes_list.end(), std::greater<std::string>());
            std::string etypes = std::accumulate(etypes_list.begin(), etypes_list.end(), std::string(""));
            
            
            cycle_distribution[c][dists[i]][std::forward_as_tuple(vtypes, etypes)]++;
            
        }
                
    }
    
    return cycle_distribution;
    
    
}

    
    
#endif //ALPHACOMPLEX_HPP