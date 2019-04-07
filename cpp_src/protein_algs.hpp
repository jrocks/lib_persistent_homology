#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <queue>
#include <tuple> 
#include <algorithm>

#include "eigen_macros.hpp"
#include "cell_complex.hpp"
#include "filtration.hpp"
#include "search.hpp"
#include "persistence_simplification.hpp"
#include "deformation.hpp"


#include <pybind11/pybind11.h>
namespace py = pybind11;



std::unordered_set<int> get_restricted_neighbors(int start, CellComplex &comp, std::unordered_set<int> &restricted, int max_dist) {
    
    auto cmp = [](const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) {
        return lhs > rhs;
    };
    
    std::priority_queue< std::pair<int, int>, 
        std::vector<std::pair<int, int> >, decltype(cmp)> PQ(cmp);
    
    PQ.emplace(0, start);
    
    std::unordered_map<int, int> dists;
    
    while(!PQ.empty()) {
        
        auto top = PQ.top();
        PQ.pop();
        
        int d = top.first;
        int vi = top.second;
        
        if(!dists.count(vi) || dists[vi] > d) {
            dists[vi] = d;
        } else {
            continue;
        }
        
        for(auto ei: comp.get_cofacets(vi)) {
            if(restricted.count(ei)) {
                for(auto vj: comp.get_facets(ei)) {
                    if(vj != vi && d < max_dist) {
                        PQ.emplace(d+1, vj);
                    }
                }
            }
        }
    }
    
    std::unordered_set<int> neigh;
    for(auto pair: dists) {
        neigh.insert(pair.first);
    }
    return neigh;
    
}



// template <int DIM> CellComplex shrink_alpha_complex(CellComplex &comp, Filtration &filt, int max_dist, double threshold=0.0, bool verbose=false) {
    
    
//     auto cmp = [&filt](const int lhs, const int rhs) {
//         return filt.get_order(lhs) < filt.get_order(rhs);
//     };
        
//     std::vector<int> cell_order(comp.ndcells[1]);
//     std::iota(cell_order.begin(), cell_order.end(), comp.dcell_range[1].first);
//     std::sort(cell_order.begin(), cell_order.end(), cmp);
    
    
//     std::unordered_set<int> restricted;
//     std::unordered_set<int> incomplete;
    
//     for(int i = 0; i < comp.ndcells[0]; i++) {
//         incomplete.insert(i);
//     }
    
//     int c;
//     for(int i = 0; i < comp.ndcells[1]; i++) {
//         c = cell_order[i];
        
        
//         if(verbose && restricted.size() % 1000 == 0) {
//             py::print(restricted.size(), incomplete.size(), filt.get_func(c), py::arg("flush")=true);
//         }
        
//         restricted.insert(c);
        
//         if(filt.get_func(c) < threshold) {
//             continue;
//         }
                
//         std::vector<int> complete;
//         for(int vi: incomplete) {
                        
//             auto neigh = get_restricted_neighbors(vi, comp, restricted, max_dist);
//             // Includes vi
//             if(neigh.size() >= DIM+1) {
//                 complete.push_back(vi);
//             }
            
//         }
        
//         for(int vi: complete) {
//             incomplete.erase(vi);
//         }
        
//         if(incomplete.empty() && filt.get_func(c) > threshold) {
//             break;
//         }
        
//     }
    
//     std::vector<int> rem_cells;
//     for(int i = comp.dcell_range[1].first; i < comp.ncells; i++) {
//         if(filt.get_order(i) > filt.get_order(c)) {
//             rem_cells.push_back(i);
//         }
//     }
    
//     if(verbose) {
//         py::print("Final alpha", filt.get_func(c), py::arg("flush")=true);
//     }
    
    
//     return prune_cell_complex(rem_cells, comp);
    

// }


template <int DIM> CellComplex shrink_alpha_complex(CellComplex &comp, Embedding<DIM> &embed) {
    
    std::vector<int> edge_order(comp.ndcells[1]);
    std::iota(edge_order.begin(), edge_order.end(), comp.dcell_range[1].first);
    
    std::vector<double> eq_lengths(comp.ndcells[1]);
    for(int c = comp.dcell_range[1].first; c < comp.dcell_range[1].second; c++) {
        
        auto facets = comp.get_facets(c);
        int vi = facets[0];
        int vj = facets[1];
        
        DVec bvec = embed.get_diff(embed.get_vpos(vj), embed.get_vpos(vi));
        
//         py::print(c, comp.get_label(c), bvec.norm());
        
        eq_lengths[comp.get_label(c)] = bvec.norm();
        
    }
    
    auto cmp = [&comp, &eq_lengths](const int &lhs, const int &rhs) {
        return eq_lengths[comp.get_label(lhs)] > eq_lengths[comp.get_label(rhs)];
    };
    
    std::sort(edge_order.begin(), edge_order.end(), cmp);
    
//     for(int i = 0; i < comp.ndcells[1]; i++) {
//         py::print(i, edge_order[i], eq_lengths[comp.get_label(edge_order[i])]);
//     }
        
    std::vector<int> rem_cells;
    for(int i = 0; i < comp.ndcells[1]; i++) {
        int c = edge_order[i];
        
        rem_cells.push_back(c);
        
        CellComplex comp_tmp = prune_cell_complex(rem_cells, comp);
        
        auto betti = calc_betti_numbers(comp_tmp);
        
        py::print(i, c, eq_lengths[comp.get_label(c)], betti);
        
        if(betti[0] > 1 || betti[comp.dim-1] > 0) {
            rem_cells.pop_back();
            break;
        }
        
        
    }
    
    
    return prune_cell_complex(rem_cells, comp);
    
    
}




template <int DIM> std::tuple<std::map<std::vector<int>, std::vector<int> >, 
std::map<std::vector<int>, std::vector<std::vector<int> > > > get_neighbor_grid(double grid_size, Embedding<DIM> &embed) {



    std::map<std::vector<int>, std::vector<int> > grid;
    
    DVec L = embed.box_mat.diagonal();    
    
    DiVec nbins;
    for(int d = 0; d < DIM; d++) {
        nbins(d) = int(L(d) / grid_size);
    }
    
    for(int vi = 0; vi < embed.NV; vi++) {
        
        std::vector<int> index(DIM);
        DVec pos = embed.get_pos(vi);
        for(int d = 0; d < DIM; d++) {            
            index[d] = int(pos(d) / L(d) * nbins(d));     
        }
        
        grid[index].push_back(vi);
        
    }
    
    
        
    std::map<std::vector<int>, std::vector<std::vector<int> > > grid_connections;
    for(auto gi: grid) {
        for(auto gj: grid) {
            
            DiVec diff(DIM);
            for(int d = 0; d < DIM; d++) {
                diff(d) = gi.first[d] - gj.first[d];
            }
            
            diff = diff.cwiseAbs();
            
            if(diff.maxCoeff() <= 1) {
                grid_connections[gi.first].push_back(gj.first);
            }
            
        }
    }
    
    return std::make_tuple(grid, grid_connections);
    
    
}



template <int DIM> XVec transform_coords(RXVec X, RDMat F, RDVec u) {
    
    XVec x = XVec::Zero(X.size());
    
    for(int vi = 0; vi < X.size() / DIM; vi++) {
        x.segment<DIM>(DIM*vi) = F * X.segment<DIM>(DIM*vi) + u;
    }
    
    return x;
    
}


template <int DIM> std::tuple<XVec, XVec> calc_neighborhood_D2min_strain(RXVec disp, Embedding<DIM> &embed, double max_dist, bool linear=true) {
    
    
    
    std::map<std::vector<int>, std::vector<int> > grid;
    std::map<std::vector<int>, std::vector<std::vector<int> > > grid_connections;
    
    std::tie(grid, grid_connections) = get_neighbor_grid<DIM>(max_dist, embed);
    
    
    
//     DVec L = embed.box_mat.diagonal();    
    
//     DiVec nbins;
//     for(int d = 0; d < DIM; d++) {
//         nbins(d) = int(L(d) / max_dist);
//     }
    
//     for(int vi = 0; vi < embed.NV; vi++) {
        
//         std::vector<int> index(DIM);
//         DVec pos = embed.get_pos(vi);
//         for(int d = 0; d < DIM; d++) {            
//             index[d] = int(pos(d) / L(d) * nbins(d));     
//         }
        
//         grid[index].push_back(vi);
        
//     }
    
    
        
//     std::map<std::vector<int>, std::vector<std::vector<int> > > grid_connections;
//     for(auto gi: grid) {
//         for(auto gj: grid) {
            
//             DiVec diff(DIM);
//             for(int d = 0; d < DIM; d++) {
//                 diff(d) = gi.first[d] - gj.first[d];
//             }
            
//             diff = diff.cwiseAbs();
            
//             if(diff.maxCoeff() <= 1) {
//                 grid_connections[gi.first].push_back(gj.first);
//             }
            
//         }
//     }
    
    
    DVec L = embed.box_mat.diagonal();    
    
    DiVec nbins;
    for(int d = 0; d < DIM; d++) {
        nbins(d) = int(L(d) / max_dist);
    }

    XVec D2min = XVec::Zero(embed.NV);
    XVec strains = XVec::Zero(embed.NV);

    for(int vi = 0; vi < embed.NV; vi++) {

        std::vector<int> index(DIM);
        DVec pos = embed.get_pos(vi);
        for(int d = 0; d < DIM; d++) {
            index[d] = int(pos(d) / L(d) * nbins(d));
        }
        
        
        DVec vposi = embed.get_vpos(vi);
        
        std::vector<int> verts;
        
        // Find neighbors
        for(auto gi: grid_connections[index]) {
            for(int vj: grid[gi]) {
                DVec bvec = embed.get_diff(vposi, embed.get_vpos(vj));
                
                if(bvec.norm() <= max_dist) {
                    verts.push_back(vj);
                }
            }
        }
        
        verts.erase(std::remove(verts.begin(), verts.end(), vi), verts.end());
        verts.insert(verts.begin(), vi);
        
        DMat F;
        std::tie(F, D2min[vi]) = calc_def_grad(verts, disp, embed, true);
        DMat eps  = def_grad_to_strain(F, linear);
        
        strains[vi] = eps.norm();
        
    }

    return std::make_tuple(D2min, strains);


}

template <int DIM> CellComplex get_contact_network(Embedding<DIM> &embed, RXVec rad) {
    
    
    double max_dist = 2.5*rad.maxCoeff();
    
    std::map<std::vector<int>, std::vector<int> > grid;
    std::map<std::vector<int>, std::vector<std::vector<int> > > grid_connections;
    
    std::tie(grid, grid_connections) = get_neighbor_grid<DIM>(max_dist, embed);
    
    
        
    DVec L = embed.box_mat.diagonal();    
    
    DiVec nbins;
    for(int d = 0; d < DIM; d++) {
        nbins(d) = int(L(d) / max_dist);
    }
    
    
    
    CellComplex comp(1, true, false);
    for(int vi = 0; vi < embed.NV; vi++) {
        std::vector<int> facets;
        std::vector<int> coeffs;
        comp.add_cell(vi, 0, facets, coeffs);
    }
    
    for(int vi = 0; vi < embed.NV; vi++) {
        
        std::vector<int> index(DIM);
        DVec pos = embed.get_pos(vi);
        for(int d = 0; d < DIM; d++) {
            index[d] = int(pos(d) / L(d) * nbins(d));
        }
        
        DVec vposi = embed.get_vpos(vi);
        
        // Find neighbors
        for(auto gi: grid_connections[index]) {
            for(int vj: grid[gi]) {
                
                if(vj <= vi) {
                    continue;
                }
                
                DVec bvec = embed.get_diff(vposi, embed.get_vpos(vj));
                
                if(bvec.norm() <= rad[vi] + rad[vj]) {
                    std::vector<int> facets;
                    facets.push_back(vi);
                    facets.push_back(vj);
                    
                    std::vector<int> coeffs;
                    comp.add_cell(comp.ndcells[1], 1, facets, coeffs);
                    
                }
            }
        }
        
    }
    
    comp.construct_cofacets();
    comp.make_compressed();
    
    return comp;
    
    
}



void merge_basins(std::vector<int> &saddles, RXiVec V, RXiVec coV, Filtration &filt, CellComplex &comp) {
    
    for(auto s: saddles) {
     
//             py::print("Saddle", i, j, s, filt.get_func(s), py::arg("flush")=true);

        auto morse_boundary = find_morse_boundary(s, V, comp, false, comp.oriented);

        std::vector<int> mverts;

        for(auto trip: morse_boundary) {
            int c, k;
            std::tie(c, k, std::ignore) = trip;

            if(k % 2 == 0) {
                py::print("Error", c, k);
            }

            mverts.push_back(c);

//                 py::print("vert", c, filt.get_func(c));

        }

        if(mverts.size() != 2) {
            py::print("Error", mverts, mverts.size());
        }


        auto pers_cmp = [s, &filt](const int &lhs, const int &rhs) {

            return std::abs(filt.get_func(s)-filt.get_func(lhs)) 
                < std::abs(filt.get_func(s)-filt.get_func(rhs));
        };

        int min_vert = *std::min_element(mverts.begin(), mverts.end(), pers_cmp);

//             py::print(filt.get_func(min_vert), filt.get_func(s), std::abs(filt.get_func(s)-filt.get_func(min_vert)));

        std::pair<int, int> cpair(min_vert, s);
        cancel_close_pair(cpair, V, coV, comp);

    }

    
    
    
}





std::vector<std::pair<int, int> > find_hinge_persistence_pairs(std::vector<std::pair<int, int> > &pairs, RXiVec V, RXiVec coV, Filtration &filt, CellComplex &comp, int n_basins=2, int min_size=1, bool reset=false, bool verbose=false) {
    
    
    XiVec V_tmp = V;
    XiVec coV_tmp = coV;
    
    auto cmp = [](const std::pair<double, std::pair<int, int> > &lhs, const std::pair<double, std::pair<int, int> > &rhs) {
        return lhs > rhs;
    };
    
    std::priority_queue<std::pair<double, std::pair<int, int> >, 
        std::vector<std::pair<double, std::pair<int, int> > >, decltype(cmp)> cancel_pairs(cmp);
    
            
    std::set<std::pair<int, int> > unsorted;
    for(auto pair: pairs) {
        unsorted.insert(pair);
    }
    
    
    if(verbose) {
        py::print("Pass:", 0, py::arg("flush")=true);
        py::print("Unsorted Pairs:", unsorted.size(), py::arg("flush")=true);
    }
    
    
    
    
    
    
     // Find morse basins
    auto mcomp = construct_morse_complex(V_tmp, comp);
    auto mbasins = find_morse_basins(coV_tmp, mcomp, comp);

    if(verbose) {
        py::print("Basins", mbasins.size(), py::arg("flush")=true);
    }  

    // Calculate basin sizes
    
    std::vector<int> basins;
    std::vector<int> sizes;
    std::vector<int> size_sort(mbasins.size());
    
    auto scmp = [&sizes] (const int lhs, const int rhs) {
        return sizes[lhs] > sizes[rhs];
    };
    
    for(auto pair: mbasins) {
        basins.push_back(pair.first);
        sizes.push_back(pair.second.size());
    }
    
    std::iota(size_sort.begin(), size_sort.end(), 0);
    std::sort(size_sort.begin(), size_sort.end(), scmp);
    
    
    


    bool analyze = true;
    
    std::unordered_set<int> hinge_basins;
    
    std::pair<int, int> prev_pair;
    
    for(int n = 1; ; n++) {
        
        
        if(analyze) {
            
            analyze = false;
            
            if(verbose) {
                py::print("Basins", mbasins.size(), py::arg("flush")=true);
            }  
            
            if(verbose) {
                if(mbasins.size() >= 5) {
                    py::print("Basin sizes", sizes[size_sort[0]], sizes[size_sort[1]], sizes[size_sort[2]], 
                              sizes[size_sort[3]], sizes[size_sort[4]], py::arg("flush")=true);
                } else if(mbasins.size() == 4) {
                    py::print("Basin sizes", sizes[size_sort[0]], sizes[size_sort[1]], sizes[size_sort[2]], 
                              sizes[size_sort[3]], py::arg("flush")=true);
                } else if(mbasins.size() == 3) {
                    py::print("Basin sizes", sizes[size_sort[0]], sizes[size_sort[1]], sizes[size_sort[2]], 
                              py::arg("flush")=true);
                } else if(mbasins.size() == 2) {
                    py::print("Basin sizes", sizes[size_sort[0]], sizes[size_sort[1]], py::arg("flush")=true);
                }
            }
            
            bool recorded = false;
            // If enough basins
            if((int)mbasins.size() >= n_basins) {
                
                // If nth basin is larger than min_size 
                // Or haven't recorded basins yet and there are only n left
                if(sizes[size_sort[n_basins-1]] >= min_size 
                   || (hinge_basins.empty() && (int)mbasins.size() == n_basins)) {
                    
                    // Record basins
                    hinge_basins.clear();
                    for(int i = 0; i < n_basins; i++) {
                        hinge_basins.insert(basins[size_sort[i]]);
                    }
                    
                    recorded = true;
                    
                    py::print("Recorded", hinge_basins);
                    
                }

            }
            
            // If a pair was previously cancelled (n > 1)
            // And basins have not been recorded this pass
            if(n > 1 && !recorded) {
                
                // Follow simplification path of basins and replace (but don't delete)
                if(hinge_basins.count(prev_pair.first)) {
                    
                    auto morse_boundary = find_morse_boundary(prev_pair.second, V_tmp, comp, false, comp.oriented);
                    for(auto trip: morse_boundary) {
                        int c = std::get<0>(trip);
                        if(!hinge_basins.count(c)) {
                            hinge_basins.erase(prev_pair.first);
                            hinge_basins.insert(c);
                            
                            py::print("Replaced", prev_pair.first, "with", c);
                            py::print(hinge_basins);
                        }
                    }
                    
                }   
                
            }
            
            
            
        }
        
        // Reset pairs with each pass
        if(reset) {
            // Clear priority queue of cancellable pairs
            while(!cancel_pairs.empty()) {
                auto top = cancel_pairs.top();
                cancel_pairs.pop();
                auto cpair = top.second;
                unsorted.insert(cpair);
            }  
            
        }

        // Identify cancellable pairs
        std::vector<std::pair<int, int> > remove;
        for(auto cpair: unsorted) {
        
            auto morse_boundary = find_morse_boundary(cpair.second, V_tmp, comp, false, comp.oriented);
                     
            for(auto trip: morse_boundary) {
                int c = std::get<0>(trip);

                if(c == cpair.first) {
                    
                    cancel_pairs.emplace(std::abs(filt.get_func(cpair.second)-filt.get_func(cpair.first)), cpair);
                    remove.push_back(cpair);
                    break;
                }
            }
        }
        
        for(auto cpair: remove) {
            unsorted.erase(cpair);
        }
    
        
         if(verbose) {
            py::print("Pass:", n, py::arg("flush")=true);
            py::print("Cancellable Pairs:", cancel_pairs.size(), py::arg("flush")=true);
            py::print("Unsorted Pairs:", unsorted.size(), py::arg("flush")=true);
             
        }
    
        if(!cancel_pairs.empty()) {
            auto top = cancel_pairs.top();
            cancel_pairs.pop();

            auto t = top.first;
            auto cpair = top.second;
            
            cancel_close_pair(cpair, V_tmp, coV_tmp, comp);
            
            if(verbose) {
                py::print("Threshold:", t, cpair, py::arg("flush")=true);
            }
            
            
            if(comp.get_dim(cpair.first) == 0) {
               
                analyze = true;
                prev_pair = cpair;
                
                
                auto morse_boundary = find_morse_boundary(cpair.second, V_tmp, comp, false, comp.oriented);
                int c = std::get<0>(morse_boundary[0]);
                mbasins[c].insert(mbasins[cpair.first].begin(), mbasins[cpair.first].end());
                mbasins.erase(cpair.first);
                
                if(verbose) {
                    py::print("Joining", cpair.first, c, py::arg("flush")=true);
                }
                
                basins.clear();
                sizes.clear();
                size_sort.resize(mbasins.size());
                for(auto pair: mbasins) {
                    basins.push_back(pair.first);
                    sizes.push_back(pair.second.size());
                }

                std::iota(size_sort.begin(), size_sort.end(), 0);
                std::sort(size_sort.begin(), size_sort.end(), scmp);
                
                
                
            }
            
        } else {
            break;
        }
        
    }
    
    
    std::vector<std::pair<int, int> > hinge_pairs;
    for(auto pair: pairs) {
        if(hinge_basins.count(pair.first)) {
            hinge_pairs.push_back(pair);
        }
    }

    
    return hinge_pairs;
    
}





void simplify_morse_complex(std::vector<std::pair<int, int> > &pairs, RXiVec V, RXiVec coV, Filtration &filt, CellComplex &comp, bool reset=false, bool verbose=false) {
    
    
    
    auto cmp = [](const std::pair<double, std::pair<int, int> > &lhs, const std::pair<double, std::pair<int, int> > &rhs) {
        return lhs > rhs;
    };
    
    std::priority_queue<std::pair<double, std::pair<int, int> >, 
        std::vector<std::pair<double, std::pair<int, int> > >, decltype(cmp)> cancel_pairs(cmp);
    
            
    std::set<std::pair<int, int> > unsorted;
    for(auto pair: pairs) {
        unsorted.insert(pair);
    }
    
    
    if(verbose) {
        py::print("Pass:", 0, py::arg("flush")=true);
        py::print("Unsorted Pairs:", unsorted.size(), py::arg("flush")=true);
    }

    std::pair<int, int> hinge_pair;
    
    for(int n = 1; ; n++) {
        
        // Reset pairs with each pass
        if(reset) {
            // Clear priority queue of cancellable pairs
            while(!cancel_pairs.empty()) {
                auto top = cancel_pairs.top();
                cancel_pairs.pop();
                auto cpair = top.second;
                unsorted.insert(cpair);
            }  
            
        }

        
        std::vector<std::pair<int, int> > remove;
        for(auto cpair: unsorted) {
        
            std::vector<std::tuple<int, int, int> > morse_boundary = find_morse_boundary(cpair.second, V, comp, false, comp.oriented);
            
            for(auto trip: morse_boundary) {
                int c, k;
                std::tie(c, k, std::ignore) = trip;

                if(c == cpair.first) {
                    
                    cancel_pairs.emplace(std::abs(filt.get_func(cpair.second)-filt.get_func(cpair.first)), cpair);
                    remove.push_back(cpair);
                    break;
                }
            }
        }
        
        for(auto cpair: remove) {
            unsorted.erase(cpair);
        }
    
        
        if(verbose) {
            py::print("Pass:", n, py::arg("flush")=true);
            py::print("Cancellable Pairs:", cancel_pairs.size(), py::arg("flush")=true);
            py::print("Unsorted Pairs:", unsorted.size(), py::arg("flush")=true);
             
        }
    
        if(!cancel_pairs.empty()) {
            auto top = cancel_pairs.top();
            cancel_pairs.pop();

            auto t = top.first;
            auto cpair = top.second;
            
            
            if(verbose) {
                py::print("Cancelling:", cpair, t, py::arg("flush")=true);
            }
            
            
            cancel_close_pair(cpair, V, coV, comp);
            
        } else {
            break;
        }
        
    }
    
}




#endif // PROTEIN_HPP