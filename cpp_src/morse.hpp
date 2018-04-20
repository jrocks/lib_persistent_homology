#ifndef MORSE_HPP
#define MORSE_HPP
   
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
    
namespace py = pybind11;
    
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <numeric>
#include <utility>
#include <math.h>
#include <time.h>
    
#include "cell_complex.hpp"
#include "filtration.hpp"
    

std::unordered_set<int> get_lower_star(int alpha, bool co, StarFiltration &filt, CellComplex &comp, int target_dim) {
        
    int star_index = filt.get_filt_cell(alpha);
    
    std::unordered_set<int> star;
    
    std::queue<int> Q;
    Q.push(alpha);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        if(star_index == filt.get_filt_cell(a)) {
            
            // These are purposely backwards
            // Explore cofacets for star and facets for costar
            auto range = co ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto it = range.first; it != range.second; it++) {
                Q.push(*it);
            }

            if(target_dim == -1 || comp.get_dim(a) == target_dim) {
                star.insert(a);
            }
        }
    }
    
    return star;
 
}



std::tuple<py::array_t<int>, py::array_t<int>>  construct_discrete_gradient(StarFiltration &filt, CellComplex &comp) {     
        
    auto cmp = [&filt](const int& lhs, const int &rhs) {
        // Always pop largest value
        // or Always pop lowest value
        return (filt.co && filt(lhs, rhs)) || (!filt.co && filt(rhs, lhs));
        
    };
        
    py::array_t<int> V;
    V.resize({comp.ncells});    
    std::fill(V.mutable_data(), V.mutable_data() + V.size(), -1);
    
    auto Vbuf = V.mutable_unchecked<1>();
    
    
    for(int i = 0; i < comp.ncells; i++) {
        if((filt.co && comp.get_dim(i) != comp.dim) || (!filt.co && comp.get_dim(i) != 0)) {
            continue;
        }
        
        // get lower/upper star/costar
        std::unordered_set<int> star = get_lower_star(i, filt.co, filt, comp, -1);
        
        // record all unpaired facets/cofacets of each cell in star
        std::unordered_map<int, std::unordered_set<int> > unpaired;
        for(auto alpha: star) {
            unpaired[alpha] = std::unordered_set<int>();
            
            auto range = filt.co ? comp.get_cofacet_range(alpha) : comp.get_facet_range(alpha);
            for(auto it = range.first; it != range.second; it++) {
                int beta = *it;
                if(star.count(beta)) {
                    unpaired[alpha].insert(beta);
                }
            }
            
        }
        
        std::priority_queue<int, std::vector<int>, decltype(cmp)> PQone(cmp);
        std::priority_queue<int, std::vector<int>, decltype(cmp)> PQzero(cmp);
        
        for(auto alpha: star) {
            if(unpaired[alpha].empty()) {
                PQzero.push(alpha);
            } else if(unpaired[alpha].size() == 1) {
                PQone.push(alpha);
            }
        }
        
        while(!PQone.empty() || !PQzero.empty()) {
            while(!PQone.empty()) {
                int alpha = PQone.top();
                PQone.pop();
                
                if(unpaired[alpha].empty()) {
                    PQzero.push(alpha);
                    continue;
                }
                
                int beta = *(unpaired[alpha].begin());
                if(filt.co) {
                    Vbuf(alpha) = beta;
                } else {
                    Vbuf(beta)= alpha;
                }
                
                star.erase(alpha);
                star.erase(beta);
                unpaired.erase(alpha);
                
                for(auto gamma: star) {
                    if(unpaired[gamma].count(alpha) || unpaired[gamma].count(beta)) {
                        unpaired[gamma].erase(alpha);
                        unpaired[gamma].erase(beta);
                        
                        if(unpaired[gamma].size() == 1) {
                            PQone.push(gamma);
                        }
                    }
                }
            }
            
            
            if(!PQzero.empty()) {
                int alpha = PQzero.top();
                PQzero.pop();
                
                if(star.count(alpha)) {
                    Vbuf(alpha) = alpha;
                    
                    star.erase(alpha);
                    unpaired.erase(alpha);
                    
                    for(auto gamma: star) {
                        if(unpaired[gamma].count(alpha)) {
                            unpaired[gamma].erase(alpha);

                            if(unpaired[gamma].size() == 1) {
                                PQone.push(gamma);
                            }
                        }
                    }
                    
                }
            }
            
        }
                
    }
    
    py::array_t<int> coV;
    coV.resize({Vbuf.size()});    
    std::fill(coV.mutable_data(), coV.mutable_data() + coV.size(), -1);
    
    auto coVbuf = coV.mutable_unchecked<1>();
    
    for(int i = 0; i < Vbuf.size(); i++) {
        if(Vbuf(i) != -1) {
            coVbuf(Vbuf(i)) = i;
        }
    }
    
    return std::make_tuple(V, coV);
    
}



std::vector<std::tuple<int, int, int> > traverse_flow(int s, py::array_t<int> V, 
                                                       CellComplex &comp, bool co, bool coordinated) {
    
    auto Vbuf = V.unchecked<1>();
    
    std::unordered_map<int, int> nrIn;
    std::unordered_map<int, int> k;
    
    if(coordinated) {
        std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, V, comp, co, false);
        for(auto trip: traversal) {
            int c = std::get<2>(trip);
            nrIn[c] += 1;
        }
    }
    
    std::queue<int> Q;
    Q.push(s);
    
    std::unordered_set<int> seen;
    seen.insert(s);
    
    std::vector<std::tuple<int, int, int> > traversal;
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        auto range = co ? comp.get_cofacet_range(a) : comp.get_facet_range(a);
        for(auto it = range.first; it != range.second; it++) {
            int b = *it;
            if(Vbuf(b) != -1 && Vbuf(b) != a) {
                int c = Vbuf(b);
                traversal.emplace_back(a, b, c);
                
                if(!seen.count(c) && b != c) {
                    if(coordinated) {
                        k[c] += 1;
                        
                        if(k[c] != nrIn[c]) {
                            continue;
                        }
                        
                    }
                    
                    
                    Q.push(c);
                    seen.insert(c);
                    
                }   
            }
        }
        
    }
    
    return traversal;
    
}


std::vector<std::tuple<int, int, int> > calc_morse_boundary(int s, py::array_t<int> V, 
                                                             CellComplex &comp, bool co, bool oriented) {
    
    std::unordered_map<int, int> counts;
    counts[s] = 1;
    
    std::unordered_set<int> boundary;
    
    std::unordered_map<int, int> mult;
    if(oriented) {
        mult[s] = 1;
    }
    
    std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, V, comp, co, true);
    for(auto trip: traversal) {
        int a, b, c;
        std::tie(a, b, c) = trip;
        
        counts[c] += counts[a];
        
        if(b == c) {
            boundary.insert(c);
            
            // if oriented:
//                 mult[c] = mult.get(c, 0) + mult[a] * ( -coeff(a)[b] )
            
//         elif oriented:
//             mult[c] = mult.get(c, 0) + mult[a] * ( -coeff(a)[b] * coeff(c)[b] )
        }
    }
    
    std::vector<std::tuple<int, int, int> > morse_boundary;
    for(auto c: boundary) {
        if(oriented) {
            morse_boundary.emplace_back(c, counts[c], mult[c]);
        } else {
            morse_boundary.emplace_back(c, counts[c], counts[c] % 2);
        }
    }
    
    return morse_boundary;
    
}



CellComplex construct_morse_complex(py::array_t<int> V, CellComplex &comp, bool oriented) {
    
    auto Vbuf = V.unchecked<1>();
    
    CellComplex mcomp(comp.dim, false, oriented);
    
    
    // Morse complex is labeled according to corresponding cell in original complex
    // To get label of cell, consult original complex
    std::unordered_map<int, int> label_to_cell;
    std::vector<int> cell_to_label;
    
    int index = 0;
    for(int s = 0; s < comp.ncells; s++) {
        if(Vbuf(s) == s) {
            
            label_to_cell[s] = index;
            cell_to_label.push_back(s);
            
            index++;
            
        }
    }

    
    for(auto s: cell_to_label) {
        
        std::vector<int> facets;
        std::vector<int> coeffs;

        std::vector<std::tuple<int, int, int> > morse_boundary = calc_morse_boundary(s, V, comp, false, oriented);
        for(auto trip: morse_boundary) {
            int c, k, m;
            std::tie(c, k, m) = trip;

            facets.push_back(label_to_cell[c]);
            coeffs.push_back(m);
            
            if(m > 1 || m < -1) {
                py::print("Large Coefficient:", k, m, comp.get_dim(c), comp.get_dim(s));
            }
        }
        
        mcomp.add_cell(s, comp.get_dim(s), facets, coeffs);
    }
    
    mcomp.construct_cofacets();
    mcomp.make_compressed();
    return mcomp;
    
}


std::vector<std::tuple<int, int, int> > find_connections(int s, int t, py::array_t<int> V,  py::array_t<int> coV, 
                                                          CellComplex &comp) {
    
    std::unordered_set<int> active;
    active.insert(t);
        
    std::vector<std::tuple<int, int, int> > traversal = traverse_flow(t, coV, comp, true, false);
    for(auto trip: traversal) {
        int c = std::get<2>(trip);
        active.insert(c);
    }
        
    std::vector<std::tuple<int, int, int> > connections;
    
    traversal = traverse_flow(s, V, comp, false, false);
    for(auto trip: traversal) {
        int b = std::get<1>(trip);
        if(active.count(b)) {
            connections.push_back(trip);
        }
    }
        
    return connections;
    
}



std::unordered_set<int> convert_morse_to_real(std::unordered_set<int> &mfeature, py::array_t<int> V, py::array_t<int> coV,
                                                      CellComplex &comp, bool follow_bounds = true) {
    
    std::unordered_set<int> feature;
    
    for(auto s: mfeature) {
        
        bool co = !follow_bounds;
        
        if(comp.get_dim(s) == 0) {
            co = true;
        } else if(comp.get_dim(s) == comp.dim) {
            co = false;
        }
        
        feature.insert(s);
        std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, (co ? coV : V), comp, co, false);
        for(auto trip: traversal) {
            int b, c;
            std::tie(std::ignore, b, c) = trip;
            if(b != c) {
                feature.insert(c);
            }
        }
    }
    
    return feature;
    
}


// Update this to look more like extract_persistence_feature
// Also change so that it immediately spits out morse feature by default
std::unordered_set<int> extract_morse_feature(int i, int j, CellComplex &mcomp, StarFiltration &filt, int target_dim=-1, bool complement=false) {
        
    
    bool co = (mcomp.get_dim(i) != 0);
    
    if(target_dim == -1) {
        target_dim = co ? mcomp.dim : 0;
    }
    
    std::unordered_set<int> seen;
    std::queue<int> Q;
    if(!co) {
        seen.insert(i);
        Q.push(i);
    } else {
        seen.insert(j);
        Q.push(j);
    }
  
    int orderi = filt.get_total_order(mcomp.get_label(i));
    int orderj = filt.get_total_order(mcomp.get_label(j));
         
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        // py::print("a", a);
        
        for(auto b: get_star(a, co, mcomp, -1)) {
            
            // py::print("b", b);
                        
            if((!co && filt.get_total_order(mcomp.get_label(b)) >= orderj)
              || (co && filt.get_total_order(mcomp.get_label(b)) <= orderi)) {
                continue;
            }

            for(auto c: get_star(b, !co, mcomp, -1)) {
                // py::print("c", c);
                
                if((!co && filt.get_total_order(mcomp.get_label(c)) <= orderi)
                  || (co && filt.get_total_order(mcomp.get_label(c)) >= orderj)) {
                    continue;
                }
                
                if(!seen.count(c) && c != a) {
                    Q.push(c);
                    seen.insert(c);
                }
                
            }
        }
    }
    
    if(complement) {
        
        std::unordered_set<int> comp_seen;
        if(!co) {
            seen.insert(j);
            Q.push(j);
        } else {
            seen.insert(i);
            Q.push(i);
        }
        
        while(!Q.empty()) {
            int a = Q.front();
            Q.pop();

            // py::print("a", a);

            for(auto b: get_star(a, !co, mcomp, -1)) {

                // py::print("b", b);

                for(auto c: get_star(b, co, mcomp, -1)) {
                    // py::print("c", c);
                    
                    if((!co && filt.get_total_order(mcomp.get_label(c)) >= orderj)
                      || (co && filt.get_total_order(mcomp.get_label(c)) <= orderi)) {
                        continue;
                    }
                    
                    if(!seen.count(c) && !comp_seen.count(c) && c != a) {
                        Q.push(c);
                        comp_seen.insert(c);
                    }

                }
            }
        }
        
        seen = comp_seen;
        
    }
    
    std::unordered_set<int> feature;
    
    for(auto s: seen) {
        if(mcomp.get_dim(s) == target_dim) {
            feature.insert(mcomp.get_label(s));
        }
    }
    
    return feature;
    
}

    

/***************************************************
If partition is true:

Case 1: 2-cells are representative cells

Ex. Convert vertices to 2-cells

Each 2-cell should map to exactly one vertex
But each vertex can map to none or multiple 2-cells

Each 2-cell maps to lowest vertex in inclusion set
or if no vertices, then just to lowest boundary vertex

1. For each vertex find all 2-cells that are cofaces
2. For each 2-cell, find vertices in inclusion set
3. If nonzero vertices in inclusion set, check if current vertex is lowest
4. Otherwise check if current vertex is lowest boundary

This method ensure that vertices corresponding to minima always map to 2-cell


Case 2: vertices are representative cells

Ex. Convert 2-cells to vertices

Each vertex should map to exactly one 2-cell
But each 2-cell can map to none or multiple vertices

Each vertex maps to highest 2-cell in inclusion set
or if no 2-cells, then just to highest coboundary 2-cell

1. For each 2-cells find all vertices that are faces
2. For each vertex, find 2-cells in inclusion set
3. If nonzero 2-cells in inclusion set, check if current 2-cell is highest
4. Otherwise check if current 2-cell is highest cobounday

This method ensure that 2-cells corresponding to maxima always map to a vertex

If partition is false: 
Converts each cell to cells of new dimension that are in same inclusion set

***************************************************/
                                
std::unordered_set<int> change_feature_dim(std::unordered_set<int> &feature, int target_dim, StarFiltration &filt, CellComplex &comp, bool minimal = true) {
        
    std::unordered_set<int> new_feature;
    
    auto cmp = [&filt](const int &lhs, const int &rhs) {
        
        return (filt.co && filt(lhs, rhs)) || (!filt.co && filt(rhs, lhs));
        
    };
    
    for(auto s: feature) {
        
        bool co = (comp.get_dim(s) > target_dim);
        
        // Find all cofaces/faces of target_dim
        // Star is cofaces, costar is faces
        std::unordered_set<int> cofaces = get_star(s, co, comp, target_dim);
        
        // For each coface / face, check if it should be added to new feature
        for(auto c: cofaces) {
            
            // Don't uniquely partition the cofaces among the cells of the feature
            if(!minimal) {
                // Just add coface if is in same inclusion set
                if (filt.get_filt_cell(s) == filt.get_filt_cell(c)) {
                    new_feature.insert(c);
                }
                
            // Uniquely partition cofaces among the cells of the feature
            } else {

                // Get the star and upper costar of c
                std::unordered_set<int> star = get_star(c, !co, comp, comp.get_dim(s));
                std::unordered_set<int> ucostar;
                for(auto a: star) {
                    if(filt.get_filt_cell(c) == filt.get_filt_cell(a)) {
                        ucostar.insert(a);
                    }
                }

                // Upper costar is not empty and s is the minimum in upper costar
                // Or if upper costar is empty and s is the minimum in the star
                if( (!ucostar.empty() && (s == *std::min_element(ucostar.begin(), ucostar.end(), cmp)))
                  || (ucostar.empty() && (s == *std::min_element(star.begin(), star.end(), cmp))) ) {
                    new_feature.insert(c);
                }
                
            }
        }
                
    }
    
    
    return new_feature;
    
}




std::unordered_set<int> extract_morse_feature_to_real(int i, int j, CellComplex &mcomp, py::array_t<int> V, py::array_t<int> coV,
                                                      CellComplex &comp, StarFiltration &filt, bool complement=false, int target_dim=-1) {
    
    
    if(target_dim == -1) {
        target_dim = filt.fdim;
    }
    
    std::unordered_set<int> feature = extract_morse_feature(i, j, mcomp, filt, -1, complement);
    feature = convert_morse_to_real(feature, V, coV, comp);    
    feature = change_feature_dim(feature, target_dim, filt, comp, true);
    feature = comp.get_labels(feature);
    
    return feature;
    
}

/***************************************************
Finds the morse basins and converts to target_dim

***************************************************/
std::unordered_map<int, std::unordered_set<int>> find_morse_basins(CellComplex &mcomp, py::array_t<int> V, py::array_t<int> coV, 
                                                                   StarFiltration &filt, CellComplex &comp, int target_dim=-1) {
        
    std::unordered_map<int, std::unordered_set<int> > basins;
    
    if(target_dim == -1) {
        target_dim = filt.fdim;
    }
        
    for(int c = 0; c < mcomp.ncells; c++) {
        if(mcomp.get_dim(c) == 0) {
            // Cell in original complex
            int s = mcomp.get_label(c);
            
            std::unordered_set<int> mfeature;
            mfeature.insert(s);
            
            // Find the corresponding cells in real complex
            std::unordered_set<int> rfeature = convert_morse_to_real(mfeature, V, coV, comp);
                               
            // Change dimension of features to representative dimension
            std::unordered_set<int> feature = change_feature_dim(rfeature, target_dim, filt, comp, true);
            
            std::unordered_set<int> feature_labels = comp.get_labels(feature);
            
            // Find cell representing minimum in real complex
            // Since s is a minimum, this should just be the star cell
            int star_cell = comp.get_label(filt.get_filt_cell(s));
            basins.emplace(std::piecewise_construct, std::forward_as_tuple(star_cell), 
                           std::forward_as_tuple(feature_labels.begin(), feature_labels.end()));
                     
        }
    }
    
    return basins;
}


std::unordered_set<int> find_border(std::unordered_set<int> &feature, CellComplex &comp) {

    std::unordered_set<int> border;
    
    for(auto c: feature) {
        
        auto rangea = comp.get_cofacet_range(c);
        
        for(auto ita = rangea.first; ita != rangea.second; ita++) {
            
            int a = *ita;
            
            auto rangeb = comp.get_facet_range(a);
            for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
            
                int b = *itb;
                
                if(!feature.count(b)) {
                    border.insert(a);
                }
                
            
            }
            
        }
        
    }
    

    return border;
    
}

std::unordered_set<int> find_morse_basin_borders(CellComplex &mcomp, py::array_t<int> V, py::array_t<int> coV, 
                                                                   StarFiltration &filt, CellComplex &comp, int target_dim=-1) {
        
    std::unordered_set<int> basin_borders;
    
    if(target_dim == -1) {
        target_dim = filt.fdim;
    }
        
    for(int c = 0; c < mcomp.ncells; c++) {
        if(mcomp.get_dim(c) == 0) {
            // Cell in original complex
            int s = mcomp.get_label(c);
            
            std::unordered_set<int> mfeature;
            mfeature.insert(s);
            
            // Find the corresponding cells in real complex
            std::unordered_set<int> rfeature = convert_morse_to_real(mfeature, V, coV, comp);
                               
            std::unordered_set<int> rborder = find_border(rfeature, comp);
            
             // Change dimension of features to representative dimension
            std::unordered_set<int> border = change_feature_dim(rborder, target_dim, filt, comp, false);
            
            std::unordered_set<int> feature_labels = comp.get_labels(border);

            basin_borders.insert(feature_labels.begin(), feature_labels.end());
                     
        }
    }
    
    return basin_borders;
}

// Finds the morse skeleton in dimension sdim
std::unordered_set<int> find_morse_skeleton(CellComplex &mcomp, int sdim, py::array_t<int> V, py::array_t<int> coV,
                                            StarFiltration &filt, CellComplex &comp) {
        
    std::unordered_set<int> skeleton;
    
    for(int c = 0; c < mcomp.ncells; c++) {
        if(mcomp.get_dim(c) == sdim) {
            int s = mcomp.get_label(c);
                        
            std::unordered_set<int> mfeature;
            mfeature.insert(s);
            
            // Find the corresponding cells in real complex
            std::unordered_set<int> rfeature = convert_morse_to_real(mfeature, V, coV, comp);
            
            // Change dimension of features to representative dimension
            std::unordered_set<int> feature = change_feature_dim(rfeature, filt.fdim, filt, comp, false);
            
            std::unordered_set<int> feature_labels = comp.get_labels(feature);
            
            skeleton.insert(feature_labels.begin(), feature_labels.end());
        }
    }
    
    return skeleton;
}




std::pair<int, int> find_cancel_pair(int s, py::array_t<int> V, py::array_t<int> coV, 
                           StarFiltration &filt, CellComplex &comp) {
    
    if(comp.get_dim(s) == 0) {
        return std::make_pair(-1, -1);
    }

    bool single_Vpath = false;
    int close_alpha = 0;
    int close_alpha_time = 0;
    std::vector<std::tuple<int, int, int> > morse_boundary = calc_morse_boundary(s, V, comp, false, comp.oriented);
    for(auto trip: morse_boundary) {
        int c, k;
        std::tie(c, k, std::ignore) = trip;
        
        // if(s == 1073) {
        //     py::print("alpha", trip, filt.get_total_order(c));
        // }
        
        // If alpha is latest boundary so far, then record time and cell        
        if(k % 2 != 0 && (close_alpha == -1 || filt.get_total_order(c) > close_alpha_time)) {
            
            close_alpha = c;
            close_alpha_time = filt.get_total_order(c);
            
            // If there is only one V-path, then a cancellable close pair may exist
            single_Vpath = (k == 1);
        }
    }

    // if(s == 1073) {
    //     py::print("close_alpha", close_alpha, "single_path", single_Vpath);
    // }
        
    if(close_alpha == -1 || !single_Vpath) {
        return std::make_pair(-1, -1);
    }

    int close_beta = -1;
    int close_beta_time = 0;
    morse_boundary = calc_morse_boundary(close_alpha, coV, comp, true, comp.oriented);
    for(auto trip: morse_boundary) {
        int c, k;
        std::tie(c, k, std::ignore) = trip;


        if(k % 2 != 0 && (close_beta == -1 || filt.get_total_order(c) < close_beta_time)) {
            close_beta = c;
            close_beta_time = filt.get_total_order(c);
                        
        }
        
        
        
    }
        
    if(s == close_beta) {
        return std::make_pair(close_alpha, close_beta);   
    } else {
        return std::make_pair(-1, -1);
    }
    
    
}

void cancel_close_pair(std::pair<int, int> &pair, py::array_t<int> V, py::array_t<int> coV, CellComplex &comp) {
    
    auto Vbuf = V.mutable_unchecked<1>();
    auto coVbuf = coV.mutable_unchecked<1>();   
    
    int s = pair.second;
    int t = pair.first;

    std::vector<std::pair<int, int> > reverse_pairs;
    
    std::vector<std::tuple<int, int, int> > connections = find_connections(s, t, V, coV, comp);
    for(auto trip: connections) {
        int a, b;
        std::tie(a, b, std::ignore) = trip;
        reverse_pairs.emplace_back(a, b);
    }
    
    for(auto pair: reverse_pairs) {
        int a = pair.first;
        int b = pair.second;

        Vbuf(b) = a;
        coVbuf(a) = b;
        Vbuf(a) = -1;
        coVbuf(b) = -1;
    }
    
    
}

double find_join_threshold(std::vector<int> &verts, py::array_t<int> V, py::array_t<int> coV, 
                           StarFiltration &filt, CellComplex &comp, bool verbose) {
    
    
    py::array_t<int> V_tmp;
    py::array_t<int> coV_tmp;
    
    V_tmp.resize({comp.ncells});   
    coV_tmp.resize({comp.ncells});   
    
    std::copy(V.mutable_data(), V.mutable_data() + V.size(), V_tmp.mutable_data());
    std::copy(coV.mutable_data(), coV.mutable_data() + coV.size(), coV_tmp.mutable_data());
    
    

    auto Vbuf = V_tmp.mutable_unchecked<1>();
    std::vector<int> unpaired_crit_cells;
    for(int s = 0; s < Vbuf.size(); s++) {
        if(Vbuf(s) == s) {
            unpaired_crit_cells.push_back(s);
        }
    }
    

    auto cmp = [](const std::pair<double, std::pair<int, int> > &lhs, const std::pair<double, std::pair<int, int> > &rhs) {
        return lhs > rhs;
    };
    
    std::priority_queue<std::pair<double, std::pair<int, int> >, 
        std::vector<std::pair<double, std::pair<int, int> > >, decltype(cmp)> cancel_pairs(cmp);
    
    double threshold = 0.0;
    
    for(int n = 1; ; n++) {
        
        if(verbose) {
            py::print("Pass:", n, py::arg("flush")=true);
            py::print("Unpaired Critical Cells:", unpaired_crit_cells.size(), py::arg("flush")=true);
            py::print("Cancellable Pairs:", cancel_pairs.size(), py::arg("flush")=true);
        }
        
                
        // First check how many basins the vertices are divided into
        std::unordered_set<int> basins;
        for(auto s: verts) {
            
            // Check if vertex is critical
            if(Vbuf(s) == s) {
                basins.insert(s);
                continue;
            }
            
            // Need to make sure only goes one direction when starting from edge
        
            // If not, then start from adjacent edge and flow to critical vertex
            std::vector<std::tuple<int, int, int> > traversal = traverse_flow(Vbuf(s), V_tmp, comp, false, false);
            for(auto trip: traversal) {
                int b, c;
                std::tie(std::ignore, b, c) = trip;
                if(b == c) {
                    basins.insert(c);
                }
            }
        }
        
        
        if(basins.size() == 1) {
            return threshold;
        }

                
        // Otherwise pass through all unpaired critical cells and find cancellable pairs
        std::vector<int> remove;
        for(auto s: unpaired_crit_cells) {
                        
            auto cpair = find_cancel_pair(s, V_tmp, coV_tmp, filt, comp);
            if(cpair.first != -1) {
            
                // py::print(filt.get_time(cpair.second)-filt.get_time(cpair.first), cpair);
                
                cancel_pairs.emplace(std::abs(filt.get_time(cpair.second)-filt.get_time(cpair.first)), cpair);
                
                remove.push_back(cpair.first);
                remove.push_back(cpair.second);
                
            }   
        }
                
        for(auto s: remove) {
            unpaired_crit_cells.erase(std::find(unpaired_crit_cells.begin(), unpaired_crit_cells.end(), s));
        }        
                
        // Cancel critical pair with lowest persistence
        
        if(cancel_pairs.empty()) {
            break;
        }
        
        auto top = cancel_pairs.top();
        cancel_pairs.pop();
        auto cpair = top.second;
                
        if(top.first > threshold) {
            threshold = top.first;
        }
        
                
        cancel_close_pair(cpair, V_tmp, coV_tmp, comp);
                
    } 
    
    return -1.0;
    
}


double find_cancel_threshold(std::pair<int, int> pair, py::array_t<int> V, py::array_t<int> coV, 
                           StarFiltration &filt, CellComplex &comp, bool verbose) {
    
    
    py::array_t<int> V_tmp;
    py::array_t<int> coV_tmp;
    
    V_tmp.resize({comp.ncells});   
    coV_tmp.resize({comp.ncells});   
    
    std::copy(V.mutable_data(), V.mutable_data() + V.size(), V_tmp.mutable_data());
    std::copy(coV.mutable_data(), coV.mutable_data() + coV.size(), coV_tmp.mutable_data());
    
    

    auto Vbuf = V_tmp.mutable_unchecked<1>();
    std::vector<int> unpaired_crit_cells;
    for(int s = 0; s < Vbuf.size(); s++) {
        if(Vbuf(s) == s) {
            unpaired_crit_cells.push_back(s);
        }
    }
    

    auto cmp = [](const std::pair<double, std::pair<int, int> > &lhs, const std::pair<double, std::pair<int, int> > &rhs) {
        return lhs > rhs;
    };
    
    std::priority_queue<std::pair<double, std::pair<int, int> >, 
        std::vector<std::pair<double, std::pair<int, int> > >, decltype(cmp)> cancel_pairs(cmp);
    
    double threshold = 0.0;
    
    for(int n = 1; ; n++) {
        
        if(verbose) {
            py::print("Pass:", n, py::arg("flush")=true);
            py::print("Unpaired Critical Cells:", unpaired_crit_cells.size(), py::arg("flush")=true);
            py::print("Cancellable Pairs:", cancel_pairs.size(), py::arg("flush")=true);
        }
        
                
        // First check if pair is already a cancellable pair
        std::pair<int, int> cpair = find_cancel_pair(pair.second, V_tmp, coV_tmp, filt, comp);
        if(cpair == pair) {
            
            double p = std::abs(filt.get_time(cpair.second)-filt.get_time(cpair.first));
            
            return p > threshold ? p : threshold;
        }
                
        // Otherwise pass through all unpaired critical cells and find cancellable pairs
        std::vector<int> remove;
        for(auto s: unpaired_crit_cells) {
                        
            cpair = find_cancel_pair(s, V_tmp, coV_tmp, filt, comp);
            if(cpair.first != -1) {
            
                // py::print(filt.get_time(cpair.second)-filt.get_time(cpair.first), cpair);
                
                cancel_pairs.emplace(std::abs(filt.get_time(cpair.second)-filt.get_time(cpair.first)), cpair);
                
                remove.push_back(cpair.first);
                remove.push_back(cpair.second);
                
            }   
        }
                
        for(auto s: remove) {
            unpaired_crit_cells.erase(std::find(unpaired_crit_cells.begin(), unpaired_crit_cells.end(), s));
        }        
                
        // Cancel critical pair with lowest persistence
        
        if(cancel_pairs.empty()) {
            break;
        }
        
        auto top = cancel_pairs.top();
        cancel_pairs.pop();
        cpair = top.second;
        
        // py::print("top", top, py::arg("flush")=true);
        
        if(top.first > threshold) {
            threshold = top.first;
        }
        
                
        cancel_close_pair(cpair, V_tmp, coV_tmp, comp);
                
    } 
    
    return -1.0;
    
}

void simplify_morse_complex(double threshold, py::array_t<int> V, py::array_t<int> coV, StarFiltration &filt, CellComplex &comp, bool leq=true, bool verbose=false) {
    
    auto Vbuf = V.mutable_unchecked<1>();
    std::vector<int> unpaired_crit_cells;
    for(int s = 0; s < Vbuf.size(); s++) {
        if(Vbuf(s) == s) {
            unpaired_crit_cells.push_back(s);
        }
    }
    
    auto cmp = [](const std::pair<double, std::pair<int, int> > &lhs, const std::pair<double, std::pair<int, int> > &rhs) {
        return lhs > rhs;
    };
    
    std::priority_queue<std::pair<double, std::pair<int, int> >, 
        std::vector<std::pair<double, std::pair<int, int> > >, decltype(cmp)> cancel_pairs(cmp);
    
    for(int n = 1; ; n++) {
        
        if(verbose) {
            py::print("Pass:", n, py::arg("flush")=true);
            py::print("Critical Cells:", unpaired_crit_cells.size() + 2*cancel_pairs.size(), py::arg("flush")=true);
        }
        
        // Pass through all unpaired critical cells and find cancellable pairs
        std::vector<int> remove;
        for(auto s: unpaired_crit_cells) {
            
            auto cpair = find_cancel_pair(s, V, coV, filt, comp);            
            
            if(cpair.first != -1) {
                
                cancel_pairs.emplace(std::abs(filt.get_time(cpair.second)-filt.get_time(cpair.first)), cpair);
                
                remove.push_back(cpair.first);
                remove.push_back(cpair.second);
                
            }   
        }
        
        for(auto s: remove) {
            
            
            unpaired_crit_cells.erase(std::find(unpaired_crit_cells.begin(), unpaired_crit_cells.end(), s));
            
        }       
        
        
        // Cancel critical pairs with persistence below threshold
        
        if(verbose) {
            py::print("Cancellable Pairs:", cancel_pairs.size(), py::arg("flush")=true);
        }
        
        int n_cancel = 0;
        
        while(!cancel_pairs.empty() && ((leq &&  cancel_pairs.top().first <= threshold)
                  || (!leq && cancel_pairs.top().first < threshold))) {
            
            auto top = cancel_pairs.top();
            cancel_pairs.pop();
            auto cpair = top.second;
            cancel_close_pair(cpair, V, coV, comp);
                              
            n_cancel++;
            
        }
        
        if(verbose) {
            py::print("Cancelled Pairs:", n_cancel, py::arg("flush")=true);
            py::print("Remaining Critical Cells:", unpaired_crit_cells.size() + 2*cancel_pairs.size(), py::arg("flush")=true);
        }
        
        if(n_cancel == 0) {
            return;
        }
        
    } 
        
}


// void simplify_morse_complex(double threshold, py::array_t<int> V, py::array_t<int> coV, StarFiltration &filt, CellComplex &comp, bool leq=true, bool verbose=false) {
    
//     auto Vbuf = V.mutable_unchecked<1>();
//     std::vector<int> unpaired_crit_cells;
//     for(int s = 0; s < Vbuf.size(); s++) {
//         if(Vbuf(s) == s) {
//             unpaired_crit_cells.push_back(s);
//         }
//     }
    
//     auto cmp = [](const std::pair<double, std::pair<int, int> > &lhs, const std::pair<double, std::pair<int, int> > &rhs) {
//         return lhs > rhs;
//     };
    
//     std::priority_queue<std::pair<double, std::pair<int, int> >, 
//         std::vector<std::pair<double, std::pair<int, int> > >, decltype(cmp)> cancel_pairs(cmp);
    
//     for(int n = 1; ; n++) {
        
//         if(verbose) {
//             py::print("Pass:", n, py::arg("flush")=true);
//             py::print("Critical Cells:", unpaired_crit_cells.size() + 2*cancel_pairs.size(), py::arg("flush")=true);
//         }
        
//         // Pass through all unpaired critical cells and find cancellable pairs
//         std::vector<int> remove;
//         for(auto s: unpaired_crit_cells) {
            
//             auto cpair = find_cancel_pair(s, V, coV, filt, comp);            
            
//             if(cpair.first != -1) {
                
//                 py::print("Cancel", filt.get_time(cpair.second)-filt.get_time(cpair.first), cpair, comp.get_dim(cpair.first), comp.get_dim(cpair.second), py::arg("flush")=true);
            
//                 cancel_pairs.emplace(filt.get_time(cpair.second)-filt.get_time(cpair.first), cpair);
                
//                 remove.push_back(cpair.first);
//                 remove.push_back(cpair.second);
                
//             }   
//         }
        
//         py::print(remove, py::arg("flush")=true);
//         py::print(unpaired_crit_cells, py::arg("flush")=true);

        
//         for(auto s: remove) {
            
//             if(n == 3) {
//                 py::print(s, py::arg("flush")=true);
//             }
            
//             unpaired_crit_cells.erase(std::find(unpaired_crit_cells.begin(), unpaired_crit_cells.end(), s));
            
//         }       
        
        
//         // Cancel critical pairs with persistence below threshold
        
//         if(verbose) {
//             py::print("Cancellable Pairs:", cancel_pairs.size(), py::arg("flush")=true);
//         }
        
//         int n_cancel = 0;
        
//         while(!cancel_pairs.empty() && ((leq &&  cancel_pairs.top().first <= threshold)
//                   || (!leq && cancel_pairs.top().first < threshold))) {
            
//             auto top = cancel_pairs.top();
//             cancel_pairs.pop();
//             auto cpair = top.second;
//             cancel_close_pair(cpair, V, coV, comp);
                        
//             n_cancel++;
            
//         }
        
//         if(verbose) {
//             py::print("Cancelled Pairs:", n_cancel, py::arg("flush")=true);
//             py::print("Remaining Critical Cells:", unpaired_crit_cells.size() + 2*cancel_pairs.size(), py::arg("flush")=true);
//         }
        
//         if(n_cancel == 0) {
//             return;
//         }
        
//     } 
        
// }


    
#endif // MORSE_HPP