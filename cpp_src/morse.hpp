#ifndef MORSE
#define MORSE
   
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
    
#include "cell_complex.hpp"
#include "filtration.hpp"
    

std::unordered_set<int> get_lower_star(int alpha, bool co, Filtration &filt, CellComplex &comp, int target_dim) {
        
    int star_index = filt.get_star(alpha);
    
    std::unordered_set<int> star;
    
    std::queue<int> Q;
    Q.push(alpha);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        if(star_index == filt.get_star(a)) {
            
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



std::tuple<py::array_t<int>, py::array_t<int>>  construct_discrete_gradient(Filtration &filt, CellComplex &comp, bool dual) {     
    
    auto cmp = [&filt, dual](const int& lhs, const int &rhs) {
        if(dual) {                    
            return (-filt.get_filtration_order(lhs) > -filt.get_filtration_order(rhs));
        } else {
            return (filt.get_filtration_order(lhs) > filt.get_filtration_order(rhs));
        }
    };
        
    py::array_t<int> V;
    V.resize({comp.ncells});    
    std::fill(V.mutable_data(), V.mutable_data() + V.size(), -1);
    
    auto Vbuf = V.mutable_unchecked<1>();
    
    
    for(int i = 0; i < comp.ncells; i++) {
        if((dual && comp.get_dim(i) != comp.dim) || (!dual && comp.get_dim(i) != 0)) {
            continue;
        }
        
        // get lower/upper star/costar
        std::unordered_set<int> star = get_lower_star(i, dual, filt, comp, -1);
        
        // record all unpaired facets/cofacets of each cell in star
        std::unordered_map<int, std::unordered_set<int> > unpaired;
        for(auto alpha: star) {
            unpaired[alpha] = std::unordered_set<int>();
            
            auto range = dual ? comp.get_cofacet_range(alpha) : comp.get_facet_range(alpha);
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
                if(dual) {
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


void simplify_morse_complex(double threshold, py::array_t<int> V, py::array_t<int> coV, Filtration &filt, CellComplex &comp, bool verbose=false) {
    
    auto Vbuf = V.mutable_unchecked<1>();
    auto coVbuf = coV.mutable_unchecked<1>();    
    std::vector<int> crit_cells;
    for(int s = 0; s < Vbuf.size(); s++) {
        if(Vbuf(s) == s) {
            crit_cells.push_back(s);
        }
    }
    
    std::sort(crit_cells.begin(), crit_cells.end(),
       [&filt](int lhs, int rhs) {return filt.get_filtration_order(lhs) < filt.get_filtration_order(rhs);});
    
    for(int n = 1; ; n++) {
        
        
        if(verbose) {
            py::print("Pass:", n);
            py::print("Critical Cells:", crit_cells.size());
        }
        
        int n_cancel = 0;
        std::vector<std::pair<int, int> >  cancel_pairs;
        for(auto s: crit_cells) {
            if(comp.get_dim(s) == 0) {
                continue;
            }
            
            int close_alpha = -1;
            double close_alpha_time = 0.0;
            std::vector<std::tuple<int, int, int> > morse_boundary = calc_morse_boundary(s, V, comp, false, comp.oriented);
            for(auto trip: morse_boundary) {
                int c, k;
                std::tie(c, k, std::ignore) = trip;
                
                if(k == 1) {
                    if(close_alpha == -1 || filt.get_filtration_order(c) > close_alpha_time) {
                        close_alpha = c;
                        close_alpha_time = filt.get_filtration_order(c);
                    }
                }
            }
            
            if(close_alpha == -1) {
                continue;
            }
            
            int close_beta = -1;
            double close_beta_time = 0.0;
            morse_boundary = calc_morse_boundary(close_alpha, coV, comp, true, comp.oriented);
            for(auto trip: morse_boundary) {
                int c, k;
                std::tie(c, k, std::ignore) = trip;
                
                if(k == 1) {
                    if(close_beta == -1 || filt.get_filtration_order(c) < close_beta_time) {
                        close_beta = c;
                        close_beta_time = filt.get_filtration_order(c);
                    }
                }
            }
            
            if(s == close_beta) {
                n_cancel++;
                if(filt.get_time(close_beta)-filt.get_time(close_alpha) <= threshold) {
                    cancel_pairs.emplace_back(close_alpha, close_beta); 
                }
            }
                
        }
        
        
        if(verbose) {
            py::print("Cancellable Pairs:", n_cancel);
        }
        
        for(auto pair: cancel_pairs) {
            
            // s is the higher dimensional cell
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
            
            crit_cells.erase(std::find(crit_cells.begin(), crit_cells.end(), s));
            crit_cells.erase(std::find(crit_cells.begin(), crit_cells.end(), t));
            
        }
        
        if(verbose) {
            py::print("Cancelled Pairs:", cancel_pairs.size());
            py::print("Remaining Critical Cells:", crit_cells.size());
        }
        
        if(cancel_pairs.empty()) {
            break;
        }
        
    } 
    
}



std::unordered_set<int> convert_morse_to_real_complex(std::unordered_set<int> &mfeature, py::array_t<int> V,
                                                      CellComplex &comp, bool co) {
    
    std::unordered_set<int> feature;
    
    for(auto s: mfeature) {
        feature.insert(s);
        std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, V, comp, co, false);
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


std::unordered_map<int, std::unordered_set<int>> find_basins(CellComplex &mcomp, py::array_t<int> coV, Filtration &filt, CellComplex &comp, bool dual) {
        
    std::unordered_map<int, std::unordered_set<int> > basins;
    
    for(int c = 0; c < mcomp.ncells; c++) {
        if(mcomp.get_dim(c) == 0) {
            // Cell in original complex
            int s = mcomp.get_label(c);
                       
            
            std::unordered_set<int> mfeature;
            mfeature.insert(s);
            
            std::unordered_set<int> feature = convert_morse_to_real_complex(mfeature, coV, comp, true);
                        
            std::unordered_set<int> pixels = convert_to_pixels(feature, filt, comp, dual);
                        
            int pix = comp.get_label((dual ? filt.get_star(s) : s)); 
            basins.emplace(std::piecewise_construct, std::forward_as_tuple(pix), 
                           std::forward_as_tuple(pixels.begin(), pixels.end()));
                     
        }
    }
    
    return basins;
}


std::unordered_set<int> find_morse_skeleton(CellComplex &mcomp, int d, py::array_t<int> V, Filtration &filt, CellComplex &comp, bool dual) {
        
    std::unordered_set<int> skeleton;
    
    for(int c = 0; c < mcomp.ncells; c++) {
        if(mcomp.get_dim(c) == d) {
            int s = mcomp.get_label(c);
                        
            std::unordered_set<int> mfeature;
            mfeature.insert(s);
            
            std::unordered_set<int> feature = convert_morse_to_real_complex(mfeature, V, comp, false);
            
            std::unordered_set<int> pixels = convert_to_pixels(feature, filt, comp, dual);
            
            skeleton.insert(pixels.begin(), pixels.end());
            
        }
    }
    
    
    return skeleton;
}


std::unordered_set<int> extract_persistence_feature(int i, int j, CellComplex &mcomp, py::array_t<int> V, py::array_t<int> coV, Filtration &filt, CellComplex &comp) {
        
    
    std::unordered_set<int> seen;
    std::queue<int> Q;
    if(mcomp.get_dim(i) == 0) {
        seen.insert(i);
        Q.push(i);
    } else {
        seen.insert(j);
        Q.push(j);
    }
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        auto rangea = (mcomp.get_dim(i) == 0) ? mcomp.get_cofacet_range(a) : mcomp.get_facet_range(a);
        for(auto ita = rangea.first; ita != rangea.second; ita++) {

            int b = *ita;
            
            // maybe change to get_filtration_order instead?
            if((mcomp.get_dim(i) == 0 && filt.get_time(mcomp.get_label(b)) >= filt.get_time(mcomp.get_label(j))) 
              ||(mcomp.get_dim(i) != 0 && filt.get_time(mcomp.get_label(b)) <= filt.get_time(mcomp.get_label(i)))) {
                continue;
            }

            auto rangeb = (mcomp.get_dim(i) == 0) ? mcomp.get_facet_range(b) : mcomp.get_cofacet_range(b);
            for(auto itb = rangeb.first; itb != rangeb.second; itb++) {

                int c = *itb;
                if(!seen.count(c) && c != a) {
                    Q.push(c);
                    seen.insert(c);
                }
                
            }
        }
    }
    
    std::unordered_set<int> mfeature;
    for(auto s: seen) {
        mfeature.insert(mcomp.get_label(s));
    }
    
    if(mcomp.get_dim(i) == 0) {
        return convert_morse_to_real_complex(mfeature, coV, comp, true);     
    } else {
        return convert_morse_to_real_complex(mfeature, V, comp, false);
    }
        
}


    
#endif // MORSE