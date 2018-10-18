#ifndef SIMP_HPP
#define SIMP_HPP

#include "eigen_macros.hpp"
#include "cell_complex.hpp"
#include "filtration.hpp"
#include "morse_complex.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <numeric>
#include <utility>
#include <math.h>
#include <time.h>
   
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
    
namespace py = pybind11;
    

    



std::pair<int, int> find_cancel_pair(int s, RXiVec V, RXiVec coV, Filtration &filt, CellComplex &comp) {
    
    if(comp.get_dim(s) == 0) {
        return std::make_pair(-1, -1);
    }

    bool single_Vpath = false;
    int close_alpha = 0;
    int close_alpha_time = 0;
    std::vector<std::tuple<int, int, int> > morse_boundary = find_morse_boundary(s, V, comp, false, comp.oriented);
    for(auto trip: morse_boundary) {
        int c, k;
        std::tie(c, k, std::ignore) = trip;
        
        // If alpha is latest boundary so far, then record time and cell        
        if(k % 2 != 0 && (close_alpha == -1 || filt.get_order(c) > close_alpha_time)) {
            
            close_alpha = c;
            close_alpha_time = filt.get_order(c);
            
            // If there is only one V-path, then a cancellable close pair may exist
            single_Vpath = (k == 1);
        }
    }

    if(close_alpha == -1 || !single_Vpath) {
        return std::make_pair(-1, -1);
    }

    int close_beta = -1;
    int close_beta_time = 0;
    morse_boundary = find_morse_boundary(close_alpha, coV, comp, true, comp.oriented);
    for(auto trip: morse_boundary) {
        int c, k;
        std::tie(c, k, std::ignore) = trip;


        if(k % 2 != 0 && (close_beta == -1 || filt.get_order(c) < close_beta_time)) {
            close_beta = c;
            close_beta_time = filt.get_order(c);
                        
        }
        
        
        
    }
        
    if(s == close_beta) {
        return std::make_pair(close_alpha, close_beta);   
    } else {
        return std::make_pair(-1, -1);
    }
    
    
}

void cancel_close_pair(std::pair<int, int> &pair, RXiVec V, RXiVec coV, CellComplex &comp) {

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

        V(b) = a;
        coV(a) = b;
        V(a) = -1;
        coV(b) = -1;
    }
    
    
}

double find_join_threshold(std::vector<int> &verts, RXiVec V, RXiVec coV, 
                           Filtration &filt, CellComplex &comp, bool verbose) {
    
    
    XiVec V_tmp = V;
    XiVec coV_tmp = coV;

    std::vector<int> unpaired_crit_cells;
    for(int s = 0; s < V.size(); s++) {
        if(V(s) == s) {
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
            if(V_tmp(s) == s) {
                basins.insert(s);
                continue;
            }
                    
            // If not, then start from adjacent edge and flow to critical vertex
            std::vector<std::tuple<int, int, int> > traversal = traverse_flow(V_tmp(s), V_tmp, comp, false, false);
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
            
                // py::print(filt.get_func(cpair.second)-filt.get_func(cpair.first), cpair);
                
                cancel_pairs.emplace(std::abs(filt.get_func(cpair.second)-filt.get_func(cpair.first)), cpair);
                
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


double find_cancel_threshold(std::pair<int, int> pair, RXiVec V, RXiVec coV, 
                           Filtration &filt, CellComplex &comp, bool verbose) {
    
    
    XiVec V_tmp = V;
    XiVec coV_tmp = coV;

    std::vector<int> unpaired_crit_cells;
    for(int s = 0; s < V.size(); s++) {
        if(V(s) == s) {
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
            
            double p = std::abs(filt.get_func(cpair.second)-filt.get_func(cpair.first));
            
            return p > threshold ? p : threshold;
        }
                
        // Otherwise pass through all unpaired critical cells and find cancellable pairs
        std::vector<int> remove;
        for(auto s: unpaired_crit_cells) {
                        
            cpair = find_cancel_pair(s, V_tmp, coV_tmp, filt, comp);
            if(cpair.first != -1) {
            
                // py::print(filt.get_func(cpair.second)-filt.get_func(cpair.first), cpair);
                
                cancel_pairs.emplace(std::abs(filt.get_func(cpair.second)-filt.get_func(cpair.first)), cpair);
                
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

void simplify_morse_complex(double threshold, RXiVec V, RXiVec coV, Filtration &filt, CellComplex &comp, bool leq=true, bool verbose=false) {
    
    std::vector<int> unpaired_crit_cells;
    for(int s = 0; s < V.size(); s++) {
        if(V(s) == s) {
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
                
                cancel_pairs.emplace(std::abs(filt.get_func(cpair.second)-filt.get_func(cpair.first)), cpair);
                
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

#endif // SIMP_HPP