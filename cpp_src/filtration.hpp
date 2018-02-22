#ifndef FILTRATION
#define FILTRATION
    
    
#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <numeric>
    
#include "cell_complex.hpp"
    
class Filtration {
    
public:
    
    int ncells;
    int nprimary_cells;
    
private:
    
    
    
    
    
    // Assumes primary cells are first, when that is not necessarily true
    // Need to set primary
    
    // The problem is that time and primary order are indexed from 0 to nprimary_cell, while actual cells may have larger labels
    
    // Mapping of each cell to primary cell index corresponding to star_cell
    // Use this to access time, primary_order and star_cell
    std::vector<int> primary_index;
    
    // Filtration time of primary cells
    std::vector<double> time;
    // Strict ordering index of primary cells
    std::vector<int> primary_order;
    // Primary cell whose lower star/upper costar each cell belongs to
    std::vector<int> star_cell;
    
    // Full filtration ordering
    std::vector<int> filtration_order;
    
    
public:
    
    Filtration(CellComplex &comp, std::vector<double> &time, bool dual) : 
        ncells(comp.ncells), nprimary_cells(time.size()), primary_index(comp.ncells),
        time(time), primary_order(time.size()), star_cell(time.size()),
        filtration_order(ncells) {
            
            int target_dim = dual ? comp.dim : 0;
            int pindex = 0;
            for(int i = 0; i < ncells; i++) {
                if(comp.get_dim(i) == target_dim) {
                    primary_index[i] = pindex;
                    star_cell[pindex] = i;
                    pindex++;
                }
            }
            
        }
    
    void set_primary_order(int alpha, int order) {
        primary_order[primary_index[alpha]] = order;
    }
    
    void set_star(int alpha, int beta) {
        primary_index[alpha] = primary_index[beta];
    }
    
    void set_filtration_order(int alpha, int order) {
        filtration_order[alpha] = order;
    }
    
    
    double get_time(int alpha) {
        return time[primary_index[alpha]];
    }
    
    int get_primary_order(int alpha) {
        return primary_order[primary_index[alpha]];
    }
    
    int get_star(int alpha) {
        return star_cell[primary_index[alpha]];
    }
    
    int get_filtration_order(int alpha) {
        return filtration_order[alpha];
    }
    
    std::vector<int> get_filtration() {
        std::vector<int> filtration(ncells);
        std::iota(filtration.begin(), filtration.end(), 0);
        std::sort(filtration.begin(), filtration.end(), 
                 [this](const int &lhs, const int &rhs) {return this->filtration_order[lhs] < this->filtration_order[rhs];});
        
        return filtration;
    }
    
};


void perform_watershed_transform(Filtration &filt, CellComplex &comp, bool dual) {    
    
    std::vector<bool> submerged(filt.nprimary_cells, false);
    
    std::vector<int> pcell_argsort;
    int target_dim = dual ? comp.dim : 0;
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) == target_dim) {
            pcell_argsort.push_back(i);
        }
    }
    std::sort(pcell_argsort.begin(), pcell_argsort.end(),
       [&filt](int lhs, int rhs) {return filt.get_time(lhs) <= filt.get_time(rhs);});
    
    unsigned int ti = 0;
    while(ti < pcell_argsort.size()) {
        
        std::unordered_set<int> plateau;
        
        unsigned int tj = ti;
        while(true) {
            int vi = pcell_argsort[tj];
            double t = filt.get_time(vi);
            plateau.insert(vi);
            
            if((tj+1 == pcell_argsort.size()) || (t != filt.get_time(pcell_argsort[tj+1]))) {
                break;
            }
            
            tj++;
        }
        
        if(plateau.size() == 1) {
            int vi = *(plateau.begin());
            submerged[vi] = true;
            filt.set_primary_order(vi, ti);
            ti++;
            
            continue;
        }
        
        std::unordered_map<int, double> dist;
        std::queue<int> Q;
        for(auto a: plateau) {
            
            auto rangea = dual ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto ita = rangea.first; ita != rangea.second; ita++) {
                
                int b = *ita;
                
                auto rangeb = dual ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
                    
                    int c = *itb;
                    
                    if(submerged[c]) {
                        
                        double new_dist = 1.0;
                        
                        if(!dist.count(a) || new_dist < dist[a]) {
                            dist[a] = new_dist;
                            Q.push(a);
                        }
                        
                    }

                }
                
            }
            
        }
        
        while(!Q.empty()) {
            int a = Q.front();
            Q.pop();
            
            // This differs from the python version;
            double current_dist = dist[a];
            
            if(!plateau.count(a)) {
                continue;
            }
            
            
            plateau.erase(a);
            submerged[a] = true;
            filt.set_primary_order(a, ti);
            ti++;
            
            auto rangea = dual ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto ita = rangea.first; ita != rangea.second; ita++) {
                
                int b = *ita;
                
                auto rangeb = dual ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
                    
                    int c = *itb;
                    
                    if(plateau.count(c)) {
                        
                        double new_dist = current_dist + 1.0;
                        
                        if(!dist.count(c) || new_dist < dist[c]) {
                            dist[c] = new_dist;
                            Q.push(c);
                        }
                    }
                    
                }
                
            }
            
            
            
        }
        
        
        while(!plateau.empty()) {
            int s = *(plateau.begin());
            plateau.erase(s);
            
            Q.push(s);
            
            while(!Q.empty()) {
                int a = Q.front();
                Q.pop();

                // This differs from the python version;
                double current_dist = dist[a];

                plateau.erase(a);
                submerged[a] = true;
                filt.set_primary_order(a, ti);
                ti++;



                auto rangea = dual ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
                for(auto ita = rangea.first; ita != rangea.second; ita++) {

                    int b = *ita;

                    std::vector<int>::iterator beginb;
                    std::vector<int>::iterator endb;

                    auto rangeb = dual ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                    for(auto itb = rangeb.first; itb != rangeb.second; itb++) {

                        int c = *itb;

                        if(plateau.count(c)) {
                            double new_dist = current_dist + 1.0;

                            if(!dist.count(c) || new_dist < dist[c]) {
                                dist[c] = new_dist;
                                Q.push(c);
                            }
                        }

                    }
                }
            }
            
        }
        
    }    
    
}




void construct_filtration(Filtration &filt, CellComplex &comp, bool dual) {
    
    std::vector<int> cell_to_star(comp.ncells);
    std::vector<std::vector<int> > lex_val(comp.ncells, std::vector<int>());
    
    int target_dim = dual ? comp.dim : 0;
    for(int i = 0; i < comp.ncells; i++) {
                  
        std::unordered_set<int> star = get_star(i, !dual, comp, target_dim);
        int star_cell = *(star.begin());
                
        for(auto alpha: star) {
                                    
            lex_val[i].push_back(filt.get_primary_order(alpha));
            
            if((dual && filt.get_primary_order(alpha) < filt.get_primary_order(star_cell))
               || (!dual && filt.get_primary_order(alpha) > filt.get_primary_order(star_cell))) {
                
                star_cell = alpha;
            }
        }
                
        cell_to_star[i] = star_cell;
        
        // Sort from high to low
        std::sort(lex_val[i].begin(), lex_val[i].end(), std::greater<int>());  
        
    }    
    
    auto cmp = [&comp, &lex_val, dual](const int &lhs, const int &rhs) {
        // If dual, first compare upper costars (smallest goes first)
        if(dual && lex_val[lhs].back() != lex_val[rhs].back()) {
            return lex_val[lhs].back() < lex_val[rhs].back();
        // If not dual, first compare lower stars (smallest goes first)
        } else if(!dual && lex_val[lhs].front() != lex_val[rhs].front()) {
            return lex_val[lhs].front() < lex_val[rhs].front();
        // Second compare dimensions (smallest goes first)
        } else if(comp.get_dim(lhs) != comp.get_dim(rhs)) {
            return comp.get_dim(lhs) < comp.get_dim(rhs);
        // Third compare lexicographic values
        } else if(lex_val[lhs] != lex_val[rhs]) {
            return lex_val[lhs] < lex_val[rhs];
            
        // Finally, if cells have identical lexicographic orderings, 
        // then sort by raw cell index
        } else {
            return lhs < rhs;
        }
            
            
            
//         // Third compare length of lexicographic values (smallest goes first)
//         } else if(lex_val[lhs].size() != lex_val[rhs].size()) {
//             return lex_val[lhs].size() < lex_val[rhs].size();
//         // Finally compare lexicographic values (smallest goes first)
//         // This will be inverted when actually constructing gradient
//         } else {
//             for(unsigned int i = 0; i < lex_val[lhs].size(); i++) {
//                 if(lex_val[lhs][i] != lex_val[rhs][i]) {
//                     return lex_val[lhs][i] < lex_val[rhs][i];
//                 }
//             }

//             // If cells have identical lexicographic orderings, 
//             // then sort by raw cell index
//             return (lhs < rhs);
//         } 
    };
    
    std::vector<int> cells(comp.ncells);
    std::iota(cells.begin(), cells.end(), 0);
    std::sort(cells.begin(), cells.end(), cmp);
    
    for(int i = 0; i < comp.ncells; i++) {
        int c = cells[i];
        int star = cell_to_star[c];
        
        filt.set_star(c, star);
        filt.set_filtration_order(c, i);
             
    }
            
    
}


Filtration reduce_filtration(Filtration &full_filt, CellComplex &full_comp, CellComplex &red_comp, bool dual) {
    
    // Reduce filtration to only cover specified complex
    int target_dim = dual ? red_comp.dim : 0;
    
    std::vector<int> pcells;
    std::vector<double> time;
    for(int i = 0; i < full_comp.ncells; i++) {
        // Full complex must have label -1 for all cells not to be included in reduced complex below the reduced complex dimension
        if(full_comp.get_dim(i) == target_dim &&  full_comp.get_label(i) != -1) {
            time.push_back(full_filt.get_time(i));
            pcells.push_back(full_comp.get_label(i));
        }
    }
    
    Filtration red_filt(red_comp, time, dual);

    std::sort(pcells.begin(), pcells.end(),
             [&full_filt] (const int &lhs, const int &rhs) { 
                 return full_filt.get_filtration_order(lhs) < full_filt.get_filtration_order(rhs);
             });
    
    for(int ti = 0; ti < red_filt.nprimary_cells; ti++) {
        red_filt.set_primary_order(pcells[ti], ti);
    }
    
    
    construct_filtration(red_filt, red_comp, dual);
    
    return red_filt;
    
    // 1. Use full filtration ordering to get primary order of edges
    // 2. How to choose ordering of vertices?
        // can either redo vertex ordering by constructing a new filtration, or preserve ordering from full filtration
        // Might want to reorder
    // I don't think the ordering within the upper costar matters as long as lower dimension comes before higher dimension
    // In fact, the paper says it is so, so might as well run construct_filtration again to get the proper star cells
    

    
}
    
#endif // FILTRATION