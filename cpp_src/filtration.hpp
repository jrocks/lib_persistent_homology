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
    
    // Number of cells
    const int ncells;
    // Dimension of cells defining filtration
    const int fdim;
    // Number of steps in filtration (number of cells of dimension fdim)
    const int nsteps;
    // Ascending or descending filtration?
    const bool ascend;
    // star or costar filtration?
    const bool co;
    
private:
 
    // Cell which defines each subcomplex
    // Cells on which time is explicitly defined
    std::vector<int> filt_cell;
    
    // Filtration time of subcomplexes
    std::vector<double> time;
    // Ordering of subcomplexes
    std::vector<int> subcomplex_order;
       
    // Which subcomplex each cells first appears
    std::vector<int> subcomplex;
    // Total filtration orderering of cells
    std::vector<int> total_order;
    
    
public:
    
    Filtration(std::vector<double> &time, std::vector<int> &subcomplex_order, CellComplex &comp, bool ascend = true, bool co = false) : 
        ncells(comp.ncells), fdim(co ? comp.dim : 0), nsteps(time.size()), ascend(ascend), co(co),
        filt_cell(time.size()), time(time), subcomplex_order(subcomplex_order), 
        subcomplex(comp.ncells), total_order(comp.ncells) {
            
            for(int i = 0, index = 0; i < ncells; i++) {
                if(comp.get_dim(i) == fdim) {
                    subcomplex[i] = index;
                    filt_cell[index] = i;
                    index++;
                }
            }
                        
        }
    
    // Get filtration cell defining subcomplex which contains cell alpha
    int get_filt_cell(int alpha) {
        return filt_cell[subcomplex[alpha]];
    }
    
    // Get insertion time of the subcomplex of cell alpha
    double get_time(int alpha) {
        return time[subcomplex[alpha]];
    }
    
    // Get insertion order of subcomplex containing alpha
    int get_subcomplex_order(int alpha) {
        return subcomplex_order[subcomplex[alpha]];
    }
    
    // Add cell alpha to the subcomplex containing cell beta
    void add_to_subcomplex(int alpha, int beta) {
        subcomplex[alpha] = subcomplex[beta];
    }
    
    // Set the total ordering of cell alpha
    void set_total_order(int alpha, int order) {
        total_order[alpha] = order;
    }
    
    // Get the total ordering of cell alpha
    int get_total_order(int alpha) {
        return total_order[alpha];
    }
    
    std::vector<int> get_filtration() {
        std::vector<int> filtration(ncells);
        std::iota(filtration.begin(), filtration.end(), 0);
        std::sort(filtration.begin(), filtration.end(), 
                 [this](const int &lhs, const int &rhs) {return this->total_order[lhs] < this->total_order[rhs];});
        
        return filtration;
    }
    
};


std::vector<int> perform_watershed_transform(std::vector<double> &time, CellComplex &comp, bool ascend = true, bool co = false) {    
    
    int fdim = co ? comp.dim : 0;
    
    std::vector<int> subcomplex_order(time.size());
    
    std::unordered_map<int, int> cell_to_index;
    std::vector<bool> submerged(time.size(), false);
            
    // Find primary cells and sort according to insertion time
    std::vector<int> filt_cell_argsort;
    for(int i = 0, index = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) == fdim) {
            filt_cell_argsort.push_back(i);
            cell_to_index[i] = index;
            index++;
        }
    }
    std::sort(filt_cell_argsort.begin(), filt_cell_argsort.end(),
       [&time, &cell_to_index, ascend](int lhs, int rhs) {
           if(ascend) {
               // <=?
               return time[cell_to_index[lhs]] < time[cell_to_index[rhs]];
           } else {
               // >=?
               return time[cell_to_index[lhs]] > time[cell_to_index[rhs]];
           }
       });
                
    // Iterate through each level set
    unsigned int ti = 0;
    while(ti < filt_cell_argsort.size()) {
        
        std::unordered_set<int> level;
        
        unsigned int tj = ti;
        while(true) {
            int ci = filt_cell_argsort[tj];
            double t = time[cell_to_index[ci]];
            level.insert(ci);
            
            if((tj+1 == filt_cell_argsort.size()) || (t != time[cell_to_index[filt_cell_argsort[tj+1]]])) {
                break;
            }
            
            tj++;
        }
                
        if(level.size() == 1) {
            int ci = *(level.begin());
            submerged[cell_to_index[ci]] = true;
            subcomplex_order[cell_to_index[ci]] = ti;
            ti++;
            
            continue;
        }
        
        std::unordered_map<int, double> dist;
        std::queue<int> Q;
        for(auto a: level) {
            
            auto rangea = co ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto ita = rangea.first; ita != rangea.second; ita++) {
                
                int b = *ita;
                
                auto rangeb = co ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
                    
                    int c = *itb;
                    
                    if(submerged[cell_to_index[c]]) {
                        
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
            
            if(!level.count(a)) {
                continue;
            }
            
            
            level.erase(a);
            submerged[cell_to_index[a]] = true;
            subcomplex_order[cell_to_index[a]] = ti;
            ti++;
            
            auto rangea = co ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto ita = rangea.first; ita != rangea.second; ita++) {
                
                int b = *ita;
                
                auto rangeb = co ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
                    
                    int c = *itb;
                    
                    if(level.count(c)) {
                        
                        double new_dist = current_dist + 1.0;
                        
                        if(!dist.count(c) || new_dist < dist[c]) {
                            dist[c] = new_dist;
                            Q.push(c);
                        }
                    }
                    
                }
                
            }
            
            
            
        }
        
        
        while(!level.empty()) {
            int s = *(level.begin());
            level.erase(s);
            
            Q.push(s);
            
            while(!Q.empty()) {
                int a = Q.front();
                Q.pop();

                // This differs from the python version;
                double current_dist = dist[a];

                level.erase(a);
                submerged[cell_to_index[a]] = true;
                subcomplex_order[cell_to_index[a]] = ti;
                ti++;



                auto rangea = co ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
                for(auto ita = rangea.first; ita != rangea.second; ita++) {

                    int b = *ita;

                    std::vector<int>::iterator beginb;
                    std::vector<int>::iterator endb;

                    auto rangeb = co ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                    for(auto itb = rangeb.first; itb != rangeb.second; itb++) {

                        int c = *itb;

                        if(level.count(c)) {
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
    
    
    return subcomplex_order;
    
}





Filtration construct_filtration(std::vector<double> &time, std::vector<int> &subcomplex_order, CellComplex &comp, bool ascend = true , bool co = false) {
    
    Filtration filt(time, subcomplex_order, comp, ascend, co);
    
    auto star_cmp = [&filt, ascend, co](const int &lhs, const int &rhs) {
        if((ascend && co) || (!ascend && !co)) {
            return filt.get_subcomplex_order(lhs) < filt.get_subcomplex_order(rhs); 
        } else {
            return filt.get_subcomplex_order(lhs) > filt.get_subcomplex_order(rhs); 
        }
    };
    
    
    // Calculate lexicographic value of each cell and assign to subcomplex
    std::vector<std::vector<int> > lex_val(comp.ncells);
    for(int i = 0; i < comp.ncells; i++) {
                            
        std::unordered_set<int> star = get_star(i, !co, comp, filt.fdim);
        
        for(auto alpha: star) {
            lex_val[i].push_back(filt.get_subcomplex_order(alpha));
        }
        
        int filt_cell = *std::min_element(star.begin(), star.end(), star_cmp);
        filt.add_to_subcomplex(i, filt_cell);
                        
        // If ascend filtration sort from high to low
        // If upper filtration sort from low to high
        if(ascend) {
            std::sort(lex_val[i].begin(), lex_val[i].end(), std::greater<int>()); 
        } else {
            std::sort(lex_val[i].begin(), lex_val[i].end(), std::less<int>()); 
        }
            
    }    
    
    
        
    auto lex_cmp = [&comp, &lex_val, ascend, co](const int &lhs, const int &rhs) {
        
        // Ascending filtraiton
        if(ascend) {
            // Upper costar filtration
            // Compare upper costar subcomplex value (smallest goes first) 
            if(co && lex_val[lhs].back() != lex_val[rhs].back()) {
                return lex_val[lhs].back() < lex_val[rhs].back();
            // Lower star filtration
            // Compare lower star subcomplex value (smallest goes first) 
            } else if(!co && lex_val[lhs].front() != lex_val[rhs].front()) {
                return lex_val[lhs].front() < lex_val[rhs].front();
            // Compare dimensions (smallest goes first)
            } else if(comp.get_dim(lhs) != comp.get_dim(rhs)) {
                return comp.get_dim(lhs) < comp.get_dim(rhs);
            // Compare lexicographic values (smallest goes first)
            } else if(lex_val[lhs] != lex_val[rhs]) {
                return lex_val[lhs] < lex_val[rhs];
            // Finally, if cells have identical lexicographic orderings, 
            // then sort by raw cell index
            } else {
                // py::print("Breaking tie with indices", lhs, rhs);
                return lhs < rhs;
            }
        } else {
            // Lower costar filtration
            // Compare lower costar subcomplex value (largest goes first) 
            if(co && lex_val[lhs].back() != lex_val[rhs].back()) {
                return lex_val[lhs].back() > lex_val[rhs].back();
            // Upper star filtration
            // Compare upper star subcomplex value (largest goes first) 
            } else if(!co && lex_val[lhs].front() != lex_val[rhs].front()) {
                return lex_val[lhs].front() > lex_val[rhs].front();
            // Compare dimensions (smallest goes first)
            } else if(comp.get_dim(lhs) != comp.get_dim(rhs)) {
                return comp.get_dim(lhs) < comp.get_dim(rhs);
            // Compare lexicographic values (largest goes first)
            } else if(lex_val[lhs] != lex_val[rhs]) {
                return lex_val[lhs] > lex_val[rhs];
            // Finally, if cells have identical lexicographic orderings, 
            // then sort by raw cell index
            } else {
                // py::print("Breaking tie with indices", lhs, rhs);
                return lhs > rhs;
            }
        }
            
    };
    
    std::vector<int> cells(comp.ncells);
    std::iota(cells.begin(), cells.end(), 0);
    std::sort(cells.begin(), cells.end(), lex_cmp);
    
    for(int i = 0; i < comp.ncells; i++) {
        filt.set_total_order(cells[i], i);
    }
            
    
    return filt;
    
    
}


Filtration reduce_filtration(Filtration &full_filt, CellComplex &full_comp, CellComplex &red_comp) {
    
    // This other method mixes subcomplexes up because it relabels them
    // Edges are still added in correct order
    
//     // Reduce filtration to only cover specified complex
//     int red_fdim = full_filt.co ? red_comp.dim : 0;
    
//     // Get map of filtration cells to index in reduced filtration
//     std::unordered_map<int, int> filt_cells_to_index;
//     std::vector<int> filt_cells;
//     for(int i = 0; i < full_comp.ncells; i++) {
//         // Full complex must have label -1 for all cells not to be included in reduced complex below the reduced complex dimension
//         if(full_comp.get_dim(i) == red_fdim && full_comp.get_label(i) != -1) {
//             filt_cells_to_index[i] = filt_cells_to_index.size();
//             filt_cells.push_back(i);
//         }
//     }
    
//     auto cmp = [&full_filt] (const int &lhs, const int &rhs) { 
//                   return full_filt.get_total_order(lhs) < full_filt.get_total_order(rhs);
//               };
    
    
//     // Sort filtration cells
//     std::sort(filt_cells.begin(), filt_cells.end(), cmp);
    
//     // Get new subcomplex order and time
//     std::vector<double> time(filt_cells.size());
//     std::vector<int> subcomplex_order(filt_cells.size());
//     for(unsigned int i = 0; i < filt_cells.size(); i++) {
//         time[filt_cells_to_index[filt_cells[i]]] = full_filt.get_time(filt_cells[i]);
//         subcomplex_order[filt_cells_to_index[filt_cells[i]]] = i;
         
//     }
    
    
    // This second method also mixes subcomplexes, but can just pretend that the corners represent the true subcomplexes
    // Hence the fact that the raw original subcomplex order numbers are used
    // However, each vertex is still mapped to a specific edge
    
    // Reduce filtration to only cover specified complex
    int red_fdim = full_filt.co ? red_comp.dim : 0;
    
    // Get subcomplex order and time    
    std::vector<double> time;
    std::vector<int> subcomplex_order;
    for(int i = 0; i < full_comp.ncells; i++) {
        // Full complex must have label -1 for all cells not to be included in reduced complex below the reduced complex dimension
        if(full_comp.get_dim(i) == red_fdim && full_comp.get_label(i) != -1) {
            time.push_back(full_filt.get_time(i));
            subcomplex_order.push_back(full_filt.get_subcomplex_order(i));
        }
    }
    
    Filtration red_filt(time, subcomplex_order, red_comp, full_filt.ascend, full_filt.co);
    
    
    std::vector<int> red_cells;
    for(int i = 0; i < full_comp.ncells; i++) {
        if(full_comp.get_dim(i) <= red_comp.dim && full_comp.get_label(i) != -1) {
            red_cells.push_back(i);
        }
    }
    
    auto cmp = [&full_filt] (const int &lhs, const int &rhs) { 
                  return full_filt.get_total_order(lhs) < full_filt.get_total_order(rhs);
              };
    
    std::sort(red_cells.begin(), red_cells.end(), cmp);
    
    
    for(int i = 0; i < red_comp.ncells; i++) {
        int ci = red_cells[i];
        
        red_filt.set_total_order(full_comp.get_label(ci), i);
        
        std::unordered_set<int> star = get_star(i, !red_filt.co, red_comp, red_filt.fdim);
        int filt_cell = *std::min_element(star.begin(), star.end(), cmp);
        red_filt.add_to_subcomplex(i, filt_cell);
        
    }
                         
    return red_filt;
    
}
    
#endif // FILTRATION