#ifndef PERSIST
#define PERSIST
    
#include <pybind11/pybind11.h>
namespace py = pybind11;
    
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <utility>
    
#include "cell_complex.hpp"
#include "filtration.hpp"
 
// Note: This algorithm does not work for a discrete Morse complex
// Maxima and minima must both be the same type of cell for this to work (aka normal cell complex or Morse-Smale complex)
std::tuple<std::vector<std::pair<int, int> >, std::vector<std::pair<int, int> >, std::vector<std::pair<int, int> > >
    calc_extended_persistence(Filtration &filt_asc, Filtration &filt_desc, CellComplex &comp) {
    
    std::vector<std::vector<int> > columns(2*comp.ncells);
    
    std::vector<int> cell_to_col_asc(comp.ncells);
    std::vector<int> cell_to_col_desc(comp.ncells);
    
    std::vector<int> col_to_cell(2*comp.ncells);
    
    std::vector<int> filtration(comp.ncells);
    std::iota(filtration.begin(), filtration.end(), 0);
        
 
    std::sort(filtration.begin(), filtration.end(), 
         [&filt_asc](const int &lhs, const int &rhs) {return filt_asc.get_total_order(lhs) < filt_asc.get_total_order(rhs);});
        
    int icol = 0;
    for(auto ci: filtration) {
                        
        if(comp.regular) {
            auto facet_range = comp.get_facet_range(ci);
            for(auto it = facet_range.first; it != facet_range.second; it++) {
                columns[icol].push_back(cell_to_col_asc[*it]);
            }
        } else {
            auto facet_range = comp.get_facet_range(ci);
            auto coeff_range = comp.get_coeff_range(ci);
            for(auto itf = facet_range.first, itc = coeff_range.first; itf != facet_range.second; itf++, itc++) {
                if((*itc) % 2!= 0) {
                    columns[icol].push_back(cell_to_col_asc[*itf]);
                }
            }
        }
        
        columns[icol].shrink_to_fit();
        std::sort(columns[icol].begin(), columns[icol].end());
        cell_to_col_asc[ci] = icol;
        col_to_cell[icol] = ci;
        icol++;
    }
        
    std::sort(filtration.begin(), filtration.end(), 
         [&filt_desc](const int &lhs, const int &rhs) {return filt_desc.get_total_order(lhs) < filt_desc.get_total_order(rhs);});
        
    for(auto ci: filtration) {
                
        columns[icol].push_back(cell_to_col_asc[ci]);
        
        if(comp.regular) {
            auto facet_range = comp.get_facet_range(ci);
            for(auto it = facet_range.first; it != facet_range.second; it++) {
                columns[icol].push_back(cell_to_col_desc[*it]);
            }
        } else {
            auto facet_range = comp.get_facet_range(ci);
            auto coeff_range = comp.get_coeff_range(ci);
            for(auto itf = facet_range.first, itc = coeff_range.first; itf != facet_range.second; itf++, itc++) {
                if((*itc) % 2 != 0) {
                    columns[icol].push_back(cell_to_col_desc[*itf]);
                }
            }
        }
        
        columns[icol].shrink_to_fit();
        std::sort(columns[icol].begin(), columns[icol].end());
        cell_to_col_desc[ci] = icol;
        col_to_cell[icol] = ci;
        icol++;
    }
        
    // py::print(columns);
    
    // pivot row of each column if it has one
    std::unordered_map<int, int> pivot_row;
    // row to reduced column with pivot in that row
    std::unordered_map<int, int> pivot_col;
    
    for(unsigned int j = 0; j < columns.size(); j++) {
         
        if(columns[j].size()) {
            pivot_row[j] = *std::max_element(columns[j].begin(), columns[j].end());
        }
        
        while(columns[j].size() && pivot_col.count(pivot_row[j])) {

            int l = pivot_col[pivot_row[j]];
            
            // py::print(j, columns[j]);
            // py::print(l, columns[l]);
            
            // Symmetric difference
            // columns[j] ^= columns[l]
            std::vector<int> col_sum(columns[j].size() + columns[l].size());
            
            auto it = std::set_symmetric_difference (columns[j].begin(), columns[j].end(), 
                                                     columns[l].begin(), columns[l].end(), 
                                                     col_sum.begin());
            col_sum.resize(it-col_sum.begin());
            
            columns[j].assign(col_sum.begin(), col_sum.end());
            
            // py::print(j, "+", l, col_sum);
            
            if(columns[j].size()) {
                pivot_row[j] = *std::max_element(columns[j].begin(), columns[j].end());
            } else {
                pivot_row.erase(j);
            }
                        
        }

        if(columns[j].size()) {
            pivot_col[pivot_row[j]] = j;
        }
                
        columns[j].shrink_to_fit();
                            
    }
        
    std::vector<std::pair<int, int> > ord_pairs; 
    std::vector<std::pair<int, int> > rel_pairs;
    std::vector<std::pair<int, int> > ext_pairs;
    for(auto pair: pivot_row) {
        int i = pair.second;
        int j = pair.first;
        int ci = col_to_cell[i];
        int cj = col_to_cell[j];
        
        if(j < comp.ncells) {
            ord_pairs.emplace_back(ci, cj);
        } else if(i < comp.ncells) {
            ext_pairs.emplace_back(ci, cj);
        } else {
            rel_pairs.emplace_back(ci, cj);
        }
        
        
    }
    
    return std::forward_as_tuple(ord_pairs, rel_pairs, ext_pairs);
    
    
}
    
#endif // PERSIST