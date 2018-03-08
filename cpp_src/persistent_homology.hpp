#ifndef PERSIST_HPP
#define PERSIST_HPP
    
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
    
    
std::tuple<std::vector<std::vector<int> >, std::vector<int>, std::vector<int> > 
    calc_boundary_mat(Filtration &filt, CellComplex &comp) {
    
    std::vector<std::vector<int> > columns(comp.ncells);
    std::vector<int> cell_to_col(comp.ncells);    
    std::vector<int> col_to_cell(comp.ncells);
    
    std::vector<int> cell_order = filt.get_filtration();
    for(std::size_t i = 0; i < cell_order.size(); i++) {
        
        int ci = cell_order[i];

        if(comp.regular) {
            auto facet_range = comp.get_facet_range(ci);
            for(auto it = facet_range.first; it != facet_range.second; it++) {
                columns[i].push_back(cell_to_col[*it]);
            }
        } else {
            auto facet_range = comp.get_facet_range(ci);
            auto coeff_range = comp.get_coeff_range(ci);
            for(auto itf = facet_range.first, itc = coeff_range.first; itf != facet_range.second; itf++, itc++) {
                if((*itc) % 2!= 0) {
                    columns[i].push_back(cell_to_col[*itf]);
                }
            }
        }
        
        columns[i].shrink_to_fit();
        std::sort(columns[i].begin(), columns[i].end());
        cell_to_col[ci] = i;
        col_to_cell[i] = ci;
        
    }
    
    return std::forward_as_tuple(columns, cell_to_col, col_to_cell);
    
}


std::vector<int> add_cols_Z2(std::vector<int> &col1, std::vector<int> &col2) {
    
    std::vector<int> col_sum(col1.size() + col2.size());
    auto it = std::set_symmetric_difference (col1.begin(), col1.end(), 
                                             col2.begin(), col2.end(), 
                                             col_sum.begin());
    col_sum.resize(it-col_sum.begin());
    
    return col_sum;
    
}


std::vector<std::vector<int> > reduce_smith_normal_form(std::vector<std::vector<int> > &columns, bool birth_cycles=false) {
    
    // pivot row of each column if it has one
    std::unordered_map<int, int> pivot_row;
    // row to reduced column with pivot in that row
    std::unordered_map<int, int> pivot_col;
     
    std::vector<std::vector<int> > g;
    
    for(unsigned int j = 0; j < columns.size(); j++) {
        
        if(birth_cycles) {
            g.emplace_back(1, j);
        }
        
        while(columns[j].size()) {
            
            int pivot_row = columns[j].back();
            
            if(!pivot_col.count(pivot_row)) {
                break;
            }

            int l = pivot_col[pivot_row];
            
            // py::print(j, columns[j]);
            // py::print(l, columns[l]);
            
            std::vector<int> col_sum = add_cols_Z2(columns[j], columns[l]);
            columns[j].assign(col_sum.begin(), col_sum.end());
            
            // py::print(j, "+", l, col_sum);
            
            if(birth_cycles) {
                std::vector<int> g_sum = add_cols_Z2(g[j], g[l]);
                g[j].assign(g_sum.begin(), g_sum.end());
            }
       
        }

        if(columns[j].size()) {
            pivot_col[columns[j].back()] = j;
        }
                
        columns[j].shrink_to_fit();
                            
    }
        
    return g;
    
}

// // Note: This algorithm does not work for a discrete Morse complex
// // Maxima and minima must both be the same type of cell for this to work (e.g. standard cell complex or Morse-Smale complex)
std::tuple<std::vector<std::pair<int, int> >, std::vector<std::pair<int, int> >, std::vector<std::pair<int, int> > >
    calc_extended_persistence(Filtration &filt_asc, Filtration &filt_desc, CellComplex &comp) {

    // Get boundary matrices for ascending and descending filtrations
    auto bound_mat_asc = calc_boundary_mat(filt_asc, comp);
    auto bound_mat_desc = calc_boundary_mat(filt_desc, comp);

    // Initialize columns and column maps with ascending filtration
    std::vector<std::vector<int> > columns = std::get<0>(bound_mat_asc);
    std::vector<int> col_to_cell = std::get<2>(bound_mat_asc);

    
    // Create extended boundary matrix
    columns.resize(2*comp.ncells);
    col_to_cell.resize(2*comp.ncells);
    for(int i = 0; i < comp.ncells; i++) {
        
        // Get cell in descending filtration
        int ci = std::get<2>(bound_mat_desc)[i];
        // Get row of cell in ascending ascending boundary matrix
        int irow = std::get<1>(bound_mat_asc)[ci];
        // Add cell to permutation matrix
        columns[i+comp.ncells].push_back(irow);
        
        // Add column to extended boundary matrix at an offset
        std::vector<int> col = std::get<0>(bound_mat_desc)[i];
        for(auto cj: col) {
            columns[i+comp.ncells].push_back(cj+comp.ncells);
        }
        columns[i+comp.ncells].shrink_to_fit();
        
        col_to_cell[i+comp.ncells] = ci;
        
    }
    
    py::print(columns);
        
    reduce_smith_normal_form(columns);
    
    py::print(columns);    
    
    std::vector<std::pair<int, int> > ord_pairs; 
    std::vector<std::pair<int, int> > rel_pairs;
    std::vector<std::pair<int, int> > ext_pairs;
    for(std::size_t j = 0; j < columns.size(); j++) {
        if(!columns[j].size()) {
            continue;
        }
        
        int i = columns[j].back();
        int ci = col_to_cell[i];
        int cj = col_to_cell[j];
        
        if(j < (std::size_t)comp.ncells) {
            ord_pairs.emplace_back(ci, cj);
        } else if(i < comp.ncells) {
            ext_pairs.emplace_back(ci, cj);
        } else {
            rel_pairs.emplace_back(ci, cj);
        }
        
        
    }
    
    return std::forward_as_tuple(ord_pairs, rel_pairs, ext_pairs);
    
}

std::unordered_map<int, std::vector<int> > calc_birth_cycles(Filtration &filt, CellComplex &comp) {
    
    auto bound_mat = calc_boundary_mat(filt, comp);
    
    std::vector<std::vector<int> > columns = std::get<0>(bound_mat);
    std::vector<int> col_to_cell = std::get<2>(bound_mat);
    
    std::vector<std::vector<int> > g = reduce_smith_normal_form(columns, true);
    
    py::print(g);
    
    
    std::unordered_map<int, std::vector<int> > cycles;
    
    for(std::size_t j = 0; j < columns.size(); j++) {
        if(columns[j].size()) {
            continue;
        }
        
        int cj = col_to_cell[j];
        cycles[cj];
        for(auto gi: g[j]) {
            cycles[cj].push_back(col_to_cell[gi]);
        }
    }
    
    return cycles;
}


// std::tuple<std::vector<std::pair<int, int> >, std::vector<std::pair<int, int> >, std::vector<std::pair<int, int> > >
//     calc_extended_persistence(Filtration &filt_asc, Filtration &filt_desc, CellComplex &comp) {
    
//     std::vector<std::vector<int> > columns(2*comp.ncells);
    
//     std::vector<int> cell_to_col_asc(comp.ncells);
//     std::vector<int> cell_to_col_desc(comp.ncells);
    
//     std::vector<int> col_to_cell(2*comp.ncells);
            
//     int icol = 0;
//     for(auto ci: filt_asc.get_filtration()) {
                        
//         if(comp.regular) {
//             auto facet_range = comp.get_facet_range(ci);
//             for(auto it = facet_range.first; it != facet_range.second; it++) {
//                 columns[icol].push_back(cell_to_col_asc[*it]);
//             }
//         } else {
//             auto facet_range = comp.get_facet_range(ci);
//             auto coeff_range = comp.get_coeff_range(ci);
//             for(auto itf = facet_range.first, itc = coeff_range.first; itf != facet_range.second; itf++, itc++) {
//                 if((*itc) % 2!= 0) {
//                     columns[icol].push_back(cell_to_col_asc[*itf]);
//                 }
//             }
//         }
        
//         columns[icol].shrink_to_fit();
//         std::sort(columns[icol].begin(), columns[icol].end());
//         cell_to_col_asc[ci] = icol;
//         col_to_cell[icol] = ci;
//         icol++;
//     }
          
//     for(auto ci: filt_desc.get_filtration()) {
                
//         columns[icol].push_back(cell_to_col_asc[ci]);
        
//         if(comp.regular) {
//             auto facet_range = comp.get_facet_range(ci);
//             for(auto it = facet_range.first; it != facet_range.second; it++) {
//                 columns[icol].push_back(cell_to_col_desc[*it]);
//             }
//         } else {
//             auto facet_range = comp.get_facet_range(ci);
//             auto coeff_range = comp.get_coeff_range(ci);
//             for(auto itf = facet_range.first, itc = coeff_range.first; itf != facet_range.second; itf++, itc++) {
//                 if((*itc) % 2 != 0) {
//                     columns[icol].push_back(cell_to_col_desc[*itf]);
//                 }
//             }
//         }
        
//         columns[icol].shrink_to_fit();
//         std::sort(columns[icol].begin(), columns[icol].end());
//         cell_to_col_desc[ci] = icol;
//         col_to_cell[icol] = ci;
//         icol++;
//     }
        
//     // py::print(columns);
    
//     // pivot row of each column if it has one
//     std::unordered_map<int, int> pivot_row;
//     // row to reduced column with pivot in that row
//     std::unordered_map<int, int> pivot_col;
    
//     for(unsigned int j = 0; j < columns.size(); j++) {
         
//         if(columns[j].size()) {
//             pivot_row[j] = *std::max_element(columns[j].begin(), columns[j].end());
//         }
        
//         while(columns[j].size() && pivot_col.count(pivot_row[j])) {

//             int l = pivot_col[pivot_row[j]];
            
// //             py::print(j, columns[j]);
// //             py::print(l, columns[l]);
            
//             // Symmetric difference
//             std::vector<int> col_sum(columns[j].size() + columns[l].size());
//             auto it = std::set_symmetric_difference (columns[j].begin(), columns[j].end(), 
//                                                      columns[l].begin(), columns[l].end(), 
//                                                      col_sum.begin());
//             col_sum.resize(it-col_sum.begin());
            
//             columns[j].assign(col_sum.begin(), col_sum.end());
            
//             // py::print(j, "+", l, col_sum);
            
//             if(columns[j].size()) {
//                 pivot_row[j] = *std::max_element(columns[j].begin(), columns[j].end());
//             } else {
//                 pivot_row.erase(j);
//             }
                        
//         }

//         if(columns[j].size()) {
//             pivot_col[pivot_row[j]] = j;
//         }
                
//         columns[j].shrink_to_fit();
                            
//     }
    
    
//     // py::print(columns);
        
//     std::vector<std::pair<int, int> > ord_pairs; 
//     std::vector<std::pair<int, int> > rel_pairs;
//     std::vector<std::pair<int, int> > ext_pairs;
//     for(auto pair: pivot_row) {
//         int i = pair.second;
//         int j = pair.first;
//         int ci = col_to_cell[i];
//         int cj = col_to_cell[j];
        
//         if(j < comp.ncells) {
//             ord_pairs.emplace_back(ci, cj);
//         } else if(i < comp.ncells) {
//             ext_pairs.emplace_back(ci, cj);
//         } else {
//             rel_pairs.emplace_back(ci, cj);
//         }
        
        
//     }
    
//     return std::forward_as_tuple(ord_pairs, rel_pairs, ext_pairs);
    
    
// }


    
#endif // PERSIST_HPP