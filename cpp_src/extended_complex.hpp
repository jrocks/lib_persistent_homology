#ifndef EXTCOMPLEX_HPP
#define EXTCOMPLEX_HPP

#include "eigen_macros.hpp"
#include "cell_complex.hpp"
#include "filtration.hpp"
#include "morse_complex.hpp"

#include <vector>

#include <pybind11/pybind11.h>    
namespace py = pybind11;


CellComplex extend_complex(CellComplex &comp) {
    
    
    CellComplex ext_comp(comp.dim+1, comp.regular, comp.oriented);

    // Sort cells accoring to dimension
    std::vector<std::vector<int> > cells(comp.dim+1);
    for(int c = 0; c < comp.ncells; c++) {
        cells[comp.get_dim(c)].push_back(c);
    }
    
    // First add dummy vertex
    std::vector<int> facets;
    std::vector<int> coeffs;
    ext_comp.add_cell(-1, 0, facets, coeffs);
    
    // Map of cells in original complex to same cells in extended complex
    std::vector<int> orig_map(comp.ncells);
    // Map of cells in original complex to analogous cells of one dimension higher in extended complex
    std::vector<int> ext_map(comp.ncells);
    
    // Copy cells from original complex of dimension d
    for(int d = 0; d <= comp.dim; d++) {
        
        for(auto c: cells[d]) {
        
            orig_map[c] = ext_comp.ncells;

            std::vector<int> facets;
            auto facet_range = comp.get_facet_range(c);
            for(auto it = facet_range.first; it != facet_range.second; it++) {
                facets.push_back(orig_map[*it]);
            }
            
            auto coeff_range = comp.get_coeff_range(c);
            std::vector<int> coeffs(coeff_range.first, coeff_range.second);
            ext_comp.add_cell(c, d, facets, coeffs);
            
        }
    }
    
    // Add cell of dimension d+1 for each cell of dimension d in original complex
    for(int d = 0; d <= comp.dim; d++) {
        
        for(auto c: cells[d]) {
            
            ext_map[c] = ext_comp.ncells;
                        
            std::vector<int> facets;
            // Extended cell will have analog of itself from lower dimension as face
            facets.push_back(orig_map[c]);
            // If cell was originally a vertex, then add dummy vertex as face
            if(d == 0) {
                facets.push_back(0);
            // Otherwise add analogous faces of original cell, but one dimension higher
            } else {
                auto facet_range = comp.get_facet_range(c);
                for(auto it = facet_range.first; it != facet_range.second; it++) {
                    facets.push_back(ext_map[*it]);
                }   
            }
            
            auto coeff_range = comp.get_coeff_range(c);
            std::vector<int> coeffs(coeff_range.first, coeff_range.second);
            
            ext_comp.add_cell(c, d+1, facets, coeffs);
            
        }
        
    }
    
    ext_comp.construct_cofacets();
    ext_comp.make_compressed(); 
    
    return ext_comp;
    
}


Filtration extend_filtration(Filtration &filt_asc, Filtration &filt_desc, CellComplex &ext_comp) {
    
    
    // Number of original cells
    int N = (ext_comp.ncells-1)/2;
    
    std::vector<double> func(ext_comp.ncells);  
    std::vector<int> digi_func(ext_comp.ncells);
    std::vector<int> order(ext_comp.ncells);
    
    // First add dummy vertex
    func[0] = 0.0; // This value doesn't matter because this component is born first and never dies
    digi_func[0] = -1; // digi_func typically starts at zero, so dummy vertex doesn't have same value as any other vertex
    order[0] = -1; // same with order
    
    // Process ascending filtration
    for(int c = 1; c < 1 + N; c++) {
        func[c] = filt_asc.get_func(ext_comp.get_label(c));
        digi_func[c] = filt_asc.get_digi_func(ext_comp.get_label(c));
        order[c] = filt_asc.get_digi_func(ext_comp.get_label(c));
    }
    
    // Process descending filtration
    for(int c = 1 + N; c < ext_comp.ncells; c++) {
        func[c] = filt_desc.get_func(ext_comp.get_label(c));
        digi_func[c] = filt_desc.get_digi_func(ext_comp.get_label(c))+N;
        order[c] = filt_desc.get_digi_func(ext_comp.get_label(c))+N;
    }
    
    return Filtration(ext_comp, func, digi_func, order, true);
    
    
}

std::tuple<XiVec, XiVec> extend_discrete_gradient(RXiVec V_asc, RXiVec V_desc, CellComplex &ext_comp) {
    
     // Number of original cells
    int N = (ext_comp.ncells-1)/2;
    
    // Map of cells in original complex to same cells in extended complex
    std::vector<int> orig_map(N);
    for(int c = 1; c < 1 + N; c++) {
        orig_map[ext_comp.get_label(c)] = c;
    }
    
    // Map of cells in original complex to analogous cells of one dimension higher in extended complex
    std::vector<int> ext_map(N);
    for(int c = 1 + N; c < ext_comp.ncells; c++) {
        ext_map[ext_comp.get_label(c)] = c;
    }
    
    XiVec V = XiVec::Constant(ext_comp.ncells, -1);
    XiVec coV = XiVec::Constant(ext_comp.ncells, -1);
    
    // Dummy vertex is critical cell
    V[0] = 0;
    coV[0] = 0;
    
    // Process ascending gradient
    for(int c = 1; c < 1 + N; c++) {
        V[c] = orig_map[V_asc[ext_comp.get_label(c)]];
        coV[V[c]] = c;
    }
    
    // Process descending filtration
    for(int c = 1 + N; c < ext_comp.ncells; c++) {
        V[c] = ext_map[V_desc[ext_comp.get_label(c)]];
        coV[V[c]] = c;
    }
    
    return std::make_tuple(V, coV);
    
}



#endif // EXTCOMPLEX_HPP