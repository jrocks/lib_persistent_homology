#ifndef CUBICALCOMPLEX
#define CUBICALCOMPLEX
 
    
#include <vector>
#include "cell_complex.hpp"
#include "filtration.hpp"
    

CellComplex construct_cubical_complex(std::vector<int> &shape, bool oriented, bool dual) {
    
    int dim = shape.size();
    
    CellComplex comp(dim, true, oriented);
    
    // cubical complexes are labeled according to corresponding pixel
    
    if(dim == 2) {
        int nrows = dual ? shape[0] + 1: shape[0];
        int ncols = dual ? shape[1] + 1: shape[1];
        
        int nverts = nrows*ncols;
        int nhedges = nrows*(ncols-1);

        
        // vertices
        for(int i = 0; i < nrows*ncols; i++) {
            std::vector<int> facets;
            std::vector<int> coeffs;
            comp.add_cell(i, 0, facets, coeffs);
        }
        
        // horizontal edges
        for(int i = 0; i < nrows; i++) {
            for(int j = 0; j < ncols-1; j++) {
                std::vector<int> facets;                
                facets.push_back(ncols*i + j);
                facets.push_back(ncols*i + j+1);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell((ncols-1)*i + j, 1, facets, coeffs);
                
            }
        }
        
        // vertical edges
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols; j++) {
                std::vector<int> facets;                
                facets.push_back(ncols*i + j);
                facets.push_back(ncols*(i+1) + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(nhedges + ncols*i + j, 1, facets, coeffs);
                
            }
        }
        
        // faces
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols-1; j++) {
                std::vector<int> facets;                
                facets.push_back(nverts + (ncols-1)*i + j);
                facets.push_back(nverts + nhedges + ncols*i + j+1);
                facets.push_back(nverts + (ncols-1)*(i+1) + j);
                facets.push_back(nverts + nhedges + ncols*i + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(1);
                    coeffs.push_back(1);
                    coeffs.push_back(-1);
                    coeffs.push_back(-1);
                }
                
                comp.add_cell((ncols-1)*i + j, 2, facets, coeffs);
                
            }
        }
    }
    
    comp.construct_cofacets();
    comp.make_compressed(); 
        
    return comp;
    
}


CellComplex construct_masked_cubical_complex(std::vector<bool> &mask, std::vector<int> &shape, bool oriented, bool dual) {
      
    CellComplex complete_comp = construct_cubical_complex(shape, oriented, dual);
    
    std::vector<bool> include(complete_comp.ncells, false);
    
    for(int i = 0; i < complete_comp.ncells; i++) {
        
        if((dual && complete_comp.get_dim(i) == complete_comp.dim) 
          || (!dual && complete_comp.get_dim(i) == 0)) {
            
            include[i] = !mask[complete_comp.get_label(i)];

            
        }
    }
    
    if(dual) {
        for(int d = complete_comp.dim-1; d >= 0; d--) {
            for(int i = 0; i < complete_comp.ncells; i++) {
                if(complete_comp.get_dim(i) == d) {
                                        
                    auto range = complete_comp.get_cofacet_range(i);
                    for(auto it = range.first; it != range.second; it++) {
                        
                        if(include[*it]) {
                            include[i] = true;
                        }
                        
                    }
                }
            }
        }
    } else {
        for(int d = 1; d <= complete_comp.dim; d++) {
            for(int i = 0; i < complete_comp.ncells; i++) {
                if(complete_comp.get_dim(i) == d) {
                                        
                    auto range = complete_comp.get_facet_range(i);
                    for(auto it = range.first; it != range.second; it++) {
                        
                        if(include[*it]) {
                            mask[i] = true;
                        }
                        
                    }
                }
            }
        }
    }
      
    
    std::vector<int> complete_to_mask(complete_comp.ncells);
    int index = 0;
    for(int i = 0; i < complete_comp.ncells; i++) {
        if(include[i]) {
            complete_to_mask[i] = index;
            index++;
        }
    }
    
    CellComplex mask_comp(complete_comp.dim, true, oriented);
    
    for(int i = 0; i < complete_comp.ncells; i++) {
        if(include[i]) {
            
            std::vector<int> facets;
            auto facet_range = complete_comp.get_facet_range(i);
            for(auto it = facet_range.first; it != facet_range.second; it++) {
                facets.push_back(complete_to_mask[*it]);
            }
            std::vector<int> coeffs;
            auto coeff_range = complete_comp.get_coeff_range(i);
            coeffs.assign(coeff_range.first, coeff_range.second);
            
            mask_comp.add_cell(complete_comp.get_label(i), complete_comp.get_dim(i), facets, coeffs);
        }
    }
       
    mask_comp.construct_cofacets();
    mask_comp.make_compressed(); 
        
    return mask_comp;
    
}



// CellComplex construct_masked_cubical_complex(std::vector<bool> &mask, std::vector<int> &shape, bool oriented, bool dual) {
      
//     CellComplex complete_comp = construct_cubical_complex(shape, oriented, dual);
    
//     std::vector<bool> include(complete_comp.ncells, false);
    
//     for(int i = 0; i < complete_comp.ncells; i++) {
        
//         if((dual && complete_comp.get_dim(i) == complete_comp.dim) 
//           || (!dual && complete_comp.get_dim(i) == 0)) {
            
//             include[i] = !mask[complete_comp.get_label(i)];
//         }
//     }
    
//     if(dual) {
//         for(int d = complete_comp.dim-1; d >= 0; d--) {
//             for(int i = 0; i < complete_comp.ncells; i++) {
//                 if(complete_comp.get_dim(i) == d) {
                                        
//                     auto range = complete_comp.get_cofacet_range(i);
//                     for(auto it = range.first; it != range.second; it++) {
                        
//                         if(include[*it]) {
//                             include[i] = true;
//                         }
                        
//                     }
//                 }
//             }
//         }
//     } else {
//         for(int d = 1; d <= complete_comp.dim; d++) {
//             for(int i = 0; i < complete_comp.ncells; i++) {
//                 if(complete_comp.get_dim(i) == d) {
                                        
//                     auto range = complete_comp.get_facet_range(i);
//                     for(auto it = range.first; it != range.second; it++) {
                        
//                         if(include[*it]) {
//                             mask[i] = true;
//                         }
                        
//                     }
//                 }
//             }
//         }
//     }
      
    
//     std::vector<int> complete_to_mask(complete_comp.ncells);
//     int index = 0;
//     for(int i = 0; i < complete_comp.ncells; i++) {
//         if(include[i]) {
//             complete_to_mask[i] = index;
//             index++;
//         }
//     }
    
//     CellComplex mask_comp(complete_comp.dim, true, oriented);
    
//     for(int i = 0; i < complete_comp.ncells; i++) {
//         if(include[i]) {
            
//             std::vector<int> facets;
//             auto facet_range = complete_comp.get_facet_range(i);
//             for(auto it = facet_range.first; it != facet_range.second; it++) {
//                 facets.push_back(complete_to_mask[*it]);
//             }
//             std::vector<int> coeffs;
//             auto coeff_range = complete_comp.get_coeff_range(i);
//             coeffs.assign(coeff_range.first, coeff_range.second);
            
//             mask_comp.add_cell(complete_comp.get_label(i), complete_comp.get_dim(i), facets, coeffs);
//         }
//     }
       
//     mask_comp.construct_cofacets();
//     mask_comp.make_compressed(); 
        
//     return mask_comp;
    
// }


std::unordered_set<int> get_boundary_pixels(std::unordered_set<int> &pixels, std::vector<int> &shape) {
 
    std::unordered_set<int> boundary;
        
    int nrows = shape[0];
    int ncols = shape[1];
    
    for(auto pix: pixels) {
        int col = pix % ncols;
        int row = (pix - col) / ncols;
                
        if(row == 0 || row == nrows-1 || col == 0 || col == ncols-1) {
            boundary.insert(pix);
        } else if(!pixels.count(ncols*(row-1)+col) || !pixels.count(ncols*(row-1)+col+1) 
                 || !pixels.count(ncols*row+col+1) || !pixels.count(ncols*(row+1)+col+1)
                 || !pixels.count(ncols*(row+1)+col) || !pixels.count(ncols*(row+1)+col-1)
                 || !pixels.count(ncols*row+col-1) || !pixels.count(ncols*(row-1)+col-1)) {
            boundary.insert(pix);
        }
    }
    
    return boundary;
}

    
#endif // CUBICALCOMPLEX