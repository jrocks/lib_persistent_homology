#ifndef CUBICALCOMPLEX_HPP
#define CUBICALCOMPLEX_HPP
 
#include <map>
#include <algorithm>
#include <vector>   
    
#include "cell_complex.hpp"
#include "filtration.hpp"
    
    

    

CellComplex construct_cubical_complex(std::vector<int> &shape, bool oriented, bool dual) {
    
//     int dim = shape.size();
    
//     CellComplex comp(dim, true, oriented);
    
    int dim = shape.size();
    
    CellComplex comp(dim, true, false);
    
    // In d dimensions there are (d choose k) cells of dimension k at each vertex. 
    // Each dimension k cell as 2^k vertices.
    
    
    int size = 1;
    std::vector<int> lengths(dim);
    for(int d = 0; d < dim; d++) {
        lengths[d] = dual ? shape[d] + 1: shape[d];
        size *= lengths[d];
        
    }
    
    auto coords_to_index = [&lengths, dim](std::vector<int> multi_index) {
        
        int index = multi_index[0];
        for(int d = 1; d < dim; d++) {
            index *= lengths[d];
            index += multi_index[d];
        }
        
        return index;
        
    };
    
    int ncells = 0;
        
    // Map of cell dimension to map of cell vertex list to cell index
    std::vector<std::map<std::vector<int>, int> > verts_to_cell(dim+2, std::map<std::vector<int>, int>());
    
    std::vector<int> multi_index(dim, 0);
    for(int k = 0; k <= dim; k++) {
        
        // py::print("k", k, py::arg("flush")=true);
        
//         for(auto l: verts_to_cell[k]) {
        
//             py::print(l, py::arg("flush")=true);
//         }
        
        for(int i = 0; i < size; i++) {
            
//             py::print("i", i, py::arg("flush")=true);
            
//             py::print(i, multi_index, py::arg("flush")=true);
            
            // For dimension k cells, iterate through every combination of k coordinates
            std::vector<bool> mask(k, true);
            mask.resize(dim, false);
            do { 
                
                // Find coords to consider
                bool skip = false;
                std::vector<int> coords;
                for(int m = 0; m < dim; m++) {
                    if(mask[m]) {
                        
                        if(multi_index[m] == lengths[m] - 1) {
                            skip = true;
                            break;
                        }
                        
                        coords.push_back(m);
                    }
                }
                
                // py::print("coords", coords, py::arg("flush")=true);
                
                if(skip) {
                    continue;
                }

                // List of vertices comprising cell
                std::vector<int> verts;
                          
                // this is not correct and i don't use the coords list
                // Go through every subset of coordinates and add vertex
                for(int p = 0; p < (int) pow(2, k); p++) {
                    
                    std::vector<int> tmp_index = multi_index;
                    for(int q = 0; q < k; q++) {
                        if ((p & (1 << q)) != 0) {
                            tmp_index[coords[q]]++;
                        }
                    }
                    
                    verts.push_back(coords_to_index(tmp_index));
                    
                }
                
                // py::print("verts", verts, py::arg("flush")=true);
                
                
                // Find facets
                std::vector<int> facets;
                
                std::vector<bool> facet_mask((int) pow(2, k-1), true);
                facet_mask.resize(verts.size(), false);
                do { 
                    
                    // List of vertices comprising facet
                    std::vector<int> facet_verts;

                    // Go through every subset of coordinates and add vertex
                    for(std::size_t p = 0; p < verts.size(); p++) {
                        if(facet_mask[p]) {
                            facet_verts.push_back(verts[p]);
                        }
                    }
                    
                    
                    // py::print("facet_verts", facet_verts, py::arg("flush")=true);
                    
                    if(verts_to_cell[k].count(facet_verts)) {
                        facets.push_back(verts_to_cell[k][facet_verts]);
                    }
                    
                } while(std::prev_permutation(facet_mask.begin(), facet_mask.end()));
                
                
                // py::print("facets", facets, py::arg("flush")=true);
                
                
                int label = verts_to_cell[k+1].size();
                
                std::vector<int> coeffs;
                
                
                // py::print(label, k, facets, coeffs, py::arg("flush")=true);
                
                
                comp.add_cell(label, k, facets, coeffs);
                
                // Store list of cells
                verts_to_cell[k+1][verts] = ncells;
                ncells++;
                
            
            } while(std::prev_permutation(mask.begin(), mask.end()));
            
            
            
            for(int d = dim-1; d >= 0; d--) {
                if(multi_index[d] == lengths[d] - 1) {
                    multi_index[d] = 0;
                } else {
                    multi_index[d]++;
                    break;
                }
                
            }
            
            
        }
        
        
        verts_to_cell[k].clear();
        
    }
    

    
    
    // cubical complexes are labeled according to corresponding pixel
    
    
//     if(dim == 2) {
//         int nrows = dual ? shape[0] + 1: shape[0];
//         int ncols = dual ? shape[1] + 1: shape[1];
        
//         int nverts = nrows*ncols;
//         int nhedges = nrows*(ncols-1);

        
//         // vertices
//         for(int i = 0; i < nrows*ncols; i++) {
//             std::vector<int> facets;
//             std::vector<int> coeffs;
//             comp.add_cell(i, 0, facets, coeffs);
//         }
        
//         // horizontal edges
//         for(int i = 0; i < nrows; i++) {
//             for(int j = 0; j < ncols-1; j++) {
//                 std::vector<int> facets;                
//                 facets.push_back(ncols*i + j);
//                 facets.push_back(ncols*i + j+1);
                
//                 std::vector<int> coeffs;
//                 if(oriented) {
//                     coeffs.push_back(-1);
//                     coeffs.push_back(1);
//                 }
                
//                 comp.add_cell((ncols-1)*i + j, 1, facets, coeffs);
                
//             }
//         }
        
//         // vertical edges
//         for(int i = 0; i < nrows-1; i++) {
//             for(int j = 0; j < ncols; j++) {
//                 std::vector<int> facets;                
//                 facets.push_back(ncols*i + j);
//                 facets.push_back(ncols*(i+1) + j);
                
//                 std::vector<int> coeffs;
//                 if(oriented) {
//                     coeffs.push_back(-1);
//                     coeffs.push_back(1);
//                 }
                
//                 comp.add_cell(nhedges + ncols*i + j, 1, facets, coeffs);
                
//             }
//         }
        
//         // faces
//         for(int i = 0; i < nrows-1; i++) {
//             for(int j = 0; j < ncols-1; j++) {
//                 std::vector<int> facets;                
//                 facets.push_back(nverts + (ncols-1)*i + j);
//                 facets.push_back(nverts + nhedges + ncols*i + j+1);
//                 facets.push_back(nverts + (ncols-1)*(i+1) + j);
//                 facets.push_back(nverts + nhedges + ncols*i + j);
                
//                 std::vector<int> coeffs;
//                 if(oriented) {
//                     coeffs.push_back(1);
//                     coeffs.push_back(1);
//                     coeffs.push_back(-1);
//                     coeffs.push_back(-1);
//                 }
                
//                 comp.add_cell((ncols-1)*i + j, 2, facets, coeffs);
                
//             }
//         }
//     } else if(dim == 3) {
        
//         int size = 1;
//         std::vector<int> lengths(dim);
//         for(int d = 0; d < dim; d++) {
//             lengths[d] = dual ? shape[d] + 1: shape[d];
//             size *= lengths[d];
//         }
        
//         // vertices
//         for(int i = 0; i < size; i++) {
//             std::vector<int> facets;
//             std::vector<int> coeffs;
//             comp.add_cell(i, 0, facets, coeffs);
//         }
        
//         // first edge batch
//         for(int i = 0; i < lengths[0]; i++) {
//             for(int j = 0; j < lengths[1]; j++) {
//                 for(int k = 0; k < lengths[2]-1; k++){
//                     std::vector<int> facets; 
//                     facets.push_back(lengths[1]*lengths[2]*i + lengths[2]*j + k);
//                     facets.push_back(lengths[1]*lengths[2]*i + lengths[2]*j + k+1);

//                     std::vector<int> coeffs;
//                     if(oriented) {
//                         coeffs.push_back(-1);
//                         coeffs.push_back(1);
//                     }

//                     comp.add_cell(lengths[1]*(lengths[2]-1)*i + (lengths[2]-1)*j + k, 1, facets, coeffs);
                    
//                 }
//             }
//         }
        
//         // second edge batch
//         int nedges = lengths[0]*lengths[1]*(lengths[2]-1);
//         for(int i = 0; i < lengths[0]; i++) {
//             for(int j = 0; j < lengths[1]-1; j++) {
//                 for(int k = 0; k < lengths[2]; k++){
//                     std::vector<int> facets; 
//                     facets.push_back(lengths[1]*lengths[2]*i + lengths[2]*j + k);
//                     facets.push_back(lengths[1]*lengths[2]*i + lengths[2]*(j+1) + k);

//                     std::vector<int> coeffs;
//                     if(oriented) {
//                         coeffs.push_back(-1);
//                         coeffs.push_back(1);
//                     }

//                     comp.add_cell(nedges + (lengths[1]-1)*lengths[2]*i + lengths[2]*j + k, 1, facets, coeffs);
//                 }
//             }
//         }
        
//         // third edge batch
//         int nedges += lengths[0]*(lengths[1]-1)*lengths[2];
//         for(int i = 0; i < lengths[0]-1; i++) {
//             for(int j = 0; j < lengths[1]; j++) {
//                 for(int k = 0; k < lengths[2]; k++){
//                     std::vector<int> facets; 
//                     facets.push_back(lengths[1]*lengths[2]*i + lengths[2]*j + k);
//                     facets.push_back(lengths[1]*lengths[2]*(i+1) + lengths[2]*j + k);

//                     std::vector<int> coeffs;
//                     if(oriented) {
//                         coeffs.push_back(-1);
//                         coeffs.push_back(1);
//                     }

//                     comp.add_cell(nedges + lengths[1]*lengths[2]*i + lengths[2]*j + k, 1, facets, coeffs);
//                 }
//             }
//         }
        
        
//         int nverts = lengths[0]*lengths[1]*lengths[2];
//         // first face batch
//         for(int i = 0; i < lengths[0]-1; i++) {
//             for(int j = 0; j < lengths[1]-1; j++) {
//                 for(int k = 0; k < lengths[2]; k++){

//                     std::vector<int> facets;                
//                     facets.push_back(nverts + (ncols-1)*i + j);
//                     facets.push_back(nverts + nhedges + ncols*i + j+1);
//                     facets.push_back(nverts + (ncols-1)*(i+1) + j);
//                     facets.push_back(nverts + nhedges + ncols*i + j);

//                     std::vector<int> coeffs;
//                     if(oriented) {
//                         coeffs.push_back(1);
//                         coeffs.push_back(1);
//                         coeffs.push_back(-1);
//                         coeffs.push_back(-1);
//                     }

//                     comp.add_cell((ncols-1)*i + j, 2, facets, coeffs);
                    
//                 }
//             }
//         }
        
//     }
    
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

    
#endif // CUBICALCOMPLEX_HPP