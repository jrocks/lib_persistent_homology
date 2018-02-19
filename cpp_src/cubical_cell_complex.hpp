#ifndef CUBICALCOMPLEX
#define CUBICALCOMPLEX
 
    
#include <vector>
#include "cell_complex.hpp"
#include "filtration.hpp"
    


CellComplex construct_cubical_complex(std::vector<int> &shape, bool oriented, bool dual) {
    
    int dim = shape.size();
    
    CellComplex comp(dim, true, oriented);
    
    if(dim == 2 and !dual) {
        int nrows = shape[0];
        int ncols = shape[1];
        
        int ilabel = 0;
        // vertices
        for(int i = 0; i < nrows*ncols; i++, ilabel++) {
            std::vector<int> facets;
            std::vector<int> coeffs;
            comp.add_cell(ilabel, 0, facets, coeffs);
        }
        
        // horizontal edges
        for(int i = 0; i < nrows; i++) {
            for(int j = 0; j < ncols-1; j++, ilabel++) {
                std::vector<int> facets;                
                facets.push_back(ncols*i + j);
                facets.push_back(ncols*i + j+1);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(ilabel, 1, facets, coeffs);
                
            }
        }
        
        // vertical edges
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols; j++, ilabel++) {
                std::vector<int> facets;                
                facets.push_back(ncols*i + j);
                facets.push_back(ncols*(i+1) + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(ilabel, 1, facets, coeffs);
                
            }
        }
        
        // faces
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols-1; j++, ilabel++) {
                std::vector<int> facets;                
                facets.push_back(nrows*ncols + (ncols-1)*i + j);
                facets.push_back(nrows*ncols + (ncols-1)*nrows + ncols*i + j+1);
                facets.push_back(nrows*ncols + (ncols-1)*(i+1) + j);
                facets.push_back(nrows*ncols + (ncols-1)*nrows + ncols*i + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(1);
                    coeffs.push_back(1);
                    coeffs.push_back(-1);
                    coeffs.push_back(-1);
                }
                
                comp.add_cell(ilabel, 2, facets, coeffs);
                
            }
        }
    } else if(dim == 2 and dual) {
        int nrows = shape[0]+1;
        int ncols = shape[1]+1;
        
        int nhedges = nrows*(ncols-1);
        int nvedges = (nrows-1)*ncols;
        int nfaces = (nrows-1)*(ncols-1);
        
        int ilabel = 0;
        // faces
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols-1; j++, ilabel++) {
                std::vector<int> facets;                
                facets.push_back(nfaces + (ncols-1)*i + j);
                facets.push_back(nfaces + nhedges + ncols*i + j+1);
                facets.push_back(nfaces + (ncols-1)*(i+1) + j);
                facets.push_back(nfaces + nhedges + ncols*i + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(1);
                    coeffs.push_back(1);
                    coeffs.push_back(-1);
                    coeffs.push_back(-1);
                }
                
                comp.add_cell(ilabel, 2, facets, coeffs);
                
            }
        }
        
        // horizontal edges
        for(int i = 0; i < nrows; i++) {
            for(int j = 0; j < ncols-1; j++, ilabel++) {
                std::vector<int> facets;                
                facets.push_back(nfaces + nhedges + nvedges + ncols*i + j);
                facets.push_back(nfaces + nhedges + nvedges + ncols*i + j+1);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(ilabel, 1, facets, coeffs);
                
            }
        }
        
        // vertical edges
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols; j++, ilabel++) {
                std::vector<int> facets;                
                facets.push_back(nfaces + nhedges + nvedges + ncols*i + j);
                facets.push_back(nfaces + nhedges + nvedges + ncols*(i+1) + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(ilabel, 1, facets, coeffs);
                
            }
        }
        
        // vertices
        for(int i = 0; i < nrows*ncols; i++, ilabel++) {
            std::vector<int> facets;
            std::vector<int> coeffs;
            comp.add_cell(ilabel, 0, facets, coeffs);
        }

    }
    
    comp.construct_cofacets();
    comp.make_compressed(); 
        
    return comp;
    
}



// CellComplex construct_masked_cubical_complex(std::vector<bool> &mask, std::vector<int> &shape, bool oriented, bool dual) {
    
//     int dim = shape.size();
    
//     int ncells = 0;
//     std::vector<int> dims;
//     std::vector<int> labels;
//     std::vector<std::vector<int> > facet_list;
//     std::vector<std::vector<int> > coeff_list;
    
//     if(dim == 2 and dual) {
//         int nrows = shape[0]+1;
//         int ncols = shape[1]+1;
        
//         int nverts = nrows*ncols;
//         int nhedges = nrows*(ncols-1);
//         int nvedges = (nrows-1)*ncols;
//         int nfaces = (nrows-1)*(ncols-1);
        
//         ncells = nverts+nhedges+nvedges+nfaces;
//         dims.resize(ncells);
//         labels.resize(ncells);
//         facet_list.resize(ncells);
//         if(!oriented) {
//             coeff_list.resize(ncells);
//         }
        
//         int icell = 0;
//         // faces
//         for(int i = 0; i < nrows-1; i++) {
//             for(int j = 0; j < ncols-1; j++, icell++) {
                
//                 labels[icell] = icell;
//                 dims[icell] = 2;
                
//                 facet_list[icell].push_back(nfaces + (ncols-1)*i + j);
//                 facet_list[icell].push_back(nfaces + nhedges + ncols*i + j+1);
//                 facet_list[icell].push_back(nfaces + (ncols-1)*(i+1) + j);
//                 facet_list[icell].push_back(nfaces + nhedges + ncols*i + j);
                
//                 if(oriented) {
//                     coeff_list[icell].push_back(1);
//                     coeff_list[icell].push_back(1);
//                     coeff_list[icell].push_back(-1);
//                     coeff_list[icell].push_back(-1);
//                 }
                
//             }
//         }
        
//         // horizontal edges
//         for(int i = 0; i < nrows; i++) {
//             for(int j = 0; j < ncols-1; j++, icell++) {
                
//                 labels[icell] = icell;
//                 dims[icell] = 1;
                
//                 facet_list[icell].push_back(nfaces + nhedges + nvedges + ncols*i + j);
//                 facet_list[icell].push_back(nfaces + nhedges + nvedges + ncols*i + j+1);
                
//                 if(oriented) {
//                     coeff_list[icell].push_back(-1);
//                     coeff_list[icell].push_back(1);
//                 }
                
//             }
//         }
        
//         // vertical edges
//         for(int i = 0; i < nrows-1; i++) {
//             for(int j = 0; j < ncols; j++, icell++) {
                
//                 labels[icell] = icell;
//                 dims[icell] = 1;
                
//                 facet_list[icell].push_back(nfaces + nhedges + nvedges + ncols*i + j);
//                 facet_list[icell].push_back(nfaces + nhedges + nvedges + ncols*(i+1) + j);
                
//                 if(oriented) {
//                     coeff_list[icell].push_back(-1);
//                     coeff_list[icell].push_back(1);
//                 }
                
//             }
//         }
        
//         // vertices
//         for(int i = 0; i < nrows*ncols; i++, icell++) {
            
//             labels[icell] = icell;
//             dims[icell] = 0;
            
//         }

//     }
    
        
//     std::vector<bool> include(ncells, false);
//     for(int i = 0; i < ncells; i++) {
        
//         if((dual && dims[i] == dim) 
//           || (!dual && dims[i] == 0)) {
            
//             include[i] = !mask[labels[i]];
//         }
//     }
    
//     if(dual) {
//         // Cells are already sorted from highest to lowest dimension
//         for(int i = 0; i < ncells; i++) {
//             if(include[i]) {             
//                 for(auto alpha: facet_list[i]) {
//                     include[alpha] = true;
//                 }
//             }
//         }
//     } else {

//     }
    
//     std::vector<int> complete_to_mask(ncells);
//     int index = 0;
//     for(int i = 0; i < ncells; i++) {
//         if(include[i]) {
//             complete_to_mask[i] = index;
//             index++;
//         }
//     }
    
//     CellComplex mask_comp(dim, true, oriented);
    
//     for(int i = 0; i < ncells; i++) {
//         if(include[i]) {
            
//             std::vector<int> facets;
//             for(auto alpha: facet_list[i]) {
//                 facets.push_back(complete_to_mask[alpha]);
//             }
            
//             std::vector<int> coeffs;
//             if(oriented) {
//                 for(auto alpha: coeff_list[i]) {
//                     coeffs.push_back(complete_to_mask[alpha]);
//                 }
//             }
            
//             mask_comp.add_cell(labels[i], dims[i], facets, coeffs);
//         }
        
//         facet_list[i].clear();
//         if(oriented) {
//             coeff_list[i].clear();
//         }
        
//     }
       
//     mask_comp.construct_cofacets();
//     mask_comp.make_compressed(); 
        
//     return mask_comp;
    
// }



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
    

std::unordered_set<int> convert_to_pixels(std::unordered_set<int> &feature, Filtration &filt, CellComplex &comp, bool dual) {
        
    std::unordered_set<int> pixels;
    
    for(auto s: feature) {
        if(comp.get_dim(s) == 0) {
            if(!dual) {
                pixels.insert(comp.get_label(s));
            } else {
                
                std::unordered_set<int> cofaces = get_star(s, false, comp, comp.dim);
                for(auto c: cofaces) {
                    
                    std::unordered_set<int> verts = get_star(c, true, comp, 0);
                    
                    std::unordered_set<int> ucostar_verts;
                    for(auto v: verts) {
                        if(filt.get_star(c) == filt.get_star(v)) {
                            ucostar_verts.insert(v);
                        }
                    }
                    
                    if(ucostar_verts.empty()) {
                        int min_vert = *(verts.begin());
                        for(auto v: verts) {
                            if(filt.get_filtration_order(v) < filt.get_filtration_order(min_vert)) {
                                min_vert = v;
                            }
                        }
                        
                        if(feature.count(min_vert)) {
                            pixels.insert(comp.get_label(c));
                        }
                        
                    } else {
                        
                        bool issubset = true;
                        for(auto v: ucostar_verts) {
                            if(feature.count(v) == 0) {
                                issubset = false;
                                break;
                            }
                            
                        }
                        
                        if(issubset) {
                            
                            pixels.insert(comp.get_label(c));
                            
                        } else {
                            int min_vert = *(ucostar_verts.begin());
                            for(auto v: ucostar_verts) {
                                if(filt.get_filtration_order(v) < filt.get_filtration_order(min_vert)) {
                                    min_vert = v;
                                }
                            }
                            
                            if(feature.count(min_vert)) {
                                pixels.insert(comp.get_label(c));
                            }
                            
                        }
                        
                    }
                                            
                }
                
            }
        } else if(comp.get_dim(s) == 1) {
            if(!dual) {
                auto range = comp.get_facet_range(s);
                for(auto it = range.first; it != range.second; it++) {
                    pixels.insert(comp.get_label(*it));
                }
            } else {
                int ucostar = filt.get_star(s);
                pixels.insert(comp.get_label(ucostar));
            }
            
        } else if(comp.get_dim(s) == 2) {
            
            if(!dual) {  
                auto rangei = comp.get_facet_range(s);
                for(auto iti = rangei.first; iti != rangei.second; iti++) {
                    auto rangej = comp.get_facet_range(*iti);
                    for(auto itj = rangej.first; itj != rangej.second; itj++) {
                        pixels.insert(comp.get_label(*itj));
                    }
                }
                
            } else {
                pixels.insert(comp.get_label(s));
            }
            
        }
    }
    
    
    return pixels;
    
}

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