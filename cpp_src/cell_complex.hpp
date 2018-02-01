#ifndef CELLCOMPLEX
#define CELLCOMPLEX
    
#include <pybind11/pybind11.h>
namespace py = pybind11;
    
#include <vector>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <queue>
        
class CellComplex {
    
public:
    
    int dim;    
    int ncells;
    int nfacets;
    bool regular;
    bool oriented;
    bool ordered;
            
    std::unordered_map<int,int> cell_to_index;
    
    std::vector<int> dims;
    std::vector<int> facet_ind;
    std::vector<int> facets;
    std::vector<int> coeffs;
    
    std::vector<int> cofacet_ind;
    std::vector<int> cofacets;
    
    CellComplex(int dim, bool regular, bool oriented, bool ordered) :
        dim(dim), ncells(0), nfacets(0),
        regular(regular), oriented(oriented), ordered(ordered) {
            
        facet_ind.push_back(0);            
    }
    
    void add_cell(int dim, std::vector<int> &cell_facets, int label, std::vector<int> &cell_coeffs) {
       
        if(!ordered) {
            cell_to_index[label] = ncells;
        }
        
        dims.push_back(dim);
        facet_ind.push_back(facet_ind[ncells] + cell_facets.size());
        facets.insert(facets.end(), cell_facets.begin(), cell_facets.end());
        if(!regular || oriented) {
            coeffs.insert(coeffs.end(), cell_coeffs.begin(), cell_coeffs.end());
        }
        
        ncells++;
        nfacets += cell_facets.size();
        
    }
    
    void make_compressed() {
        dims.shrink_to_fit();
        facet_ind.shrink_to_fit();
        facets.shrink_to_fit();
        coeffs.shrink_to_fit();
        cofacet_ind.shrink_to_fit();
        cofacets.shrink_to_fit();
    }
    
    int get_dim(int alpha) {
        if(ordered) {
            return dims[alpha];
        } else {
            return dims[cell_to_index[alpha]];
        }
    }
    
    std::vector<int> get_facets(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        return std::vector<int>(facets.begin()+facet_ind[index], facets.begin()+facet_ind[index+1]);
    }
                      
    std::unordered_map<int, int> get_coeffs(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        std::unordered_map<int, int> coeff_map;
                
        for(int i = facet_ind[index]; i < facet_ind[index+1]; i++) {
            coeff_map[facets[i]] = coeffs[i];
        }
                
        return coeff_map;
                
    }
    
    std::vector<int> get_cofacets(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        return std::vector<int>(cofacets.begin()+cofacet_ind[index], cofacets.begin()+cofacet_ind[index+1]);
    }
    
    
    
    std::pair<std::vector<int>::iterator, std::vector<int>::iterator> get_facet_range(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        return std::make_pair(facets.begin()+facet_ind[index], facets.begin()+facet_ind[index+1]);
    }
    
    std::pair<std::vector<int>::iterator, std::vector<int>::iterator> get_cofacet_range(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        return std::make_pair(cofacets.begin()+cofacet_ind[index], cofacets.begin()+cofacet_ind[index+1]);
    }
    
    
    
    std::vector<int> get_cells() {
        std::vector<int> cells;
        
        // convert to ranged based for loop
        if(ordered) {
            for(int i = 0; i < ncells; i++) {
                cells.push_back(i);
            }
        } else {
            for(auto it=cell_to_index.begin(); it!=cell_to_index.end(); ++it) {
                cells.push_back(it->first);
            }
        }
        
        return cells;
    }
    
    void construct_cofacets() {
                 
        std::vector<std::vector<int> > cell_list(ncells, std::vector<int>());
        if(ordered) {
            for(int i = 0; i < ncells; i++) {
                
                for(auto j=facets.begin()+facet_ind[i]; j!=facets.begin()+facet_ind[i+1]; j++) {
                    cell_list[*j].push_back(i);
                }
            }
        } else {
            // for(auto i=cell_to_index.begin(); i!=cell_to_index.end(); ++i) {
            for(auto i: cell_to_index) {
                
                // for(auto j=facets.begin()+facet_ind[i->second]; j!=facets.begin()+facet_ind[(i->second)+1]; j++) {
                //     cell_list[cell_to_index[*j]].push_back(i->first);
                // }
                
                auto range = get_facet_range(i.first);
                for(auto itj = range.first; itj != range.second; itj++) {
                    cell_list[cell_to_index[*itj]].push_back(i.first);
                }
            }
        }
        
        cofacet_ind.push_back(0);
        for(int i = 0; i < ncells; i++) {
            cofacet_ind.push_back(cofacet_ind[i] + cell_list[i].size());
            cofacets.insert(cofacets.end(), cell_list[i].begin(), cell_list[i].end());
        }
        
    }
       
};


CellComplex construct_cubical_complex(std::vector<int> &shape, bool oriented, bool dual) {
    
    int dim = shape.size();
    
    CellComplex comp(dim, true, oriented, true);
    
    if(dim == 2 and !dual) {
        int nrows = shape[0];
        int ncols = shape[1];
        
        // vertices
        for(int i = 0; i < nrows*ncols; i++) {
            std::vector<int> facets;
            std::vector<int> coeffs;
            comp.add_cell(0, facets, -1, coeffs);
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
                
                comp.add_cell(1, facets, -1, coeffs);
                
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
                
                comp.add_cell(1, facets, -1, coeffs);
                
            }
        }
        
        // faces
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols-1; j++) {
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
                
                comp.add_cell(2, facets, -1, coeffs);
                
            }
        }
    } else if(dim == 2 and dual) {
        int nrows = shape[0]+1;
        int ncols = shape[1]+1;
        
        int nhedges = nrows*(ncols-1);
        int nvedges = (nrows-1)*ncols;
        int nfaces = (nrows-1)*(ncols-1);
        
        // faces
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols-1; j++) {
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
                
                comp.add_cell(2, facets, -1, coeffs);
                
            }
        }
        
        // horizontal edges
        for(int i = 0; i < nrows; i++) {
            for(int j = 0; j < ncols-1; j++) {
                std::vector<int> facets;                
                facets.push_back(nfaces + nhedges + nvedges + ncols*i + j);
                facets.push_back(nfaces + nhedges + nvedges + ncols*i + j+1);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(1, facets, -1, coeffs);
                
            }
        }
        
        // vertical edges
        for(int i = 0; i < nrows-1; i++) {
            for(int j = 0; j < ncols; j++) {
                std::vector<int> facets;                
                facets.push_back(nfaces + nhedges + nvedges + ncols*i + j);
                facets.push_back(nfaces + nhedges + nvedges + ncols*(i+1) + j);
                
                std::vector<int> coeffs;
                if(oriented) {
                    coeffs.push_back(-1);
                    coeffs.push_back(1);
                }
                
                comp.add_cell(1, facets, -1, coeffs);
                
            }
        }
        
        // vertices
        for(int i = 0; i < nrows*ncols; i++) {
            std::vector<int> facets;
            std::vector<int> coeffs;
            comp.add_cell(0, facets, -1, coeffs);
        }

    }
    
    comp.construct_cofacets();
    comp.make_compressed(); 
        
    return comp;
    
}



    
bool check_boundary_op(CellComplex &comp) {
    bool valid = true;
    
    for(auto i: comp.get_cells()) {
        if(comp.get_dim(i) > 0) {
            if(!comp.regular || comp.oriented) {
                std::unordered_map<int, int> sub_face_coeffs;
                
                std::unordered_map<int, int> coeffsi = comp.get_coeffs(i);
                for(auto j: coeffsi) {
                    std::unordered_map<int, int> coeffsj = comp.get_coeffs(j.first);
                    for(auto k: coeffsj) {
                        if(sub_face_coeffs.count(k.first)) {
                            sub_face_coeffs[k.first] += coeffsi[j.first] * coeffsj[k.first];
                        } else {
                            sub_face_coeffs[k.first] = coeffsi[j.first] * coeffsj[k.first];
                        }
                    }
                    
                }
                
                for(auto j: sub_face_coeffs) {
                    if((comp.oriented && j.second != 0) || (!comp.oriented && j.second % 2 != 0) ) {
                        py::print("Error:", i, j.first, j.second);
                        valid = false;
                    }
                }

                
            } else {
                std::unordered_set<int> sub_faces;
                
                auto rangei = comp.get_facet_range(i);
                
                for(auto iti = rangei.first; iti != rangei.second; iti++) {
                    
                    auto rangej = comp.get_facet_range(*iti);
                                    
                    for(auto itj = rangej.first; itj != rangej.second; itj++) {
                        if(sub_faces.count(*itj)) {
                            sub_faces.erase(*itj);
                        } else {
                            sub_faces.insert(*itj);
                        }
                    }
                    
                }
                
                
                if(!sub_faces.empty()) {
                    py::print("Error:", i);
                    valid = false;
                    
                }
                
                
            }
            
        }
    }
        
    return valid;
    
}



std::unordered_set<int> get_star(int alpha, bool co, CellComplex &comp, int target_dim) {
    
    std::unordered_set<int> star;
    
    std::queue<int> Q;
    Q.push(alpha);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
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
    
    return star;
 
}



#endif // CELLCOMPLEX