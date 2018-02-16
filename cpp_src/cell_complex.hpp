#ifndef CELLCOMPLEX
#define CELLCOMPLEX
    
#include <pybind11/pybind11.h>
namespace py = pybind11;
    
#include <vector>
#include <unordered_map>
#include <utility>
#include <queue>
        
class CellComplex {
    
public:
    
    int dim;    
    int ncells;
    int nfacets;
    bool regular;
    bool oriented;
                
    std::vector<int> labels;
    std::vector<int> dims;
    std::vector<int> facet_ind;
    std::vector<int> facets;
    std::vector<int> coeffs;
    
    std::vector<int> cofacet_ind;
    std::vector<int> cofacets;
    
    CellComplex(int dim, bool regular = true, bool oriented = false) :
        dim(dim), ncells(0), nfacets(0),
        regular(regular), oriented(oriented) {
            
        facet_ind.push_back(0);            
    }
    
    
    void add_cell(int label, int dim, std::vector<int> &cell_facets, std::vector<int> &cell_coeffs) {
       
        labels.push_back(label);
        dims.push_back(dim);
        facet_ind.push_back(facet_ind[ncells] + cell_facets.size());
        facets.insert(facets.end(), cell_facets.begin(), cell_facets.end());
        if(!regular || oriented) {
            coeffs.insert(coeffs.end(), cell_coeffs.begin(), cell_coeffs.end());
        }
        
        ncells++;
        nfacets += cell_facets.size();
        
    }
    
    
    
    void add_cell(int label, int dim, 
                  std::pair<std::vector<int>::iterator, std::vector<int>::iterator> &facet_range, 
                  std::pair<std::vector<int>::iterator, std::vector<int>::iterator> &coeff_range) {
       
        labels.push_back(label);
        dims.push_back(dim);
        
        int delta_size = facets.size();
        facets.insert(facets.end(), facet_range.first, facet_range.second);
        delta_size = facets.size() - delta_size;
        if(!regular || oriented) {
            coeffs.insert(coeffs.end(), coeff_range.first, coeff_range.second);
        }
        
        facet_ind.push_back(facet_ind[ncells] + delta_size);
        
        ncells++;
        nfacets += delta_size;
        
    }
    
    
    
    
    
    int get_dim(int alpha) {
        return dims[alpha];
    }
    
    int get_label(int alpha) {
        return labels[alpha];
    }
    
    std::vector<int> get_facets(int alpha) {        
        return std::vector<int>(facets.begin()+facet_ind[alpha], facets.begin()+facet_ind[alpha+1]);
    }
          
    std::vector<int> get_cofacets(int alpha) {        
        return std::vector<int>(cofacets.begin()+cofacet_ind[alpha], cofacets.begin()+cofacet_ind[alpha+1]);
    }
    
    std::unordered_map<int, int> get_coeffs(int alpha) {
        
        std::unordered_map<int, int> coeff_map;
                
        for(int i = facet_ind[alpha]; i < facet_ind[alpha+1]; i++) {
            if(oriented || !regular) {
                coeff_map[facets[i]] = coeffs[i];
            } else {
                coeff_map[facets[i]] = 1;
            }
        }
                
        return coeff_map;
                
    }
    
    
    
    std::pair<std::vector<int>::iterator, std::vector<int>::iterator> get_facet_range(int alpha) {
        return std::make_pair(facets.begin()+facet_ind[alpha], facets.begin()+facet_ind[alpha+1]);
    }
    
    std::pair<std::vector<int>::iterator, std::vector<int>::iterator> get_cofacet_range(int alpha) {
        return std::make_pair(cofacets.begin()+cofacet_ind[alpha], cofacets.begin()+cofacet_ind[alpha+1]);
    }
    
    std::pair<std::vector<int>::iterator, std::vector<int>::iterator> get_coeff_range(int alpha) {
        if(oriented || !regular) {
            return std::make_pair(coeffs.begin()+facet_ind[alpha], coeffs.begin()+facet_ind[alpha+1]);
        } else {
            return std::make_pair(coeffs.begin(), coeffs.begin());
        }
    }
    
    
    
    
    void make_compressed() {
        labels.shrink_to_fit();
        dims.shrink_to_fit();
        facet_ind.shrink_to_fit();
        facets.shrink_to_fit();
        coeffs.shrink_to_fit();
        cofacet_ind.shrink_to_fit();
        cofacets.shrink_to_fit();
    }
    
    void construct_cofacets() {
                 
        std::vector<std::vector<int> > cell_list(ncells, std::vector<int>());
        for(int i = 0; i < ncells; i++) {
            auto range = get_facet_range(i);
            for(auto j = range.first; j != range.second; j++) {
                cell_list[*j].push_back(i);
            }
        }

        cofacet_ind.push_back(0);
        for(int i = 0; i < ncells; i++) {
            cofacet_ind.push_back(cofacet_ind[i] + cell_list[i].size());
            cofacets.insert(cofacets.end(), cell_list[i].begin(), cell_list[i].end());
        }
        
    }
       
};





    
bool check_boundary_op(CellComplex &comp) {
    bool valid = true;
    
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) > 0) {
            if(!comp.regular || comp.oriented) {
                std::unordered_map<int, int> sub_face_coeffs;
                
                std::unordered_map<int, int> coeffsi = comp.get_coeffs(i);
                for(auto j: coeffsi) {
                    std::unordered_map<int, int> coeffsj = comp.get_coeffs(j.first);
                    for(auto k: coeffsj) {
                        sub_face_coeffs[k.first] += coeffsi[j.first] * coeffsj[k.first];
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


std::unordered_set<int> get_boundary(std::unordered_set<int> &cells, CellComplex &comp) {
    std::unordered_set<int> cycle;
    if(comp.regular) {
        for(auto c: cells) {
            auto range = comp.get_facet_range(c);
            for(auto it = range.first; it != range.second; it++) {
                int a = *it;
                if(cycle.count(a)) {
                    cycle.erase(a);
                } else {
                    cycle.insert(a);
                }
            }
        }
    }
    
    return cycle;
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