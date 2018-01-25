#ifndef CELLCOMPLEX
#define CELLCOMPLEX
    
#include <vector>
#include <unordered_map>
#include <iostream>
    
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
    
    std::vector<int>* get_facets(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        return new std::vector<int>(facets.begin()+facet_ind[index], facets.begin()+facet_ind[index+1]);
    }
                      
    std::unordered_map<int, int>* get_coeffs(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        std::unordered_map<int, int>* coeff_map = new std::unordered_map<int, int>();
                
        for(int i = facet_ind[index]; i < facet_ind[index+1]; i++) {
            (*coeff_map)[facets[i]] = coeffs[i];
        }
                
        return coeff_map;
        
        // return new std::vector<int>(coeffs.begin()+facet_ind[index], coeffs.begin()+facet_ind[index+1]);
        
    }
    
    std::vector<int>* get_cofacets(int alpha) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        return new std::vector<int>(cofacets.begin()+cofacet_ind[index], cofacets.begin()+cofacet_ind[index+1]);
    }
    
    void get_facets(int alpha, std::vector<int>::iterator &begin, std::vector<int>::iterator &end) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        begin = facets.begin()+facet_ind[index];
        end = facets.begin()+facet_ind[index+1];
    }
    
    void get_cofacets(int alpha, std::vector<int>::iterator &begin, std::vector<int>::iterator &end) {
        int index = ordered ? alpha : cell_to_index[alpha];
        
        begin = cofacets.begin()+cofacet_ind[index];
        end = cofacets.begin()+cofacet_ind[index+1];
    }
    
    std::vector<int>* get_cells() {
        std::vector<int>* cells = new std::vector<int>();
        
        // convert to ranged based for loop
        if(ordered) {
            for(int i = 0; i < ncells; i++) {
                cells->push_back(i);
            }
        } else {
            for(auto it=cell_to_index.begin(); it!=cell_to_index.end(); ++it) {
                cells->push_back(it->first);
            }
        }
        
        return cells;
    }
    
    void construct_cofacets() {
         
        // std::cout << "hi" <<std::endl;
        
        // convert to ranged based for loop
        std::vector<std::vector<int> > cell_list(ncells, std::vector<int>());
        if(ordered) {
            for(int i = 0; i < ncells; i++) {
                
                for(auto j=facets.begin()+facet_ind[i]; j!=facets.begin()+facet_ind[i+1]; j++) {
                    cell_list[*j].push_back(i);
                }
            }
        } else {
            for(auto i=cell_to_index.begin(); i!=cell_to_index.end(); ++i) {
                for(auto j=facets.begin()+facet_ind[i->second]; j!=facets.begin()+facet_ind[(i->second)+1]; j++) {
                    cell_list[cell_to_index[*j]].push_back(i->first);
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


#endif // CELLCOMPLEX