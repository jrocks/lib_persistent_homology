#ifndef CELLCOMPLEXALGS
#define CELLCOMPLEXALGS
   
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
    
namespace py = pybind11;
    
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
    
#include "cell_complex.hpp"
    
    
    
bool check_boundary_op(CellComplex &comp) {
    bool valid = true;
    std::vector<int>* cells = comp.get_cells();
    
    for(auto i: *cells) {
        if(comp.get_dim(i) > 0) {
            if(!comp.regular || comp.oriented) {
                std::unordered_map<int, int> sub_face_coeffs;
                
                std::unordered_map<int, int>* coeffsi = comp.get_coeffs(i);
                for(auto j: *coeffsi) {
                    std::unordered_map<int, int>* coeffsj = comp.get_coeffs(j.first);
                    for(auto k: *coeffsj) {
                        if(sub_face_coeffs.count(k.first) == 1) {
                            sub_face_coeffs[k.first] += (*coeffsi)[j.first] * (*coeffsj)[k.first];
                        } else {
                            sub_face_coeffs[k.first] = (*coeffsi)[j.first] * (*coeffsj)[k.first];
                        }
                                            }
                    delete coeffsj;
                    
                }
                delete coeffsi;
                
                for(auto j: sub_face_coeffs) {
                    if((comp.oriented && j.second != 0) || (!comp.oriented && j.second % 2 != 0) ) {
                        std::cout << "Error: " << i << std::endl;
                        valid = false;
                        
                    }
                }

                
            } else {
                std::unordered_set<int> sub_faces;
                
                std::vector<int>::iterator begini;
                std::vector<int>::iterator endi;
                comp.get_facets(i, begini, endi);
                for(std::vector<int>::iterator iti = begini; iti != endi; iti++) {
                    std::vector<int>::iterator beginj;
                    std::vector<int>::iterator endj;
                    comp.get_facets(*iti, beginj, endj);
                    
                    for(std::vector<int>::iterator itj = beginj; itj != endj; itj++) {
                        if(sub_faces.count(*itj) == 1) {
                            sub_faces.erase(*itj);
                        } else {
                            sub_faces.insert(*itj);
                        }
                    }
                    
                }
                
                
                if(!sub_faces.empty()) {
                    std::cout << "Error: " << i << std::endl;
                    valid = false;
                    
                }
                
                
            }
            
        }
    }
    
    delete cells;
    
    return valid;
    
}
    
    
    
    
// def check_boundary_op(comp):
    
//     valid = True
//     for i in comp.get_cells():
//         if comp.get_dim(i) > 0:
            
//             if not comp.regular or comp.oriented:
//                 sub_faces = {}
//                 for j in comp.get_facets(i):
//                     for k in comp.get_facets(j):
//                         sub_faces[k] = sub_faces.get(k, 0) + comp.get_coeffs(i)[j] * comp.get_coeffs(j)[k]
                        
//                 for j in sub_faces:
//                     if (comp.oriented and sub_faces[j] != 0) or (not comp.oriented and sub_faces[j] % 2 != 0):
//                         print("Error:", i, j, sub_faces[j])
//                         valid = False
//             else:
//                 sub_faces = set()
//                 for j in comp.get_facets(i):
//                     sub_faces ^= set(comp.get_facets(j))
                   
//                 if len(sub_faces) != 0:
//                     print("Error:", i, sub_faces)
//                     valid = False
                    
//     return valid
    
    
    
    
    
    
std::unordered_set<int> *get_star(int alpha, bool co, CellComplex &comp, int target_dim) {
    
    std::unordered_set<int> *star = new std::unordered_set<int>();
    
    std::queue<int> Q;
    Q.push(alpha);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        std::vector<int>::iterator begin;
        std::vector<int>::iterator end;
        // These are purposely backwards
        // Explore cofacets for star and facets for costar
        if(co) {
            comp.get_facets(a, begin, end);
        } else {
            comp.get_cofacets(a, begin, end);
        }
        
        for(std::vector<int>::iterator it = begin; it != end; it++) {
            Q.push(*it);
        }
        
        if(target_dim == -1 || comp.get_dim(a) == target_dim) {
            star->insert(a);
        }
    }
    
    return star;
 
}


std::unordered_set<int> *get_lower_star(int alpha, bool co, CellComplex &comp, int target_dim, py::array_t<double> insert_order) {
    
    auto order = insert_order.unchecked<2>();
    
    int star_index = order(alpha, STAR);
    
    std::unordered_set<int> *star = new std::unordered_set<int>();
    
    std::queue<int> Q;
    Q.push(alpha);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        if(star_index == order(a, STAR)) {
            std::vector<int>::iterator begin;
            std::vector<int>::iterator end;
            // These are purposely backwards
            // Explore cofacets for star and facets for costar
            if(co) {
                comp.get_facets(a, begin, end);
            } else {
                comp.get_cofacets(a, begin, end);
            }
            
            for(std::vector<int>::iterator it = begin; it != end; it++) {
                Q.push(*it);
            }

            if(target_dim == -1 || comp.get_dim(a) == target_dim) {
                star->insert(a);
            }
        }
    }
    
    return star;
 
}


py::array_t<double> construct_time_of_insertion_map(CellComplex &comp, py::array_t<int> vertex_order, 
                                py::array_t<double> vertex_time, bool dual) {
    
        
    auto vorder = vertex_order.unchecked<1>();
    auto vtime = vertex_time.unchecked<1>();
    
    class Comparator {
        
        CellComplex &comp;
        std::vector<std::vector<int> > &lex_val;
        bool dual;
        public:
            Comparator(CellComplex &comp, std::vector<std::vector<int> > &lex_val, const bool& dual=false):
                        comp(comp), lex_val(lex_val), dual(dual) {}
            bool operator() (const int& lhs, const int&rhs) const {
                
                // If dual, first compare upper costars (smallest goes first)
                if(dual && lex_val[lhs].back() != lex_val[rhs].back()) {
                    return lex_val[lhs].back() < lex_val[rhs].back();
                // If not dual, first compare lower stars (smallest goes first)
                } else if(!dual && lex_val[lhs].front() != lex_val[rhs].front()) {
                    return lex_val[lhs].front() < lex_val[rhs].front();
                // Second compare dimensions (smallest goes first)
                } else if(comp.get_dim(lhs) != comp.get_dim(rhs)) {
                    return comp.get_dim(lhs) < comp.get_dim(rhs);
                // Third compare length of lexicographic values (smallest goes first)
                } else if(lex_val[lhs].size() != lex_val[rhs].size()) {
                    return lex_val[lhs].size() < lex_val[rhs].size();
                // Finally compare lexicographic values (smallest goes first)
                // This will be inverted when actually constructing gradient
                } else {
                    for(unsigned int i = 0; i < lex_val[lhs].size(); i++) {
                        if(lex_val[lhs][i] != lex_val[rhs][i]) {
                            return lex_val[lhs][i] < lex_val[rhs][i];
                        }
                    }
                    
                    // If cells have identical lexicographic orderings, 
                    // then sort by raw cell index
                    return (lhs < rhs);
                } 
                           
            }
    };
        
    std::vector<int> cell_to_star(comp.ncells);
    std::vector<std::vector<int> > lex_val(comp.ncells, std::vector<int>());
    
    int target_dim = dual ? comp.dim : 0;
    for(int i = 0; i < comp.ncells; i++) {
        
        std::unordered_set<int>* star = get_star(i, !dual, comp, target_dim);
        int star_cell = -1;
                
        for(auto alpha: *star) {
                        
            lex_val[i].push_back(vorder(alpha));
            
            if(star_cell == -1 
               || (dual && vorder(alpha) < vorder(star_cell))
               || (!dual && vorder(alpha) > vorder(star_cell))) {
                
                star_cell = alpha;
            }
        }
        
        cell_to_star[i] = star_cell;
        
        // Sort from high to low
        std::sort(lex_val[i].begin(), lex_val[i].end(), std::greater<int>());
        
        delete star;
    }
    
        
    std::vector<int> cells(comp.ncells);
    for(int i = 0; i < comp.ncells; i++) {
        cells[i] = i;
    }
    
    Comparator cmp(comp, lex_val, dual);
    std::sort(cells.begin(), cells.end(), cmp);
        
    py::array_t<double> insert_order;
    insert_order.resize({comp.ncells, 3});
    auto order = insert_order.mutable_unchecked<2>();
    
    for(int i = 0; i < comp.ncells; i++) {
        int c = cells[i];
        int star = cell_to_star[c];
        order(c, TIME) = vtime(star);
        order(c, STAR) = star;
        order(c, LEX) = i;
        
    }
        
    return insert_order;
    
    
}


py::array_t<int> construct_discrete_gradient(CellComplex &comp, py::array_t<double> insert_order, bool dual) {     
    
    class Comparator {
        bool dual;
        py::detail::unchecked_reference<double, 2> order;
        public:
            Comparator(py::array_t<double> insert_order, const bool& dual=false): 
                dual(dual), order(insert_order.unchecked<2>()) {}
            bool operator() (const int& lhs, const int&rhs) const {
                
                // These are use (>) because with default (<) priority queue will always pop the largest value
                // We want the smallest value instead
                if(dual) {                    
                    return (-order(lhs, LEX) > -order(rhs, LEX));
                } else {
                    return (order(lhs, LEX) > order(rhs, LEX));
                }
              
            }
    };
        
    py::array_t<int> V;
    V.resize({comp.ncells});    
    std::fill(V.mutable_data(), V.mutable_data() + V.size(), -1);
    
    auto Vbuf = V.mutable_unchecked<1>();
    
    
    for(int i = 0; i < comp.ncells; i++) {
        if((dual && comp.get_dim(i) != comp.dim) || (!dual && comp.get_dim(i) != 0)) {
            continue;
        }
        
        // get lower/upper star/costar
        std::unordered_set<int> *star = get_lower_star(i, dual, comp, -1, insert_order);
        
        // record all unpaired facets/cofacets of each cell in star
        std::unordered_map<int, std::unordered_set<int> > unpaired;
        for(auto alpha: *star) {
            unpaired[alpha] = std::unordered_set<int>();
            
            std::vector<int>::iterator begin;
            std::vector<int>::iterator end;
            if(dual) {
                comp.get_cofacets(alpha, begin, end);
            } else {
                comp.get_facets(alpha, begin, end);
            }
            for(std::vector<int>::iterator it = begin; it != end; it++) {
                if(star->count(*it) == 1) {
                    unpaired[alpha].insert(*it);
                }
            }
            
        }
        
        Comparator cmp(insert_order, dual);
        std::priority_queue<int, std::vector<int>, Comparator> PQone(cmp);
        std::priority_queue<int, std::vector<int>, Comparator> PQzero(cmp);
        
        for(auto alpha: *star) {
            if(unpaired[alpha].empty()) {
                PQzero.push(alpha);
            } else if(unpaired[alpha].size() == 1) {
                PQone.push(alpha);
            }
        }
        
        while(!PQone.empty() || !PQzero.empty()) {
            while(!PQone.empty()) {
                int alpha = PQone.top();
                PQone.pop();
                
                if(unpaired[alpha].empty()) {
                    PQzero.push(alpha);
                    continue;
                }
                
                int beta = *(unpaired[alpha].begin());
                if(dual) {
                    Vbuf(alpha) = beta;
                } else {
                    Vbuf(beta)= alpha;
                }
                
                star->erase(alpha);
                star->erase(beta);
                unpaired.erase(alpha);
                
                for(auto gamma: *star) {
                    if(unpaired[gamma].count(alpha) == 1 || unpaired[gamma].count(beta) == 1) {
                        unpaired[gamma].erase(alpha);
                        unpaired[gamma].erase(beta);
                        
                        if(unpaired[gamma].size() == 1) {
                            PQone.push(gamma);
                        }
                    }
                }
            }
            
            
            if(!PQzero.empty()) {
                int alpha = PQzero.top();
                PQzero.pop();
                
                if(star->count(alpha) == 1) {
                    Vbuf(alpha) = alpha;
                    
                    star->erase(alpha);
                    unpaired.erase(alpha);
                    
                    for(auto gamma: *star) {
                        if(unpaired[gamma].count(alpha) == 1) {
                            unpaired[gamma].erase(alpha);

                            if(unpaired[gamma].size() == 1) {
                                PQone.push(gamma);
                            }
                        }
                    }
                    
                }
            }
            
        }
        
        delete star;
        
    }
    
    
    return V;
    
}

























































    
#endif // CELLCOMPLEXALGS