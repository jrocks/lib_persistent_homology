#ifndef MORSE
#define MORSE
   
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
    
namespace py = pybind11;
    
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <numeric>
#include <utility>
    
#include "cell_complex.hpp"
    
    
    
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



py::array_t<int> construct_vertex_filtration_order(CellComplex &comp, py::array_t<double> vertex_time, bool dual) {
    
    auto vtime = vertex_time.unchecked<1>();
    
    
    py::array_t<int> vertex_order;
    vertex_order.resize({vtime.size()});
    auto vorder = vertex_order.mutable_unchecked<1>();
    
    
    std::vector<bool> submerged(vertex_order.size(), false);
    
    std::vector<int> vertex_argsort(vertex_order.size());
    std::iota(vertex_argsort.begin(), vertex_argsort.end(), 0);
    std::sort(vertex_argsort.begin(), vertex_argsort.end(),
       [&vtime](int lhs, int rhs) {return vtime(lhs) <= vtime(rhs);});
    
    unsigned int ti = 0;
    while(ti < vertex_argsort.size()) {
        
        std::unordered_set<int> plateau;
        
        unsigned int tj = ti;
        while(true) {
            int vi = vertex_argsort[tj];
            double t = vtime(vi);
            plateau.insert(vi);
            
            if((tj+1 == vertex_argsort.size()) || (t != vtime(vertex_argsort[tj+1]))) {
                break;
            }
            
            tj++;
        }
        
        if(plateau.size() == 1) {
            int vi = *(plateau.begin());
            submerged[vi] = true;
            vorder(vi) = ti;
            ti++;
            
            continue;
        }
        
        std::unordered_map<int, double> dist;
        std::queue<int> Q;
        for(auto a: plateau) {
            
            auto rangea = dual ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto ita = rangea.first; ita != rangea.second; ita++) {
                
                int b = *ita;
                
                auto rangeb = dual ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
                    
                    int c = *itb;
                    
                    if(submerged[c]) {
                        
                        double new_dist = 1.0;
                        
                        if(!dist.count(a) || new_dist < dist[a]) {
                            dist[a] = new_dist;
                            Q.push(a);
                        }
                        
                    }

                }
                
            }
            
        }
        
        while(!Q.empty()) {
            int a = Q.front();
            Q.pop();
            
            // This differs from the python version;
            double current_dist = dist[a];
            
            if(!plateau.count(a)) {
                continue;
            }
            
            
            plateau.erase(a);
            submerged[a] = true;
            vorder(a) = ti;
            ti++;
            
            auto rangea = dual ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
            for(auto ita = rangea.first; ita != rangea.second; ita++) {
                
                int b = *ita;
                
                auto rangeb = dual ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                for(auto itb = rangeb.first; itb != rangeb.second; itb++) {
                    
                    int c = *itb;
                    
                    if(plateau.count(c)) {
                        
                        double new_dist = current_dist + 1.0;
                        
                        if(!dist.count(c) || new_dist < dist[c]) {
                            dist[c] = new_dist;
                            Q.push(c);
                        }
                    }
                    
                }
                
            }
            
            
            
        }
        
        
        while(!plateau.empty()) {
            int s = *(plateau.begin());
            plateau.erase(s);
            
            Q.push(s);
            
            while(!Q.empty()) {
                int a = Q.front();
                Q.pop();

                // This differs from the python version;
                double current_dist = dist[a];

                plateau.erase(a);
                submerged[a] = true;
                vorder(a) = ti;
                ti++;



                auto rangea = dual ? comp.get_facet_range(a) : comp.get_cofacet_range(a);
                for(auto ita = rangea.first; ita != rangea.second; ita++) {

                    int b = *ita;

                    std::vector<int>::iterator beginb;
                    std::vector<int>::iterator endb;

                    auto rangeb = dual ? comp.get_cofacet_range(b) : comp.get_facet_range(b);
                    for(auto itb = rangeb.first; itb != rangeb.second; itb++) {

                        int c = *itb;

                        if(plateau.count(c)) {

                            double new_dist = current_dist + 1.0;

                            if(!dist.count(c) || new_dist < dist[c]) {
                                dist[c] = new_dist;
                                Q.push(c);
                            }
                        }

                    }
                }
            }
            
        }
        
    }
    
    return vertex_order;
    
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


std::unordered_set<int> get_lower_star(int alpha, bool co, CellComplex &comp, int target_dim, py::array_t<double> insert_order) {
    
    auto iorder = insert_order.unchecked<2>();
    
    int star_index = iorder(alpha, STAR);
    
    std::unordered_set<int> star;
    
    std::queue<int> Q;
    Q.push(alpha);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        if(star_index == iorder(a, STAR)) {
            
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
    }
    
    return star;
 
}


py::array_t<double> construct_time_of_insertion_map(CellComplex &comp, py::array_t<int> vertex_order, 
                                py::array_t<double> vertex_time, bool dual) {
    
        
    auto vorder = vertex_order.unchecked<1>();
    auto vtime = vertex_time.unchecked<1>();
    
    class Comparator {
        
        CellComplex &comp;
        
        
        // Use vector of tuples and then can actually use lexicographic ordering
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
        
        std::unordered_set<int> star = get_star(i, !dual, comp, target_dim);
        int star_cell = *(star.begin());
                
        for(auto alpha: star) {
                        
            lex_val[i].push_back(vorder(alpha));
            
            if((dual && vorder(alpha) < vorder(star_cell))
               || (!dual && vorder(alpha) > vorder(star_cell))) {
                
                star_cell = alpha;
            }
        }
        
        cell_to_star[i] = star_cell;
        
        // Sort from high to low
        std::sort(lex_val[i].begin(), lex_val[i].end(), std::greater<int>());        
    }
    
    std::vector<int> cells(comp.ncells);
    std::iota(cells.begin(), cells.end(), 0);
    
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


std::tuple<py::array_t<int>, py::array_t<int>>  construct_discrete_gradient(CellComplex &comp, py::array_t<double> insert_order, bool dual) {     
    
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
        std::unordered_set<int> star = get_lower_star(i, dual, comp, -1, insert_order);
        
        // record all unpaired facets/cofacets of each cell in star
        std::unordered_map<int, std::unordered_set<int> > unpaired;
        for(auto alpha: star) {
            unpaired[alpha] = std::unordered_set<int>();
            
            auto range = dual ? comp.get_cofacet_range(alpha) : comp.get_facet_range(alpha);
            for(auto it = range.first; it != range.second; it++) {
                int beta = *it;
                if(star.count(beta)) {
                    unpaired[alpha].insert(beta);
                }
            }
            
        }
        
        Comparator cmp(insert_order, dual);
        std::priority_queue<int, std::vector<int>, Comparator> PQone(cmp);
        std::priority_queue<int, std::vector<int>, Comparator> PQzero(cmp);
        
        for(auto alpha: star) {
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
                
                star.erase(alpha);
                star.erase(beta);
                unpaired.erase(alpha);
                
                for(auto gamma: star) {
                    if(unpaired[gamma].count(alpha) || unpaired[gamma].count(beta)) {
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
                
                if(star.count(alpha)) {
                    Vbuf(alpha) = alpha;
                    
                    star.erase(alpha);
                    unpaired.erase(alpha);
                    
                    for(auto gamma: star) {
                        if(unpaired[gamma].count(alpha)) {
                            unpaired[gamma].erase(alpha);

                            if(unpaired[gamma].size() == 1) {
                                PQone.push(gamma);
                            }
                        }
                    }
                    
                }
            }
            
        }
                
    }
    
    py::array_t<int> coV;
    coV.resize({Vbuf.size()});    
    std::fill(coV.mutable_data(), coV.mutable_data() + coV.size(), -1);
    
    auto coVbuf = coV.mutable_unchecked<1>();
    
    for(int i = 0; i < Vbuf.size(); i++) {
        if(Vbuf(i) != -1) {
            coVbuf(Vbuf(i)) = i;
        }
    }
    
    return std::make_tuple(V, coV);
    
}


std::vector<std::tuple<int, int, int> > traverse_flow(int s, py::array_t<int> V, 
                                                       CellComplex &comp, bool co, bool coordinated) {
    
    auto Vbuf = V.unchecked<1>();
    
    std::unordered_map<int, int> nrIn;
    std::unordered_map<int, int> k;
    
    if(coordinated) {
        std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, V, comp, co, false);
        for(auto trip: traversal) {
            int c = std::get<2>(trip);
            nrIn[c] += 1;
        }
    }
    
    std::queue<int> Q;
    Q.push(s);
    
    std::unordered_set<int> seen;
    seen.insert(s);
    
    std::vector<std::tuple<int, int, int> > traversal;
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        auto range = co ? comp.get_cofacet_range(a) : comp.get_facet_range(a);
        for(auto it = range.first; it != range.second; it++) {
            int b = *it;
            if(Vbuf(b) != -1 && Vbuf(b) != a) {
                int c = Vbuf(b);
                traversal.emplace_back(a, b, c);
                
                if(!seen.count(c) && b != c) {
                    if(coordinated) {
                        k[c] += 1;
                        
                        if(k[c] != nrIn[c]) {
                            continue;
                        }
                        
                    }
                    
                    
                    Q.push(c);
                    seen.insert(c);
                    
                }   
            }
        }
        
    }
    
    return traversal;
    
}


std::vector<std::tuple<int, int, int> > calc_morse_boundary(int s, py::array_t<int> V, 
                                                             CellComplex &comp, bool co, bool oriented) {
    
    std::unordered_map<int, int> counts;
    counts[s] = 1;
    
    std::unordered_set<int> boundary;
    
    std::unordered_map<int, int> mult;
    if(oriented) {
        mult[s] = 1;
    }
    
    std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, V, comp, co, true);
    for(auto trip: traversal) {
        int a, b, c;
        std::tie(a, b, c) = trip;
        
        counts[c] += counts[a];
        
        if(b == c) {
            boundary.insert(c);
            
            // if oriented:
//                 mult[c] = mult.get(c, 0) + mult[a] * ( -coeff(a)[b] )
            
//         elif oriented:
//             mult[c] = mult.get(c, 0) + mult[a] * ( -coeff(a)[b] * coeff(c)[b] )
        }
    }
    
    std::vector<std::tuple<int, int, int> > morse_boundary;
    for(auto c: boundary) {
        if(oriented) {
            morse_boundary.emplace_back(c, counts[c], mult[c]);
        } else {
            morse_boundary.emplace_back(c, counts[c], counts[c] % 2);
        }
    }
    
    return morse_boundary;
    
}



CellComplex construct_morse_complex(py::array_t<int> V, CellComplex &comp, bool oriented) {
    
    auto Vbuf = V.unchecked<1>();
    
    CellComplex mcomp(comp.dim, false, oriented, false);
    
    for(int s = 0; s < comp.ncells; s++) {
        if(Vbuf(s) == s) {
            std::vector<int> facets;
            std::vector<int> coeffs;
            
            std::vector<std::tuple<int, int, int> > morse_boundary = calc_morse_boundary(s, V, comp, false, oriented);
            for(auto trip: morse_boundary) {
                int c, k, m;
                std::tie(c, k, m) = trip;
                
                facets.push_back(c);
                coeffs.push_back(m);
                
                if(m > 1 || m < -1) {
                    py::print("Large Coefficient:", k, m, comp.get_dim(c), comp.get_dim(s));
                }
            }
            
            mcomp.add_cell(comp.get_dim(s), facets, s, coeffs);
        }
    }
    
    mcomp.construct_cofacets();
    mcomp.make_compressed();
    return mcomp;
    
}


std::vector<std::tuple<int, int, int> > find_connections(int s, int t, py::array_t<int> V,  py::array_t<int> coV, 
                                                          CellComplex &comp) {
    
    std::unordered_set<int> active;
    active.insert(t);
        
    std::vector<std::tuple<int, int, int> > traversal = traverse_flow(t, coV, comp, true, false);
    for(auto trip: traversal) {
        int c = std::get<2>(trip);
        active.insert(c);
    }
        
    std::vector<std::tuple<int, int, int> > connections;
    
    traversal = traverse_flow(s, V, comp, false, false);
    for(auto trip: traversal) {
        int b = std::get<1>(trip);
        if(active.count(b)) {
            connections.push_back(trip);
        }
    }
        
    return connections;
    
}


// def find_connections(s, t, V, coV, I, coI):
    
//     active = set([t])
//     for (a, b, c) in traverse_flow(t, coV, coI):
//         active.add(c)
        
//     for (a, b, c) in traverse_flow(s, V, I):
//         if b in active:
//             yield (a, b, c)
    


void simplify_morse_complex(double threshold, py::array_t<int> V, py::array_t<int> coV, CellComplex &comp,
                           py::array_t<double> insert_order, bool verbose=false) {
    
    auto Vbuf = V.mutable_unchecked<1>();
    auto coVbuf = coV.mutable_unchecked<1>();
    auto iorder = insert_order.unchecked<2>();
    
    std::vector<int> crit_cells;
    for(int s = 0; s < Vbuf.size(); s++) {
        if(Vbuf(s) == s) {
            crit_cells.push_back(s);
        }
    }
    
    std::sort(crit_cells.begin(), crit_cells.end(),
       [&iorder](int lhs, int rhs) {return iorder(lhs, LEX) < iorder(rhs, LEX);});
    
    for(int n = 1; ; n++) {
        
        
        if(verbose) {
            py::print("Pass:", n);
            py::print("Critical Cells:", crit_cells.size());
        }
        
        int n_cancel = 0;
        std::vector<std::pair<int, int> >  cancel_pairs;
        for(auto s: crit_cells) {
            if(comp.get_dim(s) == 0) {
                continue;
            }
            
            int close_alpha = -1;
            double close_alpha_time;
            std::vector<std::tuple<int, int, int> > morse_boundary = calc_morse_boundary(s, V, comp, false, comp.oriented);
            for(auto trip: morse_boundary) {
                int c, k;
                std::tie(c, k, std::ignore) = trip;
                
                if(k == 1) {
                    if(close_alpha == -1 || iorder(c, LEX) > close_alpha_time) {
                        close_alpha = c;
                        close_alpha_time = iorder(c, LEX);
                    }
                }
            }
            
            if(close_alpha == -1) {
                continue;
            }
            
            int close_beta = -1;
            double close_beta_time;
            morse_boundary = calc_morse_boundary(close_alpha, coV, comp, true, comp.oriented);
            for(auto trip: morse_boundary) {
                int c, k;
                std::tie(c, k, std::ignore) = trip;
                
                if(k == 1) {
                    if(close_beta == -1 || iorder(c, LEX) < close_beta_time) {
                        close_beta = c;
                        close_beta_time = iorder(c, LEX);
                    }
                }
            }
            
            if(s == close_beta) {
                n_cancel++;
                if(iorder(close_beta, TIME)-iorder(close_alpha, TIME) <= threshold) {
                    cancel_pairs.emplace_back(close_alpha, close_beta); 
                }
            }
                
        }
        
//         std::sort(cancel_pairs.begin(), cancel_pairs.end(),
//                   [&iorder](std::pair<int, int> lhs, std::pair<int, int> rhs) {
//                       return (iorder(lhs.second, TIME)-iorder(lhs.first, TIME)) 
//                           < (iorder(rhs.second, TIME)-iorder(rhs.first, TIME));
//                   });
        
        if(verbose) {
            py::print("Cancellable Pairs:", n_cancel);
        }
        
        for(auto pair: cancel_pairs) {
            
            // s is the higher dimensional cell
            int s = pair.second;
            int t = pair.first;
            
                                    
            std::vector<std::pair<int, int> > reverse_pairs;
            
            std::vector<std::tuple<int, int, int> > connections = find_connections(s, t, V, coV, comp);
            for(auto trip: connections) {
                int a, b;
                std::tie(a, b, std::ignore) = trip;
                reverse_pairs.emplace_back(a, b);
            }
            
            for(auto pair: reverse_pairs) {
                int a = pair.first;
                int b = pair.second;
                
                Vbuf(b) = a;
                coVbuf(a) = b;
                Vbuf(a) = -1;
                coVbuf(b) = -1;
            }
            
            crit_cells.erase(std::find(crit_cells.begin(), crit_cells.end(), s));
            crit_cells.erase(std::find(crit_cells.begin(), crit_cells.end(), t));
            
        }
        
        if(verbose) {
            py::print("Cancelled Pairs:", cancel_pairs.size());
            py::print("Remaining Critical Cells:", crit_cells.size());
        }
        
        if(cancel_pairs.empty()) {
            break;
        }
        
    } 
    
}



// def simplify_morse_complex(threshold, V, coV, comp, insert_order):
     
//     TIME = 0
//     LEX = 2
        
//     crit_cells = {s for s in range(len(V)) if V[s] == s}
        
//     n = 0
//     while True:
        
//         n += 1
                 
//         print("Pass:", n)
            
//         close_pairs = []
        
//         # print(sorted(crit_cells))
//         print("Critical Cells", len(crit_cells))
                
//         for s in crit_cells:
            
//             if comp.get_dim(s) == 0:
//                 continue
                
//             # print("s", s, insert_order[s], comp.get_dim(s))
            
//             close_alpha = None
//             close_alpha_time = None
//             for (c, k, m) in calc_morse_boundary(s, V, comp.get_facets, comp.get_coeffs, oriented=comp.oriented):

//                 if k == 1:
//                     if close_alpha is None or ((insert_order[c][TIME], insert_order[c][LEX]) > close_alpha_time):
//                         close_alpha = c
//                         close_alpha_time = (insert_order[c][TIME], insert_order[c][LEX])
                    
//             if close_alpha is None:
//                 continue
                
                                
//             close_beta = None
//             close_beta_time = None
//             for (c, k, m) in calc_morse_boundary(close_alpha, coV, comp.get_cofacets, comp.get_coeffs, oriented=comp.oriented):
//                 if k == 1:
//                     if close_beta is None or ((insert_order[c][TIME], insert_order[c][LEX]) < close_beta_time):
//                         close_beta = c
//                         close_beta_time = (insert_order[c][TIME], insert_order[c][LEX])
               
           
//             if s == close_beta:
//                 close_pairs.append((insert_order[s][TIME] - insert_order[close_alpha][TIME], (close_alpha, s)))
        
        
//         close_pairs = sorted(close_pairs)
        
//         print("Close Pairs:", len(close_pairs))
        
//         if len(close_pairs) == 0 or close_pairs[0][0] > threshold:
//             print("Cancelled Pairs:", 0)
//             break
          
//         i = 0
//         for (time, (t, s)) in close_pairs:
//             if time > threshold:
//                 break
        
//             # print(i, time, t, s, comp.get_dim(t), comp.get_dim(s))
//             i += 1
//             reverse_pairs = []
//             for (a, b, c) in find_connections(s, t, V, coV, comp.get_facets, comp.get_cofacets):
//                 reverse_pairs.append((a, b))
                
//             for (a, b) in reverse_pairs:                
//                 V[b] = a
                
//                 coV[a] = b
                
//                 V[a] = -1
//                 coV[b] = -1
                
                
//                 # if a in V:
//                 #     del V[a]
//                 # if b in coV:
//                 #     del coV[b]
                    
//             crit_cells.remove(t)
//             crit_cells.remove(s)
               
//         print("Cancelled Pairs:", i)
        
//         print("Remaining Critical Cells", len(crit_cells))
            
        
                
//     return (V, coV)


std::unordered_set<int> convert_morse_to_real_complex(std::unordered_set<int> &mfeature, py::array_t<int> V,
                                                      CellComplex &comp, bool co) {
    
    std::unordered_set<int> feature;
    
    for(auto s: mfeature) {
        feature.insert(s);
        std::vector<std::tuple<int, int, int> > traversal = traverse_flow(s, V, comp, co, false);
        for(auto trip: traversal) {
            int b, c;
            std::tie(std::ignore, b, c) = trip;
            if(b != c) {
                feature.insert(c);
            }
        }
    }
    
    return feature;
    
}


std::unordered_set<int> convert_to_pixels(std::unordered_set<int> &feature, CellComplex &comp, py::array_t<double> insert_order, bool dual) {
    
    auto iorder = insert_order.unchecked<2>();
    
    std::unordered_set<int> pixels;
    
    for(auto s: feature) {
        if(comp.get_dim(s) == 0) {
            if(!dual) {
                pixels.insert(s);
            } else {
                
                std::unordered_set<int> cofaces = get_star(s, false, comp, comp.dim);
                for(auto c: cofaces) {
                    
                    std::unordered_set<int> verts = get_star(c, true, comp, 0);
                    
                    std::unordered_set<int> ucostar_verts;
                    for(auto v: verts) {
                        if(iorder(c, STAR) == iorder(v, STAR)) {
                            ucostar_verts.insert(v);
                        }
                    }
                    
                    if(ucostar_verts.empty()) {
                        int min_vert = *(verts.begin());
                        for(auto v: verts) {
                            if(iorder(v, LEX) < iorder(min_vert, LEX)) {
                                min_vert = v;
                            }
                        }
                        
                        if(feature.count(min_vert)) {
                            pixels.insert(c);
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
                            
                            pixels.insert(c);
                            
                        } else {
                            int min_vert = *(ucostar_verts.begin());
                            for(auto v: ucostar_verts) {
                                if(iorder(v, LEX) < iorder(min_vert, LEX)) {
                                    min_vert = v;
                                }
                            }
                            
                            if(feature.count(min_vert)) {
                                pixels.insert(c);
                            }
                            
                        }
                        
                    }
                                            
                }
                
            }
        } else if(comp.get_dim(s) == 1) {
            if(!dual) {
                auto range = comp.get_facet_range(s);
                pixels.insert(range.first, range.second);
            } else {
                int ucostar = (int)iorder(s, STAR);
                pixels.insert(ucostar);
            }
            
        } else if(comp.get_dim(s) == 2) {
            
            if(!dual) {  
                auto rangei = comp.get_facet_range(s);
                for(auto iti = rangei.first; iti != rangei.second; iti++) {
                    auto rangej = comp.get_facet_range(*iti);
                    pixels.insert(rangej.first, rangej.second);
                }
                
            } else {
                pixels.insert(s);
            }
            
        }
    }
    
    
    return pixels;
    
}


std::unordered_map<int, std::unordered_set<int>> find_basins(CellComplex &mcomp, py::array_t<int> coV, CellComplex &comp, 
                                                       py::array_t<double> insert_order, bool dual) {
    
    auto iorder = insert_order.unchecked<2>();
    
    std::unordered_map<int, std::unordered_set<int> > basins;
    
    for(auto c: mcomp.get_cells()) {
        if(mcomp.get_dim(c) == 0) {
            int index = dual ? (int)iorder(c, STAR) : c;
            
            
            std::unordered_set<int> mfeature;
            mfeature.insert(c);
            
            std::unordered_set<int> feature = convert_morse_to_real_complex(mfeature, coV, comp, true);
            
            std::unordered_set<int> pixels = convert_to_pixels(feature, comp, insert_order, dual);
            
            // might ont need piecewise_construct, but it may be safer this way
            basins.emplace(std::piecewise_construct, std::forward_as_tuple(index), 
                           std::forward_as_tuple(pixels.begin(), pixels.end()));
                     
        }
    }
    
    return basins;
}


std::unordered_set<int> find_morse_skeleton(CellComplex &mcomp, py::array_t<int> V, CellComplex &comp, 
                                                       int d, py::array_t<double> insert_order, bool dual) {
    
    auto Vbuf = V.unchecked<1>();
    
    std::unordered_set<int> skeleton;
    
    for(int s = 0; s < comp.ncells; s++) {
        if(Vbuf(s) == s && mcomp.get_dim(s) == d) {
                        
            std::unordered_set<int> mfeature;
            mfeature.insert(s);
            
            std::unordered_set<int> feature = convert_morse_to_real_complex(mfeature, V, comp, false);
            
            std::unordered_set<int> pixels = convert_to_pixels(feature, comp, insert_order, dual);
            
            skeleton.insert(pixels.begin(), pixels.end());
            
        }
    }
    
    
    return skeleton;
}


std::unordered_set<int> extract_persistence_feature(int i, int j, CellComplex &mcomp, CellComplex &comp,
                                                     py::array_t<int> V, py::array_t<int> coV, 
                                                     py::array_t<double> insert_order) {
    
    auto iorder = insert_order.unchecked<2>();
    
    std::unordered_set<int> seen;
    std::queue<int> Q;
    if(comp.get_dim(i) == 0) {
        seen.insert(i);
        Q.push(i);
    } else {
        seen.insert(j);
        Q.push(j);
    }
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        auto rangea = (comp.get_dim(i) == 0) ? comp.get_cofacet_range(a) : comp.get_facet_range(a);
        for(auto ita = rangea.first; ita != rangea.second; ita++) {

            int b = *ita;
            
            if((comp.get_dim(i) == 0 && iorder(b, TIME) >= iorder(j, TIME)) 
              ||(comp.get_dim(i) != 0 && iorder(b, TIME) <= iorder(i, TIME))) {
                continue;
            }

            auto rangeb = (comp.get_dim(i) == 0) ? comp.get_facet_range(b) : comp.get_cofacet_range(b);
            for(auto itb = rangeb.first; itb != rangeb.second; itb++) {

                int c = *itb;
                if(!seen.count(c) && c != a) {
                    Q.push(c);
                    seen.insert(c);
                }
                
            }
        }
    }
    
    if(comp.get_dim(i) == 0) {
        return convert_morse_to_real_complex(seen, coV, comp, true);     
    } else {
        return convert_morse_to_real_complex(seen, V, comp, false);
    }
        
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

    
#endif // MORSE