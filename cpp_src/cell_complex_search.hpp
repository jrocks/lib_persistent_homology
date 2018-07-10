#ifndef COMPLEXSEARCH_HPP
#define COMPLEXSEARCH_HPP
    
#include "cell_complex.hpp"
    
#include <vector>
#include <unordered_set>
#include <queue>
#include "math.h"
    
typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;
typedef Eigen::VectorXi XiVec;
typedef Eigen::MatrixXi XiMat;
typedef Eigen::Ref<XVec > RXVec;
typedef Eigen::Ref<XMat > RXMat;
typedef Eigen::Ref<XiMat > RXiMat;
    
#include <pybind11/pybind11.h>
namespace py = pybind11;



std::vector<int> find_distances(int start, CellComplex &comp, int max_dist=-1) {
    
    
    std::vector<int> dist(comp.ncells, -1);
    
    dist[start] = 0;
    
    std::unordered_set<int> seen;
    seen.insert(start);
    
    std::queue<int> Q;
    Q.push(start);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        std::unordered_set<int> cofaces = get_star(a, false, comp, -1);
        for(auto c: cofaces) {
            if(!seen.count(c)) {
                seen.insert(c);
                dist[c] = dist[a] + 1;
                 
                if(max_dist == -1 || dist[c] < max_dist) {
                    Q.push(c);
                }
            }
        }
        
        std::unordered_set<int> faces = get_star(a, true, comp, -1);
        for(auto c: faces) {
            if(!seen.count(c)) {
                seen.insert(c);
                dist[c] = dist[a] + 1;
                
                if(max_dist == -1 || dist[c] < max_dist) {
                    Q.push(c);
                }
            }
        }
    }
        
    return dist;
    
}

// RXiMat find_pairwise_distances(CellComplex &comp, int target_dim=-1, int max_cutoff=-1) {
    
//     int ncells = 0;
    
//     if(target_dim == -1) {
//         ncells = comp.ncells;
//     } else {
//         for(int i = 0; i < comp.ncells; i++) {
//             if(comp.get_dim(i) == target_dim) {
//                 ncells++;
//             }
//         }
//     }
    
//     py::print(ncells, py::arg("flush")=true);
    
    
//     // Map of all cells to all other cells with dist < max_cutoff
//     std::unordered_map<int, std::unordered_map<int, int> > dist_map;
//     for(int i = 0; i < comp.ncells; i++) {
//         dist_map[i][i] = 0;
//     }
    
//     int count = 0;
    
//     for(int i = 0; i < comp.ncells; i++) {
//         if(target_dim != -1 && comp.get_dim(i) != target_dim) {
//             continue;
//         }
        
//         if(count % 100 == 0) {
//             py::print(i, count, "/", ncells, py::arg("flush")=true);
//         }
        
//         count++;
        
//         std::vector<int> dist(comp.ncells, -1);
//         dist[i] = 0;
        
//         std::unordered_set<int> seen;
//         seen.insert(i);
        
//         std::queue<int> Q;
//         Q.push(i);
        
//         while(!Q.empty()) {
//             int a = Q.front();
//             Q.pop();

//             std::unordered_set<int> cofaces = get_star(a, false, comp, -1);
//             for(auto c: cofaces) {
//                 if(!seen.count(c)) {
//                     seen.insert(c);
//                     dist[c] = dist[a] + 1;
                    
//                     if(max_cutoff == -1 || dist[c] < max_cutoff) {
//                         Q.push(c);
//                     }
                    
//                 }
//             }

//             std::unordered_set<int> faces = get_star(a, true, comp, -1);
//             for(auto c: faces) {
//                 if(!seen.count(c)) {
//                     seen.insert(c);
//                     dist[c] = dist[a] + 1;
                    
//                     if(max_cutoff == -1 || dist[c] < max_cutoff) {
//                         Q.push(c);
//                     }
                    
//                 }
//             }
//         }
                
        
        
        
// //         std::queue<int> Q;
// //         Q.push(i);
        
// //         while(!Q.empty()) {
// //             int a = Q.front();
// //             Q.pop();

// //             std::unordered_set<int> cofaces = get_star(a, false, comp, -1);
// //             for(auto c: cofaces) {
// //                 // If already seen, then skip
// //                 if(dist_map[i].count(c)) {
// //                     continue;
// //                 }
                
// //                 // Otherwise update distance
// //                 int new_dist = dist_map[i][a] + 1;
// //                 dist_map[i][c] = new_dist;
// //                 dist_map[c][i] = new_dist;
                
// //                 // If within cutoff distance, then add to queue
// //                 if(max_cutoff == -1 || new_dist < max_cutoff) {
// //                     Q.push(c);
// //                 }

// //             }

// //             std::unordered_set<int> faces = get_star(a, true, comp, -1);
// //             for(auto c: faces) {
// //                 // If already seen, then skip
// //                 if(dist_map[i].count(c)) {
// //                     continue;
// //                 }
                
// //                 // Otherwise update distance
// //                 int new_dist = dist_map[i][a] + 1;
// //                 dist_map[i][c] = new_dist;
// //                 dist_map[c][i] = new_dist;
                
// //                 // If within cutoff distance, then add to queue
// //                 if(max_cutoff == -1 || new_dist < max_cutoff) {
// //                     Q.push(c);
// //                 }
// //             }
// //         }
        
// //         py::print("seen", seen.size());
        
        
// //         if(target_dim == -1) {
// //             for(int j = 0; j < comp.ncells; j++) {
// //                 dist_mat(i, j) = dist[j];
// //             }
// //         } else {
// //             for(int j = 0; j < comp.ncells; j++) {
// //                 if(comp.get_dim(j) == target_dim) {
// //                     dist_mat(comp.get_label(i), comp.get_label(j)) = dist[j];
// //                 }
// //             }
// //         }
  
//     }
    
//     XiMat dist_mat = XiMat::Constant(ncells, ncells, -1);
    
// //     if(target_dim == -1) {
// //         for(auto pairi: dist_map) {
// //             for(auto pairj: pairi.second) {
// //                 dist_mat(pairi.first, pairj.first) = pairj.second;
// //             }
// //         }
// //     } else {
// //         for(auto pairi: dist_map) {
// //             for(auto pairj: pairi.second) {
// //                 dist_mat(comp.get_label(pairi.first), comp.get_label(pairj.first)) = pairj.second;
// //             }
// //         }
        
// //     }
        
//     return dist_mat;
    
// }


std::tuple<std::vector<int>, std::vector<XVec > > get_neighborhood(int p, double dist, int DIM, RXVec pos, RXMat L) {
    
    XVec posi = pos.segment(DIM*p, DIM);
    
    std::vector<int> neighborhood;
    std::vector<XVec > neigh_pos;
    
    for(int i = 0; i < pos.size() / DIM; i++) {
        XVec posj = pos.segment(DIM*i, DIM);
        
        XVec bvec = posj - posi;
        
        for(int d = 0; d < DIM; d++) {
            if(std::fabs(bvec[d]) > L(d, d) / 2.0) {
                
                bvec -= ((bvec[d] > 0) - (bvec[d] < 0)) * L.col(d);
                            
            }
            
        }
        
        if(bvec.norm() < dist) {
            neighborhood.push_back(i);
            neigh_pos.push_back(posi+bvec);
        }
    }
    
    
    return std::make_tuple(neighborhood, neigh_pos);
    
    
}

std::unordered_set<int> find_thresholded_component(int start, double threshold, StarFiltration &filt, CellComplex &comp) {
    
    std::unordered_set<int> component;
    
    std::unordered_set<int> seen;
    seen.insert(start);
    
    std::queue<int> Q;
    Q.push(start);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        if(comp.get_dim(a) == 0) {
            component.insert(a);
        }
        
        std::unordered_set<int> cofaces = get_star(a, false, comp, -1);
        for(auto c: cofaces) {
            if(!seen.count(c) && filt.get_time(c) < threshold) {
                Q.push(c);
                seen.insert(c);
                
            }
        }
        
        std::unordered_set<int> faces = get_star(a, true, comp, -1);
        for(auto c: faces) {
            if(!seen.count(c) && filt.get_time(c) < threshold) {
                Q.push(c);
                seen.insert(c);
            }
        }
    }
        
    return component;
    
}


#endif // COMPLEXSEARCH_HPP