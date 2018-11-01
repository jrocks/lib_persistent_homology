#ifndef SEARCH_HPP
#define SEARCH_HPP

#include "eigen_macros.hpp"
#include "cell_complex.hpp"
#include "embedding.hpp"
    
#include <vector>
#include <unordered_set>
#include <queue>
#include "math.h"

    
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Euclidean distances between all pairs
template <int DIM> XMat calc_euclid_pair_dists(std::vector<int> &verts, Embedding<DIM> &embed) {
    
    XMat dist_mat = XMat::Zero(verts.size(), verts.size());
        
    for(std::size_t i = 0; i < verts.size(); i++) {
        
        int pi = verts[i];
        
        DVec posi = embed.get_vpos(pi);
                
        for(std::size_t j = i+1; j < verts.size(); j++) {
            
            int pj = verts[j];
            
            DVec posj = embed.get_vpos(pj);
            
            DVec bvec = posj - posi;
            
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
            
            dist_mat(i, j) = (embed.box_mat*bvec).norm();
            dist_mat(j, i) = dist_mat(i, j);
        
        }
                
    }
        
    return dist_mat;
    
}


// Euclidean distances from point
template <int DIM> XVec calc_euclid_point_dists(int p, int NP, Embedding<DIM> &embed) {
    
    
    XVec dist = XVec::Zero(NP);
    
    DVec posi = embed.get_vpos(p);
        
    for(int pj = 0; pj < NP; pj++) {
                        
        XVec posj = embed.get_vpos(pj);

        XVec bvec = posj - posi;

        for(int d = 0; d < DIM; d++) {
            if(std::fabs(bvec(d)) > 0.5) {
                bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
            }
        }

        dist(pj) = (embed.box_mat*bvec).norm();        
                
    }
        
    return dist;
    
}

// Find neighborhood of points within distance of another point
template <int DIM> std::tuple<std::vector<int>, std::vector<DVec > > get_point_neighborhood(int p, double max_dist, int NP, Embedding<DIM> &embed) {
    
    
    XVec posi = embed.get_vpos(p);
    
    std::vector<int> neighborhood;
    std::vector<DVec > neigh_pos;
    
    for(int pj = 0; pj < NP; pj++) {
        XVec posj = embed.get_vpos(pj);

        XVec bvec = posj - posi;

        for(int d = 0; d < DIM; d++) {
            if(std::fabs(bvec(d)) > 0.5) {
                bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
            }
        }
                
        if((embed.box_mat*bvec).norm() < max_dist) {
            neighborhood.push_back(pj);
            neigh_pos.push_back(posi+bvec);
        }
    }
    
    
    return std::make_tuple(neighborhood, neigh_pos);
    
    
}

// Discrete distances between all pairs of vertices within cell complex
// Uses breadth-first search
XiMat calc_comp_pair_dists(std::vector<int> &verts, CellComplex &comp) {
    
    int NV = 0;
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) == 0) {
            NV++;
        }
    }
    
    XiMat full_dist_mat = XiMat::Constant(NV, NV, -1);
    
    
    for(std::size_t i = 0; i < verts.size(); i++) {
        
        int pi = verts[i];
        
        std::queue<int> Q;
        Q.push(pi);
        
        full_dist_mat(pi, pi) = 0;
        
        std::unordered_set<int> seen;
        seen.insert(pi);
        
        while(!Q.empty()) {
            int a = Q.front();
            Q.pop();
            
            auto corange = comp.get_cofacet_range(a);
            for(auto coit = corange.first; coit != corange.second; coit++) {
                
                auto range = comp.get_facet_range(*coit);
                for(auto it = range.first; it != range.second; it++) {
                    
                    int pj = *it;
     
                    if(!seen.count(pj)) {
                        seen.insert(pj);
                        Q.push(pj);
                    }
                    
                    if(full_dist_mat(pi, pj) == -1) {
                        
                        full_dist_mat(pi, pj) = full_dist_mat(pi, a) + 1;
                        full_dist_mat(pj, pi) = full_dist_mat(pi, pj);
                    }
                    
                }

            }
            
        }
        
    }
    
    
    XiMat dist_mat = XiMat::Constant(verts.size(), verts.size(), -1);
    
    for(std::size_t i = 0; i < verts.size(); i++) {
        for(std::size_t j = 0; j < verts.size(); j++) {
            dist_mat(i, j) = full_dist_mat(verts[i], verts[j]);
            dist_mat(j, i) = dist_mat(i, j);
        }
    }
    
        
    return dist_mat;
    
}

// Discrete distances from vertex within cell complex
// Uses breadth-first search
std::vector<int> calc_comp_point_dists(int p, CellComplex &comp, int max_dist=-1) {
    
    
    std::vector<int> dist(comp.ncells, -1);
    
    dist[p] = 0;
    
    std::unordered_set<int> seen;
    seen.insert(p);
    
    std::queue<int> Q;
    Q.push(p);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        std::unordered_set<int> cofaces = comp.get_cofaces(a);
        for(auto c: cofaces) {
            if(!seen.count(c)) {
                seen.insert(c);
                dist[c] = dist[a] + 1;
                 
                if(max_dist == -1 || dist[c] < max_dist) {
                    Q.push(c);
                }
            }
        }
        
        std::unordered_set<int> faces = comp.get_faces(a);
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

// Calculate discrete distance between pair of cells
int calc_cell_pair_dist(int ci, int cj, CellComplex &comp) {
    
    
    
    std::vector<int> dist(comp.ncells, -1);
    
    dist[ci] = 0;
    
    std::unordered_set<int> seen;
    seen.insert(ci);
    
    std::queue<int> Q;
    Q.push(ci);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        std::unordered_set<int> cofaces = comp.get_cofaces(a);
        for(auto c: cofaces) {
            if(!seen.count(c)) {
                seen.insert(c);
                dist[c] = dist[a] + 1;

                if(c == cj) {
                    return dist[c];
                }

                Q.push(c);
            }
        }
        
        std::unordered_set<int> faces = comp.get_faces(a);
        for(auto c: faces) {
            if(!seen.count(c)) {
                seen.insert(c);
                dist[c] = dist[a] + 1;
                                    
                if(c == cj) {
                    return dist[c];
                }

                Q.push(c);
            }
        }
    }
        
    return dist[cj];
    
    
    
}

// Discrete distance from vertex within cell complex
// Limits to cells in search zone
std::unordered_map<int, int> calc_comp_point_dists_search_zone(int p, std::unordered_set<int> &search_zone, CellComplex &comp) {
    
    
    std::unordered_map<int, int> dist;
    dist[p] = 0;
    
    std::queue<int> Q;
    Q.push(p);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        std::unordered_set<int> cofaces = comp.get_cofaces(a);
        for(auto c: cofaces) {
            if(search_zone.count(c) && !dist.count(c)) {
                dist[c] = dist[a] + 1;
                Q.push(c);
            }
        }
        
        std::unordered_set<int> faces = comp.get_faces(a);
        for(auto c: faces) {
            if(search_zone.count(c) && !dist.count(c)) {
                dist[c] = dist[a] + 1;
                Q.push(c);
            }
        }
    }
        
    return dist;
    
}

// Find cells within discrete distance max_dist to vertex p within cell complex
std::vector<std::vector<int> > find_nearest_neighbors(int p, int max_dist, int target_dim, CellComplex &comp) {
    
    
    std::vector<std::vector<int> > neighbors(max_dist+1);
    
    auto dist = calc_comp_point_dists(p, comp, max_dist);
    
    for(int c = 0; c < comp.ncells; c++) {
        if(dist[c] != -1) {
            neighbors[dist[c]].push_back(c);
        }
    }
    
    
    return neighbors;
    
    
    
}


// std::unordered_set<int> find_thresholded_component(int start, double threshold, StarFiltration &filt, CellComplex &comp) {
    
//     std::unordered_set<int> component;
    
//     std::unordered_set<int> seen;
//     seen.insert(start);
    
//     std::queue<int> Q;
//     Q.push(start);
    
//     while(!Q.empty()) {
//         int a = Q.front();
//         Q.pop();
        
//         if(comp.get_dim(a) == 0) {
//             component.insert(a);
//         }
        
//         std::unordered_set<int> cofaces = get_star(a, false, comp, -1);
//         for(auto c: cofaces) {
//             if(!seen.count(c) && filt.get_time(c) < threshold) {
//                 Q.push(c);
//                 seen.insert(c);
                
//             }
//         }
        
//         std::unordered_set<int> faces = get_star(a, true, comp, -1);
//         for(auto c: faces) {
//             if(!seen.count(c) && filt.get_time(c) < threshold) {
//                 Q.push(c);
//                 seen.insert(c);
//             }
//         }
//     }
        
//     return component;
    
// }


#endif // SEARCH_HPP