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

XMat find_pairwise_euclid_distances(std::vector<int> &particles, int DIM, RXVec pos, RXMat box_mat) {
    
    XMat dist_mat = XMat::Zero(particles.size(), particles.size());
        
    for(std::size_t i = 0; i < particles.size(); i++) {
        
        int pi = particles[i];
        
        XVec posi = pos.segment(DIM*pi, DIM);
        
        for(std::size_t j = i+1; j < particles.size(); j++) {
            
            int pj = particles[j];
            
            XVec posj = pos.segment(DIM*pj, DIM);
            
            XVec bvec = posj - posi;
            
            for(int d = 0; d < DIM; d++) {
                if(std::fabs(bvec(d)) > 0.5) {
                    bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
                }
            }
            
            dist_mat(i, j) = (box_mat*bvec).norm();
            dist_mat(j, i) = dist_mat(i, j);
        
        }
                
    }
        
    return dist_mat;
    
}

XiMat find_pairwise_tri_distances(std::vector<int> &particles, CellComplex &comp) {
    
    int NV = 0;
    for(int i = 0; i < comp.ncells; i++) {
        if(comp.get_dim(i) == 0) {
            NV++;
        }
    }
    
    XiMat full_dist_mat = XiMat::Constant(NV, NV, -1);
    
    
    for(std::size_t i = 0; i < particles.size(); i++) {
        
        int pi = particles[i];
        
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
    
    
    XiMat dist_mat = XiMat::Constant(particles.size(), particles.size(), -1);
    
    for(std::size_t i = 0; i < particles.size(); i++) {
        for(std::size_t j = 0; j < particles.size(); j++) {
            dist_mat(i, j) = full_dist_mat(particles[i], particles[j]);
            dist_mat(j, i) = dist_mat(i, j);
        }
    }
    
        
    return dist_mat;
    
}


std::vector<double> find_all_euclid_distances(int start, int NP, int DIM, RXVec pos, RXMat box_mat) {
    
    
    std::vector<double> dist(NP);
    
    XVec posi = pos.segment(DIM*start, DIM);
        
    for(int pj = 0; pj < NP; pj++) {
                        
        XVec posj = pos.segment(DIM*pj, DIM);

        XVec bvec = posj - posi;

        for(int d = 0; d < DIM; d++) {
            if(std::fabs(bvec(d)) > 0.5) {
                bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));
            }
        }

        dist[pj] = (box_mat*bvec).norm();        
                
    }
        
    return dist;
    
}

std::vector<int> find_all_tri_distances(int start, CellComplex &comp, int max_dist=-1) {
    
    
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

std::vector<std::vector<int> > find_nearest_neighbors(int start, int max_dist, int target_dim, CellComplex &comp) {
    
    std::vector<std::vector<int> > neighbors(max_dist+1);
    neighbors[0].push_back(start);
    
    std::unordered_map<int, int> dist;
    dist[start] = 0;
    
    
    std::queue<int> Q;
    Q.push(start);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        auto range = comp.get_facet_range(a);
        
        for(auto it = range.first; it!= range.second; it++) {
            int b = *it;
            
            auto corange = comp.get_cofacet_range(b);
            for(auto coit = corange.first; coit!= corange.second; coit++) {
                int c = *coit;
                
                if(!dist.count(c)) {
                    dist[c] = dist[a] + 1;
                    
                    neighbors[dist[c]].push_back(c);
                    
                    if(dist[c] < max_dist) {
                        Q.push(c);
                    }
                }
                
            }
            
        }
        
    }
    
    return neighbors;
    
    
    
}


std::unordered_map<int, int> find_search_zone_distances(int start, std::unordered_set<int> &search_zone, CellComplex &comp) {
    
    
    std::unordered_map<int, int> dist;
    dist[start] = 0;
    
    std::queue<int> Q;
    Q.push(start);
    
    while(!Q.empty()) {
        int a = Q.front();
        Q.pop();
        
        std::unordered_set<int> cofaces = get_star(a, false, comp, -1);
        for(auto c: cofaces) {
            if(search_zone.count(c) && !dist.count(c)) {
                dist[c] = dist[a] + 1;
                Q.push(c);
            }
        }
        
        std::unordered_set<int> faces = get_star(a, true, comp, -1);
        for(auto c: faces) {
            if(search_zone.count(c) && !dist.count(c)) {
                dist[c] = dist[a] + 1;
                Q.push(c);
            }
        }
    }
        
    return dist;
    
}



std::tuple<std::vector<int>, std::vector<XVec > > get_neighborhood(int p, double dist, int DIM, RXVec pos, RXMat box_mat) {
    
    XVec posi = pos.segment(DIM*p, DIM);
    
    std::vector<int> neighborhood;
    std::vector<XVec > neigh_pos;
    
    for(int i = 0; i < pos.size() / DIM; i++) {
        XVec posj = pos.segment(DIM*i, DIM);
        
        XVec bvec = posj - posi;
        
        for(int d = 0; d < DIM; d++) {
            if(std::fabs(bvec(d)) > 0.5) {
                bvec(d) -= ((bvec(d) > 0) - (bvec(d) < 0));          
            }
        }
                
        if((box_mat * bvec).norm() < dist) {
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