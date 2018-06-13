#ifndef UTILALGS_HPP
#define UTILALGS_HPP
    
#include "cell_complex.hpp"
    
#include <vector>
#include <unordered_set>
#include <queue>
#include "math.h"
    
typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;
typedef Eigen::Ref<XVec > RXVec;
typedef Eigen::Ref<XMat > RXMat;
    
#include <pybind11/pybind11.h>
namespace py = pybind11;



std::vector<int> find_distances(int start, CellComplex &comp) {
    
    
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
                Q.push(c);
                seen.insert(c);
                dist[c] = dist[a] + 1;
            }
        }
        
        std::unordered_set<int> faces = get_star(a, true, comp, -1);
        for(auto c: faces) {
            if(!seen.count(c)) {
                Q.push(c);
                seen.insert(c);
                dist[c] = dist[a] + 1;
            }
        }
    }
        
    return dist;
    
}


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


#endif // UTILALGS_HPP