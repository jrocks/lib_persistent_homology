#ifndef EMBEDDING_HPP
#define EMBEDDING_HPP

#include "eigen_macros.hpp"
  
template <int DIM> class Embedding {
 
public:
    
    // Embedding informtaion for a collection of points
        
    // dimension
    const int dim;
    
protected:
    // vertex positions (undeformed and in unit box)
    XVec vert_pos;
    
public:
    // box matrix (multiply by pos to get actual embedded position)
    DMat box_mat;
    // boundary conditions
    bool periodic;
    
    Embedding(RXVec vert_pos, RDMat box_mat, bool periodic) :
        dim(DIM), vert_pos(vert_pos), box_mat(box_mat), periodic(periodic) {}
    
    // Get virtual position
    inline DVec get_vpos(int vi) {
        return vert_pos.segment<DIM>(DIM*vi);
    }
    
    // Get real position
    inline DVec get_pos(int vi) {
        return box_mat * vert_pos.segment<DIM>(DIM*vi);
    };
    
};

    
#endif // EMBEDDING_HPP