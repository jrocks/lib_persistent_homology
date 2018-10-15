#ifndef EIGEN_MACROS_HPP
#define EIGEN_MACROS_HPP
  
// Variable length vectors and and matrices
typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;

// Reference type variable length vectors and and matrices
// Useful for pybind interfacing
typedef Eigen::Ref<XVec > RXVec;
typedef Eigen::Map<XVec > XVecMap;
    
// DIM-dimensional fixed length vectors and matrices
# define DEIGEN(DIM) \
    typedef Eigen::Matrix<double, DIM, 1> DVec; \
    typedef Eigen::Matrix<double, DIM, DIM> DMat;
    
#endif