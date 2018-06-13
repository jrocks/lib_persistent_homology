#ifndef LANDSCAPE_HPP
#define LANDSCAPE_HPP
    
#include <algorithm>
#include "math.h"
    
#include <Eigen/Core>
#include <Eigen/Geometry>
  
typedef Eigen::VectorXd XVec;
typedef Eigen::MatrixXd XMat;

typedef Eigen::Ref<const XVec > CRXVec;
typedef Eigen::Ref<const XMat > CRXMat;
    
#include <pybind11/pybind11.h>
namespace py = pybind11;



XMat calc_landscape(CRXVec birth, CRXVec death, CRXVec t, int K) {

    int NL = (K > 0) ? K : birth.size();
    
    XMat landscape = XMat::Zero(NL, t.size());
    
    XVec b = birth.cwiseMin(death);
    XVec d = birth.cwiseMax(death);
    
    for(int i = 0; i < t.size(); i++) {

        // \lambda_k(t) = max(0, min(t - b, d -t) )
        
        XVec lambdakt = ((t[i] - b.array()).min(d.array() - t[i])).max(0.0);
        
        std::sort(lambdakt.data(), lambdakt.data()+lambdakt.size(), std::greater<double>());
        
        
        landscape.block(0, i, std::min(NL, (int)lambdakt.size()), 1) = lambdakt.segment(0, std::min(NL, (int)lambdakt.size()));
        
    }
    
    
    return landscape;
    
    
    
}

double calc_landscape_norm(CRXMat landscape, CRXVec t) {
    
    XMat landscape2 = landscape.array().square();
    
    XMat dland = landscape2.block(0, 1, landscape2.rows(), landscape2.cols()-1) + landscape2.block(0, 0, landscape2.rows(), landscape2.cols()-1);
    
    XVec dt = t.segment(1, t.size()-1) - t.segment(0, t.size()-1);
    
    return sqrt(0.5 * (dland * dt).sum());
    
}


double calc_landscape_dist(CRXMat landscapei, CRXMat landscapej, CRXVec t) {
        
    XMat landscape2 = (landscapej - landscapei).array().square();
    
    XMat dland = landscape2.block(0, 1, landscape2.rows(), landscape2.cols()-1) + landscape2.block(0, 0, landscape2.rows(), landscape2.cols()-1);
    
    XVec dt = t.segment(1, t.size()-1) - t.segment(0, t.size()-1);
    
    return sqrt(0.5 * (dland * dt).sum());
    
}


double calc_landscape_norm(CRXVec birth, CRXVec death, CRXVec t) {
    
    XMat landscape2 = calc_landscape(birth, death, t, -1).array().square();
    
    XMat dland = landscape2.block(0, 1, landscape2.rows(), landscape2.cols()-1) + landscape2.block(0, 0, landscape2.rows(), landscape2.cols()-1);
    
    XVec dt = t.segment(1, t.size()-1) - t.segment(0, t.size()-1);
    
    return sqrt(0.5 * (dland * dt).sum());
    
}


double calc_landscape_dist(CRXVec birth1, CRXVec death1, CRXVec birth2, CRXVec death2, CRXVec t) {
    
    int K = std::max(birth1.size(), birth2.size());
    
    XMat landscape2 = (calc_landscape(birth2, death2, t, K) - calc_landscape(birth1, death1, t, K)).array().square();
    
    XMat dland = landscape2.block(0, 1, landscape2.rows(), landscape2.cols()-1) + landscape2.block(0, 0, landscape2.rows(), landscape2.cols()-1);
    
    XVec dt = t.segment(1, t.size()-1) - t.segment(0, t.size()-1);
    
    return sqrt(0.5 * (dland * dt).sum());
    
}



#endif // LANDSCAPE_HPP
