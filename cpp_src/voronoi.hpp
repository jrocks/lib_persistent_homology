#ifndef VORONOI_HPP
#define VORONOI_HPP

#include "alpha_complex.hpp"

//Although everything is templated with "Dim", algorithms are only promised to
//work in 2d.


template <int DIM>
class Voronoi {

private:
  Embedding<DIM> embed;
  CellComplex comp;
  XVec voronoi_vertices;
  XVec cell_centroids;
  XVec cell_areas;

public:
  Voronoi();
  Voronoi(int NV, RXVec vert_pos, RDMat box_mat, bool periodic);

  //construction subroutines
  void construct_voronoi_vertices();
  void construct_cell_areas_and_centroids();


  XVec get_cell_centroids() const; 
  Embedding<DIM> get_embedding() const;
};

template <int DIM>
Voronoi<DIM>::Voronoi() {}

template <int DIM>
Voronoi<DIM>::Voronoi(int NV, RXVec vert_pos, RDMat box_mat, bool periodic) {
  embed = Embedding<DIM>(NV, vert_pos, box_mat, periodic);
  comp = construct_alpha_complex(embed,XVec::Zero(NV), false);
  construct_voronoi_vertices();
  construct_cell_areas_and_centroids();
}

template <int DIM>
void Voronoi<DIM>::construct_voronoi_vertices() {
  voronoi_vertices = XVec::Zero(DIM*comp.ndcells[2]);
  for (int v = 0; v < comp.ndcells[2]; v++) {
    std::vector<int> neighb_cells = comp.get_cofaces(v+comp.dcell_range[2].first,0);
    voronoi_vertices.segment<DIM>(DIM*v) = calc_circumsphere(neighb_cells,
        embed, std::vector<double>(DIM+1,0.0))[1];
  }
}

//use the "shoelace formula" - need to compute vertices first.
template <int DIM>
void Voronoi<DIM>::construct_cell_areas_and_centroids() {
  cell_areas = XVec::Zero(comp.ndcells[0]);
  cell_centroids = XVec::Zero(DIM*comp.ndcells[0]);
  for (int c = 0; c < comp.ndcells[0]; c++) {
    std::vector<int> neighb_vertices = comp.get_cofaces(c+comp.dcell_range[0].first, 2);
    std::vector< DVec > vert_pos(neighb_vertices.size(), DVec::Zero(DIM));
    for (int i = 0; i < vert_pos.size(); i++) {vert_pos[i] = voronoi_vertices.segment<DIM>(DIM*neighb_vertices[i]);}

    //now sort vertices counterclockwise - definitely only works in 2d...
    auto angle_comp = [] (DVec x, DVec y) {return atan2(x[1], x[0]) < atan2(y[1], y[0]) ;};
    std::sort(vert_pos.begin(), vert_pos.end(), angle_comp);

    int n = vert_pos.size();
    for (int i = 0; i < n ; i++) {
      cell_areas[c] += 0.5*(vert_pos[i](0)*vert_pos[(i+1) % n](1)-vert_pos[(i+1) % n](0)*vert_pos[i](1));
    }

    for (int i = 0; i < vert_pos.size(); i++) {
      double temp = ((1.0)/(6.0*cell_areas[c]))*(vert_pos[i](0)*vert_pos[(i+1) % n](1)-vert_pos[(i+1) % n](0)*vert_pos[i](1));
      cell_centroids.segment<DIM>(DIM*c) += temp*(vert_pos[i] + vert_pos[(i+1)%n]);
    }


  }
}

template <int DIM>
XVec Voronoi<DIM>::get_cell_centroids() const {
  return cell_centroids;
}

template<int DIM>
Embedding<DIM> Voronoi<DIM>::get_embedding() const {
  return embed;
}

#endif
