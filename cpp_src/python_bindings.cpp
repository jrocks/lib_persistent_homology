#include "cell_complex.hpp"
#include "cubical_cell_complex.hpp"
#include "graph_cell_complex.hpp"
#include "filtration.hpp"
#include "morse.hpp"
#include "persistent_homology.hpp"

#ifdef DELAUNAY
    #include "delaunay_cell_complex.hpp"
#endif   


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
    
namespace py = pybind11;


template <int DIM> void init(py::module &m) {
    
    m.def((std::string("construct_delaunay_complex_")+std::to_string(DIM)+std::string("D")).c_str(), &construct_delaunay_complex<DIM>);
    m.def((std::string("calc_alpha_vals_")+std::to_string(DIM)+std::string("D")).c_str(), &calc_alpha_vals<DIM>);
    
}


PYBIND11_MODULE(chomology, m) {
    
    // init<1>(m);
    init<2>(m);
    // init<3>(m);
    
    
    py::class_<CellComplex>(m, "CellComplex")
        .def(py::init<int, bool, bool>(), py::arg("dim"), py::arg("regular")=true, py::arg("oriented")=false)
        .def("add_cell", (void (CellComplex::*)(int, int, std::vector<int>&, std::vector<int>&)) &CellComplex::add_cell)
        .def("make_compressed", &CellComplex::make_compressed)
        .def("get_label", &CellComplex::get_label)
        .def("get_labels", &CellComplex::get_labels)
        .def("get_dim", &CellComplex::get_dim)
        .def("get_facets", &CellComplex::get_facets, py::return_value_policy::take_ownership)
        .def("get_coeffs", &CellComplex::get_coeffs, py::return_value_policy::take_ownership)
        .def("get_cofacets", &CellComplex::get_cofacets, py::return_value_policy::take_ownership)
        .def("construct_cofacets", &CellComplex::construct_cofacets)
        .def_readonly("dim", &CellComplex::dim)
        .def_readonly("ncells", &CellComplex::ncells)
        .def_readonly("regular", &CellComplex::regular)
        .def_readonly("oriented", &CellComplex::oriented)
        .def_readonly("dims", &CellComplex::dims)
        .def_readonly("facet_ind", &CellComplex::facet_ind)
        .def_readonly("facets", &CellComplex::facets)
        .def_readonly("coeffs", &CellComplex::coeffs)
        .def_readonly("cofacet_ind", &CellComplex::cofacet_ind)
        .def_readonly("cofacets", &CellComplex::cofacets);
    
    
    m.def("construct_cubical_complex", &construct_cubical_complex, "Constructs a cublical complex.");
    m.def("construct_masked_cubical_complex", &construct_masked_cubical_complex);
    m.def("check_boundary_op", &check_boundary_op, 
          "Checks the boundary operator of a complex to ensure that \\delta_d\\delta_(d-1) = 0 for each cell.");
    m.def("get_boundary_pixels", &get_boundary_pixels);
    
    py::class_<Graph>(m, "Graph")
        .def(py::init<int, int, std::vector<int>&, std::vector<int>&>())
        .def("set_embedding", &Graph::set_embedding);
    
    m.def("construct_graph_complex", &construct_graph_complex);
    m.def("find_corners", &find_corners);
    m.def("calc_corner_strains", &calc_corner_strains);
    m.def("calc_edge_extension", &calc_edge_extension);
    m.def("perform_corner_transform", &perform_corner_transform);
    m.def("construct_corner_complex", &construct_corner_complex);
    
        
    py::class_<Filtration>(m, "Filtration")
        .def_readonly("ncells", &Filtration::ncells)
        .def(py::init<std::vector<double>&, CellComplex&, bool>(),
            py::arg("time"), py::arg("comp"), py::arg("ascend")=true)
        .def("get_time", &Filtration::get_time)
        .def("set_total_order", &Filtration::set_total_order)
        .def("get_total_order", &Filtration::get_total_order)
        .def("get_filtration", &Filtration::get_filtration);
    
    py::class_<StarFiltration>(m, "StarFiltration")
        .def_readonly("ncells", &StarFiltration::ncells)
        .def_readonly("fdim", &StarFiltration::fdim)
        .def_readonly("nsteps", &StarFiltration::nsteps)
        .def(py::init<std::vector<double>&, std::vector<int>&, CellComplex&, bool, bool>(),
            py::arg("time"), py::arg("subcomplex_order"), py::arg("comp"), py::arg("ascend")=true, py::arg("co")=false)
        .def("get_filt_cell", &StarFiltration::get_filt_cell)
        .def("get_time", &StarFiltration::get_time)
        .def("get_subcomplex_order", &StarFiltration::get_subcomplex_order)
        .def("add_to_subcomplex", &StarFiltration::add_to_subcomplex)
        .def("set_total_order", &StarFiltration::set_total_order)
        .def("get_total_order", &StarFiltration::get_total_order)
        .def("get_filtration", &StarFiltration::get_filtration);
        
    
    m.def("construct_filtration", &construct_filtration,
         py::arg("time"), py::arg("comp"), py::arg("ascend")=true);
    
    m.def("perform_watershed_transform", &perform_watershed_transform, 
          py::arg("time"), py::arg("comp"), py::arg("ascend")=true, py::arg("co")=false);
    m.def("construct_star_filtration", &construct_star_filtration,
         py::arg("time"), py::arg("subcomplex_order"), py::arg("comp"), py::arg("ascend")=true, py::arg("co")=false);
    m.def("reduce_filtration", &reduce_filtration);
        
    m.def("get_star", &get_star,
         "Finds the star/costar of a give cell.");
    m.def("get_lower_star", &get_lower_star,
         "Finds the lower star/upper costar of a cell with a given Morse function defined by a complete cell filtration order.");
    
    
    m.def("construct_discrete_gradient", &construct_discrete_gradient,
         "Constructs the discrete gradient of a Morse function given a complete cell filtration order.");
    m.def("traverse_flow", &traverse_flow);
    m.def("calc_morse_boundary", &calc_morse_boundary);
    m.def("construct_morse_complex", &construct_morse_complex);
    m.def("find_connections", &find_connections);
    m.def("simplify_morse_complex", &simplify_morse_complex, 
          py::arg("threshold"), py::arg("V"), py::arg("coV"), py::arg("comp"),
          py::arg("insert_order"), py::arg("verbose") = false);
    m.def("convert_morse_to_real", &convert_morse_to_real, 
          py::arg("mfeature"), py::arg("V"), py::arg("coV"), py::arg("comp"), py::arg("follow_bounds")=true);
    m.def("change_feature_dim", &change_feature_dim);
    m.def("find_morse_basins", &find_morse_basins);
    m.def("find_morse_skeleton", &find_morse_skeleton);
    m.def("extract_persistence_feature", &extract_persistence_feature);
    m.def("extract_morse_feature", &extract_morse_feature);
    m.def("get_boundary", &get_boundary);
    
    m.def("calc_extended_persistence", &calc_extended_persistence);
     
};