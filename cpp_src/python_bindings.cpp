#define TIME 0
#define STAR 1
#define LEX 2

#include "cell_complex.hpp"
#include "cubical_cell_complex.hpp"
#include "graph_cell_complex.hpp"
#include "filtration.hpp"
#include "morse.hpp"
#include "persistent_homology.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
    
namespace py = pybind11;


PYBIND11_MODULE(chomology, m) {
    
    py::class_<CellComplex>(m, "CellComplex")
        .def(py::init<int, bool, bool>(), py::arg("dim"), py::arg("regular")=true, py::arg("oriented")=false)
        .def("add_cell", (void (CellComplex::*)(int, int, std::vector<int>&, std::vector<int>&)) &CellComplex::add_cell)
        .def("make_compressed", &CellComplex::make_compressed)
        .def("get_label", &CellComplex::get_label)
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
    
    m.def("construct_graph_complex", &construct_graph_complex);
    m.def("find_corners", &find_corners);
    
    py::class_<Filtration>(m, "Filtration")
        .def_readonly("ncells", &Filtration::ncells)
        .def_readonly("nprimary_cells", &Filtration::nprimary_cells)
        .def(py::init<int, std::vector<double>& >())
        .def("set_primary_order", &Filtration::set_primary_order)
        .def("set_filtration_order", &Filtration::set_filtration_order)
        .def("set_star", &Filtration::set_star)
        .def("get_time", &Filtration::get_time)
        .def("get_primary_order", &Filtration::get_primary_order)
        .def("get_filtration_order", &Filtration::get_filtration_order)
        .def("get_star", &Filtration::get_star)
        .def("get_filtration", &Filtration::get_filtration);
        
    
    m.def("perform_watershed_transform", &perform_watershed_transform);
    m.def("construct_filtration", &construct_filtration);
        
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
    m.def("convert_morse_to_real_complex", &convert_morse_to_real_complex);
    m.def("convert_to_pixels", &convert_to_pixels);
    m.def("get_boundary_pixels", &get_boundary_pixels);
    m.def("find_basins", &find_basins);
    m.def("find_morse_skeleton", &find_morse_skeleton);
    m.def("extract_persistence_feature", &extract_persistence_feature);
    m.def("get_boundary", &get_boundary);
     
};