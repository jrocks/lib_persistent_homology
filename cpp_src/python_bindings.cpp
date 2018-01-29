#define TIME 0
#define STAR 1
#define LEX 2

#include "cell_complex.hpp"
#include "filtration.hpp"
#include "morse.hpp"
#include "persistent_homology.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
    
namespace py = pybind11;


PYBIND11_MODULE(phutil, m) {
    
    py::class_<CellComplex>(m, "CellComplex")
        .def(py::init<int, bool, bool, bool>())
        .def("add_cell", &CellComplex::add_cell)
        .def("make_compressed", &CellComplex::make_compressed)
        .def("get_dim", &CellComplex::get_dim)
        .def("get_facets", &CellComplex::get_facets, py::return_value_policy::take_ownership)
        .def("get_coeffs", &CellComplex::get_coeffs, py::return_value_policy::take_ownership)
        .def("get_cofacets", &CellComplex::get_cofacets, py::return_value_policy::take_ownership)
        .def("get_cells", &CellComplex::get_cells)
        .def("construct_cofacets", &CellComplex::construct_cofacets)
        .def_readonly("dim", &CellComplex::dim)
        .def_readonly("ncells", &CellComplex::ncells)
        .def_readonly("ordered", &CellComplex::ordered)
        .def_readonly("regular", &CellComplex::regular)
        .def_readonly("oriented", &CellComplex::oriented)
        .def_readonly("dims", &CellComplex::dims)
        .def_readonly("facet_ind", &CellComplex::facet_ind)
        .def_readonly("facets", &CellComplex::facets)
        .def_readonly("coeffs", &CellComplex::coeffs)
        .def_readonly("cofacet_ind", &CellComplex::cofacet_ind)
        .def_readonly("cofacets", &CellComplex::cofacets);
    
    m.def("construct_cubical_complex", &construct_cubical_complex, "Constructs a cublical complex.",
          py::return_value_policy::take_ownership);
    
    m.def("check_boundary_op", &check_boundary_op, 
          "Checks the boundary operator of a complex to ensure that \\delta_d\\delta_(d-1) = 0 for each cell.");
    
    
    
     m.def("get_star", &get_star,
         "Finds the star/costar of a give cell.");
    m.def("get_lower_star", &get_lower_star,
         "Finds the lower star/upper costar of a cell with a given Morse function defined by a complete cell filtration order.");
    
    
    
    m.def("construct_vertex_filtration_order", &construct_vertex_filtration_order,
         "Finds a valid ordering for cells with explicitly defined Morse function values in a cell complex.");
    m.def("construct_time_of_insertion_map", &construct_time_of_insertion_map,
         "Constructs a complete cell filtration order.");
    
    
    
    
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
    m.def("find_basins", &find_basins);
    m.def("find_morse_skeleton", &find_morse_skeleton);
    m.def("extract_persistence_feature", &extract_persistence_feature);
    m.def("get_boundary", &get_boundary);
     
};