#define TIME 0
#define STAR 1
#define LEX 2

#include "cell_complex.hpp"
#include "persistent_homology_algs.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
    
namespace py = pybind11;


PYBIND11_MODULE(phutil, m) {
    
    py::class_<CellComplex>(m, "CellComplex")
        .def(py::init<int, bool, bool, bool>())
        .def("add_cell", &CellComplex::add_cell)
        .def("make_compressed", &CellComplex::make_compressed)
        .def("get_dim", &CellComplex::get_dim)
        .def("get_facets", (std::vector<int>* (CellComplex::*)(int))&CellComplex::get_facets,
             py::return_value_policy::take_ownership)
        .def("get_coeffs", &CellComplex::get_coeffs, py::return_value_policy::take_ownership)
        .def("get_cofacets", (std::vector<int>* (CellComplex::*)(int))&CellComplex::get_cofacets,
             py::return_value_policy::take_ownership)
        .def("get_cells", &CellComplex::get_cells, py::return_value_policy::take_ownership)
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
    
    m.def("check_boundary_op", &check_boundary_op);
    m.def("get_star", &get_star, py::return_value_policy::take_ownership);
    m.def("get_lower_star", &get_lower_star, py::return_value_policy::take_ownership);
    m.def("construct_time_of_insertion_map", &construct_time_of_insertion_map);
    m.def("construct_discrete_gradient", &construct_discrete_gradient);
     
};