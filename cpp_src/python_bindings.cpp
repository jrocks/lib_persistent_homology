#include "cell_complex.hpp"
#include "cubical_cell_complex.hpp"
#include "filtration.hpp"
#include "morse.hpp"
#include "persistent_homology.hpp"
#include "persistence_landscape.hpp"
#include "cell_complex_search.hpp"

#ifdef GRAPH
    #include "graph_cell_complex.hpp"
#endif
    
#ifdef ALPHA
    #include "alpha_cell_complex.hpp"
#endif  
    
#ifdef OPTIMAL
    #include "optimal.hpp"
#endif


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/eigen.h>
    
namespace py = pybind11;

#ifdef GRAPH
template <int DIM> void init_graph_templates(py::module &m) {
    
    py::class_<Embedding<DIM> >(m, (std::string("Embedding")+std::to_string(DIM)+std::string("D")).c_str())
        .def(py::init<RXVec, RXVec>());
    m.def((std::string("find_corners_")+std::to_string(DIM)+std::string("D")).c_str(), &find_corners<DIM>);
    m.def((std::string("calc_corner_strains_")+std::to_string(DIM)+std::string("D")).c_str(), &calc_corner_strains<DIM>);
    m.def((std::string("calc_edge_extension_")+std::to_string(DIM)+std::string("D")).c_str(), &calc_edge_extension<DIM>, 
         py::arg("disp"), py::arg("graph"), py::arg("embed"), py::arg("is_strain") = false);
    
    m.def((std::string("construct_corner_complex_")+std::to_string(DIM)+std::string("D")).c_str(), &construct_corner_complex<DIM>);
        
}
#endif

#ifdef ALPHA
template <int DIM> void init_alpha_templates(py::module &m) {
    
    
    m.def((std::string("construct_alpha_complex_")+std::to_string(DIM)+std::string("D")).c_str(), &construct_alpha_complex<DIM>,
         py::arg("NV"), py::arg("vert_pos"), py::arg("weights"), py::arg("box_mat"), py::arg("periodic")=false);
    
     m.def((std::string("calc_alpha_vals_")+std::to_string(DIM)+std::string("D")).c_str(), &calc_alpha_vals<DIM>,
         py::arg("vert_pos"), py::arg("weights"), py::arg("comp"), py::arg("box_mat"), py::arg("periodic")=false, py::arg("alpha0")=-1.0);
    
    
    m.def((std::string("calc_strains_")+std::to_string(DIM)+std::string("D")).c_str(), &calc_strains<DIM>,
         py::arg("disp"), py::arg("vert_pos"), py::arg("comp"), py::arg("box_mat"));
    
    
    m.def((std::string("calc_voronoi_D2min_")+std::to_string(DIM)+std::string("D")).c_str(), &calc_voronoi_D2min<DIM>,
         py::arg("NP"), py::arg("disp"), py::arg("vert_pos"), py::arg("comp"), py::arg("box_mat"));
}
#endif


PYBIND11_MODULE(chomology, m) {
    
    
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
    m.def("calc_elongation", &calc_elongation);
    
#ifdef GRAPH
    
    init_graph_templates<1>(m);
    init_graph_templates<2>(m);
    

    py::class_<Graph>(m, "Graph")
        .def(py::init<int, int, std::vector<int>&, std::vector<int>&>())
        .def_readonly("NV", &Graph::NV)
        .def_readonly("NE", &Graph::NE);
    
    m.def("construct_graph_complex", &construct_graph_complex);
    // m.def("perform_corner_transform", &perform_corner_transform);
#endif
    
    
#ifdef ALPHA
    
    init_alpha_templates<2>(m);
    init_alpha_templates<3>(m);
    // init_alpha_templates<4>(m);
    
    m.def("calc_radial_gap_distribution", &calc_radial_gap_distribution,
          py::arg("cell_list"), py::arg("alpha_vals"), py::arg("comp"), py::arg("max_dist")=-1, py::arg("verbose")=false);
    
    m.def("calc_angular_gap_distribution", &calc_angular_gap_distribution,
          py::arg("cell_list"), py::arg("alpha_vals"), py::arg("comp"), py::arg("max_dist")=-1, py::arg("verbose")=false);
    
    
    m.def("calc_radial_cycle_distribution", &calc_radial_cycle_distribution,
          py::arg("cell_list"), py::arg("alpha_vals"), py::arg("comp"), py::arg("vertex_types"), py::arg("max_dist")=-1, py::arg("verbose")=false);

#endif
    
    
#ifdef OPTIMAL
    m.def("calc_optimal_cycles", &calc_optimal_cycles, 
          py::arg("filt"), py::arg("comp"), py::arg("weights"), py::arg("dim")=-1, py::arg("verbose")=false,
         py::call_guard<py::scoped_ostream_redirect,
                     py::scoped_estream_redirect>());
    
//     m.def("calc_optimal_homologous_cycles", &calc_optimal_homologous_cycles, 
//           py::arg("filt"), py::arg("comp"), py::arg("weights"), py::arg("dim")=-1, py::arg("verbose")=false,
//          py::call_guard<py::scoped_ostream_redirect,
//                      py::scoped_estream_redirect>());
#endif
        
    py::class_<Filtration>(m, "Filtration")
        .def_readonly("ncells", &Filtration::ncells)
        .def(py::init<std::vector<double>&, CellComplex&, bool>(),
            py::arg("time"), py::arg("comp"), py::arg("ascend")=true)
        .def_readonly("time", &Filtration::time)
        .def("get_time", &Filtration::get_time)
        .def("set_total_order", &Filtration::set_total_order)
        .def("get_total_order", &Filtration::get_total_order)
        .def("get_filtration", &Filtration::get_filtration);
    
    py::class_<StarFiltration, Filtration>(m, "StarFiltration")
        .def_readonly("ncells", &StarFiltration::ncells)
        .def_readonly("fdim", &StarFiltration::fdim)
        .def_readonly("nsteps", &StarFiltration::nsteps)
        .def_readonly("time", &Filtration::time)
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
    m.def("convert_morse_to_real", &convert_morse_to_real, 
          py::arg("mfeature"), py::arg("V"), py::arg("coV"), py::arg("comp"), py::arg("follow_bounds")=true);
    m.def("change_feature_dim", &change_feature_dim);
    m.def("extract_morse_basin", &extract_morse_basin, 
         py::arg("s"), py::arg("mcomp"), py::arg("V"), py::arg("coV"), 
          py::arg("comp"), py::arg("filt"), py::arg("target_dim")=-1);
    m.def("find_morse_basins", &find_morse_basins,
         py::arg("mcomp"), py::arg("V"), py::arg("coV"), py::arg("filt"), py::arg("comp"), py::arg("target_dim")=-1);
    m.def("find_morse_basin_borders", &find_morse_basin_borders,
         py::arg("mcomp"), py::arg("V"), py::arg("coV"), py::arg("filt"), py::arg("comp"), py::arg("target_dim")=-1);
    m.def("find_morse_skeleton", &find_morse_skeleton);
    m.def("extract_morse_feature", &extract_morse_feature, 
         py::arg("i"), py::arg("j"), py::arg("mcomp"), py::arg("filt"), py::arg("target_dim")=-1, py::arg("complement")=false);
    m.def("extract_morse_feature_to_real", &extract_morse_feature_to_real, 
         py::arg("i"), py::arg("j"), py::arg("mcomp"), py::arg("V"), py::arg("coV"), 
          py::arg("comp"), py::arg("filt"), py::arg("complement")=false, py::arg("target_dim")=-1);
    
    
    m.def("simplify_morse_complex", &simplify_morse_complex, 
          py::arg("threshold"), py::arg("V"), py::arg("coV"), py::arg("comp"),
          py::arg("insert_order"), py::arg("leq") = true, py::arg("verbose") = false);
    m.def("find_cancel_threshold", &find_cancel_threshold);
    m.def("find_join_threshold", &find_join_threshold);
    
    m.def("get_boundary", &get_boundary);
    
    m.def("calc_extended_persistence", (std::tuple<std::tuple<std::vector<std::pair<int, int> >, 
                                        std::vector<std::pair<int, int> >, 
                                        std::vector<std::pair<int, int> > >,
                                        std::unordered_map<int, std::vector<int> > > (*)
                                        (Filtration&, Filtration&, CellComplex&, bool, int)) &calc_extended_persistence,
                                     py::arg("filt_asc"), py::arg("filt_desc"), py::arg("comp"),
                                     py::arg("ext_cycles"), py::arg("dim")=-1);
    
    m.def("calc_extended_persistence", (std::tuple<std::vector<std::pair<int, int> >, 
                                        std::vector<std::pair<int, int> >, 
                                        std::vector<std::pair<int, int> > > (*)
                                        (Filtration&, Filtration&, CellComplex&)) &calc_extended_persistence,
                                     py::arg("filt_asc"), py::arg("filt_desc"), py::arg("comp"));
    
    
    m.def("calc_birth_cycles", &calc_birth_cycles, py::arg("filt"), py::arg("comp"), py::arg("dim")=-1);
    m.def("calc_homologous_birth_cycles", &calc_homologous_birth_cycles, py::arg("filt"), py::arg("comp"), py::arg("dim")=-1);
    
    m.def("extract_persistence_feature", &extract_persistence_feature, 
         py::arg("i"), py::arg("j"), py::arg("comp"), py::arg("filt"), py::arg("target_dim")=-1, py::arg("complement")=false);
    
    // m.def("calc_persistence_landscape", (std::vector<std::vector<double> > (*)(std::vector<double>&, std::vector<double>&, std::vector<double>&, int)) &calc_persistence_landscape);
    // m.def("calc_persistence_landscape", (std::vector<std::vector<double> > (*)(std::vector<std::pair<int, int> >&, std::vector<double>&, int, StarFiltration&)) &calc_persistence_landscape);
    
    m.def("calc_landscape", &calc_landscape, py::return_value_policy::take_ownership);
    m.def("calc_landscape_norm", (double (*)(CRXMat, CRXVec))&calc_landscape_norm);
    m.def("calc_landscape_dist", (double (*)(CRXMat, CRXMat, CRXVec))&calc_landscape_dist);
    m.def("calc_landscape_norm", (double (*)(CRXVec, CRXVec, CRXVec))&calc_landscape_norm);
    m.def("calc_landscape_dist", (double (*)(CRXVec, CRXVec, CRXVec, CRXVec, CRXVec))&calc_landscape_dist);
    
    m.def("find_all_tri_distances", &find_all_tri_distances,
         py::arg("start"), py::arg("comp"), py::arg("max_dist")=-1);
    // m.def("find_pairwise_distances", &find_pairwise_distances,
    //       py::arg("comp"), py::arg("target_dim")=-1, py::arg("max_cutoff")=-1);
    m.def("find_all_euclid_distances", &find_all_euclid_distances);
    m.def("get_neighborhood", &get_neighborhood);
    m.def("find_thresholded_component", &find_thresholded_component);
    m.def("find_pairwise_euclid_distances", &find_pairwise_euclid_distances);
    m.def("find_pairwise_tri_distances", &find_pairwise_tri_distances);
    
    m.def("find_nearest_neighbors", &find_nearest_neighbors);
     
};



