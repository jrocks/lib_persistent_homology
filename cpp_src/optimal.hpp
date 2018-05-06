#ifndef OPTIMAL_HPP
#define OPTIMAL_HPP
 
#include <vector>
#include <unordered_map>
    
#include <Eigen/Core>
#include <Eigen/Sparse>
    
typedef Eigen::MatrixXd XMat;
typedef Eigen::SparseMatrix<double> SMat;
typedef Eigen::Triplet<double> Trip;

#include <ilcplex/ilocplex.h>
    
#include "cell_complex.hpp"
#include "filtration.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;


void optimize_cycle(int j, std::vector<int> &x, std::vector<int> &y, std::vector<int> &a, 
                    SMat &g, SMat &A, SMat &z, std::vector<double> &w, bool verbose=false) {
        
    // Create map of rows in full boundary matrix to rows in boundary matrix of just cells in optimization problem
    std::unordered_map<int, int> full_to_red_index;
    for(std::size_t i = 0; i < x.size(); i++) {
        full_to_red_index[x[i]] = i;
    }
    
    // Initialize environment and problem model
    IloEnv env;
    IloNumVarArray vars(env);
    IloModel mod(env);
    
    // Set up equality constraint
    IloNumArray zj(env, x.size());
    for(SMat::InnerIterator it(g, j); it; ++it) {
        if(full_to_red_index.count(it.row())) {
            zj[full_to_red_index[it.row()]] = it.value();
        }
    }
            
    // Range is vector of equality constraints, so the lower and upper bounds are both zj
    IloRangeArray range (env, zj, zj);
    // Add range to model
    mod.add(range);
    
    // Initialize cost functions
    IloObjective cost = IloAdd(mod, IloMinimize(env));
    
    // Construct model column by column
    
    // First add identity columns for x^+
    for(std::size_t i = 0; i < x.size(); i++) {
        // Add w[i]*x[i] to objective function
        IloNumColumn col = cost(w[i]);
        // Add column with 1.0 for x[i]
        col += range[i](1.0);
        // Add x[i] to list of variables
        vars.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
        col.end();
    }
    
    // Add negative identity columns for x^-
    for(std::size_t i = 0; i < x.size(); i++) {
        // Add w[i]*x[i] to objective function
        IloNumColumn col = cost(w[i]);
        // Add column with -1.0 for x[i]
        col += range[i](-1.0);
        // Add x[i] to list of variables
        vars.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
        col.end();
    }
    
    // Add columns from boundary matrix for y^+
    for(std::size_t i = 0; i < y.size(); i++) {
        // Add 0.0*y[i] to objective function
        IloNumColumn col = cost(0.0);
        // Add column with -A_j for y[i]
        for(SMat::InnerIterator it(A, y[i]); it; ++it) {
            if(full_to_red_index.count(it.row())) {
                col += range[full_to_red_index[it.row()]](-it.value());
            }
        }
        
        // Add y[i] to list of variables
        vars.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
        col.end();
    }
    
    // Add columns from boundary matrix for y^-
    for(std::size_t i = 0; i < y.size(); i++) {
        // Add 0.0*y[i] to objective function
        IloNumColumn col = cost(0.0);
        // Add column with A_j for y[i]
        for(SMat::InnerIterator it(A, y[i]); it; ++it) {
            if(full_to_red_index.count(it.row())) {
                col += range[full_to_red_index[it.row()]](it.value());
            }
        }
        
        // Add y[i] to list of variables
        vars.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
        col.end();
    }
    
    // Add columns from boundary matrix for a^+
    for(std::size_t i = 0; i < a.size(); i++) {
        // Add 0.0*a[i] to objective function
        IloNumColumn col = cost(0.0);
        // Add column with z_j for a[i]
        for(SMat::InnerIterator it(z, a[i]); it; ++it) {
            if(full_to_red_index.count(it.row())) {
                col += range[full_to_red_index[it.row()]](it.value());
            }
        }
        
        // Add y[i] to list of variables
        vars.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
        col.end();
    }
    
    // Add columns from boundary matrix for a^-
    for(std::size_t i = 0; i < a.size(); i++) {
        // Add 0.0*a[i] to objective function
        IloNumColumn col = cost(0.0);
        // Add column with -z_j for a[i]
        for(SMat::InnerIterator it(z, a[i]); it; ++it) {
            if(full_to_red_index.count(it.row())) {
                col += range[full_to_red_index[it.row()]](-it.value());
            }
        }
        
        // Add y[i] to list of variables
        vars.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
        col.end();
    }
    
    range.end();
    
    
    // Solve model

    IloCplex cplex(mod);
    
    if(!verbose) {
        cplex.setOut(env.getNullStream());
    }
    
    
    // cplex.exportModel("optimal_cycle.lp");

    cplex.solve();
    
    if(verbose) {
        cplex.out() << "solution status = " << cplex.getStatus() << std::endl;
        cplex.out() << std::endl;
        cplex.out() << "cost   = " << cplex.getObjValue() << std::endl;
        for (int i = 0; i < vars.getSize(); i++) {
            cplex.out() << "  x" << i << " = " << cplex.getValue(vars[i]) << std::endl;
        }
    }
    
    
    for(std::size_t i = 0; i < x.size(); i++) {
        z.insert(x[i], j) = cplex.getValue(vars[i]) - cplex.getValue(vars[x.size() + i]);
    }
    
    z.prune(0.0);
    
    env.end();    
    
}


std::unordered_map<int, std::vector<int> > calc_optimal_cycles(Filtration &filt, CellComplex &comp, std::vector<double> &weights, int dim=-1, bool verbose=false) {
        
    if(!comp.oriented) {
        py::print("Cell complex does not have oriented cells");
        return std::unordered_map<int, std::vector<int> >();
    }
    
    if((int)weights.size() != comp.ncells) {
        weights.assign(comp.ncells, 1.0);
    }
    
    SMat A(comp.ncells, comp.ncells);
    std::vector<int> cell_to_col(comp.ncells);    
    std::vector<int> col_to_cell(comp.ncells);
    
    
    std::vector<Trip> trip_list;
    
    std::vector<int> cell_order = filt.get_filtration();    
    for(int j = 0; j < comp.ncells; j++) {
        
        int cj = cell_order[j];

        auto facet_range = comp.get_facet_range(cj);
        auto coeff_range = comp.get_coeff_range(cj);
        for(auto itf = facet_range.first, itc = coeff_range.first; itf != facet_range.second; itf++, itc++) {
            int ci = *itf;
            int c = *itc;
            if(c != 0) {
                trip_list.emplace_back(cell_to_col[ci], j, c);
            }
        }
        
        cell_to_col[cj] = j;
        col_to_cell[j] = cj;
        
    }
    A.setFromTriplets(trip_list.begin(), trip_list.end());
            
    // Cycle basis
    SMat g(comp.ncells, comp.ncells);
    g.setIdentity();
    
    SMat z(comp.ncells, comp.ncells);
            
    // Simplices to include in optimization
    std::vector<int> x;
    // Columns of A to include in optimization
    std::vector<int> y;
    // Columns of z to include in optimization
    std::vector<int> a;
    // Weights of simplices in x
    std::vector<double> w;
    
     // row to reduced column with pivot in that row
    std::unordered_map<int, int> pivot_col;    
    
    for(int j = 0; j < A.cols(); j++) { 
        
        if(dim == -1 || comp.get_dim(col_to_cell[j]) == dim) {
            x.push_back(j);
            w.push_back(weights[col_to_cell[j]]);
        }
        
        while(A.col(j).nonZeros()) {
            
            SMat::ReverseInnerIterator it(A,j);
            int pivot_row = it.row();
                                     
            if(!pivot_col.count(pivot_row)) {
                break;
            }

            int l = pivot_col[pivot_row];
            
            double r = it.value() / SMat::ReverseInnerIterator(A,l).value();
            
            A.col(j) = (A.col(j) - r * A.col(l)).pruned();
            
            g.col(j) = (g.col(j) - r * g.col(l)).pruned();
            
            // py::print(j, "+", l, SMat(A.col(j)));
       
        }
        
        if(!A.col(j).nonZeros() && comp.get_dim(col_to_cell[j]) == dim) {
            optimize_cycle(j, x, y, a, g, A, z, w, verbose);
            a.push_back(j);
        }

        if(A.col(j).nonZeros()) {
            pivot_col[SMat::ReverseInnerIterator(A,j).row()] = j;
            
            if(comp.get_dim(col_to_cell[j]) == dim+1) {
                y.push_back(j);
                a.erase(std::remove(a.begin(), a.end(), SMat::ReverseInnerIterator(A,j).row()), a.end());
            }
        }
                                            
    }
    
    std::unordered_map<int, std::vector<int> > cycles;
    
    for(int j = 0; j < comp.ncells; j++) {
        if(A.col(j).nonZeros()) {
            continue;
        }
        
        int cj = col_to_cell[j];
        
        if(dim != -1 && comp.get_dim(cj) != dim) {
            continue;
        }
        
        cycles[cj];
        for(SMat::InnerIterator it(z,j); it; ++it) {
            cycles[cj].push_back(col_to_cell[it.row()]);
        }
    }
        
    return cycles;
    
                
}

    
#endif // OPTIMAL_HPP