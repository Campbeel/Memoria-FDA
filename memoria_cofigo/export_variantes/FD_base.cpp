// FD_base.cpp
//
// Variante base de Feasible Diving y B&B mínimo.

#include "fda_bridge_common.h"
#include <chrono>

using namespace ibex;
using namespace std;

// Subproceso de Feasible Diving base (Algoritmo 3) sobre una caja inicial.
RunStats dive_base(const System& sys,
                   Ctc& contractor,
                   const IntervalVector& start_box,
                   int run_id,
                   double eps_box,
                   int max_iters,
                   double UB_global_in,
                   double& best_value_out)
{
    RunStats res;
    res.run_id  = run_id;
    res.variant = "base";

    auto t0 = std::chrono::high_resolution_clock::now();

    IntervalVector X = start_box;
    int depth = 0;
    long depth_accum = 0;

    double UB = UB_global_in;

    for (int it = 0; it < max_iters; ++it) {
        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contractor.contract(X);
        if (X.is_empty()) {
            break;
        }

        if (sys.goal) {
            Interval f_int = sys.goal->eval(X);
            (void)f_int;
            double f_mid = eval_at_mid(sys, X);
            if (f_mid < UB) UB = f_mid;
        }

        // parada por tamaño de caja
        if (X.max_diam() < eps_box) {
            res.reached_optimum = true;
            break;
        }

        IntervalVector left(X.size()), right(X.size());
        bisect_box(X, left, right);
        double scoreL = heur_diving_score(sys, left);
        double scoreR = heur_diving_score(sys, right);

        if (scoreL >= scoreR) {
            X = left;
        } else {
            X = right;
        }
        depth++;
    }

    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();

    res.best_value = UB;
    best_value_out = UB;
    return res;
}

RunStats run_fda_base(const System& sys,
                      Ctc& contractor,
                      const IntervalVector& root_box,
                      int run_id,
                      double eps_box,
                      int max_iters)
{
    double best_value = std::numeric_limits<double>::infinity();
    return dive_base(sys, contractor, root_box, run_id,
                     eps_box, max_iters, best_value, best_value);
}

// B&B sencillo estilo ibexsolve, sin heurística extra.
RunStats run_fda_base_bb(const System& sys,
                         Ctc& contractor,
                         const IntervalVector& root_box,
                         int run_id,
                         double eps_box,
                         int max_bb_nodes,
                         int /*max_iters_fdive*/)
{
    RunStats res;
    res.run_id  = run_id;
    res.variant = "base_bb";

    auto t0 = std::chrono::high_resolution_clock::now();

    // Delegamos en DefaultSolver de IBEX (similar a ibexsolve) para obtener un recorrido válido.
    Vector eps_x_min(sys.nb_var, eps_box);
    DefaultSolver solver(sys, eps_x_min, DefaultSolver::default_eps_x_max, true,
                         DefaultSolver::default_random_seed);
    solver.solve(sys.box);

    const CovSolverData& data = solver.get_data();
    res.nodes_visited   = static_cast<long>(data.nb_cells());
    res.best_value      = std::numeric_limits<double>::quiet_NaN(); // DefaultSolver no expone loup
    res.max_depth       = 0;
    res.avg_depth       = 0.0;
    res.reached_optimum = (data.nb_solution() > 0);

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();
    return res;
}
