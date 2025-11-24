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
                   int base_depth,
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
    int depth = base_depth;
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
                     /*base_depth=*/0,
                     eps_box, max_iters, best_value, best_value);
}

// B&B sencillo pero acumulando las iteraciones del FD interno.
RunStats run_fda_base_bb(const System& sys,
                         Ctc& contractor,
                         const IntervalVector& root_box,
                         int run_id,
                         double eps_box,
                         int max_bb_nodes,
                         int max_iters_fdive)
{
    const double lb_tol = 1e-12;

    RunStats res;
    res.run_id  = run_id;
    res.variant = "base_bb";

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<BbNode> active;
    active.push_back({root_box, 0, 0.0});

    double UB_global = std::numeric_limits<double>::infinity();
    bool   has_ub    = false;
    long depth_accum = 0;
    int  bb_nodes_processed = 0;

    const int min_force_depth = 3;

    while (!active.empty() && bb_nodes_processed < max_bb_nodes) {
        BbNode node = active.back();
        active.pop_back();

        IntervalVector X = node.box;
        int depth = node.depth;

        if (!is_valid_box(X, sys.nb_var) || X.is_empty()) {
            if (depth < min_force_depth && X.size() == sys.nb_var) {
                IntervalVector left(X.size()), right(X.size());
                bisect_box(node.box, left, right);
                active.push_back({right, depth + 1, 0.0});
                active.push_back({left, depth + 1, 0.0});
            }
            continue;
        }

        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contractor.contract(X);
        if (X.is_empty()) {
            if (depth < min_force_depth) {
                IntervalVector left(node.box.size()), right(node.box.size());
                bisect_box(node.box, left, right);
                active.push_back({right, depth + 1, 0.0});
                active.push_back({left, depth + 1, 0.0});
            }
            continue;
        }

        double ub_hint = UB_global;
        if (sys.goal) {
            Interval f_int = sys.goal->eval(X);
            if (!f_int.is_empty()) {
                double f_lb = f_int.lb();
                double f_ub = f_int.ub();
                if (has_ub && std::isfinite(f_lb) && f_lb >= UB_global - lb_tol) continue;
                if (std::isfinite(f_ub) && f_ub < ub_hint) ub_hint = f_ub;
                if (!X.is_unbounded()) {
                    double f_mid = eval_at_mid(sys, X);
                    if (std::isfinite(f_mid) && f_mid < ub_hint) ub_hint = f_mid;
                }
            }
        }

        double best_local = ub_hint;
        RunStats local_stats = dive_base(sys, contractor, X,
                                         run_id, depth,
                                         eps_box, max_iters_fdive,
                                         ub_hint, best_local);

        if (std::isfinite(best_local) && best_local < UB_global) {
            UB_global = best_local;
            has_ub = true;
        }

        if (local_stats.nodes_visited > 0) {
            long local_nodes = local_stats.nodes_visited;
            long local_depth_sum =
                static_cast<long>(local_stats.avg_depth * local_nodes);
            res.nodes_visited += local_nodes - 1; // ya contamos el nodo B&B
            depth_accum += depth * (local_nodes - 1) + local_depth_sum;
            int local_max_depth = depth + local_stats.max_depth;
            if (local_max_depth > res.max_depth) res.max_depth = local_max_depth;
        }

        ++bb_nodes_processed;

        IntervalVector left(X.size()), right(X.size());
        bisect_box(X, left, right);

        active.push_back({right, depth + 1, 0.0});
        active.push_back({left, depth + 1, 0.0});
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    }
    res.best_value = UB_global;
    res.reached_optimum = has_ub;
    return res;
}
