// FDA_depth_k_rand.cpp
//
// Variante profundidad + temperatura con aleatoriedad.

#include "fda_bridge_common.h"
#include <chrono>

using namespace ibex;
using namespace std;

// FD con profundidad límite, temperatura determinista y ruido en T y d_max.
RunStats dive_depth_k_rand(OptContext& opt,
                           const IntervalVector& start_box,
                           int base_depth,
                           int run_id,
                           double T0,
                           double k,
                           int d_max_global,
                           double eps_box,
                           int max_iters,
                           double UB_global_in,
                           double& best_value_out,
                           std::mt19937& rng)
{
    RunStats res;
    res.run_id  = run_id;
    res.variant = "depth_k_rand";

    auto t0 = std::chrono::high_resolution_clock::now();

    IntervalVector X = start_box;
    int depth = 0;
    long depth_accum = 0;

    double UB = UB_global_in;
    double UB_prev = UB;
    double T  = T0;
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    Interval goal_bounds;

    int delta_d = 10;
    std::uniform_int_distribution<int> unif_int(-delta_d, delta_d);
    int d_max_run = std::max(1, d_max_global + unif_int(rng));
    int d_budget  = std::max(0, d_max_run - base_depth);

    for (int it = 0; it < max_iters; ++it) {
        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contract_with_goal(opt, X, goal_bounds, UB);
        if (X.is_empty()) break;

        if (!goal_bounds.is_empty()) {
            double f_ub = goal_bounds.ub();
            if (std::isfinite(f_ub) && f_ub < UB) UB = f_ub;
            double f_mid = eval_at_mid(opt, X);
            if (std::isfinite(f_mid) && f_mid < UB) UB = f_mid;
        }

        if (depth >= d_budget) {
            d_budget += d_max_run; // extiende el presupuesto para seguir explorando
        }

        if (depth > base_depth && X.max_diam() < eps_box) break;

        IntervalVector left(X.size()), right(X.size());
        bisect_box(opt, X, goal_bounds, left, right);

        double k_eff = adaptive_k(k, UB_prev, UB);
        UB_prev = UB;
        double rL = unif01(rng);
        double rR = unif01(rng);
        double T_childL = k_eff * T / 2.0 * rL;
        double T_childR = k_eff * T / 2.0 * rR;

        double scoreL = heur_diving_score(opt, left, Interval())  + T_childL;
        double scoreR = heur_diving_score(opt, right, Interval()) + T_childR;

        if (scoreL >= scoreR) {
            X = left;
            T = T_childL;
            goal_bounds = Interval();
        } else {
            X = right;
            T = T_childR;
            goal_bounds = Interval();
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
    res.reached_optimum = std::isfinite(UB);
    best_value_out = UB;
    return res;
}

RunStats run_fda_depth_Tk_rand(OptContext& opt,
                               const IntervalVector& root_box,
                               int run_id,
                               double T0,
                               double k,
                               int d_max,
                               double eps_box,
                               int max_iters,
                               std::mt19937& rng)
{
    double best_value = std::numeric_limits<double>::infinity();
    return dive_depth_k_rand(opt, root_box, /*base_depth=*/0, run_id,
                             T0, k, d_max, eps_box, max_iters,
                             best_value, best_value, rng);
}

RunStats run_fda_depth_Tk_rand_bb(OptContext& opt,
                                  const IntervalVector& root_box,
                                  int run_id,
                                  double T0,
                                  double k,
                                  int d_max_fdive,
                                  double eps_box,
                                  int max_bb_nodes,
                                  int max_iters_fdive)
{
    const double lb_tol = 1e-12;

    RunStats res;
    res.run_id  = run_id;
    res.variant = "depth_k_rand_bb";

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<BbNode> active;
    active.push_back({root_box, Interval(), 0, T0});

    double UB_global = std::numeric_limits<double>::infinity();
    double UB_global_prev = UB_global;
    bool   has_ub    = false;
    long depth_accum = 0;
    int  bb_nodes_processed = 0;

    std::mt19937 rng(123456u + run_id);
    const int min_force_depth = 3;

    while (!active.empty() && bb_nodes_processed < max_bb_nodes) {
        BbNode node = active.back();
        active.pop_back();

        IntervalVector X = node.box;
        int depth = node.depth;
        Interval goal_bounds = node.goal_bounds;

        if (!is_valid_box(X, opt.sys.nb_var) || X.is_empty()) {
            if (depth < min_force_depth && X.size() == opt.sys.nb_var) {
                IntervalVector left(node.box.size()), right(node.box.size());
                bisect_box(opt, node.box, node.goal_bounds, left, right);
                active.push_back({right, Interval(), depth + 1, node.T});
                active.push_back({left, Interval(), depth + 1, node.T});
            }
            continue;
        }

        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contract_with_goal(opt, X, goal_bounds, UB_global);
        if (X.is_empty()) {
            if (depth < min_force_depth) {
                IntervalVector left(node.box.size()), right(node.box.size());
                bisect_box(opt, node.box, node.goal_bounds, left, right);
                active.push_back({right, Interval(), depth + 1, node.T});
                active.push_back({left, Interval(), depth + 1, node.T});
            }
            continue;
        }

        double ub_hint = UB_global;
        if (!goal_bounds.is_empty()) {
            double f_lb = goal_bounds.lb();
            double f_ub = goal_bounds.ub();
            if (has_ub && std::isfinite(f_lb) && f_lb >= UB_global - lb_tol) continue;
            if (std::isfinite(f_ub) && f_ub < ub_hint) ub_hint = f_ub;
            if (!X.is_unbounded()) {
                double f_mid = eval_at_mid(opt, X);
                if (std::isfinite(f_mid) && f_mid < ub_hint) ub_hint = f_mid;
            }
        }

        double best_local = ub_hint;
        RunStats local_stats = dive_depth_k_rand(opt, X,
                                                 /*base_depth=*/depth,
                                                 run_id, node.T,
                                                 k, d_max_fdive,
                                                 eps_box, max_iters_fdive,
                                                 ub_hint, best_local, rng);
        if (std::isfinite(best_local) && best_local < UB_global) {
            UB_global_prev = UB_global;
            UB_global = best_local;
            has_ub = true;
        }

        if (local_stats.nodes_visited > 0) {
            long local_nodes = local_stats.nodes_visited;
            long local_depth_sum =
                static_cast<long>(local_stats.avg_depth * local_nodes);
            res.nodes_visited += local_nodes - 1;
            depth_accum += depth * (local_nodes - 1) + local_depth_sum;
            int local_max_depth = depth + local_stats.max_depth;
            if (local_max_depth > res.max_depth) res.max_depth = local_max_depth;
        }

        ++bb_nodes_processed;

        IntervalVector left(X.size()), right(X.size());
        bisect_box(opt, X, goal_bounds, left, right);

        double k_eff = adaptive_k(k, UB_global_prev, UB_global);
        double T_child = k_eff * node.T / 2.0;
        active.push_back({right, Interval(), depth + 1, T_child});
        active.push_back({left, Interval(), depth + 1, T_child});
        UB_global_prev = UB_global;
    }

    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();
    res.best_value = UB_global;
    res.reached_optimum = has_ub;
    return res;
}
