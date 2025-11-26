// FDA_vol_k.cpp
//
// Variante volumen + temperatura determinista.

#include "fda_bridge_common.h"
#include <chrono>

using namespace ibex;
using namespace std;

// FD guiado por volumen relativo + temperatura determinista.
RunStats dive_vol_Tk(OptContext& opt,
                     const IntervalVector& start_box,
                     int run_id,
                     double T0,
                     double k,
                     double eps_V,
                     double beta,
                     double eps_box,
                     int max_iters,
                     double UB_global_in,
                     double& best_value_out)
{
    RunStats res;
    res.run_id  = run_id;
    res.variant = "vol_k";

    auto t0 = std::chrono::high_resolution_clock::now();

    bool use_volume = !start_box.is_unbounded();
    double V0 = use_volume ? start_box.volume() : 1.0;
    IntervalVector X = start_box;

    int depth = 0;
    long depth_accum = 0;

    double UB = UB_global_in;
    bool   has_ub = std::isfinite(UB);
    double UB_prev = UB;
    double T  = T0;
    Interval goal_bounds;
    int depth_vol = 0; // contador para aplicar corte por volumen tras algunos niveles

    for (int it = 0; it < max_iters; ++it) {
        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contract_with_goal(opt, X, goal_bounds, UB);
        if (X.is_empty()) break;

        if (use_volume && depth_vol >= 5) {
            double Vrel = X.volume() / V0;
            double tauV = eps_V * std::exp(-beta * depth); // τ_V(d): umbral decreciente
            if (Vrel <= tauV) {
                res.reached_optimum = std::isfinite(UB);
                break; // corte solo por volumen relativo
            }
        }

        if (!goal_bounds.is_empty()) {
            double f_ub = goal_bounds.ub();
            if (std::isfinite(f_ub) && f_ub < UB) { UB = f_ub; has_ub = true; }
            double f_mid = eval_at_mid(opt, X);
            if (std::isfinite(f_mid) && f_mid < UB) { UB = f_mid; has_ub = true; }
        }

        if (depth > 0 && X.max_diam() < eps_box) break;

        IntervalVector left(X.size()), right(X.size());
        bisect_box(opt, X, goal_bounds, left, right);

        double k_eff = adaptive_k(k, UB_prev, UB);
        UB_prev = UB;
        double T_child = k_eff * T / 2.0;
        double scoreL = heur_diving_score(opt, left, Interval())  + T_child;
        double scoreR = heur_diving_score(opt, right, Interval()) + T_child;

        if (scoreL >= scoreR) {
            X = left;
            goal_bounds = Interval();
            use_volume = !X.is_unbounded();
            if (use_volume) V0 = X.volume();
        } else {
            X = right;
            goal_bounds = Interval();
            use_volume = !X.is_unbounded();
            if (use_volume) V0 = X.volume();
        }

        depth++;
        T = T_child;
        depth_vol++;
    }

    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();

    res.best_value = UB;
    res.reached_optimum = std::isfinite(res.best_value);
    best_value_out = UB;
    return res;
}

RunStats run_fda_vol_Tk(OptContext& opt,
                        const IntervalVector& root_box,
                        int run_id,
                        double T0,
                        double k,
                        double eps_V,
                        double beta,
                        double eps_box,
                        int max_iters)
{
    double best_value = std::numeric_limits<double>::infinity();
    return dive_vol_Tk(opt, root_box, run_id,
                       T0, k, eps_V, beta, eps_box,
                       max_iters, best_value, best_value);
}

RunStats run_fda_vol_Tk_bb(OptContext& opt,
                           const IntervalVector& root_box,
                           int run_id,
                           double T0,
                           double k,
                           double eps_V,
                           double beta,
                           double eps_box,
                           int max_bb_nodes,
                           int max_iters_fdive)
{
    const double lb_tol = 1e-12;

    RunStats res;
    res.run_id  = run_id;
    res.variant = "vol_k_bb";

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<BbNode> active;
    active.push_back({root_box, Interval(), 0, T0});

    double UB_global = std::numeric_limits<double>::infinity();
    double UB_global_prev = UB_global;
    bool   has_ub    = false;
    long depth_accum = 0;
    int  bb_nodes_processed = 0;

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
            if (has_ub) res.reached_optimum = true;
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
        RunStats local_stats = dive_vol_Tk(opt, X, run_id,
                                           node.T, k, eps_V, beta,
                                           eps_box, max_iters_fdive,
                                           ub_hint, best_local);
        if (std::isfinite(best_local) && best_local < UB_global) {
            UB_global_prev = UB_global;
            UB_global = best_local;
            has_ub = true;
        }
        if (has_ub) res.reached_optimum = true;

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

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    }
    res.best_value = UB_global;
    res.reached_optimum = std::isfinite(res.best_value);
    return res;
}
