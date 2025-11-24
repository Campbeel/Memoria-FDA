// FDA_vol_k_rand.cpp
//
// Variante volumen + temperatura con aleatoriedad.

#include "fda_bridge_common.h"
#include <chrono>

using namespace ibex;
using namespace std;

RunStats dive_vol_Tk_rand(const System& sys,
                          Ctc& contractor,
                          const IntervalVector& start_box,
                          int run_id,
                          double T0,
                          double k,
                          double eps_V,
                          double beta,
                          double eps_box,
                          int max_iters,
                          double UB_global_in,
                          double& best_value_out,
                          std::mt19937& rng)
{
    RunStats res;
    res.run_id  = run_id;
    res.variant = "vol_k_rand";

    auto t0 = std::chrono::high_resolution_clock::now();

    bool use_volume = !start_box.is_unbounded();
    double V0 = 1.0;
    if (use_volume) {
        V0 = start_box.volume();
    }
    IntervalVector X = start_box;

    int depth = 0;
    long depth_accum = 0;

    double UB = UB_global_in;
    double T  = T0;
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double eps_tau = 0.3;
    std::uniform_real_distribution<double> unif_tau(1.0 - eps_tau, 1.0 + eps_tau);
    double tau_factor = unif_tau(rng);

    for (int it = 0; it < max_iters; ++it) {
        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contractor.contract(X);
        if (X.is_empty()) break;

        if (use_volume) {
            double Vrel = X.volume() / V0;
            double tauV = eps_V * std::exp(-beta * depth) * tau_factor;
            if (Vrel <= tauV) break; // corte solo por volumen relativo
        }

        if (sys.goal) {
            Interval f_int = sys.goal->eval(X);
            double f_mid = eval_at_mid(sys, X);
            if (f_mid < UB) UB = f_mid;
        }

        if (depth > 20 && X.max_diam() < eps_box) {
            res.reached_optimum = true;
            break;
        }

        IntervalVector left(X.size()), right(X.size());
        bisect_box(X, left, right);

        double rL = unif01(rng);
        double rR = unif01(rng);

        double T_childL = k * T / 2.0 * rL;
        double T_childR = k * T / 2.0 * rR;

        double scoreL = heur_diving_score(sys, left)  + T_childL;
        double scoreR = heur_diving_score(sys, right) + T_childR;

        if (scoreL >= scoreR) {
            X = left;
            T = T_childL;
        } else {
            X = right;
            T = T_childR;
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

RunStats run_fda_vol_Tk_rand(const System& sys,
                             Ctc& contractor,
                             const IntervalVector& root_box,
                             int run_id,
                             double T0,
                             double k,
                             double eps_V,
                             double beta,
                             double eps_box,
                             int max_iters,
                             std::mt19937& rng)
{
    double best_value = std::numeric_limits<double>::infinity();
    return dive_vol_Tk_rand(sys, contractor, root_box, run_id,
                            T0, k, eps_V, beta, eps_box,
                            max_iters, best_value, best_value, rng);
}

RunStats run_fda_vol_Tk_rand_bb(const System& sys,
                                Ctc& contractor,
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
    res.variant = "vol_k_rand_bb";

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<BbNode> active;
    active.push_back({root_box, 0, T0});

    double UB_global = std::numeric_limits<double>::infinity();
    bool   has_ub    = false;
    long depth_accum = 0;
    int  bb_nodes_processed = 0;

    const int prune_min_depth = 5;
    const int min_force_depth = 3;

    std::mt19937 rng(123456u + run_id);

    while (!active.empty() && bb_nodes_processed < max_bb_nodes) {
        BbNode node = active.back();
        active.pop_back();

        IntervalVector X = node.box;
        int depth = node.depth;

        res.nodes_visited++;
        depth_accum += depth;
        if (depth > res.max_depth) res.max_depth = depth;

        contractor.contract(X);
        if (X.is_empty()) {
            if (depth < min_force_depth) {
                IntervalVector left(node.box.size()), right(node.box.size());
                bisect_box(node.box, left, right);
                active.push_back({right, depth + 1, node.T});
                active.push_back({left, depth + 1, node.T});
            }
            continue;
        }

        double ub_hint = UB_global;
        if (sys.goal) {
            Interval f_int = sys.goal->eval(X);
            double f_lb = f_int.lb();
            double f_ub = f_int.ub();
            if (has_ub && std::isfinite(f_lb) &&
                depth >= prune_min_depth &&
                f_lb >= UB_global - lb_tol) continue;
            if (std::isfinite(f_ub) && f_ub < ub_hint) ub_hint = f_ub;
            if (!X.is_unbounded()) {
                double f_mid = eval_at_mid(sys, X);
                if (std::isfinite(f_mid) && f_mid < ub_hint) ub_hint = f_mid;
            }
        }

        double best_local = ub_hint;
        RunStats local_stats = dive_vol_Tk_rand(sys, contractor, X, run_id,
                                                node.T, k, eps_V, beta,
                                                eps_box, max_iters_fdive,
                                                ub_hint, best_local, rng);
        if (std::isfinite(best_local) && best_local < UB_global) {
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
        bisect_box(X, left, right);

        double T_child = k * node.T / 2.0;
        active.push_back({right, depth + 1, T_child});
        active.push_back({left, depth + 1, T_child});
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds =
        std::chrono::duration<double>(t1 - t0).count();
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    }
    res.best_value = UB_global;
    return res;
}
