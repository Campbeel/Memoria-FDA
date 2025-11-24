// fda_bridge_ibex.cpp
//
// Puente entre problemas de IBEX (.bch -> System)
// y las variantes de Feasible Diving (FDA-base, profundidad, volumen, +rand).
// Implementaciones separadas por archivo para facilitar mantenimiento.
//
// Uso típico:
//   ./fda_bridge_ibex problema.bch base          > salida.csv
//   ./fda_bridge_ibex problema.bch depth_k      >> salida.csv
//   ./fda_bridge_ibex problema.bch depth_k_rand >> salida.csv
//   ./fda_bridge_ibex problema.bch vol_k        >> salida.csv
//   ./fda_bridge_ibex problema.bch vol_k_rand   >> salida.csv
//
// Columna CSV:
// run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal

#include "fda_bridge_common.h"
#include <chrono>
#include <iostream>
#include <random>
#include <string>
using namespace ibex;
using namespace std;

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Uso: " << argv[0]
             << " problema.bch variante [num_runs]\n"
             << " variantes: base | base_bb | depth_k | depth_k_bb | depth_k_rand | depth_k_rand_bb | vol_k | vol_k_bb | vol_k_rand | vol_k_rand_bb\n";
        return 1;
    }

    std::string bch_file = argv[1];
    std::string variant  = argv[2];
    int num_runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    System sys(bch_file.c_str());
    IntervalVector root_box = sys.box;
    CtcHC4 hc4(sys);
    bool root_bounded = !root_box.is_unbounded();

    double eps_box   = 1e-9;
    int    max_iters = 100000;
    double T0        = 10.0;
    double k         = 0.5;
    int    d_max     = 100;
    double eps_V     = 1e-4;
    double beta      = 0.1;

    cout << "run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal\n";

    for (int r = 0; r < num_runs; ++r) {
        std::mt19937 rng(123456u + r);
        RunStats stats;
        if (variant == "base") {
            stats = run_fda_base(sys, hc4, root_box, r, eps_box, max_iters);
        } else if (variant == "base_bb") {
            int max_bb_nodes  = 5000000;
            int max_fd_iters  = 100000;
            stats = run_fda_base_bb(sys, hc4, root_box, r,
                                    eps_box, max_bb_nodes, max_fd_iters);
        } else if (variant == "depth_k") {
            if (!root_bounded) {
                stats = run_fda_base(sys, hc4, root_box, r, eps_box, max_iters);
                stats.variant = "base_depth_fallback";
            } else {
                stats = run_fda_depth_Tk(sys, hc4, root_box, r, T0, k,
                                         d_max, eps_box, max_iters);
            }
        } else if (variant == "depth_k_bb") {
            int d_max_fdive   = 50;
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_depth_Tk_bb(sys, hc4, root_box, r, T0, k,
                                        d_max_fdive, eps_box,
                                        max_bb_nodes, max_fd_iters);
        } else if (variant == "depth_k_rand") {
            if (!root_bounded) {
                stats = run_fda_base(sys, hc4, root_box, r, eps_box, max_iters);
                stats.variant = "base_depth_rand_fallback";
            } else {
                stats = run_fda_depth_Tk_rand(sys, hc4, root_box, r, T0, k,
                                              d_max, eps_box, max_iters, rng);
            }
        } else if (variant == "depth_k_rand_bb") {
            int d_max_fdive   = 50;
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_depth_Tk_rand_bb(sys, hc4, root_box, r, T0, k,
                                             d_max_fdive, eps_box,
                                             max_bb_nodes, max_fd_iters);
        } else if (variant == "vol_k") {
            stats = run_fda_vol_Tk(sys, hc4, root_box, r, T0, k,
                                   eps_V, beta, eps_box, max_iters);
        } else if (variant == "vol_k_bb") {
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_vol_Tk_bb(sys, hc4, root_box, r, T0, k,
                                      eps_V, beta, eps_box,
                                      max_bb_nodes, max_fd_iters);
        } else if (variant == "vol_k_rand") {
            stats = run_fda_vol_Tk_rand(sys, hc4, root_box, r, T0, k,
                                        eps_V, beta, eps_box, max_iters, rng);
        } else if (variant == "vol_k_rand_bb") {
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_vol_Tk_rand_bb(sys, hc4, root_box, r, T0, k,
                                           eps_V, beta, eps_box,
                                           max_bb_nodes, max_fd_iters);
        } else {
            cerr << "Variante desconocida: " << variant << "\n";
            return 1;
        }

        cout << stats.run_id          << ","
             << stats.variant         << ","
             << stats.best_value      << ","
             << stats.nodes_visited   << ","
             << stats.elapsed_seconds << ","
             << stats.max_depth       << ","
             << stats.avg_depth       << ","
             << (stats.reached_optimum ? 1 : 0)
             << "\n";
    }

    return 0;
}
