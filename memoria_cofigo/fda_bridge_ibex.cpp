// fda_bridge_ibex.cpp
//
// Puente entre problemas de IBEX (.bch -> System)
// y las variantes de Feasible Diving (FDA-base, profundidad, volumen, +rand).
// Implementaciones separadas por archivo para facilitar mantenimiento.
//
// Uso típico:
//   ./fda_bridge_ibex problema.bch base [num_runs] [output.csv]
//   ./fda_bridge_ibex problema.bch depth_k 5 resultados.csv
//   ./fda_bridge_ibex problema.bch depth_k_rand > salida.csv   # sigue funcionando la redirección
//   ./fda_bridge_ibex problema.bch vol_k_rand                # imprime en stdout si no se da archivo
//
// Columna CSV:
// run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal

#include "fda_bridge_common.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
using namespace ibex;
using namespace std;

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Uso: " << argv[0]
             << " problema.bch variante [num_runs] [output.csv]\n"
             << " variantes: base | base_bb | depth_k | depth_k_bb | depth_k_rand | depth_k_rand_bb | vol_k | vol_k_bb | vol_k_rand | vol_k_rand_bb\n";
        return 1;
    }

    std::string bch_file = argv[1];
    std::string variant  = argv[2];
    int num_runs = 1;
    std::string output_path;

    if (argc >= 4) {
        if (parse_int_arg(argv[3], num_runs)) {
            if (argc >= 5) output_path = argv[4];
        } else {
            output_path = argv[3];
        }
    }

    System sys(bch_file.c_str());
    OptContext opt(sys);
    IntervalVector root_box = sys.box;
    bool root_bounded = !root_box.is_unbounded();

    double eps_box   = 1e-9;
    int    max_iters = 100000;
    double T0        = 10.0;
    double k         = 0.5;
    int    d_max     = 100;
    double eps_V     = 1e-4;
    double beta      = 0.1;

    std::ofstream file_out;
    std::ostream* out = &cout;
    if (!output_path.empty()) {
        file_out.open(output_path);
        if (!file_out) {
            std::cerr << "[ERROR] No se pudo abrir " << output_path << " para escritura\n";
            return 1;
        }
        out = &file_out;
    }

    *out << "run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal\n";

    for (int r = 0; r < num_runs; ++r) {
        std::mt19937 rng(123456u + r);
        RunStats stats;
        if (variant == "base") {
            stats = run_fda_base(opt, root_box, r, eps_box, max_iters);
        } else if (variant == "base_bb") {
            int max_bb_nodes  = 5000000;
            int max_fd_iters  = 100000;
            stats = run_fda_base_bb(opt, root_box, r,
                                    eps_box, max_bb_nodes, max_fd_iters);
        } else if (variant == "depth_k") {
            if (!root_bounded) {
                stats = run_fda_base(opt, root_box, r, eps_box, max_iters);
                stats.variant = "base_depth_fallback";
            } else {
                stats = run_fda_depth_Tk(opt, root_box, r, T0, k,
                                         d_max, eps_box, max_iters);
            }
        } else if (variant == "depth_k_bb") {
            int d_max_fdive   = 50;
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_depth_Tk_bb(opt, root_box, r, T0, k,
                                        d_max_fdive, eps_box,
                                        max_bb_nodes, max_fd_iters);
        } else if (variant == "depth_k_rand") {
            if (!root_bounded) {
                stats = run_fda_base(opt, root_box, r, eps_box, max_iters);
                stats.variant = "base_depth_rand_fallback";
            } else {
                stats = run_fda_depth_Tk_rand(opt, root_box, r, T0, k,
                                              d_max, eps_box, max_iters, rng);
            }
        } else if (variant == "depth_k_rand_bb") {
            int d_max_fdive   = 50;
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_depth_Tk_rand_bb(opt, root_box, r, T0, k,
                                             d_max_fdive, eps_box,
                                             max_bb_nodes, max_fd_iters);
        } else if (variant == "vol_k") {
            stats = run_fda_vol_Tk(opt, root_box, r, T0, k,
                                   eps_V, beta, eps_box, max_iters);
        } else if (variant == "vol_k_bb") {
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_vol_Tk_bb(opt, root_box, r, T0, k,
                                      eps_V, beta, eps_box,
                                      max_bb_nodes, max_fd_iters);
        } else if (variant == "vol_k_rand") {
            stats = run_fda_vol_Tk_rand(opt, root_box, r, T0, k,
                                        eps_V, beta, eps_box, max_iters, rng);
        } else if (variant == "vol_k_rand_bb") {
            int max_bb_nodes  = 1000000;
            int max_fd_iters  = 100000;
            stats = run_fda_vol_Tk_rand_bb(opt, root_box, r, T0, k,
                                           eps_V, beta, eps_box,
                                           max_bb_nodes, max_fd_iters);
        } else {
            cerr << "Variante desconocida: " << variant << "\n";
            return 1;
        }

        *out << stats.run_id          << ","
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
