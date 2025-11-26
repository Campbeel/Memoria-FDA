// run_vol_k_main.cpp
// Ejecuta solo la variante vol_k de FDA sobre un .bch

#include "fda_bridge_common.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace ibex;
using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " problema.bch [num_runs] [output.csv]\n";
        return 1;
    }
    std::string bch_file = argv[1];
    int num_runs = 1;
    std::string output_path;

    if (argc >= 3) {
        if (parse_int_arg(argv[2], num_runs)) {
            if (argc >= 4) output_path = argv[3];
        } else {
            output_path = argv[2];
        }
    }

    try {
        System sys(bch_file.c_str());
        OptContext opt(sys);
        IntervalVector root_box = sys.box;

        double eps_box   = 1e-9;
        int    max_iters = 100000;
        double T0        = 10.0;
        double k         = 0.5;
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
            RunStats stats = run_fda_vol_Tk(opt, root_box, r, T0, k,
                                            eps_V, beta, eps_box, max_iters);
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
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "[ERROR] Excepción desconocida\n";
        return 1;
    }
}
