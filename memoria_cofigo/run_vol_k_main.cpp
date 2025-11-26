// run_vol_k_main.cpp
// Ejecuta solo la variante vol_k de FDA sobre un .bch

#include "fda_bridge_common.h"
#include <chrono>
#include <iostream>
#include <random>
#include <string>

using namespace ibex;
using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " problema.bch [num_runs]\n";
        return 1;
    }
    std::string bch_file = argv[1];
    int num_runs = (argc >= 3) ? std::atoi(argv[2]) : 1;

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

        cout << "run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal\n";
        for (int r = 0; r < num_runs; ++r) {
            RunStats stats = run_fda_vol_Tk(opt, root_box, r, T0, k,
                                            eps_V, beta, eps_box, max_iters);
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
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "[ERROR] Excepción desconocida\n";
        return 1;
    }
}
