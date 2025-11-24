#include "ibex.h"
#include <iostream>

using namespace ibex;
using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " problema.bch [--trace]\n";
        return 1;
    }

    const char* bch = argv[1];
    bool trace = (argc >= 3 && std::string(argv[2]) == "--trace");

    try {
        System sys(bch);

        // Si hay objetivo, intentamos el optimizador; si falla, caemos al solver.
        if (sys.goal) {
            double eps_rel = 1e-7;
            double eps_abs = 1e-7;
            try {
                DefaultOptimizer opt(sys, eps_rel, eps_abs);
                opt.optimize(sys.box);
                const CovSolverData& data = opt.get_data();
                cout << "status: " << data.solver_status() << "\n";
                cout << "cells: "  << data.nb_cells() << "\n";
                cout << "solutions: " << data.nb_solution() << "\n";
                cout << "best (loup): " << opt.get_loup() << "\n";
                if (data.nb_solution() > 0) {
                    const IntervalVector& sol = data.inner(0);
                    cout << "sol[0]: " << sol << "\n";
                }
                return 0;
            } catch (const ibex::Exception& e) {
                cerr << "[fallback] Optimizer falló: " << e.what()
                     << ". Probando solver de satisfacción...\n";
                // sigue abajo al solver
            }
        }

        DefaultSolver solver(sys);
        solver.solve(sys.box);

        const CovSolverData& data = solver.get_data();
        cout << "status: " << data.solver_status() << "\n";
        cout << "cells: "  << data.nb_cells() << "\n";
        cout << "solutions: " << data.nb_solution() << "\n";

        if (data.nb_solution() > 0) {
            const IntervalVector& sol = data.inner(0);
            cout << "sol[0]: " << sol << "\n";
        }
    } catch (const ibex::Exception& e) {
        cerr << "IBEX exception: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
