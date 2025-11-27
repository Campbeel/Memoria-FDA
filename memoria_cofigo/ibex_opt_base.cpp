// ibex_opt_base.cpp
// Optimización base usando Ibex pero sin loup XTaylor (evita LP/CLP).

#include "ibex.h"
#include <cstdlib>
#include <iostream>
#include <atomic>
#include <thread>
#include <chrono>
#include <string>

using namespace ibex;
using namespace std;

int main(int argc, char** argv) {
	if (argc < 2) {
		cerr << "Uso: " << argv[0] << " problema.bch [rel_eps_f] [abs_eps_f] [--lp]\n";
		return 1;
	}

	const char* filename = argv[1];
	const double rel_eps_f = (argc >= 3) ? atof(argv[2]) : OptimizerConfig::default_rel_eps_f;
	const double abs_eps_f = (argc >= 4) ? atof(argv[3]) : OptimizerConfig::default_abs_eps_f;
	const bool use_lp = (argc >= 5) && std::string(argv[4])=="--lp";

	try {
		System sys(filename);
		if (!sys.goal) {
			cerr << "[error] El archivo no es un problema de optimizacion (no hay goal).\n";
			return 1;
		}

		// Componentes:
		NormalizedSystem norm_sys(sys);
		ExtendedSystem ext_sys(sys);

		CtcHC4 hc4(ext_sys, 0.01, true);
		CtcFixPoint contractor(hc4, 0.1);

		Bsc* bisector = nullptr;
		LoupFinder* loup_finder = nullptr;
		if (use_lp) {
			// LoupFinderDefault usa XTaylor (LP) + probing; bisector sin LP para evitar crash en CLP.
			bisector = new OptimLargestFirst(ext_sys.goal_var(), /*choose_obj=*/false, OptimizerConfig::default_eps_x);
			loup_finder = new LoupFinderDefault(norm_sys, /*inHC4=*/true);
			std::cerr << "[ibex_opt_base] modo LP activado (XTaylor + bisector LF)\n";
		} else {
			bisector = new SmearSumRelative(ext_sys, OptimizerConfig::default_eps_x);
			loup_finder = new LoupFinderInHC4(norm_sys);
		}
		CellDoubleHeap buffer(ext_sys);

		Optimizer opt(sys.nb_var, contractor, *bisector, *loup_finder, buffer,
		              ext_sys.goal_var(), OptimizerConfig::default_eps_x,
		              rel_eps_f, abs_eps_f, /*stats=*/false);

		// Heartbeat: avisa cada 5s que sigue trabajando.
		std::atomic<bool> stop{false};
		std::thread ping([&]() {
			using namespace std::chrono_literals;
			auto t0 = std::chrono::steady_clock::now();
			while (!stop.load(std::memory_order_relaxed)) {
				std::this_thread::sleep_for(5s);
				if (stop.load(std::memory_order_relaxed)) break;
				auto t1 = std::chrono::steady_clock::now();
				double elapsed = std::chrono::duration<double>(t1 - t0).count();
				std::cerr << "[ibex_opt_base] trabajando... " << elapsed << " s\n";
			}
		});

		try {
			opt.optimize(sys.box);
		} catch (...) {
			stop = true;
			ping.join();
			throw;
		}
		stop = true;
		ping.join();

		cout << "================ Resultado ================\n";
		cout << "f* en [" << opt.get_uplo() << ", " << opt.get_loup() << "]\n";
		cout << "punto factible: " << opt.get_loup_point() << "\n";
		cout << "celdas: " << opt.get_nb_cells() << "\n";
		cout << "tiempo: " << opt.get_time() << " s\n";
		delete bisector;
		delete loup_finder;
		return 0;
	} catch (const exception& e) {
		cerr << "[error] " << e.what() << "\n";
		return 1;
	} catch (...) {
		cerr << "[error] Excepcion desconocida\n";
		return 1;
	}
}
