//============================================================================
//                                  I B E X
//
//                               ************
//                                  IbexOpt
//                               ************
//
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Last Update : Feb 13, 2025
//============================================================================

#include "ibex.h"
#include "parse_args.h"
#include "ibex_OptimLargestFirst.h"
#include <memory>

#include <sstream>
#include <random>
#include <cstring>
#include "temp_buffer.h"
#include <filesystem>

using namespace std;
using namespace ibex;

// Simple subclass to expose protected getters from DefaultOptimizerConfig.
struct FdConfig : public DefaultOptimizerConfig {
	using DefaultOptimizerConfig::DefaultOptimizerConfig;
	using DefaultOptimizerConfig::nb_var;
	using DefaultOptimizerConfig::get_ctc;
	using DefaultOptimizerConfig::get_bsc;
	using DefaultOptimizerConfig::get_loup_finder;
	using DefaultOptimizerConfig::get_cell_buffer;
	using DefaultOptimizerConfig::goal_var;
	using DefaultOptimizerConfig::get_ext_sys;
	using DefaultOptimizerConfig::get_norm_sys;
	using DefaultOptimizerConfig::get_eps_x;
	using DefaultOptimizerConfig::get_rel_eps_f;
	using DefaultOptimizerConfig::get_abs_eps_f;
	using DefaultOptimizerConfig::with_statistics;
	using DefaultOptimizerConfig::get_trace;
};

int main(int argc, char** argv) {
	auto exe_basename = [](const char* path) -> std::string {
		const char* slash = strrchr(path, '/');
		return slash ? std::string(slash+1) : std::string(path);
	};

#ifdef __IBEX_NO_LP_SOLVER__
	ibex_error("ibexopt requires a LP Solver (use -DLP_LIB with cmake)");
	exit(1);
#endif

	stringstream _rel_eps_f, _abs_eps_f, _eps_h, _random_seed, _eps_x;
	_rel_eps_f << "Relative precision on the objective. Default value is 1e" << round(::log10(OptimizerConfig::default_rel_eps_f)) << ".";
	_abs_eps_f << "Absolute precision on the objective. Default value is 1e" << round(::log10(OptimizerConfig::default_abs_eps_f)) << ".";
	_eps_h << "Equality relaxation value. Default value is 1e" << round(::log10(NormalizedSystem::default_eps_h)) << ".";
	_random_seed << "Random seed (useful for reproducibility). Default value is " << DefaultOptimizerConfig::default_random_seed << ".";
	_eps_x << "Precision on the variable (**Deprecated**). Default value is 0.";

	args::ArgumentParser parser("********* IbexOpt (defaultoptimizer) *********.", "Solve a Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::Flag version(parser, "version", "Display the version of this plugin (same as the version of Ibex).", {'v',"version"});
	args::ValueFlag<double> rel_eps_f(parser, "float", _rel_eps_f.str(), {'r', "rel-eps-f"});
	args::ValueFlag<double> abs_eps_f(parser, "float", _abs_eps_f.str(), {'a', "abs-eps-f"});
	args::ValueFlag<double> eps_h(parser, "float", _eps_h.str(), {"eps-h"});
	args::ValueFlag<double> timeout(parser, "float", "Timeout (time in seconds). Default value is +oo.", {'t', "timeout"});
	args::ValueFlag<double> random_seed(parser, "float", _random_seed.str(), {"random-seed"});
	args::ValueFlag<double> eps_x_arg(parser, "float", _eps_x.str(), {"eps-x"});
	args::ValueFlag<int>    simpl_level(parser, "int", "Expression simplification level. Possible values are:\n"
			"\t\t* 0:\tno simplification at all (fast).\n"
			"\t\t* 1:\tbasic simplifications (fairly fast). E.g. x+1+1 --> x+2\n"
			"\t\t* 2:\tmore advanced simplifications without developing (can be slow). E.g. x*x + x^2 --> 2x^2\n"
			"\t\t* 3:\tsimplifications with full polynomial developing (can blow up!). E.g. x*(x-1) + x --> x^2\n"
			"Default value is : 1.", {"simpl"});
	args::ValueFlag<double> initial_loup(parser, "float", "Intial \"loup\" (a priori known upper bound).", {"initial-loup"});
	args::ValueFlag<string> input_file(parser, "filename", "COV input file. The file contains "
			"optimization data in the COV (binary) format.", {'i',"input"});
	args::ValueFlag<string> output_file(parser, "filename", "COV output file. The file will contain the "
			"optimization data in the COV (binary) format. See --format", {'o',"output"});
	args::Flag rigor(parser, "rigor", "Activate rigor mode (certify feasibility of equalities).", {"rigor"});
	args::Flag kkt(parser, "kkt", "Activate contractor based on Kuhn-Tucker conditions.", {"kkt"});
	args::Flag output_no_obj(parser, "output-no-obj", "Generate a COV with domains of variables only (not objective values).", {"output-no-obj"});
	args::Flag trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});
	args::Flag stats(parser, "stats", "Enable statistics. Note: This may slightly deteriorate performances. Only statistics for the LP solver exist so far.", {"stats"});
	args::Flag format(parser, "format", "Give a description of the COV format used by IbexOpt", {"format"});
	args::ValueFlag<string> no_split_arg(parser, "vars","Prevent some variables to be bisected, separated by '+'.\nExample: --no-split=x+y",{"no-split"});
	args::Flag quiet(parser, "quiet", "Print no report on the standard output.",{'q',"quiet"});
	args::ValueFlag<string> fd_mode(parser, "string", "FD mode (depth_k). Optional; if set, overrides the bisector with a depth-oriented variant.", {"fd-mode"});

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	if (version) {
		cout << "IbexOpt Release " << _IBEX_RELEASE_ << endl;
		exit(0);
	}

	if (format) {
		cout << CovOptimData::format();
		exit(0);
	}

	if (filename.Get()=="") {
		ibex_error("no input file (try ibexopt --help)");
		exit(1);
	}

	try {

		System *sys;
		string exe_name = exe_basename(argv[0]);

		string extension = filename.Get().substr(filename.Get().find_last_of('.')+1);
		if (extension == "nl") {
			cerr << "\n\033[31mAMPL files can only be read with optimizer04 (ibex-opt-extra package).\n\n";
			exit(0);
		}
		else
			// Load a system of equations
			sys = new System(filename.Get().c_str(), simpl_level? simpl_level.Get() : ExprNode::default_simpl_level);

		FdConfig config(*sys);

		string output_cov_file; // cov output file
		bool overwitten=false;  // is it overwritten?
		string cov_copy;

		if (!sys->goal) {
			ibex_error(" input file has not goal (it is not an optimization problem).");
		}

		auto infer_mode = [&](const std::string& exe)->std::string {
			if (exe.find("depth_k_rand")!=string::npos) return "depth_k_rand";
			if (exe.find("depth_k")!=string::npos)      return "depth_k";
			if (exe.find("vol_k_rand")!=string::npos)   return "vol_k_rand";
			if (exe.find("vol_k")!=string::npos)        return "vol_k";
			return "";
		};

		std::string inferred = infer_mode(exe_name);
		std::string fd_choice = fd_mode ? fd_mode.Get() : inferred;

		bool use_fd_variant = (fd_choice=="depth_k" || fd_choice=="depth_k_rand" ||
		                       fd_choice=="vol_k"   || fd_choice=="vol_k_rand");

		bool want_fd_logging = use_fd_variant;
		if (want_fd_logging) {
			cout << "  [fd-debug] sys.box empty? " << (sys->box.is_empty() ? "yes" : "no")
			     << " unbounded? " << (sys->box.is_unbounded() ? "yes" : "no")
			     << " dim=" << sys->nb_var << "\n";
			if (!sys->box.is_empty() && !sys->box.is_unbounded()) {
				try {
					cout << "  [fd-debug] root volume: " << sys->box.volume() << "\n";
				} catch (...) {
					cout << "  [fd-debug] root volume: (error computing)\n";
				}
			}
		}

		if (!quiet) {
			cout << endl << "************************ setup ************************" << endl;
			cout << "  file loaded:\t\t" << filename.Get() << endl;
		}

		if (rel_eps_f) {
			config.set_rel_eps_f(rel_eps_f.Get());

			if (!quiet)
				cout << "  rel-eps-f:\t\t" << rel_eps_f.Get() << "\t(relative precision on objective)" << endl;
		}

		if (abs_eps_f) {
			config.set_abs_eps_f(abs_eps_f.Get());
			if (!quiet)
				cout << "  abs-eps-f:\t\t" << abs_eps_f.Get() << "\t(absolute precision on objective)" << endl;
		}

		if (eps_h) {
			config.set_eps_h(eps_h.Get());
			if (!quiet)
				cout << "  eps-h:\t\t" << eps_h.Get() << "\t(equality thickening)" << endl;
		}

		if (eps_x_arg) {
			if (!quiet)
				cout << "  eps-x:\t\t" << eps_x_arg.Get() << "\t(precision on variables domain)" << endl;
		}

		Vector eps_x(sys->nb_var, eps_x_arg ? eps_x_arg.Get() : OptimizerConfig::default_eps_x);

		if (no_split_arg) {
			if (!quiet)
				cout << "  don't split:\t\t";

			vector<const ExprNode*> no_split = parse_symbols_list(sys->args, no_split_arg.Get());

			if (!quiet) {
				for (vector<const ExprNode*>::const_iterator it=no_split.begin(); it!=no_split.end(); ++it)
					cout << **it << ' ';
				cout << endl;
			}

			if (!no_split.empty()) {
				// we use VarSet for convenience (handling of indexed symbols)
				VarSet varset(sys->f_ctrs,no_split,true);

				for (int i=0; i<varset.nb_var; i++) {
					eps_x[varset.var(i)]=POS_INFINITY;
				}
				for (vector<const ExprNode*>::iterator it=no_split.begin(); it!=no_split.end(); ++it) {
					cleanup(**it,false);
				}
			}
		}

		config.set_eps_x(eps_x);

		// This option certifies feasibility with equalities
		if (rigor) {
			config.set_rigor(rigor.Get());
			if (!quiet)
				cout << "  rigor mode:\t\tON\t(feasibility of equalities certified)" << endl;
		}

	if (kkt) {
		config.set_kkt(kkt.Get());
		if (!quiet)
			cout << "  KKT contractor:\tON" << endl;
	}

		if (simpl_level)
			cout << "  symbolic simpl level:\t" << simpl_level.Get() << "\t" << endl;

		if (initial_loup) {
			if (!quiet)
				cout << "  initial loup:\t\t" << initial_loup.Get() << " (a priori upper bound of the minimum)" << endl;
		}

		// Fix the random seed for reproducibility.
		if (random_seed) {
			config.set_random_seed(random_seed.Get());
			if (!quiet)
				cout << "  random seed:\t\t" << random_seed.Get() << endl;
		}

		if (input_file) {
			if (!quiet) {
				cout << "  input COV file:\t" << input_file.Get().c_str() << "\n";
			}
		}

		if (output_file) {
			output_cov_file = output_file.Get();
		} else {
			// got from stackoverflow.com:
			string::size_type const p(filename.Get().find_last_of('.'));
			// filename without extension
			string filename_no_ext=filename.Get().substr(0, p);
			stringstream ss;
			ss << filename_no_ext << ".cov";
			output_cov_file=ss.str();

			ifstream file;
			file.open(output_cov_file.c_str(), ios::in); // to check if it exists

			if (file.is_open()) {
				overwitten = true;
				stringstream ss;
				ss << output_cov_file << "~";
				cov_copy=ss.str();
				// got from stackoverflow.com:
				ofstream dest(cov_copy, ios::binary);

			    istreambuf_iterator<char> begin_source(file);
			    istreambuf_iterator<char> end_source;
			    ostreambuf_iterator<char> begin_dest(dest);
			    copy(begin_source, end_source, begin_dest);
			}
			file.close();
		}

		if (!quiet) {
			cout << "  output COV file:\t" << output_cov_file << "\n";
		}

		// This option limits the search time
		if (timeout) {
			if (!quiet)
				cout << "  timeout:\t\t" << timeout.Get() << "s" << endl;
			config.set_timeout(timeout.Get());
		}

		// This option prints each better feasible point when it is found
		if (trace) {
			if (!quiet)
				cout << "  trace:\t\tON" << endl;
			config.set_trace(trace.Get());
		}

		// This option enables statistics
		if (stats) {
			if (!quiet)
				cout << "  statistics:\t\tON" << endl;
			config.set_statistics(stats.Get());
		}

		// Question: is really inHC4 good?
		config.set_inHC4(true);

		if (!config.with_inHC4()) {
			cerr << "\n  \033[33mwarning: inHC4 disabled\033[0m (unimplemented operator)" << endl;
		}

		if (output_no_obj) {
			cout << "  Generates COV with:\tvariable domains only\n";
			config.set_extended_cov(false);
		}

		if (!quiet) {
			cout << "*******************************************************" << endl << endl;
		}

		// Build the optimizer (default or FD mode)
		Optimizer* opt_ptr = NULL;
		std::unique_ptr<OptimLargestFirst> fd_bisector;
		std::unique_ptr<Optimizer> opt_owner;
		std::unique_ptr<TempBuffer> temp_buffer;
		TempBuffer* temp_raw = nullptr;

		if (use_fd_variant) {
			string mode = fd_choice;
			if (mode=="depth_k" || mode=="depth_k_rand" || mode=="vol_k" || mode=="vol_k_rand") {
				// Para estas variantes, usamos un bisector OptimLargestFirst sin dividir el objetivo.
				ExtendedSystem& ext_sys = config.get_ext_sys();
				const Vector& ex = config.get_eps_x();
				double eps_scalar = ex.size()>=1 ? ex[0] : OptimizerConfig::default_eps_x;
				fd_bisector.reset(new OptimLargestFirst(ext_sys.goal_var(), /*choose_obj=*/false, eps_scalar));
				CellBufferOptim* buffer_ptr = &config.get_cell_buffer();
				double V0_ref = 1.0;
				if (!sys->box.is_empty() && !sys->box.is_unbounded()) {
					try { V0_ref = sys->box.volume(); } catch (...) { V0_ref = 1.0; }
				}
				bool is_rand = (mode.find("_rand")!=string::npos);
				bool use_depth = (mode.find("depth")!=string::npos);
				bool use_vol   = (mode.find("vol")!=string::npos);

				// Semilla base: si no dan random-seed, usamos un valor no determinista.
				uint64_t seed = random_seed ? static_cast<uint64_t>(random_seed.Get()) : static_cast<uint64_t>(std::random_device{}());

				auto u01 = [&]() {
					static thread_local std::mt19937_64 gen(seed);
					static thread_local std::uniform_real_distribution<double> dist(0.0,1.0);
					return dist(gen);
				};

				TempBuffer::Params params;
				// k muestreado por corrida (±20%) y con ruido de tie-break más alto
				double k_base = 10.0;
				double k_spread = is_rand ? 0.3 : 0.2; // más dispersión en rand
				params.k = k_base * (1.0 + k_spread * (2.0*u01()-1.0));
				params.k = std::max(5.0, std::min(15.0, params.k));
				params.bias = 1e-2; // peso térmico moderado (se ajusta más abajo si dominio enorme)
				params.T0 = 100.0;
				// Ajustes adaptativos según el volumen inicial: dominios grandes => vol_ratio más alto, depth_cut más bajo.
				double logV = 0.0;
				if (V0_ref>0) logV = std::log10(std::max(1e-12, V0_ref));
				int nv = sys->nb_var;
				auto scale_vol = [&](double base) {
					double v = base;
					if (logV > 6) v *= 1.5;   // volumen enorme -> umbral más alto
					else if (logV > 3) v *= 1.2;
					else if (logV < -2) v *= 0.7;  // volumen pequeño -> umbral más bajo
					if (nv > 20) v *= 1.2;         // más vars: subir umbral para explorar
					if (nv < 8) v *= 0.9;          // pocos vars: bajar umbral
					return v;
				};
				auto scale_depth = [&](int base) {
					int d = base;
					if (logV > 6) d = std::max(2, d-2); // dominios grandes: menos profundidad
					else if (logV > 3) d = std::max(3, d-1);
					else if (logV < -2) d = d+1;             // dominios pequeños: más profundidad
					if (nv > 20) d = std::max(2, d-1);       // muchos vars: bajar profundidad
					if (nv < 8) d = d+1;                     // pocos vars: permitir más
					return d;
				};
				// depth cut muestreado por corrida (alrededor de 2-3 para disparar)
				if (use_depth) {
					int d = static_cast<int>(std::round(2.5 * (1.0 + 0.1 * (2.0*u01()-1.0))));
					d = scale_depth(d);
					params.depth_cut = std::max(2, std::min(6, d));
					params.depth_cut_jitter = 0.1; // variación por nodo
					params.depth_hard_cut = params.depth_cut + 3; // tope duro para evitar explosiones
				} else {
					params.depth_cut = 0;
				}
				// vol cut muestreado por corrida (alrededor de 5%)
				if (use_vol) {
					// Umbral medio para que dispare en medium.
					double v = 0.1 * (1.0 + 0.1 * (2.0*u01()-1.0));
					v = scale_vol(v);
					params.vol_ratio_cut = std::max(0.07, std::min(0.13, v));
					params.vol_cut_jitter = 0.1; // variación por nodo
					params.vol_hard_ratio = 0.02; // corte duro adicional para evitar tiempos enormes
				} else {
					params.vol_ratio_cut = 0.0;
				}
				params.V0_ref = V0_ref;
				params.depth_penalty = 1e6;
				params.vol_penalty = 1e6;
				params.rand_k = is_rand;
				params.rand_seed = seed;
				params.tie_noise = 1e-3; // ruido para desempate moderado

				// Ajustes extra para dominios muy grandes: reducir ruido térmico y profundidades.
				if (logV > 6) {
					params.bias *= 0.1;          // menos peso térmico
					params.tie_noise *= 0.1;     // menos ruido
					params.vol_hard_ratio = 0.03; // corte duro moderado en volumen relativo
					if (params.vol_ratio_cut > 0.15) params.vol_ratio_cut = 0.15;
					if (use_depth) {
						params.depth_cut = std::max(2, std::min(4, params.depth_cut));
						params.depth_hard_cut = params.depth_cut + 2;
					}
					if (is_rand) {
						// bajar la dispersión de k en dominios enormes para no ir a rutas muy malas
						double kspread = 0.15;
						params.k = k_base * (1.0 + kspread * (2.0*u01()-1.0));
						params.k = std::max(5.0, std::min(15.0, params.k));
					}
				} else if (logV > 3) {
					params.bias *= 0.3;
					params.tie_noise *= 0.3;
					params.vol_hard_ratio = 0.025;
					if (params.vol_ratio_cut > 0.12) params.vol_ratio_cut = 0.12;
				}

				// Aplicamos TempBuffer a las variantes FD (depth/vol) con parámetros moderados.
				if (use_depth || use_vol) {
					temp_buffer.reset(new TempBuffer(ext_sys, ext_sys.goal_var(), params));
					temp_raw = temp_buffer.get();
					buffer_ptr = temp_raw;
				}

				opt_owner.reset(new Optimizer(
					config.nb_var(),
					config.get_ctc(),
					*fd_bisector,
					config.get_loup_finder(),
					*buffer_ptr,
					config.goal_var(),
					OptimizerConfig::default_eps_x,
					config.get_rel_eps_f(),
					config.get_abs_eps_f(),
					config.with_statistics()));
				opt_ptr = opt_owner.get();
				if (!quiet) cout << "  fd-mode:\t\t" << mode << " (bisector OptimLargestFirst + TempBuffer)\n";
			} else {
				if (!quiet) cerr << "  [warning] fd-mode '" << mode << "' no soportado; usando modo base.\n";
			}
		}

		if (!opt_ptr) {
			opt_owner.reset(new Optimizer(config));
			opt_ptr = opt_owner.get();
		}
		Optimizer& o = *opt_ptr;

		// display solutions with up to 12 decimals
		cout.precision(12);

		if (!quiet)
			cout << "running............" << endl << endl;

		// Search for the optimum
		// Get the solutions
		if (input_file)
			if (initial_loup)
				o.optimize(input_file.Get().c_str(), initial_loup.Get());
			else
				o.optimize(input_file.Get().c_str());
		else
			if (initial_loup)
				o.optimize(sys->box, initial_loup.Get());
			else
				o.optimize(sys->box);

		if (!quiet) {
			cout << "nodes (cells):\t\t" << o.get_nb_cells() << endl;
		}

		if (trace) cout << endl;

		// Report some information (computation time, etc.)

		if (!quiet)
			o.report(); // will include statistics if they are enabled

		if (temp_raw && use_fd_variant) {
			// Línea fácil de parsear para el menú/CSV.
			cout << "fd_triggers:" << temp_raw->trigger_count() << endl;
			cout << "triggers_csv:" << temp_raw->trigger_count() << endl;
		}

		o.get_data().save(output_cov_file.c_str());

		if (!quiet) {
			cout << " results written in " << output_cov_file << "\n";
			if (overwitten)
				cout << " (old file saved in " << cov_copy << ")\n";
		}

		delete sys;
		return 0;

	}
	catch(ibex::UnknownFileException& e) {
		cerr << "Error: cannot read file '" << filename.Get() << "'" << endl;
	}
	catch(ibex::SyntaxError& e) {
		cout << e << endl;
	}
}
