// fda_bridge_common.h
//
// Declaraciones compartidas para las variantes de FDA + IBEX.

#pragma once

#include "ibex.h"
#include <limits>
#include <random>
#include <string>

// Contexto simplificado: usa solo HC4 sobre el sistema original (sin ibexopt extendido).
struct OptContext {
    const ibex::System& sys;
    ibex::CtcHC4 contractor;

    explicit OptContext(const ibex::System& s);
};

// Estadísticas de ejecución (compatibles con el CSV de salida).
struct RunStats {
    int    run_id = 0;
    std::string variant;
    double best_value = std::numeric_limits<double>::infinity();
    long   nodes_visited = 0;
    double elapsed_seconds = 0.0;
    int    max_depth = 0;
    double avg_depth = 0.0;
    bool   reached_optimum = false;
};

// Nodo básico usado por el B&B sencillo.
struct BbNode {
    ibex::IntervalVector box;
    ibex::Interval       goal_bounds; // cota [lb,ub] de f(x) tras contracción
    int depth;
    double T = 0.0; // temperatura heredada (solo aplica a variantes con T)
};

// Ajusta k cuando hay mejora en la mejor cota (determinista).
double adaptive_k(double k_base, double ub_prev, double ub_new);
double heur_diving_score(const OptContext& opt,
                         const ibex::IntervalVector& box,
                         const ibex::Interval& goal_bounds);
bool contract_with_goal(OptContext& opt,
                        ibex::IntervalVector& box,
                        ibex::Interval& goal_bounds,
                        double current_ub = std::numeric_limits<double>::infinity());
void bisect_box(OptContext& opt,
                const ibex::IntervalVector& parent,
                const ibex::Interval& goal_bounds,
                ibex::IntervalVector& left,
                ibex::IntervalVector& right);
double eval_at_mid(const OptContext& opt, const ibex::IntervalVector& box);
bool is_valid_box(const ibex::IntervalVector& box, int expected_dim);

// Interfaces públicas de cada variante.
RunStats run_fda_base(OptContext& opt,
                      const ibex::IntervalVector& root_box,
                      int run_id,
                      double eps_box,
                      int max_iters);

RunStats run_fda_base_bb(OptContext& opt,
                         const ibex::IntervalVector& root_box,
                         int run_id,
                         double eps_box,
                         int max_bb_nodes,
                         int max_iters_fdive);

RunStats run_fda_depth_Tk(OptContext& opt,
                          const ibex::IntervalVector& root_box,
                          int run_id,
                          double T0,
                          double k,
                          int d_max,
                          double eps_box,
                          int max_iters);

RunStats run_fda_depth_Tk_bb(OptContext& opt,
                             const ibex::IntervalVector& root_box,
                             int run_id,
                             double T0,
                             double k,
                             int d_max_fdive,
                             double eps_box,
                             int max_bb_nodes,
                             int max_iters_fdive);

RunStats run_fda_depth_Tk_rand(OptContext& opt,
                               const ibex::IntervalVector& root_box,
                               int run_id,
                               double T0,
                               double k,
                               int d_max,
                               double eps_box,
                               int max_iters,
                               std::mt19937& rng);

RunStats run_fda_depth_Tk_rand_bb(OptContext& opt,
                                  const ibex::IntervalVector& root_box,
                                  int run_id,
                                  double T0,
                                  double k,
                                  int d_max_fdive,
                                  double eps_box,
                                  int max_bb_nodes,
                                  int max_iters_fdive);

RunStats run_fda_vol_Tk(OptContext& opt,
                        const ibex::IntervalVector& root_box,
                        int run_id,
                        double T0,
                        double k,
                        double eps_V,
                        double beta,
                        double eps_box,
                        int max_iters);

RunStats run_fda_vol_Tk_bb(OptContext& opt,
                           const ibex::IntervalVector& root_box,
                           int run_id,
                           double T0,
                           double k,
                           double eps_V,
                           double beta,
                           double eps_box,
                           int max_bb_nodes,
                           int max_iters_fdive);

RunStats run_fda_vol_Tk_rand(OptContext& opt,
                             const ibex::IntervalVector& root_box,
                             int run_id,
                             double T0,
                             double k,
                             double eps_V,
                             double beta,
                             double eps_box,
                             int max_iters,
                             std::mt19937& rng);

RunStats run_fda_vol_Tk_rand_bb(OptContext& opt,
                                const ibex::IntervalVector& root_box,
                                int run_id,
                                double T0,
                                double k,
                                double eps_V,
                                double beta,
                                double eps_box,
                                int max_bb_nodes,
                                int max_iters_fdive);
