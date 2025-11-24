// fda_bridge_common.h
//
// Declaraciones compartidas para las variantes de FDA + IBEX.
// Cada implementación vive en su propio .cpp para aislar cambios.

#pragma once

#include "ibex.h"
#include <limits>
#include <random>
#include <string>

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
    int depth;
    double T = 0.0; // temperatura heredada (solo aplica a variantes con T)
};

// Utilidades generales sobre el problema IBEX.
double heur_diving_score(const ibex::System& sys, const ibex::IntervalVector& box);
int widest_var(const ibex::IntervalVector& box);
void bisect_box(const ibex::IntervalVector& parent,
                ibex::IntervalVector& left,
                ibex::IntervalVector& right);
double eval_at_mid(const ibex::System& sys, const ibex::IntervalVector& box);
bool is_valid_box(const ibex::IntervalVector& box, int expected_dim);

// Interfaces públicas de cada variante.
RunStats run_fda_base(const ibex::System& sys,
                      ibex::Ctc& contractor,
                      const ibex::IntervalVector& root_box,
                      int run_id,
                      double eps_box,
                      int max_iters);

RunStats run_fda_base_bb(const ibex::System& sys,
                         ibex::Ctc& contractor,
                         const ibex::IntervalVector& root_box,
                         int run_id,
                         double eps_box,
                         int max_bb_nodes,
                         int max_iters_fdive);

RunStats run_fda_depth_Tk(const ibex::System& sys,
                          ibex::Ctc& contractor,
                          const ibex::IntervalVector& root_box,
                          int run_id,
                          double T0,
                          double k,
                          int d_max,
                          double eps_box,
                          int max_iters);

RunStats run_fda_depth_Tk_bb(const ibex::System& sys,
                             ibex::Ctc& contractor,
                             const ibex::IntervalVector& root_box,
                             int run_id,
                             double T0,
                             double k,
                             int d_max_fdive,
                             double eps_box,
                             int max_bb_nodes,
                             int max_iters_fdive);

RunStats run_fda_depth_Tk_rand(const ibex::System& sys,
                               ibex::Ctc& contractor,
                               const ibex::IntervalVector& root_box,
                               int run_id,
                               double T0,
                               double k,
                               int d_max,
                               double eps_box,
                               int max_iters,
                               std::mt19937& rng);

RunStats run_fda_depth_Tk_rand_bb(const ibex::System& sys,
                                  ibex::Ctc& contractor,
                                  const ibex::IntervalVector& root_box,
                                  int run_id,
                                  double T0,
                                  double k,
                                  int d_max_fdive,
                                  double eps_box,
                                  int max_bb_nodes,
                                  int max_iters_fdive);

RunStats run_fda_vol_Tk(const ibex::System& sys,
                        ibex::Ctc& contractor,
                        const ibex::IntervalVector& root_box,
                        int run_id,
                        double T0,
                        double k,
                        double eps_V,
                        double beta,
                        double eps_box,
                        int max_iters);

RunStats run_fda_vol_Tk_bb(const ibex::System& sys,
                           ibex::Ctc& contractor,
                           const ibex::IntervalVector& root_box,
                           int run_id,
                           double T0,
                           double k,
                           double eps_V,
                           double beta,
                           double eps_box,
                           int max_bb_nodes,
                           int max_iters_fdive);

RunStats run_fda_vol_Tk_rand(const ibex::System& sys,
                             ibex::Ctc& contractor,
                             const ibex::IntervalVector& root_box,
                             int run_id,
                             double T0,
                             double k,
                             double eps_V,
                             double beta,
                             double eps_box,
                             int max_iters,
                             std::mt19937& rng);

RunStats run_fda_vol_Tk_rand_bb(const ibex::System& sys,
                                ibex::Ctc& contractor,
                                const ibex::IntervalVector& root_box,
                                int run_id,
                                double T0,
                                double k,
                                double eps_V,
                                double beta,
                                double eps_box,
                                int max_bb_nodes,
                                int max_iters_fdive);
