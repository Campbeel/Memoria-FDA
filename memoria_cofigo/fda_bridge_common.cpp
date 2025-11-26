// fda_bridge_common.cpp
//
// Implementación de utilidades compartidas.

#include "fda_bridge_common.h"
#include <cstdlib>

using namespace ibex;
using namespace std;

namespace {
int widest_var(const IntervalVector& box) {
    int idx = 0;
    double wmax = box[0].diam();
    for (int i = 1; i < box.size(); ++i) {
        double w = box[i].diam();
        if (w > wmax) { wmax = w; idx = i; }
    }
    return idx;
}
}

OptContext::OptContext(const System& s)
    : sys(s),
      contractor(s) {} // HC4 contractor sobre el sistema original

double adaptive_k(double k_base, double ub_prev, double ub_new) {
    if (!std::isfinite(ub_prev) || !std::isfinite(ub_new) || ub_new >= ub_prev) {
        return k_base;
    }
    double rel_impr = (ub_prev - ub_new) / (std::fabs(ub_prev) + 1e-9);
    double alpha = 0.5;
    return k_base * (1.0 + alpha * rel_impr);
}

double heur_diving_score(const OptContext& opt,
                         const IntervalVector& box,
                         const Interval& goal_bounds) {
    if (!opt.sys.goal) return 0.0;
    if (!goal_bounds.is_empty()) {
        double lb = goal_bounds.lb();
        if (std::isfinite(lb)) return -lb;
    }
    Interval f_int = opt.sys.goal->eval(box);
    return -f_int.lb();
}

bool contract_with_goal(OptContext& opt,
                        IntervalVector& box,
                        Interval& goal_bounds,
                        double /*current_ub*/) {
    opt.contractor.contract(box);
    if (box.is_empty()) {
        goal_bounds.set_empty();
        return false;
    }
    if (opt.sys.goal) {
        goal_bounds = opt.sys.goal->eval(box);
    } else {
        goal_bounds = Interval();
    }
    return true;
}

void bisect_box(OptContext& /*opt*/,
                const IntervalVector& parent,
                const Interval& /*goal_bounds*/,
                IntervalVector& left,
                IntervalVector& right) {
    left  = parent;
    right = parent;
    int var = widest_var(parent);
    double mid = 0.5 * (parent[var].lb() + parent[var].ub());
    left[var]  = Interval(parent[var].lb(), mid);
    right[var] = Interval(mid, parent[var].ub());
}

double eval_at_mid(const OptContext& opt, const IntervalVector& box) {
    if (!opt.sys.goal) return std::numeric_limits<double>::infinity();
    if (box.is_unbounded()) {
        return std::numeric_limits<double>::infinity();
    }
    Vector mid = box.mid();
    Interval val = opt.sys.goal->eval(mid);
    return 0.5 * (val.lb() + val.ub());
}

bool is_valid_box(const IntervalVector& box, int expected_dim) {
    if (box.size() != expected_dim) return false;
    for (int i = 0; i < box.size(); ++i) {
        double lb = box[i].lb();
        double ub = box[i].ub();
        if (!std::isfinite(lb) || !std::isfinite(ub) || lb > ub) {
            return false;
        }
    }
    return true;
}

bool parse_int_arg(const char* text, int& value) {
    if (!text) return false;
    char* endptr = nullptr;
    long parsed = std::strtol(text, &endptr, 10);
    if (endptr == text || *endptr != '\0') {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}
