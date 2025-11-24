// fda_bridge_common.cpp
//
// Implementación de utilidades compartidas.

#include "fda_bridge_common.h"

using namespace ibex;
using namespace std;

double heur_diving_score(const System& sys, const IntervalVector& box) {
    if (!sys.goal) return 0.0;
    Interval f_int = sys.goal->eval(box);
    return -f_int.lb();
}

int widest_var(const IntervalVector& box) {
    int idx = 0;
    double wmax = box[0].diam();
    for (int i = 1; i < box.size(); ++i) {
        double w = box[i].diam();
        if (w > wmax) { wmax = w; idx = i; }
    }
    return idx;
}

void bisect_box(const IntervalVector& parent,
                IntervalVector& left,
                IntervalVector& right) {
    left  = parent;
    right = parent;
    int var = widest_var(parent);
    double mid = 0.5 * (parent[var].lb() + parent[var].ub());
    left[var]  = Interval(parent[var].lb(), mid);
    right[var] = Interval(mid, parent[var].ub());
}

double eval_at_mid(const System& sys, const IntervalVector& box) {
    if (!sys.goal) return std::numeric_limits<double>::infinity();
    if (box.is_unbounded()) {
        return std::numeric_limits<double>::infinity();
    }
    Vector mid = box.mid();
    Interval val = sys.goal->eval(mid);
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
