// temp_buffer.cpp

#include "temp_buffer.h"
#include <cmath>
#include <cstdint>

using namespace ibex;

namespace {
// Pequeño generador determinístico (xorshift-like) para evitar RNG global.
inline double deterministic_noise(uint64_t seed, uint64_t depth) {
    uint64_t z = (depth + 1) ^ (seed + 0x9e3779b97f4a7c15ULL);
    z ^= (z >> 30); z *= 0xbf58476d1ce4e5b9ULL;
    z ^= (z >> 27); z *= 0x94d049bb133111ebULL;
    z ^= (z >> 31);
    const double norm = static_cast<double>((uint64_t(1) << 53) - 1);
    return 0.5 + (z & ((uint64_t(1) << 53) - 1)) / norm; // [0.5 , 1.5)
}

inline double jitter(double base, double rel, uint64_t seed, uint64_t depth) {
    if (rel <= 0.0) return base;
    double n = deterministic_noise(seed, depth) - 1.0; // [-0.5,0.5]
    return base * (1.0 + rel * n);
}
}

TempBuffer::TempCost::TempCost(const ExtendedSystem& sys, int goal_var, const Params& p)
: CellCostFunc(sys, false), goal_var(goal_var), params(p) {}

double TempBuffer::TempCost::cost(const Cell& c) const {
    double lb = c.box[goal_var].lb();
    double ub = c.box[goal_var].ub();
    if (!std::isfinite(lb)) lb = ub;
    if (!std::isfinite(ub)) ub = lb;
    // Usamos lb como coste principal; ub sólo si lb no está definido.
    double bound = std::isfinite(lb) ? lb : ub;

    double k_eff = params.k;
    if (params.rand_k) {
        k_eff *= deterministic_noise(params.rand_seed, c.depth);
    }
    double T = c.temperature;
    if (!(std::isfinite(T) && T>0.0)) {
        T = params.T0 * std::pow(k_eff/2.0, static_cast<double>(c.depth));
    }

    // Coste híbrido: lb principal + tie-break térmico (menor score para T alta).
    double temp_term = params.bias * (1.0 / (1.0 + T));
    double score = bound + temp_term;

    int depth_limit = params.depth_cut;
    if (depth_limit > 0 && params.depth_cut_jitter > 0.0) {
        double jd = jitter(static_cast<double>(depth_limit), params.depth_cut_jitter, params.rand_seed ^ 0x44444444ULL, c.depth);
        depth_limit = std::max(1, static_cast<int>(std::round(jd)));
    }

    bool depth_trigger = false;
    if (depth_limit > 0 && static_cast<int>(c.depth) >= depth_limit && c.depth>0) {
        depth_trigger = true;
    }

    if (params.vol_ratio_cut > 0.0 && !c.box.is_unbounded()) {
        double V = c.box.volume();
        if (std::isfinite(V) && params.V0_ref > 0.0) {
            double vr = params.vol_ratio_cut;
            if (params.vol_cut_jitter > 0.0) {
                vr = jitter(vr, params.vol_cut_jitter, params.rand_seed ^ 0x33333333ULL, c.depth);
            }
            double thresh = vr * params.V0_ref;
            if (V <= thresh && c.depth>0) {
                // solo contar, sin penalizar
            }
        }
    }

    if (params.tie_noise > 0.0) {
        double tie = params.tie_noise * deterministic_noise(params.rand_seed ^ 0x22222222ULL, c.depth);
        score += tie;
    }

    return score;
}

TempBuffer::TempBuffer(const ExtendedSystem& sys,
                       int goal_var,
                       const Params& params)
: cost_obj_(sys, goal_var, params),
  heap_(cost_obj_),
  last_min_(POS_INFINITY) {
    // Configura parámetros de temperatura en la celda para herencia T_hijo = (T_padre/2)*k
    Cell::set_temp_params(params.k, params.rand_k, params.rand_seed, params.T0);
    debug_triggers_ = std::getenv("FD_TRIGGER_DEBUG") != nullptr;
}

TempBuffer::~TempBuffer() {
    flush();
}

void TempBuffer::flush() {
    heap_.flush();
    last_min_ = POS_INFINITY;
}
unsigned int TempBuffer::size() const      { return static_cast<unsigned int>(heap_.size()); }
bool TempBuffer::empty() const             { return heap_.empty(); }
double TempBuffer::minimum() const {
    return heap_.empty() ? POS_INFINITY : heap_.minimum();
}
void TempBuffer::contract(double new_loup) { heap_.contract(new_loup); }

void TempBuffer::push(Cell* cell) {
    if (!cell) return;
    // Contar trigger solo en nodos no raíz (hijos) y sin modificar la cola.
    if (!reinserting_ && cell->depth>0) {
        int depth_limit = cost_obj_.params.depth_cut;
        if (depth_limit > 0 && cost_obj_.params.depth_cut_jitter > 0.0) {
            double jd = jitter(static_cast<double>(depth_limit), cost_obj_.params.depth_cut_jitter, cost_obj_.params.rand_seed ^ 0x44444444ULL, cell->depth);
            depth_limit = std::max(1, static_cast<int>(std::round(jd)));
        }
        bool depth_trigger = (depth_limit > 0 && static_cast<int>(cell->depth) >= depth_limit);
        bool vol_trigger = false;
        if (cost_obj_.params.vol_ratio_cut > 0.0 && !cell->box.is_unbounded()) {
            double V = cell->box.volume();
            if (std::isfinite(V) && cost_obj_.params.V0_ref > 0.0) {
                double vr = cost_obj_.params.vol_ratio_cut;
                if (cost_obj_.params.vol_cut_jitter > 0.0) {
                    vr = jitter(vr, cost_obj_.params.vol_cut_jitter, cost_obj_.params.rand_seed ^ 0x33333333ULL, cell->depth);
                }
                double thresh = vr * cost_obj_.params.V0_ref;
                vol_trigger = (V <= thresh);
            }
        }
        if (depth_trigger || vol_trigger) {
            trigger_count_++;
            if (debug_triggers_ && debug_shown_ < 5) {
                std::cout << "[fd-debug] trigger depth=" << cell->depth
                          << " depth_cut=" << depth_limit
                          << " vol_trigger=" << vol_trigger
                          << " V0=" << cost_obj_.params.V0_ref
                          << " vr=" << cost_obj_.params.vol_ratio_cut << std::endl;
                debug_shown_++;
            }
        }
    }
    heap_.push(cell);
}

std::ostream& TempBuffer::print(std::ostream& os) const {
    os << "TempBuffer(size=" << size() << ")";
    return os;
}

Cell* TempBuffer::pop() {
    if (heap_.empty()) return NULL;
    return heap_.pop();
}

Cell* TempBuffer::top() const {
    if (heap_.empty()) return NULL;
    return heap_.top();
}
