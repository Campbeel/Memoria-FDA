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
}

TempBuffer::TempCost::TempCost(const ExtendedSystem& sys, int goal_var, const Params& p)
: CellCostFunc(sys, false), goal_var(goal_var), params(p) {}

double TempBuffer::TempCost::cost(const Cell& c) const {
    double lb = c.box[goal_var].lb();
    if (!std::isfinite(lb)) lb = 0.0;

    double k_eff = params.k;
    if (params.rand_k) {
        k_eff *= deterministic_noise(params.rand_seed, c.depth);
    }
    double T = params.T0 * std::pow(k_eff/2.0, static_cast<double>(c.depth));

    double score = lb - params.bias * T;

    if (params.depth_cut > 0 && static_cast<int>(c.depth) >= params.depth_cut) {
        score += params.depth_penalty * (1.0 + static_cast<double>(c.depth - params.depth_cut));
    }

    if (params.vol_ratio_cut > 0.0 && !c.box.is_unbounded()) {
        double V = c.box.volume();
        if (std::isfinite(V) && params.V0_ref > 0.0) {
            double thresh = params.vol_ratio_cut * params.V0_ref;
            if (V <= thresh) {
                double gap = (thresh - V) / params.V0_ref;
                score += params.vol_penalty * (1.0 + gap);
            }
        }
    }

    return score;
}

TempBuffer::TempBuffer(const ExtendedSystem& sys,
                       int goal_var,
                       const Params& params)
: cost_obj_(sys, goal_var, params),
  heap_(cost_obj_) {}

TempBuffer::~TempBuffer() {
    flush();
}

void TempBuffer::flush()                   { heap_.flush(); }
unsigned int TempBuffer::size() const      { return static_cast<unsigned int>(heap_.size()); }
bool TempBuffer::empty() const             { return heap_.empty(); }
Cell* TempBuffer::pop()                    { return heap_.pop(); }
Cell* TempBuffer::top() const              { return heap_.top(); }
double TempBuffer::minimum() const         { return heap_.minimum(); }
void TempBuffer::contract(double new_loup) { heap_.contract(new_loup); }

void TempBuffer::push(Cell* cell) {
    if (!cell) return;
    heap_.push(cell);
}

std::ostream& TempBuffer::print(std::ostream& os) const {
    os << "TempBuffer(size=" << size() << ")";
    return os;
}
