// temp_buffer.cpp
// Cola con buckets por lb. Tie-break configurable (profundidad o volumen), sin alterar la cota inferior.

#include "temp_buffer.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>

using namespace ibex;

TempBuffer::TempBuffer(const ExtendedSystem& sys,
                       int goal_var,
                       const Params& params)
: goal_var_(goal_var), params_(params) {
    Cell::set_temp_params(params.k, params.rand_k, params.rand_seed, params.T0);
    debug_triggers_ = std::getenv("FD_TRIGGER_DEBUG") != nullptr;
}

TempBuffer::~TempBuffer() {
    flush();
}

void TempBuffer::flush() {
    buckets_.clear();
    size_ = 0;
    trigger_count_ = depth_trigger_count_ = vol_trigger_count_ = 0;
    vol_eval_count_ = vol_nonfinite_count_ = 0;
}

unsigned int TempBuffer::size() const { return static_cast<unsigned int>(size_); }
bool TempBuffer::empty() const { return size_==0; }
double TempBuffer::minimum() const {
    if (buckets_.empty()) return POS_INFINITY;
    return buckets_.begin()->first;
}
void TempBuffer::contract(double loup) {
    auto it = buckets_.upper_bound(loup);
    while (it != buckets_.end()) {
        size_ -= it->second.size();
        it = buckets_.erase(it);
    }
}

double TempBuffer::log_volume(const IntervalVector& box) const {
    double acc = 0.0;
    for (int i=0;i<box.size();++i) {
        double w = box[i].diam();
        if (!std::isfinite(w)) return POS_INFINITY;
        w = std::max(1e-300, std::min(1e300, w));
        acc += std::log10(w);
    }
    return acc;
}

double TempBuffer::temp_value(const Cell& c) const {
    double T = c.temperature;
    if (!(std::isfinite(T) && T>0.0)) {
        T = params_.T0 * std::pow(params_.k/2.0, static_cast<double>(c.depth));
    }
    return T;
}

size_t TempBuffer::choose_index(const std::vector<Cell*>& vec) const {
    if (vec.empty()) return 0;
    size_t best = 0;
    if (params_.tie_break_mode==1) {
        int best_d = vec[0]->depth;
        double best_t = temp_value(*vec[0]);
        for (size_t i=1;i<vec.size();++i) {
            int d = vec[i]->depth;
            double t = temp_value(*vec[i]);
            if (d < best_d || (d==best_d && t<best_t)) {
                best = i; best_d = d; best_t = t;
            }
        }
    } else if (params_.tie_break_mode==2) {
        double best_v = log_volume(vec[0]->box);
        double best_t = temp_value(*vec[0]);
        for (size_t i=1;i<vec.size();++i) {
            double v = log_volume(vec[i]->box);
            double t = temp_value(*vec[i]);
            if (v < best_v || (v==best_v && t<best_t)) {
                best = i; best_v = v; best_t = t;
            }
        }
    }
    return best;
}

void TempBuffer::push(Cell* cell) {
    if (!cell) return;
    if (cell->depth>0) {
        int depth_limit = params_.depth_cut;
        if (depth_limit > 0 && params_.depth_cut_jitter > 0.0) {
            double n = ((double)std::rand()/RAND_MAX - 0.5);
            depth_limit = std::max(1, static_cast<int>(std::round(depth_limit * (1.0 + params_.depth_cut_jitter * n))));
        }
        bool depth_trigger = (depth_limit > 0 && static_cast<int>(cell->depth) >= depth_limit);
        bool vol_trigger = false;
        if (!cell->box.is_unbounded()) {
            vol_eval_count_++;
            double V = cell->box.volume();
            double logV = log_volume(cell->box);
            bool use_log = params_.use_log_volume || !std::isfinite(V) || !std::isfinite(params_.V0_ref);
            double vr = params_.vol_ratio_cut;
            if (params_.vol_cut_jitter > 0.0) {
                double n = ((double)std::rand()/RAND_MAX - 0.5);
                vr = vr * (1.0 + params_.vol_cut_jitter * n);
            }
            if (use_log) {
                if (std::isfinite(logV) && std::isfinite(params_.log_V0_ref) && vr>0.0) {
                    double log_thresh = params_.log_V0_ref + std::log10(vr);
                    vol_trigger = (logV <= log_thresh);
                } else {
                    vol_nonfinite_count_++;
                }
            } else if (std::isfinite(V) && params_.V0_ref > 0.0 && vr>0.0) {
                double thresh = vr * params_.V0_ref;
                vol_trigger = (V <= thresh);
            } else {
                vol_nonfinite_count_++;
            }
        }
        if (depth_trigger || vol_trigger) {
            trigger_count_++;
            if (depth_trigger) depth_trigger_count_++;
            if (vol_trigger)   vol_trigger_count_++;
            if (debug_triggers_ && debug_shown_ < 5) {
                std::cout << "[fd-debug] trigger depth=" << cell->depth
                          << " depth_cut=" << depth_limit
                          << " vol_trigger=" << vol_trigger
                          << " V0=" << params_.V0_ref
                          << " vr=" << params_.vol_ratio_cut << std::endl;
                debug_shown_++;
            }
        }
    }
    double lb = cell->box[goal_var_].lb();
    if (!std::isfinite(lb)) lb = cell->box[goal_var_].ub();
    if (!std::isfinite(lb)) lb = POS_INFINITY;
    buckets_[lb].push_back(cell);
    size_++;
}

Cell* TempBuffer::pop() {
    if (buckets_.empty()) return NULL;
    auto it = buckets_.begin();
    auto& vec = it->second;
    if (vec.empty()) {
        buckets_.erase(it);
        return pop();
    }
    size_t idx = choose_index(vec);
    Cell* res = vec[idx];
    vec[idx] = vec.back();
    vec.pop_back();
    if (vec.empty()) buckets_.erase(it);
    size_--;
    return res;
}

Cell* TempBuffer::top() const {
    if (buckets_.empty()) return NULL;
    auto it = buckets_.begin();
    const auto& vec = it->second;
    if (vec.empty()) return NULL;
    size_t idx = choose_index(vec);
    return vec[idx];
}

std::ostream& TempBuffer::print(std::ostream& os) const {
    os << "TempBuffer(size=" << size() << ")";
    return os;
}
