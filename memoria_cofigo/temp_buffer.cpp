// temp_buffer.cpp
// Buffer con score térmico para seleccionar nodos sin romper uplo/contracción.

#include "temp_buffer.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <climits>

using namespace ibex;

TempBuffer::TempBuffer(const ExtendedSystem& sys,
                       int goal_var,
                       const Params& params,
                       CellBufferOptim& delegate)
: delegate_(delegate), goal_var_(goal_var), params_(params) {
    (void) sys;
    Cell::set_temp_params(params.k, params.rand_k, params.rand_seed, params.T0);
    debug_triggers_ = std::getenv("FD_TRIGGER_DEBUG") != nullptr;
    uint64_t seed = params.rand_seed ? params.rand_seed : static_cast<uint64_t>(std::random_device{}());
    rng_.seed(seed);
    current_logV_ref_ = params.log_V0_ref;
}

TempBuffer::~TempBuffer() {
    flush();
}

void TempBuffer::flush() {
    // Limpia nuestras estructuras y el delegado por si acumula estado.
    for (auto& ptr : items_) {
        if (ptr && ptr->alive && ptr->cell) {
            delete ptr->cell;
        }
    }
    items_.clear();
    score_heap_ = {};
    lb_heap_ = {};
    alive_count_ = 0;
    depth_floor_ = 0;
    current_logV_ref_ = params.log_V0_ref;
    trigger_count_ = depth_trigger_count_ = vol_trigger_count_ = 0;
    vol_eval_count_ = vol_nonfinite_count_ = 0;
    delegate_.flush();
}

unsigned int TempBuffer::size() const { return static_cast<unsigned int>(alive_count_); }
bool TempBuffer::empty() const { return alive_count_ == 0; }

double TempBuffer::log_volume(const IntervalVector& box) const {
    double acc = 0.0;
    for (int i = 0; i < box.size(); ++i) {
        double w = box[i].diam();
        if (!std::isfinite(w) || w <= 0.0) return POS_INFINITY;
        w = std::max(1e-300, std::min(1e300, w));
        acc += std::log10(w);
    }
    return acc;
}

double TempBuffer::compute_score(const Item& item) {
    double score = item.lb;

    // Sesgo térmico: preferir nodos "calientes".
    if (std::isfinite(item.cell->temperature) && params_.bias != 0.0) {
        score -= params_.bias * item.cell->temperature;
    }

    if (params_.depth_penalty > 0.0) {
        score += params_.depth_penalty * static_cast<double>(item.depth);
    }

    if (params_.vol_penalty > 0.0 && std::isfinite(item.vol_ratio)) {
        double vterm = params_.use_log_volume
                           ? (item.log_volume - params_.log_V0_ref)
                           : std::log(std::max(1e-300, item.vol_ratio));
        score += params_.vol_penalty * vterm;
    }

    if (params_.tie_noise > 0.0) {
        score += params_.tie_noise * rand_unit();
    }

    if (!std::isfinite(score)) score = item.lb;
    return score;
}

template <typename Heap>
void TempBuffer::prune_heap(Heap& h) const {
    while (!h.empty()) {
        size_t idx = h.top().second;
        if (idx < items_.size() && items_[idx] && items_[idx]->alive) break;
        h.pop();
    }
}

void TempBuffer::prune_heaps() const {
    prune_heap(score_heap_);
    prune_heap(lb_heap_);
}

void TempBuffer::drop_item(Item& item) {
    if (!item.alive) return;
    item.alive = false;
    if (item.cell) delete item.cell;
    if (alive_count_ > 0) alive_count_--;
}

void TempBuffer::push(Cell* cell) {
    if (!cell) return;

    Item item;
    item.cell = cell;
    item.depth = cell->depth;
    item.lb = cell->box[goal_var_].lb();
    if (std::isnan(item.lb)) item.lb = POS_INFINITY;

    bool depth_trigger = false;
    bool vol_trigger = false;

    if (cell->depth > 0) {
        int depth_limit = params_.depth_cut;
        if (depth_limit > 0 && params_.depth_cut_jitter > 0.0) {
            double n = rand_unit();
            depth_limit = std::max(1, static_cast<int>(std::round(depth_limit * (1.0 + params_.depth_cut_jitter * n))));
        }
        depth_trigger = (depth_limit > 0 && static_cast<int>(cell->depth) >= depth_limit);

        if (!cell->box.is_unbounded()) {
            vol_eval_count_++;
            // Referencia local: el último nodo seleccionado (current_logV_ref_).
            double log_ref = current_logV_ref_;
            if (!std::isfinite(log_ref)) log_ref = params_.log_V0_ref;
            if (!std::isfinite(log_ref)) log_ref = item.log_volume;
            item.log_volume = log_volume(cell->box);
            double vr = params_.vol_ratio_cut;
            if (params_.vol_cut_jitter > 0.0) {
                double n = rand_unit();
                vr = vr * (1.0 + params_.vol_cut_jitter * n);
            }
            bool use_log = true;
            if (use_log) {
                if (std::isfinite(item.log_volume) && std::isfinite(log_ref) && vr > 0.0) {
                    double log_thresh = log_ref + std::log10(vr);
                    vol_trigger = (item.log_volume <= log_thresh);
                } else {
                    vol_nonfinite_count_++;
                }
            }
        }
    }

    if (depth_trigger || vol_trigger) {
        trigger_count_++;
        if (depth_trigger) depth_trigger_count_++;
        if (vol_trigger) vol_trigger_count_++;
        if (debug_triggers_ && debug_shown_ < 5) {
            std::cout << "[fd-debug] trigger depth=" << cell->depth
                      << " depth_cut=" << params_.depth_cut
                      << " vol_trigger=" << vol_trigger
                      << " V0=" << params_.V0_ref
                      << " vr=" << params_.vol_ratio_cut << std::endl;
            debug_shown_++;
        }
    }

    item.score = compute_score(item);

    size_t idx = items_.size();
    items_.push_back(std::make_unique<Item>(item));
    score_heap_.push({item.score, idx});
    lb_heap_.push({item.lb, idx});
    alive_count_++;
}

Cell* TempBuffer::pop() {
    prune_heap(score_heap_);
    if (score_heap_.empty()) return NULL;

    double depth_limit = params_.depth_cut > 0 ? static_cast<double>(depth_floor_ + params_.depth_cut) : POS_INFINITY;
    std::vector<HeapEntry> deferred;
    size_t attempts = 0;
    size_t selected_idx = static_cast<size_t>(-1);

    auto accept = [&](size_t idx)->bool {
        if (idx >= items_.size() || !items_[idx] || !items_[idx]->alive) return false;
        if (params_.depth_cut>0 && items_[idx]->depth > depth_limit) return false;
        return true;
    };

    while (!score_heap_.empty() && attempts < alive_count_) {
        auto top = score_heap_.top(); score_heap_.pop();
        size_t idx = top.second;
        attempts++;
        if (accept(idx)) { selected_idx = idx; break; }
        if (idx < items_.size() && items_[idx] && items_[idx]->alive) deferred.push_back(top);
    }

    if (selected_idx == static_cast<size_t>(-1)) {
        // Todos superan depth_limit: elevar el piso y aceptar el mejor diferido.
        unsigned int min_depth = UINT_MAX;
        for (auto& e : deferred) {
            size_t idx = e.second;
            if (idx < items_.size() && items_[idx]) {
                min_depth = std::min(min_depth, items_[idx]->depth);
            }
        }
        if (min_depth!=UINT_MAX) depth_floor_ = min_depth;
        if (!deferred.empty()) {
            selected_idx = deferred.front().second;
            deferred.erase(deferred.begin());
        } else if (!score_heap_.empty()) {
            selected_idx = score_heap_.top().second;
            score_heap_.pop();
        }
    }

    // Reinsertar diferidos que no fueron elegidos
    for (auto& e : deferred) score_heap_.push(e);

    if (selected_idx == static_cast<size_t>(-1) || selected_idx >= items_.size() || !items_[selected_idx]) return NULL;
    Item& item = *items_[selected_idx];
    if (!item.alive) return NULL;
    item.alive = false;
    prune_heap(lb_heap_);
    if (alive_count_ > 0) alive_count_--;
    // Actualizar referencia de volumen para hijos de este nodo.
    current_logV_ref_ = item.log_volume;
    return item.cell;
}

Cell* TempBuffer::top() const {
    prune_heap(score_heap_);
    if (score_heap_.empty()) return NULL;
    size_t idx = score_heap_.top().second;
    if (idx >= items_.size() || !items_[idx]) return NULL;
    return items_[idx]->cell;
}

double TempBuffer::minimum() const {
    if (alive_count_ == 0) return POS_INFINITY;
    prune_heap(lb_heap_);
    if (lb_heap_.empty()) return POS_INFINITY;
    return lb_heap_.top().first;
}

void TempBuffer::contract(double loup) {
    for (auto& ptr : items_) {
        if (!ptr || !ptr->alive) continue;
        if (ptr->lb > loup) {
            drop_item(*ptr);
        }
    }
    prune_heaps();
}

std::ostream& TempBuffer::print(std::ostream& os) const {
    os << "TempBuffer(size=" << size() << ")";
    return os;
}
