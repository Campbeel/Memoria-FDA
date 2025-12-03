// temp_buffer.h
// Buffer con selección térmica: balancea la exploración usando la temperatura
// del nodo, pero conserva los límites inferiores reales para uplo/contracción.

#pragma once

#include "ibex.h"
#include <cstdint>
#include <functional>
#include <memory>
#include <queue>
#include <random>
#include <utility>
#include <vector>

class TempBuffer : public ibex::CellBufferOptim {
public:
    struct Params {
        double k = 10.0;
        double bias = 1e-3;
        double T0 = 100.0;
        int depth_cut = 0;
        double vol_ratio_cut = 0.0;
        double V0_ref = 1.0;
        double log_V0_ref = 0.0;
        bool use_log_volume = false;
        double depth_penalty = 0.0;
        double vol_penalty = 0.0;
        double depth_hard_cut = 0.0;
        double vol_hard_ratio = 0.0;
        bool rand_k = false;
        uint64_t rand_seed = 1;
        double tie_noise = 0.0;
        double depth_cut_jitter = 0.0;
        double vol_cut_jitter = 0.0;
        int tie_break_mode = 0; // se ignora; se mantiene por compatibilidad
    };

    TempBuffer(const ibex::ExtendedSystem& sys,
               int goal_var,
               const Params& params,
               ibex::CellBufferOptim& delegate);
    virtual ~TempBuffer();

    void flush() override;
    unsigned int size() const override;
    bool empty() const override;
    void push(ibex::Cell* cell) override;
    ibex::Cell* pop() override;
    ibex::Cell* top() const override;
    double minimum() const override;
    void contract(double loup) override;
    std::ostream& print(std::ostream& os) const override;

private:
    ibex::CellBufferOptim& delegate_;
    int goal_var_;
    Params params_;
    std::vector<std::unique_ptr<struct Item>> items_;
    using HeapEntry = std::pair<double, size_t>;
    mutable std::priority_queue<HeapEntry, std::vector<HeapEntry>, std::greater<HeapEntry>> score_heap_;
    mutable std::priority_queue<HeapEntry, std::vector<HeapEntry>, std::greater<HeapEntry>> lb_heap_;
    size_t alive_count_ = 0;
    unsigned int depth_floor_ = 0;
    double current_logV_ref_ = ibex::POS_INFINITY;
    size_t trigger_count_ = 0;
    size_t depth_trigger_count_ = 0;
    size_t vol_trigger_count_ = 0;
    size_t vol_eval_count_ = 0;
    size_t vol_nonfinite_count_ = 0;
    bool debug_triggers_ = false;
    mutable size_t debug_shown_ = 0;
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> noise_dist_{-0.5, 0.5};

    struct Item {
        ibex::Cell* cell = nullptr;
        double lb = ibex::POS_INFINITY;
        double score = ibex::POS_INFINITY;
        double vol_ratio = ibex::POS_INFINITY;
        double log_volume = ibex::POS_INFINITY;
        unsigned int depth = 0;
        bool alive = true;
    };

    double rand_unit() { return noise_dist_(rng_); }
    double log_volume(const ibex::IntervalVector& box) const;
    double compute_score(const Item& item);
    template <typename Heap>
    void prune_heap(Heap& h) const;
    void prune_heaps() const;
    void drop_item(Item& item);

public:
    size_t trigger_count() const { return trigger_count_; }
    size_t depth_trigger_count() const { return depth_trigger_count_; }
    size_t vol_trigger_count() const { return vol_trigger_count_; }
    size_t vol_eval_count() const { return vol_eval_count_; }
    size_t vol_nonfinite_count() const { return vol_nonfinite_count_; }
};
