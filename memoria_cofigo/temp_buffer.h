// temp_buffer.h
// Cola simple con buckets por lb: cuenta triggers y desempata por profundidad/volumen dentro del mismo lb.
// minimum() devuelve el lb mínimo, contract() elimina buckets > loup. No toca la cota inferior.

#pragma once

#include "ibex.h"
#include <cstdint>
#include <map>
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
        int tie_break_mode = 0; // 0 ninguno, 1 profundidad asc+T, 2 volumen asc+T
    };

    TempBuffer(const ibex::ExtendedSystem& sys,
               int goal_var,
               const Params& params);
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
    int goal_var_;
    Params params_;
    std::map<double, std::vector<ibex::Cell*>> buckets_; // lb -> celdas
    size_t size_ = 0;
    size_t trigger_count_ = 0;
    size_t depth_trigger_count_ = 0;
    size_t vol_trigger_count_ = 0;
    size_t vol_eval_count_ = 0;
    size_t vol_nonfinite_count_ = 0;
    bool debug_triggers_ = false;
    mutable size_t debug_shown_ = 0;

    double log_volume(const ibex::IntervalVector& box) const;
    double temp_value(const ibex::Cell& c) const;
    size_t choose_index(const std::vector<ibex::Cell*>& vec) const;

public:
    size_t trigger_count() const { return trigger_count_; }
    size_t depth_trigger_count() const { return depth_trigger_count_; }
    size_t vol_trigger_count() const { return vol_trigger_count_; }
    size_t vol_eval_count() const { return vol_eval_count_; }
    size_t vol_nonfinite_count() const { return vol_nonfinite_count_; }
};
