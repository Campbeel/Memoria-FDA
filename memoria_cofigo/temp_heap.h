// temp_heap.h
// Heap/Buffer con coste térmico y triggers de corte (profundidad/volumen).
#pragma once

#include "ibex.h"
#include <memory>

class TempHeap : public ibex::Heap<ibex::Cell>, public ibex::CellBufferOptim {
public:
    TempHeap(const ibex::ExtendedSystem& sys,
             int goal_var,
             double k_factor,
             double bias,
             int depth_cut,
             double vol_ratio_cut,
             double V0_ref);
    virtual ~TempHeap();

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
    struct TempCost : public ibex::CellCostFunc {
        TempCost(const ibex::ExtendedSystem& sys, int goal_var, double k, double bias)
            : ibex::CellCostFunc(sys, false), goal_var(goal_var), k(k), bias(bias), T0(100.0) {}
        double cost(const ibex::Cell& c) const override {
            double lb = c.box[goal_var].lb();
            if (!std::isfinite(lb)) lb = 0.0;
            double T = T0 * std::pow(k/2.0, static_cast<double>(c.depth));
            return lb - bias * T;
        }
        int goal_var;
        double k;
        double bias;
        double T0;
    };

    int goal_var_;
    int depth_cut_;
    double vol_ratio_cut_;
    double V0_;
    std::unique_ptr<TempCost> cost_owner_;
};
