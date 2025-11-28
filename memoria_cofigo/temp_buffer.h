// temp_buffer.h
// Buffer personalizado con coste térmico y triggers de corte (profundidad/volumen).

#pragma once

#include "ibex.h"
#include <cstdint>

class TempBuffer : public ibex::CellBufferOptim {
public:
    struct Params {
        double k = 10.0;                 // factor k en T_hijo = k*T_padre/2
        double bias = 1e-3;              // cuánto pesa la temperatura en el coste
        double T0 = 100.0;               // temperatura de la raíz
        int depth_cut = 0;               // si >0 aplica penalización por profundidad
        double vol_ratio_cut = 0.0;      // si >0 aplica penalización por volumen relativo
        double V0_ref = 1.0;             // volumen de referencia
        double depth_penalty = 1e6;      // penalización al superar el corte de profundidad
        double vol_penalty = 1e6;        // penalización al superar el corte de volumen
        bool rand_k = false;             // si true, k se perturba con ruido determinístico
        uint64_t rand_seed = 1;          // semilla del ruido
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
    struct TempCost : public ibex::CellCostFunc {
        TempCost(const ibex::ExtendedSystem& sys, int goal_var, const Params& p);
        double cost(const ibex::Cell& c) const override;

        int goal_var;
        Params params;
    };

    TempCost cost_obj_;
    ibex::Heap<ibex::Cell> heap_;
};
