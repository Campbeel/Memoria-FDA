// temp_buffer.h
// Buffer personalizado con coste térmico y triggers de corte (profundidad/volumen).

#pragma once

#include "ibex.h"
#include <cstdint>
#include <vector>

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
        double depth_hard_cut = 0.0;     // si >0 descarta al superar esta profundidad
        double vol_hard_ratio = 0.0;     // si >0 descarta si vol <= vol_hard_ratio * V0_ref
        bool rand_k = false;             // si true, k se perturba con ruido determinístico
        uint64_t rand_seed = 1;          // semilla del ruido
        double tie_noise = 0.0;          // ruido pequeño para desempate en el coste
        double depth_cut_jitter = 0.0;   // variación relativa por nodo en depth_cut
        double vol_cut_jitter = 0.0;     // variación relativa por nodo en vol_ratio_cut
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
    bool reinserting_ = false;
    mutable double last_min_; // para evitar new_uplo descendente (ya no se usa)
    size_t trigger_count_ = 0; // para estadísticas
    bool debug_triggers_ = false;
    mutable size_t debug_shown_ = 0;

public:
    size_t trigger_count() const { return trigger_count_; }
};
