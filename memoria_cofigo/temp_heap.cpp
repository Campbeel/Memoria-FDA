// temp_heap.cpp

#include "temp_heap.h"
#include <cmath>

using namespace ibex;

TempHeap::TempHeap(const ExtendedSystem& sys,
                   int goal_var,
                   double k_factor,
                   double bias,
                   int depth_cut,
                   double vol_ratio_cut,
                   double V0_ref)
: Heap<Cell>(*(new TempCost(sys, goal_var, k_factor, bias))),
  CellBufferOptim(),
  goal_var_(goal_var),
  depth_cut_(depth_cut),
  vol_ratio_cut_(vol_ratio_cut),
  V0_(V0_ref > 0 ? V0_ref : 1.0) {
    cost_owner_.reset((TempCost*)&costf);
}

TempHeap::~TempHeap() {
    flush();
    TempCost* c = cost_owner_.release();
    delete c;
}

void TempHeap::flush()                   { Heap<Cell>::flush(); }
unsigned int TempHeap::size() const      { return Heap<Cell>::size(); }
bool TempHeap::empty() const             { return Heap<Cell>::empty(); }
Cell* TempHeap::pop()                    { return Heap<Cell>::pop(); }
Cell* TempHeap::top() const              { return Heap<Cell>::top(); }
double TempHeap::minimum() const         { return Heap<Cell>::minimum(); }
void TempHeap::contract(double new_loup) { Heap<Cell>::contract(new_loup); }

void TempHeap::push(Cell* cell) {
    if (!cell) return;
    // Triggers: si supera profundidad o volumen pequeño relativo a V0, se descarta.
    if (depth_cut_ > 0 && static_cast<int>(cell->depth) >= depth_cut_) {
        delete cell;
        return;
    }
    if (vol_ratio_cut_ > 0.0 && !cell->box.is_unbounded()) {
        double V = cell->box.volume();
        if (std::isfinite(V) && V <= vol_ratio_cut_ * V0_) {
            delete cell;
            return;
        }
    }
    Heap<Cell>::push(cell);
}

std::ostream& TempHeap::print(std::ostream& os) const {
    os << "TempHeap(size=" << size() << ")";
    return os;
}
