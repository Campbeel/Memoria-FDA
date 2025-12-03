//============================================================================
//                                  I B E X                                   
// File        : ibex_Cell.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 10, 2012
// Last Update : Jun 07, 2018
//============================================================================

#include "ibex_Cell.h"
#include "ibex_Bsc.h"
#include <limits.h>
#include "ibex_Bxp.h"
#include "ibex_Bxp.h"
#include <cmath>

using namespace std;

namespace ibex {

static double deterministic_noise(uint64_t seed, uint64_t depth) {
	uint64_t z = (depth + 1) ^ (seed + 0x9e3779b97f4a7c15ULL);
	z ^= (z >> 30); z *= 0xbf58476d1ce4e5b9ULL;
	z ^= (z >> 27); z *= 0x94d049bb133111ebULL;
	z ^= (z >> 31);
	const double norm = static_cast<double>((uint64_t(1) << 53) - 1);
	return 0.5 + (z & ((uint64_t(1) << 53) - 1)) / norm; // [0.5 , 1.5)
}

Cell::Cell(const IntervalVector& box, int var, unsigned int depth) :
	box(box), prop(this->box), bisected_var(var), depth(depth), temperature(temp_T0_) {

}

Cell::Cell(const Cell& e) :
	box(e.box), prop(this->box, e.prop), bisected_var(e.bisected_var), depth(e.depth), temperature(e.temperature) {

}

pair<Cell*,Cell*> Cell::bisect(const BisectionPoint& pt) const {

	Cell* cleft;
	Cell* cright;

	if (pt.rel_pos) {
		pair<IntervalVector,IntervalVector> boxes=box.bisect(pt.var,pt.pos);
		cleft = new Cell(boxes.first, pt.var, depth+1);
		cright = new Cell(boxes.second, pt.var, depth+1);
	} else {
		IntervalVector b1(box);
		IntervalVector b2(box);
		b1[pt.var]=Interval(box[pt.var].lb(), pt.pos);
		b2[pt.var]=Interval(pt.pos, box[pt.var].ub());
		cleft = new Cell(b1, pt.var, depth+1);
		cright = new Cell(b2, pt.var, depth+1);
	}

	double noise = temp_rand_ ? deterministic_noise(temp_seed_, depth+1) : 1.0;
	double child_temp = (temperature/2.0) * temp_k_ * noise;
	cleft->temperature = child_temp;
	cright->temperature = child_temp;

	prop.update_bisect(Bisection(box, pt, cleft->box, cright->box), cleft->prop, cright->prop);

	return pair<Cell*,Cell*>(cleft,cright);
}

Cell::~Cell() {

}

std::ostream& operator<<(std::ostream& os, const Cell& c) {
	os << c.box;
	return os;
}

} // end namespace ibex
