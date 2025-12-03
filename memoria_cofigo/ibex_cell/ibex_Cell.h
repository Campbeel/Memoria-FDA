//============================================================================
//                                  I B E X                                   
// File        : ibex_Cell.h
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 10, 2012
// Last Update : Jun 07, 2018
//============================================================================

#ifndef __IBEX_CELL_H__
#define __IBEX_CELL_H__

#include "ibex_IntervalVector.h"
#include "ibex_BoxProperties.h"
#include "ibex_BisectionPoint.h"
#include "ibex_Map.h"
#include <cstdint>

namespace ibex {

/**
 * \defgroup strategy Strategies
 */

/**
 * \ingroup strategy
 *
 * \brief Node in an interval exploration binary tree.
 *
 * This representation includes default data (current box) and data required by
 * specific operators (contractors, bisectors, cell buffers, ...), stored in
 * a "property list" (see #ibex::BoxProperties). A different cell is associated
 * to each node and the way the properties are inherited from a cell to its children
 * can be controlled thanks to the "update_bisect" function of Bxp. (see #ibex::Bxp).
 *
 * The cell on its own contains the minimum of information associated to the actual
 * search space: the current box and the last bisected variable. Other fields might
 * be added with future releases.
 *
 */
class Cell {
public:
	static inline double temp_k_ = 10.0;
	static inline double temp_T0_ = 100.0;
	static inline bool   temp_rand_ = false;
	static inline uint64_t temp_seed_ = 1;

	static inline void set_temp_params(double k, bool rand_flag, uint64_t seed, double T0) {
		temp_k_ = k;
		temp_rand_ = rand_flag;
		temp_seed_ = seed;
		temp_T0_ = T0;
	}
	/**
	 * \brief Create the root cell.
	 *
	 * \param box - Box (passed by copy).
	 */
	explicit Cell(const IntervalVector& box, int bisected_var=-1, unsigned int depth=0);

	/**
	 * \brief Constructor by copy.
	 */
	explicit Cell(const Cell& e);

	/**
	 * \brief Bisect this cell.
	 *
	 * The box of the first (resp. second) cell is \a left (resp. \a right).
	 * Each sub-cell inherits from the properties of this cell via the
	 * \link #ibex::Bxp::update(const Bisection&, const BoxProperties&) update_bisect \endlink
	 * function.
	 *
	 * <p>
	 * This function is called by the bisector of boxes (see #ibex::Bsc).
	 */
	std::pair<Cell*,Cell*> bisect(const BisectionPoint& b) const;

	/**
	 * \brief Delete *this.
	 */
	virtual ~Cell();

	/**
	 * \brief The box
	 */
	IntervalVector box;

	/**
	 * \brief Properties bound to the box.
	 *
	 * Box properties are transmitted to operators (contractors, etc.)
	 * through the "trust chain" principle (see documentation).
	 */
	BoxProperties prop;

	/**
	 * Last bisected variable (-1 if root node).
	 */
	int bisected_var;

	/**
	 * Cell depth (0 if root node).
	 */
	unsigned int depth;

	/**
	 * Temperature carried by the node (para variantes FDA).
	 */
	double temperature;

private:
};

/**
 * \brief Print the cell.
 */
std::ostream& operator<<(std::ostream& os, const Cell& c);

} // end namespace ibex

#endif // __IBEX_CELL_H__
