////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2015-2017 Clemson University.
//
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file hydro_helper.h
//! \brief A header file to define useful hydro helper routines.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef MY_HYDRO_HELPER_H
#define MY_HYDRO_HELPER_H

#include <boost/format.hpp>

#include <boost/program_options.hpp>

#include "nnt/iter.h"
#include "nnt/string_defs.h"

#include "user/evolve.h"
#include "user/hydro_helper.h"

#define S_DELTA_TRAJ    "delta"
#define S_RHO_1         "rho_1"
#define S_RHO_2         "rho_2"
#define S_ROOT_FACTOR   "root_factor"

namespace po = boost::program_options;

/**
 * @brief A namespace for user-defined functions.
 */
namespace my_user
{

typedef std::map<std::string, boost::any> param_map_t;

typedef std::vector< double > state_type;

//##############################################################################
// Prototypes.
//##############################################################################

void
get_user_defined_descriptions( po::options_description& );

void
set_user_defined_options( po::variables_map&, param_map_t& );

void
initialize_state( param_map_t&, state_type& );

double
acceleration( param_map_t&, nnt::Zone&, const state_type&, const double );

double rho_function( param_map_t&, const state_type& );

double t9_function( nnt::Zone& zone, param_map_t&, Libnucnet__NetView * );

void
observer_function( nnt::Zone&,
  const state_type&,
  const state_type&,
  const double
);

} // namespace my_user

#endif // MY_HYDRO_HELPER_H
