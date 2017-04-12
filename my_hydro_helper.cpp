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
//! \file hydro_helper.cpp
//! \brief A file to define useful hydro helper routines.
//!
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include "my_hydro_helper.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace my_user
{

typedef std::vector< double > state_type;

//##############################################################################
// get_user_defined_descriptions().
//##############################################################################

void
get_user_defined_descriptions( po::options_description& user )
{

  try
  {

    user.add_options()

      ( nnt::s_T9_0, po::value<double>()->default_value( 10., "10." ),
        "Initial T (in 10^9 K)"
      )

      ( nnt::s_RHO_0, po::value<double>()->default_value( 1.e8, "1.e8" ),
        "Initial density (g/cc)"
      )
      
      ( S_RHO_1, po::value<double>()->default_value( 9.e7, "9.e7" ),
        "rho_1 density (g/cc)"
      )

      ( nnt::s_TAU, po::value<double>()->default_value( 0.1, "0.1" ),
        "Expansion timescale (s)"
      )
      
      ( S_DELTA_TRAJ, po::value<double>()->default_value( 0.1, "0.1" ),
        "Cutoff time (s)"
      )

      ( S_ROOT_FACTOR, po::value<double>()->default_value( 1.001, "1.001" ),
        "Root expansion factor"
      )

    ;

// Add checks on input.

  }
  catch( std::exception& e )
  {
    std::cerr << "Error: " << e.what() << "\n";
    exit( EXIT_FAILURE );
  }
  catch(...)
  {
    std::cerr << "Exception of unknown type!\n";
    exit( EXIT_FAILURE );
  }
}

//##############################################################################
// set_user_defined_options().
//##############################################################################

void
set_user_defined_options( po::variables_map& vmap, param_map_t& param_map )
{

  param_map[nnt::s_T9_0] = vmap[nnt::s_T9_0].as<double>();
  param_map[nnt::s_RHO_0] = vmap[nnt::s_RHO_0].as<double>();
  param_map[S_RHO_1] = vmap[S_RHO_1].as<double>();
  param_map[nnt::s_TAU] = vmap[nnt::s_TAU].as<double>();
  param_map[S_DELTA_TRAJ] = vmap[S_DELTA_TRAJ].as<double>();
  param_map[S_ROOT_FACTOR] = vmap[S_ROOT_FACTOR].as<double>();
  
  if(
    boost::any_cast<double>( param_map[S_RHO_1] ) >
    boost::any_cast<double>( param_map[nnt::s_RHO_0] )
  )
  {
    std::cerr << "rho_1 must be less than rho_0." << std::endl;
    exit( EXIT_FAILURE );
  }
  
  param_map[S_RHO_2] =
    boost::any_cast<double>( param_map[nnt::s_RHO_0] ) -
    boost::any_cast<double>( param_map[S_RHO_1] );
 
}

//##############################################################################
// initialize_state().
//##############################################################################

void
initialize_state(
  param_map_t& param_map,
  state_type& x
)
{

  x[0] = 1.;
  x[1] = 
     GSL_POW_4( x[0] ) *
     (
       (
         boost::any_cast<double>( param_map[S_RHO_1] ) /
         boost::any_cast<double>( nnt::s_TAU )
       )
       +
       (
         2. * boost::any_cast<double>( param_map[S_RHO_2] ) /
         boost::any_cast<double>( param_map[S_DELTA_TRAJ] )
       )
     )
     /
     ( 3. * boost::any_cast<double>( nnt::s_RHO_0 ) );

}

//##############################################################################
// acceleration().
//##############################################################################

double
acceleration(
  param_map_t& param_map,
  nnt::Zone & zone,
  const state_type& x,
  const double time
)
{
  x[2] =
     GSL_POW_3( x[0] ) *
     (
       (4. * x[1]*
         (
           (
             (
               boost::any_cast<double>( param_map[S_RHO_1] ) /
               boost::any_cast<double>( nnt::s_TAU )
             )
             *
             exp( -time / boost::any_cast<double>( param_map[nnt::s_TAU] )
           )
           +
           (
             2. * boost::any_cast<double>( param_map[S_RHO_2] ) /
             boost::any_cast<double>( param_map[S_DELTA_TRAJ] ) 
           )
         )
       )
       -
       ( x[0] *
         (
           ( boost::any_cast<double>( param_map[S_RHO_1] ) /
             GSL_POW_2( boost::any_cast<double>( nnt::s_TAU )
           )
           *
           (
             exp( -time / boost::any_cast<double>( param_map[nnt::s_TAU] )
           )
           +
           (
             6. * boost::any_cast<double>( param_map[S_RHO_2] ) /
             GSL_POW_2(boost::any_cast<double>( param_map[S_DELTA_TRAJ] )
           )
         )
       )
     )
     /
     ( 3. * boost::any_cast<double>( nnt::s_RHO_0 ) );
   
  return
    x[1] / ( 3. * boost::any_cast<double>( param_map[nnt::s_TAU] ) );

}

//##############################################################################
// rho_function().
//##############################################################################

double rho_function( param_map_t& param_map, const state_type& x )
{

  return
    boost::any_cast<double>( param_map[nnt::s_RHO_0] ) / gsl_pow_3( x[0] );

}

//##############################################################################
// t9_function().
//##############################################################################

double t9_function(
  nnt::Zone& zone,
  param_map_t& param_map,
  Libnucnet__NetView * p_view )
{

  double t9 =
    nnt::compute_1d_root(
      boost::bind(
        user::t9_from_entropy_root,
        _1,
        boost::ref( zone ),
        p_view
      ),
      zone.getProperty<double>( nnt::s_T9 ),
      boost::any_cast<double>( param_map[S_ROOT_FACTOR] )
    );

  return t9;

}

//##############################################################################
// observer_function().
//##############################################################################

void
observer_function(
  nnt::Zone& zone,
  const state_type& x,
  const state_type& dxdt,
  const double d_t
)
{

  double d_dt =
    d_t - zone.getProperty<double>( nnt::s_TIME );

  std::cout <<
    boost::format( "t = %.5e dt = %.5e\n" ) %
    d_t %
    d_dt;

  std::cout <<
    boost::format( "x = {%.5e, %.5e, %.5e}\n" ) %
    x[0] %
    x[1] %
    x[2]; 

  std::cout <<
    boost::format( "dxdt = {%.5e, %.5e, %.5e}\n\n" ) %
    dxdt[0] %
    dxdt[1] %
    dxdt[2];

}

}  // namespace my_user
