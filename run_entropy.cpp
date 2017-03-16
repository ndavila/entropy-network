////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Example code for running a network calculation with entropy
//!        generation.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>
#include <iostream>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>
#include <boost/assign.hpp>

#include <Libnucnet.h>

#include "nnt/two_d_weak_rates.h"
#include "user/remove_duplicate.h"
#include "user/user_rate_functions.h"
#include "user/network_limiter.h"
#include "user/flow_utilities.h"
#include "user/hydro_helper.h"

#include "my_hydro_helper.h"

typedef my_user::state_type my_state_type;

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_X_REG_T      0.15    /* x change regulator for dt update */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */
#define D_LIM_CUTOFF   1.e-25  /* Cutoff abundance for network limiter */

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "no"

//##############################################################################
// Strings.
//##############################################################################


#define S_X            "x"

#define S_ACCELERATION_FUNCTION  "acceleration function"
#define S_ENTROPY_FUNCTION  "entropy function"
#define S_ENTROPY_GENERATION_FUNCTION  "entropy generation function"
#define S_EVOLVE_FUNCTION "evolution function"
#define S_NUCNET       "nucnet"
#define S_NUC_XPATH       "nuc_xpath"
#define S_OBSERVE  "observe"
#define S_OBSERVER_FUNCTION  "observer function"
#define S_PARTICLE     nnt::s_TOTAL
#define S_PROGRAM_OPTIONS  "program_options"
#define S_REAC_XPATH       "reac_xpath"
#define S_RESPONSE_FILE    "response_file"
#define S_RHO_FUNCTION  "rho function"
#define S_SOLVER       nnt::s_ARROW // Solver type: ARROW or GSL
#define S_SDOT_NUC_XPATH  "sdot_nuc_xpath"
#define S_SDOT_REAC_XPATH  "sdot_reac_xpath"
#define S_T9_FUNCTION   "t9 function"
#define S_T9_GUESS   "t9_guess"
#define S_VIEW         "view"


#define B_OUTPUT_EVERY_TIME_DUMP    false  // Change to true to write to xml
                                           // every time dump.  False just
                                           // writes output at end of
                                           // calculation.
                                           //

namespace po = boost::program_options;

//##############################################################################
// at_option_parser().
//##############################################################################

std::pair<std::string,std::string> at_option_parser( std::string const& s )
{
  if ( '@' == s[0] )
    return std::make_pair( std::string( S_RESPONSE_FILE ), s.substr( 1 ) );
  else
    return std::pair<std::string,std::string>();
}

//##############################################################################
// response_file().
//##############################################################################

po::variables_map
response_file(
  po::variables_map& vm,
  po::options_description& all
)
{
  // Load the file and tokenize it
  std::ifstream ifs( vm[S_RESPONSE_FILE].as<std::string>().c_str() );
  if( !ifs )
  {
    std::cout << "Could not open the response file\n";
    exit( EXIT_FAILURE );
  }

  // Read the whole file into a string
  std::stringstream ss;
  ss << ifs.rdbuf();

  // Split the file content
  std::string sep1("");
  std::string sep2(" \n\r");
  std::string sep3("\"");
  std::string sstr = ss.str();

  boost::escaped_list_separator<char> els( sep1, sep2, sep3 );
  boost::tokenizer< boost::escaped_list_separator<char> > tok( sstr, els );

  std::vector<std::string> args;
  std::copy( tok.begin(), tok.end(), std::back_inserter( args ) );

  // Parse the file and store the options
  store( po::command_line_parser( args ).options( all ).run(), vm );

  return vm;
}

//##############################################################################
// entropy_generation_rhs.
//##############################################################################

class entropy_generation_rhs
{
  nnt::Zone& zone;
  Libnucnet__NetView *pView;

  public:
    entropy_generation_rhs(
      nnt::Zone& _zone,
      Libnucnet__NetView * p_view
    ) : zone( _zone ), pView( p_view ) {}

    void operator()(
      const my_state_type &x, my_state_type &dxdt, const double d_t
    )
    {

      double d_dt, d_entropy_generation, d_t9_old;
      gsl_vector * p_abundances, * p_abundance_changes;
      boost::function<double( const my_state_type& )> rho_func;

      p_abundances = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

      p_abundance_changes =
        Libnucnet__Zone__getAbundanceChanges( zone.getNucnetZone() );
  
      d_t9_old = zone.getProperty<double>( nnt::s_T9 );

      d_dt =
        d_t 
        -
        zone.getProperty<double>( nnt::s_TIME );

      zone.updateProperty(
        nnt::s_DTIME,
        d_dt
      );

      zone.updateProperty( nnt::s_ENTROPY_PER_NUCLEON, x[2] );

      zone.updateProperty(
        nnt::s_RHO,
        boost::any_cast<boost::function<double( const my_state_type& )>
        >(
          zone.getFunction( S_RHO_FUNCTION )
        )( x )
      );

      zone.updateProperty(
        nnt::s_T9,
        boost::any_cast<
          boost::function<double( Libnucnet__NetView * )>
        >(
          zone.getFunction( S_T9_FUNCTION )
        )( pView )
      );

      boost::any_cast<
        boost::function<void( Libnucnet__NetView *, const double )>
      >(
        zone.getFunction( S_EVOLVE_FUNCTION )
      )( pView, d_dt );

      dxdt[0] = x[1];

      dxdt[1] = 
        boost::any_cast<
          boost::function<double( const my_state_type&, const double )>
        >(
          zone.getFunction( S_ACCELERATION_FUNCTION )
        )( x, d_t );

      d_entropy_generation =
        boost::any_cast<boost::function<double( Libnucnet__NetView * )> >(
          zone.getFunction( S_ENTROPY_GENERATION_FUNCTION )
        )( pView );

      dxdt[2] = d_entropy_generation; // - d_energy_loss;

      if( zone.hasFunction( S_OBSERVER_FUNCTION ) )
      {
        boost::any_cast<
          boost::function<
            void(
              const my_state_type&,
              const my_state_type&,
              const double
            )
          >
        >( zone.getFunction( S_OBSERVER_FUNCTION ) )( x, dxdt, d_t );
      }

      Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );
      
      Libnucnet__Zone__updateAbundanceChanges(
        zone.getNucnetZone(),
        p_abundance_changes
      );

      zone.updateProperty( nnt::s_T9, d_t9_old );

  }
      
}; 

//##############################################################################
// program_options().
//##############################################################################

void
program_options(
  po::variables_map& vmap,
  po::options_description& help,
  po::options_description& general,
  po::options_description& network,
  po::options_description& user,
  po::options_description& all
)
{

  const std::string& s = vmap[S_PROGRAM_OPTIONS].as<std::string>();

  if( s == "help" )
  {
    std::cout << help << std::endl;
  }
  else if( s == "general" )
  {
    std::cout << general << std::endl;
  }
  else if( s == "network" )
  {
    std::cout << network << std::endl;
  }
  else if( s == "user" )
  {
    std::cout << user << std::endl;
  }
  else if( s == "all" )
  {
    std::cout << all << std::endl;
  }
  else
  {
    std::cout << "\nUnknown options_description '" << s << "' in the "
    "--program_options option\n\n";

    exit( EXIT_FAILURE );
  }

  exit( EXIT_SUCCESS );

}

//##############################################################################
// get_input().
//##############################################################################

my_user::param_map_t
get_input( int argc, char **argv )
{

  try
  {

    my_user::param_map_t param_map;
    po::variables_map vm;        
    std::string s_nuc_xpath = "", s_reac_xpath = "";
    std::string s_sdot_nuc_xpath = "", s_sdot_reac_xpath = "";

    std::string s_purpose = "\nPurpose: run a network calculation with entropy generation for the input xml_file for the selected nuclei and reactions and for the selected nuclei and reactions for entropy generation.";

    po::options_description help( "\nHelp Options" );
    help.add_options()

      ( "help", "print out usage statement and exit\n" )

      ( "example", "print out example usage and exit\n" )

      ( "program_options", po::value<std::string>(),
        "print out list of program options (help, general, network,"
        " user, or all)"
        " and exit"
      )
    ;

    po::options_description general("\nGeneral options");
    general.add_options()
      (
       nnt::s_TIME,
       po::value<double>()->default_value( 0., "0." ),
       "Initial time (s)"
      )
      (
       nnt::s_DTIME,
       po::value<double>()->default_value( 1.e-15, "1.e-15" ),
       "Initial time step (s)"
      )
      (
       nnt::s_TEND,
       po::value<double>()->default_value( 10., "10." ),
       "End time (s)"
      )
      (
       nnt::s_STEPS,
       po::value<size_t>()->default_value( 20 ),
       "Frequency of time step dump"
      )
      (
       nnt::s_MU_NUE_KT,
       po::value<std::string>()->default_value( "-inf" ),
       "Electron neutrino chemical potential / kT"
      )
      (
       S_T9_GUESS,
       po::value<std::string>()->default_value( "yes" ),
       "Guess next T9"
      )
      (
       S_OBSERVE,
       po::value<std::string>()->default_value( "no" ),
       "Observe steps"
      )

      ( S_RESPONSE_FILE, po::value<std::string>(),
        "can be specified with '@name', too\n"
      )

    ;

    po::options_description network("\nNetwork options");
    network.add_options()
      (
       S_NUC_XPATH,
       po::value<std::vector<std::string> >()->multitoken()->composing(),
       "XPath to select nuclei (default: all nuclides)\n"
      )
      (
       S_REAC_XPATH,
       po::value<std::vector<std::string> >()->multitoken()->composing(),
       "XPath to select reactions (default: all reactions)\n"
      )
      (
       S_SDOT_NUC_XPATH,
       po::value<std::vector<std::string> >()->multitoken()->composing(),
       "XPath to select nuclides for entropy generation (default: a step's evolution network nuclides)"
      )
      (
       S_SDOT_REAC_XPATH,
       po::value<std::vector<std::string> >()->multitoken()->composing(),
       "XPath to select reactions for entropy generation (default: a step's evolution network reactions)"
      )
      (
       nnt::s_USE_SCREENING,
       po::value<std::string>()->default_value( "no" ),
       "Use screening"
      )
      (
       nnt::s_USE_NSE_CORRECTION,
       po::value<std::string>()->default_value( "no" ),
       "Use NSE correction"
      )

    ;
 
    // Read in user-defined options
    po::options_description user("\nUser-defined options");
    my_user::get_user_defined_descriptions( user );

    po::options_description all( "\nAll Allowed Options" );
    all.add( help ).add( general ).add( network ).add( user );

    store(
      po::command_line_parser( argc, argv ).
      options( all ).
      extra_parser( at_option_parser ).
      run(),
      vm
    );

    if( vm.count( "example" ) )
    {
      std::cerr <<
        "\n" << argv[0] <<
        " ../../data_pub/my_net.xml ../../data/my_zone.xml output.xml " <<
        " --nuc_xpath \"[z <= 30]\"\n " << std::endl;
      exit( EXIT_SUCCESS );
    }

    if( argc == 1 || vm.count("help") == 1 )
    {
      std::cerr <<
        "\nUsage: " << argv[0] << " net_xml zone_xml output_xml [options]" <<
        std::endl;
      std::cerr << s_purpose << std::endl;
      std::cout << help << "\n";
      exit( EXIT_FAILURE );
    }

    if( vm.count( S_PROGRAM_OPTIONS ) )
    {
      program_options(
        vm,
        help,
        general,
        network,
        user,
        all
      );
    }

    if( vm.count( S_RESPONSE_FILE ) )
      vm = response_file( vm, general );

    //==========================================================================
    // XPath strings.
    //==========================================================================

    if( vm.count( S_NUC_XPATH ) )
    {
      BOOST_FOREACH(
        std::string s,
        vm[S_NUC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_nuc_xpath += s + " ";
      }
    }

    if( vm.count( S_REAC_XPATH ) )
    {
      BOOST_FOREACH(
        std::string s,
        vm[S_REAC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_reac_xpath += s + " ";
      }
    }

    if( vm.count( S_SDOT_NUC_XPATH ) )
    {
      BOOST_FOREACH(
        std::string s,
        vm[S_SDOT_NUC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_sdot_nuc_xpath += s + " ";
      }
    }

    if( vm.count( S_SDOT_REAC_XPATH ) )
    {
      BOOST_FOREACH(
        std::string s,
        vm[S_SDOT_REAC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_sdot_reac_xpath += s + " ";
      }
    }

    //==========================================================================
    // Validate input file.
    //==========================================================================

    if( strcmp( VALIDATE, "yes" ) == 0 )
    {
      if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
        fprintf( stderr, "Not valid libnucnet input!\n" );
        exit( EXIT_FAILURE );
      }
    }

    //==========================================================================
    // Get network and view.
    //==========================================================================

    param_map[S_NUCNET] = Libnucnet__new();

    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( boost::any_cast<Libnucnet *>( param_map[S_NUCNET] ) ),
      argv[1],
      s_nuc_xpath.c_str(),
      s_reac_xpath.c_str()
    );

    Libnucnet__assignZoneDataFromXml(
      boost::any_cast<Libnucnet *>( param_map[S_NUCNET] ),
      argv[2],
      ""
    );

    if( vm.count(S_SDOT_NUC_XPATH) == 1 || vm.count(S_SDOT_REAC_XPATH) == 1 )
    {
      param_map[S_VIEW] = 
        Libnucnet__NetView__new(
          Libnucnet__getNet(
            boost::any_cast<Libnucnet *>( param_map[S_NUCNET] )
          ),
          s_sdot_nuc_xpath.c_str(),
          s_sdot_reac_xpath.c_str()
        );
    }

    //==========================================================================
    // Get other data.
    //==========================================================================

    param_map[nnt::s_TIME] = vm[nnt::s_TIME].as<double>();
    param_map[nnt::s_DTIME] = vm[nnt::s_DTIME].as<double>();
    param_map[nnt::s_TEND] = vm[nnt::s_TEND].as<double>();
    param_map[nnt::s_STEPS] = vm[nnt::s_STEPS].as<size_t>();
    param_map[nnt::s_USE_SCREENING] = vm[nnt::s_USE_SCREENING].as<std::string>();
    param_map[nnt::s_USE_NSE_CORRECTION] =
      vm[nnt::s_USE_NSE_CORRECTION].as<std::string>();
    param_map[S_T9_GUESS] = vm[S_T9_GUESS].as<std::string>();
    param_map[S_OBSERVE] = vm[S_OBSERVE].as<std::string>();
    param_map[nnt::s_MU_NUE_KT] = vm[nnt::s_MU_NUE_KT].as<std::string>();

    // Set user-defined options
    my_user::set_user_defined_options( vm, param_map );

    return param_map;

  }
  catch( std::exception& e )
  {
    std::cerr << "error: " << e.what() << "\n";
    exit( EXIT_FAILURE );
  }
  catch(...)
  {
    std::cerr << "Exception of unknown type!\n";
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int k = 0;
  size_t i_step = 0;
  double d_t, d_dt, d_t9_old, d_dt9dt;
  my_user::param_map_t param_map;
  Libnucnet * p_my_nucnet, * p_my_output;
  Libnucnet__NetView * p_view = NULL;
  nnt::Zone zone;
  char s_property[32];
  std::set<std::string> isolated_species_set;

  my_state_type
    x(3), xold(3), x_lim = boost::assign::list_of(1.e-10)(1.)(1.e-5);

  //============================================================================
  // Check input.
  //============================================================================

  param_map = get_input( argc, argv );

  p_my_nucnet = boost::any_cast<Libnucnet *>( param_map[S_NUCNET] );

  if( param_map.find( S_VIEW ) != param_map.end() )
  {
    p_view = boost::any_cast<Libnucnet__NetView *>( param_map[S_VIEW] );
  }

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Set the zone.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) &&
    zone.getProperty<std::string>( nnt::s_USE_APPROXIMATE_WEAK_RATES ) == "yes"
  )
  {
      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
  }

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Set screening, Coulomb correction, nue/kT, and rate update functions.
  //============================================================================

  if( boost::any_cast<std::string>( param_map[nnt::s_USE_SCREENING] ) == "yes" )
  {
    user::set_screening_function( zone );
  }

  if(
    boost::any_cast<std::string>(
      param_map[nnt::s_USE_NSE_CORRECTION] ) == "yes"
    )
  {
    user::set_nse_correction_function( zone );
  }

  user::set_rate_data_update_function( zone );

  zone.updateProperty(
    nnt::s_MU_NUE_KT,
    boost::any_cast<std::string>( param_map[nnt::s_MU_NUE_KT] )
  );

  //============================================================================
  // Remove isolated species if desired.
  //============================================================================

  if(
    Libnucnet__Reac__getNumberOfReactions(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
    ) != 0
  )
  {

    isolated_species_set =
      user::get_isolated_species(
        Libnucnet__getNet( p_my_nucnet ),
        "",
        ""
      );

    BOOST_FOREACH( std::string s_species, isolated_species_set )
    {

//    Careful that you don't remove a species with non-zero abundance!

      std::cout << s_species << std::endl;

      Libnucnet__Nuc__removeSpecies(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__getNet( p_my_nucnet )
          ),
          s_species.c_str()
        )
      );

    }

  }
  
  //============================================================================
  // Set the acceleration function.
  // Must return the scaled acceleration for the given available data.
  //============================================================================

  zone.updateFunction(
    S_ACCELERATION_FUNCTION,
    static_cast<boost::function<double( const my_state_type&, const double )> >(
      boost::bind(
        my_user::acceleration,
        boost::ref( param_map ),
        boost::ref( zone ),
        _1,
        _2
      )
    )
  );

  //============================================================================
  // Set the rho function.
  // Must return the density for the given available data.
  //============================================================================

  zone.updateFunction(
    S_RHO_FUNCTION,
    static_cast<boost::function<double( const my_state_type& )> >(
      boost::bind(
        my_user::rho_function,
        boost::ref( param_map ),
        _1
      )
    )
  );

  //============================================================================
  // Set the t9 function.
  // Must return the t9 for the given available data.
  //============================================================================

  zone.updateFunction(
    S_T9_FUNCTION,
    static_cast<boost::function<double( Libnucnet__NetView * )> >(
      boost::bind(
        my_user::t9_function,
        boost::ref( zone ),
        boost::ref( param_map ),
        _1
      )
    )
  );

  //============================================================================
  // Set the entropy function.
  //============================================================================

  zone.updateFunction(
    S_ENTROPY_FUNCTION,
    static_cast<boost::function<double( )> >(
      boost::bind(
        user::compute_entropy,
        boost::ref( zone )
      )
    )
  );

  //============================================================================
  // Set the entropy generation function.
  //============================================================================

  zone.updateFunction(
    S_ENTROPY_GENERATION_FUNCTION,
    static_cast<boost::function<double( Libnucnet__NetView * )> >(
      boost::bind(
        user::compute_entropy_generation_rate,
        boost::ref( zone ),
        _1
      )
    )
  );

  //============================================================================
  // Set the abundance evolver with basic prototype
  // Must evolve the abundances over a given input timestep.
  //============================================================================

  zone.updateFunction(
    S_EVOLVE_FUNCTION,
    static_cast<boost::function<void( Libnucnet__NetView *, const double )> >(
      boost::bind(
        user::evolve_function,
        boost::ref( zone ),
        _1,
        _2
      )
    )
  );

  //============================================================================
  // Set the observer function with basic prototype
  //   void( const my_state_type& x, const my_state_type& dxdt, const double t )
  //============================================================================

  if( boost::any_cast<std::string>( param_map[S_OBSERVE] ) == "yes" )
  {
    zone.updateFunction(
      S_OBSERVER_FUNCTION,
      static_cast<
        boost::function<
          void(
            const my_state_type&,
            const my_state_type&,
            const double
          )
        >
      >( boost::bind(
           my_user::observer_function,
           boost::ref( zone ),
           _1,
           _2,
           _3
         )
      )
    );
  }
            
  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( S_SOLVER == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "3" );

  }

  //============================================================================
  // Create output.
  //============================================================================

  p_my_output = nnt::create_network_copy( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  //============================================================================
  // Initialize the system.
  //============================================================================

  zone.updateProperty(
    nnt::s_T9,
    boost::any_cast<double>( param_map[nnt::s_T9_0] )
  );

  zone.updateProperty(
    nnt::s_RHO,
    boost::any_cast<double>( param_map[nnt::s_RHO_0] )
  );

  d_t9_old = zone.getProperty<double>( nnt::s_T9 );

  d_dt9dt = 0;
  
  zone.updateProperty( nnt::s_PARTICLE, S_PARTICLE);

  d_dt = boost::any_cast<double>( param_map[nnt::s_DTIME] );

  d_t = boost::any_cast<double>( param_map[nnt::s_TIME] );

  my_user::initialize_state( param_map, x );

  x[2] =
    boost::any_cast< boost::function<double( )> >(
      zone.getFunction( S_ENTROPY_FUNCTION )
    )( );

  user::limit_evolution_network( zone, D_LIM_CUTOFF );

  //============================================================================
  // Choose the stepper.
  //============================================================================

  boost::numeric::odeint::adams_bashforth<4, my_state_type > stepper;

  //============================================================================
  // Evolve network while t < final t. 
  //============================================================================

  while ( d_t < boost::any_cast<double>( param_map[nnt::s_TEND] ) )
  {

  //============================================================================
  // Set time.
  //============================================================================

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );  

  //============================================================================
  // Save old values.
  //============================================================================

    std::copy( x.begin(), x.end(), xold.begin() );

  //============================================================================
  // Evolve step.
  //============================================================================

    Libnucnet__NetView * p_sdot_view;

    if( p_view )
      p_sdot_view = p_view;
    else
      p_sdot_view = zone.getNetView( EVOLUTION_NETWORK );

    entropy_generation_rhs my_rhs( zone, p_sdot_view );

    stepper.do_step( my_rhs, x, d_t, d_dt );

  //============================================================================
  // Update properties.
  //============================================================================

    d_t += d_dt;

    zone.updateProperty( nnt::s_DTIME, d_dt );

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

    zone.updateProperty(
      nnt::s_RHO,
      my_user::rho_function( param_map, x )
    );

    zone.updateProperty(
      nnt::s_ENTROPY_PER_NUCLEON,
      x[2]
    );

    if( boost::any_cast<std::string>( param_map[S_T9_GUESS] ) == "yes" )
    {
      zone.updateProperty(
        nnt::s_T9,
        d_t9_old + d_dt9dt * d_dt
      );
    }

std::cout << "Now update: " << std::endl;
    zone.updateProperty(
      nnt::s_T9,
      my_user::t9_function(
        zone,
        param_map,
        zone.getNetView( EVOLUTION_NETWORK )
      )
    );

    if( boost::any_cast<std::string>( param_map[S_T9_GUESS] ) == "yes" )
    {
      d_dt9dt = ( zone.getProperty<double>( nnt::s_T9 ) - d_t9_old ) / d_dt;
      d_t9_old = zone.getProperty<double>( nnt::s_T9 );
    }

    boost::any_cast<
      boost::function<void( Libnucnet__NetView *, const double )>
    >(
      zone.getFunction( S_EVOLVE_FUNCTION )
    )( zone.getNetView( EVOLUTION_NETWORK ), d_dt );

    zone.updateProperty( S_X, "0", x[0] );

    zone.updateProperty( S_X, "1", x[1] );

  //============================================================================
  // Output step data.
  //============================================================================

    if( boost::any_cast<std::string>( param_map[S_OBSERVE] ) == "yes" )
    {
      std::cout <<
        boost::format( "t = %g, x = {%g, %g, %g}\n\n" ) %
        d_t % x[0] % x[1] % x[2];
      std::cout << boost::format( "-----------\n\n" );
    }

  //============================================================================
  // Print out abundances.
  //============================================================================

    if( i_step++ % boost::any_cast<size_t>( param_map[nnt::s_STEPS] ) == 0 ||
        d_t >= boost::any_cast<double>( param_map[nnt::s_TEND] )
    )
    {
      sprintf( s_property, "%d", ++k );
      Libnucnet__relabelZone(
        p_my_nucnet,
        zone.getNucnetZone(),
        s_property,
        NULL,
        NULL
      );
      nnt::print_zone_abundances( zone );
      nnt::write_xml( p_my_output, zone.getNucnetZone() );
      if( B_OUTPUT_EVERY_TIME_DUMP )
      {
        Libnucnet__writeToXmlFile( p_my_output, argv[3] );
      }
    }

  //============================================================================
  // Limit network.
  //============================================================================

  user::limit_evolution_network( zone, D_LIM_CUTOFF );

  //============================================================================
  // Update timestep.
  //============================================================================

    double d_h = 1.e99;
    for( size_t i = 0; i < x.size(); i++ )
    {
      double delta = fabs( ( x[i] - xold[i] ) / x[i] );
      if( delta > 0 && fabs( x[i] ) > x_lim[i] )
        d_h = GSL_MIN( d_h, D_X_REG_T * d_dt / delta );
    }

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( d_dt > d_h ) d_dt = d_h;

    if ( d_t + d_dt > boost::any_cast<double>( param_map[nnt::s_TEND] ) )
    {
      d_dt = boost::any_cast<double>( param_map[nnt::s_TEND] ) - d_t;
    }

  }  

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__writeToXmlFile( p_my_output, argv[3] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
