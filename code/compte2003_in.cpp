/*
 *  compte2003_in.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "compte2003_in.h"
#include "nest_names.h"

#ifdef HAVE_GSL_1_11

#include "universal_data_logger_impl.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdio>

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::compte2003_in > nest::compte2003_in::recordablesMap_;

namespace nest
{
/*
 * template specialization must be placed in namespace
 *
 * Override the create() method with one call to RecordablesMap::insert_()
 * for each quantity to be recorded.
 */
template <>
void
RecordablesMap< compte2003_in >::create()
{
 
  // soma
  insert_( names::V_m_s, &compte2003_in::get_y_elem_< compte2003_in::State_::V_M_S > );
  insert_( names::Na_h, &compte2003_in::get_y_elem_< compte2003_in::State_::NA_H > );
  insert_( names::K_m, &compte2003_in::get_y_elem_< compte2003_in::State_::K_M > );
  insert_( names::g_ampa, &compte2003_in::get_y_elem_< compte2003_in::State_::G_AMPA > );
  insert_( names::g_nmda_fast, &compte2003_in::get_y_elem_< compte2003_in::State_::G_NMDA_FAST > );
  insert_( names::g_nmda_slow, &compte2003_in::get_y_elem_< compte2003_in::State_::G_NMDA_SLOW > );
  insert_( names::g_gaba, &compte2003_in::get_y_elem_< compte2003_in::State_::G_GABA > );
}
}


extern "C" int
nest::compte2003_in_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef nest::compte2003_in::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::compte2003_in& node = *( reinterpret_cast< nest::compte2003_in* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].
 
  
  // membrane potential
  const double_t& V = y[ S::V_M_S ];
  
  //-----inhibitory synaptic current-----
  double g_syn_inh = y[S::G_GABA]; 
  
  double I_inh = g_syn_inh*( node.P_.E_in-V);

  //-----excitatory synaptic currents-----
  double g_syn_exc = y[S::G_AMPA]+(y[S::G_NMDA_SLOW]-y[S::G_NMDA_FAST]);

  double I_exc = g_syn_exc*( node.P_.E_ex - V);
    
  //-----fast Na+ current-----
  double Na_m_inf; 
  double Na_h_inf; 
  double Na_h_tau;
  double a,b;

  //activation
  a=0.5*(V+35.)/(1.-exp(-(V+35.)/10.));
  b=20.*exp(-(V+60.)/18.);
  Na_m_inf=a/(a+b); 
     
  //inactivation   
  a=0.35*exp(-(V+58.)/20.);
  b=5./(1.+exp(-(V+28.)/10.));
  Na_h_inf=a/(a+b); 
  Na_h_tau=1./(a+b);    //[ms]
     
  double I_Na = node.P_.g_Na*pow(Na_m_inf,3)*y[S::NA_H]*(node.P_.E_Na-V); //[pA]
     
  //-----fast K+ current
  double K_m_inf; 
  double K_m_tau;

  a=0.05*(V+34.)/(1-exp(-(V+34.)/10.));
  b=0.625*exp(-(V+44.)/80.);
  K_m_inf=a/(a+b); 
  K_m_tau=1./(a+b);     //[ms]
     
  double I_K = node.P_.g_K*pow(y[S::K_M],4)*(node.P_.E_K-V);  //[pA]

  //-----leak current
  double I_L = node.P_.g_L*(node.P_.E_L-V);


  //-----derivatives
  // membrane potential
  f[S::V_M_S] = ( I_L + I_Na + I_K + I_inh + I_exc + node.B_.I_stim_ + node.P_.I_e ) / node.P_.C_m_s;     
     
  // excitatory conductance
  f[S::G_AMPA] = -y[S::G_AMPA] / node.P_.tau_syn_ampa; // Synaptic Conductance (nS)
  f[S::G_NMDA_FAST] = -y[S::G_NMDA_FAST] / node.P_.tau_syn_nmda_fast; // Synaptic Conductance (nS)
  f[S::G_NMDA_SLOW] = -y[S::G_NMDA_SLOW] / node.P_.tau_syn_nmda_slow; // Synaptic Conductance (nS)
    
  // inhibitory conductance
  f[S::G_GABA] = -y[S::G_GABA] / node.P_.tau_syn_gaba; // Synaptic Conductance (nS)
     
  // fast Na+ inactivation
  f[S::NA_H] = (Na_h_inf - y[S::NA_H]) / Na_h_tau;  
  // fast K+ activation
  f[S::K_M]  = (K_m_inf  - y[S::K_M])  / K_m_tau;               
       
  return GSL_SUCCESS;
   
   
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::compte2003_in::Parameters_::Parameters_()
  : t_ref   (2.0 ), // ms
    E_ex    (0.0),  // mV
    E_in    (-85.0),  // mV
    E_L     (-70.0),  // mV
    E_Na     (50.0),  // mV
    E_K     (-90.0),  // mV
    
    tau_syn_ampa  (3.),  // ms
    tau_syn_nmda_fast  (5.), // ms
    tau_syn_nmda_slow  (50.), // ms
    tau_syn_gaba        (3.),  // ms

    g_Na    (0.0),  // nS
    g_K     (0.0),  // nS
    g_L     (10.0),  // nS
    C_m_s   (150.0),  // pF        
    I_e     (0.0)  // pA


{
}


nest::compte2003_in::State_::State_(const Parameters_ & p)
  : r_(0)
{
  // for simplicity, we first initialize all values to 0,
  // then set the membrane potentials for each compartment
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = 0;

  double V=-80.;
   
  //Fast Na+ current inactivation
  double a=0.35*exp(-(V+58.)/20.);
  double b=5./(1.+exp(-(V+28.)/10.));
  double Na_h_inf=a/(a+b); 

  //Fast K+ current
  a=0.05*(V+34.)/(1-exp(-(V+34.)/10.));
  b=0.625*exp(-(V+44.)/80.);
  double K_m_inf=a/(a+b); 
      
  y_[0] = V;    //V_m_s    
  y_[5] = Na_h_inf; //Na_h
  y_[6]  = K_m_inf; //K_m
 
}


nest::compte2003_in::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
}

nest::compte2003_in::State_& nest::compte2003_in::State_::operator=( const State_& s )
{
  assert( this != &s ); // would be bad logical error in program

  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
  r_ = s.r_;
  return *this;
}

/* ----------------------------------------------------------------
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::compte2003_in::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::t_ref,      t_ref);
  def<double>(d,names::E_ex,       E_ex);
  def<double>(d,names::E_in,       E_in);
  def<double>(d,names::E_L,        E_L); 
  def<double>(d,names::E_Na,       E_Na);
  def<double>(d,names::E_K ,       E_K );


  def<double>(d,names::tau_syn_ampa, tau_syn_ampa);
  def<double>(d,names::tau_syn_nmda_fast, tau_syn_nmda_fast);
  def<double>(d,names::tau_syn_nmda_slow, tau_syn_nmda_slow);
  def<double>(d,names::tau_syn_gaba, tau_syn_gaba);
    
  def<double>(d, names::g_L,        g_L);
  def<double>(d, names::g_Na,       g_Na);
  def<double>(d, names::g_K,        g_K);
  def<double>(d, names::C_m_s,      C_m_s);
  def<double>(d, names::I_e,        I_e);


  
}



void nest::compte2003_in::Parameters_::set(const DictionaryDatum &d)
{
  // allow setting the membrane potential
  updateValue<double>(d,names::t_ref,     t_ref);
  updateValue<double>(d,names::E_ex,       E_ex);
  updateValue<double>(d,names::E_in,       E_in);
  updateValue<double>(d,names::E_L,        E_L); 
  updateValue<double>(d,names::E_Na,       E_Na);
  updateValue<double>(d,names::E_K ,       E_K );


  updateValue<double>(d,names::tau_syn_ampa, tau_syn_ampa);
  updateValue<double>(d,names::tau_syn_nmda_fast, tau_syn_nmda_fast);
  updateValue<double>(d,names::tau_syn_nmda_slow, tau_syn_nmda_slow);
  updateValue<double>(d,names::tau_syn_gaba, tau_syn_gaba);

  updateValue<double>(d, names::g_L,        g_L);
  updateValue<double>(d, names::g_Na,       g_Na);
  updateValue<double>(d, names::g_K,        g_K);
  updateValue<double>(d, names::C_m_s,      C_m_s);
  updateValue<double>(d, names::I_e,        I_e);


  // apply checks 
  if ( C_m_s <= 0 )
      throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref < 0 )
      throw BadProperty("Refractory time cannot be negative.");
          
    if ( tau_syn_ampa <= 0 || tau_syn_nmda_fast <= 0 || tau_syn_nmda_slow <= 0 || tau_syn_gaba <= 0 )
      throw BadProperty("All time constants must be strictly positive.");
  
}

void nest::compte2003_in::State_::get(DictionaryDatum &d) const
{

  // we assume here that State_::get() always is called after Parameters_::get(),
  // so that the per-compartment dictionaries exist
  def<double>(d, names::V_m_s, y_[V_M_S]); // Membrane potential
  def<double>(d, names::Na_h, y_[NA_H]); // 
  def<double>(d, names::K_m , y_[K_M]); // 
  def<double>(d, names::g_ampa , y_[G_AMPA]); // 
  def<double>(d, names::g_nmda_slow , y_[G_NMDA_SLOW]); // 
  def<double>(d, names::g_nmda_fast , y_[G_NMDA_FAST]); // 
  def<double>(d, names::g_gaba , y_[G_GABA]); // 
 
}

void nest::compte2003_in::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  // extract from sub-dictionaries
  updateValue<double>(d, names::V_m_s, y_[V_M_S]);
  updateValue<double>(d, names::Na_h, y_[NA_H]);
  updateValue<double>(d, names::K_m , y_[K_M]);
  updateValue<double>(d, names::g_ampa , y_[G_AMPA]);
  updateValue<double>(d, names::g_nmda_fast , y_[G_NMDA_FAST]);
  updateValue<double>(d, names::g_nmda_slow , y_[G_NMDA_SLOW]);
  updateValue<double>(d, names::g_gaba , y_[G_GABA]);

  if ( y_[ G_AMPA ] < 0 || y_[ G_NMDA_FAST ] < 0 || y_[ G_NMDA_SLOW ] < 0 || y_[ G_GABA ] < 0  )
    throw BadProperty( "Conductances must not be negative." );

}


nest::compte2003_in::Buffers_::Buffers_( compte2003_in& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::compte2003_in::Buffers_::Buffers_( const Buffers_&, compte2003_in& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::compte2003_in::compte2003_in()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

nest::compte2003_in::compte2003_in( const compte2003_in& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

nest::compte2003_in::~compte2003_in()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
    gsl_odeiv2_step_free( B_.s_ );
  if ( B_.c_ )
    gsl_odeiv2_control_free( B_.c_ );
  if ( B_.e_ )
    gsl_odeiv2_evolve_free( B_.e_ );
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::compte2003_in::init_state_( const Node& proto )
{
  const compte2003_in& pr = downcast< compte2003_in >( proto );
  S_ = pr.S_;
}

void
nest::compte2003_in::init_buffers_()
{
  B_.spikes_.resize( NUM_SPIKE_RECEPTORS );
  for ( size_t n = 0; n < NUM_SPIKE_RECEPTORS; ++n )
    B_.spikes_[ n ].clear(); // includes resize

  B_.currents_.resize( NUM_CURR_RECEPTORS );
  for ( size_t n = 0; n < NUM_CURR_RECEPTORS; ++n )
    B_.currents_[ n ].clear(); // includes resize

  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv2_step_alloc( gsl_odeiv2_step_rkf45, State_::STATE_VEC_SIZE );
  else
    gsl_odeiv2_step_reset( B_.s_ );

  if ( B_.c_ == 0 )
    B_.c_ = gsl_odeiv2_control_yp_new( 1e-6, 1e-6 );
  else
    gsl_odeiv2_control_init( B_.c_, 1e-6, 1e-6, 0.0, 1.0 );

  if ( B_.e_ == 0 )
    B_.e_ = gsl_odeiv2_evolve_alloc( State_::STATE_VEC_SIZE );
  else
    gsl_odeiv2_evolve_reset( B_.e_ );

  B_.sys_.function = compte2003_in_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::compte2003_in::calibrate()
{
  B_.logger_.init(); // ensures initialization in case mm connected after Simulate
  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref ) ).get_steps();
  assert( V_.RefractoryCounts_ >= 0 ); // since t_ref >= 0, this can only fail in error
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::compte2003_in::update( const Time& origin, const long_t from, const long_t to )
{
  assert( to >= 0 && ( delay ) from < Scheduler::get_min_delay() );
  assert( from < to );

  for ( long_t lag = from; lag < to; ++lag )
  {
    double t = 0.0;
    const double V_old = S_.y_[State_::V_M_S]; 
    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv2_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t
    while ( t < B_.step_ )
    {
        
      const int status = gsl_odeiv2_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &t,                   // from t
        B_.step_,             // to t <= step
        &B_.IntegrationStep_, // integration step size
        S_.y_ );              // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure( get_name(), status );
    }
    
    
    // refractoriness and spiking
    if ( S_.r_ )
      { // neuron is absolute refractory and cannot spike
        --S_.r_;
      }
    else if  ( S_.y_[State_::V_M_S] >= 0 && V_old > S_.y_[State_::V_M_S])
        //spike is emitted if Vm > 0 and starts to decay

	  {
	    S_.r_ = V_.RefractoryCounts_;
	  
	    set_spiketime(Time::step(origin.get_steps()+lag+1));
	  
	    SpikeEvent se;
	    network()->send(*this, se, lag);
	  }
    
      


    //---manage synaptic input
    double  G_syn;
    //excitatory inputs
    S_.y_[State_::G_AMPA] += B_.spikes_[0].get_value(lag);
	G_syn = B_.spikes_[1].get_value(lag);
    S_.y_[State_::G_NMDA_FAST] += G_syn;
    S_.y_[State_::G_NMDA_SLOW] += G_syn;
      
    //inhibitory inputs
    S_.y_[State_::G_GABA] += B_.spikes_[3].get_value(lag);

    // set new input current
    B_.I_stim_ = B_.currents_[0].get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}


void
nest::compte2003_in::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );
  assert( 0 <= e.get_rport() && e.get_rport() < 5 );

  B_.spikes_[ e.get_rport() ].add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
}

void
nest::compte2003_in::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  // add weighted current; HEP 2002-10-04
  B_.currents_[ e.get_rport() ].add_value(
    e.get_rel_delivery_steps( network()->get_slice_origin() ), e.get_weight() * e.get_current() );
}



void
nest::compte2003_in::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL_1_11
