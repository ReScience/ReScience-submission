/*
 *  compte2003_ex.cpp
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

#include "compte2003_ex.h"
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

nest::RecordablesMap< nest::compte2003_ex > nest::compte2003_ex::recordablesMap_;

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
RecordablesMap< compte2003_ex >::create()
{
 
  // soma
  insert_( names::V_m_s, &compte2003_ex::get_y_elem_< compte2003_ex::State_::V_M_S > );
  insert_( names::n_Na, &compte2003_ex::get_y_elem_< compte2003_ex::State_::N_NA > );
  insert_( names::Na_h, &compte2003_ex::get_y_elem_< compte2003_ex::State_::NA_H > );
  insert_( names::K_m, &compte2003_ex::get_y_elem_< compte2003_ex::State_::K_M > );
  insert_( names::K_A_h, &compte2003_ex::get_y_elem_< compte2003_ex::State_::K_A_H > );
  insert_( names::K_s_m, &compte2003_ex::get_y_elem_< compte2003_ex::State_::K_S_M > );
  insert_( names::g_ampa, &compte2003_ex::get_y_elem_< compte2003_ex::State_::G_AMPA > );
  insert_( names::g_nmda_fast, &compte2003_ex::get_y_elem_< compte2003_ex::State_::G_NMDA_FAST > );
  insert_( names::g_nmda_slow, &compte2003_ex::get_y_elem_< compte2003_ex::State_::G_NMDA_SLOW > );
  
  //dendr
  insert_( names::V_m_d, &compte2003_ex::get_y_elem_< compte2003_ex::State_::V_M_D > );
  insert_( names::n_Ca, &compte2003_ex::get_y_elem_< compte2003_ex::State_::N_CA > );
  insert_( names::g_gaba, &compte2003_ex::get_y_elem_< compte2003_ex::State_::G_GABA > );
}
}


extern "C" int
nest::compte2003_ex_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef nest::compte2003_ex::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::compte2003_ex& node = *( reinterpret_cast< nest::compte2003_ex* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].
 
  
  // membrane potential
  const double_t& V_s = y[ S::V_M_S ];
  const double_t& V_d = y[ S::V_M_D ];

  
  
  //-----somatic currents
  
  
  //-----inhibitory synaptic current-----
  double g_syn_inh = y[S::G_GABA]; 
  
  double I_inh = g_syn_inh*( node.P_.E_in-V_s);
  
  //-----fast Na+ current-----
  double Na_m_inf; 
  double Na_h_inf; 
  double Na_h_tau;
  double a,b;
  //activation
  a=0.1*(V_s+33.)/(1.-exp(-(V_s+33.)/10.)); 
  b=4.*exp(-(V_s+53.7)/12.);
  Na_m_inf=a/(a+b); 
 
  //inactivation   
  a=0.07*exp(-(V_s+50.)/10.);
  b=1./(1.+exp(-(V_s+20.)/10.));
  Na_h_inf=a/(a+b); 
  Na_h_tau=0.25/(a+b);    //[ms]
      
  double I_Na = node.P_.g_Na*pow(Na_m_inf,3)*y[S::NA_H]*(node.P_.E_Na-V_s); //[pA]

  //-----fast K+ current
  double K_m_inf; 
  double K_m_tau;
  a=0.01*(V_s+34.)/(1-exp(-(V_s+34.)/10.));
  b=0.125*exp(-(V_s+44.)/25.);
  K_m_inf=a/(a+b); 
  K_m_tau=0.25/(a+b);     //[ms]

  double I_K = node.P_.g_K*pow(y[S::K_M],4)*(node.P_.E_K-V_s);  //[pA]

  //-----K+ A-current-----
  //activation
  double K_A_m=1./(1.+exp(-(V_s+50.)/20.));
  //inactivation  
  double K_A_h_inf=1./(1.+exp((V_s+80.)/6.));
  double K_A_h_tau=15.;                  //[ms]
     
  double I_K_A = node.P_.g_K_A*pow(K_A_m,3)*y[S::K_A_H]*(node.P_.E_K-V_s);  //[pA]    
  
  //-----K+ slow non-inactivating current
  double K_s_m_inf=1./(1.+exp(-(V_s+34.)/6.5));
  double K_s_m_tau=8./(exp(-(V_s+55.)/30.)+exp((V_s+55.)/30.));    //[ms]
     
  double I_K_s = node.P_.g_K_s*y[S::K_S_M]*(node.P_.E_K-V_s);   //[pA]  
  
  //-----K(Na) current
  double n_Na=y[S::N_NA];
  double K_Na_m_inf=0.37/(1.+pow(38.7/n_Na,3.5));
  double I_K_Na = node.P_.g_K_Na*K_Na_m_inf*(node.P_.E_K-V_s);     //[pA]  
  
  //-----leak current
  double I_L = node.P_.g_L*(node.P_.E_L-V_s);
     
  
     
  //-----dendritic currents
  
        
  //-----excitatory synaptic currents-----
  double g_syn_exc = y[S::G_AMPA]+(y[S::G_NMDA_SLOW]-y[S::G_NMDA_FAST]);

  double I_exc = g_syn_exc*( node.P_.E_ex - V_d);
              
  //-----Na+ persistent current
  double Na_p_m_inf=1./(1.+exp(-(V_d+55.7)/7.7));
    
  double I_Na_p = node.P_.g_Na_p*pow(Na_p_m_inf,3)*(node.P_.E_Na-V_d);
     
  //-----K+ inward-rectifying current
  double K_AR_m_inf=1./(1.+exp((V_d+75.)/4.));
    
  double I_K_AR = node.P_.g_K_AR*K_AR_m_inf*(node.P_.E_K-V_d); //[pA]
                    
  //-----Ca2+ current   
  double Ca_m_inf=1./(1.+exp(-(V_d+20.)/9.));
    
  double I_Ca = node.P_.g_Ca*pow(Ca_m_inf,2)*(node.P_.E_Ca-V_d); //[pA]
     
  //-----K(Ca) current
  double n_Ca=y[S::N_CA];
  double K_Ca_m_inf=n_Ca/(n_Ca+30.);
  double I_K_Ca = node.P_.g_K_Ca*K_Ca_m_inf*(node.P_.E_K-V_d);     //[pA]


  //-----coupling currents between segments-----
  double I_conn = node.P_.g_conn * ( V_d - V_s ); 


  //-----derivatives
  // membrane potential
  f[S::V_M_S] = ( I_L + I_Na + I_K+ I_K_s + I_K_A + I_K_Na + I_inh + I_conn + node.B_.I_stim_ + node.P_.I_e ) / node.P_.C_m_s;    
  f[S::V_M_D] = ( I_Na_p + I_K_AR+ I_Ca + I_K_Ca + I_exc - I_conn ) / node.P_.C_m_d;     
  /*  
  if (y[S::V_M_S] < -100.)
    {printf("V_m_s=%f",y[S::V_M_S]);
     printf("Na_m_inf=%f",Na_m_inf);
     printf("Na_h_inf=%f",Na_h_inf);
     printf("Na_h_tau=%f",Na_h_tau);
     printf("K_m_inf=%f" ,K_m_inf);
     printf("K_m_tau=%f" ,K_m_tau); 
     printf("f_V_s=%f\n" ,f[S::V_M_S]);   
    }
  */  
  if (y[S::V_M_S] < -200.)
    {printf("V_m_s=%f ",y[S::V_M_S]);
     printf("f_V_s=%f ",f[S::V_M_S]);
     printf("V_m_d=%f ",y[S::V_M_D]);
     printf("f_V_d=%f ",f[S::V_M_D]);
     printf("I_L=%f ",I_L);
     printf("I_Na=%f ",I_Na);
     printf("I_K=%f ",I_K);
     printf("I_K_s=%f ",I_K_s);
     printf("I_K_Na=%f " ,I_K_Na); 
     printf("I_conn=%f " ,I_conn);  
     printf("I_Na_p=%f ",I_Na_p);
     printf("I_K_AR=%f ",I_K_AR);
     printf("I_Ca=%f ",I_Ca);
     printf("I_K_Ca=%f\n",I_K_Ca);
    }
  // Ca2+ amount
  double alpha_Ca=(0.005)*1E-3;   //[uM/pA/ms] - 1E-3 is a transition factor from original nA to present pA
  double tau_Ca=150.;             //[ms]       - Ca concentration decay time
  f[S::N_CA] = alpha_Ca*I_Ca-n_Ca/tau_Ca;
     
  // Na+ amount
  double alpha_Na=(0.01)*1E-3;   //[mM/pA/ms] - 1E-3 is a transition factor from original nA to present pA
  double R_pump=0.018;           //[mM/ms]
  double Na_eq=9.5;              //[mM]
  f[S::N_NA] = alpha_Na*(I_Na+I_Na_p)-R_pump*(pow(n_Na,3)/(pow(n_Na,3)+3375.)-pow(Na_eq,3)/(pow(Na_eq,3)+3375));
    
     
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
  // K+ A-current inactivation
  f[S::K_A_H]  = (K_A_h_inf  - y[S::K_A_H])  / K_A_h_tau;      
  // K+ slow non-inactivating activation
  f[S::K_S_M]   = (K_s_m_inf   - y[S::K_S_M])   / K_s_m_tau; 
           
       
   return GSL_SUCCESS;
   
   
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::compte2003_ex::Parameters_::Parameters_()
  : t_ref   (2.0 ), // ms
    E_ex    (0.0),  // mV
    E_in    (-85.0),  // mV
    E_L     (-70.0),  // mV
    E_Na     (50.0),  // mV
    E_K     (-90.0),  // mV
    E_Ca    (140.0),  // mV
    
    tau_syn_ampa  (3.),  // ms
    tau_syn_nmda_fast  (5.), // ms
    tau_syn_nmda_slow  (50.), // ms
    tau_syn_gaba        (3.),  // ms

    g_conn  (2.5),  // nS,  conductances between compartments
    
    //soma
    g_Na    (0.0),  // nS
    g_K     (0.0),  // nS
    g_K_s   (0.0),  // nS
    g_K_A   (0.0),  // nS
    g_K_Na  (0.0),  // nS
    g_L     (10.0),  // nS
    C_m_s   (150.0),  // pF        
    I_e     (0.0),  // pA
    
    //dendr    
    g_Na_p  (0.0),  // nS      
    g_K_AR  (0.0),  // nS
    g_Ca    (0.0),  // nS
    g_K_Ca  (0.0),  // nS
    C_m_d   (150.0)  // pF


{
}


nest::compte2003_ex::State_::State_(const Parameters_ & p)
  : r_(0)
{
  // for simplicity, we first initialize all values to 0,
  // then set the membrane potentials for each compartment
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = 0;

  double V=-80.;
   
  //Fast Na+ current inactivation
  double a=0.07*exp(-(V+50.)/10.);
  double b=1./(1.+exp(-(V+20.)/10.));
  double Na_h_inf=a/(a+b); 

  //Fast K+ current
  a=0.01*(V+34.)/(1-exp(-(V+34.)/10.));
  b=0.125*exp(-(V+44.)/25.);
  double K_m_inf=a/(a+b); 
   
  //K+ A-current inactivation
  double K_A_h_inf=1./(1.+exp((V+80.)/6.));

  //K+ slow non-inactivating current
  double K_s_m_inf=1./(1.+exp(-(V+34.)/6.5));
   
  y_[0] = V;    //V_m_s
  y_[1] = V;    //V_m_d
  y_[6]= 9.5;   //n_Ca
    
  y_[8] = Na_h_inf; //Na_h
  y_[9]  = K_m_inf; //K_m
  y_[10] = K_s_m_inf; //K_s_m   
  y_[11] = K_A_h_inf; //K_A_h
 
}


nest::compte2003_ex::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
}

nest::compte2003_ex::State_& nest::compte2003_ex::State_::operator=( const State_& s )
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

void nest::compte2003_ex::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::t_ref,      t_ref);
  def<double>(d,names::E_ex,       E_ex);
  def<double>(d,names::E_in,       E_in);
  def<double>(d,names::E_L,        E_L); 
  def<double>(d,names::E_Na,       E_Na);
  def<double>(d,names::E_K ,       E_K );
  def<double>(d,names::E_Ca,       E_Ca);


  def<double>(d,names::tau_syn_ampa, tau_syn_ampa);
  def<double>(d,names::tau_syn_nmda_fast, tau_syn_nmda_fast);
  def<double>(d,names::tau_syn_nmda_slow, tau_syn_nmda_slow);
  def<double>(d,names::tau_syn_gaba, tau_syn_gaba);
    
  def<double>(d,names::g_conn, g_conn);
  //soma
  def<double>(d, names::g_L,        g_L);
  def<double>(d, names::g_Na,       g_Na);
  def<double>(d, names::g_K,        g_K);
  def<double>(d, names::g_K_s ,     g_K_s);
  def<double>(d, names::g_K_A ,     g_K_A);
  def<double>(d, names::g_K_Na ,    g_K_Ca);
  def<double>(d, names::C_m_s,      C_m_s);
  def<double>(d, names::I_e,        I_e);
  //dendr  
  def<double>(d, names::g_Na_p,     g_Na_p);
  def<double>(d, names::g_K_AR,     g_K_AR);
  def<double>(d, names::g_Ca ,      g_Ca);
  def<double>(d, names::g_K_Ca ,    g_K_Ca);
  def<double>(d, names::C_m_d,      C_m_d);

  
}



void nest::compte2003_ex::Parameters_::set(const DictionaryDatum &d)
{
  // allow setting the membrane potential
  updateValue<double>(d,names::t_ref,     t_ref);
  updateValue<double>(d,names::E_ex,       E_ex);
  updateValue<double>(d,names::E_in,       E_in);
  updateValue<double>(d,names::E_L,        E_L); 
  updateValue<double>(d,names::E_Na,       E_Na);
  updateValue<double>(d,names::E_K ,       E_K );
  updateValue<double>(d,names::E_Ca ,      E_Ca );


  updateValue<double>(d,names::tau_syn_ampa, tau_syn_ampa);
  updateValue<double>(d,names::tau_syn_nmda_fast, tau_syn_nmda_fast);
  updateValue<double>(d,names::tau_syn_nmda_slow, tau_syn_nmda_slow);
  updateValue<double>(d,names::tau_syn_gaba, tau_syn_gaba);
           
  updateValue<double>(d,names::g_conn, g_conn);
  //soma
  updateValue<double>(d, names::g_L,        g_L);
  updateValue<double>(d, names::g_Na,       g_Na);
  updateValue<double>(d, names::g_K,        g_K);
  updateValue<double>(d, names::g_K_s,      g_K_s);
  updateValue<double>(d, names::g_K_A,      g_K_A);
  updateValue<double>(d, names::g_K_Na,     g_K_Na);
  updateValue<double>(d, names::C_m_s,      C_m_s);
  updateValue<double>(d, names::I_e,        I_e);
  //dendr  
  updateValue<double>(d, names::g_Na_p,     g_Na_p);
  updateValue<double>(d, names::g_K_AR,     g_K_AR);
  updateValue<double>(d, names::g_Ca,       g_Ca);
  updateValue<double>(d, names::g_K_Ca,     g_K_Ca);
  updateValue<double>(d, names::C_m_d,        C_m_d);


  // apply checks 
  if ( C_m_s <= 0  || C_m_d <=0 )
      throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref < 0 )
      throw BadProperty("Refractory time cannot be negative.");
          
    if ( tau_syn_ampa <= 0 || tau_syn_nmda_fast <= 0 || tau_syn_nmda_slow <= 0 || tau_syn_gaba <= 0 )
      throw BadProperty("All time constants must be strictly positive.");
  
}

void nest::compte2003_ex::State_::get(DictionaryDatum &d) const
{

  // we assume here that State_::get() always is called after Parameters_::get(),
  // so that the per-compartment dictionaries exist
  def<double>(d, names::V_m_s, y_[V_M_S]); // Membrane potential
  def<double>(d, names::V_m_d, y_[V_M_D]); // Membrane potential
  def<double>(d, names::n_Ca, y_[N_CA]); // 
  def<double>(d, names::n_Na, y_[N_NA]); // 
  def<double>(d, names::Na_h, y_[NA_H]); // 
  def<double>(d, names::K_m , y_[K_M]); // 
  def<double>(d, names::K_s_m, y_[K_S_M]); // 
  def<double>(d, names::K_A_h , y_[K_A_H]); // 
  def<double>(d, names::g_ampa , y_[G_AMPA]); // 
  def<double>(d, names::g_nmda_slow , y_[G_NMDA_SLOW]); // 
  def<double>(d, names::g_nmda_fast , y_[G_NMDA_FAST]); // 
  def<double>(d, names::g_gaba , y_[G_GABA]); // 
 
}

void nest::compte2003_ex::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  // extract from sub-dictionaries
  updateValue<double>(d, names::V_m_s, y_[V_M_S]);
  updateValue<double>(d, names::V_m_d, y_[V_M_D]);
  updateValue<double>(d, names::n_Ca, y_[N_CA]);
  updateValue<double>(d, names::n_Na, y_[N_NA]);
  updateValue<double>(d, names::Na_h, y_[NA_H]);
  updateValue<double>(d, names::K_m , y_[K_M]);
  updateValue<double>(d, names::K_s_m, y_[K_S_M]);
  updateValue<double>(d, names::K_A_h , y_[K_A_H]);
  updateValue<double>(d, names::g_ampa , y_[G_AMPA]);
  updateValue<double>(d, names::g_nmda_fast , y_[G_NMDA_FAST]);
  updateValue<double>(d, names::g_nmda_slow , y_[G_NMDA_SLOW]);
  updateValue<double>(d, names::g_gaba , y_[G_GABA]);

  if ( y_[ G_AMPA ] < 0 || y_[ G_NMDA_FAST ] < 0 || y_[ G_NMDA_SLOW ] < 0 || y_[ G_GABA ] < 0  )
    throw BadProperty( "Conductances must not be negative." );

}


nest::compte2003_ex::Buffers_::Buffers_( compte2003_ex& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::compte2003_ex::Buffers_::Buffers_( const Buffers_&, compte2003_ex& n )
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

nest::compte2003_ex::compte2003_ex()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

nest::compte2003_ex::compte2003_ex( const compte2003_ex& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

nest::compte2003_ex::~compte2003_ex()
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
nest::compte2003_ex::init_state_( const Node& proto )
{
  const compte2003_ex& pr = downcast< compte2003_ex >( proto );
  S_ = pr.S_;
}

void
nest::compte2003_ex::init_buffers_()
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

  B_.sys_.function = compte2003_ex_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::compte2003_ex::calibrate()
{
  B_.logger_.init(); // ensures initialization in case mm connected after Simulate
  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref ) ).get_steps();
  assert( V_.RefractoryCounts_ >= 0 ); // since t_ref >= 0, this can only fail in error
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::compte2003_ex::update( const Time& origin, const long_t from, const long_t to )
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
nest::compte2003_ex::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );
  assert( 0 <= e.get_rport() && e.get_rport() < 5 );

  B_.spikes_[ e.get_rport() ].add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
}

void
nest::compte2003_ex::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  // add weighted current; HEP 2002-10-04
  B_.currents_[ e.get_rport() ].add_value(
    e.get_rel_delivery_steps( network()->get_slice_origin() ), e.get_weight() * e.get_current() );
}



void
nest::compte2003_ex::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL_1_11
