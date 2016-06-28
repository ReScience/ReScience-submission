/* -----------------------------------------------------------------------------
 * Distributed under the GNU General Public License.
 *
 * Contributors: Andrei Maksimov (maksimov.andrei7@gmail.com)
 *-----------------------------------------------------------------------------
 * References:
 *
 * Cellular and network mechanisms of slow oscillatory activity (<1 Hz) 
 * and wave propagations in a cortical network model*, A. Compte, 
 * M.V. Sanchez-Vives, D.A. McCormick, X.-J. Wang, 
 * Journal of Neurophysiology, 2707--2725, 2003"
 *------------------------------------------------------------------------------
 *
 * This model is based on the "iaf_cond_alpha_mc" model in NEST
 * written by Hans Ekkehard Plesser
 *------------------------------------------------------------------------------
*/ 

#include "compte2003_ex.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <iomanip>
#include <iostream>
#include <cstdio>

/* ----------------------------------------------------------------
 * Compartment name list
 * ---------------------------------------------------------------- */

/* Harold Gutch reported some static destruction problems on OSX 10.4.
   He pointed out that the problem is avoided by defining the comp_names_
   vector with its final size. See also #348.
*/
std::vector< Name > nest::compte2003_ex::comp_names_( NCOMP );

/* ----------------------------------------------------------------
 * Receptor dictionary
 * ---------------------------------------------------------------- */

// leads to seg fault on exit, see #328
// DictionaryDatum nest::compte2003_ex::receptor_dict_ = new Dictionary();

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::compte2003_ex > nest::compte2003_ex::recordablesMap_;

namespace nest
{
// specialization must be place in namespace

template <>
void
RecordablesMap< compte2003_ex >::create()
{
//compartment 0
  insert_(Name("V_m.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::V_M,compte2003_ex::comp_0>);
  insert_(Name("n_Ca.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::n_Ca,compte2003_ex::comp_0>);
  insert_(Name("n_Na.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::n_Na,compte2003_ex::comp_0>);
  insert_(Name("Na_h.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::Na_h,compte2003_ex::comp_0>);
  insert_(Name("K_m.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::K_m,compte2003_ex::comp_0>);
  insert_(Name("K_s_m.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::K_s_m,compte2003_ex::comp_0>);
  insert_(Name("K_A_h.0"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::K_A_h,compte2003_ex::comp_0>);
  insert_(Name("g_ampa.0"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_AMPA,compte2003_ex::comp_0>);
  insert_(Name("g_nmda_fast.0"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_NMDA_FAST,compte2003_ex::comp_0>);
  insert_(Name("g_nmda_slow.0"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_NMDA_SLOW,compte2003_ex::comp_0>);
  insert_(Name("g_gaba.0"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_GABA,compte2003_ex::comp_0>);

  //compartment 1
  insert_(Name("V_m.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::V_M,compte2003_ex::comp_1>);
  insert_(Name("n_Ca.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::n_Ca,compte2003_ex::comp_1>);
  insert_(Name("n_Na.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::n_Na,compte2003_ex::comp_1>);
  insert_(Name("Na_h.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::Na_h,compte2003_ex::comp_1>);
  insert_(Name("K_m.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::K_m,compte2003_ex::comp_1>);
  insert_(Name("K_s_m.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::K_s_m,compte2003_ex::comp_1>);
  insert_(Name("K_A_h.1"), &compte2003_ex::get_y_elem_<compte2003_ex::State_::K_A_h,compte2003_ex::comp_1>);
  insert_(Name("g_ampa.1"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_AMPA,compte2003_ex::comp_1>);
  insert_(Name("g_nmda_fast.1"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_NMDA_FAST,compte2003_ex::comp_1>);
  insert_(Name("g_nmda_slow.1"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_NMDA_SLOW,compte2003_ex::comp_1>);
  insert_(Name("g_gaba.1"),&compte2003_ex::get_y_elem_<compte2003_ex::State_::G_GABA,compte2003_ex::comp_1>);

}
}

/* ----------------------------------------------------------------
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C" int
nest::compte2003_ex_dynamics( double, const double y[], double f[], void* pnode )
{
  // some shorthands
  typedef nest::compte2003_ex N;
  typedef nest::compte2003_ex::State_ S;

  // get access to node so we can work almost as in a member function
  assert( pnode );
  const nest::compte2003_ex& node = *( reinterpret_cast< nest::compte2003_ex* >( pnode ) );

  // compute dynamics for each compartment
  // computations written quite explicitly for clarity, assume compile
  // will optimized most stuff away ...
   
   
  //for ( size_t n = 0; n < N::NCOMP; ++n )
  //for ( int n = 0; n<=N::NCOMP-1 ;  ++n )
  for ( int n = N::NCOMP-1 ; n >=0 ; --n )
   {
     // membrane potential for current compartment
     const double V = y[S::idx(n, S::V_M)];

     //-----excitatory synaptic currents-----
     const double g_syn_exc = y[S::idx(n, S::G_AMPA)]+(y[S::idx(n, S::G_NMDA_SLOW)]-y[S::idx(n, S::G_NMDA_FAST)]);
     
     //-----inhibitory synaptic current-----
     const double g_syn_inh = y[S::idx(n, S::G_GABA)];


     //-----coupling currents between segments-----
     double I_conn =
       (n>0               ? node.P_.g_conn[n-1] * ( V - y[S::idx(n-1, S::V_M)] ) : 0. ) //connect proximal neighbor segment
      +(n<(N::NCOMP-1)    ? node.P_.g_conn[n] * ( V - y[S::idx(n+1, S::V_M)] ) : 0. ); //connect distal neighbor segment
     
     //-----fast Na+ current-----
     double Na_m_inf; 
     double Na_h_inf; 
     double Na_h_tau;
     double a,b;
     //activation
     a=0.1*(V+33.)/(1.-exp(-(V+33.)/10.));
     b=4.*exp(-(V+53.7)/12.);
     Na_m_inf=a/(a+b); 
     
     //inactivation   
     a=0.07*exp(-(V+50.)/10.);
     b=1./(1.+exp(-(V+20.)/10.));
     Na_h_inf=a/(a+b); 
     Na_h_tau=0.25/(a+b);    //[ms]
     
     double I_Na = node.P_.g_Na[n]*pow(Na_m_inf,3)*y[S::idx(n, S::Na_h)]*(node.P_.E_Na-V); //[pA]
     
     //-----fast K+ current
     double K_m_inf; 
     double K_m_tau;
     a=0.01*(V+34.)/(1-exp(-(V+34.)/10.));
     b=0.125*exp(-(V+44.)/25.);
     K_m_inf=a/(a+b); 
     K_m_tau=0.25/(a+b);     //[ms]

     double I_K = node.P_.g_K[n]*pow(y[S::idx(n, S::K_m)],4)*(node.P_.E_K-V);  //[pA]

     //-----K+ A-current-----
     //activation
     double K_A_m=1./(1.+exp(-(V+50.)/20.));
     //inactivation  
     double K_A_h_inf=1./(1.+exp((V+80.)/6.));
     double K_A_h_tau=15.;                  //[ms]
     
     double I_K_A = node.P_.g_K_A[n]*pow(K_A_m,3)*y[S::idx(n, S::K_A_h)]*(node.P_.E_K-V);  //[pA]    
          
     //-----Na+ persistent current
     double Na_p_m_inf=1./(1+exp(-(V+55.7)/7.7));
    
     double I_Na_p = node.P_.g_Na_p[n]*pow(Na_p_m_inf,3)*(node.P_.E_Na-V);
     
     //-----K+ slow non-inactivating current
     double K_s_m_inf=1./(1.+exp(-(V+34.)/6.5));
     double K_s_m_tau=8./(exp(-(V+55.)/30.)+exp((V+55.)/30.));    //[ms]
     
     double I_K_s = node.P_.g_K_s[n]*y[S::idx(n, S::K_s_m)]*(node.P_.E_K-V);   //[pA]

     //-----K+ inward-rectifying current
     double K_AR_m_inf=1./(1+exp((V+75.)/4.));
    
     double I_K_AR = node.P_.g_K_AR[n]*K_AR_m_inf*(node.P_.E_K-V); //[pA]
          
          
     //-----Ca2+ current   
     double Ca_m_inf=1./(1.+exp(-(V+20.)/9.));
    
     double I_Ca = node.P_.g_Ca[n]*pow(Ca_m_inf,2)*(node.P_.E_Ca-V); //[pA]

     
     //-----K(Ca) current
     double n_Ca=y[S::idx(n, S::n_Ca)];
     double K_Ca_m_inf=n_Ca/(n_Ca+30.);
     double I_K_Ca = node.P_.g_K_Ca[n]*K_Ca_m_inf*(node.P_.E_K-V);     //[pA]

     //-----K(Na) current
     double n_Na=y[S::idx(n, S::n_Na)];
     double K_Na_m_inf=0.37/(1.+pow(38.7/n_Na,3.5));
     double I_K_Na = node.P_.g_K_Na[n]*K_Na_m_inf*(node.P_.E_K-V);     //[pA]

     //-----leak current
     double I_L = node.P_.g_L[n]*(node.P_.E_L-V);

     // derivatives
     // membrane potential
     f[S::idx(n, S::V_M)] = ( I_L + I_Na + I_K+ I_Na_p + I_K_s + I_K_AR+ I_K_A+I_Ca+I_K_Ca+I_K_Na+g_syn_exc*( node.P_.E_ex-V) + g_syn_inh*( node.P_.E_in-V) - I_conn + node.B_.I_stim_[n] + node.P_.I_e[n] ) / node.P_.C_m[n];     
     
     // Ca2+ amount
     double alpha_Ca=(0.005)*1E-3;   //[uM/pA/ms] - 1E-3 is a transition factor from original nA to present pA
     double tau_Ca=150.;             //[ms]       - Ca concentration decay time
     f[S::idx(n, S::n_Ca)] = alpha_Ca*I_Ca-n_Ca/tau_Ca;
     
     // Na+ amount
     double alpha_Na=(0.01)*1E-3;   //[mM/pA/ms] - 1E-3 is a transition factor from original nA to present pA
     double R_pump=0.018;           //[mM/ms]
     double Na_eq=9.5;              //[mM]
     f[S::idx(n, S::n_Na)] = alpha_Na*(I_Na+I_Na_p)-R_pump*(pow(n_Na,3)/(pow(n_Na,3)+3375.)-pow(Na_eq,3)/(pow(Na_eq,3)+3375));
    
     
     // excitatory conductance
     f[S::idx(n, S::G_AMPA)] = -y[S::idx(n, S::G_AMPA)] / node.P_.tau_syn_ampa; // Synaptic Conductance (nS)
     f[S::idx(n, S::G_NMDA_FAST)] = -y[S::idx(n, S::G_NMDA_FAST)] / node.P_.tau_syn_nmda_fast; // Synaptic Conductance (nS)
     f[S::idx(n, S::G_NMDA_SLOW)] = -y[S::idx(n, S::G_NMDA_SLOW)] / node.P_.tau_syn_nmda_slow; // Synaptic Conductance (nS)
    
     // inhibitory conductance
     f[S::idx(n, S::G_GABA)] = -y[S::idx(n, S::G_GABA)] / node.P_.tau_syn_gaba; // Synaptic Conductance (nS)
     
     // fast Na+ inactivation
     f[S::idx(n,S::Na_h)] = (Na_h_inf - y[S::idx(n, S::Na_h)]) / Na_h_tau;  
     // fast K+ activation
     f[S::idx(n,S::K_m)]  = (K_m_inf  - y[S::idx(n, S::K_m)])  / K_m_tau;    
     // K+ A-current inactivation
     f[S::idx(n,S::K_A_h)]  = (K_A_h_inf  - y[S::idx(n, S::K_A_h)])  / K_A_h_tau;      
     // K+ slow non-inactivating activation
     f[S::idx(n,S::K_s_m)]   = (K_s_m_inf   - y[S::idx(n, S::K_s_m)])   / K_s_m_tau; 
           
       
       
   }
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
    tau_syn_gaba        (3.)  // ms

{
    // conductances between compartments
  for ( size_t n = 0 ; n < NCOMP-1 ; ++n )
    g_conn[n] = 2.5;  // nS, soma-proximal


    // compartment parameters
  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    g_Na     [n] =  0.0;  // nS
    g_K      [n] =  0.0;  // nS
    g_Na_p   [n] =  0.0;  // nS      
    g_K_s    [n] =  0.0;  // nS
    g_K_A    [n] =  0.0;  // nS
    g_K_AR   [n] =  0.0;  // nS
    g_Ca     [n] =  0.0;  // nS
    g_K_Ca   [n] =  0.0;  // nS
    g_K_Na   [n] =  0.0;  // nS
    g_L      [n] =  10.0;  // nS
    C_m      [n] = 150.0;  // pF
    I_e      [n] =   0.0;  // pA
  }


}

nest::compte2003_ex::Parameters_::Parameters_(const Parameters_& p)
  : t_ref (p.t_ref),
    E_ex  (p.E_ex),
    E_in  (p.E_in),
    E_Na  (p.E_Na),
    E_K   (p.E_K),     
    E_L   (p.E_L),   
    E_Ca  (p.E_Ca),   

    tau_syn_ampa (p.tau_syn_ampa),
    tau_syn_nmda_fast (p.tau_syn_nmda_fast),
    tau_syn_nmda_slow (p.tau_syn_nmda_slow),
    tau_syn_gaba      (p.tau_syn_gaba)
{
  // copy C-arrays
  for ( size_t n = 0 ; n < NCOMP-1 ; ++n )
    g_conn[n] = p.g_conn[n];

  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    g_Na    [n] = p.g_Na    [n];
    g_K     [n] = p.g_K     [n];
    g_Na_p  [n] = p.g_Na_p  [n];
    g_K_s   [n] = p.g_K_s   [n];
    g_K_A   [n] = p.g_K_A   [n];
    g_K_AR  [n] = p.g_K_AR  [n];
    g_Ca    [n] = p.g_Ca    [n];
    g_K_Ca  [n] = p.g_K_Ca  [n];
    g_K_Na  [n] = p.g_K_Na  [n];
    g_L     [n] = p.g_L     [n];
    C_m     [n] = p.C_m     [n];
    I_e     [n] = p.I_e     [n];
  }
}


nest::compte2003_ex::Parameters_& nest::compte2003_ex::Parameters_::operator=(const Parameters_& p)
{
  assert(this != &p);  // would be bad logical error in program

  t_ref   = p.t_ref;
  E_ex    = p.E_ex;
  E_in    = p.E_in;
  E_Na    = p.E_Na;
  E_K     = p.E_K ;  
  E_Ca    = p.E_Ca;
  E_L     = p.E_L ;

  tau_syn_ampa = p.tau_syn_ampa;
  tau_syn_nmda_fast = p.tau_syn_nmda_fast;
  tau_syn_nmda_slow = p.tau_syn_nmda_slow;
  tau_syn_gaba      = p.tau_syn_gaba;
       
  // copy C-arrays
  for ( size_t n = 0 ; n < NCOMP-1 ; ++n )
    g_conn[n] = p.g_conn[n];

  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    g_Na    [n] = p.g_Na    [n];
    g_K     [n] = p.g_K     [n];
    g_Na_p  [n] = p.g_Na_p  [n];
    g_K_s   [n] = p.g_K_s   [n];
    g_K_A   [n] = p.g_K_A   [n];
    g_K_AR  [n] = p.g_K_AR  [n];
    g_Ca    [n] = p.g_Ca    [n];
    g_K_Ca  [n] = p.g_K_Ca  [n];
    g_K_Na  [n] = p.g_K_Na  [n];
    g_L     [n] = p.g_L     [n];
    C_m     [n] = p.C_m     [n];
    I_e     [n] = p.I_e     [n];
  }

  return *this;
}


nest::compte2003_ex::State_::State_(const Parameters_& p)
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

     
     
  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    y_[idx(n, V_M)] = V;
    y_[idx(n, n_Ca)]= 9.5;
    
    y_[idx(n, Na_h)] = Na_h_inf;
    y_[idx(n, K_m)] = K_m_inf;
    y_[idx(n, K_s_m)] = K_s_m_inf;   
    y_[idx(n, K_A_h)] = K_A_h_inf;   
  }
}


nest::compte2003_ex::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
}

nest::compte2003_ex::State_& nest::compte2003_ex::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

nest::compte2003_ex::Buffers_::Buffers_(compte2003_ex& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::compte2003_ex::Buffers_::Buffers_(const Buffers_&, compte2003_ex& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}
/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::compte2003_ex::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::t_ref,     t_ref);
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
    
  for ( size_t n = 0 ; n < NCOMP-1 ; ++n )  
    def<double>(d,Name(String::compose("g_%1",n)), g_conn[n]);

  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    DictionaryDatum dd = new Dictionary();
    def<double>(dd, names::g_L,        g_L[n]);
    def<double>(dd, names::g_Na,       g_Na[n]);
    def<double>(dd, names::g_K,        g_K[n]);
    def<double>(dd, names::g_Na_p,     g_Na_p[n]);
    def<double>(dd, names::g_K_s ,     g_K_s[n]);
    def<double>(dd, names::g_K_A ,     g_K_A[n]);
    def<double>(dd, names::g_K_AR,     g_K_AR[n]);
    def<double>(dd, names::g_Ca ,      g_Ca[n]);
    def<double>(dd, names::g_K_Ca ,    g_K_Ca[n]);
    def<double>(dd, names::g_K_Na ,    g_K_Ca[n]);
    def<double>(dd, names::C_m,        C_m[n]);
    def<double>(dd, names::I_e,        I_e[n]);

    (*d)[comp_names_[n]] = dd;
  }
}

void nest::compte2003_ex::Parameters_::set(const DictionaryDatum& d)
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

                  
  for ( size_t n = 0 ; n < NCOMP-1 ; ++n )  
    updateValue<double>(d,Name(String::compose("g_%1",n)), g_conn[n]);

  // extract from sub-dictionaries
  for ( size_t n = 0 ; n < NCOMP ; ++n )
    if ( d->known(comp_names_[n]) )
    {
      DictionaryDatum dd = getValue<DictionaryDatum>(d, comp_names_[n]);
      updateValue<double>(dd, names::g_L,        g_L[n]);
      updateValue<double>(dd, names::g_Na,       g_Na[n]);
      updateValue<double>(dd, names::g_K,        g_K[n]);
      updateValue<double>(dd, names::g_Na_p,     g_Na_p[n]);
      updateValue<double>(dd, names::g_K_s,      g_K_s[n]);
      updateValue<double>(dd, names::g_K_A,      g_K_A[n]);
      updateValue<double>(dd, names::g_K_AR,     g_K_AR[n]);
      updateValue<double>(dd, names::g_Ca,       g_Ca[n]);
      updateValue<double>(dd, names::g_K_Ca,     g_K_Ca[n]);
      updateValue<double>(dd, names::g_K_Na,     g_K_Na[n]);
      updateValue<double>(dd, names::C_m,        C_m[n]);
      updateValue<double>(dd, names::I_e,        I_e[n]);


    }

 // apply checks compartment-wise
  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    if ( C_m[n] <= 0 )
      throw BadProperty("Capacitance (" 
			+ comp_names_[n].toString() + ") must be strictly positive.");
    
    if ( t_ref < 0 )
      throw BadProperty("Refractory time cannot be negative.");
          
    if ( tau_syn_ampa <= 0 || tau_syn_nmda_fast <= 0 || tau_syn_nmda_slow <= 0 || tau_syn_gaba <= 0 )
      throw BadProperty("All time constants ("
			+ comp_names_[n].toString() + ") must be strictly positive.");
  }
}

void nest::compte2003_ex::State_::get(DictionaryDatum &d) const
{

  // we assume here that State_::get() always is called after Parameters_::get(),
  // so that the per-compartment dictionaries exist
  for ( size_t n = 0 ; n < NCOMP ; ++n )
  {
    assert(d->known(comp_names_[n]));
    DictionaryDatum dd = getValue<DictionaryDatum>(d, comp_names_[n]);

    def<double>(dd, names::V_m , y_[idx(n, V_M)]); // Membrane potential
    def<double>(dd, names::n_Ca, y_[idx(n, n_Ca)]); // 
    def<double>(dd, names::n_Na, y_[idx(n, n_Na)]); // 
    def<double>(dd, names::Na_h, y_[idx(n, Na_h)]); // 
    def<double>(dd, names::K_m , y_[idx(n, K_m )]); // 
    def<double>(dd, names::K_s_m, y_[idx(n, K_s_m)]); // 
    def<double>(dd, names::K_A_h , y_[idx(n, K_A_h)]); // 
  }
}

void nest::compte2003_ex::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  // extract from sub-dictionaries
  for ( size_t n = 0 ; n < NCOMP ; ++n )
    if ( d->known(comp_names_[n]) )
    {
      DictionaryDatum dd = getValue<DictionaryDatum>(d, comp_names_[n]);
      updateValue<double>(dd, names::V_m, y_[idx(n, V_M)]);
      updateValue<double>(dd, names::n_Ca, y_[idx(n, n_Ca)]);
      updateValue<double>(dd, names::n_Na, y_[idx(n, n_Na)]);
      updateValue<double>(dd, names::Na_h, y_[idx(n, Na_h)]);
      updateValue<double>(dd, names::K_m , y_[idx(n, K_m )]);
      updateValue<double>(dd, names::K_s_m, y_[idx(n, K_s_m)]);
      updateValue<double>(dd, names::K_A_h , y_[idx(n, K_A_h )]);
 
    }
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

  // set up table of compartment names
  // comp_names_.resize(NCOMP); --- Fixed size, see comment on definition
  for ( size_t n = 0 ; n < NCOMP ; ++n )
    comp_names_[n] = Name(String::compose("comp_%1",n));
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
    gsl_odeiv_step_free( B_.s_ );
  if ( B_.c_ )
    gsl_odeiv_control_free( B_.c_ );
  if ( B_.e_ )
    gsl_odeiv_evolve_free( B_.e_ );
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

  B_.logger_.reset();
  Archiving_Node::clear_history();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_step_reset( B_.s_ );

  if ( B_.c_ == 0 )
    B_.c_ = gsl_odeiv_control_y_new( 1e-6, 0.0 );
  else
    gsl_odeiv_control_init( B_.c_, 1e-6, 0.0, 1.0, 0.0 );

  if ( B_.e_ == 0 )
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_evolve_reset( B_.e_ );

  B_.sys_.function = compte2003_ex_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  for ( size_t n = 0; n < NCOMP; ++n )
    B_.I_stim_[ n ] = 0.0;
}

void nest::compte2003_ex::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate
  V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref >= 0, this can only fail in error
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::compte2003_ex::update( Time const& origin, const long_t from, const long_t to )
{

  assert( to >= 0 && ( delay ) from < Scheduler::get_min_delay() );
  assert( from < to );

  for ( long_t lag = from; lag < to; ++lag )
  {

    double t = 0.0;
    //previous step value of membrane potential. Used to detect negative derivative V_m - V_old <0
    const double V_old = S_.y_[State_::V_M]; 

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
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

    // add incoming spikes at end of interval
    // exploit here that spike buffers are compartment for compartment,
    // alternating between excitatory and inhibitory
    double  G_syn;
    for ( size_t n = 0 ; n < NCOMP ; ++n )
    {
      //excitatory inputs
      S_.y_[n*State_::STATE_VEC_COMPS + State_::G_AMPA] 
	+= B_.spikes_[4*n].get_value(lag);
	  G_syn = B_.spikes_[4*n+1].get_value(lag);
      S_.y_[n*State_::STATE_VEC_COMPS + State_::G_NMDA_FAST] 
	+= G_syn;
      S_.y_[n*State_::STATE_VEC_COMPS + State_::G_NMDA_SLOW] 
	+= G_syn;
      
      //inhibitory inputs
      S_.y_[n*State_::STATE_VEC_COMPS + State_::G_GABA] 
	+= B_.spikes_[4*n+3].get_value(lag);

    }

    // refractoriness and spiking
    // exploit here that plain offset enum value V_M indexes soma V_M
    if ( S_.r_ )
    { // neuron is absolute refractory and cannot spike
      --S_.r_;
    }
    else if  ( S_.y_[State_::V_M] >= 0 && V_old > S_.y_[State_::V_M])
      //spike is emitted if Vm > 0 and starts to decay

	  {
	    S_.r_ = V_.RefractoryCounts_;
	  
	    set_spiketime(Time::step(origin.get_steps()+lag+1));
	  
	    SpikeEvent se;
	    network()->send(*this, se, lag);
	    }

    // set new input currents
    for ( size_t n = 0; n < NCOMP; ++n )
      B_.I_stim_[ n ] = B_.currents_[ n ].get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::compte2003_ex::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );
  assert( 0 <= e.get_rport() && e.get_rport() < 5 * NCOMP );

  B_.spikes_[ e.get_rport() ].add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
}

void
nest::compte2003_ex::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );
  assert( 0 <= e.get_rport() && e.get_rport() < NCOMP ); // not 100% clean, should look at MIN, SUP

  // add weighted current; HEP 2002-10-04
  B_.currents_[ e.get_rport() ].add_value(
    e.get_rel_delivery_steps( network()->get_slice_origin() ), e.get_weight() * e.get_current() );
}

void
nest::compte2003_ex::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL
