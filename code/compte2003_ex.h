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
 * This model is based on the "aeif_cond_exp" model in NEST 
 *------------------------------------------------------------------------------
*/ 
#ifndef COMPTE2003_EX_H
#define COMPTE2003_EX_H

#include "config.h"

#ifdef HAVE_GSL_1_11

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include "dictdatum.h"
#include "name.h"
#include <vector>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


namespace nest
{
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int compte2003_ex_dynamics( double, const double*, double*, void* );

class compte2003_ex : public Archiving_Node
{

public:
  compte2003_ex();
  compte2003_ex( const compte2003_ex& );
  ~compte2003_ex();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool );

  void handle( SpikeEvent& );
  void handle( CurrentEvent& );
  void handle( DataLoggingRequest& );

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();
  void update( const Time&, const long_t, const long_t );

  // Enumerations and constants specifying structure and properties ----
  
  /**
   * Minimal spike receptor type.
   * @note Start with 1 so we can forbid port 0 to avoid accidental
   *       creation of connections with no receptor type set.
   */
  static const port MIN_SPIKE_RECEPTOR = 1;

  /**
   * Spike receptors.
   */
  enum SpikeSynapseTypes { 
       AMPA=MIN_SPIKE_RECEPTOR, NMDA_FAST, NMDA_SLOW, GABA, SUP_SPIKE_RECEPTOR };

  static const size_t NUM_SPIKE_RECEPTORS = SUP_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR;
  /**
   * Minimal current receptor type.
   *  @note Start with SUP_SPIKE_RECEPTOR to avoid any overlap and
   *        accidental mix-ups.
   */
  static const port MIN_CURR_RECEPTOR = SUP_SPIKE_RECEPTOR;

  /**
   * Current receptors.
   */
  enum CurrentSynapseTypes { I=MIN_CURR_RECEPTOR, SUP_CURR_RECEPTOR };

  static const size_t NUM_CURR_RECEPTORS = SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR;

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int compte2003_ex_dynamics( double, const double*, double*, void* );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< compte2003_ex >;
  friend class UniversalDataLogger< compte2003_ex >;

private:
  //! Independent parameters
    struct Parameters_ {

      double_t t_ref;      //!< Refractory period in ms
      double_t E_L;        //!< Leak reversal Potential (aka resting potential) in mV
      double_t E_Na;       //!<   
      double_t E_K;        //!<
      double_t E_Ca;       //!<         
     
      double_t g_conn;    //!< Conductances connecting compartments, in nS
      
      double_t tau_syn_ampa;    //!< Synaptic Time Constant Excitatory Synapse in ms
      double_t tau_syn_nmda_fast;         //
      double_t tau_syn_nmda_slow;         //
      double_t E_ex;        //!< Excitatory reversal Potential in mV
      double_t E_in;        //!< Inhibitory reversal Potential in mV
      double_t tau_syn_gaba; //!< Synaptic Time Constant for Inhibitory Synapse in ms

      //soma
      double_t C_m_s;         //!< Somatic membrane Capacitance in pF
      double_t g_L;         //!< Leak Conductance in nS
      double_t g_Na;        //!< 
      double_t g_K;         //!< 
      double_t g_K_s;       //!< 
      double_t g_K_A;       //!<       
      double_t g_K_Na;       //!< 
      double_t I_e;         //!< Constant Current to soma in pA
            
      //dendr            
      double_t C_m_d;         //!< Dendr membrane Capacitance in pF
      double_t g_Na_p;      //!< 
      double_t g_K_AR;       //!< 
      double_t g_K_Ca;       //!< 
      double_t g_Ca;       //!< 
      
    Parameters_();                                //!< Sets default parameter values
    //!!!Parameters_( const Parameters_& );            //!< needed to copy C-arrays
    //!!!Parameters_& operator=( const Parameters_& ); //!< needed to copy C-arrays

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };
  

public:
  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor and assignment operator required because
   *       of C-style array.
   */
  struct State_
  {
    /**
     * Enumeration identifying elements in state array State_::y_.
     * The state vector must be passed to GSL as a C array. This enum
     * identifies the elements of the vector. It must be public to be
     * accessible from the iteration function.
     */
    enum StateVecElems { 
         V_M_S = 0, 
         V_M_D, //1
         G_AMPA,//2 
         G_NMDA_FAST,//3 
         G_NMDA_SLOW,//4 
         G_GABA,//5
         N_CA,//6
         N_NA,//7
         NA_H,//8 
         K_M, //9
         K_S_M,//10 
         K_A_H,//11 
         STATE_VEC_SIZE };


    double_t y_[ STATE_VEC_SIZE ]; //!< neuron state, must be C-array for GSL solver
    int_t r_;                      //!< number of refractory steps remaining

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );
    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_& );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( compte2003_ex& );                  //!<Sets buffer pointers to 0
    Buffers_( const Buffers_&, compte2003_ex& ); //!<Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< compte2003_ex > logger_;

    /** buffers and sums up incoming spikes/currents
     *  @note Using STL vectors here to ensure initialization.
     */
    std::vector< RingBuffer > spikes_;
    std::vector< RingBuffer > currents_;

    /** GSL ODE stuff */
    gsl_odeiv2_step* s_;    //!< stepping function
    gsl_odeiv2_control* c_; //!< adaptive stepsize control function
    gsl_odeiv2_evolve* e_;  //!< evolution function
    gsl_odeiv2_system sys_; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double_t step_;          //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    /**
     * Input current injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double_t I_stim_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    int_t RefractoryCounts_;
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger
  template < State_::StateVecElems elem >
  double_t
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  // ----------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< compte2003_ex > recordablesMap_;
};

inline port
compte2003_ex::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
compte2003_ex::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type < MIN_SPIKE_RECEPTOR || receptor_type >= SUP_SPIKE_RECEPTOR )
  {
    if ( receptor_type < 0 || receptor_type >= SUP_CURR_RECEPTOR )
      throw UnknownReceptorType( receptor_type, get_name() );
    else
      throw IncompatibleReceptorType( receptor_type, get_name(), "SpikeEvent" );
  }
  return receptor_type - MIN_SPIKE_RECEPTOR;
}

inline port
compte2003_ex::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type < MIN_CURR_RECEPTOR || receptor_type >= SUP_CURR_RECEPTOR )
  {
    if ( receptor_type >= 0 && receptor_type < MIN_CURR_RECEPTOR )
      throw IncompatibleReceptorType( receptor_type, get_name(), "CurrentEvent" );
    else
      throw UnknownReceptorType( receptor_type, get_name() );
  }
  return receptor_type - MIN_CURR_RECEPTOR;
}

inline port
compte2003_ex::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    if ( receptor_type < 0 || receptor_type >= SUP_CURR_RECEPTOR )
      throw UnknownReceptorType( receptor_type, get_name() );
    else
      throw IncompatibleReceptorType( receptor_type, get_name(), "DataLoggingRequest" );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}


inline void
compte2003_ex::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();

  /**
   * @todo dictionary construction should be done only once for
   * static member in default c'tor, but this leads to
   * a seg fault on exit, see #328
   */
  DictionaryDatum receptor_dict_ = new Dictionary();
  (*receptor_dict_)[Name("ampa")]  = AMPA;
  (*receptor_dict_)[Name("nmda_fast")]  = NMDA_FAST;
  (*receptor_dict_)[Name("nmda_slow")]  = NMDA_SLOW;
  (*receptor_dict_)[Name("gaba")]  = GABA;
  (*receptor_dict_)[Name("curr")]  = I;


  ( *d )[ names::receptor_types ] = receptor_dict_;
}


inline void
compte2003_ex::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp );   // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif // HAVE_GSL_1_11
#endif // COMPTE2003_E_H
