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

#ifndef compte2003_in_H
#define compte2003_in_H

#include "config.h"

#ifdef HAVE_GSL

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
#include <gsl/gsl_odeiv.h>

namespace nest
{
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 */
extern "C" int compte2003_in_dynamics( double, const double*, double*, void* );

/**
 * @note All parameters that occur for both compartments
 *       and dendrite are stored as C arrays, with index 0 being soma.
 */
class compte2003_in : public Archiving_Node
{

  // Boilerplate function declarations --------------------------------

public:
  compte2003_in();
  compte2003_in( const compte2003_in& );
  ~compte2003_in();

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
  void update( Time const&, const long_t, const long_t );

  // Enumerations and constants specifying structure and properties ----

  //! Compartments, NCOMP is number
        enum Compartments_ {comp_0 = 0, comp_1,NCOMP };

  /**
   * Minimal spike receptor type.
   * @note Start with 1 so we can forbid port 0 to avoid accidental
   *       creation of connections with no receptor type set.
   */
  static const port MIN_SPIKE_RECEPTOR = 1;

  /**
   * Spike receptors.
   */
        enum SpikeSynapseTypes { C0_AMPA=MIN_SPIKE_RECEPTOR,C0_NMDA_FAST,C0_NMDA_SLOW,C0_GABA,
				C1_AMPA,C1_NMDA_FAST,C1_NMDA_SLOW,C1_GABA,
			       SUP_SPIKE_RECEPTOR };

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
        enum CurrentSynapseTypes { C0_I=MIN_CURR_RECEPTOR,C1_I, 
   			       SUP_CURR_RECEPTOR };

  static const size_t NUM_CURR_RECEPTORS = SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR;

  // Friends --------------------------------------------------------

  friend int compte2003_in_dynamics( double, const double*, double*, void* );

  friend class RecordablesMap< compte2003_in >;
  friend class UniversalDataLogger< compte2003_in >;


  // Parameters ------------------------------------------------------

  /**
   * Independent parameters of the model.
   * These parameters must be passed to the iteration function that
   * is passed to the GSL ODE solvers. Since the iteration function
   * is a C++ function with C linkage, the parameters can be stored
   * in a C++ struct with member functions, as long as we just pass
   * it by void* from C++ to C++ function. The struct must be public,
   * though, since the iteration function is a function with C-linkage,
   * whence it cannot be a member function of compte2003_in.
   * @note One could achieve proper encapsulation by an extra level
   *       of indirection: Define the iteration function as a member
   *       function, plus an additional wrapper function with C linkage.
   *       Then pass a struct containing a pointer to the node and a
   *       pointer-to-member-function to the iteration function as void*
   *       to the wrapper function. The wrapper function can then invoke
   *       the iteration function on the node (Stroustrup, p 418). But
   *       this appears to involved, and the extra indirections cost.
   */
    struct Parameters_ {

      double_t t_ref;      //!< Refractory period in ms
      double_t E_L;        //!< Leak reversal Potential (aka resting potential) in mV
      double_t E_Na;       //!<   
      double_t E_K;        //!<
      double_t E_Ca;       //!<         
     
      double_t g_conn[NCOMP-1];    //!< Conductances connecting compartments, in nS
      double_t g_L[NCOMP];         //!< Leak Conductance in nS
      double_t g_Na[NCOMP];        //!< 
      double_t g_K[NCOMP];         //!< 
      double_t g_Na_p[NCOMP];      //!< 
      double_t g_K_s[NCOMP];       //!< 
      double_t g_K_A[NCOMP];       //!< 
      double_t g_K_AR[NCOMP];       //!< 
      double_t g_K_Ca[NCOMP];       //!< 
      double_t g_K_Na[NCOMP];       //!< 
      double_t g_Ca[NCOMP];       //!< 
      double_t C_m[NCOMP];         //!< Membrane Capacitance in pF

      double_t tau_syn_ampa;    //!< Synaptic Time Constant Excitatory Synapse in ms
      double_t tau_syn_nmda_fast;         //
      double_t tau_syn_nmda_slow;         //
      double_t E_ex;        //!< Excitatory reversal Potential in mV
      double_t E_in;        //!< Inhibitory reversal Potential in mV
      double_t tau_syn_gaba;    //!< Synaptic Time Constant for Inhibitory Synapse in ms

      double_t I_e[NCOMP];         //!< Constant Current in pA

    Parameters_();                                //!< Sets default parameter values
    Parameters_( const Parameters_& );            //!< needed to copy C-arrays
    Parameters_& operator=( const Parameters_& ); //!< needed to copy C-arrays

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };


  // State variables  ------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor and assignment operator required because
   *       of C-style array.
   */
public:
  struct State_
  {

    /**
     * Elements of state vector.
     * For the multicompartmental case here, these are offset values.
     * The state variables are stored in contiguous blocks for each
     * compartment, beginning with the soma.
     */
      enum StateVecElems_ { V_M = 0, G_AMPA, G_NMDA_FAST,G_NMDA_SLOW,
			    G_GABA,n_Ca,n_Na,Na_h, K_m, K_s_m, K_A_h,STATE_VEC_COMPS };

    //! total size of state vector
    static const size_t STATE_VEC_SIZE = STATE_VEC_COMPS * NCOMP;

    //! neuron state, must be C-array for GSL solver
    double_t y_[ STATE_VEC_SIZE ];
    int_t r_; //!< number of refractory steps remaining

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );
    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_& );

    /**
     * Compute linear index into state array from compartment and element.
     * @param comp compartment index
     * @param elem elemet index
     * @note compartment argument is not of type Compartments_, since looping
     *       over enumerations does not work.
     */
    static size_t
    idx( size_t comp, StateVecElems_ elem )
    {
      return comp * STATE_VEC_COMPS + elem;
    }
  };

private:
  // Internal buffers --------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( compte2003_in& );                  //!<Sets buffer pointers to 0
    Buffers_( const Buffers_&, compte2003_in& ); //!<Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< compte2003_in > logger_;

    /** buffers and sums up incoming spikes/currents
     *  @note Using STL vectors here to ensure initialization.
     */
    std::vector< RingBuffer > spikes_;
    std::vector< RingBuffer > currents_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double_t step_;          //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    /**
     * Input currents injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double_t I_stim_[ NCOMP ]; //!< External Stimulus in pA
  };

  // Internal variables ---------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    /** initial value to normalise excitatory synaptic conductance */
    double_t PSConInit_E_[ NCOMP ];

    /** initial value to normalise inhibitory synaptic conductance */
    double_t PSConInit_I_[ NCOMP ];

    int_t RefractoryCounts_;
  };

  // Access functions for UniversalDataLogger -------------------------------

  /**
   * Read out state vector elements, used by UniversalDataLogger
   * First template argument is component "name", second compartment "name".
   */
  template < State_::StateVecElems_ elem, Compartments_ comp >
  double_t
  get_y_elem_() const
  {
    return S_.y_[ S_.idx( comp, elem ) ];
  }

  //! Read out number of refractory steps, used by UniversalDataLogger
  double_t
  get_r_() const
  {
    return Time::get_resolution().get_ms() * S_.r_;
  }

  // Data members ----------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Table of compartment names
  static std::vector< Name > comp_names_;

  //! Dictionary of receptor types, leads to seg fault on exit, see #328
  // static DictionaryDatum receptor_dict_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< compte2003_in > recordablesMap_;
};

inline port
compte2003_in::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
compte2003_in::handles_test_event( SpikeEvent&, rport receptor_type )
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
compte2003_in::handles_test_event( CurrentEvent&, rport receptor_type )
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
compte2003_in::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
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
compte2003_in::get_status( DictionaryDatum& d ) const
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
  (*receptor_dict_)[Name("C0_ampa")]  = C0_AMPA;
  (*receptor_dict_)[Name("C0_nmda_fast")]  = C0_NMDA_FAST;
  (*receptor_dict_)[Name("C0_nmda_slow")]  = C0_NMDA_SLOW;
  (*receptor_dict_)[Name("C0_gaba")]  = C0_GABA;
  (*receptor_dict_)[Name("C0_curr")]  = C0_I;
  //Compartment 1
  (*receptor_dict_)[Name("C1_ampa")]  = C1_AMPA;
  (*receptor_dict_)[Name("C1_nmda_fast")]  = C1_NMDA_FAST;
  (*receptor_dict_)[Name("C1_nmda_slow")]  = C1_NMDA_SLOW;
  (*receptor_dict_)[Name("C1_gaba")]  = C1_GABA;
  (*receptor_dict_)[Name("C1_curr")]  = C1_I;

  ( *d )[ names::receptor_types ] = receptor_dict_;
}

inline void
compte2003_in::set_status( const DictionaryDatum& d )
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


#endif // HAVE_GSL
#endif // compte2003_in_H
