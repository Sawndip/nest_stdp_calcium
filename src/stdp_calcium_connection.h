/*
 *  stdp_connection_calcium.h
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

#ifndef STDP_CONNECTION_CALCIUM
#define STDP_CONNECTION_CALCIUM

/* BeginDocumentation
  Name: stdp_connection_calcium - Synapse type for spike-timing dependent
   plasticity, based on calcium-history, with stability transfer function
   (individual synapse is assumed to be bistable)

  Description:
   stdp_connection_calcium is a connector to create synapses with a spike time
   dependent plasticity based on a calcium history.
   The model is a modified version of [1].
   It integrates a phenomenological dynamics for a phosphorylation variable, on which the synaptic weight directly depends
   (see Graupner & Brunel 2007, http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030221)
   
   The model relies on the following equations:
   tau_ca * d[ca]/dt = -ca + sum[C_pre*dirac[t_pre] + c_post*dirac[t_post]]
   tau_rho * d[rho]/dt = gamma_p*(rho_max - rho)*(ca > theta_p) - gamma_d*rho*(ca > theta_d)
   tau_w * d[weight]/dt = transfer(rho) - weight

  Variables:
   weight    double    The synaptic weight, between 0 and 1
   rho       double    A phosphorylation proxy variable that directly determines the weight
   ca        double    The current calcium level at the synapse

  Parameters:
   tau_ca    double    The time constant of calcium decay at the synapse
   C_pre     double    The amount of calcium added by a presynaptic spike
   C_post    double    The amount of calcium added by a postsynaptic spike
   theta_p   double    The threshold that the calcium level must exceed for potentiation to take place
   theta_d   double    The threshold that the calcium level must exceed for depression to take place
   gamma_p   double    Controls the speed of potentiation of the phosphorylation variable
   gamma_d   double    Controls the speed of depression of the phosphorylation variable
   tau_rho   double    Time constant for phosphorylation
   S_attr    double    Phosphorylation level for the saddle node bifurcation that determines bistability of the phosphorylation var at basal calcium level
   sigma     double    Noise level
   rho_max   double    Maximum phosphorylation level
   p_active_ bool      Active potentiation flag (ca > theta_p)?
   d_active_ bool      Active depression flag (ca > theta_d)?

  Transmits: SpikeEvent

  References:
   [1]  Graupner M, Brunel N (2012)
        Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location.
        Proc Natl Acad Sci U S A. 2012 Mar 6; 109(10): 3991–3996.
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3309784/

   [2]  Graupner M, Brunel N (2007)
        STDP in a Bistable Synapse Model Based on CaMKII and Associated Signaling Pathways.
        PLoS Comput Biol 3(11): e221.
        https://doi.org/10.1371/journal.pcbi.0030221

  FirstVersion: August 2018
  Author: Simon Lebastard
*/

// C++ includes:
#include <cmath>
// Including boost for inverse erfc function
#include <boost/math/special_functions/erf.hpp>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

// Include stdpcalcium namespace
#include "stdp_calcium_names.h"


namespace stdpcalcium
{

  class STDPCalciumCommonProperties : public nest::CommonSynapseProperties
  {
  
      template < typename targetidentifierT >
      friend class STDPCalciumConnection;
    
    public:
      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      STDPCalciumCommonProperties();
    
      /**
       * Get all properties and put them into a dictionary.
       */
      void get_status( DictionaryDatum& d ) const;
    
      /**
       * Set properties from the values given in dictionary.
       */
      void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );
    
      // data members common to all connections
      double tau_ca_;
      double C_pre_;
      double C_post_;
      double theta_p_;
      double theta_d_;
      double gamma_p_;
      double gamma_d_;
      double tau_rho_;
      double S_attr_;
      double sigma_;
      double rho_max_;
  };



  template < typename targetidentifierT >
  class STDPCalciumConnection : public nest::Connection< targetidentifierT >
  {
    public:

      typedef STDPCalciumCommonProperties CommonPropertiesType;
      typedef nest::Connection< targetidentifierT > ConnectionBase;

      STDPCalciumConnection ();
      STDPCalciumConnection ( const STDPCalciumConnection& );
      ~STDPCalciumConnection (){}

      using ConnectionBase::get_delay_steps;
      using ConnectionBase::get_delay;
      using ConnectionBase::get_rport;
      using ConnectionBase::get_target;

      /**
       * Get all properties of this connection and put them into a dictionary.
       */
      void get_status( DictionaryDatum& d ) const;

      /**
       * Set properties of this connection from the values given in dictionary.
       */
      void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

      /**
       * Send an event to the receiver of this connection.
       * \param e The event to send
       * \param cp common properties of all synapses (empty).
       */
  	  void send( nest::Event& e,
    	nest::thread t,
    	double t_lastspike,
    	const STDPCalciumCommonProperties& cp );

      class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
      {
      public:
        // Ensure proper overriding of overloaded virtual functions.
        // Return values from functions are ignored.
        using nest::ConnTestDummyNodeBase::handles_test_event;
        nest::port
        handles_test_event( nest::SpikeEvent&, nest::rport )
        {
          return nest::invalid_port_;
        }
      };


      void
      check_connection( nest::Node& s,
        nest::Node& t,
        nest::rport receptor_type,
        double t_lastspike,
        const CommonPropertiesType& )
      {
        ConnTestDummyNode dummy_target;
    
        ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
    
        t.register_stdp_connection( t_lastspike - get_delay() );
      }

     /*
	  * set_weight is required by nest infrastructure, but is not required here
	  * the phosphorylation variable will NOT be set to the corresponding value,
	  * and the calcium history will remain untouched
	  * In particular, setting weight to zero may not guarantee that the weight directly
	  * remains at zero since there could be a high local calcium concentration,
	  * triggering instantaneous plasticity
      */
      void set_weight( double w )
      {
          weight_ = w;
      }


    private:
      double weight_;
      double rho_;
      double ca_;

      double t_lastspike;
      bool p_active_;
      bool d_active_;

      double transfer_( const STDPCalciumCommonProperties cp, double rho );

      double invtransfer_( const STDPCalciumCommonProperties cp, double w );

      void calcium_update_prepost_( const STDPCalciumCommonProperties cp, double dt );

      void calcium_update_postpre_( const STDPCalciumCommonProperties cp, double dt );

      void update_rho_( const STDPCalciumCommonProperties cp, double t_close_p, double t_close_d );
  };

  template < typename targetidentifierT >
  double STDPCalciumConnection< targetidentifierT >::transfer_( const STDPCalciumCommonProperties cp, double rho )
  {
    return erfc((cp.S_attr_ - rho)/cp.sigma_);
  }

  template < typename targetidentifierT >
  double STDPCalciumConnection< targetidentifierT >::invtransfer_( const STDPCalciumCommonProperties cp, double w )
  {
  	using namespace boost::math;
    return std::min(std::max(0.0, cp.S_attr_ - 2*std::pow(cp.sigma_,2.0)*erfc_inv(w)), cp.rho_max_);
  }
  
  template < typename targetidentifierT >
  void STDPCalciumConnection< targetidentifierT >::calcium_update_prepost_( const STDPCalciumCommonProperties cp, double dt )
  {
    ca_ = ca_*std::exp(-dt/cp.tau_ca_) + cp.C_post_;
  }
  
  template < typename targetidentifierT >
  void STDPCalciumConnection< targetidentifierT >::calcium_update_postpre_( const STDPCalciumCommonProperties cp, double dt )
  {
    ca_ = ca_*std::exp(-dt/cp.tau_ca_) + cp.C_pre_;
  }
  
  template < typename targetidentifierT >
  void STDPCalciumConnection< targetidentifierT >::update_rho_( const STDPCalciumCommonProperties cp, double t_close_p, double t_close_d )
  {
    double t_pure_d = t_close_d - t_close_p;
    assert( t_pure_d >= 0 );
  
    rho_ = ((rho_ - cp.rho_max_*cp.gamma_p_/(cp.gamma_p_+cp.gamma_d_))*std::exp(-(t_close_p)*(cp.gamma_p_+cp.gamma_d_)/cp.tau_rho_) + cp.rho_max_*cp.gamma_p_/(cp.gamma_p_+cp.gamma_d_));
    rho_ = rho_* std::exp(-cp.gamma_d_*(t_pure_d)/cp.tau_rho_);
  
    weight_ = transfer_( cp, rho_ );
  }

  template < typename targetidentifierT >
  void STDPCalciumConnection< targetidentifierT >::send( nest::Event& e,
    nest::thread t,
    double t_lastspike,
    const STDPCalciumCommonProperties& cp )
  {
    // synapse STDP depressing/facilitation dynamics
    double t_spike = e.get_stamp().get_ms();
  
    // use accessor functions (inherited from Connection< >) to obtain delay and
    // target
    nest::Node* target = get_target( t );
    double dendritic_delay = get_delay();
  
    // get spike history in relevant range (t1, t2] from post-synaptic neuron
    std::deque< nest::histentry >::iterator start;
    std::deque< nest::histentry >::iterator finish;
  
    // For a new synapse, t_lastspike contains the point in time of the last
    // spike. So we initially read the
    // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
    // which increases the access counter for these entries.
    // At registration, all entries' access counters of
    // history[0, ..., t_last_spike - dendritic_delay] have been
    // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
    // details.
    target->get_history( t_lastspike - dendritic_delay,
      t_spike - dendritic_delay,
      &start,
      &finish );
    // facilitation due to post-synaptic spikes since last pre-synaptic spike
    double post_spike;
    double post_previousspike;
    double post_lastspike;
    double ca_delay;

    double t_max_p_duration;
    double t_close_p;
    double t_max_d_duration;
    double t_close_d;

    post_lastspike = finish->t_;
    post_previousspike = t_lastspike;
  
    while ( start != finish )
    {
      post_spike = start->t_;
      ++start;
      ca_delay = post_spike - post_previousspike;
      post_previousspike = post_spike;
  
      t_max_p_duration = p_active_ * cp.tau_ca_*std::log(ca_/cp.theta_p_);
      t_close_p = std::min(p_active_ * (post_spike - post_previousspike), t_max_p_duration);
  
      t_max_d_duration = d_active_ * cp.tau_ca_*std::log(ca_/cp.theta_d_);
      t_close_d = std::min(d_active_ * (post_spike - post_previousspike), t_max_d_duration);
  
      calcium_update_prepost_( cp, ca_delay );
      p_active_ = (ca_ >= cp.theta_p_);
      d_active_ = (ca_ >= cp.theta_d_);
  
      update_rho_( cp, t_close_p, t_close_d );
    }
  
    t_max_p_duration = p_active_ * cp.tau_ca_*std::log(ca_/cp.theta_p_);
    t_close_p = std::min(p_active_ * (t_spike - post_previousspike), t_max_p_duration);
  
    t_max_d_duration = d_active_ * cp.tau_ca_*std::log(ca_/cp.theta_d_);
    t_close_d = std::min(d_active_ * (t_spike - post_previousspike), t_max_d_duration);
  
    calcium_update_postpre_( cp, t_spike - post_lastspike );
    p_active_ = (ca_ >= cp.theta_p_);
    d_active_ = (ca_ >= cp.theta_d_);
    update_rho_( cp, t_close_p, t_close_d );
  
  
    e.set_receiver( *target );
    e.set_weight( weight_ );
    // use accessor functions (inherited from Connection< >) to obtain delay in
    // steps and rport
    e.set_delay( get_delay_steps() );
    e.set_rport( get_rport() );
    e();
  }
  
  
  
  template < typename targetidentifierT >
  STDPCalciumConnection< targetidentifierT >::STDPCalciumConnection()
    : ConnectionBase()
    , weight_(0.0)
    , rho_(0.0)
    , ca_(0.0)
    , p_active_( false )
    , d_active_( false )
  {
  }
  
  template < typename targetidentifierT >
  STDPCalciumConnection< targetidentifierT >::STDPCalciumConnection(
    const STDPCalciumConnection< targetidentifierT >& rhs )
    : ConnectionBase( rhs )
    , weight_( rhs.weight_ )
    , rho_( rhs.rho_ )
    , ca_( rhs.ca_ )
    , p_active_( rhs.p_active_ )
    , d_active_( rhs.d_active_ )
  {
  }
  
  template < typename targetidentifierT >
  void
  STDPCalciumConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
  {
    ConnectionBase::get_status( d );
    def< double >( d, names::weight, weight_ );
    def< double >( d, names::rho, rho_ );
    def< double >( d, names::ca, ca_ );
    def< bool >( d, names::p_active, p_active_ ); 
    def< bool >( d, names::d_active, d_active_ ); 
  }
  
  template < typename targetidentifierT >
  void
  STDPCalciumConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
    nest::ConnectorModel& cm )
  {
    ConnectionBase::set_status( d, cm );
    updateValue< double >( d, names::weight, weight_ );
    updateValue< double >( d, names::rho, rho_ );
    updateValue< double >( d, names::ca, ca_ );
    updateValue< bool >( d, names::p_active, p_active_ ); 
    updateValue< bool >( d, names::d_active, d_active_ ); 
  }

} // namespace nest

#endif /* #ifndef STDP_CONNECTION_CALCIUM */