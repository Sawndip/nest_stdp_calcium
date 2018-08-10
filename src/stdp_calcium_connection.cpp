/*
 *  stdp_connection_calcium.cpp
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

#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_calcium_connection.h"
#include "event.h"

namespace stdpcalcium
{
//
// Implementation of class STDPCalciumCommonProperties.
//
  
  STDPCalciumCommonProperties::STDPCalciumCommonProperties()
    : CommonSynapseProperties()
    , tau_ca_(80.0)
    , C_pre_(0.40)
    , C_post_(0.84)
    , theta_p_(1.08)
    , theta_d_(1.00)
    , gamma_p_(120.0)
    , gamma_d_(200.0)
    , tau_rho_(100000.0)
    , S_attr_(40.0)
    , sigma_(25.0)
    , rho_max_(200.0)
  {
  }
  
  void
  STDPCalciumCommonProperties::get_status( DictionaryDatum& d ) const
  {
    CommonSynapseProperties::get_status( d );
  
    def< double >( d, names::tau_ca, tau_ca_ );
    def< double >( d, names::C_pre, C_pre_ );
    def< double >( d, names::C_post, C_post_ );
    def< double >( d, names::theta_p, theta_p_ );
    def< double >( d, names::theta_d, theta_d_ );
    def< double >( d, names::gamma_p, gamma_p_ );
    def< double >( d, names::gamma_d, gamma_d_ );
    def< double >( d, names::tau_rho, tau_rho_ );
    def< double >( d, names::S_attr, S_attr_ );
    def< double >( d, names::sigma, sigma_ );
    def< double >( d, names::rho_max, rho_max_ );
  }
  
  void
  STDPCalciumCommonProperties::set_status( const DictionaryDatum& d,
    nest::ConnectorModel& cm )
  {
    CommonSynapseProperties::set_status( d, cm );
  
    updateValue< double >( d, names::tau_ca, tau_ca_ );
    updateValue< double >( d, names::C_pre, C_pre_ );
    updateValue< double >( d, names::C_post, C_post_ );
    updateValue< double >( d, names::theta_p, theta_p_ );
    updateValue< double >( d, names::theta_d, theta_d_ );
    updateValue< double >( d, names::gamma_p, gamma_p_ );
    updateValue< double >( d, names::gamma_d, gamma_d_ );
    updateValue< double >( d, names::tau_rho, tau_rho_ );
    updateValue< double >( d, names::S_attr, S_attr_ );
    updateValue< double >( d, names::sigma, sigma_ );
    updateValue< double >( d, names::rho_max, rho_max_ );
  }
  
} // of namespace nest
