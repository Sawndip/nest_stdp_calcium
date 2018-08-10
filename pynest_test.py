#TEST_STDPCalciumSynapse

import pylab
import nest
import numpy as np

nest.Install("stdp_calcium_module")

# First let us set the default parameters of model "iaf_psc_alpha" as we see fit
ndict = {"I_e": 200.0, "tau_m": 20.0}
nest.SetDefaults("iaf_psc_alpha", ndict)

# We can clone the model "iaf_psc_alpha" and modify the default parameters of the clone model with
edict = {"I_e": 200.0, "tau_m": 20.0}
nest.CopyModel("iaf_psc_alpha", "exc_iaf_neuron", params=edict)

idict = {"I_e": 300.0}
nest.CopyModel("iaf_psc_alpha", "inh_iaf_neuron", params=idict)

# Now creating populations of excitatory and inhibitory neurons with their parameters previously set
epop1 = nest.Create("exc_iaf_neuron", 100)
epop2 = nest.Create("exc_iaf_neuron", 100)
ipop1 = nest.Create("inh_iaf_neuron", 30)
ipop2 = nest.Create("inh_iaf_neuron", 30)

# Parameterization before initialization is not possible when random parameters are involved
# One can then loop on the neurons in the population:
Vth=-55.
Vrest=-70.

Vms = Vrest+(Vth-Vrest)*np.random.rand(len(epop1))
nest.SetStatus(epop1, "V_m", Vms)
# Note that this method of randomization only works for randomizing a single parameter
# Other more general (but less efficient) methods are available in tutorial 2

conn_dict = {"rule": "fixed_indegree", "indegree": 15}
syn_dict = {"model": "stdp_calcium_synapse"}
nest.Connect(epop1, epop2, conn_dict, syn_dict)

nest.GetConnections(epop1)