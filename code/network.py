"""
Definition of the model with the neuronal simulator ANNarchy.
Original model by Clopath et al. (2010).
The implementation based on the Matlab code, available on modeldb
(https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=144566).
Note: homeostatic mechanism are here working as it is needed for the
connectivity experiments
"""

from ANNarchy import *
""" Neuron Model for V1-Layer, after Clopath et al.(2010) """

# Neuron parameters
params = """
gL = 30.0       :population
DeltaT = 2.0    :population
tauw = 144.0    :population
a = 4.0         :population
b = 0.0805      :population
EL = -70.6      :population
C = 281.0       :population
tauz = 40.0     :population
tauVT= 50.0     :population
Isp = 400.0     :population
VTMax = -30.8   :population
VTrest = -50.4  :population
taux = 15.0     :population
tauLTD = 10.0   :population
tauLTP= 7.0     :population
taumean = 1200.0:population
tau_gExc = 1.0  :population
inter_vm = -49.5
"""
"""
 With the ':population' keyword,
 ANNarchy knows that the parameter set is equal for the complete population
 and do not create a set for each neuron to save memory.
"""

"""
 Equations for the adaptive exponential integrate-and-fire neuron
 note: the changes in the membrane potential (vm) are based on the
 modeldb source code:
 first step after a spike: vm is set 29.0+3.462 = 32.463 mV,
 second step after spike: vm is set to EL (vm - (32.463+70.6)) and then increased
 with hyperpolarisation depending variables and inputs.
 After second step: normal behavior."""

neuron_eqs = """
dvm/dt = if state>=2:+3.462 else: if state==1:-vm + inter_vm +1/C*( - wad+b + z)+g_Exc else:1/C * ( -gL * (vm - EL) + gL * DeltaT * exp((vm - VT) / DeltaT) - wad + z ) + g_Exc : init = -70.6
dvmean/dt = (pos(vm - EL)**2 - vmean)/taumean    :init = 0.0
dumeanLTD/dt = (vm - umeanLTD)/tauLTD : init=-70.0
dumeanLTP/dt = (vm - umeanLTP)/tauLTP : init =-70.0
dxtrace /dt = (- xtrace )/taux
dwad/dt = if state==1:+b else: (a * (vm - EL) - wad)/tauw : init = 0.0
dz/dt = if state==1:-z+Isp-10 else:-z/tauz  : init = 0.0
dVT/dt =if state==2: 0 else: if state==1: -VT +VTMax else:(VTrest - VT)/tauVT  : init=-50.4
dg_Exc/dt = -g_Exc/tau_gExc
state = if state > 0: state-1 else:0
Spike = 0.0
dresetvar / dt = 1/(1.0) * (-resetvar)
           """

"""
 Create the 'Neuron' object with the parameters and equations for the
 adaptive exponential integrate-and-fire (AdEx) neuron model.
 Define the spiking conditions (in 'spike') and what happens after a spike ('reset').
"""
AdExNeuron = Neuron(parameters = params,equations=neuron_eqs,spike="""(vm>VT) and (state==0)""",
                         reset="""vm = 29.4
                                  state = 2.0
                                  Spike = 1.0
                                  xtrace+= 1/taux""")

# Synapses29.4
inputSynapse =  Synapse(
    parameters = "",
    equations = "",
    pre_spike = """
        g_target += w
                """
)

""" STDP Synapses after Clopath et. al(2010)for Input- to Exitatory- Layer """
equatSTDP = """
    ltdTerm_fix = if w>wMin : (aLTD*(vmean_fix/urefsquare)*pre.Spike * pos(post.umeanLTD - thetaLTD)) else : 0.0
    ltdTerm = if w>wMin : (aLTD*(post.vmean/urefsquare)*pre.Spike * pos(post.umeanLTD - thetaLTD)) else : 0.0
    ltpTerm = if w<wMax : (aLTP * pos(post.vm - thetaLTP) *(pre.xtrace)* pos(post.umeanLTP - thetaLTD)) else : 0.0
      deltaW = if set_fix==1: ltpTerm - ltdTerm_fix else: ltpTerm - ltdTerm 
        dw/dt = deltaW :min=0.0, max=wMax,explicite"""
#deltaW = if set_fix==1: ( -ltdTerm_fix + ltpTerm) else: ( -ltdTerm + ltpTerm)
"""
 Definition of the parameters for the synapses,
 with the ':projection' keyword knows ANNarchy, that the parameter set is equal
 for all the synapses in the project object and do not create a set for every synapse
 to save memory.
"""

parameterFF="""
vmean_fix = 60.0    :projection
urefsquare = 60.0   :projection
thetaLTD = -70.6    :projection
thetaLTP = -45.3    :projection
aLTD = 0.00014      :projection
aLTP = 0.00008      :projection
wMin = 0.0          :projection
wMax =3.0           :projection
transmit = 0.0      :projection
set_fix = 0.0       :projection
"""
# Notice the additional transmit variable (which should be 1 or 0)
# to transmit or not an EPSP to the postsynaptic neuron

"""
 Create the synapse object with the parameters ('parameters'),
 the learning equation ('equations')
 and what happend, if a presynaptic spike occurs.
"""
ffSyn = Synapse( parameters = parameterFF,
    equations= equatSTDP,
    pre_spike='''g_target += w*transmit''')
