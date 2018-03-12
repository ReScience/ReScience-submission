"""
Definition of the model with the neuronal simulator ANNarchy.
Original model by Clopath et al. 2010.
The implementation based on the matlab code, available on modeldb
(https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=144566).
Note: The homeostatic mechanism are fixed, as it is needed for the
plasticity experiments
"""

from ANNarchy import *
## Neuron Model for V1-Layer, after Clopath et al.(2010) ##
# neuron parameters
params = """
gL = 30.0
DeltaT = 2.0
tauw = 144.0
a = 4.0
b = 0.0805
EL = -70.6
C = 281.0
tauz = 40.0
tauVT= 50.0
Isp = 400.0
VTMax = -30.4
VTrest = -50.4
taux = 15.0
tauLTD = 10.0
tauLTP= 7.0
taumean = 2400.0
tau_gExc = 1.0
"""

# equations for the adaptive exponential integrate-and-fire neuron
# note: the changes in the membrane potential (vm) are based on the
# modeldb sourcecode:
# first step after a spike: vm is set 29.0+3.462 = 32.463 mV
# second step after spike: vm is set to EL (vm - (32.463+70.6)) and then increased
# with an value of [] (see matlab code for more information) + hyperpolarisation depending variables and inputs
# after second step: normal behaviour
neuron_eqs = """
dvm/dt = if state>=2:+3.462 else: if state==1:-(vm+51.75)+1/C*(Isp - (wad+b))+g_Exc-g_Inh else:1/C * ( -gL * (vm - EL) + gL * DeltaT * exp((vm - VT) / DeltaT) - wad + z ) + g_Exc: init = -70.6
dvmean/dt = (pos(vm - EL)**2 - vmean)/taumean    :init = 0.0
dumeanLTD/dt = (vm - umeanLTD)/tauLTD : init=-70.0
dumeanLTP/dt = (vm - umeanLTP)/tauLTP : init =-70.0
dxtrace /dt = (- xtrace )/taux
dwad/dt = if state ==2:0 else:if state==1:+b/tauw else: (a * (vm - EL) - wad)/tauw : init = 0.0
dz/dt = if state==1:-z+Isp-10 else:-z/tauz  : init = 0.0
dVT/dt =if state==1: +(VTMax - VT)-0.4 else:(VTrest - VT)/tauVT  : init=-50.4
dg_Exc/dt = -g_Exc/tau_gExc
state = if state > 0: state-1 else:0
Spike = 0.0
dresetvar / dt = 1/(1.0) * (-resetvar)
           """
spkNeurV1 = Neuron(parameters = params,equations=neuron_eqs,spike="""(vm>VT) and (state==0)""",
                         reset="""vm = 29.0
                                  state = 2.0
                                  VT = VTMax
                                  Spike = 1.0
                                  resetvar = 1.0
                                  xtrace+= 1/taux""")

#----------------------------------synapse definitions----------------------

#----- Synapse from Poisson to Input-Layer -----#
inputSynapse =  Synapse(
    parameters = "",
    equations = "",
    pre_spike = """
        g_target += w
                """
)

#--STDP Synapses after Clopath et. al(2008)for Input- to Exitatory- Layer--#
equatSTDP = """
    ltdTerm = if w>wMin : (aLTD*(vmean/urefsquare)*pre.Spike * pos(post.umeanLTD - thetaLTD)) else : 0.0
    ltpTerm = if w<wMax : (aLTP * pos(post.vm - thetaLTP) *(pre.xtrace)* pos(post.umeanLTP - thetaLTD)) else : 0.0
      deltaW = ( -ltdTerm + ltpTerm)
        dw/dt = deltaW :min=0.0,explicite"""

# 100/urefsquare for rateCode
parameterFF="""
vmean = 60.0
urefsquare = 60.0
thetaLTD = -70.6
thetaLTP = -45.3
aLTD = 0.00014
aLTP = 0.00008
wMin = 0.0
wMax =3.0
transmit = 0.0
"""
# notice the additional transmit variable (wich shoudl be 1 or 0)
# to transmit or not an EPSP to the postsynaptic neuron

#post.vmean
#*(post.vmean/urefsquare)
ffSyn = Synapse( parameters = parameterFF,
    equations= equatSTDP,
    pre_spike='''g_target += w*transmit''')
