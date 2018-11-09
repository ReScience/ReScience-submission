TITLE HVA-channel channel
: Model from Tabak et al. 2011 J Neurosci 31(46)
: Implemented by G.Halnes August 2018
  
NEURON {
    SUFFIX ihvat
    USEION Ca WRITE iCa VALENCE 2
    RANGE ghvat
}
     
PARAMETER {        
        ghvat = 3.2e-4 (mho/cm2)
        eCa = 60 (mV)		: must be explicitely def. in hoc
        sm = 12   (mV)
        vm = -20   (mV)
        v (mV)
}

                                     
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}


ASSIGNED {
        iCa (mA/cm2)
        finf
}

BREAKPOINT {
        finf = 1/(1+exp((vm -v)/sm) )
        iCa = ghvat*finf*(v-eCa)
}

INITIAL {
	finf = 0.01
}


