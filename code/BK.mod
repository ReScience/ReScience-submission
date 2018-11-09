TITLE BK channel
: Model from Tabak et al. 2011 J Neurosci 31(46)
: Used as BK channel for gonandotrope cell in pituitary gland of medaka.
: Implemented by G.Halnes March 2018

NEURON {
    SUFFIX bk
    USEION k READ ek WRITE ik
    GLOBAL finf
    RANGE gbk, ftau
}

PARAMETER {
        gbk = 0.00 (mho/cm2)
        ek (mV)		: must be explicitely def. in hoc
        sf = 2   (mV)   : TODO changed this from 1
        vf = -20   (mV) : TODO changed this from 15
        v (mV)
        ftau = 5 (ms)
}


UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}


ASSIGNED {
        ik (mA/cm2)
        finf
}

STATE {f}

BREAKPOINT {
        SOLVE states METHOD cnexp
        ik = gbk*f*(v-ek)
}

INITIAL {
	trates(v)
}

DERIVATIVE states {
        trates(v)
        f' = (finf-f)/ftau
}

PROCEDURE trates(v (mV)) { :callable from hoc
        finf = 1/(1+exp((vf -v)/sf) )
}
