TITLE K-DR channel
: From Tabak et al. 2011
: G.Halnes 2018

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
        v (mV)
        ek (mV)		: must be explicitely def. in hoc
        celsius		(degC)
        gkdrbar= 4.8e-4 (mho/cm2)
        vlj = 0 (mV) : Added for medaka-cell (shift threshold higher by 0 mV)
        sn = 10 (mV)
        vn = -5 (mV)
        taun = 30 (ms)
: Noter at jeg hastet denne opp til 10 ms tidligere
}


NEURON {
	SUFFIX kdrt
	USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar
	GLOBAL ninf,taun
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkdr
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkdrbar*n*(v-ek)

}

INITIAL {
	rates(v+vlj, vn, sn)
	n=0.1
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v+vlj,vn, sn)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v,vn,sn) { :callable from hoc
        ninf = 1/(1+exp((vn-v)/sn))
}




















