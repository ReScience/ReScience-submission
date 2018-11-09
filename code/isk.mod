TITLE Medium duration Ca-dependent potassium current
:
:   Ca++ dependent K+ current IC responsible for medium duration AHP
:
:   SK-model from Tabak et al. 2011
:   Written by Geir Halnes 2018

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX sk
    USEION k READ ek WRITE ik VALENCE 1
    USEION Ca READ Cai VALENCE 2
    RANGE gskbar
    GLOBAL sinf
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v (mV)
	ek = -75 (mV)
	Cai 	= 1e-4 (mM)			: Initial [Ca]i = 50 nM (Cai is simulated by separate mod-file)
	gskbar	= 3.2e-4	(mho/cm2)	: Conductance (can be reset in hoc-file)
	ks	= 4e-4(mM)		: Middle point of sinf fcn
}


STATE {
	s
}


ASSIGNED {
    ik 	(mA/cm2)
    g
    sinf
}


BREAKPOINT { 
    sinf = Cai^2/(Cai^2+ks^2)
    ik = gskbar*sinf*(v - ek)
}


UNITSOFF
INITIAL {
:  activation kinetics are assumed to be at 22 deg. C
:  Q10 is assumed to be 3

	VERBATIM
	Cai = _ion_Cai;
	ENDVERBATIM

    sinf = Cai^2/(Cai^2+ks^2)
    
}

UNITSON
