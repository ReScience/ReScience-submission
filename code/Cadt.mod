TITLE decay of internal calcium concentration
:
: Simple extrusion mechanism for internal calium dynamics
: Model by Tabak et al. 2011
: Written by Geir Halnes, Norwegian Life Science University of Life Sciences, Aug 2018


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Cadt
	USEION Ca READ iCa, Cai WRITE Cai VALENCE 2
	RANGE fc,kc,alpha
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}


PARAMETER {
    fc = 0.01 (1): Fraction of Ca++ in cytoplasm
    alpha = 0.00945 (mM*cm^2/muC) : Converted from Tabak 2011, Membrane area-dependent
    kc = 0.12   (1/ms)  : From Tabak 2011
    Cainit  = 1e-4 (mM)	: Initial Ca-level
}


STATE {
	Cai		(mM) <1e-8> : to have tolerance of .01nM
}


INITIAL {
	Cai = Cainit
}


ASSIGNED {
	iCa		(mA/cm2)
}

	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 
    Cai' = - fc*(alpha*iCa + kc*Cai)
}
