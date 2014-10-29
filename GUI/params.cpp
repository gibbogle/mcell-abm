#include <qstring.h>
#include "params.h"

Params::Params()
{
    PARAM_SET params[] = {

{"NCIRC", 15, 0, 50,
"Ncirc",
"Number of circumferential cells"},

{"NLONG", 60, 0, 50,
"Nlong",
"Number of longitudinal cells"},

{"NHOURS", 20.0, 0.0, 0.0,
"Number of hours",
"Length of the simulation.(hours)"},

{"DELTA_T", 3, 0, 0,
"Time step",
"Time step (mins)."},

{"DELTA_X", 1, 0, 0,
"Cell size",
"Cell size."},

{"PRESSURE", 0.0, 0, 0,
"Pressure",
"Internal pressure (cardiac jelly)"},

{"TENSION", 0.0, 0, 0,
"Tension",
"Tension at free vertices"},

{"FALPHA_AXIAL", 1, 0, 0,
"Axial force constant",
"Axial force constant"},

{"FALPHA_SHEAR", 0.01, 0, 0,
"Shear force constant",
"Shear force constant"},

{"FALPHA_BEND", 0.001, 0, 0,
"Bending force constant",
"Bening force constant"},

{"SOLVER", 0, 0, 0,
"Solver #",
"Solver number"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 1, 0, 0,
"NT_VTK",
"Animation interval."},

{"CELLML_FILE", 0, 0, 0,
"",
"CellML file for the cell growth model"},

// Entries after this point are QMyLabel dummies, to enable display of explanatory info  - no input data is transmitted


};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}
