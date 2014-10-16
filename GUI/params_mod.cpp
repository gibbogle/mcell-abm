#include <qstring.h>
#include "params.h"

Params::Params()
{
    PARAM_SET params[] = {


{"NX", 100, 0, 0,
"Lattice size",
"Dimension of the lattice (number of sites in X,Y and Z directions).  Typically 5*BLOB_RADIUS is OK."},

{"INITIAL_COUNT", 5000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},



{"DIVIDE_TIME_MEDIAN", 18, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},


{"NDAYS", 0.1, 0, 0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 300, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.\n\
[days]"},

{"NT_CONC", 10, 0, 0,
"Number of ODE solver sub-steps.",
"The number of subdivisions of the major time step, for the ODE diffusion-reaction solver.\n\
[days]"},

{"TEST_CASE", 0, 0, 0,
"Test case #",
"Number of the test case to run.  The default value of 0 is for a normal run"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 3, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 20, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"INPUT_FILE", 0, 0, 0,
"spheroid_fixed.inpdata",
"The auxiliary input file contains data that (almost!) never changes"}

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
