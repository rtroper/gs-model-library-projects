
#include "pch.h"
#include "TimeSeries.h"

//  Global enumerations, useful for C-style implementations 

//  Every external function is called by GoldSim with specific requests including:
//    - initialization
//    - performing calculations
//    - obtaining the version number
//    - obtaining the number of input and output arguments
//    - "cleaning up" after the simulation (e.g. releasing memory and closing files)
//
//  XFMethodID (below) identifies the method types used to identify the phase of the simulation that is currently in progress.
//
//    XF_INITIALIZE    - Called after DLL is loaded and before each realization (no arguments are passed on this call). 
//    XF_CALCULATE     - Called during the simulation, each time the inputs change. 
//    XF_REP_VERSION   - Called after DLL load to report the external fcn version number. 
//    XF_REP_ARGUMENTS - Called after DLL load to report the number of input and output arguments. 
//    XF_CLEANUP       - Called before the DLL is unloaded. 

enum XFMethodID
{
	XF_INITIALIZE = 0,
	XF_CALCULATE = 1,
	XF_REP_VERSION = 2,
	XF_REP_ARGUMENTS = 3,
	XF_CLEANUP = 99
};

// Declare global variables that persist between calls from GoldSim to the DLL
TimeSeries* ts1;		// Use only for function GetTimeSeriesStatistics
TimeSeries* ts2;		// Use only for function GetTimeSeriesCorrelation
TimeSeries* ts3;		// Use only for function GetTimeSeriesCorrelation
int initialized1 = -1;	// Set to a value >= 0 on the first call to GetTimeSeriesStatistics
int initialized2 = -1;	// Set to a value >= 0 on the first call to GetTimeSeriesCorrelation
double time_index1;		// Use in functions where time values do not matter (use for GetTimeSeriesStatistics)
double time_index2;		// Use in functions where time values do not matter (use for GetTimeSeriesCorrelation)

//  XFStatusID (below) identifies the return codes for external functions. 
//      
//    XF_SUCCESS          – Call completed successfully. 
//    XF_CLEANUP_NOW      - Call was successful, but GoldSim should clean up and unload the DLL immediately. 
//    XF_FAILURE          - Failure (no error information returned).  
//    XF_FAILURE_WITH_MSG – Failure, with DLL-supplied error message available. Address of error message is returned in the first element of the output arguments array. 
//    XF_INCREASE_MEMORY  - Failed because the memory allocated for output arguments is too small.  GoldSim will increase the size of the output argument array and try again. 

// Return codes:
//   0 indicates OK, continue GoldSim
//   >0 and <99 indicates to terminate GoldSim
//   99 indicates OK, unload the DLL
//
// The following codes can only be used for XF_CALCULATE:
//   -1 indicates fatal error and an error message pointer is returned
//   -2 indicates more result memory is required; total amount in doubles is returned in outargs[0]
enum XFStatusID
{
	XF_SUCCESS = 0,
	XF_FAILURE = 1,
	XF_CLEANUP_NOW = 99,
	XF_FAILURE_WITH_MSG = -1,
	XF_INCREASE_MEMORY = -2
};

//  When calling methods from a DLL, GoldSim always expects the following C/C++ function signature:
//     extern "C" void __declspec(dllexport) MyExternalFcn(int XFMethod, int* XFState, double* inargs, double* outargs)
//
//  Notes:
//    - The function name and argument names can be different from those shown above
//    - "C" specifies C language style linkage between GoldSim and the DLL
//    - __declspec(dllexport) makes the function visible outside the DLL
//    - XFMethod: specifies the action that the external function must perform
//    - XFState: a code indicating the status of the external function
//    - inargs: array of input arguments
//    - outargs: this array returns different information for different XFMethod values

// Calculate basic statistics for values in a time series
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) GetTimeSeriesStatistics(int methodID, int* status, double* inargs, double* outargs)
{
	*status = XF_SUCCESS;
	//double current_x;

	switch (methodID)
	{
		// Initialize (called at beginning of each realization)
		// For example, allocate memory or open files
		// No arguments are passed on this call
	case  XF_INITIALIZE:
		// Initialize global variables
		break;	// no initialization in this example

	// The external function reports its version
	// No input arguments are passed on this call.
	// outargs[0] is set equal to the external fcn version
	case  XF_REP_VERSION:
		outargs[0] = 1.01;
		break;

		// External fcn reports the number of input and output arguments
		// outargs[0] is set equal to the # of inputs arguments
		// outargs[1] is set equal to the # of output arguments
	case  XF_REP_ARGUMENTS:
		// Two values from GoldSim are expected on each call to this function
		// The second value is designated to specify a numerical code for missing values in the time series
		outargs[0] = 2.0;

		// For now, there are two return values, the mean and standard deviation of the time series
		outargs[1] = 2.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize global variables
		if (initialized1 < 0)
		{
			// Initialize time index
			time_index1 = 0.0;

			// Set 'initialized1' >= 0 so that this block is only executed on the first call to the DLL
			initialized1 = 1;

			// Initialize the time series
			ts1 = new TimeSeries();
		}

		// Add a new time point to the time series if it is not a missing value
		if (inargs[0] > inargs[1] + 1.0e-6)
		{
			ts1->addTimepoint(time_index1, inargs[0]);
		}

		// Calculate and return the mean and standard deviation of the stored time series values
		outargs[0] = ts1->mean();
		outargs[1] = ts1->stdev();

		// Increment the time series index
		time_index1 += 1.0;

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		// Free memory allocated for time series object(s)
		delete ts1;
		ts1 = 0;
		break;

		// Error if this point is reached
		// This means the switch statement did not provide the cases that GoldSim expected.
		// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}

// Calculate the correlation between two time series
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) GetTimeSeriesCorrelation(int methodID, int* status, double* inargs, double* outargs)
{
	*status = XF_SUCCESS;
	//double current_x;

	switch (methodID)
	{
		// Initialize (called at beginning of each realization)
		// For example, allocate memory or open files
		// No arguments are passed on this call
	case  XF_INITIALIZE:
		// Initialize global variables
		break;	// no initialization in this example

	// The external function reports its version
	// No input arguments are passed on this call.
	// outargs[0] is set equal to the external fcn version
	case  XF_REP_VERSION:
		outargs[0] = 1.01;
		break;

		// External fcn reports the number of input and output arguments
		// outargs[0] is set equal to the # of inputs arguments
		// outargs[1] is set equal to the # of output arguments
	case  XF_REP_ARGUMENTS:
		// Three values from GoldSim are expected on each call to this function
		// The third value is designated to specify a numerical code for missing values in the time series
		outargs[0] = 3.0;

		// Return the value of the correlation between the two time series 
		// as well as the slope and intercept of the linear regression
		outargs[1] = 3.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize global variables
		if (initialized2 < 0)
		{
			// Initialize time index
			time_index2 = 0.0;

			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			initialized2 = 1;

			// Initialize the two time series
			ts2 = new TimeSeries();
			ts3 = new TimeSeries();
		}

		// Add new time points to the two time series only if neither has a missing value at the current time point
		if ( (inargs[0] > inargs[2] + 1.0e-6) && (inargs[1] > inargs[2] + 1.0e-6) )
		{
			ts2->addTimepoint(time_index2, inargs[0]);
			ts3->addTimepoint(time_index2, inargs[1]);
		}

		// Calculate and return the correlation of values in the two time series
		outargs[0] = ts2->correlation(*ts3);
		outargs[1] = ts2->linear_regression_slope(*ts3);
		outargs[2] = ts2->linear_regression_intercept(*ts3);

		// Increment the time series index
		time_index2 += 1.0;

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		// Free memory allocated for time series object(s)
		delete ts2;
		delete ts3;
		ts2 = ts3 = 0;
		break;

		// Error if this point is reached
		// This means the switch statement did not provide the cases that GoldSim expected.
		// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}
