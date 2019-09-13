
#include "pch.h"
#include "TimeSeries.h"

//  Global enumerations, useful for C-style implementations 

//  XFMethodID (below) identifies the method types used to identify the phase of the simulation that is currently in progress.
//
//    XF_INITIALIZE    - Called after DLL is loaded and before each realization (no arguments are passed on this call). 
//    XF_CALCULATE     - Called during the simulation, each time the inputs change. 
//    XF_REP_VERSION   - Called after DLL load to report the external function version number. 
//    XF_REP_ARGUMENTS - Called after DLL load to report the number of input and output arguments. 
//    XF_CLEANUP       - Called before the DLL is unloaded to "clean up" after the simulation (e.g. release memory and close files)

enum XFMethodID
{
	XF_INITIALIZE = 0,
	XF_CALCULATE = 1,
	XF_REP_VERSION = 2,
	XF_REP_ARGUMENTS = 3,
	XF_CLEANUP = 99
};

// Declare global variables that persist between calls from GoldSim to the DLL
enum TSDefinitionConstants
{
	TS_START_IDX = 0,
	TS_SIZE_IDX = 7,
	TS_DATA_START_IDX = 8
};

TimeSeries* ts4;		// Use only for function GetTimeSeriesAutoCorrelation
TimeSeries* ts5;		// Use only for function GetTimeSeriesAutoCorrelation

int initialized3 = -1;	// Set to a value >= 0 on the first call to GetTimeSeriesAutoCorrelation

double time_index3;		// Use in functions where time values do not matter (use for GetTimeSeriesAutoCorrelation)
int counter;			// Use only in GetTimeSeriesAutoCorrelation

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

// Reserved for function GetTimeSeriesStatistics
namespace TSStatistics
{
	TimeSeries* ts;			// Stores series times and values
	int initialized = -1;	// Set to a value >= 0 on the first call to the function
}

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
		// Two values from GoldSim are expected when this function is called
		// The first value is designated to specify a numerical code for missing values in the time series
		// The second value is the first value in the time series definition
		outargs[0] = 2.0;

		// For now, there are two return values, the mean and standard deviation of the time series
		outargs[1] = 2.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize global variables
		if (TSStatistics::initialized < 0)
		{
			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			TSStatistics::initialized = 1;

			// Initialize the time series and specify missing value code
			TSStatistics::ts = new TimeSeries();
			TSStatistics::ts->setMissingValueCode(inargs[0]);

			// Confirm that GoldSim is providing a time series definition (if so, load times and values of the time series definition)
			// Note that a value of 20 is a necessary but insufficient requirement (i.e. a minimum requirement)
			// for this to be a time series definition
			if (int(inargs[TS_START_IDX+1]) == 20)
			{
				// Assume that this is a time series definition and get the number of data points
				int number_of_data_points = int(inargs[TS_SIZE_IDX+1]);

				// Load time and data values
				for (int i = 0; i < number_of_data_points; i++)
				{
					TSStatistics::ts->addTimepoint(inargs[TS_DATA_START_IDX + 1 + i], inargs[TS_DATA_START_IDX + 1 + number_of_data_points + i]);
				}
			}
		}

		// Calculate and return the mean and standard deviation of the stored time series values
		outargs[0] = TSStatistics::ts->mean();
		outargs[1] = TSStatistics::ts->stdev();

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		// Free memory allocated for time series object(s)
		delete TSStatistics::ts;
		TSStatistics::ts = 0;
		break;

		// Error if this point is reached
		// This means the switch statement did not provide the cases that GoldSim expected.
		// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}

// Reserved for function GetTimeSeriesCorrelation
namespace TSCorrelation
{
	TimeSeries* ts1;		// Stores first series times and values
	TimeSeries* ts2;		// Stores second series times and values
	int initialized = -1;	// Set to a value >= 0 on the first call to the function
}

// Calculate the correlation between two time series
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) GetTimeSeriesCorrelation(int methodID, int* status, double* inargs, double* outargs)
{
	*status = XF_SUCCESS;
	std::pair<double, double> linear_regression_coefs;

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
		if (TSCorrelation::initialized < 0)
		{
			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			TSCorrelation::initialized = 1;

			// Initialize the two time series
			TSCorrelation::ts1 = new TimeSeries();
			TSCorrelation::ts2 = new TimeSeries();

			// Specify missing value codes
			TSCorrelation::ts1->setMissingValueCode(inargs[0]);
			TSCorrelation::ts2->setMissingValueCode(inargs[0]);

			// Confirm that GoldSim is providing a time series definition (if so, load times and values of the time series definition)
			// Note that a value of 20 is a necessary but insufficient requirement (i.e. a minimum requirement)
			// for this to be a time series definition
			if (int(inargs[TS_START_IDX + 1]) == 20)
			{
				// Assume that this is a time series definition and get the number of data points
				int number_of_data_points1 = int(inargs[TS_SIZE_IDX + 1]);

				// Load time and data values
				for (int i = 0; i < number_of_data_points1; i++)
				{
					TSCorrelation::ts1->addTimepoint(inargs[TS_DATA_START_IDX + 1 + i], inargs[TS_DATA_START_IDX + 1 + number_of_data_points1 + i]);
				}

				// Now, confirm that the next index after the first time series definition also contains 20
				// If so, assume that another time series definition follows and load times and values
				if (int(inargs[TS_DATA_START_IDX + 1 + 2 * number_of_data_points1]) == 20)
				{
					// Calculate indexes for retrieving information about the second time series definition
					int ts2_size_idx = TS_DATA_START_IDX + 1 + 2 * number_of_data_points1 + TS_SIZE_IDX;
					int number_of_data_points2 = int(inargs[ts2_size_idx]);
					int ts2_data_start_idx = ts2_size_idx + 1;

					// Load time and data values
					for (int i = 0; i < number_of_data_points2; i++)
					{
						TSCorrelation::ts2->addTimepoint(inargs[ts2_data_start_idx + i], inargs[ts2_data_start_idx + number_of_data_points2 + i]);
					}
				}
			}
		}

		// Calculate and return the correlation of values in the two time series
		outargs[0] = TSCorrelation::ts1->correlation(*TSCorrelation::ts2);

		// Calculate and return the linear regression slope and intercept for the two time series
		linear_regression_coefs = TSCorrelation::ts1->linear_regression_coefs(*TSCorrelation::ts2);
		outargs[1] = linear_regression_coefs.first;
		outargs[2] = linear_regression_coefs.second;

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		// Free memory allocated for time series object(s)
		delete TSCorrelation::ts1;
		delete TSCorrelation::ts2;
		TSCorrelation::ts1 = TSCorrelation::ts2 = 0;
		break;

		// Error if this point is reached
		// This means the switch statement did not provide the cases that GoldSim expected.
		// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}

// Calculate the auto-correlation for a given time series and specified shift or offset
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) GetTimeSeriesAutoCorrelation(int methodID, int* status, double* inargs, double* outargs)
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
		// The first value is the current value of the time series
		// The second value is the shift for the auto-correlation calculation
		// The third value is designated to specify a numerical code for missing values in the time series
		outargs[0] = 3.0;

		// Return the value of the auto-correlation of the time series
		outargs[1] = 1.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize global variables
		if (initialized3 < 0)
		{
			// Initialize time index
			time_index3 = 0.0;

			// Initialize counter
			counter = 0;

			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			initialized3 = 1;

			// Initialize the two time series
			ts4 = new TimeSeries();
			ts5 = new TimeSeries();
		}

		// Add new time point to the time series only if it is not a missing value at the current time point
		if ( inargs[0] > inargs[2] + 1.0e-6 )
		{
			// Add time point and increment counter
			ts4->addTimepoint(time_index3, inargs[0]);
			counter++;

			// Add time point from ts4 to ts5 if counter is >= auto-correlation shift 
			if (counter > int(inargs[1]))
			{
				ts5->addTimepoint(time_index3, inargs[0]);
			}
		}

		// Calculate and return the correlation of values in the two time series
		outargs[0] = ts4->correlation(*ts5);

		// Increment the time series index
		time_index3 += 1.0;

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		// Free memory allocated for time series object(s)
		delete ts4;
		delete ts5;
		ts4 = ts5 = 0;
		break;

		// Error if this point is reached
		// This means the switch statement did not provide the cases that GoldSim expected.
		// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}
