
#include "pch.h"
#include "SplineInterpolator.h"

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
SplineInterpolator* interpolator;
int initialized = -1;		// This should be set to a value >= 0 on the first call to the DLL
int DATA_BUFFER_SIZE = 50;	// Max size of x and y data

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

// Calculate the spline-interpolated values (non-GSL) for given x and y input vectors of data
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) Interpolate(int methodID, int* status, double* inargs, double* outargs)
{
	*status = XF_SUCCESS;
	double current_x;

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
		// The first argument from GS is the actual data size (e.g. if it is 20, then x is size 20 and y is size 20).
		// The next 2 * DATA_BUFFER_SIZE values are reserved for x and y so that DATA_BUFFER_SIZE is the maximum size
		// of x and y. The last argument from GS is the x value for which to calculate the interpolated value.
		outargs[0] = 2.0 * DATA_BUFFER_SIZE + 2.0;

		// The first return argument is the interpolated value and the second is the derivative of the fitted spline
		outargs[1] = 2.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize the spline interpolator
		if (initialized < 0)
		{
			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			initialized = 1;

			// Initialize x and y vectors to store input data
			int size_data = min(static_cast<int>(inargs[0]), DATA_BUFFER_SIZE);
			std::vector<double> x(size_data, 0.0), y(size_data, 0.0);

			// Load input data into x and y
			for (int i = 0; i < size_data; i++)
			{
				x[i] = inargs[i + 1];
				y[i] = inargs[i + 1 + DATA_BUFFER_SIZE];
			}

			// Initialize the spline interpolator
			interpolator = new SplineInterpolator(x, y);
		}

		// Calculate and return the interpolated value and derivative
		current_x = inargs[2 * DATA_BUFFER_SIZE + 1];
		outargs[0] = interpolator->interpolate(current_x);
		outargs[1] = interpolator->interpolate_derivative(current_x);

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		break;	// No clean-up required

	// Error if this point is reached
	// This means the switch statement did not provide the cases that GoldSim expected.
	// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}

// Calculate the spline-interpolated values (GSL cubic spline) for given x and y input vectors of data
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) Interpolate_CSpline(int methodID, int* status, double* inargs, double* outargs)
{
	*status = XF_SUCCESS;
	double current_x;

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
		// The first argument from GS is the actual data size (e.g. if it is 20, then x is size 20 and y is size 20).
		// The next 2 * DATA_BUFFER_SIZE values are reserved for x and y so that DATA_BUFFER_SIZE is the maximum size
		// of x and y. The last argument from GS is the x value for which to calculate the interpolated value.
		outargs[0] = 2.0 * DATA_BUFFER_SIZE + 2.0;

		// The first return argument is the interpolated value and the second is the derivative of the fitted spline
		outargs[1] = 2.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize the spline interpolator
		if (initialized < 0)
		{
			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			initialized = 1;

			// Initialize x and y vectors to store input data
			int size_data = min(static_cast<int>(inargs[0]), DATA_BUFFER_SIZE);
			std::vector<double> x(size_data, 0.0), y(size_data, 0.0);

			// Load input data into x and y
			for (int i = 0; i < size_data; i++)
			{
				x[i] = inargs[i + 1];
				y[i] = inargs[i + 1 + DATA_BUFFER_SIZE];
			}

			// Initialize the spline interpolator
			interpolator = new SplineInterpolator(x, y);
		}

		// Calculate and return the interpolated value and derivative
		current_x = inargs[2 * DATA_BUFFER_SIZE + 1];
		outargs[0] = interpolator->interpolate_cspline(current_x);
		outargs[1] = interpolator->interpolate_cspline_derivative(current_x);

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		break;	// No clean-up required

	// Error if this point is reached
	// This means the switch statement did not provide the cases that GoldSim expected.
	// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}

// Calculate the spline-interpolated values (GSL Steffen spline) for given x and y input vectors of data
//-----------------------------------------------------------------------------------------------
extern "C" void __declspec(dllexport) Interpolate_SteffenSpline(int methodID, int* status, double* inargs, double* outargs)
{
	*status = XF_SUCCESS;
	double current_x;

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
		// The first argument from GS is the actual data size (e.g. if it is 20, then x is size 20 and y is size 20).
		// The next 2 * DATA_BUFFER_SIZE values are reserved for x and y so that DATA_BUFFER_SIZE is the maximum size
		// of x and y. The last argument from GS is the x value for which to calculate the interpolated value.
		outargs[0] = 2.0 * DATA_BUFFER_SIZE + 2.0;

		// The first return argument is the interpolated value and the second is the derivative of the fitted spline
		outargs[1] = 2.0;
		break;

		// Normal calculation.
		// Results are returned as outarg[0], outarg[1], etc. depending on number of outputs
	case  XF_CALCULATE:
		// Initialize the spline interpolator
		if (initialized < 0)
		{
			// Set 'initialized' >= 0 so that this block is only executed on the first call to the DLL
			initialized = 1;

			// Initialize x and y vectors to store input data
			int size_data = min(static_cast<int>(inargs[0]), DATA_BUFFER_SIZE);
			std::vector<double> x(size_data, 0.0), y(size_data, 0.0);

			// Load input data into x and y
			for (int i = 0; i < size_data; i++)
			{
				x[i] = inargs[i + 1];
				y[i] = inargs[i + 1 + DATA_BUFFER_SIZE];
			}

			// Initialize the spline interpolator
			interpolator = new SplineInterpolator(x, y);
		}

		// Calculate and return the interpolated value and derivative
		current_x = inargs[2 * DATA_BUFFER_SIZE + 1];
		outargs[0] = interpolator->interpolate_steffen(current_x);
		outargs[1] = interpolator->interpolate_steffen_derivative(current_x);

		break;

		// Close any open files
		// Optionally release any memory that's been allocated
		// No arguments are passed on this call.
	case  XF_CLEANUP:
		break;	// No clean-up required

	// Error if this point is reached
	// This means the switch statement did not provide the cases that GoldSim expected.
	// The external function must have cases 0, 1, 2, 3 and 99
	default:
		*status = XF_FAILURE;
		break;
	}
}
