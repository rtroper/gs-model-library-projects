#include "pch.h"
#include "SplineInterpolator.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

bool SplineInterpolator::store_input_data(const std::vector<double>& x_in, const std::vector<double>& y_in)
{
	// Default return value is true
	bool success = true;

	// Store data up to the maximum available (for x and y)
	size_data = static_cast<size_t>(std::min(x_in.size(), y_in.size()));

	// Store source data in x and y and calculate h
	for (size_t i = 0; i < size_data; i++)
	{
		x.push_back(x_in[i]);
		y.push_back(y_in[i]);

		if (i < size_data - 1)
		{
			// Store difference between x[i+1] and x[i]
			h.push_back(x_in[i + 1] - x_in[i]);
		}
	}

	return success;
}

// Left-hand-side matrix (lhs) and right-hand-side vector (rhs) must not be initialized until
// store_input_data has been called since y and h must first be initialized
bool SplineInterpolator::initialize_lhs_and_rhs(gsl_matrix* lhs, gsl_vector* rhs)
{
	// Default return value is true
	bool success = true;

	// NOTE: lhs and rhs are used to solve for first derivatives D0, D1, ..., DN at data points
	// (x0, y0), (x1, y1), ... , (xN, yN), where N = size-1, and second derivatives D0' and DN' at 
	// (x0, y0) and (xN, yN). Thus, lhs, rhs, and the resulting solution are of size 'size + 2',
	// where 'size' is the number of data points through which the spline function passes. In the
	// resulting solution vector, the first 'size' values are D0, D1, ... , DN and the last 2 values
	// are D0' and DN'.

	// Initialize the first 'size_data' rows of lhs and rhs
	for (size_t i = 0; i < size_data; i++)
	{
		if (i == 0)
		{
			// Set lhs values
			gsl_matrix_set(lhs, 0,    0,          2.0 * h[0] );
			gsl_matrix_set(lhs, 0,    1,                h[0] );
			gsl_matrix_set(lhs, 0, size_data, 0.5 * pow(h[0], 2.0));

			// Set rhs value
			gsl_vector_set(rhs, 0, 3.0 * (y[1] - y[0]));
		}
		else if (i == size_data - 1)
		{
			// Set lhs values
			gsl_matrix_set(lhs, size_data - 1, size_data - 1,      2.0 * h[size_data - 2]);
			gsl_matrix_set(lhs, size_data - 1, size_data - 2,            h[size_data - 2]);
			gsl_matrix_set(lhs, size_data - 1, size_data + 1, -0.5 * pow(h[size_data - 2], 2.0));

			// Set rhs value
			gsl_vector_set(rhs, size_data - 1, 3.0 * (y[size_data - 1] - y[size_data - 2]));
		}
		else
		{
			// Set lhs values
			gsl_matrix_set(lhs, i, i - 1,                           h[i - 1]);
			gsl_matrix_set(lhs, i, i    , 2.0 * h[i - 1] * (1 + pow(h[i - 1] / h[i], 2.0)) );
			gsl_matrix_set(lhs, i, i + 1,       h[i - 1] *      pow(h[i - 1] / h[i], 2.0)  );

			// Set rhs value
			gsl_vector_set(rhs, i, 3.0 * (y[i + 1] - y[i]) * pow(h[i - 1] / h[i], 2.0) + 3.0 * (y[i] - y[i - 1]));
		}
	}

	// Set values in last 2 rows of lhs to specify boundary conditions. The default is for
	// boundary first derivatives to be equal to each other and boundary second derivatives to be
	// equal to each other. Note that for these conditions, the last 2 rows of rhs are zero and,
	// therefore, do not need to be modified. TO-DO: Create a separate function that can be used
	// to set different boundary conditions.
	gsl_matrix_set(lhs,     size_data,             0,  1.0);
	gsl_matrix_set(lhs,     size_data, size_data - 1, -1.0);
	gsl_matrix_set(lhs, size_data + 1,     size_data,  1.0);
	gsl_matrix_set(lhs, size_data + 1, size_data + 1, -1.0);

	return success;
}

bool SplineInterpolator::calculate_derivatives(gsl_vector* D)
{
	// Default return value is true
	bool success = true;

	// Define integer for error returns
	int status = 0;

	// Initialize LU matrix and right-hand-side (rhs) vector
	gsl_matrix* LU = gsl_matrix_calloc(size_data + 2, size_data + 2);
	gsl_vector* rhs = gsl_vector_calloc(size_data + 2);
	success = initialize_lhs_and_rhs(LU, rhs);
	if (!success) return false;	// Add code for proper error handling/messaging

	// Create permutation vector to store permutation
	gsl_permutation* p = gsl_permutation_calloc(size_data + 2);

	// Create pointer to int for third argument to LU decomposition function
	int* signum = new int(0);

	// Get the LU decomposition
	status = gsl_linalg_LU_decomp(LU, p, signum);
	if (status) return false;	// Add code for proper error handling/messaging

	// Solve for D
	const gsl_vector* rhs_const = rhs;
	status = gsl_linalg_LU_solve(LU, p, rhs, D);
	if (status) return false;	// Add code for proper error handling/messaging

	// Free memory used for solving matrix equations
	gsl_matrix_free(LU);
	gsl_vector_free(rhs);
	gsl_permutation_free(p);
	delete signum;

	return success;
}

bool SplineInterpolator::calculate_spline_parameters()
{
	// Default return value is true
	bool success = true;

	// Initialize vector for derivatives
	gsl_vector* D = gsl_vector_calloc(size_data + 2);
	success = calculate_derivatives(D);
	if (!success) return false;	// Add code for proper error handling/messaging

	// Load calculated values into a, b, c, and d
	double D1 = 0.0, D2 = 0.0;
	for (size_t i = 0; i < size_data - 1; i++)
	{
		// Get D values needed for calculations
		D1 = gsl_vector_get(D,     i);
		D2 = gsl_vector_get(D, i + 1);

		// Calculate a, b, c, and d values
		a.push_back(y[i]);
		b.push_back(h[i] * D1);
		c.push_back(3.0 * (y[i + 1] - y[i]) - h[i] * (2.0 * D1 + D2));
		d.push_back(2.0 * (y[i] - y[i + 1]) + h[i] * (D1 + D2));
	}

	// Free memory that is no longer needed
	gsl_vector_free(D);

	return success;
}

bool SplineInterpolator::initialize_default_interpolator()
{
	// Default return value is true
	bool success = true;

	// Calculate spline parameters a, b, c, and d
	success = calculate_spline_parameters();
	if (!success) return false;	// Add code for proper error handling/messaging

	return success;
}

bool SplineInterpolator::initialize_cspline_interpolator(SplineBoundaryCondition boundary)
{
	// Default return value is true
	bool success = true;

	// Define integer for error returns
	int status = 0;

	// Initialize the cubic spline interpolator
	switch (boundary)
	{
	case(SplineBoundaryCondition::Natural):
		gsl_interpolator = gsl_interp_alloc(gsl_interp_cspline, size_data);
		break;
	case(SplineBoundaryCondition::Periodic):
		gsl_interpolator = gsl_interp_alloc(gsl_interp_cspline_periodic, size_data);
		break;
	default:
		gsl_interpolator = gsl_interp_alloc(gsl_interp_cspline, size_data);
		break;
	}
	status = gsl_interp_init(gsl_interpolator, x.data(), y.data(), size_data);
	if (status) return false;	// Add code for proper error handling/messaging

	// Initialize the accelerator for cubic spline interpolation
	gsl_accelerator = gsl_interp_accel_alloc();

	return success;
}

bool SplineInterpolator::initialize_steffen_interpolator()
{
	// Default return value is true
	bool success = true;

	// Define integer for error returns
	int status = 0;

	// Initialize the Steffen interpolator
	gsl_interpolator = gsl_interp_alloc(gsl_interp_steffen, size_data);
	status = gsl_interp_init(gsl_interpolator, x.data(), y.data(), size_data);
	if (status) return false;	// Add code for proper error handling/messaging

	// Initialize the accelerator for Steffen interpolation
	gsl_accelerator = gsl_interp_accel_alloc();

	return success;
}

bool SplineInterpolator::initialize_akima_interpolator(SplineBoundaryCondition boundary)
{
	// Default return value is true
	bool success = true;

	// Define integer for error returns
	int status = 0;

	// Initialize the Akima interpolator
	switch (boundary)
	{
	case(SplineBoundaryCondition::Natural):
		gsl_interpolator = gsl_interp_alloc(gsl_interp_akima, size_data);
		break;
	case(SplineBoundaryCondition::Periodic):
		gsl_interpolator = gsl_interp_alloc(gsl_interp_akima_periodic, size_data);
		break;
	default:
		gsl_interpolator = gsl_interp_alloc(gsl_interp_akima, size_data);
		break;
	}
	status = gsl_interp_init(gsl_interpolator, x.data(), y.data(), size_data);
	if (status) return false;	// Add code for proper error handling/messaging

	// Initialize the accelerator for Akima interpolation
	gsl_accelerator = gsl_interp_accel_alloc();

	return success;
}

bool SplineInterpolator::initialize_interpolator(const std::vector<double>& x_in, const std::vector<double>& y_in)
{
	// Default return value is true
	bool success = true;

	// NOTE: store_input_data must be called before anything else, since it
	// initializes the value of size_data

	// Store data and calculate intervals h (also initializes size_data)
	success = store_input_data(x_in, y_in);
	if (!success) return false;	// Add code for proper error handling/messaging

	// Switch on interpolation method
	switch (method)
	{
	case(InterpolationMethod::BasicNonGSL):
		success = initialize_default_interpolator();
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	case(InterpolationMethod::GSLCubic):
		success = initialize_cspline_interpolator(SplineBoundaryCondition::Natural);
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	case(InterpolationMethod::GSLCubicPeriodic):
		success = initialize_cspline_interpolator(SplineBoundaryCondition::Periodic);
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	case(InterpolationMethod::GSLSteffen):
		success = initialize_steffen_interpolator();
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	case(InterpolationMethod::GSLAkima):
		success = initialize_akima_interpolator(SplineBoundaryCondition::Natural);
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	case(InterpolationMethod::GSLAkimaPeriodic):
		success = initialize_akima_interpolator(SplineBoundaryCondition::Periodic);
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	default:
		success = initialize_default_interpolator();
		if (!success) return false;	// Add code for proper error handling/messaging
		break;
	}

	return success;
}

std::pair<int, double> SplineInterpolator::get_i_and_t_values(double x_value)
{
	// Initialize return value
	std::pair<int, double> i_and_t = std::pair<int, double>();

	// First mod the input x value with the maximum x so that values outside x range are wrapped
	double x_adjusted = std::fmod(x_value, x[size_data - 1]);

	// Loop over x to find i such that x[i] < x_value < x[i+1]
	int i = 0;
	for (size_t idx = 1; idx < size_data; idx++)
	{
		if (x_adjusted > x[idx]) i++;
		else break;
	}

	// Calculate t value i.e. (x_value - x[i]) / (x[i+1] - x[i])
	double t = (x_adjusted - x[i]) / h[i];

	// Store i and t values to return
	i_and_t.first = i;
	i_and_t.second = t;

	return i_and_t;
}

SplineInterpolator::SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source)
	: method(InterpolationMethod::BasicNonGSL), size_data(0), gsl_interpolator(0), gsl_accelerator(0)
{
	// Call initialization function
	initialize_interpolator(x_source, y_source);
}

SplineInterpolator::SplineInterpolator(InterpolationMethod method, const std::vector<double>& x_source, const std::vector<double>& y_source)
	: method(method), size_data(0), gsl_interpolator(0), gsl_accelerator(0)
{
	// Call initialization function
	initialize_interpolator(x_source, y_source);
}

SplineInterpolator::~SplineInterpolator()
{
	// Free memory allocated for GSL interpolation objects
	if (gsl_interpolator) gsl_interp_free(gsl_interpolator);
	if (gsl_accelerator) gsl_interp_accel_free(gsl_accelerator);
}

double SplineInterpolator::interpolate(double x_value)
{
	double y_interp = 0.0;

	if (method == InterpolationMethod::BasicNonGSL)
	{
		// Get i and t values that correspond to the input x value
		std::pair<int, double> i_and_t = get_i_and_t_values(x_value);
		int i = i_and_t.first;
		double t = i_and_t.second;

		// Calculate the interpolated y value
		y_interp = a[i] + b[i] * t + c[i] * std::pow(t, 2.0) + d[i] * std::pow(t, 3.0);
	}
	else
	{
		// First mod the input x value with the maximum x so that values outside x range are wrapped
		double x_adjusted = std::fmod(x_value, x[size_data - 1]);

		// Get the interpolated value
		y_interp = gsl_interp_eval(gsl_interpolator, x.data(), y.data(), x_adjusted, gsl_accelerator);
	}

	return y_interp;
}

double SplineInterpolator::interpolate_derivative(double x_value)
{
	double deriv_interp = 0.0;

	if (method == InterpolationMethod::BasicNonGSL)
	{
		// Get i and t values that correspond to the input x value
		std::pair<int, double> i_and_t = get_i_and_t_values(x_value);
		int i = i_and_t.first;
		double t = i_and_t.second;

		// Calculate the interpolated y value
		deriv_interp = (b[i] + 2.0 * c[i] * t + 3.0 * d[i] * std::pow(t, 2.0)) / h[i];
	}
	else
	{
		// First mod the input x value with the maximum x so that values outside x range are wrapped
		double x_adjusted = std::fmod(x_value, x[size_data - 1]);

		// Get the interpolated derivative
		deriv_interp = gsl_interp_eval_deriv(gsl_interpolator, x.data(), y.data(), x_adjusted, gsl_accelerator);
	}

	return deriv_interp;
}

