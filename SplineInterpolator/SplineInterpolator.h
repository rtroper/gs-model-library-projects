#pragma once

#include <vector>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>

class SplineInterpolator
{
public:
	// Enum to specify the interpolation method
	enum InterpolationMethod
	{
		BasicNonGSL = 0,
		GSLCubic,
		GSLCubicPeriodic,
		GSLSteffen,
		GSLAkima,
		GSLAkimaPeriodic
	};

	// Enum to specify spline boundary conditions
	enum SplineBoundaryCondition
	{
		Natural = 0,
		Periodic
	};

private:
	InterpolationMethod method;
	size_t size_data;				// Number of data points to which to fit spline
	std::vector<double> x, y, h;	// Independent (x) and dependent (y) variable data and x-variable intervals (h)
	std::vector<double> a, b, c, d; // Coefficients for spline interpolation

	// Variables below are for spline interpolation using GSL's built-in spline interpolation functions
	gsl_interp* gsl_interpolator;
	gsl_interp_accel* gsl_accelerator;

	// Store input data in member variables (return true if success)
	bool store_input_data(const std::vector<double>& x_in, const std::vector<double>& y_in);

	// Initialize left-hand-side matrix and right-hand-side vector for derivative calculations (return true if success)
	bool initialize_lhs_and_rhs(gsl_matrix* lhs, gsl_vector* rhs);

	// Calculate first derivatives at all data points and second derivatives at bounding points
	bool calculate_derivatives(gsl_vector* D);

	// Calculate spline parameter vectors a, b, c, and d
	bool calculate_spline_parameters();

	// Initialize lhs, rhs, a, b, c, d, etc used in default interpolation method
	bool initialize_default_interpolator();

	// Initialize GSL cspline interpolator
	bool initialize_cspline_interpolator(SplineBoundaryCondition boundary);

	// Initialize GSL Steffen interpolator
	bool initialize_steffen_interpolator();

	// Initialize GSL Akima interpolator
	bool initialize_akima_interpolator(SplineBoundaryCondition boundary);

	// Store input data, initialize lhs and rhs and calculate derivatives
	bool initialize_interpolator(const std::vector<double>& x_in, const std::vector<double>& y_in);

	// Calculate i and t values for a specified x value
	std::pair<int, double> get_i_and_t_values(double x_value);

public:
	// Constructor that uses default interpolator (InterpolationMethod::BasicNonGSL)
	SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source);

	// Constructor that allows specification of the interpolation method
	SplineInterpolator(InterpolationMethod method, const std::vector<double>& x_source, const std::vector<double>& y_source);

	// Destructor to free memory allocated by GSL
	~SplineInterpolator();

	// Calculate spline interpolation value
	double interpolate(double x_value);

	// Calculate spline interpolation derivative
	double interpolate_derivative(double x_value);
};

