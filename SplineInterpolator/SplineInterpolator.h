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
		GSLSteffen,
		GSLAkima
	};

private:
	InterpolationMethod method;
	size_t size_data;				// Number of data points to which to fit spline
	gsl_vector* x, * y;				// Independent and dependent variable data to which to fit spline functions

	// Variables below are for spline interpolation without GSL's built-in spline interpolation functions

	gsl_matrix* lhs;				// Left-hand-side matrix in equation to solve for derivatives
	gsl_vector* rhs;				// Right-hand-side vector in equation to solve for derivatives
	gsl_vector* a, * b, * c, * d;	// Coefficients in spline functions
	gsl_vector* D, * h;				// Derivatives (D) and x-variable intervals (h)

	// Variables below are for spline interpolation using GSL's built-in spline interpolation functions

	gsl_interp* gsl_interpolator;
	gsl_interp_accel* gsl_accelerator;

	// Store input data in member variables (return true if success)
	bool store_input_data(const double* x_source, int x_size, const double* y_source, int y_size);

	// Initialize left-hand-side matrix and right-hand-side vector for derivative calculations (return true if success)
	bool initialize_lhs_and_rhs();

	// Calculate first derivatives at all data points and second derivatives at bounding points
	bool calculate_derivatives();

	// Calculate spline parameter vectors a, b, c, and d
	bool calculate_spline_parameters();

	// Initialize lhs, rhs, a, b, c, d, etc used in default interpolation method
	bool initialize_default_interpolator();

	// Initialize GSL cspline interpolator
	bool initialize_cspline_interpolator();

	// Initialize GSL Steffen interpolator
	bool initialize_steffen_interpolator();

	// Initialize GSL Akima interpolator
	bool initialize_akima_interpolator();

	// Store input data, initialize lhs and rhs and calculate derivatives
	bool initialize_interpolator(const double* x_source, int x_size, const double* y_source, int y_size);

	// Calculate i and t values for a specified x value
	std::pair<int, double> get_i_and_t_values(double x_value);

public:
	// Constructors that use default interpolator (InterpolationMethod::BasicNonGSL)
	SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source);
	SplineInterpolator(const double* x_source, int x_size, const double* y_source, int y_size);

	// Constructors that allow specification of the interpolator
	SplineInterpolator(InterpolationMethod method, const std::vector<double>& x_source, const std::vector<double>& y_source);
	SplineInterpolator(InterpolationMethod method, const double* x_source, int x_size, const double* y_source, int y_size);

	// Destructor to free memory allocated by GSL
	~SplineInterpolator();

	// Calculate spline interpolation value
	double interpolate(double x_value);

	// Calculate spline interpolation derivative
	double interpolate_derivative(double x_value);

	// Convenience functions to write vector and matrix values to check for correctness
	void write_vector_values(std::string file_name);
	void write_matrix_values(std::string file_name);
};

