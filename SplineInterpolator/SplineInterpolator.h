#pragma once

#include <vector>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>

class SplineInterpolator
{
private:
	size_t size_data;				// Number of data points to which to fit spline
	gsl_vector* x, * y;				// Independent and dependent variable data to which to fit spline functions

	// Variables below are for spline interpolation without GSL's built-in spline interpolation functions

	gsl_matrix* lhs;				// Left-hand-side matrix in equation to solve for derivatives
	gsl_vector* rhs;				// Right-hand-side vector in equation to solve for derivatives
	gsl_vector* a, * b, * c, * d;	// Coefficients in spline functions
	gsl_vector* D, * h;				// Derivatives (D) and x-variable intervals (h)

	// Variables below are for spline interpolation using GSL's built-in spline interpolation functions

	gsl_interp* cspline_interp, * steffen_interp, * akima_interp;		// Interpolation objects for cubic spline and Steffen interpolation
	gsl_interp_accel* cspline_accel, * steffen_accel, * akima_accel;	// Accelerator objects used for cubic spline and Steffen interpolation

	// Store input data in member variables (return true if success)
	bool store_input_data(const double* x_source, int x_size, const double* y_source, int y_size);

	// Initialize left-hand-side matrix and right-hand-side vector for derivative calculations (return true if success)
	bool initialize_lhs_and_rhs();

	// Calculate first derivatives at all data points and second derivatives at bounding points
	bool calculate_derivatives();

	// Calculate spline parameter vectors a, b, c, and d
	bool calculate_spline_parameters();

	// Initialize GSL cspline interpolator
	bool initialize_cspline_interpolator();

	// Initialize GSL Steffen interpolator
	bool initialize_steffen_interpolator();

	// Initialize GSL Akima interpolator
	bool initialize_akima_interpolator();

	// Store input data, initialize lhs and rhs and calculate derivatives
	bool initialize_interpolators(const double* x_source, int x_size, const double* y_source, int y_size);

	// Calculate i and t values for a specified x value
	std::pair<int, double> get_i_and_t_values(double x_value);

public:
	// Constructors
	SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source);
	SplineInterpolator(const double* x_source, int x_size, const double* y_source, int y_size);

	// Calculate spline interpolation value
	double interpolate(double x_value);

	// Calculate spline interpolation derivative
	double interpolate_derivative(double x_value);

	// Calculate spline interpolation value using GSL's cubic spline
	double interpolate_cspline(double x_value);

	// Calculate spline interpolation derivative using GSL's cubic spline derivative
	double interpolate_cspline_derivative(double x_value);

	// Calculate spline interpolation value using GSL's Steffen spline
	double interpolate_steffen(double x_value);

	// Calculate spline interpolation derivative using GSL's Steffen spline derivative
	double interpolate_steffen_derivative(double x_value);

	// Calculate spline interpolation value using GSL's Akima spline
	double interpolate_akima(double x_value);

	// Calculate spline interpolation derivative using GSL's Akima spline derivative
	double interpolate_akima_derivative(double x_value);

	// Convenience functions to write vector and matrix values to check for correctness
	void write_vector_values(std::string file_name);
	void write_matrix_values(std::string file_name);
};

