#pragma once

#include <vector>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class SplineInterpolator
{
private:
	size_t size_data;				// Number of data points to which to fit spline
	gsl_vector* x, * y;				// Independent and dependent variable data to which to fit spline functions
	gsl_matrix* lhs;				// Left-hand-side matrix in equation to solve for derivatives
	gsl_vector* rhs;				// Right-hand-side vector in equation to solve for derivatives
	gsl_vector* a, * b, * c, * d;	// Coefficients in spline functions
	gsl_vector* D, * h;				// Derivatives (D) and x-variable intervals (h)

	// Store input data in member variables (return true if success)
	bool store_input_data(const double* x_source, int x_size, const double* y_source, int y_size);

	// Initialize left-hand-side matrix and right-hand-side vector for derivative calculations (return true if success)
	bool initialize_lhs_and_rhs();

	// Calculate first derivatives at all data points and second derivatives at bounding points
	bool calculate_derivatives();

	// Store input data, initialize lhs and rhs and calculate derivatives
	bool initialize_interpolator(const double* x_source, int x_size, const double* y_source, int y_size);

public:
	SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source);
	SplineInterpolator(const double* x_source, int x_size, const double* y_source, int y_size);

	// Convenience function to check for correct derivative values
	void write_derivative_values(std::string file_name);
};

