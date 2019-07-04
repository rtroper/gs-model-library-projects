#pragma once

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

//struct gsl_vector;

class SplineInterpolator
{
private:
	size_t size_data;					// Number of data points to which to fit spline
	gsl_vector* x, * y;				// Independent and dependent variable data to which to fit spline functions
	gsl_matrix* lhs;				// Left-hand-side matrix in equation to solve for derivatives
	gsl_vector* rhs;				// Right-hand-side vector in equation to solve for derivatives
	gsl_vector* a, * b, * c, * d;	// Coefficients in spline functions
	gsl_vector* D, * h;				// Derivatives (D) and x-variable intervals (h)

	// Store input data in member variables (return true if success)
	bool store_input_data(const double* x_source, int x_size, const double* y_source, int y_size);

	// Initialize left-hand-side matrix and right-hand-side vector for derivative calculations (return true if success)
	bool initialize_lhs_and_rhs();

public:
	SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source);
	SplineInterpolator(const double* x_source, int x_size, const double* y_source, int y_size);

};

