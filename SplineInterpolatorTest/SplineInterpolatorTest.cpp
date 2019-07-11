// SplineInterpolatorTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>

#include "SplineInterpolator.h"

int main()
{
	// Specify data to fit spline functions to
	std::vector<double> x({0.0, 31.0, 59.0, 90.0, 120.0, 151.0, 181.0, 212.0, 243.0, 273.0, 304.0, 334.0, 365.25});
	std::vector<double> y({ 0.0, 93.0, 219.0, 265.5, 331.5, 424.5, 628.5, 1040.8, 1319.8, 1496.8, 1614.6, 1671.6, 1702.85 });

	// Create/initialize interpolators
	SplineInterpolator interp_default(x, y);
	SplineInterpolator interp_cspline(SplineInterpolator::GSLCubic, x, y);
	SplineInterpolator interp_steffen(SplineInterpolator::GSLSteffen, x, y);
	SplineInterpolator interp_akima(SplineInterpolator::GSLAkima, x, y);

	// Write to console interpolated value and derivative at a specified x value
	double x_in = 1072.0;
	std::cout << "Interpolated value (non-GSL) for " << x_in << " is " << interp_default.interpolate(x_in) << std::endl;
	std::cout << "Derivative is " << interp_default.interpolate_derivative(x_in) << std::endl;

	std::cout << "\nInterpolated value (cubic spline) for " << x_in << " is " << interp_cspline.interpolate(x_in) << std::endl;
	std::cout << "Derivative is " << interp_cspline.interpolate_derivative(x_in) << std::endl;

	std::cout << "\nInterpolated value (Steffen spline) for " << x_in << " is " << interp_steffen.interpolate(x_in) << std::endl;
	std::cout << "Derivative is " << interp_steffen.interpolate_derivative(x_in) << std::endl;

	std::cout << "\nInterpolated value (Akima spline) for " << x_in << " is " << interp_akima.interpolate(x_in) << std::endl;
	std::cout << "Derivative is " << interp_akima.interpolate_derivative(x_in) << std::endl;

	// Print vector and matrix data for default interpolator
	interp_default.write_vector_values("vector_values.txt");
	interp_default.write_matrix_values("matrix_values.txt");
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
