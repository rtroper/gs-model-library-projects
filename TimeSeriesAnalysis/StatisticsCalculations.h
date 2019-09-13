#pragma once
#include <vector>

// Static functions for calculating statistics
class StatisticsCalculations
{
public:
	// Functions for calculating mean of a vector of values
	static double mean(std::vector<double> &values);
	static double mean(std::vector<double>& values, size_t N);	// Get mean of at most the first N values

	// Function for calculating standard deviation of a vector of values
	static double stdev(std::vector<double> &values);
	static double stdev(std::vector<double>& values, size_t N);	// Get standard deviation of at most the first N values

	// Get covariance (use number of values equal to min of sizes of each input vector)
	static double covar(std::vector<double> &v1, std::vector<double> &v2);

	// Get correlation (use number of values equal to min of sizes of each input vector)
	static double correlation(std::vector<double>& v1, std::vector<double>& v2);

	// Get linear regression coefficients, slope and intercept (use number of values equal to min of sizes of input vectors)
	static std::pair<double, double> linear_regression_coefs(std::vector<double>& v1, std::vector<double>& v2);
};

