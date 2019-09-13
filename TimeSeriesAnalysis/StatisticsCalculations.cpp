#include "pch.h"
#include "StatisticsCalculations.h"

double StatisticsCalculations::mean(std::vector<double>& values)
{
	// Return the mean on all values in the vector
	return mean(values, values.size());
}

// Get mean of at most the first N values
double StatisticsCalculations::mean(std::vector<double>& values, size_t N)
{
	// Initialize count and sum of values
	size_t count = std::min(values.size(), N);
	double sum = 0.0;

	// Return 0.0 if there are no values
	if (count == 0) return 0.0;

	// Add up the values
	for (size_t i = 0; i < count; i++) sum += values[i];

	// Return the sum
	return sum / count;
}

double StatisticsCalculations::stdev(std::vector<double>& values)
{
	// Return the standard deviation on all values in the vector
	return stdev(values, values.size());
}

// Get standard deviation of at most the first N values
double StatisticsCalculations::stdev(std::vector<double>& values, size_t N)
{
	// Initialize count and sum of squared deviations from mean of the values
	size_t count = std::min(values.size(), N);
	double sum_squared_deviations = 0.0;

	// Return 0.0 if there are insufficient values (need at least 2 to avoid divide by zero)
	if (count < 2) return 0.0;

	// Get the mean
	double mean = StatisticsCalculations::mean(values);

	// Add up the squared deviations from mean
	for (size_t i = 0; i < count; i++) sum_squared_deviations += pow(values[i] - mean, 2.0);

	return sqrt(sum_squared_deviations / (count - 1));
}

// Get covariance (use number of values equal to min of sizes of each input vector)
double StatisticsCalculations::covar(std::vector<double>& v1, std::vector<double>& v2)
{
	// Initialize count and sum of product of deviations from means
	size_t count = std::min(v1.size(), v2.size());
	double sum_product_dev_from_means = 0.0;

	// Return 0.0 if there are insufficient values (need at least 2 to avoid divide by zero)
	if (count < 2) return 0.0;

	// Get the means
	double mean1 = StatisticsCalculations::mean(v1, count);
	double mean2 = StatisticsCalculations::mean(v2, count);

	// Add up product of deviations from the means
	for (size_t i = 0; i < count; i++) sum_product_dev_from_means += (v1[i] - mean1) * (v2[i] - mean2);

	return sum_product_dev_from_means / (count - 1);
}

// Get correlation (use number of values equal to min of sizes of each input vector)
double StatisticsCalculations::correlation(std::vector<double>& v1, std::vector<double>& v2)
{
	// Initialize count
	size_t count = std::min(v1.size(), v2.size());

	// Get standard deviations
	double sd1 = stdev(v1, count);
	double sd2 = stdev(v2, count);

	// If either standard deviation is zero, return 0.0
	if ((sd1 <= 0.0) || (sd2 <= 0.0)) return 0.0;

	// Get the covariance
	double covar = StatisticsCalculations::covar(v1, v2);

	return covar / (sd1 * sd2);
}

// Get linear regression coefficients, slope and intercept (use number of values equal to min of sizes of input vectors)
std::pair<double, double> StatisticsCalculations::linear_regression_coefs(std::vector<double>& v1, std::vector<double>& v2)
{
	// Initialize count
	size_t count = std::min(v1.size(), v2.size());

	// Initialize return values
	std::pair<double, double> coefs;
	coefs.first = coefs.second = 0.0;

	// Note, there is redundancy in the calculations below (e.g. stdev and covar require calculation of means)
	// But, for sake of keeping the API clean, there is no attempt to reduce duplication of calculations
	// If performance ever became an issue, this could be optimized (e.g. by having stdev return means in addition
	// to standard deviation so that the means could be used in other calculations)

	// Get standard deviation for v1
	double sd1 = StatisticsCalculations::stdev(v1, count);

	// If standard deviation is zero, return 0.0 values for coefficients
	if ( sd1 <= 0.0 ) return coefs;

	// Get means for both v1 and v2
	double mean1 = StatisticsCalculations::mean(v1, count);
	double mean2 = StatisticsCalculations::mean(v2, count);

	// Get covariance
	double covar = StatisticsCalculations::covar(v1, v2);

	// Calculate slope
	coefs.first = covar / pow(sd1, 2.0);

	// Calculate intercept
	coefs.second = mean2 - coefs.first * mean1;

	return coefs;
}


