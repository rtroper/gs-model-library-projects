#include "pch.h"
#include "TimeSeries.h"

// Default constructor
TimeSeries::TimeSeries()
{
}

// Add a new time point and value to the time series
bool TimeSeries::addTimepoint(double time, double value)
{
	// Push new time and value
	times.push_back(time);
	values.push_back(value);

	// Return true if this point is reached
	// (For now, always assume success)
	return true;
}

// Get value at the specified index location (regardless of time)
double TimeSeries::getValueByIndex(int i)
{
	// Check that the index is in the valid range
	if (i >= 0 && size_t(i) < values.size())
	{
		return values[i];
	}
	// Otherwise, return a zero value
	else
	{
		return 0.0;
	}
}

// Get the size of the time series
int TimeSeries::size()
{
	// Assume that size of time and value vectors are equal
	// Note: addTimepoint function doesn't allow these to be unequal
	return values.size();
}

// Mean of all values in the time series
double TimeSeries::mean()
{
	return mean(values.size());
}

// Get the mean of a maximum of the first n values
double TimeSeries::mean(int n)
{
	double sum = 0.0;

	// Number of values to use to calculate the mean
	int number_of_values_to_use = min(values.size(), size_t(n));

	if ( number_of_values_to_use > 0 )
	{
		// Calculate sum of values and return the mean
		for (int i = 0; i < number_of_values_to_use; i++) sum += values[i];
		return sum / values.size();
	}
	else
	{
		// Return 0 if there is no data
		return 0.0;
	}
}

// Standard deviation of all values in the time series
double TimeSeries::stdev()
{
	return stdev(values.size());
}

// Get the standard deviation of a maximum of the first n values
double TimeSeries::stdev(int n)
{
	// Sum of squared deviations
	double sum_of_squared_deviations = 0.0;

	// Number of values to use to calculate the standard deviation
	int number_of_values_to_use = min(values.size(), size_t(n));

	if ( number_of_values_to_use > 0 )
	{
		// Get mean of values
		double mean = TimeSeries::mean(number_of_values_to_use);

		// Calculate sum of squared deviations from mean and return standard deviation
		for (int i = 0; i < number_of_values_to_use; i++) sum_of_squared_deviations += pow(values[i] - mean, 2.0);
		return sqrt(sum_of_squared_deviations / number_of_values_to_use);
	}
	else
	{
		// Return 0 if there is no data
		return 0.0;
	}
}

// Get correlation of values in 'this' TimeSeries with respect to TimeSeries ts
double TimeSeries::correlation(TimeSeries& ts)
{
	// Sum of the product of the normalized deviations from means
	double sum_product_deviations = 0.0;

	// Get the maximum number of values that can be used in calculating the correlation
	int number_of_values = min(TimeSeries::size(), ts.size());

	// At least two values are required to calculate correlation
	if ( number_of_values >= 2 )
	{
		// Get the means of the two time series
		double mean1 = mean(number_of_values);
		double mean2 = ts.mean(number_of_values);

		// Get the standard deviations of the two time series
		double stdev1 = stdev(number_of_values);
		double stdev2 = ts.stdev(number_of_values);

		// Calculate the sum of the products of normalized deviations from means and return correlation
		for (int i = 0; i < number_of_values; i++)
		{
			sum_product_deviations += ((values[i] - mean1) / stdev1 * (ts.getValueByIndex(i) - mean2) / stdev2);
		}
		return sum_product_deviations / (number_of_values);// -1.0);	// Without the -1, results are identical to Excel
	}
	else
	{
		// Return 0 (no correlation) if there are insufficient values
		return 0.0;
	}
}

// Get the slope linear regression coefficient with respect to another time series
double TimeSeries::linear_regression_slope(TimeSeries& ts)
{
	// Needed variables
	double sum1(0.0), sum2(0.0), sum_products(0.0), sum_squares(0.0);

	// Get the maximum number of values that can be used in calculating the correlation
	int number_of_values = min(TimeSeries::size(), ts.size());

	// At least two values are required to calculate slope
	if (number_of_values >= 2)
	{
		// Calculate necessary values
		for (int i = 0; i < number_of_values; i++)
		{
			sum1 += values[i];
			sum2 += ts.getValueByIndex(i);
			sum_products += (values[i] * ts.getValueByIndex(i));
			sum_squares += pow(values[i], 2.0);
		}
		return (number_of_values*sum_products - sum1*sum2) / (number_of_values*sum_squares - pow(sum1, 2.0));
	}
	else
	{
		// Return 0 if there are insufficient values
		return 0.0;
	}
}

// Get the intercept linear regression coefficient with respect to another time series
double TimeSeries::linear_regression_intercept(TimeSeries& ts)
{
	// Needed variables
	double sum1(0.0), sum2(0.0);

	// Get the maximum number of values that can be used in calculating the correlation
	int number_of_values = min(TimeSeries::size(), ts.size());

	// At least two values are required to calculate slope
	if (number_of_values >= 2)
	{
		// Calculate necessary values
		for (int i = 0; i < number_of_values; i++)
		{
			sum1 += values[i];
			sum2 += ts.getValueByIndex(i);
		}
		return (sum2 - TimeSeries::linear_regression_slope(ts) * sum1) / number_of_values;
	}
	else
	{
		// Return 0 if there are insufficient values
		return 0.0;
	}
}
