#include "pch.h"
#include "TimeSeries.h"
#include "Constants.h"
#include "StatisticsCalculations.h"

std::pair<double, double> TimeSeries::getTimeBounds()
{
	// Initialize return values
	std::pair<double, double> bounds(0.0, 0.0);

	// Get the number of time points in the series
	size_t size = TimeSeries::size();

	// Only attempt to get bounds if size >= 1;
	if (size > 0)
	{
		bounds.first = times[0];
		bounds.second = times[size-1];
	}

	return bounds;
}

std::pair<std::vector<double>, std::vector<double>> TimeSeries::getAlignedValues(TimeSeries& ts)
{
	// Initialize two vectors to store values that will be used in correlation calculation
	std::vector<double> v1;
	std::vector<double> v2;
	
	// Get time bounds of the two series
	std::pair<double, double> ts1_bounds = getTimeBounds();
	std::pair<double, double> ts2_bounds = ts.getTimeBounds();

	// Get min and max range over which to calculate correlation, which can only be calculated on overlapping portions
	double minTime = std::max(ts1_bounds.first, ts2_bounds.first);
	double maxTime = std::min(ts1_bounds.second, ts2_bounds.second);

	// Get non-missing values within the min and max time range (i.e. the overlapping portion of the two series)
	for (size_t i = 0; i < times.size(); i++)
	{
		// If already beyond max time, break out of loop
		if (times[i] > maxTime) break;

		// If value at current time is missing, skip to next iteration in the loop
		if (values[i] < missing + TSConstants::SMALL_CONSTANT) continue;

		// If current time is within the time range, get values for both series
		if (times[i] >= minTime)
		{
			// Get value (interpolated, if necessary) for ts
			double ts_interp_value = ts.getValueByTime(times[i]);

			// If the value from ts is missing, skip to next iteration in the loop
			if (ts_interp_value < missing + TSConstants::SMALL_CONSTANT) continue;

			// Store values to be returned
			v1.push_back(values[i]);
			v2.push_back(ts_interp_value);
		}
	}

	// Initialize pair and store and return vectors
	std::pair<std::vector<double>, std::vector<double>> aligned_values;
	aligned_values.first = v1;
	aligned_values.second = v2;

	return aligned_values;
}

// Default constructor
TimeSeries::TimeSeries() : missing(-999.0) // Default missing value code
{
}

// Add a new time point and value to the time series
bool TimeSeries::addTimepoint(double time, double value)
{
	// Push new time and value
	times.push_back(time);
	values.push_back(value);

	// Return true if this point is reached (For now, always assume success)
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

// Get value that corresponds to the specified time (linearly interpolates, if necessary)
double TimeSeries::getValueByTime(double time)
{
	// Initialize return value to the missing value
	double value = missing;
	
	// Return missing if time value is less than min or more than max time value
	if (time < times[0] || time > times[size() - 1]) return missing;

	// Find index of highest time value < the specified time
	size_t idx_lower = 0;
	while (idx_lower < (size() - 1) && time > times[idx_lower + 1])
		idx_lower += 1;

	// Return missing if value at lower or upper index is missing
	double adj_missing = missing + TSConstants::SMALL_CONSTANT;
	if ((times[idx_lower] < adj_missing) || (times[idx_lower + 1] < adj_missing)) return missing;

	// Get linearly interpolated value
	double fractional_idx = (time - times[idx_lower]) / (times[idx_lower + 1] - times[idx_lower]);
	value = values[idx_lower] + (values[idx_lower + 1] - values[idx_lower]) * fractional_idx;

	return value;
}

// Mean of all non-missing values in the time series
double TimeSeries::mean()
{
	// Generate a vector with only non-missing value
	std::vector<double> v;
	for (size_t i = 0; i < values.size(); i++)
	{
		if (values[i] > missing + TSConstants::SMALL_CONSTANT) v.push_back(values[i]);
	}

	return StatisticsCalculations::mean(v);
}

// Standard deviation of all non-missing values in the time series
double TimeSeries::stdev()
{
	// Generate a vector with only non-missing value
	std::vector<double> v;
	for (size_t i = 0; i < values.size(); i++)
	{
		if (values[i] > missing + TSConstants::SMALL_CONSTANT) v.push_back(values[i]);
	}

	return StatisticsCalculations::stdev(v);
}

// Get correlation of values with respect to TimeSeries ts
double TimeSeries::correlation(TimeSeries& ts)
{
	// Get aligned values for the two series
	std::pair<std::vector<double>, std::vector<double>> aligned_values;
	aligned_values = getAlignedValues(ts);

	// Return the correlation
	return StatisticsCalculations::correlation(aligned_values.first, aligned_values.second);
}

// Get linear regression coefficients (slope and intercept) with respect to TimeSeries ts
std::pair<double, double> TimeSeries::linear_regression_coefs(TimeSeries& ts)
{
	// Get aligned values for the two series
	std::pair<std::vector<double>, std::vector<double>> aligned_values;
	aligned_values = getAlignedValues(ts);

	// Return the linear regression coefficients (slope and intercept)
	return StatisticsCalculations::linear_regression_coefs(aligned_values.first, aligned_values.second);
}

double TimeSeries::autocorrelation(double time_shift)
{
	// Get time bounds
	std::pair<double, double> bounds = getTimeBounds();

	// Initialize vectors to store non-missing values and shifted values (not necessary to store time values)
	std::vector<double> values1;
	std::vector<double> values2;

	// Store time-shifted values
	size_t time_idx = 0;
	while ((time_idx < times.size()) && (times[time_idx] + time_shift <= bounds.second))
	{
		// Only get time-shifted value if current value is not missing
		if (values[time_idx] > missing + TSConstants::SMALL_CONSTANT)
		{
			double time_shifted_value = getValueByTime(times[time_idx] + time_shift);

			// Only store value and time-shifted value if the latter is not missing
			if (time_shifted_value > missing + TSConstants::SMALL_CONSTANT)
			{
				values1.push_back(values[time_idx]);
				values2.push_back(time_shifted_value);
			}
		}

		// Increment time index to avoid infinite loop
		time_idx++;
	}
	
	return StatisticsCalculations::correlation(values1, values2);
}
