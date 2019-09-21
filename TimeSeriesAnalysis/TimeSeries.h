#pragma once

#include <vector>

class TimeSeries
{
private:
	// Stores time values (assumed to be in seconds)
	std::vector<double> times;
	
	// Stores data values
	std::vector<double> values;

	// Missing value code (values <= to this will be treated as 'missing')
	double missing;

	// Get the minimum and maximum times in the series
	std::pair<double, double> getTimeBounds();

	// Get aligned non-missing values (based on time) between 'this' and another time series 
	std::pair<std::vector<double>, std::vector<double>> getAlignedValues(TimeSeries &ts);

public:
	// Basic constructor requires no arguments
	TimeSeries();

	// Set the missing value code
	void setMissingValueCode(double value) { missing = value; };

	// Set the missing value code
	double getMissingValueCode() { return missing; };

	// Add a new time and value (return true if successfully added)
	bool addTimepoint(double time, double value);

	// Get value at the specified index location (regardless of time)
	double getValueByIndex(int i);

	// Get value that corresponds to the specified time (linearly interpolates, if necessary)
	double getValueByTime(double time);

	// Get the size of the time series (assume that size of time and value vectors are equal)
	// Note: addTimepoint function doesn't allow these to be unequal
	size_t size() { return values.size(); };

	// Mean of all non-missing values in the time series
	double mean();

	// Standard deviation of all non-missing values in the time series
	double stdev();

	// Get the correlation value with respect to another time series
	double correlation(TimeSeries &ts);

	// Get linear regression coefficients (slope and intercept) with respect to another time series
	std::pair<double, double> linear_regression_coefs(TimeSeries& ts);

	// Get autocorrelation for a specified time shift
	double autocorrelation(double time_shift);
};

