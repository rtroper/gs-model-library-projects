#pragma once

#include <vector>

class TimeSeries
{
private:
	// Stores time values (assumed to be in seconds)
	std::vector<double> times;
	
	// Stores data values
	std::vector<double> values;

public:
	// Basic constructor requires no arguments
	TimeSeries();

	// Add a new time and value (return true if successfully added)
	bool addTimepoint(double time, double value);

	// Get value at the specified index location (regardless of time)
	double getValueByIndex(int i);

	// Get the size of the time series
	int size();

	// Get the mean of the data
	double mean();

	// Get the mean of a maximum of the first n values
	double mean(int n);

	// Get the standard deviation of the data
	double stdev();

	// Get the standard deviation of a maximum of the first n values
	double stdev(int n);

	// Get the correlation value with respect to another time series
	double correlation(TimeSeries &ts);

	// Get the slope linear regression coefficient with respect to another time series
	double linear_regression_slope(TimeSeries& ts);

	// Get the intercept linear regression coefficient with respect to another time series
	double linear_regression_intercept(TimeSeries& ts);
};

