#include "pch.h"
#include "SplineInterpolator.h"

bool SplineInterpolator::store_input_data(const double* x_source, int x_size, const double* y_source, int y_size)
{
	// Default return value is true
	bool success = true;

	// Store data up to the maximum available (for x and y)
	size_data = static_cast<size_t>(min(x_size, y_size));

	// Initialize x and y vectors to store input data
	x = gsl_vector_calloc(size_data);
	y = gsl_vector_calloc(size_data);

	// Initialize h to store x intervals
	h = gsl_vector_calloc(size_data - 1);

	// Store source data in x and y and calculate h
	for (size_t i = 0; i < size_data; i++)
	{
		gsl_vector_set(x, i, x_source[i]);
		gsl_vector_set(y, i, y_source[i]);

		if (i < size_data - 1)
		{
			// Store difference between x[i+1] and x[i]
			gsl_vector_set(h, i, x_source[i + 1] - x_source[i]);
		}
	}

	return success;
}

// Left-hand-side matrix (lhs) and right-hand-side vector (rhs) must not be initialized until
// store_input_data has been called since y and h must first be initialized
bool SplineInterpolator::initialize_lhs_and_rhs()
{
	// Default return value is true
	bool success = true;

	// Cast data size to correct type to initialize lhs and rhs
	const size_t size = static_cast<size_t>(size_data);

	// Initialize lhs and rhs with 0 values
	lhs = gsl_matrix_calloc(size + 2, size + 2);
	rhs = gsl_vector_calloc(size + 2);

	// NOTE: lhs and rhs are used to solve for first derivatives D0, D1, ..., DN at data points
	// (x0, y0), (x1, y1), ... , (xN, yN), where N = size-1, and second derivatives D0' and DN' at 
	// (x0, y0) and (xN, yN). Thus, lhs, rhs, and the resulting solution are of size 'size + 2',
	// where 'size' is the number of data points through which the spline function passes. In the
	// resulting solution vector, the first 'size' values are D0, D1, ... , DN and the last 2 values
	// are D0' and DN'.

	// Initialize the first 'size' rows of lhs and rhs
	double h1 = 0.0, h2 = 0.0, y1 = 0.0, y2 = 0.0, y3 = 0.0;
	for (size_t i = 0; i < size; i++)
	{
		if (i == 0)
		{
			// Get h and y values needed for calculations
			h1 = gsl_vector_get(h, 0);
			y1 = gsl_vector_get(y, 0);
			y2 = gsl_vector_get(y, 1);

			// Set lhs values
			gsl_matrix_set(lhs, 0,    0,          2.0 * h1 );
			gsl_matrix_set(lhs, 0,    1,                h1 );
			gsl_matrix_set(lhs, 0, size, 0.5 * pow(h1, 2.0));

			// Set rhs value
			gsl_vector_set(rhs, 0, 3.0 * (y2 - y1));
		}
		else if (i == size - 1)
		{
			// Get h and y values needed for calculations
			h1 = gsl_vector_get(h, size - 2);
			y1 = gsl_vector_get(y, size - 2);
			y2 = gsl_vector_get(y, size - 1);

			// Set lhs values
			gsl_matrix_set(lhs, size - 1, size - 1,           2.0 * h1 );
			gsl_matrix_set(lhs, size - 1, size - 2,                 h1 );
			gsl_matrix_set(lhs, size - 1, size + 1, -0.5 * pow(h1, 2.0));

			// Set rhs value
			gsl_vector_set(rhs, size - 1, 3.0 * (y2 - y1));
		}
		else
		{
			// Get h and y values needed for calculations
			h1 = gsl_vector_get(h, i - 1);
			h2 = gsl_vector_get(h, i    );
			y1 = gsl_vector_get(y, i - 1);
			y2 = gsl_vector_get(y, i    );
			y3 = gsl_vector_get(y, i + 1);

			// Set lhs values
			gsl_matrix_set(lhs, i, i - 1,                                h1  );
			gsl_matrix_set(lhs, i, i    , 2.0 * h1 * (1 + pow(h1 / h2, 2.0)) );
			gsl_matrix_set(lhs, i, i + 1,       h1 *      pow(h1 / h2, 2.0)  );

			// Set rhs value
			gsl_vector_set(rhs, i, 3.0 * (y3 - y2) * pow(h1 / h2, 2.0) + 3.0 * (y2 - y1));
		}
	}

	// Set values in last 2 rows of lhs to specify boundary conditions. The default is for
	// boundary first derivatives to be equal to each other and boundary second derivatives to be
	// equal to each other. Note that for these conditions, the last 2 rows of rhs are zero and,
	// therefore, do not need to be modified. TO-DO: Create a separate function that can be used
	// to set different boundary conditions.
	gsl_matrix_set(lhs,     size,        0,  1.0);
	gsl_matrix_set(lhs,     size, size - 1, -1.0);
	gsl_matrix_set(lhs, size + 1,     size,  1.0);
	gsl_matrix_set(lhs, size + 1, size + 1, -1.0);

	return success;
}

bool SplineInterpolator::calculate_derivatives()
{
	// Default return value is true
	bool success = true;

	// Define integer for error returns
	int status = 0;

	// Create and initialize LU, a matrix to store LU decomposition when solving for derivatives
	gsl_matrix* LU = gsl_matrix_alloc(size_data + 2, size_data + 2);
	status = gsl_matrix_memcpy(LU, lhs);
	if (status) return false;	// Add code for proper error handling/messaging

	// Create permutation vector to store permutation
	gsl_permutation* p = gsl_permutation_calloc(size_data + 2);

	// Create pointer to int for third argument to LU decomposition function
	int* signum = new int(0);

	// Get the LU decomposition
	status = gsl_linalg_LU_decomp(LU, p, signum);
	if (status) return false;	// Add code for proper error handling/messaging

	// Initialize D (to store results of LU solve)
	D = gsl_vector_calloc(size_data + 2);

	// Solve for D
	const gsl_vector* rhs_const = rhs;
	status = gsl_linalg_LU_solve(LU, p, rhs, D);
	if (status) return false;	// Add code for proper error handling/messaging

	return success;
}

bool SplineInterpolator::calculate_spline_parameters()
{
	// Default return value is true
	bool success = true;

	// Initialize a, b, c, and d
	a = gsl_vector_calloc(size_data - 1);
	b = gsl_vector_calloc(size_data - 1);
	c = gsl_vector_calloc(size_data - 1);
	d = gsl_vector_calloc(size_data - 1);

	// Load calculated values into a, b, c, and d
	double h1 = 0.0, y1 = 0.0, y2 = 0.0, D1 = 0.0, D2 = 0.0;
	for (size_t i = 0; i < size_data - 1; i++)
	{
		// Get h, y and D values needed for calculations
		h1 = gsl_vector_get(h,     i);
		y1 = gsl_vector_get(y,     i);
		y2 = gsl_vector_get(y, i + 1);
		D1 = gsl_vector_get(D,     i);
		D2 = gsl_vector_get(D, i + 1);

		// Calculate a, b, c, and d values
		gsl_vector_set(a, i,                                    y1 );
		gsl_vector_set(b, i,                               h1 * D1 );
		gsl_vector_set(c, i, 3.0 * (y2 - y1) - h1 * (2.0 * D1 + D2));
		gsl_vector_set(d, i,       2.0 * (y1 - y2) + h1 * (D1 + D2));
	}

	return success;
}

bool SplineInterpolator::initialize_interpolator(const double* x_source, int x_size, const double* y_source, int y_size)
{
	// Default return value is true
	bool success = true;

	// Store data and calculate intervals h
	success = store_input_data(x_source, x_size, y_source, y_size);
	if (!success) return success;	// Add code for proper error handling/messaging

	// Initialize left-hand-side matrix and right-hand-side vector used to calculate derivatives
	success = initialize_lhs_and_rhs();
	if (!success) return success;	// Add code for proper error handling/messaging

	// Calculate derivatives at all data points (to be used in calculating a, b, c, and d parameters)
	success = calculate_derivatives();
	if (!success) return success;	// Add code for proper error handling/messaging

	// Calculate spline parameters a, b, c, and d
	success = calculate_spline_parameters();
	if (!success) return success;	// Add code for proper error handling/messaging

	return success;
}

std::pair<int, double> SplineInterpolator::get_i_and_t_values(double x_value)
{
	// Initialize return value
	std::pair<int, double> i_and_t = std::pair<int, double>();

	// First mod the input x value with the maximum x so that values
	// outside x range are wrapped
	double x_adjusted = std::fmod(x_value, gsl_vector_get(x, size_data - 1));

	// Loop over x to find i such that x[i] < x_value < x[i+1]
	int i = 0;
	for (size_t idx = 1; idx < size_data; idx++)
	{
		if (x_adjusted > gsl_vector_get(x, idx)) i++;
		else break;
	}

	// Calculate t value i.e. (x_value - x[i]) / (x[i+1] - x[i])
	double t = (x_adjusted - gsl_vector_get(x, i)) / gsl_vector_get(h, i);

	// Store i and t values to return
	i_and_t.first = i;
	i_and_t.second = t;

	return i_and_t;
}

SplineInterpolator::SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source)
{
	// Get sizes of source data
	int x_size = static_cast<int>(x_source.size());
	int y_size = static_cast<int>(y_source.size());

	// Call initialization function that uses C-style arrays
	initialize_interpolator(x_source.data(), x_size, y_source.data(), y_size);
}

SplineInterpolator::SplineInterpolator(const double* x_source, int x_size, const double* y_source, int y_size)
{
	// Call initialization function that uses C-style arrays
	initialize_interpolator(x_source, x_size, y_source, y_size);
}

double SplineInterpolator::interpolate(double x_value)
{
	double y_interp = 0.0;

	// Get i and t values that correspond to the input x value
	std::pair<int, double> i_and_t = get_i_and_t_values(x_value);
	int i = i_and_t.first;
	double t = i_and_t.second;

	// Get parameter values needed for interpoloation
	double a_val = gsl_vector_get(a, i);
	double b_val = gsl_vector_get(b, i);
	double c_val = gsl_vector_get(c, i);
	double d_val = gsl_vector_get(d, i);
	
	// Calculate the interpolated y value
	y_interp = a_val + b_val * t + c_val * std::pow(t, 2.0) + d_val * std::pow(t, 3.0);

	return y_interp;
}

double SplineInterpolator::interpolate_derivative(double x_value)
{
	double deriv_interp = 0.0;

	// Get i and t values that correspond to the input x value
	std::pair<int, double> i_and_t = get_i_and_t_values(x_value);
	int i = i_and_t.first;
	double t = i_and_t.second;

	// Get parameter values needed for calculation
	double a_val = gsl_vector_get(a, i);
	double b_val = gsl_vector_get(b, i);
	double c_val = gsl_vector_get(c, i);
	double d_val = gsl_vector_get(d, i);

	// Get x interval value needed for calculation
	double h_val = gsl_vector_get(h, i);

	// Calculate the interpolated y value
	deriv_interp = (b_val + 2.0 * c_val * t + 3.0 * d_val * std::pow(t, 2.0)) / h_val;

	return deriv_interp;
}

// Generic convenience function to print values of all vector-type member data
void SplineInterpolator::write_vector_values(std::string file_name)
{
	// Open output file
	FILE* file;
	errno_t err = fopen_s(&file, file_name.data(), "w");
	
	// If file is not null, write all vector data
	if (file)
	{
		// Write x
		fprintf_s(file, "Input data independent variable (x) values:\n\n");
		if (x)
		{
			gsl_vector_fprintf(file, x, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write y
		fprintf_s(file, "Input data dependent variable (y) values:\n\n");
		if (y)
		{
			gsl_vector_fprintf(file, y, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write h
		fprintf_s(file, "Independent variable interval (h) values:\n\n");
		if (h)
		{
			gsl_vector_fprintf(file, h, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write rhs
		fprintf_s(file, "Right-hand-side (rhs) of equation to solve for derivatives:\n\n");
		if (rhs)
		{
			gsl_vector_fprintf(file, rhs, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write D
		fprintf_s(file, "Derivative (D) values:\n\n");
		if (D)
		{
			gsl_vector_fprintf(file, D, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write a
		fprintf_s(file, "Spline parameter 'a' values:\n\n");
		if (a)
		{
			gsl_vector_fprintf(file, a, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write b
		fprintf_s(file, "Spline parameter 'b' values:\n\n");
		if (b)
		{
			gsl_vector_fprintf(file, b, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write c
		fprintf_s(file, "Spline parameter 'c' values:\n\n");
		if (c)
		{
			gsl_vector_fprintf(file, c, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");

		// Write d
		fprintf_s(file, "Spline parameter 'd' values:\n\n");
		if (d)
		{
			gsl_vector_fprintf(file, d, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");
	}

	// Close the file
	if (file) fclose(file);
}

// Generic convenience function to print values of all matrix-type member data
void SplineInterpolator::write_matrix_values(std::string file_name)
{
	// Open output file
	FILE* file;
	errno_t err = fopen_s(&file, file_name.data(), "w");
	
	// If file is not null, write all matrix data
	if (file)
	{
		// Write lhs
		fprintf_s(file, "Left-hand-side (lhs) of equation to solve for derivatives:\n\n");
		if (lhs)
		{
			gsl_matrix_fprintf(file, lhs, "%f");
			fprintf_s(file, "\n");
		}
		else fprintf_s(file, "NULL DATA!!\n\n");
	}

	// Close the file
	if (file) fclose(file);
}
