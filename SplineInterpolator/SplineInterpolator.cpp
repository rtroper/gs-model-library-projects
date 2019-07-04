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
	for (int i = 0; i < size_data; i++)
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

SplineInterpolator::SplineInterpolator(const std::vector<double>& x_source, const std::vector<double>& y_source)
{
	// Get sizes of source data
	int x_size = static_cast<int>(x_source.size());
	int y_size = static_cast<int>(y_source.size());

	// Simply call constructor that uses C-style arrays
	SplineInterpolator(x_source.data(), x_size, y_source.data(), y_size);
}

SplineInterpolator::SplineInterpolator(const double* x_source, int x_size, const double* y_source, int y_size)
{
	// Store data and calculate intervals h
	store_input_data(x_source, x_size, y_source, y_size);

	// Initialize left-hand-side matrix and right-hand-side vector used to calculate derivatives
	initialize_lhs_and_rhs();
}
