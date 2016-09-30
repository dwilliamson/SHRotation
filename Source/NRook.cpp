
#include "NRook.h"
#include <sdla/Maths.h>
#include <cstdlib>
#include <cstring>

using namespace sdla;


// Method of uniformly distributing samples over a sphere
// From ACM Journal of Graphics Tools (http://www.acm.org/jgt)
// Thanks for Robin Green for suggestion


namespace
{
	double drand48(void)
	{
		return ((double)rand() / RAND_MAX);
	}

	int lrand48(void)
	{
		return (rand());
	}

	void MultiStageNRooks(const int size, int* cells)
	{
		if (size == 1)
			return;

		int size1 = size >> 1;
		int size2 = size >> 1;

		if (size & 1)
		{
			if (drand48() > 0.5)
				size1++;
			else
				size2++;
		}

		int* upper_cells = new int[size1];
		int* lower_cells = new int[size2];

		int i = 0, j = 0;
		for ( ; i < size - 1; i += 2, j++)
		{
			if (lrand48() & 1)
			{
				upper_cells[j] = cells[i];
				lower_cells[j] = cells[i + 1];
			}
			else
			{
				upper_cells[j] = cells[i + 1];
				lower_cells[j] = cells[i];
			}
		}

		if (size1 != size2)
		{
			if (size1 > size2)
				upper_cells[j] = cells[i];
			else
				lower_cells[j] = cells[i];
		}

		MultiStageNRooks(size1, upper_cells);
		memcpy(cells, upper_cells, size1 * sizeof(int));
		delete [] upper_cells;

		MultiStageNRooks(size2, lower_cells);
		memcpy(cells + size1, lower_cells, size2 * sizeof(int));
		delete [] lower_cells;
	}
}


cSphericalFunction::Sample* NRook::GenerateSamples(const int nb_samples)
{
	// Generate nrook cells
	int* cells = new int[nb_samples];
	for (int i = 0; i < nb_samples; i++)
		cells[i] = i;
	MultiStageNRooks(nb_samples, cells);

	// Allocate the samples list
	cSphericalFunction::Sample* samples = new cSphericalFunction::Sample[nb_samples];

	for (int i = 0; i < nb_samples; i++)
	{
		// Generate uniform random sample
		double x = (i + drand48()) / nb_samples;
		double y = (cells[i] + drand48()) / nb_samples;

		// Generate spherical/cartesian co-ordinates
		samples[i] = cSphericalFunction::Sample(2 * acos(sqrt(1 - x)), 2 * PI * y);
	}

	delete [] cells;

	return (samples);
}