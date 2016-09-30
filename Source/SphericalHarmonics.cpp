
#include "SphericalHarmonics.h"
#include <cmath>
#include <sdla/Exception.h>
#include <sdla/Maths.h>

using namespace sdla;


struct FactorialGen
{
	FactorialGen(const int count) :
		factorials(0),
		nb_factorials(count)
	{
		factorials = new int[nb_factorials];

		// Build a factorial LUT
		factorials[0] = 1;
		for (int i = 1; i < nb_factorials; i++)
			factorials[i] = i * factorials[i - 1];
	}

	~FactorialGen(void)
	{
		delete [] factorials;
	}

	int operator () (const int i) const
	{
		// Check input for debug builds
#ifdef	_DEBUG
		ASSERT(i >= 0 && i < nb_factorials);
#endif
		return (factorials[i]);
	}

	int*	factorials;
	int		nb_factorials;
};


double SphericalHarmonics::P(const int l, const int m, const double x)
{
	// Start with P(0,0) at 1
	double pmm = 1;

	// First calculate P(m,m) since that is the only rule that requires results
	// from previous bands
	// No need to check for m>0 since SH function always gives positive m

	// Precalculate (1 - x^2)^0.5
	double somx2 = sqrt(1 - x * x);

	// This calculates P(m,m). There are three terms in rule 2 that are being iteratively multiplied:
	//
	// 0: -1^m
	// 1: (2m-1)!!
	// 2: (1-x^2)^(m/2)
	//
	// Term 2 has been partly precalculated and the iterative multiplication by itself m times
	// completes the term.
	// The result of 2m-1 is always odd so the double factorial calculation multiplies every odd
	// number below 2m-1 together. So, term 3 is calculated using the 'fact' variable.
	double fact = 1;
	for (int i = 1; i <= m; i++)
	{
		pmm *= -1 * fact * somx2;
		fact += 2;
	}

	// No need to go any further, rule 2 is satisfied
	if (l == m)
		return (pmm);

	// Since m<l in all remaining cases, all that is left is to raise the band until the required
	// l is found

	// Rule 3, use result of P(m,m) to calculate P(m,m+1)
	double pmmp1 = x * (2 * m + 1) * pmm;

	// Is rule 3 satisfied?
	if (l == m + 1)
		return (pmmp1);

	// Finally, use rule 1 to calculate any remaining cases
	double pll = 0;
	for (int ll = m + 2; ll <= l; ll++)
	{
		// Use result of two previous bands
		pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);

		// Shift the previous two bands up
		pmm = pmmp1;
		pmmp1 = pll;
	}

	return (pll);
}


// Re-normalisation constant for SH function
namespace {
double K(const int l, const int m)
{
	// Generate lut the first time this function is called
	// Note that this method doesn't interact well with explicit init/shutdown memory managers
	// Additionally, 32-bits is not enough to store any higher bands (7 max) !!!
	static FactorialGen factorial(14);

	// Note that |m| is not used here as the SH function always passes positive m
	return (sqrt(((2 * l + 1) * factorial(l - m)) / (4 * PI * factorial(l + m))));
}
}	// End namespace


double SphericalHarmonics::y(const int l, const int m, const double theta, const double phi)
{
	if (m == 0)
		return (K(l, 0) * P(l, 0, cos(theta)));

	if (m > 0)
		return (sqrt(2.0) * K(l, m) * cos(m * phi) * P(l, m, cos(theta)));

	// m < 0, m is negated in call to K
	return (sqrt(2.0) * K(l, -m) * sin(-m * phi) * P(l, -m, cos(theta)));
}


void SphericalHarmonics::Project(const cSphericalFunction& func, const cSphericalFunction::Sample* samples, const cCoeffs* coeffs, const int nb_samples, SphericalHarmonics::cCoeffs& dest)
{
	int nb_coeffs = dest.GetSize();

	// Check input
	ASSERT(samples && coeffs);
	ASSERT(nb_samples >= 0);
	ASSERT(coeffs[0].GetSize() == nb_coeffs);

	// Clear out sums
	for (int i = 0; i < nb_coeffs; i++)
		dest(i) = 0;

	for (int i = 0; i < nb_samples; i++)
	{
		// Take the sample at this point on the sphere
		const cSphericalFunction::Sample& s = samples[i];
		const cCoeffs& c = coeffs[i];
		double sample = func.ToReal(s);

		// Sum the projection of this sample onto each SH basis
		if (sample)
			for (int j = 0; j < nb_coeffs; j++)
				dest(j) += sample * c(j);
	}

	// Divide each coefficient by the number of samples and multiply by weights
	for (int i = 0; i < nb_coeffs; i++)
		dest(i) = dest(i) * ((4 * PI) / nb_samples);
}