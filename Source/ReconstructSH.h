
#pragma once


#include <sdla/Colour.h>
#include "SphericalHarmonics.h"


struct cReconstructSH : public cSphericalFunction
{
	// Construct from SH co-efficients
	cReconstructSH(const SphericalHarmonics::cCoeffs& _coeffs) : coeffs(_coeffs) { }

	double ToReal(const Sample& s) const
	{
		int nb_bands = coeffs.GetNbBands();

		// Calculate a coeff-weighted sum of all the SH bases
		double value = 0;
		for (int l = 0, i = 0; l < nb_bands; l++)
			for (int m = -l; m <= l; m++, i++)
				value += coeffs(i) * SphericalHarmonics::y(l, m, s.theta, s.phi);

		return (value);
	}

	sdla::cColour ToColour(const Sample& s) const
	{
		return (sdla::cColour((float)ToReal(s)));
	}

	// Copy of the SH list
	SphericalHarmonics::cCoeffs	coeffs;
};
