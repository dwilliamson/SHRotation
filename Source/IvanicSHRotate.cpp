
// Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion
// Joseph Ivanic and Klaus Ruedenberg
// J. Phys. Chem. 1996, 100, 6342-5347
// 
// Additions and Corrections (to previous paper)
// Joseph Ivanic and Klaus Ruedenberg
// J. Phys. Chem. A, 1998, Vol. 102, No. 45, 9099


#include "IvanicSHRotate.h"
#include "SHRotateMatrix.h"
#include <sdla/Matrix.h>
#include <sdla/Exception.h>
#include <cmath>
#include <cstring>

using namespace sdla;


namespace
{
	inline double delta(const int m, const int n)
	{
		// Kronecker Delta
		return (m == n ? 1 : 0);
	}


	void uvw(const int l, const int m, const int n, double& u, double& v, double& w)
	{
		// Pre-calculate simple reusable terms
		double d = delta(m, 0);
		int abs_m = abs(m);

		// Only calculate the required denominator once
		double denom;
		if (abs(n) == l)
			denom = (2 * l) * (2 * l - 1);

		else
			denom = (l + n) * (l - n);

		// Now just calculate the scalars
		u = sqrt((l + m) * (l - m) / denom);
		v = 0.5f * sqrt((1 + d) * (l + abs_m - 1) * (l + abs_m) / denom) * (1 - 2 * d);
		w = -0.5f * sqrt((l - abs_m - 1) * (l - abs_m) / denom) * (1 - d);
	}


	double P(const int i, const int l, const int a, const int b, const cSHRotateMatrix& M)
	{
		// Rather than passing the permuted rotation matrix around grab it directly from the first
		// rotation band which is never modified
		double ri1 = M(1, i, 1);
		double rim1 = M(1, i, -1);
		double ri0 = M(1, i, 0);

		if (b == -l)
			return (ri1 * M(l - 1, a, -l + 1) + rim1 * M(l - 1, a, l - 1));

		else if (b == l)
			return (ri1 * M(l - 1, a, l - 1) - rim1 * M(l - 1, a, -l + 1));

		else // |b|<l
			return (ri0 * M(l - 1, a, b));
	}


	double U(const int l, const int m, const int n, const cSHRotateMatrix& M)
	{
		// All cases fall correctly through here
		return (P(0, l, m, n, M));
	}


	double V(const int l, const int m, const int n, const cSHRotateMatrix& M)
	{
		if (m == 0)
		{
			double p0 = P(1, l, 1, n, M);
			double p1 = P(-1, l, -1, n, M);
			return (p0 + p1);
		}

		else if (m > 0)
		{
			double d = delta(m, 1);
			double p0 = P(1, l, m - 1, n, M);
			double p1 = P(-1, l, -m + 1, n, M);
			return (p0 * sqrt(1 + d) - p1 * (1 - d));
		}

		else // m < 0
		{
			double d = delta(m, -1);
			double p0 = P(1, l, m + 1, n, M);
			double p1 = P(-1, l, -m - 1, n, M);
			return (p0 * (1 - d) + p1 * sqrt(1 + d));
		}
	}


	double W(const int l, const int m, const int n, const cSHRotateMatrix& M)
	{
		if (m == 0)
		{
			// Never gets called as kd=0
			ASSERT(false);
			return (0);
		}

		else if (m > 0)
		{
			double p0 = P(1, l, m + 1, n, M);
			double p1 = P(-1, l, -m - 1, n, M);
			return (p0 + p1);
		}

		else // m < 0
		{
			double p0 = P(1, l, m - 1, n, M);
			double p1 = P(-1, l, -m + 1, n, M);
			return (p0 - p1);
		}
	}


	double M(const int l, const int m, const int n, const cSHRotateMatrix& M)
	{
		// First get the scalars
		double u, v, w;
		uvw(l, m, n, u, v, w);

		// Scale by their functions
		if (u)
			u *= U(l, m, n, M);
		if (v)
			v *= V(l, m, n, M);
		if (w)
			w *= W(l, m, n, M);

		return (u + v + w);
	}
}


void IvanicSHRotate(cSHRotateMatrix& shrm, const cMatrix& rotation)
{
	// To order required in paper
	cMatrix r = cMatrix::Transpose(rotation);

	// The first band is rotated by a permutation of the original matrix
	shrm(1, -1, -1) = r.e[1][1];
	shrm(1, -1,  0) = r.e[1][2];
	shrm(1, -1, +1) = r.e[1][0];
	shrm(1,  0, -1) = r.e[2][1];
	shrm(1,  0,  0) = r.e[2][2];
	shrm(1,  0, +1) = r.e[2][0];
	shrm(1, +1, -1) = r.e[0][1];
	shrm(1, +1,  0) = r.e[0][2];
	shrm(1, +1, +1) = r.e[0][0];

	// Calculate each block of the rotation matrix for each subsequent band
	for (int band = 2; band < shrm.GetNbBands(); band++)
	{
		for (int m = -band; m <= band; m++)
			for (int n = -band; n <= band; n++)
				shrm(band, m, n) = M(band, m, n, shrm);
	}
}