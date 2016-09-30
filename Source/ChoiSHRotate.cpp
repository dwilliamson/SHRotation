
// Rapid and stable determination of rotation matrices between spherical harmonics by direct recursion
// Cheol Ho Choi, Joseph Ivanic, Mark S. Gordon, and Klaus Ruedenberg
// J. Chem. Phys, 1999, Vol. 111, No. 19


#include "ChoiSHRotate.h"
#include <sdla/Matrix.h>

using namespace sdla;


namespace
{
	enum
	{
		x = 0,
		y = 1,
		z = 2
	};


	struct cChoiSHRotate
	{
		cChoiSHRotate(const cMatrix& R_, cSHRotateMatrix& F_, cSHRotateMatrix& G_) :
			R(R_),
			F(F_),
			G(G_)
		{
			// 1 / sqrt(2)
			const double s = 0.70710678118655;

			// Construct D(1), the first order Wigner D matrix, directly from the input rotation

			// Eq. 5.4
			F(1, -1, -1) = 0.5f * (R.e[y][y] + R.e[x][x]);
			F(1, -1,  0) = s * R.e[x][z];
			F(1, -1,  1) = 0.5f * (R.e[y][y] - R.e[x][x]);
			F(1,  0, -1) = s * R.e[z][x];
			F(1,  0,  0) = R.e[z][z];
			F(1,  0,  1) = -F(1, 0, -1);
			F(1,  1, -1) = F(1, -1, 1);
			F(1,  1,  0) = -F(1, -1, 0);
			F(1,  1,  1) = F(1, -1, -1);

			// Eq. 5.5
			G(1, -1, -1) = 0.5f * (R.e[y][x] - R.e[x][y]);
			G(1, -1,  0) = s * R.e[y][z];
			G(1, -1,  1) = -0.5f * (R.e[y][x] + R.e[x][y]);
			G(1,  0, -1) = -s * R.e[z][y];
			G(1,  0,  0) = 0;
			G(1,  0,  1) = G(1, 0, -1);
			G(1,  1, -1) = -G(1, 0, -1);
			G(1,  1,  0) = G(1, -1, 0);
			G(1,  1,  1) = -G(1, -1, -1);

			// Recursively construct higher orders using D(1) and D(l - 1)
			for (int l = 2; l < F.GetNbBands(); l++)
			{
				// Wrong when m=2,n=1
				for (int m = -l; m <= l; m++)
				{
					// NOTE: Eq. 5.8 can be used to simplify the calculation of these
					// Only need to calculate the upper triangle

					// Calculate for n=-l
					FG_1(l, m, -l);

					// Calculate for n=l
					FG_2(l, m, l);

					// Calculate the rest when n!=abs(l)
					for (int n = -l + 1; n < l; n++)
						FG_0(l, m, n);
				}
			}
		}

		void ComplexToReal(cSHRotateMatrix& Rn)
		{
			// Now calculate the real rotation matrix R(n) from F(n),G(n) of D(n)

			// First note that R(1) is merely a permutation of R
			Rn(1, -1, -1) = R.e[y][y];
			Rn(1, -1,  0) = R.e[y][z];
			Rn(1, -1, +1) = R.e[y][x];
			Rn(1,  0, -1) = R.e[z][y];
			Rn(1,  0,  0) = R.e[z][z];
			Rn(1,  0, +1) = R.e[z][x];
			Rn(1, +1, -1) = R.e[x][y];
			Rn(1, +1,  0) = R.e[x][z];
			Rn(1, +1, +1) = R.e[x][x];

			// Then calculate the higher orders
			for (int l = 2; l < Rn.GetNbBands(); l++)
			{
				for (int m = -l; m <= l; m++)
				{
					for (int n = -l; n <= l; n++)
						CalculateReal(l, m, n, Rn);
				}

				/*for (int m = 1; m <= l; m++)
				{
					for (int n = 1; n <= m; n++)
					{
						F(l, -m, -n) = pow(-1, m + n) * F(l, m, n);
						G(l, -m, -n) = -pow(-1, m + n) * G(l, m, n);
					}
				}*/
			}
		}

		double a(const int l, const int m, const int n) const
		{
			// Eq. 6.4
			if (l == abs(m))
				return (0);

			// Eq. 6.2
			return (sqrt(double((l + m) * (l - m)) / double((l + n) * (l - n))));
		}

		double b(const int l, const int m, const int n) const
		{
			// Eq. 6.5
			if (m == -l || m == -l + 1)
				return (0);

			// Eq. 6.3
			return (sqrt(double((l + m) * (l + m - 1)) / double(2 * (l + n) * (l - n))));
		}

		double c(const int l, const int m, const int n) const
		{
			// Eq. 6.12
			if (l == abs(m))
				return (0);

			// Eq. 6.10
			return (sqrt(double(2 * (l + m) * (l - m)) / double((l + n) * (l + n - 1))));
		}

		double d(const int l, const int m, const int n) const
		{
			// Eq. 6.13
			if (m == -l || m == -l + 1)
				return (0);

			// Eq. 6.11
			return (sqrt(double((l + m) * (l + m - 1)) / double((l + n) * (l + n - 1))));
		}

		double H(const int l, const int m, const int n, const int i, const int j) const
		{
			// Eq. 7.1
			return (F(1, i, j) * F(l - 1, m, n) - G(1, i, j) * G(l - 1, m, n));
		}

		double K(const int l, const int m, const int n, const int i, const int j) const
		{
			// Eq. 7.2
			return (F(1, i, j) * G(l - 1, m, n) + G(1, i, j) * F(l - 1, m, n));
		}

		// For |n|!=l
		void FG_0(const int l, const int m, const int n)
		{
			// Coverage for equations 7.3 & 7.4
			// Calculation split so as not to enter H or K when a or b are not defined

			double& f = F(l, m, n) = 0;
			double& g = G(l, m, n) = 0;

			if (double k = a(l, m, n))
			{
				f += k * H(l, m, n, 0, 0);
				g += k * K(l, m, n, 0, 0);
			}

			if (double k = b(l, m, n))
			{
				f += k * H(l, m - 1, n, 1, 0);
				g += k * K(l, m - 1, n, 1, 0);
			}

			if (double k = b(l, -m, n))
			{
				f += k * H(l, m + 1, n, -1, 0);
				g += k * K(l, m + 1, n, -1, 0);
			}
		}

		// For n=-l
		void FG_1(const int l, const int m, const int n)
		{
			// Coverage for equations 7.5 & 7.6
			// Calculation split so as not to enter H or K when c or d are not defined
			
			double& f = F(l, m, n) = 0;
			double& g = G(l, m, n) = 0;

			if (double k = c(l, m, -n))
			{
				f += k * H(l, m, n + 1, 0, -1);
				g += k * K(l, m, n + 1, 0, -1);
			}

			if (double k = d(l, m, -n))
			{
				f += k * H(l, m - 1, n + 1, 1, -1);
				g += k * K(l, m - 1, n + 1, 1, -1);
			}

			if (double k = d(l, -m, -n))
			{
				f += k * H(l, m + 1, n + 1, -1, -1);
				g += k * K(l, m + 1, n + 1, -1, -1);
			}
		}

		// For n=l
		void FG_2(const int l, const int m, const int n)
		{
			// Coverage for equations 7.7 & 7.8
			// Calculation split so as not to enter H or K when c or d are not defined

			double& f = F(l, m, n) = 0;
			double& g = G(l, m, n) = 0;

			if (double k = c(l, m, n))
			{
				f += k * H(l, m, n - 1, 0, 1);
				g += k * K(l, m, n - 1, 0, 1);
			}

			if (double k = d(l, m, n))
			{
				f += k * H(l, m - 1, n - 1, 1, 1);
				g += k * K(l, m - 1, n - 1, 1, 1);
			}

			if (double k = d(l, -m, n))
			{
				f += k * H(l, m + 1, n - 1, -1, 1);
				g += k * K(l, m + 1, n - 1, -1, 1);
			}
		}

		int delta(const int m, const int n) const
		{
			// Kronecker delta
			return (m == n);
		}

		double alpha(const int m) const
		{
			//return (pow(sqrt(1 + delta(m, 0)) * -1, (abs(m) + m) / 2));
			return (sqrt(1.0 + delta(m, 0)) * pow(-1.0, (abs(m) + m) / 2));
			//return (sqrt(pow(-1 * (1 + delta(m, 0)), abs(m) + m)));
		}

		// Notes from Choi:
		//
		// eq (8.10) seems to have some typos.
		// I think the alpha should read as alpha(m)= sqrt(1 + delta(m, 0)) * (-1) ^ ((abs(m) + m) / 2)
		// For all the eq. (8.11) to (8.14), I think 2 is missing in front of R, just like the 2 in eq. (8.8) and (8.9).
 

		double beta(const int m) const
		{
			// sign(s)
			int s = (m > 0) - (m < 0);

			// Eq. 8.10 (b)
			return (s * (1 - delta(m, 0)) * alpha(m));
		}

		void CalculateReal(const int l, const int m, const int n, cSHRotateMatrix& Rn)
		{
			double& rn = Rn(l, m, n);

			// Use of different expressions for the 4 quadrants in the current band rotation matrix
			// This piece of code can be well optimised by taking advantage of the fact that
			// pairs of quadrants share common terms (8.11 and 8.12 only differ by add/sub)

			if (m >= 0 && n >= 0)
			{
				// Eq. 8.11
				rn = alpha(m) * alpha(n) * F(l, m, n) - beta(m) * beta(-n) * F(l, m, -n);

				//if (m == 0 && n == 0)
				//	rn /= 2;
			}

			// CORRECT
			else if (m < 0 && n < 0)
			{
				// Eq. 8.12
				rn = alpha(m) * alpha(n) * F(l, m, n) + beta(m) * beta(-n) * F(l, m, -n);
			}

			else if (n < 0)
			{
				// Eq. 8.13
				// NOTE: Added negative to first term
				double g0 = G(l, m, n);
				double g1 = G(l, m, -n);
				double a[] =
				{
					alpha(m), alpha(-m),
					alpha(n), alpha(-n),
					beta(m), beta(-m),
					beta(n), beta(-n)
				};
				rn = -alpha(m) * beta(n) * G(l, m, n) - beta(m) * alpha(-n) * G(l, m, -n);
			}

			// CORRECT
			else // m < 0
			{
				// Eq. 8.14
				// NOTE: No negative on first term
				rn = alpha(m) * beta(n) * G(l, m, n) - beta(m) * alpha(-n) * G(l, m, -n);
			}
		}

		// 3x3 source rotation matrix
		const cMatrix&	R;

		// Real and imaginary parts of D(n)
		cSHRotateMatrix&	F;
		cSHRotateMatrix&	G;
	};


	// Complex rotation matrix is: D=F+Gi
	//
	// Recurrence relations for D(l) in terms of D(l-1) and D(1) are:
	//
	// (-l+1) <= n <= (l-1)		[R.0]
	// F(l,m,n) = a(l,m,n).H(l,m,n,0,0) + b(l,m,n).H(l,m-1,n,+,0) + b(l,-m,n).H(l,m+1,n,-,0)	Eq. 7.3
	// G(l,m,n) = a(l,m,n).K(l,m,n,0,0) + b(l,m,n).K(l,m-1,n,+,0) + b(l,-m,n).K(l,m+1,n,-,0)	Eq. 7.4
	//
	// -l <= n <= (l-2)			[R.1]
	// F(l,m,n) = c(l,m,-n).H(l,m,n+1,0,-) + d(l,m,-n).H(l,m-1,n+1,+,-) + d(l,-m,-n).H(l,m+1,n+1,-,-)	Eq. 7.5
	// G(l,m,n) = c(l,m,-n).K(l,m,n+1,0,-) + d(l,m,-n).K(l,m-1,n+1,+,-) + d(l,-m,-n).K(l,m+1,n+1,-,-)	Eq. 7.6
	//
	// (-l+2) <= n <= l			[R.2]
	// F(l,m,n) = c(l,m,n).H(l,m,n-1,0,+) + d(l,m,n).H(l,m-1,n-1,+,+) + d(l,-m,n).H(l,m+1,n-1,-,+)	Eq. 7.7
	// G(l,m,n) = c(l,m,n).K(l,m,n-1,0,+) + d(l,m,n).K(l,m-1,n-1,+,+) + d(l,-m,n).K(l,m+1,n-1,-,+)	Eq. 7.8
	//
	// where:
	//
	// H(l,m,n,i,j) = F(i,j).F(l-1,m,n) - G(i,j).G(l-1,m,n)		Eq. 7.1
	// K(l,m,n,i,j) = F(i,j).G(l-1,m,n) + G(i,j).F(l-1,m,n)		Eq. 7.2
	//
	// a(l,m,n) = A(m,l) / A(n,l)
	// b(l,m,n) = B(m,l) / sqrt(2).A(n,l)
	// c(l,m,n) = sqrt(2).A(m,l) / B(n,l)
	// d(l,m,n) = B(m,l) / B(n,l)
	//
	// and:
	//
	// A(m,l) = sqrt[(l + m)(l - m) / (2l + 1)(2l - 1)]
	// B(m,l) = sqrt[(l + m)(l + m - 1) / (2l + 1)(2l - 1)]
	//
	// Note the overlap of the ranges for which the equations are valid. This is reduced to using R.0 for when
	// n!=abs(l), using R.1 when n=-l, and using R.2 when n=l.
	//
	// The constants a,b,c,d are all of the form (a'.b'/c'.d') where a',b',c',d are all square roots
	// of integers not larger than 2L+1, where is L is the highest order you'll use. Choi's CLM table is the
	// result of this observation but I've yet to figure out how an entry in CLM relates to a',b',c',d'.
}


void ChoiSHRotate(const cMatrix& rotation, cSHRotateMatrix& F, cSHRotateMatrix& G)
{
	cChoiSHRotate(cMatrix::Transpose(rotation), F, G);
}


void ChoiSHRotate(const cMatrix& r, cSHRotateMatrix& dest)
{
	cSHRotateMatrix F(dest.GetNbBands()), G(dest.GetNbBands());
	cChoiSHRotate shr(cMatrix::Transpose(r), F, G);
	shr.ComplexToReal(dest);
}