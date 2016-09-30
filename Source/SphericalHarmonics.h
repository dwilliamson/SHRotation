
#pragma once


#include <sdla/StackAllocator.h>
#include "SphericalFunction.h"


namespace SphericalHarmonics
{
	class cCoeffs
	{
	public:
		// Allocator for co-efficients
		typedef sdla::tStackAllocator<double>	Allocator;

		cCoeffs(void) : m_Coeffs(0), m_NbBands(0), m_NbCoeffs(0) { }

		cCoeffs(const int nb_bands) :
			m_Coeffs(0),
			m_NbBands(nb_bands),
			m_NbCoeffs(m_NbBands * m_NbBands),
			m_FromAllocator(false)
		{
			// Allocate the coeff list
			m_Coeffs = new double[m_NbCoeffs];
		}

		cCoeffs(const int nb_bands, Allocator& alloc) :
			m_Coeffs(0),
			m_NbCoeffs(nb_bands * nb_bands),
			m_FromAllocator(true)
		{
			// Allocate the coeff list
			m_Coeffs = new (alloc) double[m_NbCoeffs];
		}

		~cCoeffs(void)
		{
			// Stack allocator cleans itself up
			if (!m_FromAllocator)
				delete [] m_Coeffs;
		}

		cCoeffs(const cCoeffs& other) :
			m_Coeffs(0),
			m_NbBands(other.m_NbBands),
			m_NbCoeffs(other.m_NbCoeffs),
			m_FromAllocator(other.m_FromAllocator)
		{
			// Just copy the pointer in this case
			if (m_FromAllocator)
				m_Coeffs = other.m_Coeffs;

			else
			{
				// Make a copy when not using the stack allocator
				m_Coeffs = new double[m_NbCoeffs];
				for (int i = 0; i < m_NbCoeffs; i++)
					m_Coeffs[i] = other.m_Coeffs[i];
			}
		}

		// Const/non-const accessors for band/arg
		double operator () (const int l, const int m) const
		{
			return (m_Coeffs[Check(l * (l + 1) + m)]);
		}
		double& operator () (const int l, const int m)
		{
			return (m_Coeffs[Check(l * (l + 1) + m)]);
		}

		// Const/non-const accessors by index
		double operator () (const int i) const
		{
			return (m_Coeffs[Check(i)]);
		}
		double& operator () (const int i)
		{
			return (m_Coeffs[Check(i)]);
		}

		int	GetNbBands(void) const
		{
			return (m_NbBands);
		}

		int GetSize(void) const
		{
			return (m_NbCoeffs);
		}

	private:
		inline int Check(const int index) const
		{
			// Check bounds in debug build
			#ifdef	_DEBUG
			ASSERT(index >= 0 && index < m_NbCoeffs);
			#endif

			return (index);
		}

		// List of SH co-efficients
		double*	m_Coeffs;

		// Number of bands used
		int		m_NbBands;

		// Number of coefficients (=bands^2)
		int		m_NbCoeffs;

		// Allocated from a special-purpose allocator?
		bool	m_FromAllocator;
	};

	// Evaluate Associated Legendre Polynomial
		double	P(const int l, const int m, const double x);

	// Evaluate Real Spherical Harmonic
	double	y(const int l, const int m, const double theta, const double phi);

	// Project a spherical function onto the spherical harmonic bases
	//
	// a. Spherical Function to project
	// b. List of sample locations on the sphere
	// c. List of y(l,m,theta,phi) for each sample
	// d. Number of samples on the sphere
	// e. Destination projection
	void Project(const cSphericalFunction& func,
				 const cSphericalFunction::Sample* samples,
				 const cCoeffs* coeffs,
				 const int nb_samples, 
				 cCoeffs& dest);
}