
#pragma once


#include <new>
#include <sdla/Exception.h>
#include "SphericalHarmonics.h"


class cSHRotateMatrix
{
public:
	class SubMatrix
	{
	public:
		// Construct that needs previously allocated element array
		SubMatrix(double* elements, const int width)  :

			m_Elements(elements),
			m_Width(width),
			m_Size(width * width),
			m_Shift((width - 1) / 2)

		{
		}


		// Element accessors with (+,-) lookup
		double& operator () (const int m, const int n)
		{
			return (m_Elements[Index(m, n)]);
		}
		double operator () (const int m, const int n) const
		{
			return (m_Elements[Index(m, n)]);
		}


		// Linear element accessors
		double& operator () (const int i)
		{
			return (m_Elements[Index(i)]);
		}
		double operator () (const int i) const
		{
			return (m_Elements[Index(i)]);
		}


		int	GetWidth(void) const
		{
			return (m_Width);
		}

		int	GetShift(void) const
		{
			return (m_Shift);
		}


	private:
		int Index(const int m, const int n) const
		{
			// Check bounds in debug build
			#ifdef	_DEBUG
			ASSERT(m >= -m_Shift && m <= m_Shift);
			ASSERT(n >= -m_Shift && n <= m_Shift);
			#endif

			return ((m + m_Shift) * m_Width + (n + m_Shift));
		}

		int Index(const int i) const
		{
			// Check bounds in debug build
			#ifdef	_DEBUG
			ASSERT(i >= 0 && i < m_Size);
			#endif

			return (i);
		}

		// Element array
		double*	m_Elements;

		// Width of the sub-matrix
		int		m_Width;

		// Total size of the sub-matrix in multiples of doubles
		int		m_Size;

		// Value to shift incoming matrix indices by
		int		m_Shift;
	};


	cSHRotateMatrix(const int nb_bands) :

		m_NbBands(nb_bands),
		m_Matrices(0),
		m_Elements(0)

	{
		// Allocate the sub-matrix array
		m_Matrices = (SubMatrix*)operator new (m_NbBands * sizeof(SubMatrix));

		int size = 0;
		for (int i = 0; i < m_NbBands; i++)
			size += (i * 2 + 1) * (i * 2 + 1);

		// Allocate the entire element array
		m_Elements = new double[size];

		// Construct each sub-matrix
		for (int i = 0, j = 0; i < m_NbBands; i++)
		{
			int w = i * 2 + 1;
			new (&m_Matrices[i]) SubMatrix(m_Elements + j, w);
			j += (w * w);
		}

		// The first 1x1 sub-matrix is always 1 and so doesn't need to be stored in the rotation matrix
		// It's stored simply to make the indexing logic as easy as possible
		m_Elements[0] = 1;
	}


	~cSHRotateMatrix(void)
	{
		if (m_Elements)
			delete [] m_Elements;
		if (m_Matrices)
			delete [] m_Matrices;
	}


	// Element accessors which take the band index and (+,-) index relative to the centre
	// of the sub-matrix
	double&	operator () (const int l, const int m, const int n)
	{
		return (m_Matrices[l](m, n));
	}
	double operator () (const int l, const int m, const int n) const
	{
		return (m_Matrices[l](m, n));
	}


	int	GetNbBands(void) const
	{
		return (m_NbBands);
	}


	void Transform(const SphericalHarmonics::cCoeffs& source, SphericalHarmonics::cCoeffs& dest) const
	{
		// Check number of bands match
		int nb_bands = source.GetNbBands();
		ASSERT(nb_bands == dest.GetNbBands());

		// Band 0 is always multiplied by 1 so it stays untouched
		dest(0) = source(0);

		// Loop through each band
		for (int l = 1; l < nb_bands; l++)
		{
			SubMatrix& M = m_Matrices[l];

			// Calculate band offset into coeff list
			int band_offset = l * (l + 1);

			// Now through each argument of the destination coeff-list
			for (int mo = -l; mo <= l; mo++)
			{
				// Clear destination
				double& d = dest(band_offset + mo);
				d = 0;

				// Pre-calculate of mo's lookup in the sub-matrix, plus mi's shift
				int p_mo = (mo + M.GetShift()) * M.GetWidth() + M.GetShift();

				// Multiply-add with each argument of the source coeff-list
				for (int mi = -l; mi <= l; mi++)
					d += source(band_offset + mi) * M(p_mo + mi);
			}
		}
	}

private:
	// Number of bands stored
	int		m_NbBands;

	// Sub-matrix for each band
	SubMatrix*	m_Matrices;

	// Element array for each sub-matrix
	double*	m_Elements;
};