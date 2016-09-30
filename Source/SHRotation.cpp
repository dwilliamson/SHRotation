
#include "SHRotation.h"
#include "NRook.h"
#include "IvanicSHRotate.h"
#include "ChoiSHRotate.h"
#include <sdla/Maths.h>

using namespace sdla;


namespace
{
	struct TestLight : public cSphericalFunction
	{
		static double max(double a, double b)
		{
			return (a < b ? b : a);
		}

		double ToReal(const Sample& s) const
		{
			return (max(0, 5 * cos(s.theta) - 4) + max(0, -4 * sin(s.theta - PI) * cos(s.phi - 2.5) - 3));
		}

		cColour ToColour(const Sample& s) const
		{
			return (cColour((float)ToReal(s)));
		}
	};
}


cSHRotation::cSHRotation(const int width, const int height) :

	sdla::cApplication(width, height, false),
	m_NbSamples(100),
	m_NbBands(4),
	m_Samples(0),
	m_CoeffsAllocator(0),
	m_SHValues(0),
	m_InitialCoeffs(m_NbBands),
	m_Function(0),
	m_SHRotateMatrix(m_NbBands)

{
	// Generate the sphere samples
	m_Samples = NRook::GenerateSamples(m_NbSamples);

	// Create an allocator with enough space to store y(l,m) over each band for each sample
	m_CoeffsAllocator = new SphericalHarmonics::cCoeffs::Allocator(m_NbSamples * m_NbBands * m_NbBands);

	// Allocate a set of empty coeff lists
	m_SHValues = new SphericalHarmonics::cCoeffs[m_NbSamples];

	for (int i = 0; i < m_NbSamples; i++)
	{
		// De-ref current sample
		cSphericalFunction::Sample& s = m_Samples[i];

		// Create an SH co-efficient to store the required number of bands
		SphericalHarmonics::cCoeffs& sh = m_SHValues[i] = SphericalHarmonics::cCoeffs(m_NbBands, *m_CoeffsAllocator);

		// Record all y(l,m) for this sample
		for (int l = 0; l < m_NbBands; l++)
			for (int m = -l; m <= l; m++)
				sh(l, m) = SphericalHarmonics::y(l, m, s.theta, s.phi);
	}

	// Project the test light onto SH co-efficients
	SphericalHarmonics::Project(TestLight(), m_Samples, m_SHValues, m_NbSamples, m_InitialCoeffs);

	// Create the 3D representation of the reconstructed test light function
	m_LightList = glGenLists(2);
	m_Function = new cReconstructSH(m_InitialCoeffs);
	m_Function->Draw(m_LightList);
	m_Function->Draw(m_LightList + 1);

	glDisable(GL_LIGHTING);
}


cSHRotation::~cSHRotation(void)
{
	delete m_Function;

	glDeleteLists(m_LightList, 2);

	if (m_SHValues)
		delete [] m_SHValues;
	
	delete m_CoeffsAllocator;

	if (m_Samples)
		delete m_Samples;
}


void PrintSHRotation(const cSHRotateMatrix& m)
{
	printf("-------------------------------------------\n");
	for (int i = 0; i < m.GetNbBands(); i++)
	{
		for (int y = -i; y <= i; y++)
		{
			for (int x = -i; x <= i; x++)
			{
				printf("%4.4f ", m(i, x, y));
			}
			printf("\n");
		}
	}
}


void PrintMatrix(const cMatrix& m)
{
	printf("-------------------------------------------\n");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			printf("%4.4f ", m.e[i][j]);
		}
		printf("\n");
	}
}


int minus1pow(const int p)
{
	return ((p & 1) ? -1 : 1);
}


void RealToComplexMatrix(cSHRotateMatrix& F, cSHRotateMatrix& G)
{
	ASSERT(F.GetNbBands() == G.GetNbBands());

	const double oor2 = 0.7071067811865475244008443621052;

	for (int i = 0; i < F.GetNbBands(); i++)
	{
		// Eq. 8.5
		F(i, 0, 0) = 1;
		G(i, 0, 0) = 0;

		for (int x = -i; x <= i; x++)
		{
			for (int y = -i; y <= i; y++)
			{
				// Set to zero by default
				F(i, x, y) = 0;
				G(i, x, y) = 0;

				// Only the diagonals have any value
				if (abs(x) == abs(y))
				{
					// Eq. 8.6

					if (x >= 0 && y >= 0)
						F(i, x, y) = minus1pow(abs(x)) * oor2;

					else if (x < 0 && y < 0)
						G(i, x, y) = -oor2;

					else if (x < 0)
						G(i, x, y) = minus1pow(abs(x)) * oor2;

					else // y < 0
						F(i, x, y) = oor2;
				}
			}
		}
	}
}


void RealToComplexCoeffs(const SphericalHarmonics::cCoeffs& r, SphericalHarmonics::cCoeffs& cf, SphericalHarmonics::cCoeffs& cg)
{
	ASSERT(r.GetNbBands() == cf.GetNbBands());
	ASSERT(r.GetNbBands() == cg.GetNbBands());

	// Generate W
	cSHRotateMatrix F(r.GetNbBands()), G(r.GetNbBands());
	RealToComplexMatrix(F, G);

	// Generate the complex co-efficients
	F.Transform(r, cf);
	G.Transform(r, cg);
}


void PrintCoeffs(const SphericalHarmonics::cCoeffs& c)
{
	for (int i = 0; i < c.GetSize(); i++)
		printf("%5.5f ", c(i));
	printf("\n");
}


bool cSHRotation::ProcessFrame(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	Enable3D();
	glMatrixMode(GL_MODELVIEW);

	static float s = 0, t = 0;
	cMatrix r(cMatrix::RotateY(t) * cMatrix::RotateX(s));
	//PrintMatrix(rot.ToMatrix());

	// Draw OpenGL-rotated function
	{
		glLoadMatrix(cMatrix::Translate(cVector(-0.7f, 0, 3)));
		glMultMatrix(cMatrix::RotateY(t));
		glMultMatrix(cMatrix::RotateX(s));
		glCallList(m_LightList);
	}

	// If space is pressed, match the OpenGL rotation
	if (m_Keys[SDLK_SPACE])
	{
#if	DEBUG_CHOI

		// Rotate from m_InitialCoeffs into m_Function->coeffs
		ChoiSHRotate(r, m_SHRotateMatrix);
		PrintSHRotation(m_SHRotateMatrix);
		m_SHRotateMatrix.Transform(m_InitialCoeffs, m_Function->coeffs);

		// Make the rotated coefficients complex
		SphericalHarmonics::cCoeffs cf0(m_InitialCoeffs.GetNbBands()), cg0(m_InitialCoeffs.GetNbBands());
		RealToComplexCoeffs(m_Function->coeffs, cf0, cg0);

		// Now generate complex coeffs from the initial coeffs
		SphericalHarmonics::cCoeffs cf1(m_InitialCoeffs.GetNbBands()), cg1(m_InitialCoeffs.GetNbBands());
		SphericalHarmonics::cCoeffs dcf1(m_InitialCoeffs.GetNbBands()), dcg1(m_InitialCoeffs.GetNbBands());
		RealToComplexCoeffs(m_InitialCoeffs, cf1, cg1);

		// Rotate the initial complex coeffs by the complex matrices
		cSHRotateMatrix F(m_SHRotateMatrix.GetNbBands()), G(m_SHRotateMatrix.GetNbBands());
		ChoiSHRotate(r, F, G);
		F.Transform(cf1, dcf1);
		G.Transform(cg1, dcg1);

		// cf0 and dcf1 should be equal
		// cg0 and dcg1 should be equal
		// If not, the error is in the complex->real calcs in the choi rotator

		printf("----------------------------------------\n");
		PrintCoeffs(cf0);
		PrintCoeffs(dcf1);
		PrintCoeffs(cg0);
		PrintCoeffs(dcg1);

#else

		if (m_Keys['c'])
			ChoiSHRotate(r, m_SHRotateMatrix);
		else
			IvanicSHRotate(m_SHRotateMatrix, r);

		m_SHRotateMatrix.Transform(m_InitialCoeffs, m_Function->coeffs);

#endif

		m_Function->Draw(m_LightList + 1);
	}

	// Draw SH-rotated function
	glLoadMatrix(cMatrix::Translate(cVector(0.7f, 0, 3)));
	glCallList(m_LightList + 1);

	t += 0.1f;
	s += 0.05f;

	return (true);
}


void cSHRotation::BeforeSwitch(void)
{
}


void cSHRotation::AfterSwitch(void)
{
}