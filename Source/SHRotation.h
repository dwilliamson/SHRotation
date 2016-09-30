
#pragma once


#include <sdla/Application.h>
#include "SphericalHarmonics.h"
#include "ReconstructSH.h"
#include "SHRotateMatrix.h"


class cSHRotation : public sdla::cApplication
{
public:
	// Constructor/destructor
	cSHRotation(const int width, const int height);
	~cSHRotation(void);

private:
	// cApplication implementations
	bool	ProcessFrame(void);
	void	BeforeSwitch(void);
	void	AfterSwitch(void);

	// Number of samples on the sphere
	int		m_NbSamples;

	// Number of SH bands
	int		m_NbBands;

	// Uniformly distributed samples on the sphere
	cSphericalFunction::Sample*	m_Samples;

	// Allocator for y(l,m) values
	SphericalHarmonics::cCoeffs::Allocator*	m_CoeffsAllocator;

	// y(l,m) stored for each sample over m_NbBands
	SphericalHarmonics::cCoeffs*	m_SHValues;

	// 3D representation of test light
	int		m_LightList;

	// Unrotated SH co-efficients
	SphericalHarmonics::cCoeffs	m_InitialCoeffs;

	// Function containing coeff rotation destination
	cReconstructSH*	m_Function;

	cSHRotateMatrix	m_SHRotateMatrix;
};