
#pragma once


#include <cmath>
#include <sdla/Vector.h>


namespace sdla
{
	struct cColour;
};


struct cSphericalFunction
{
	struct Sample
	{
		// Default constructor
		Sample(void) { };

		// Construct from spherical co-ordinates
		Sample(const double u, const double v) :
			theta(u),
			phi(v)
		{
			x = sin(theta) * cos(phi);
			y = sin(theta) * sin(phi);
			z = cos(theta);
		}

		// Spherical co-ordinates
		double	theta, phi;

		// Cartesian co-ordinates
		double	x, y, z;
	};

	// Single valued evaluation
	virtual double	ToReal(const Sample& s) const = 0;

	// Triple valued colour evaluation
	virtual sdla::cColour	ToColour(const Sample& s) const = 0;

	// Map from spherical co-ordinates onto the sphere using the unsigned sample as the radius
	sdla::cVector ToVector(const Sample& s) const
	{
		double r = fabs(ToReal(s));
		return (sdla::cVector((float)(r * s.x), (float)(r * s.y), float(r * s.z)));
	}

	// Helper for constructing sample type
	sdla::cVector ToVector(const double u, const double v) const
	{
		return (ToVector(Sample(u, v)));
	}

	// OpenGL draw method
	void	Draw(void) const;

	// OpenGL draw method - draws into a display list
	void	Draw(const int name) const;
};