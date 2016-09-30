
#include "SphericalFunction.h"
#include <sdla/Maths.h>
#include <sdla/SDLWrap.h>

using namespace sdla;


namespace
{
	cVector CalculateNormal(const cVector& a, const cVector& b, const cVector&c )
	{
		// Get edges
		cVector e0 = b - a;
		cVector e1 = c - a;

		// Get perpendicular
		return ((e0 ^ e1).Normalise());
	}
}


void cSphericalFunction::Draw(void) const
{
	const int	NB_SAMPLES = 50;

	// Steps along the sphere
	double du = PI / NB_SAMPLES;
	double dv = 2 * PI / NB_SAMPLES;

	glBegin(GL_QUADS);
	for (double u = 0; u < PI; u += du)
	{
		for (double v = 0; v < 2 * PI; v += dv)
		{
			// Sample 4 points in the neighbourhood
			cVector p[4] =
			{
				ToVector(u, v),
				ToVector(u + du, v),
				ToVector(u + du, v + dv),
				ToVector(u, v + dv)
			};

			double ddu = du / 10;
			double ddv = dv / 10;

			// Sample the normals as even smaller patches in the neighbourhood
			cVector n[4] =
			{
				CalculateNormal(p[0], ToVector(u + ddu, v), ToVector(u, v + ddv)),
				CalculateNormal(p[1], ToVector(u + du + ddu, v), ToVector(u + du, v + ddv)),
				CalculateNormal(p[2], ToVector(u + du + ddu, v + dv), ToVector(u + du, v + dv + ddv)),
				CalculateNormal(p[3], ToVector(u + ddu, v + dv), ToVector(u, v + dv + ddv))
			};

			// Positive SH samples are green, negative are red
			float c[3] = { 0, 0, 0 };
			if (ToReal(Sample(u, v)) < 0)
				c[0] = 1;
			else
				c[1] = 1;
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);

			// Draw the quad
			glNormal3fv(n[0]);
			glVertex3fv(p[0]);
			glNormal3fv(n[1]);
			glVertex3fv(p[1]);
			glNormal3fv(n[2]);
			glVertex3fv(p[2]);
			glNormal3fv(n[3]);
			glVertex3fv(p[3]);
		}
	}
	glEnd();
}


void cSphericalFunction::Draw(const int name) const
{
	// Draw into a display list
	glNewList(name, GL_COMPILE);
	Draw();
	glEndList();
}