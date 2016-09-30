
#pragma once


#include "SHRotateMatrix.h"


namespace sdla
{
	struct cMatrix;
}


void ChoiSHRotate(const sdla::cMatrix& rotation, cSHRotateMatrix& dest);

void ChoiSHRotate(const sdla::cMatrix& rotation, cSHRotateMatrix& F, cSHRotateMatrix& G);