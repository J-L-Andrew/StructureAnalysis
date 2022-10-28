#include "Vector.h"
int CVector::Normalize(double dv[3])
{
	double len = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
	if (len < 1.0E-15) return 0;
	if (len > 0.99999999 && len < 1.00000001) return 1;
	len = sqrt(len);
	dv[0] /= len;
	dv[1] /= len;
	dv[2] /= len;
	return 1;
}