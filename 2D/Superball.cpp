#include "Superball.h"

CPolyhedron::CPolyhedron(void)
{
	int i;
	for (i = 0; i < 3; i++)	pA[i] = 1.0;
	center.SetZero();
	e[0].Set(1.0, 0.0);
	e[1].Set(0.0, 1.0);
	rout = 0.0;
	rin = 0.0;
	area = 0.0;
}
CPolyhedron::CPolyhedron(double ppara[3], double rrr[2])
{
	int i;
	for (i = 0; i < 3; i++)	pA[i] = ppara[i];
	center.SetZero();
	e[0].Set(1.0, 0.0);
	e[1].Set(0.0, 1.0);

	rout = rrr[0];
	rin = rrr[1];
	double a1 = tgamma(0.5 / ppara[2]);
	double a2 = tgamma(1.0 / ppara[2]);
	area = ppara[0] * ppara[1] * a1 * a1 / (a2 * ppara[2]);
}
CPolyhedron::~CPolyhedron(void)
{
}
double CPolyhedron::sgnp(double x)
{
	double t;
	if (x == 0) t = 0.0;
	else if (x > 0) t = 1.0;
	else t = -1.0;
	return t;
}

void CPolyhedron::Copyfrom(CPolyhedron *pp)
{
	int i;
	for (i = 0; i < 3; i++)	pA[i] = pp->pA[i];
	center = pp->center;
	for (i = 0; i < 2; i++) e[i] = pp->e[i];
	rout = pp->rout;
	rin = pp->rin;
}

void CPolyhedron::ChangeSize(double ratio)
{
	for (int i = 0; i < 2; i++) pA[i] *= ratio;
	rin *= ratio; rout *= ratio;
	area = area * ratio * ratio;
}

void CPolyhedron::Update_EtoL(double L_inv[3][3])//欧拉坐标 to 拉格朗日相对坐标
{
	int i;
	for (i = 0; i < 3; i++)
		center_L[i] = L_inv[i][0] * center[0] + L_inv[i][1] * center[1] + L_inv[i][2] * center[2];
}
void CPolyhedron::Update_LtoE(double lambda[3][3])
{
	int i;
	for (i = 0; i < 3; i++)
		center[i] = lambda[i][0] * center_L[0] + lambda[i][1] * center_L[1] + lambda[i][2] * center_L[2];
}
void CPolyhedron::Jump(CVector TV)
{
	center += TV;
}

