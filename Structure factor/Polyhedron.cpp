#include "Polyhedron.h"

CPolyhedron::CPolyhedron(void)
{
	int i;
	for (i = 0; i < 3; i++)	pA[i] = 1.0;
	for (i = 0; i < 2; i++) pP[i] = 1.0;	
	center.SetZero();
	e[0].Set(1.0, 0.0, 0.0);
	e[1].Set(0.0, 1.0, 0.0);
	e[2].Set(0.0, 0.0, 1.0);	
	rout = 0.0;
	rin = 0.0;
	volume = 0.0;
}
CPolyhedron::CPolyhedron(double ppara[5],double rrr[2])
{
	int i;
	for (i = 0; i < 3; i++)	pA[i] = ppara[i];
	for (i = 0; i < 2; i++) pP[i] = ppara[i+3];	
	center.SetZero();	
	e[0].Set(1.0, 0.0, 0.0);
	e[1].Set(0.0, 1.0, 0.0);
	e[2].Set(0.0, 0.0, 1.0);

	rout = rrr[0];
	rin = rrr[1];
	double a1, a2, b1, b2, b3;
	a1 = tgamma(0.5 / pP[0]);
	a2 = tgamma(1.0 / pP[0]);
	b1 = tgamma(0.5 / pP[1]);
	b2 = tgamma(1.0 / pP[1]);
	b3 = tgamma(1.5 / pP[1]);
	volume = 2.0*pA[0] * pA[1] * pA[2] / 3.0 / pP[0] / pP[1] * a1*a1*b1*b2 / a2 / b3;
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
void CPolyhedron::Getroutrin()
{
	int i, j, k;
	CVector ver[26];
	double vrr, vrrmax, vrrmin;
	double x, y, z, v, u, cos1, sin1, cos2, sin2;

	k = 0;
	ver[k].Set(0.0, 0.0, pA[2]); k++;
	ver[k].Set(0.0, 0.0, -pA[2]); k++;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 8; j++)
		{
			v = PI*(-0.25 + i*0.25);
			cos2 = cos(v);
			sin2 = sin(v);
			u = PI*(-1.0 + j*0.25);
			cos1 = cos(u);
			sin1 = sin(u);
			x = pA[0] * sgnp(cos1)*pow(fabs(cos1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
			y = pA[1] * sgnp(sin1)*pow(fabs(sin1), 1.0 / pP[0])*pow(fabs(cos2), 1.0 / pP[1]);
			z = pA[2] * sgnp(sin2)*pow(fabs(sin2), 1.0 / pP[1]);
			ver[k].Set(x, y, z); k++;
		}
	vrrmax = ver[0].Length(); vrrmin = vrrmax;
	for (i = 1; i < 26; i++)
	{
		//cout << i << "\t" << ver[i][0] << "\t" << ver[i][1] << "\t" << ver[i][2] << endl;
		vrr = ver[i].Length();
		if (vrr > vrrmax) vrrmax = vrr;
		if (vrr < vrrmin) vrrmin = vrr;
	}
	rout = vrrmax;
	rin = vrrmin;
}

void CPolyhedron::ChangeSize(double ratio)
{
	for (int i = 0; i < 3; i++)	pA[i] = pA[i] * ratio;
	rout = rout*ratio;
	rin = rin*ratio;
	volume = volume*ratio*ratio*ratio;
}
void CPolyhedron::ChangePosations(CVector center0)
{
	center = center0;
}
void CPolyhedron::ChangeOritations(CVector e0[3])
{
	int i;
	for (i = 0; i < 3; i++) e[i] = e0[i];
}
void CPolyhedron::ChangePosAndOris(CVector center0, CVector e0[3])
{
	int i;
	center = center0;
	for (i = 0; i < 3; i++) e[i] = e0[i];
}
void CPolyhedron::Copyfrom(CPolyhedron *pp)
{
	int i;
	for (i = 0; i < 3; i++)	pA[i] = pp->pA[i];
	for (i = 0; i < 2; i++)	pP[i] = pp->pP[i];
	center = pp->center;
	for (i = 0; i < 3; i++) e[i] = pp->e[i];
	rout = pp->rout;
	rin = pp->rin;
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
void CPolyhedron::Rotate(CVector V, double a)
{
	CVector u;
	double sina, cosa;
	u = V;
	u.Normalize();
	sina = sin(a);
	cosa = cos(a);
	e[0] = e[0] * cosa + u*(e[0] * u)*(1.0 - cosa) + u.Cross(e[0])*sina;
	e[1] = e[1] * cosa + u*(e[1] * u)*(1.0 - cosa) + u.Cross(e[1])*sina;
	e[2] = e[0].Cross(e[1]);
}
void CPolyhedron::Randomize(CVector a[3])
{
	double b, cx, cy, cz;
	CVector X, V;

	cx = (rand() + 1.0) / (RAND_MAX + 2.0), cy = (rand() + 1.0) / (RAND_MAX + 2.0), cz = (rand() + 1.0) / (RAND_MAX + 2.0);
	center = a[0] * cx + a[1] * cy + a[2] * cz;
	do
	{
		X.Set(2.0*rand() / (RAND_MAX + 1.0) - 1.0, 2.0*rand() / (RAND_MAX + 1.0) - 1.0, 2.0*rand() / (RAND_MAX + 1.0) - 1.0);
		b = X.Length();
	} while (b<0.000001 || b>0.9999999);
	X.Normalize();
	V = e[0].Cross(X);
	b = e[0] * X;
	if (b > 1.0) b = 1.0;
	else if (b < -1.0) b = -1.0;
	Rotate(V, acos(b));
	Rotate(X, 2.0*PI*rand() / (RAND_MAX + 1.0));
}
void CPolyhedron::RandomizeCenter(CVector a[3])
{
	double cx, cy, cz;
	CVector X, V;

	cx = (rand() + 1.0) / (RAND_MAX + 2.0), cy = (rand() + 1.0) / (RAND_MAX + 2.0), cz = (rand() + 1.0) / (RAND_MAX + 2.0);
	center = a[0] * cx + a[1] * cy + a[2] * cz;
}
void CPolyhedron::RandomizeOrient()
{
	double b;
	CVector X, V;

	do
	{
		X.Set(2.0*rand() / (RAND_MAX + 1.0) - 1.0, 2.0*rand() / (RAND_MAX + 1.0) - 1.0, 2.0*rand() / (RAND_MAX + 1.0) - 1.0);
		b = X.Length();
	} while (b<0.000001 || b>0.9999999);
	X.Normalize();
	V = e[0].Cross(X);
	b = e[0] * X;
	if (b > 1.0) b = 1.0;
	else if (b < -1.0) b = -1.0;
	Rotate(V, acos(b));
	Rotate(X, 2.0*PI*rand() / (RAND_MAX + 1.0));
}

void CPolyhedron::TranslateRandom(double a, double lambda[3][3])
{
	CVector TV;
	center_L_old = center_L;
	TV.Set(1.0*rand() / double(RAND_MAX) - 0.5, 1.0*rand() / double(RAND_MAX) - 0.5, 1.0*rand() / double(RAND_MAX) - 0.5);
	center_L[0] += TV[0] * a;
	center_L[1] += TV[1] * a;
	center_L[2] += TV[2] * a;
	for (int i = 0; i < 3; i++)
	{
		while (center_L[i] >= 1.0)
			center_L[i] -= 1.0;
		while (center_L[i] < 0.0)
			center_L[i] += 1.0;
	}
	Update_LtoE(lambda);
}
void CPolyhedron::RotateRandom(double a)
{
	CVector RV;
	for (int i = 0; i < 3; i++) e_old[i] = e[i];
	double b;
	do
	{
		RV.Set(2.0*rand() / (RAND_MAX + 1.0) - 1.0, 2.0*rand() / (RAND_MAX + 1.0) - 1.0, 2.0*rand() / (RAND_MAX + 1.0) - 1.0);
		b = RV.Length();
	} while (b<0.000001 || b>0.9999999);
	RV.Normalize();
	Rotate(RV, a*2.0*PI*rand() / (RAND_MAX + 1.0));
}
void CPolyhedron::RetainTrans(double lambda[3][3])
{
	center_L = center_L_old;
	Update_LtoE(lambda);
}
void CPolyhedron::RetainRotate()
{
	for (int i = 0; i < 3; i++) e[i] = e_old[i];
}