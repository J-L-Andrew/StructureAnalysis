#pragma once

#pragma once
#include "Vector.h"
#include "Node.h"
#include <cmath>
#define PI 3.14159265359

class CParticle
{
public:
	CParticle(void);
	CParticle(double);
	CParticle(double, double, double, double);
	~CParticle(void);

	void CParticle(void);
public:
	double p; // x^2p + y^2p + z^2p = 1
	double r_scale[3];
	bool isSphere;
public:
	int ID;
	CVector center; //center coordinates
	CVector e[3];   // axis directions -> rotation matrix

	double dens; //density
	double inertia[3]; //constant in principal coordinates
	double quaternion[4];
	double mass;
	double vol; // volume
	double bound_D; //bounding sphere diameter

	double d_eq;

	//double stiff;


public:
	int ix[2], iy[2], iz[2];
	int dix, diy, diz;
	int nob;
	int* bl, * hl;
	CNode** nl;

public:
	bool isRat;
	int nol;
	CVector tv0, tv1;

	CVector force;//translation vector / force
	CVector torque;//rotation vector / torque
	double f_mag; //force magnitude

	CVector veloc;
	CVector ang_veloc;

	//CVector mid_veloc;
	//CVector mid_ang_veloc;
	CVector mid_mome;

public:

	void Copyfrom(CSuperball* psp);
	void Update();

	void Boundary(double X[2], double Y[2], double Z[2]);// calculate spheropolyhedron boundaries
	void Jump(double PBC[3]);//process periodic boundary condition
	//Translate() & Rotate() change the particle center and axises e
	//To calculate other informations, employ Update()
	void Translate();
	void Translate(double a);
	void Translate(CVector);
	void Rotate(CVector V, double a);//以V方向旋转角度a
	void Rotate(double a);//以RV方向旋转角度RV.length()*a
	void Randomize(double a);
	void Randomize2(double a, CVector X, double betra);
	void Set(double a, int i);
	void RandomizeOrientation();
};

