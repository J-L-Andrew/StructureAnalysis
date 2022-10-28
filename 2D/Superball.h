#pragma once
#include "math.h"
#include "Vector.h"
#include <cstdlib>
#define PI 3.14159265358979324

class CPolyhedron//超椭球类
{
public:
	CPolyhedron(void);
	CPolyhedron(double ppara[3], double rrr[2]);
	~CPolyhedron(void);

public:
	double pA[3];//x,y,z轴长度，a=b;  p:指数 |x/a|^2p+|y/b|^2p+|z/c|^2p=1;
	double area;
	double rout, rin;//外接球，内切球半径
	CVector center;	//质心（欧拉坐标系）
	CVector center_L;//质心	拉格朗日坐标系）
	CVector e[2];//三轴方向，e[2]为主轴Z方向，另外两个为垂直于主轴的平面内的互相垂直的任意两个向量。

	CVector center_L_old;
	CVector e_old[2];
	int ix[2], icell;
	int kind;

public:
	double sgnp(double x);
	void Getroutrin();
	void ChangeSize(double ratio);
	void ChangePosations(CVector center0);
	void ChangeOritations(CVector e0[2]);
	void ChangePosAndOris(CVector center0, CVector e0[2]);

	void Update_EtoL(double L_inv[3][3]);
	void Update_LtoE(double lambda[3][3]);
	void Copyfrom(CPolyhedron * pp);
	void Jump(CVector TV);//process periodic boundary condition	


	void TranslateRandom(double a, double lambda[3][3]);
	void RotateRandom(double a);
	void RetainTrans(double lambda[3][3]);
	void RetainRotate();
};

