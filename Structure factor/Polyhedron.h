#pragma once
#include "math.h"
#include "Vector.h"
#include <cstdlib>
#define PI 3.14159265358979324

class CPolyhedron//超椭球类
{
public:
	CPolyhedron(void);
	CPolyhedron(double ppara[5], double rrr[2]);
	~CPolyhedron(void);

public:
	double pA[3], pP[2];//x,y,z轴长度，a=b;  p:指数 |x/a|^2p+|y/b|^2p+|z/c|^2p=1;
	double volume;
	double rout, rin;//外接球，内切球半径
	CVector center;	//质心（欧拉坐标系）
	CVector center_L;//质心	拉格朗日坐标系）
	CVector e[3];//三轴方向，e[2]为主轴Z方向，另外两个为垂直于主轴的平面内的互相垂直的任意两个向量。

	CVector center_L_old;
	CVector e_old[3];
	int ix[3], icell;
	int kind;

public:
	double sgnp(double x);
	void Getroutrin();
	void ChangeSize(double ratio);
	void ChangePosations(CVector center0);
	void ChangeOritations(CVector e0[3]);
	void ChangePosAndOris(CVector center0, CVector e0[3]);

	void Update_EtoL(double L_inv[3][3]);
	void Update_LtoE(double lambda[3][3]);
	void Copyfrom(CPolyhedron * pp);
	void Jump(CVector TV);//process periodic boundary condition	

	void Rotate(CVector V, double a);//以V方向旋转角度a
	void Randomize(CVector a[3]);
	void RandomizeCenter(CVector a[3]);
	void RandomizeOrient();

	void TranslateRandom(double a, double lambda[3][3]);
	void RotateRandom(double a);
	void RetainTrans(double lambda[3][3]);
	void RetainRotate();
};

