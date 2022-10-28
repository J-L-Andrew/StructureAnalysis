#pragma once
#include "Vector.h"


void MatrixCopy(double M1[2][2], double M2[2][2])//前面copy到后面
{
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			M2[i][j] = M1[i][j];
}


void MatrixAdd(double M1[2][2], double M2[2][2])
{
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			M1[i][j] = M1[i][j] + M2[i][j];
}


double MatrixDet(double m[2][2])
{
	return -m[0][1] * m[1][0] + m[0][0] * m[1][1];
}


double MatrixMax(double m[2][2])//输出绝对值最大的矩阵元
{
	double maxv;
	maxv = 0.0;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
		{
			if (maxv < fabs(m[i][j]))
				maxv = fabs(m[i][j]);
		}
	if (maxv == 0.0)
		maxv = 1.0;
	return maxv;
}

void MatrixMult(double M1[2][2], double a, double M2[2][2])//矩阵数乘得到另一个矩阵
{
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			M2[i][j] = M1[i][j] * a;
}
void MatrixMult(double M[2][2], CVector X, CVector &res)//矩阵乘向量得到另一个矩阵
{
	for (int i = 0; i < 2; i++)
	{
		res[i] = 0;
		for (int j = 0; j < 2; j++)
		{
			res[i] += (M[i][j] * X[j]);
		}
	}
}
void MatrixMult(double M1[2][2], double M2[2][2], double M3[2][2])//矩阵乘矩阵得到另一个矩阵
{
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
		{
			M3[i][j] = 0;
			for (int k = 0; k < 2; ++k)
			{
				M3[i][j] += (M1[i][k] * M2[k][j]);
			}
		}
}

void Dyadic(CVector X1, CVector X2, double M[2][2])//向量相乘得到矩阵
{
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			M[i][j] = X1[i] * X2[j];
		}
	}
}

void MatrixInv(double M[2][2], double M1[2][2])
{
	double det = 1.0 / MatrixDet(M); //must check first
	M1[0][0] = M[1][1];
	M1[0][1] = -M[0][1];
	M1[1][0] = -M[1][0];
	M1[1][1] = M[0][0];

	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
		{
			M1[i][j] *= det;
		}
}
void MatrixInv(double M[2][2], double M1[2][2], double DetM)//一开始是伴随矩阵求逆，现在二维简单化
{
	double det = 1.0 / DetM; //must check first
	M1[0][0] = M[1][1];
	M1[0][1] = -M[0][1];
	M1[1][0] = -M[1][0];
	M1[1][1] = M[0][0];

	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
		{
			M1[i][j] *= det;
		}
}
