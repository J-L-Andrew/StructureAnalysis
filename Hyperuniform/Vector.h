#pragma once
#if !defined(_VECTOR_H)
#define _VECTOR_H
#include "math.h"
class CVector
{
public:
	double V[3];
public:
	CVector();
	CVector(const double v[3]);
	CVector(const CVector& copy);
	CVector(const double x, double y, double z);
	void Set(const double v[3]);
	void Set(const  CVector& copy);
	void Set(const double x, double y, double z);
	operator double* () { return V; }
	double& operator[](const int i) { return V[i]; }
	CVector operator -() const;//每个变量取负
	CVector& operator=(const double v[3]);	//等于
	CVector& operator +=(const double A[3]);//加等于
	CVector& operator -=(const double A[3]);//减等于
	CVector& operator *=(const double A);//数乘等于
	CVector& operator /=(const double A);//数除等于
	CVector operator +(const double A[3]);//加
	CVector operator -(const double A[3]);//减	
	CVector operator *(const double A);//数乘
	CVector operator /(const double A);//数除
	double operator *(const double A[3]);//点乘
	CVector Cross(const double A[3]);//叉乘	
	void SetZero();
	double Length();
	double LengthSquare();
	int Normalize();
public:
	int Normalize(double dv[3]);
};

inline CVector::CVector()
{
	V[0] = V[1] = V[2] = 0.0;
}
inline CVector::CVector(const double v[3])
{
	V[0] = v[0]; V[1] = v[1]; V[2] = v[2];
}
inline CVector::CVector(const CVector& copy)
{
	V[0] = copy.V[0];	V[1] = copy.V[1];	V[2] = copy.V[2];
}
inline CVector::CVector(const double x, double y, double z)
{
	V[0] = x; V[1] = y; V[2] = z;
}
inline void CVector::Set(const double v[3])
{
	V[0] = v[0]; V[1] = v[1]; V[2] = v[2];
}
inline void CVector::Set(const  CVector& copy)
{
	V[0] = copy.V[0]; V[1] = copy.V[1]; V[2] = copy.V[2];
}
inline void CVector::Set(const double x, double y, double z)
{
	V[0] = x; V[1] = y; V[2] = z;
}
inline CVector CVector::operator -() const
{
	double A[3];
	A[0] = -V[0]; A[1] = -V[1]; A[2] = -V[2];
	return A;
}
inline CVector& CVector::operator=(const double v[3])
{
	V[0] = v[0]; V[1] = v[1]; V[2] = v[2];
	return *this;
}
inline CVector& CVector::operator +=(const double A[3])
{
	V[0] += A[0]; V[1] += A[1]; V[2] += A[2];
	return *this;
}
inline CVector& CVector::operator -=(const double A[3])
{
	V[0] -= A[0]; V[1] -= A[1]; V[2] -= A[2];
	return *this;
}
inline CVector& CVector::operator *=(const double A)
{
	V[0] *= A; V[1] *= A; V[2] *= A;
	return *this;
}
inline CVector& CVector::operator /=(const double A)
{
	V[0] /= A; V[1] /= A; V[2] /= A;
	return *this;
}
inline CVector CVector::operator +(const double A[3])
{
	double B[3];
	B[0] = V[0] + A[0];	B[1] = V[1] + A[1];	B[2] = V[2] + A[2];
	return B;
}
inline CVector CVector::operator -(const double A[3])
{
	double B[3];
	B[0] = V[0] - A[0];	B[1] = V[1] - A[1];	B[2] = V[2] - A[2];
	return B;
}
inline CVector CVector::operator *(const double A)
{
	double B[3];
	B[0] = V[0] * A; B[1] = V[1] * A; B[2] = V[2] * A;
	return B;
}
inline CVector CVector::operator /(const double A)
{
	double B[3];
	B[0] = V[0] / A; B[1] = V[1] / A; B[2] = V[2] / A;
	return B;
}
inline double CVector::operator *(const double A[3])
{
	return (V[0] * A[0] + V[1] * A[1] + V[2] * A[2]);
}
inline CVector CVector::Cross(const double A[3])
{
	double d[3];
	d[0] = V[1] * A[2] - V[2] * A[1];
	d[1] = V[2] * A[0] - V[0] * A[2];
	d[2] = V[0] * A[1] - V[1] * A[0];
	return d;
}
inline void CVector::SetZero()
{
	V[0] = V[1] = V[2] = 0.0;
}
inline double CVector::Length()
{
	return sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
}
inline double CVector::LengthSquare()
{
	return (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
}
inline int CVector::Normalize()
{
	return Normalize(V);
}
#endif


