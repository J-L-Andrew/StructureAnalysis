#pragma once
#if !defined(_VECTOR_H)
#define _VECTOR_H
#include "math.h"
class CVector
{
public:
	double V[2];
public:
	CVector();
	CVector(const double v[2]);
	CVector(const CVector & copy);
	CVector(const double x, double y);
	void Set(const double v[2]);
	void Set(const  CVector & copy);
	void Set(const double x, double y);
	operator double * () { return V; }
	double & operator[](const int i) { return V[i]; }
	CVector operator -() const;//每个变量取负
	CVector & operator=(const double v[2]);	//等于
	CVector & operator +=(const double A[2]);//加等于
	CVector & operator -=(const double A[2]);//减等于
	CVector & operator *=(const double A);//数乘等于
	CVector & operator /=(const double A);//数除等于
	CVector operator +(const double A[2]);//加
	CVector operator -(const double A[2]);//减	
	CVector operator *(const double A);//数乘
	CVector operator /(const double A);//数除
	double operator *(const double A[2]);//点乘
	CVector Cross(const double A[2]);//叉乘	
	void SetZero();
	double Length();
	double LengthSquare();
	int Normalize();
public:
	int Normalize(double dv[2]);
};

inline CVector::CVector()
{
	V[0] = V[1] = 0.0;
}
inline CVector::CVector(const double v[2])
{
	V[0] = v[0]; V[1] = v[1];
}
inline CVector::CVector(const CVector & copy)
{
	V[0] = copy.V[0];	V[1] = copy.V[1];
}
inline CVector::CVector(const double x, double y)
{
	V[0] = x; V[1] = y;
}
inline void CVector::Set(const double v[2])
{
	V[0] = v[0]; V[1] = v[1];
}
inline void CVector::Set(const  CVector & copy)
{
	V[0] = copy.V[0]; V[1] = copy.V[1];
}
inline void CVector::Set(const double x, double y)
{
	V[0] = x; V[1] = y;
}
inline CVector CVector::operator -() const
{
	double A[2];
	A[0] = -V[0]; A[1] = -V[1];
	return A;
}
inline CVector & CVector::operator=(const double v[2])
{
	V[0] = v[0]; V[1] = v[1];
	return *this;
}
inline CVector & CVector::operator +=(const double A[2])
{
	V[0] += A[0]; V[1] += A[1];
	return *this;
}
inline CVector & CVector::operator -=(const double A[2])
{
	V[0] -= A[0]; V[1] -= A[1];
	return *this;
}
inline CVector & CVector::operator *=(const double A)
{
	V[0] *= A; V[1] *= A;
	return *this;
}
inline CVector & CVector::operator /=(const double A)
{
	V[0] /= A; V[1] /= A;
	return *this;
}
inline CVector CVector::operator +(const double A[2])
{
	double B[3];
	B[0] = V[0] + A[0];	B[1] = V[1] + A[1];
	return B;
}
inline CVector CVector::operator -(const double A[2])
{
	double B[2];
	B[0] = V[0] - A[0];	B[1] = V[1] - A[1];
	return B;
}
inline CVector CVector::operator *(const double A)
{
	double B[2];
	B[0] = V[0] * A; B[1] = V[1] * A;
	return B;
}
inline CVector CVector::operator /(const double A)
{
	double B[2];
	B[0] = V[0] / A; B[1] = V[1] / A;
	return B;
}
inline double CVector::operator *(const double A[2])
{
	return (V[0] * A[0] + V[1] * A[1]);
}
inline CVector CVector::Cross(const double A[2])
{
	double d[1];
	d[0] = V[0] * A[1] - V[1] * A[0];
	return d;
}
inline void CVector::SetZero()
{
	V[0] = V[1] = 0.0;
}
inline double CVector::Length()
{
	return sqrt(V[0] * V[0] + V[1] * V[1]);
}
inline double CVector::LengthSquare()
{
	return (V[0] * V[0] + V[1] * V[1]);
}
inline int CVector::Normalize()//归一化函数
{
	return Normalize(V);
}
#endif


