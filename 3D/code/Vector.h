#ifndef _VECTOR_H
#define _VECTOR_H

#include "math.h"

class CVector {
 public:
  double V[3];

 public:
  CVector();
  CVector(const CVector &copy);
  CVector(const double v[3]);
  CVector(const double from[3], const double to[3]);
  CVector(const float v[3]);
  CVector(const double x, double y, double z);
  CVector &operator=(const CVector &copy);
  void SetZero();
  void CopyTo(double v[3]);
  void CopyTo(float v[3]);
  double &operator[](const int i) { return V[i]; }
  operator double *() { return V; }
  CVector &operator=(const double v[3]);
  void Inverse();
  CVector &operator=(const float v[3]);
  CVector MatrixMult(const double M[3][3]);
  CVector &operator*=(const double M[3][3]);
  CVector &operator+=(const double A[3]);
  CVector &operator-=(const double A[3]);
  double operator*(const double A[3]);
  CVector operator+(const double A[3]);
  CVector operator-(const double A[3]);
  CVector operator-() const;
  CVector &operator*=(const double A);
  CVector &operator/=(const double A) { return operator*=(1.0 / A); }
  CVector operator*(const double A);
  CVector operator/(const double A) { return operator*(1.0 / A); }
  CVector Cross(const double A[3]);
  double Dot(CVector that) { return V[0] * that.V[0] + V[1] * that.V[1] + V[2] * that.V[2];}
  void CrossSelf(const double A[3]);
  void AssignDir(const double from[3], const double to[3]);
  double Length();
  double LengthSquare() { return V[0] * V[0] + V[1] * V[1] + V[2] * V[2]; }
  int Normalize();
  int CalPlaneNor(double p1[3], double p2[3], double p3[3]);
  void GenerateOCS(CVector &ocsx, CVector &ocsy, CVector &ocsz,
                   double *sharpAngWithDir = NULL);
  void GenerateOCS(double ocsx[3], double ocsy[3], double ocsz[3],
                   double *sharpAngWithDir = NULL);
  // 判断矢量是否与整体坐标x,y,z重叠,返回0:不重叠
  // 返回1,2,3:分别表示与x,y,z轴重叠
  int AlignWithGlobal();
  void Set(const double x, double y, double z);
  void Set(const double v[3]);

 public:
  static int Normalize(double dv[3]);
  static int Normalize(float dv[3]);
  // 当inIsOnLeftOfM=0时,计算{out} = [M]*{in}  否则计算 {out}={in}*[M]
  static void MatrixMult(const double M[3][3], const double in[3],
                         double out[3], const int inIsOnLeftOfM = 0);
  static void Copy(const double fromSrc[3], double toDest[3]);
  static void Copy(const double fromSrc[3], float toDest[3]);
  //  plusv1v2 = v1 + v2
  static void Plus(const double v1[3], const double v2[3], double plusv1v2[3]);
  //  result = subed - sub
  static void Subtract(const double subed[3], const double sub[3],
                       double result[3]);
  //  result = A*B
  static void Cross(const double a[3], const double b[3], double result[3]);
  // Value = A . B
  static double Dot(const double a[3], const double b[3]);
  // 判断两个方向是否相等
  static int DirEqual(const double dir1[3], const double dir2[3]);
  static double Length(const double a[3]);
  static double LengthSquare(const double a[3]) {
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  }
  static float Length(const float a[3]);
};

inline CVector::CVector() { SetZero(); }
inline CVector::CVector(const CVector &copy) {
  V[0] = copy.V[0];
  V[1] = copy.V[1];
  V[2] = copy.V[2];
}
inline CVector::CVector(const double v[3]) {
  V[0] = v[0];
  V[1] = v[1];
  V[2] = v[2];
}
inline CVector::CVector(const float v[3]) {
  V[0] = v[0];
  V[1] = v[1];
  V[2] = v[2];
}
inline CVector::CVector(const double from[3], const double to[3]) {
  V[0] = to[0] - from[0];
  V[1] = to[1] - from[1];
  V[2] = to[2] - from[2];
}
inline CVector::CVector(const double x, double y, double z) {
  V[0] = x;
  V[1] = y;
  V[2] = z;
}
inline CVector &CVector::operator=(const CVector &copy) {
  V[0] = copy.V[0];
  V[1] = copy.V[1];
  V[2] = copy.V[2];
  return *this;
}
inline void CVector::SetZero() { V[0] = V[1] = V[2] = 0.0; }
inline void CVector::CopyTo(double v[3]) {
  v[0] = V[0];
  v[1] = V[1];
  v[2] = V[2];
}
inline void CVector::CopyTo(float v[3]) {
  v[0] = (float)V[0];
  v[1] = (float)V[1];
  v[2] = (float)V[2];
}
inline void CVector::Inverse() {
  V[0] = -V[0];
  V[1] = -V[1];
  V[2] = -V[2];
}
inline CVector &CVector::operator=(const float v[3]) {
  V[0] = v[0];
  V[1] = v[1];
  V[2] = v[2];
  return *this;
}
inline CVector &CVector::operator=(const double v[3]) {
  V[0] = v[0];
  V[1] = v[1];
  V[2] = v[2];
  return *this;
}
// Normalize the vector
inline int CVector::Normalize() { return Normalize(V); }
inline void CVector::Copy(const double fromSrc[3], double toDest[3]) {
  toDest[0] = fromSrc[0];
  toDest[1] = fromSrc[1];
  toDest[2] = fromSrc[2];
}
inline void CVector::Copy(const double fromSrc[3], float toDest[3]) {
  toDest[0] = (float)fromSrc[0];
  toDest[1] = (float)fromSrc[1];
  toDest[2] = (float)fromSrc[2];
}
inline void CVector::Plus(const double v1[3], const double v2[3],
                          double plusv1v2[3]) {
  plusv1v2[0] = v1[0] + v2[0];
  plusv1v2[1] = v1[1] + v2[1];
  plusv1v2[2] = v1[2] + v2[2];
}
inline void CVector::Subtract(const double subed[3], const double sub[3],
                              double result[3]) {
  result[0] = subed[0] - sub[0];
  result[1] = subed[1] - sub[1];
  result[2] = subed[2] - sub[2];
}
inline double CVector::Dot(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline void CVector::Set(const double x, double y, double z) {
  V[0] = x;
  V[1] = y;
  V[2] = z;
}
inline void CVector::Set(const double v[3]) {
  V[0] = v[0];
  V[1] = v[1];
  V[2] = v[2];
}
#endif