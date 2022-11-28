#include "Vector.h"

#include <memory.h>
#include <stdio.h>

double Sqrt(double d) {
  if (d < 1.0E-20) return 0.0;
  return sqrt(d);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
void CVector::MatrixMult(const double M[3][3], const double in[3],
                         double out[3], const int inIsOnLeftOfM) {
  int i;
  double B[3];
  if (inIsOnLeftOfM) {
    for (i = 0; i < 3; i++) {
      B[i] = M[0][i] * in[0] + M[1][i] * in[1] + M[2][i] * in[2];
    }
  } else {
    for (i = 0; i < 3; i++) {
      B[i] = M[i][0] * in[0] + M[i][1] * in[1] + M[i][2] * in[2];
    }
  }
  out[0] = B[0];
  out[1] = B[1];
  out[2] = B[2];
}
// 计算M*vector
inline CVector CVector::MatrixMult(const double M[3][3]) {
  double B[3];
  for (int i = 0; i < 3; i++) {
    B[i] = M[i][0] * V[0] + M[i][1] * V[1] + M[i][2] * V[2];
  }
  return B;
}
// 计算vector*M
CVector &CVector::operator*=(const double M[3][3]) {
  double B[3];
  for (int i = 0; i < 3; i++) {
    B[i] = V[0] * M[0][i] + V[1] * M[1][i] + V[2] * M[2][i];
  }
  memcpy(V, B, sizeof(double) * 3);
  return *this;
}

CVector &CVector::operator+=(const double A[3]) {
  for (int i = 0; i < 3; i++) V[i] += A[i];
  return *this;
}
CVector &CVector::operator-=(const double A[3]) {
  for (int i = 0; i < 3; i++) V[i] -= A[i];
  return *this;
}
// 两个向量作点积
double CVector::operator*(const double A[3]) {
  return (V[0] * A[0] + V[1] * A[1] + V[2] * A[2]);
}
CVector CVector::operator+(const double A[3]) {
  double B[3];
  for (int i = 0; i < 3; i++) B[i] = V[i] + A[i];
  return B;
}
CVector CVector::operator-(const double A[3]) {
  double B[3];
  for (int i = 0; i < 3; i++) B[i] = V[i] - A[i];
  return B;
}
CVector CVector::operator-() const {
  double A[3];
  for (int i = 0; i < 3; i++) A[i] = -V[i];
  return A;
}
CVector &CVector::operator*=(const double A) {
  for (int i = 0; i < 3; i++) V[i] *= A;
  return *this;
}
CVector CVector::operator*(const double A) {
  double B[3];
  for (int i = 0; i < 3; i++) B[i] = V[i] * A;
  return B;
}

//  返回vector与A的叉乘
CVector CVector::Cross(const double A[3]) {
  double d[3];
  d[0] = V[1] * A[2] - V[2] * A[1];
  d[1] = V[2] * A[0] - V[0] * A[2];
  d[2] = V[0] * A[1] - V[1] * A[0];
  return d;
}
// vector与A的叉乘放在vector内
void CVector::CrossSelf(const double A[3]) {
  double d[3];
  d[0] = V[1] * A[2] - V[2] * A[1];
  d[1] = V[2] * A[0] - V[0] * A[2];
  d[2] = V[0] * A[1] - V[1] * A[0];
  V[0] = d[0];
  V[1] = d[1];
  V[2] = d[2];
}
// Vector=B-A;
void CVector::AssignDir(const double from[3], const double to[3]) {
  for (int i = 0; i < 3; i++) V[i] = to[i] - from[i];
}

// Get the length of a vector
double CVector::Length() {
  double len = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
  // if( len > 1.0E-28 ) return (Sqrt(len) );
  // return 0.0;
  return sqrt(len);
}
int CVector::Normalize(double dv[3]) {
  double len = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
  // if( len < 1.0E-20 ) return 0;
  // if( len>0.99999999 && len<1.00000001 ) return 1;
  // len = Sqrt(len);
  len = sqrt(len);
  dv[0] /= len;
  dv[1] /= len;
  dv[2] /= len;
  return 1;
}
int CVector::Normalize(float dv[3]) {
  float len = (float)(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
  if (len < 1.0E-10) return 0;
  if (len > 0.99999999 && len < 1.00000001) return 1;
  len = (float)Sqrt(len);
  dv[0] /= len;
  dv[1] /= len;
  dv[2] /= len;
  return 1;
}
int CVector::CalPlaneNor(double p1[3], double p2[3], double p3[3]) {
  CVector p13;
  AssignDir(p1, p2);
  p13.AssignDir(p1, p3);
  CrossSelf(p13);
  if (!Normalize()) return 0;
  return 1;
}

void CVector::GenerateOCS(double ocsx[3], double ocsy[3], double ocsz[3],
                          double *sharpAngWithDir) {
  CVector vx, vy, vz;
  GenerateOCS(vx, vy, vz, sharpAngWithDir);
  vx.CopyTo(ocsx);
  vy.CopyTo(ocsy);
  vz.CopyTo(ocsz);
}
int CVector::AlignWithGlobal() {
  Normalize();
  if (fabs(V[0]) > 0.999) return 1;
  if (fabs(V[1]) > 0.999) return 2;
  if (fabs(V[2]) > 0.999) return 3;
  return 0;
}
// 由当前的观察方向(向量方向指向眼睛)计算局部坐标轴的方向
// 一般情况下整体坐标系下的Z轴处在局部坐标系下的zy平面内
// sharpAngWithDir!=NULL: ocsz必须与sharpAngWithDir之间的夹角为锐角
void CVector::GenerateOCS(CVector &ocsx, CVector &ocsy, CVector &ocsz,
                          double *sharpAngWithDir) {
  CVector xyz;

  ocsz = *this;
  ocsz.Normalize();
  if (sharpAngWithDir && ocsz * sharpAngWithDir < 0.0) {
    ocsz.Inverse();
  }
  if (fabs(ocsz[0]) < 1.0 / 64.0 && fabs(ocsz[1]) < 1.0 / 64.0) {
    if (ocsz[2] > 0.0)
      xyz[1] = 1.0;
    else
      xyz[1] = -1.0;
  } else {
    if (ocsz[2] >= -1.0E-10)
      xyz[2] = 1.0;
    else
      xyz[2] = -1.0;
  }
  ocsx = xyz.Cross(ocsz);
  ocsx.Normalize();
  ocsy = ocsz.Cross(ocsx);
}
void CVector::Cross(const double a[3], const double b[3], double result[3]) {
  double d[3];
  d[0] = a[1] * b[2] - a[2] * b[1];
  d[1] = a[2] * b[0] - a[0] * b[2];
  d[2] = a[0] * b[1] - a[1] * b[0];
  for (int i = 0; i < 3; i++) {
    result[i] = d[i];
  }
}
int CVector::DirEqual(const double dir1[3], const double dir2[3]) {
  if (fabs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) > 1.0E-3) return 0;
  if (fabs(dir1[1] * dir2[2] - dir1[2] * dir2[1]) > 1.0E-3) return 0;
  if (fabs(dir1[0] * dir2[2] - dir1[2] * dir2[0]) > 1.0E-3) return 0;
  return 1;
}
// Get the length of a vector
double CVector::Length(const double a[3]) {
  double len = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  // if( len > 1.0E-15 ) return (Sqrt(len) );
  // return 0.0;
  return sqrt(len);
}
float CVector::Length(const float a[3]) {
  double len = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  if (len > 1.0E-15) return ((float)Sqrt(len));
  return 0.0;
}
