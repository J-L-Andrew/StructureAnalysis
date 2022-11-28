#ifndef _MATRIX_H
#define _MATRIX_H

#include "Vector.h"

inline void MatrixCopy(double M1[3][3], double M2[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M2[i][j] = M1[i][j];
}

inline void MatrixTrans(double M1[3][3], double M2[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      M2[i][j] = M1[j][i];
    }
  }
}

inline void MatrixAdd(double M1[3][3], double M2[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) M1[i][j] = M1[i][j] + M2[i][j];
}

inline void MatrixMult(double M1[3][3], double a, double M2[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M2[i][j] = M1[i][j] * a;
}

inline void MatrixMult(double M[3][3], CVector X, CVector &res) {
  for (int i = 0; i < 3; ++i) {
    res[i] = 0;
    for (int j = 0; j < 3; ++j) {
      res[i] += (M[i][j] * X[j]);
    }
  }
}
inline void MatrixMult(double M1[3][3], double M2[3][3], double M3[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      M3[i][j] = 0;
      for (int k = 0; k < 3; ++k) {
        M3[i][j] += (M1[i][k] * M2[k][j]);
      }
    }
}

inline void Dyadic(CVector X1, CVector X2, double M[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      M[i][j] = X1[i] * X2[j];
    }
  }
}

inline double MatrixDet(double m[3][3]) {
  return -m[0][2] * m[1][1] * m[2][0] + m[0][1] * m[1][2] * m[2][0] +
         m[0][2] * m[1][0] * m[2][1] - m[0][0] * m[1][2] * m[2][1] -
         m[0][1] * m[1][0] * m[2][2] + m[0][0] * m[1][1] * m[2][2];
}

inline void MatrixInv(double M[3][3], double M1[3][3]) {
  double det = 1.0 / MatrixDet(M);                    // must check first
  M1[0][0] = -M[1][2] * M[2][1] + M[1][1] * M[2][2];  //-m12 m21 + m11 m22
  M1[0][1] = M[0][2] * M[2][1] - M[0][1] * M[2][2];   // m02 m21 - m01 m22
  M1[0][2] = -M[0][2] * M[1][1] + M[0][1] * M[1][2];  //-m02 m11 + m01 m12

  M1[1][0] = M[1][2] * M[2][0] - M[1][0] * M[2][2];   // m12 m20 - m10 m22
  M1[1][1] = -M[0][2] * M[2][0] + M[0][0] * M[2][2];  //-m02 m20 + m00 m22
  M1[1][2] = M[0][2] * M[1][0] - M[0][0] * M[1][2];   // m02 m10 - m00 m12

  M1[2][0] = -M[1][1] * M[2][0] + M[1][0] * M[2][1];  //-m11 m20 + m10 m21
  M1[2][1] = M[0][1] * M[2][0] - M[0][0] * M[2][1];   // m01 m20 - m00 m21
  M1[2][2] = -M[0][1] * M[1][0] + M[0][0] * M[1][1];  //-m01 m10 + m00 m11

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      M1[i][j] *= det;
    }
}

#endif
