#ifndef _VOXEL_H
#define _VOXEL_H

#include "/mnt/Edisk/andrew/JammedPacking/3D/code/Matrix.h"
#include "/mnt/Edisk/andrew/JammedPacking/3D/code/Superball.h"
#include "/mnt/Edisk/andrew/JammedPacking/3D/code/Vector.h"

class Voxel {
 public:
  Voxel();
  Voxel(double box, double voxel_l, int nx);
  ~Voxel();

  CVector center;  // center coordinates
  CVector r; // for indicator function

  double box_l;
  double dl;
  double area;

  int n;  // number of voxel in each side
  int num_total;
  bool state;  // (1: filled, 0: empty)

 public:
  bool isContain(CSuperball *p);
};

inline int sgn(double x) {
  int t;
  if (x == 0)
    t = 0;
  else if (x > 0)
    t = 1;
  else
    t = -1;
  return t;
}

inline double SuperballFunction(double p, double r[3], CVector X) {
  return exp((double)2.0 * p * log(fabs(X[0]) / r[0])) +
         exp((double)2.0 * p * log(fabs(X[1]) / r[1])) +
         exp((double)2.0 * p * log(fabs(X[2]) / r[2]));
}

void SuperballFunctionD(double p, double r[3], CVector X, CVector &xn) {
  xn[0] = 2.0 * p * exp((2.0 * p - 1.0) * log(fabs(X[0]) / r[0])) *
          (double)sgn(X[0]) / r[0];
  xn[1] = 2.0 * p * exp((2.0 * p - 1.0) * log(fabs(X[1]) / r[1])) *
          (double)sgn(X[1]) / r[1];
  xn[2] = 2.0 * p * exp((2.0 * p - 1.0) * log(fabs(X[2]) / r[2])) *
          (double)sgn(X[2]) / r[2];
}

// Discretizing the space around a given particle into pixels
// Calculate the Fourier transform of the indicator function for particle j
void FFT_if(CSuperball* p, double resolution, CVector k, double &Re, double &Im);


int initial_voxel(Voxel** s, double box_l, double length);



#endif
