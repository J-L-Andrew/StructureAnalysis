#ifndef _VOXEL_H
#define _VOXEL_H

#include "/mnt/Edisk/andrew/JammedPacking/3D/code/Matrix.h"
#include "/mnt/Edisk/andrew/JammedPacking/3D/code/Superball.h"
#include "/mnt/Edisk/andrew/JammedPacking/3D/code/Vector.h"

class Voxel {
 public:
  Voxel();
  Voxel(double box, double voxel_l, int nx, int i, int j);
  ~Voxel();

  CVector center;  // center coordinates

  double box_l;
  double dl;
  double area;

  int n;  // number of voxel in each side
  int num_total;
  bool state;  // (1: filled, 0: empty)

 public:
  bool isContain(CSuperball *p);
};

inline double SuperballFunction(double p, double r[3], CVector X) {
  return exp((double)2.0 * p * log(fabs(X[0]) / r[0])) +
         exp((double)2.0 * p * log(fabs(X[1]) / r[1])) +
         exp((double)2.0 * p * log(fabs(X[2]) / r[2]));
}

int initial_voxel(Voxel** s, double box_l, double length);

#endif
