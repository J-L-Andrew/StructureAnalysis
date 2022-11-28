#ifndef _SUPERBALL_H
#define _SUPERBALL_H

#include <cmath>

#include "Node.h"
#include "Vector.h"
#define PI 3.14159265359

class CSuperball {
 public:
  CSuperball();
  CSuperball(double); // p
  CSuperball(double, double, double, double); // (p, a, b, c)
  ~CSuperball();

  void Destructor(void);

 public:
  double p;  // x^2p + y^2p + z^2p = 1
  double r_scale[3];

 public:
  int ID;
  double d_v;
  CVector center;  // center coordinates
  CVector e[3];    // axis directions
  double vol;      // volume
  double bound_D;  // bounding sphere diameter

 public:
  int ix[2], iy[2], iz[2];
  int dix, diy, diz;
  int nob;
  int *bl, *hl;
  CNode** nl;

 public:
  bool isRat;
  int nol;
  CVector tv0, tv1;
  CVector TV;    // translation vector / force
  CVector RV;    // rotation vector / torque
  double f_mag;  // force magnitude

 public:
  void Copyfrom(CSuperball* psp);
  void Update();

  void Boundary(double X[2], double Y[2],
                double Z[2]);  // calculate spheropolyhedron boundaries
  void Jump(double PBC[3]);    // process periodic boundary condition
  // Translate() & Rotate() change the particle center and axises e
  // To calculate other informations, employ Update()
  void Translate();
  void Translate(double a);
  void Translate(CVector);
  void Rotate(CVector V, double a);  //以V方向旋转角度a
  void Rotate(double a);             //以RV方向旋转角度RV.length()*a
  void Randomize(double a);
  void Randomize2(double a, CVector X, double betra);
  void Set(double a, int i);
  void RandomizeOrientation();
};

#endif
