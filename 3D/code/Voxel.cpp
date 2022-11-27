#include "Voxel.h"

Voxel::Voxel() { state = 0; }

// state=0 by default
Voxel::Voxel(double box, double voxel_l, int nx, int i, int j) {
  box_l = box;

  dl = voxel_l;
  area = voxel_l * voxel_l;

  n = nx;
  num_total = n * n;

  x = double(i) * dl;
  y = double(j) * dl;

  state = 1;
}

Voxel::~Voxel() {}

bool Voxel::isContain(CSuperball *pa) {
  CVector r = pa->center;
  double p = pa->p;
  double rs[3];
  rs[0] = pa->r_scale[0];
  rs[1] = pa->r_scale[1];
  rs[2] = pa->r_scale[2];

  double A[3][3], AT[3][3];  // transformation matrix
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      // Aa, Ab as rotational matrix of pa, pb
      A[j][i] = pa->e[i][j];
      AT[i][j] = pa->e[i][j];
    }

  CVector rt = center.operator-(r);
  CVector r_l;
  MatrixMult(AT, rt, r_l);
  double xa = SuperballFunction(p, rs, r_l);
  double ta = exp((double)(1.0 / p) * log(xa));

  if (ta < 0.0)
    return 1;
  else
    return 0;
}

int initial_voxel(Voxel** s, double box_l, double length) {
  int nx = int(box_l / length);
  int ntotal = nx * nx;
  double dx = box_l / double(nx);

  for (int j = 0; j < nx; ++j) {
    for (int i = 0; i < nx; ++i) {
      int id = j * nx + i;
      s[id] = new Voxel(box_l, dx, nx, i, j);
    }
  }
  cout << "Voxel initialization=" << ntotal << endl;
  return ntotal;
}
