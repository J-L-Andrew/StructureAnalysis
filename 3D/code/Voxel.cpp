#include "Voxel.h"

Voxel::Voxel() { state = 0; }

// state=0 by default
Voxel::Voxel(double box, double voxel_l, int nx) {
  box_l = box;

  dl = voxel_l;
  area = voxel_l * voxel_l;

  n = nx;
  num_total = n * n;

  state = 0;
}

Voxel::~Voxel() {}

bool Voxel::isContain(CSuperball* pa) {
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

void FFT_if(CSuperball* p, double resolution, double L, CVector K_t, double& Re,
            double& Im) {
  CVector center_min, center_max;
  CVector bound;
  bound.Set(0.5 * p->bound_D, 0.5 * p->bound_D, 0.5 * p->bound_D);
  center_min = p->center.operator-(bound);
  center_max = p->center.operator+(bound);

  // Voxel initialization
  int Nvoxel = int(pow(p->bound_D / resolution, 3.0) * 1.1);
  Voxel** vox;
  vox = new Voxel*[Nvoxel];

  int nx = int(p->bound_D / resolution);
  if (nx % 2 == 0) nx += 1;  // make sure that nx is an even integer
  int ntotal = nx * nx * nx;
  double dx = p->bound_D / double(nx);

  int mid = (nx - 1) / 2;

  Re = Im = 0.0;
  // Discretizing the space around a given particle into pixels
  for (int k = 0; k < nx; ++k) {
    for (int j = 0; j < nx; ++j) {
      for (int i = 0; i < nx; ++i) {
        int id = k * nx * nx + j * nx + i;
        vox[id] = new Voxel(p->bound_D, dx, nx);

        vox[id]->r.Set(double(i - mid) * dx, double(j - mid) * dx,
                       double(k - mid) * dx);
        vox[id]->center = vox[id]->r + p->center;

        if (vox[id]->isContain(p)) {
          vox[id]->state = 1;
          double temp = -K_t.Dot(vox[id]->r);
          Re += cos(temp);
          Im += sin(temp);
        }
      }
    }
  }
}
