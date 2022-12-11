#include "analysis.h"

int NUM_PARTICLE;
double L;  // Box length

double sumvol;
double PD;

int Nkind;
int *NpEach;
double **Para;  // (p, a, b, c)

vector<CSuperball *> particles;

ifstream myfile;

int main() {
  string file = "config1.pac";
  readfile(file);

  // POV_superball(1);

  // OutputBoundary("../vorofile/boundary.txt", 1);
  // OutputVoroPoint(20, 1);

  VoroVolume(20, 1);

  // spectral_density();

  releasespace();
}

void readfile(string filename) {
  myfile.open(filename);

  if (!myfile.is_open()) {
    cout << "Unable to open myfile" << endl;
    exit(1);
  }

  myfile >> NUM_PARTICLE;
  myfile >> Nkind;
  myfile >> L;

  askforspace();
  for (int i = 0; i < Nkind; ++i) {
    myfile >> NpEach[i] >> Para[i][0] >> Para[i][1] >> Para[i][2] >> Para[i][3];
  }

  double p, r_scale[3];
  sumvol = 0.0;
  for (int i = 0; i < NUM_PARTICLE; ++i) {
    myfile >> p >> r_scale[0] >> r_scale[1] >> r_scale[2];
    particles.push_back(new CSuperball(p, r_scale[0], r_scale[1], r_scale[2]));
    particles[i]->ID = i;

    myfile >> particles[i]->center[0] >> particles[i]->center[1] >>
        particles[i]->center[2];

    for (int j = 0; j < 3; ++j)
      myfile >> particles[i]->e[j][0] >> particles[i]->e[j][1] >>
          particles[i]->e[j][2];

    sumvol += particles[i]->vol;
  }

  PD = sumvol / L / L / L;
}

// ----------------
// Ask and release space
// ----------------
void askforspace() {
  NpEach = new int[Nkind];
  Para = new double *[Nkind];
  for (int i = 0; i < Nkind; ++i) Para[i] = new double[4];
}

void releasespace() {
  delete[] NpEach;
  for (int i = 0; i < Nkind; ++i) delete Para[i];
  delete[] Para;
}

// ----------------
// Hyperuniform Analysis
// ----------------
void Structure_factor(int species, bool global = true) {
  int nk = int(L / sqrt(3.0));  // a rough estimate
  double thickness = 2.0 * PI / L;

  double Sf[nk];
  int num[nk];

  for (int i = 0; i < nk; ++i) {
    Sf[i] = 0.0;
    num[i] = 0;
  }

  double Re, Im;

  for (int i = -nk; i < nk; ++i) {
    for (int j = -nk; j < nk; ++j) {
      for (int k = -nk; k < nk; ++k) {
        CVector K;
        K.Set(thickness * i, thickness * j, thickness * k);

        Re = Im = 0.0;
        double sf;
        if (global) {
          for (int n = 0; n < NUM_PARTICLE; ++n) {
            double temp = K.Dot(particles[n]->center);
            Re += cos(temp);
            Im += sin(temp);
          }
          sf = (Re * Re + Im * Im) / NUM_PARTICLE;
        } else {
          int id = 0;
          for (int i = 0; i < species - 1; ++i) id += NpEach[i];

          for (int n = 0; n < NpEach[species]; ++n) {
            double temp = K.Dot(particles[id + n]->center);
            Re += cos(temp);
            Im += sin(temp);
          }
          sf = (Re * Re + Im * Im) / NpEach[species];
        }

        int index = int(K.Length() / thickness);
        Sf[index] += sf;
        num[index] += 1;
      }
    }
  }

  ofstream output("SF.txt");
  output << "species\t" << species << endl;
  output << global << endl;
  for (int i = 0; i < nk; ++i) {
    if (num[i] == 0) continue;
    Sf[i] /= num[i];
    output << (i + 0.5) * thickness << "\t" << Sf[i] << endl;
  }
  output.close();
}

void spectral_density(int species, bool global = true) {
  int nk = int(L / sqrt(3.0));  // a rough estimate
  double thickness = 2.0 * PI / L;
  double density[nk];

  for (int i = 0; i < nk; ++i) density[i] = 0;

  // _m: the Fourier transform of the indicator function
  double Re, Im, Re_m, Im_m;
  for (int i = -nk; i < nk; ++i) {
    for (int j = -nk; j < nk; ++j) {
      for (int k = -nk; k < nk; ++k) {
        CVector K;
        K.Set(thickness * i, thickness * j, thickness * k);

        Re = Im = 0.0;
        if (global) {
          for (int n = 0; n < NUM_PARTICLE; ++n) {
            double temp = -K.Dot(particles[n]->center);
            FFT_if(particles[n], 0.01, K, &Re_m, &Im_m);

            Re += cos(temp) * Re_m - sin(temp) * Im_m;
            Im += cos(temp) * Im_m + sin(temp) * Re_m;
          }
        } else {
          int id = 0;
          for (int i = 0; i < species - 1; ++i) id += NpEach[i];

          for (int n = 0; n < NpEach[species]; ++n) {
            double temp = -K.Dot(particles[id + n]->center);
            FFT_if(particles[id + n], 0.01, K, &Re_m, &Im_m);
            cout << n << endl;

            Re += cos(temp) * Re_m - sin(temp) * Im_m;
            Im += cos(temp) * Im_m + sin(temp) * Re_m;
          }
        }

        double sd = (Re * Re + Im * Im) / pow(L, 3.0);
        int index = int(K.Length() / thickness);

        density[index] += sd;
      }
    }
  }

  ofstream output("SD.txt");
  output << "species\t" << species << endl;
  output << global << endl;
  for (int i = 0; i < nk; ++i) {
    double Vshell =
        4.0 / 3.0 * PI * (3 * i * i + 3 * i + 1) * pow(thickness, 3.0);
    density[i] /= Vshell;

    output << (i + 0.5) * thickness << "\t" << density[i] << endl;
  }

  output.close();
}

// ----------------
// Voronoi Analysis
// ----------------
// Discretize surface for a given species under local coordinates
void getsurpoint(double para[], int num, CVector *SurfPot) {
  // Discretize a superellipsoidal particle surface along the longitude
  // direction into 20 (Yuan: NvA(s1,s2)) equal parts (40 (2NvA(s1,s2)) parts
  // for the latitude direction, Nv = 12)
  double p = para[0];

  int id;
  double dv = PI / (double)num;

  int NSuP = 2 * (num * num - num + 1);  // 40*20-20*2+2
  // SurfPot = new CVector[NSuP];

  double *sin_fai, *cos_fai, *sin_cita, *cos_cita;
  sin_cita = new double[num - 1];  // avoid replication
  cos_cita = new double[num - 1];
  sin_fai = new double[2 * num];
  cos_fai = new double[2 * num];

  // sin(cita)cos(fai), sin(cita)sin(fai), cos(cita)

  // v=cita (0, PI)
  // start from 1 to avoid replication
  for (int i = 0; i < num - 1; ++i) {
    double v = dv * double(i + 1);
    sin_cita[i] = sin(v);
    cos_cita[i] = cos(v);
  }
  // u=fai [0, 2*PI)
  for (int j = 0; j < 2 * num; ++j) {
    double u = dv * double(j);
    sin_fai[j] = sin(u);
    cos_fai[j] = cos(u);
  }
  SurfPot[0].Set(0.0, 0.0, para[3]);

  for (int i = 0; i < num - 1; ++i) {
    for (int j = 0; j < 2 * num; ++j) {
      id = i * 2 * num + j + 1;
      SurfPot[id].Set(para[1] * pow(sin_cita[i] * cos_fai[j], 1.0 / p),
                      para[2] * pow(sin_cita[i] * sin_fai[j], 1.0 / p),
                      para[3] * pow(cos_cita[i], 1.0 / p));
    }
  }

  id += 1;
  SurfPot[id].Set(0.0, 0.0, -para[3]);

  delete[] sin_cita;
  delete[] cos_cita;
  delete[] sin_fai;
  delete[] cos_fai;
}

// num: Precision
void OutputVoroPoint(int num, int replica) {
  char cha[100];

  sprintf(cha, "../vorofile/Surfpoint%d.txt", replica);
  ofstream scr(cha);
  scr << fixed << setprecision(10);

  CVector center, point_t, point;
  double A[3][3];

  int NSuP = 2 * (num * num - num + 1);  // 40*20-20*2+2

  int id = -1;
  for (int i = 0; i < Nkind; ++i) {
    CVector *SurfPot;
    SurfPot = new CVector[NSuP];
    getsurpoint(Para[i], num, SurfPot);

    for (int j = 0; j < NpEach[i]; ++j) {
      id++;
      center = particles[id]->center;
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n) {
          A[n][m] = particles[id]->e[m][n];
        }
      for (int t = 0; t < NSuP; ++t) {
        point_t = SurfPot[t];  // temp
        MatrixMult(A, point_t, point);
        point += center;
        PeriodicCheck(point);
        int index = id * NSuP + t;
        scr << index << "\t" << point[0] << "\t" << point[1] << "\t" << point[2]
            << endl;
      }
    }

    delete[] SurfPot;
  }
  scr.close();
}

void OutputBoundary(string filename, int replica) {
  ofstream output;
  output.open(filename, ios::app);

  output << replica << "\t" << L << endl;

  output.close();
}

void VoroVolume(int num, int replica) {
  char cha[100];
  sprintf(cha, "../vorofile/V_voro%d.txt", replica);
  ifstream voro(cha);
  if (!voro.good()) {
    cout << "Unable to open file!" << endl;
    exit(1);
  }

  sprintf(cha, "../vorofile/V_analysis%d.txt", replica);
  ofstream scr(cha);
  scr << fixed << setprecision(15);

  double *Vlocal, *PDlocal;
  double *ave_v, *ave_pd, *sigma_v, *sigma_pd;

  Vlocal = new double[NUM_PARTICLE];
  PDlocal = new double[NUM_PARTICLE];

  ave_v = new double[Nkind + 1];
  ave_pd = new double[Nkind + 1];
  sigma_v = new double[Nkind + 1];
  sigma_pd = new double[Nkind + 1];

  int NSuP = 2 * (num * num - num + 1);

  for (int i = 0; i < NUM_PARTICLE; ++i) Vlocal[i] = 0.0;

  int n, id;
  double x, y, z, v;
  voro >> n >> x >> y >> z >> v;

  while (1) {
    if (voro.eof()) break;
    // scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
    id = int(n / NSuP);
    Vlocal[id] += v;
    voro >> n >> x >> y >> z >> v;
  }

  scr << "index\tkind\tVorov\tPDlocal" << endl;
  id = -1;
  for (int i = 0; i < Nkind; ++i) {
    for (int j = 0; j < NpEach[i]; ++j) {
      id++;
      Vlocal[id] /= particles[id]->vol;
      PDlocal[id] = 1.0 / Vlocal[id];
      scr << id << "\t" << i << "\t" << Vlocal[id] << "\t" << PDlocal[id]
          << endl;
    }
  }

  id = -1;
  ave_v[Nkind] = ave_pd[Nkind] = 0.0;
  for (int i = 0; i < Nkind; ++i) {
    ave_v[i] = ave_pd[i] = 0.0;
    for (int j = 0; j < NpEach[i]; ++j) {
      id++;
      ave_v[i] += Vlocal[id];
      ave_pd[i] += PDlocal[id];
    }
    ave_v[Nkind] += ave_v[i];
    ave_pd[Nkind] += ave_pd[i];

    if (NpEach[i] > 0) {
      ave_v[i] /= double(NpEach[i]);
      ave_pd[i] /= double(NpEach[i]);
    }
  }
  ave_v[Nkind] /= double(NUM_PARTICLE);
  ave_pd[Nkind] /= double(NUM_PARTICLE);

  id = -1;
  sigma_v[Nkind] = sigma_pd[Nkind] = 0.0;
  for (int i = 0; i < Nkind; ++i) {
    sigma_v[i] = sigma_pd[i] = 0.0;
    for (int j = 0; j < NpEach[i]; ++j) {
      id++;
      sigma_v[Nkind] += pow(Vlocal[id] - ave_v[Nkind], 2.0);
      sigma_pd[Nkind] += pow(PDlocal[id] - ave_pd[Nkind], 2.0);
      sigma_v[i] += pow(Vlocal[id] - ave_v[i], 2.0);
      sigma_pd[i] += pow(PDlocal[id] - ave_pd[i], 2.0);
    }
    if (NpEach[i] > 0) {
      sigma_v[i] = sqrt(sigma_v[i] / double(NpEach[i] - 1));
      sigma_pd[i] = sqrt(sigma_pd[i] / double(NpEach[i] - 1));
    }
  }
  sigma_v[Nkind] = sqrt(sigma_v[Nkind] / double(NUM_PARTICLE - 1));
  sigma_pd[Nkind] = sqrt(sigma_pd[Nkind] / double(NUM_PARTICLE - 1));

  scr << "Sum:" << endl;
  scr << ave_v[Nkind] << "\t" << sigma_v[Nkind] << "\t" << ave_pd[Nkind] << "\t"
      << sigma_pd[Nkind] << endl;

  for (int i = 0; i < Nkind; ++i) {
    scr << NpEach[i] << "\t" << ave_v[i] << "\t" << sigma_v[i] << "\t"
        << ave_pd[i] << "\t" << sigma_pd[i] << endl;
  }

  voro.close();
  scr.close();
  delete[] Vlocal;
  delete[] PDlocal;
  delete[] ave_v;
  delete[] ave_pd;
  delete[] sigma_v;
  delete[] sigma_pd;
}

// ----------------
// Other utils
// ----------------
void PeriodicCheck(CVector &point) {
  for (int i = 0; i < 3; ++i) {
    while (point[i] >= L) {
      point[i] -= L;
    }
    while (point[i] < 0) {
      point[i] += L;
    }
  }
}

void POV_superball(int num) {
  // green, blue
  double color[3][3] = {
      {0.484, 0.984, 0}, {0, 0.75, 1.0}, {0.574, 0.4375, 0.855}};
  // p < 1.0 red p>1.0 green
  double len = 1.5 * L;

  double cell[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) cell[i][j] = 0;
  cell[0][0] = cell[1][1] = cell[2][2] = L;

  CSuperball *ps;
  char name[40];
  sprintf(name, "polysuperball%d.pov", num);
  ofstream pov(name);
  pov << "#include \"colors.inc\"" << endl;
  pov << "#include \"textures.inc\"" << endl;
  pov << "camera{ location <" << len << "," << len << "," << len
      << "> look_at <0,0,0> }" << endl;
  pov << "background {rgb 1}" << endl;
  pov << "light_source{ <60,60,10> color rgb <1,1,1> }" << endl;
  pov << "light_source{ <10,60,60> color rgb <1,1,1> }" << endl;

  double br = 0.05;
  CVector p_sc;
  CVector cs_, se_, cs, se;
  int bound_flag[3];
  cs.SetZero();
  for (int i = 0; i < 3; i++) {
    bound_flag[0] = bound_flag[1] = bound_flag[2] = 0;
    bound_flag[i] = 1;
    se_.Set(bound_flag[0], bound_flag[1], bound_flag[2]);
    MatrixMult(cell, se_, se);
    pov << "cylinder { <0,0,0> , <" << se[0] << "," << se[1] << "," << se[2]
        << "> ," << br << " texture{ pigment{\t";
    pov << "color rgb <0.625,0.125,0.9375> } finish{ phong 1.0 reflection 0.0 "
           "} } }"
        << endl;  // purple
  }
  cs_.Set(1, 1, 1);
  MatrixMult(cell, cs_, cs);
  for (int i = 0; i < 3; i++) {
    bound_flag[0] = bound_flag[1] = bound_flag[2] = 1;
    bound_flag[i] = 0;
    se_.Set(bound_flag[0], bound_flag[1], bound_flag[2]);
    MatrixMult(cell, se_, se);
    pov << "cylinder { <" << cs[0] << "," << cs[1] << "," << cs[2] << "> , <"
        << se[0] << "," << se[1] << "," << se[2] << "> ," << br
        << " texture{ pigment{\t";
    pov << "color rgb <0.625,0.125,0.9375> } finish{ phong 1.0 reflection 0.0 "
           "} } }"
        << endl;  // purple
  }
  for (int i = 0; i < 3; i++) {
    bound_flag[0] = bound_flag[1] = bound_flag[2] = 0;
    bound_flag[i] = 1;
    cs_.Set(bound_flag[0], bound_flag[1], bound_flag[2]);
    MatrixMult(cell, cs_, cs);
    for (int j = 0; j < 3; j++) {
      if (j == i) continue;
      bound_flag[0] = bound_flag[1] = bound_flag[2] = 1;
      bound_flag[j] = 0;
      se_.Set(bound_flag[0], bound_flag[1], bound_flag[2]);
      MatrixMult(cell, se_, se);
      pov << "cylinder { <" << cs[0] << "," << cs[1] << "," << cs[2] << "> , <"
          << se[0] << "," << se[1] << "," << se[2] << "> ," << br
          << " texture{ pigment{\t";
      pov << "color rgb <0.625,0.125,0.9375> } finish{ phong 1.0 reflection "
             "0.0 } } }"
          << endl;  // purple
    }
  }

  int id = -1;
  for (int m = 0; m < Nkind; ++m) {
    double c[3];
    c[0] = color[m][0];
    c[1] = color[m][1];
    c[2] = color[m][2];

    for (int i = 0; i < NpEach[m]; i++) {
      id++;
      ps = particles[id];

      pov << "object { superellipsoid { <" << 1.0 / ps->p << "," << 1.0 / ps->p
          << "> scale <" << ps->r_scale[0] << "," << ps->r_scale[1] << ","
          << ps->r_scale[2] << "> texture{ pigment{\t";
      pov << "color rgb <" << c[0] << "," << c[1] << "," << c[2]
          << "> } finish{ phong 1.0 reflection 0.0 } } } matrix<";
      for (int j = 0; j < 3; j++)
        pov << ps->e[j][0] << "," << ps->e[j][1] << "," << ps->e[j][2] << ",";
      pov << ps->center[0] << "," << ps->center[1] << "," << ps->center[2]
          << "> }" << endl;
    }
  }
  pov.close();
}
