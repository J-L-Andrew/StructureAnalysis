#include "analysis.h"

int NUM_PARTICLE;
double L;  // Box length

double sumvol;
double PD;

vector<CSuperball *> particles;

ifstream myfile;

int main() {
  string file = "config0.pac";
  readfile(file);

  spectral_density();
}

void readfile(string filename) {
  myfile.open(filename);

  if (!myfile.is_open())  // unable to open file
  {
    cout << "Unable to open myfile" << endl;
    exit(1);
  }

  myfile >> NUM_PARTICLE;
  myfile >> L;

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

void spectral_density() {
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
        for (int n = 0; n < NUM_PARTICLE; ++n) {
          double temp = -K.Dot(particles[n]->center);
          FFT_if(particles[n], 0.01, K, &Re_m, &Im_m);
          cout << n << endl;

          Re += cos(temp) * Re_m - sin(temp) * Im_m;
          Im += cos(temp) * Im_m + sin(temp) * Re_m;
        }

        double sd = (Re * Re + Im * Im) / pow(L, 3.0);
        int id = int(K.Length() / thickness);

        density[id] += sd;
      }
    }
  }
  ofstream output("SD.txt");
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
void getsurpoint(double p, double r_scale[], CVector *SurfPot) {
  // Discretize a superellipsoidal particle surface along the longitude
  // direction into 20 (Yuan: NvA(s1,s2)) equal parts (40 (2NvA(s1,s2)) parts
  // for the latitude direction, Nv = 12)

  int num = 20, id;
  double dv = PI / (double)num;

  int NSuP = 2 * (num * num - num + 1);  // 40*20-20*2+2
  SurfPot = new CVector[NSuP];

  double *sin_fai, *cos_fai, *sin_cita, *cos_cita;
  sin_cita = new double[num - 1];  // avoid replication
  cos_cita = new double[num - 1];
  sin_fai = new double[2 * num];
  cos_fai = new double[2 * num];

  // sin(cita)cos(fai), sin(cita)sin(fai), cos(cita)

  // v=cita (0, PI)
  // start from 1 to avoid replication
  for (int i = 0; i < num - 1; ++i) {
    double v = dv * double(i);
    sin_cita[i] = sin(v);
    cos_cita[i] = cos(v);
  }
  // u=fai [0, 2*PI)
  for (int j = 0; j < 2 * num; ++j) {
    double u = dv * double(j);
    sin_fai[j] = sin(u);
    cos_fai[j] = cos(u);
  }
  SurfPot[0].Set(0.0, 0.0, r_scale[2]);

  for (int i = 0; i < num - 1; ++i) {
    for (int j = 0; j < 2 * num; ++j) {
      id = i * 2 * num + j + 1;
      SurfPot[id].Set(r_scale[0] * pow(sin_cita[i] * cos_fai[j], 1.0 / p),
                      r_scale[1] * pow(sin_cita[i] * sin_fai[j], 1.0 / p),
                      r_scale[2] * pow(cos_cita[i], 1.0 / p));
    }
  }

  id += 1;
  SurfPot[id].Set(0.0, 0.0, -r_scale[2]);

  delete[] sin_cita;
  delete[] cos_cita;
  delete[] sin_fai;
  delete[] cos_fai;

  /*cout << numi << "\t" << NSuP << endl;
  double radius = particle_id*0.01;
  ofstream scr("surfacepoint.scr");
  scr << fixed << setprecision(15);
  for (i = 0; i <NSuP; i++)
  {
  if (radius > 0)
  {
  scr << "sphere " << SurfPot[i][0] << "," << SurfPot[i][1] << "," <<
  SurfPot[i][2] << " " << radius << endl;
  }
  }
  scr << "grid off\nVSCURRENT c\nzoom e\n";
  scr.close();*/
}

void OutputVoroPoint(int num, int replica)  //输出用于voronoi剖分的表面离散点
{
  char cha[100];
  sprintf(cha, "surfpoint%d.txt", replica);
  ofstream scr(cha);
  scr << fixed << setprecision(10);

  int i, j, m, n, ii, jj;
  CVector surfpl0, surfpe;
  CVector center;
  double A[3][3];

  //剖分相关
  int numi;  //剖分精度
  double dv, v, u, cos1, sin1, cos2, sin2;
  double *pcos1, *psin1, *pcos2, *psin2;
  double *sgncos1, *sgnsin1, *sgnsin2;

  dv = PI / (double)num;
  NSuP = 2 * (num * num - num + 1);
  cout << "num and NSuP:\t" << num << "\t" << NSuP << endl;
  SurfPot = new CVector[NSuP];
  pcos1 = new double[2 * num];
  psin1 = new double[2 * num];
  pcos2 = new double[num - 1];
  psin2 = new double[num - 1];
  sgncos1 = new double[2 * num];
  sgnsin1 = new double[2 * num];
  sgnsin2 = new double[num - 1];

  m = -1;
  for (ii = 0; ii < NKind; ii++) {
    //表面点离散//v,cita (-0.5*PI,0.5*PI)
    for (i = 1; i < num; i++) {
      v = -0.5 * PI + dv * i;
      cos2 = cos(v);
      sin2 = sin(v);
      pcos2[i - 1] = pow(fabs(cos2), 1.0 / PPara[ii][4]);
      psin2[i - 1] = pow(fabs(sin2), 1.0 / PPara[ii][4]);
      sgnsin2[i - 1] = sgn(sin2);
    }
    // u,fai [-PI,PI)
    for (j = 0; j < 2 * num; j++) {
      u = -PI + dv * j;
      cos1 = cos(u);
      sin1 = sin(u);
      pcos1[j] = pow(fabs(cos1), 1.0 / PPara[ii][3]);
      psin1[j] = pow(fabs(sin1), 1.0 / PPara[ii][3]);
      sgncos1[j] = sgn(cos1);
      sgnsin1[j] = sgn(sin1);
    }
    SurfPot[0][0] = 0.0;
    SurfPot[0][1] = 0.0;
    SurfPot[0][2] = -PPara[ii][2];
    numi = 1;
    for (i = 1; i < num; i++)
      for (j = 0; j < 2 * num; j++) {
        SurfPot[numi][0] = PPara[ii][0] * sgncos1[j] * pcos1[j] * pcos2[i - 1];
        SurfPot[numi][1] = PPara[ii][1] * sgnsin1[j] * psin1[j] * pcos2[i - 1];
        SurfPot[numi][2] = PPara[ii][2] * sgnsin2[i - 1] * psin2[i - 1];
        numi++;
      }
    SurfPot[numi][0] = 0.0;
    SurfPot[numi][1] = 0.0;
    SurfPot[numi][2] = PPara[ii][2];
    //表面点输出
    for (jj = 0; jj < NpEach[ii]; jj++) {
      m++;
      center = polyhedra[m]->center;
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
          A[j][i] = polyhedra[m]->e[i][j];  // a局部坐标的旋转矩阵
        }
      for (i = 0; i < NSuP; i++) {
        surfpl0 = SurfPot[i];
        MatrixMult(A, surfpl0, surfpe);
        surfpe = surfpe + center;
        PeriodCheckP(surfpe);
        n = m * NSuP + i;
        scr << n << "\t" << surfpe[0] << "\t" << surfpe[1] << "\t" << surfpe[2]
            << endl;
      }
    }
  }
  scr.close();
  delete[] pcos1;
  delete[] psin1;
  delete[] pcos2;
  delete[] psin2;
  delete[] sgncos1;
  delete[] sgnsin1;
  delete[] sgnsin2;
  delete[] SurfPot;
}

int VoroVolume000(int num)  //统计vorovolume
{
  char cha[100];
  sprintf_s(cha, "vorovol %d.txt", Counterfile);
  ifstream voro;
  voro.open(cha);
  if (!voro.good()) {
    cout << cha << "不存在" << endl;
    return 1;
  }

  sprintf_s(cha, "vororesults %d.txt", Counterfile);
  ofstream scr(cha);
  scr << fixed << setprecision(15);

  int i, n;
  double x, y, z, v;
  double *Vlocal, *PDlocal;
  double avVlocal, avPDlocal, detaVlocal, detaPDlocal;

  NSuP = 2 * (num * num - num + 1);
  PDlocal = new double[Npolyhedron];
  Vlocal = new double[Npolyhedron];
  for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
  voro >> n >> x >> y >> z >> v;
  while (1) {
    if (voro.eof()) break;
    // scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
    i = int(n / NSuP);
    Vlocal[i] += v;
    voro >> n >> x >> y >> z >> v;
  }

  avVlocal = 0.0;
  avPDlocal = 0.0;
  scr << "i\tVorov\tPDlocal" << endl;
  for (i = 0; i < Npolyhedron; i++) {
    Vlocal[i] /= polyhedra[i]->volume;
    avVlocal += Vlocal[i];
    PDlocal[i] = 1.0 / Vlocal[i];
    avPDlocal += PDlocal[i];
    scr << i << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
  }

  avVlocal /= double(Npolyhedron);
  avPDlocal /= double(Npolyhedron);

  detaVlocal = 0.0;
  detaPDlocal = 0.0;
  for (i = 0; i < Npolyhedron; i++) {
    v = Vlocal[i] - avVlocal;
    detaVlocal += v * v;
    v = PDlocal[i] - avPDlocal;
    detaPDlocal += v * v;
  }
  detaVlocal = sqrt(detaVlocal / double(Npolyhedron));
  detaPDlocal = sqrt(detaPDlocal / double(Npolyhedron));

  result << avVlocal << "\t" << avPDlocal << "\t" << detaVlocal << "\t"
         << detaPDlocal << "\t";
  voro.close();
  scr.close();
  delete[] Vlocal;
  delete[] PDlocal;
  return 0;
}
int VoroVolume(int num)  //统计vorovolume
{
  char cha[100];
  sprintf_s(cha, "vorovol %d.txt", Counterfile);
  ifstream voro;
  voro.open(cha);
  if (!voro.good()) {
    cout << cha << "不存在" << endl;
    return 1;
  }

  sprintf_s(cha, "vororesults %d.txt", Counterfile);
  ofstream scr(cha);
  scr << fixed << setprecision(15);

  int i, n, ii, jj;
  double x, y, z, v;
  double *Vlocal, *PDlocal;
  double *avVlocal, *avPDlocal, *detaVlocal, *detaPDlocal;
  avVlocal = new double[NKind + 1];
  avPDlocal = new double[NKind + 1];
  detaVlocal = new double[NKind + 1];
  detaPDlocal = new double[NKind + 1];

  NSuP = 2 * (num * num - num + 1);
  cout << "num and NSuP:\t" << num << "\t" << NSuP << endl;
  result << num << "\t" << NSuP << "\t";
  PDlocal = new double[Npolyhedron];
  Vlocal = new double[Npolyhedron];
  for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
  voro >> n >> x >> y >> z >> v;
  while (1) {
    if (voro.eof()) break;
    // scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
    i = int(n / NSuP);
    Vlocal[i] += v;
    voro >> n >> x >> y >> z >> v;
  }
  scr << "i\tkind\tVorov\tPDlocal" << endl;
  for (i = 0; i < Npolyhedron; i++) {
    Vlocal[i] /= polyhedra[i]->volume;
    PDlocal[i] = 1.0 / Vlocal[i];
    scr << i << "\t" << polyhedra[i]->kind << "\t" << Vlocal[i] << "\t"
        << PDlocal[i] << endl;
  }

  i = -1;
  avVlocal[NKind] = 0.0;
  avPDlocal[NKind] = 0.0;
  for (ii = 0; ii < NKind; ii++)  //种类数
  {
    avVlocal[ii] = 0.0;
    avPDlocal[ii] = 0.0;
    for (jj = 0; jj < NpEach[ii]; jj++)  //每种的个数
    {
      i++;
      avVlocal[ii] += Vlocal[i];
      avPDlocal[ii] += PDlocal[i];
    }
    avVlocal[NKind] += avVlocal[ii];
    avPDlocal[NKind] += avPDlocal[ii];
    if (NpEach[ii] > 0) {
      avVlocal[ii] /= double(NpEach[ii]);
      avPDlocal[ii] /= double(NpEach[ii]);
    }
  }
  avVlocal[NKind] /= double(Npolyhedron);
  avPDlocal[NKind] /= double(Npolyhedron);

  i = -1;
  detaVlocal[NKind] = 0.0;
  detaPDlocal[NKind] = 0.0;
  for (ii = 0; ii < NKind; ii++)  //种类数
  {
    detaVlocal[ii] = 0.0;
    detaPDlocal[ii] = 0.0;
    for (jj = 0; jj < NpEach[ii]; jj++)  //每种的个数
    {
      i++;
      v = Vlocal[i] - avVlocal[NKind];
      detaVlocal[NKind] += v * v;
      v = PDlocal[i] - avPDlocal[NKind];
      detaPDlocal[NKind] += v * v;
      v = Vlocal[i] - avVlocal[ii];
      detaVlocal[ii] += v * v;
      v = PDlocal[i] - avPDlocal[ii];
      detaPDlocal[ii] += v * v;
    }
    if (NpEach[ii] > 0) {
      detaVlocal[ii] = sqrt(detaVlocal[ii] / double(NpEach[ii]));
      detaPDlocal[ii] = sqrt(detaPDlocal[ii] / double(NpEach[ii]));
    }
  }
  detaVlocal[NKind] = sqrt(detaVlocal[NKind] / double(Npolyhedron));
  detaPDlocal[NKind] = sqrt(detaPDlocal[NKind] / double(Npolyhedron));

  result << avVlocal[NKind] << "\t" << detaVlocal[NKind] << "\t"
         << avPDlocal[NKind] << "\t" << detaPDlocal[NKind] << "\t";
  for (ii = 0; ii < NKind; ii++) result << NpEach[ii] << "\t";
  for (ii = 0; ii < NKind; ii++)
    result << avVlocal[ii] << "\t" << detaVlocal[ii] << "\t";
  for (ii = 0; ii < NKind; ii++)
    result << avPDlocal[ii] << "\t" << detaPDlocal[ii] << "\t";

  voro.close();
  scr.close();
  delete[] Vlocal;
  delete[] PDlocal;
  delete[] avVlocal;
  delete[] avPDlocal;
  delete[] detaVlocal;
  delete[] detaPDlocal;
  return 0;
}
void OutputVoroPointC()  //输出用于voronoi剖分的表面离散点
{
  char cha[100];
  sprintf_s(cha, "surfpoint %d.txt", Counterfile);
  ofstream scr(cha);
  scr << fixed << setprecision(15);

  int m;
  CVector center;
  for (m = 0; m < Npolyhedron; m++) {
    center = polyhedra[m]->center;
    scr << m << "\t" << center[0] << "\t" << center[1] << "\t" << center[2]
        << endl;
  }
  scr.close();
}
int VoroVolumeC()  //统计vorovolume
{
  char cha[100];
  sprintf_s(cha, "vorovol %d.txt", Counterfile);
  ifstream voro;
  voro.open(cha);
  if (!voro.good()) {
    cout << cha << "不存在" << endl;
    return 1;
  }

  sprintf_s(cha, "vororesults %d.txt", Counterfile);
  ofstream scr(cha);
  scr << fixed << setprecision(15);

  int i, n;
  double x, y, z, v;
  double *Vlocal, *PDlocal;
  double avVlocal, avPDlocal, detaVlocal, detaPDlocal;

  PDlocal = new double[Npolyhedron];
  Vlocal = new double[Npolyhedron];
  for (i = 0; i < Npolyhedron; i++) Vlocal[i] = 0.0;
  voro >> n >> x >> y >> z >> v;
  while (1) {
    if (voro.eof()) break;
    // scr << n << "\t" << x << "\t" << y << "\t" << z << "\t" << v << endl;
    // i = int(n / NSuP);
    i = n;
    Vlocal[i] += v;
    voro >> n >> x >> y >> z >> v;
  }

  avVlocal = 0.0;
  avPDlocal = 0.0;
  scr << "i\tVorov\tPDlocal" << endl;
  for (i = 0; i < Npolyhedron; i++) {
    Vlocal[i] /= polyhedra[i]->volume;
    avVlocal += Vlocal[i];
    PDlocal[i] = 1.0 / Vlocal[i];
    avPDlocal += PDlocal[i];
    scr << i << "\t" << Vlocal[i] << "\t" << PDlocal[i] << endl;
  }

  avVlocal /= double(Npolyhedron);
  avPDlocal /= double(Npolyhedron);

  detaVlocal = 0.0;
  detaPDlocal = 0.0;
  for (i = 0; i < Npolyhedron; i++) {
    v = Vlocal[i] - avVlocal;
    detaVlocal += v * v;
    v = PDlocal[i] - avPDlocal;
    detaPDlocal += v * v;
  }
  detaVlocal = sqrt(detaVlocal / double(Npolyhedron));
  detaPDlocal = sqrt(detaPDlocal / double(Npolyhedron));

  result << avVlocal << "\t" << avPDlocal << "\t" << detaVlocal << "\t"
         << detaPDlocal << "\t";
  voro.close();
  scr.close();
  delete[] Vlocal;
  delete[] PDlocal;
  return 0;
}
