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
          cout<<n<<endl;

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