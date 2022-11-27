#include "analysis.h"
using namespace std;

int NUM_PARTICLE;
double L;  // Box length

double sumvol;
double PD;

vector<CSuperball *> particles;

ifstream myfile;

void readfile() {
  int i, j;
  string line;
  double nonesize;

  myfile >> NUM_PARTICLE;
  myfile >> L;

  double p, r_scale[3];
  sumvol = 0.0;
  for (int i = 0; i < NUM_PARTICLE; i++) {
    myfile >> p >> r_scale[0] >> r_scale[1] >> r_scale[2];
    particles.push_back(new CSuperball(p, r_scale[0], r_scale[1], r_scale[2]));
    particles[i]->ID = i;

    myfile >> particles[i]->center[0] >> particles[i]->center[1] >>
        particles[i]->center[2];

    for (int j = 0; j < 3; j++)
      myfile >> particles[i]->e[j][0] >> particles[i]->e[j][1] >>
          particles[i]->e[j][2];

    sumvol += particles[i]->vol;
  }

  PD = sumvol / L / L / L;
}
