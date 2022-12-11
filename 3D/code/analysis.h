#include <iostream>
#include <istream>
#include <sstream>
#include <vector>

#include "Voxel.h"
#include "fstream"
#include "iomanip"
#include "string"

using namespace std;

void readfile(string filename);

void askforspace();
void releasespace();

void spectral_density();


void getsurpoint(double para[], int num, CVector *SurfPot);
void OutputVoroPoint(int num, int replica);
void VoroVolume(int num, int replica);
void OutputBoundary(string filename, int replica);

void PeriodicCheck(CVector &point);

void POV_superball(int replica);


