#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Coord.h"
#include "Dir.h"
#include "Photon.h"
#include "RandomGenerator.h"
#include <vector>

using namespace photon_IS;


void diffuseLocal(Dir &u, Generator &G);

Dir convertToGlobal(Dir &u, Coord &r, double &R);

void traceSingleDiffuse(Photon &p, double &R, double &rho, size_t &lim, Generator &G);

void traceCollectionDiffuse(std::vector<Photon> &P, double &R, double &rho, size_t &lim, Generator &G);

double meanNumberCollision(std::vector<Photon> &P);

double stddevNumberCollision(std::vector<Photon> &P, double &mean);

double meanPathLength(std::vector<Photon> &P);

double stddevPathLength(std::vector<Photon> &P, double &mean);

#endif
