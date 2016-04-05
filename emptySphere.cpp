#include "Coord.h"
#include "Dir.h"
#include "Photon.h"
#include "RandomGenerator.h"
#include "functions.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>


using namespace photon_IS;
using namespace std;

int main(){

	Generator G(1);
	double R = 1.0;
	double rho = 0.95;
	double p = 1e-6;
	size_t collision_limit = 900;
	double max_s = 1000.0;
	double start = sqrt(R/3.0);
	size_t N_photons = 1e7;
	size_t bins = 1000;

	vector<Photon> photons(N_photons);

	for(size_t i=0; i<N_photons; i++){
		photons[i].r.x = start;
		photons[i].r.y = start;
		photons[i].r.z = start;
		diffuseLocal(photons[i].u, G);
		photons[i].u = convertToGlobal(photons[i].u, photons[i].r, R);
		traceSingleDiffuse(photons[i], R, rho, collision_limit, G);
	}
	cout << "finished trace" << endl;

	histCollisions(photons, "collisionDistr.dat");
	histPathLength(photons, "pathLengthDistr1.dat", 0.0, max_s, bins);

	double meanCollision = meanNumberCollision(photons);
	cout << "mean collision: " << meanCollision << endl;
	cout << "stddev collision: " << stddevNumberCollision(photons, meanCollision) << endl; 
	double meanPath = meanPathLength(photons);
	cout << "mean path: " << meanPath << endl;
	cout << "stddev path: " << stddevPathLength(photons, meanPath) << endl;

	{
		using namespace arma;
		vector<double> cols(N_photons);
		vector<double> path(N_photons);
		for(size_t i=0; i<N_photons; i++){
			cols[i] = (double) photons[i].N_w;
			path[i] = photons[i].s_tot;
		}
		vec collisions = conv_to< vec >::from(cols);
		vec paths = conv_to< vec >::from(path);
		uvec h_collisions = hist(collisions, linspace<vec>(1,200,200));
		h_collisions.save("h_collisions.dat", raw_ascii);
		uvec h_paths = hist(paths,1000);
		h_paths.save("h_paths.dat", raw_ascii);
	}
	return 0;
}	
