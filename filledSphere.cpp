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

	size_t seed = 100;
	Generator G(seed);
	double R = 1.0;
	double a_s = 0.105/2.0;
	double a_p = 0.125/2.0;
	double b_s = 0.105/2.0;
	double b_p = 0.125/2.0;
	double D_e = 2.1*b_p;
	double alpha = 1.2786;
	double rho = 0.95;
	double n = 1.0;
	double NA = 0.22;
	double cos_theta0 = sqrt(1.0 - (NA*NA)/(n*n));
	double z_s = sqrt(R*R - 0.25*D_e*D_e);
	size_t collision_limit = 5000;
	double start = R*sqrt(1.0/3.0);
	size_t N_photons = 1e6;
	size_t bins = 1000;

	vector<Photon> photons(N_photons);
	initPhotons(photons, R, z_s, cos_theta0, a_s, b_s, G);

	for(size_t i=0; i<N_photons; i++){
		photons[i].r.x = start;
		photons[i].r.y = start;
		photons[i].r.z = start;
		diffuseLocal(photons[i].u, G);
		photons[i].u = convertToGlobal(photons[i].u, photons[i].r, R);
		//traceSingleDiffuse(photons[i], R, rho, collision_limit, G);
		traceSingleDiffuseAbsorption(photons[i], R, rho, alpha, collision_limit, G);
		//traceSingleDiffuseInvariant(photons[i], R, rho, collision_limit, G);
	}
	cout << "finished trace" << endl;
	//savePhotons(photons, "data/filledSphere/filledSphere_Alpha1.2786_R1.0_Rho0.95.dat", R, rho, N_photons, seed);

	//histCollisions(photons, "collisionDistr.dat");
	//histPathLength(photons, "pathLengthDistr1.dat", 0.0, max_s, bins);

	double Ps = psTheory(alpha,R);
	double meanCollision = meanNumberCollision(photons);
	cout << "mean collision: " << meanCollision << endl;
	cout << "mean collision theoretical: " << (Ps/(1.0-Ps*rho)) << endl;
	cout << "stddev collision: " << stddevNumberCollision(photons, meanCollision) << endl; 
	double meanPath = meanPathLength(photons);
	cout << "mean path: " << meanPath << endl;
	cout << "mean path theoretical: " << ((1.0/alpha)*((1.0-Ps)/(1.0-rho*Ps))) << endl;
	cout << "stddev path: " << stddevPathLength(photons, meanPath) << endl;
	return 0;
}
	
