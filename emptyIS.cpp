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

	size_t seed = 2;
	Generator G(seed);
	double R = 0.7;
	double a_s = 0.105/2.0;
	double a_p = 0.200/2.0;
	double b_s = 0.125/2.0;
	double b_p = 0.225/2.0;
	double r_e = 2.0*b_p;
	double alpha = 1.2786;
	double rho = 0.99;
	double n = 1.0;
	double NA = 0.5;
	double cos_theta0 = sqrt(1.0 - (NA*NA)/(n*n));
	double z_s = sqrt(R*R - r_e*r_e);
	size_t collision_limit = 50000;
	double start = R*sqrt(1.0/3.0);
	size_t N_photons = 3e6;
	double N_d = (double) N_photons;
	size_t bins = 1000;

	vector<Photon> photons(N_photons);
	initPhotons(photons, R, z_s, cos_theta0, a_s, b_s, G);

	for(size_t i=0; i<N_photons; i++){
		//traceSingleDiffuse(photons[i], R, rho, collision_limit, G);
		//traceSingleDiffuseAbsorption(photons[i], R, rho, alpha, collision_limit, G);
		//traceSingleDiffuseInvariant(photons[i], R, rho, collision_limit, G);
		tracePhotonEmptyIS(photons[i], R, rho, z_s, cos_theta0, a_p, b_p, collision_limit, G);
	}

	double collected = 0.0;
	double absorbedEntrance = 0.0;
	double absorbedWall = 0.0;
	for(size_t i=0; i<N_photons; i++){
		if(photons[i].isCollected){
			collected += 1.0;
		}else if(photons[i].isAbsorbedEntrance){
			absorbedEntrance += 1.0;
		}else if(photons[i].isAbsorbedWall){ 
			absorbedWall += 1.0;
		}
	}

	cout << "entrance fraction: " << (r_e*r_e/(4.0*R*R)) << endl;
	cout << "A_core/A_entrance: " << (a_p*a_p/(r_e*r_e)) << endl;
	cout << "Rho->inf: " << (a_p*a_p/(r_e*r_e))*(1.0-cos_theta0) << endl;
	cout << "Collected: " << collected/N_d << endl;
	cout << "entrance: " << absorbedEntrance/N_d << endl;
	cout << "wall: " << absorbedWall/N_d << endl;
	cout << "sum: " << collected/N_d + absorbedEntrance/N_d + absorbedWall/N_d << endl;
	

	double Ps = psTheory(alpha,R);
	//double meanCollision = meanNumberCollision(photons);
	//cout << "mean collision: " << meanCollision << endl;
	//cout << "mean collision theoretical: " << (Ps/(1.0-Ps*rho)) << endl;
	//cout << "stddev collision: " << stddevNumberCollision(photons, meanCollision) << endl; 
	double meanPath = meanPathLengthCollected(photons);
	cout << "mean path: " << meanPath << endl;
	//cout << "mean path theoretical: " << ((1.0/alpha)*((1.0-Ps)/(1.0-rho*Ps))) << endl;
	cout << "stddev path: " << stddevPathLengthCollected(photons, meanPath) << endl;
	return 0;
}
	
