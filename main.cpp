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
	/*
	Coord c(0.1, 0.6, 123.2);
	Coord d(c);
	Coord e;
	e = d;
	std::cout << c.x << " " << c.y << " " << c.z << std::endl;
	std::cout << "second\n" << d.x << " " << d.y << " " << d.z << std::endl;
	std::cout << "third\n" << e.x << " " << e.y << " " << e.z << std::endl;
	
	Dir n;
	std::cout << n.u_x << " " << n.u_y << " " << n.u_z << std::endl;


	Photon p;
	Photon p2(e,n);

	std::cout << p2.u.u_z << "\n";

	Generator G(3);
	std::cout << G.uniform() << "\n";
	std::cout << G.count << "\n";

	diffuseLocal(n,G);
	std::cout << n.u_x << " " << n.u_y << " " << n.u_z << std::endl;
	*/

	Generator G(3);
	double R = 1.0;
	double rho = 0.95;
	size_t collision_limit = 500;
	double start = sqrt(R/3.0);
	Coord r0(start,start,start);
	Dir n0;
	Dir n1;
	size_t N = 100;
	vector<double> X(N);
	vector<double> Y(N);
	vector<double> Z(N);
	X[0] = sqrt(1.0/3.0);
	Y[0] = sqrt(1.0/3.0);
	Z[0] = sqrt(1.0/3.0);

	ofstream outfile("sphereReflections.dat");
	outfile << X[0] << " " << Y[0] << " " << Z[0] << "\n";

	for(int i=1; i<N; i++){
		diffuseLocal(n0,G);
		n1 = convertToGlobal(n0,r0,R);

		cout << "n0 squared: " << (n0.u_x*n0.u_x + n0.u_y*n0.u_y + n0.u_z*n0.u_z) << endl;
		std::cout << "n0: " << n0.u_x << " " << n0.u_y << " " << n0.u_z << std::endl;
		cout << "n1 squared: " << (n1.u_x*n1.u_x + n1.u_y*n1.u_y + n1.u_z*n1.u_z) << endl;
		std::cout << "n1: " << std::setprecision(10) <<   n1.u_x << " " << n1.u_y << " " << n1.u_z << std::endl;

		double s0 = - 2.0*(n1.u_x*r0.x + n1.u_y*r0.y + n1.u_z*r0.z);
		cout << "s0: " << s0 << endl;
		r0.x = r0.x + s0*n1.u_x;
		r0.y = r0.y + s0*n1.u_y;
		r0.z = r0.z + s0*n1.u_z;
		X[i] = r0.x;
		Y[i] = r0.y;
		Z[i] = r0.z;
		cout << "R squared:" << std::setprecision(10) << (X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i]) << endl;
		outfile << X[i] << " " << Y[i] << " " << Z[i] << "\n";

	}
	outfile.close();


	//TESTING SIMPLE SPHERE. Diffuse reflection. Tracing actual paths, not really necessary, but gives indiciation of compute time.
	Photon p;
	p.r.x = start;
	p.r.y = start;
	p.r.z = start;

	diffuseLocal(p.u, G);
	p.u = convertToGlobal(p.u, p.r, R);
	traceSingleDiffuse(p, R, rho, collision_limit, G);
	cout << "photon: " << p.s_tot << endl;

	size_t N_photons = 1e6;
	vector<Photon> photons(N_photons);
	ofstream colFile("pos.dat");
	for(size_t i=0; i<N_photons; i++){
		photons[i].r.x = start;
		photons[i].r.y = start;
		photons[i].r.z = start;
		diffuseLocal(photons[i].u, G);
		photons[i].u = convertToGlobal(photons[i].u, photons[i].r, R);
		traceSingleDiffuse(photons[i], R, rho, collision_limit, G);
	}

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
		uvec h_collisions = hist(collisions, 100);
		h_collisions.save("h_collisions.dat", raw_ascii);
		uvec h_paths = hist(paths,1000);
		h_paths.save("h_paths.dat", raw_ascii);
	}
		
		






	














	
















	return 0;
}
