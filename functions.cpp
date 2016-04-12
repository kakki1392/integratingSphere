#include "functions.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <armadillo>
#include <fstream>
#include <iomanip>

using namespace photon_IS;
using namespace std;

void diffuseLocal(Dir &u, Generator &G){
	double chi_1 = G.uniform();
	double chi_2 = G.uniform();
	double c = std::cos(2.0*M_PI*chi_1);
	double chi_2_sqrt = std::sqrt(chi_2);
	u.u_x = std::cos(2.0*M_PI*chi_1)*chi_2_sqrt;
	u.u_y = std::sin(2.0*M_PI*chi_1)*chi_2_sqrt; //make faster with sqrt and sign
	u.u_z = std::sqrt(1.0-chi_2);
}


Dir convertToGlobal(Dir &u, Coord &r, double &R){
	Dir u_global;
	double cos_theta = r.z/R;
	double sin_theta = std::sqrt(1.0-cos_theta*cos_theta);
	if(sin_theta < 0.0001){
		std::cout << "Warning, sin_theta is small" << std::endl;
	}
	double cos_phi = r.x/(R*sin_theta);  //WATCH OUT FOR LIMITS FOR THETA
	double sin_phi = r.y/(R*sin_theta);
	u_global.u_x = -u.u_x*cos_phi*cos_theta - u.u_y*sin_phi - u.u_z*cos_phi*sin_theta;
	u_global.u_y = -u.u_x*sin_phi*cos_theta + u.u_y*cos_phi - u.u_z*sin_phi*sin_theta;
	u_global.u_z =  u.u_x*sin_theta                         - u.u_z*cos_theta;

	//Force normalization? I think necessary. Small error blows up.
	double norm = std::sqrt(u_global.u_x*u_global.u_x + u_global.u_y*u_global.u_y + u_global.u_z*u_global.u_z);
	u_global.u_x = u_global.u_x/norm;
	u_global.u_y = u_global.u_y/norm;
	u_global.u_z = u_global.u_z/norm;

	return u_global;
}

double drawPathLength(double &alpha, Generator &G){
	return (-log(G.uniform())/alpha);
}

double psTheory(double &alpha, double &R){
	return ((1.0/(2.0*alpha*alpha*R*R))*(1.0 - (1.0 + 2.0*alpha*R)*exp(-2.0*alpha*R)));
}

	

void traceSingleDiffuse(Photon &p, double &R, double &rho, size_t &lim, Generator &G){
	double s0 = 0.0;
	//Move from initial position and direction
	s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
	p.s_tot += s0;
	p.N_w ++;

	while(p.N_w < lim){
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		if(G.uniform() < rho){
			p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			p.s_tot += s0;
		}
		else{
			break;
		}
	}
}

void traceSingleDiffuseAbsorption(Photon &p, double &R, double &rho, double &alpha, size_t &lim, Generator &G){
	double s0 = 0.0;
	double path = drawPathLength(alpha, G);
	//Move from initial position and direction
	s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
	if(s0 < path){
		p.s_tot += s0;
		p.N_w ++;
		path -= s0;
	}
	else{
		p.s_tot = path;
		return;
	}

	while(p.N_w < lim){
		//Move to wall
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		//Reflect?
		if(G.uniform() < rho){
			//Yes
		//	p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			//Can photon propagate distance s0?
			if(s0 <= path){
				p.N_w ++;
				p.s_tot += s0;
				path -= s0;
			}
			else{
				p.s_tot += path;
				break;
			}
		}
		else{
			//No
			break;
		}
	}
}

void traceSingleDiffuseInvariant(Photon &p, double &R, double &rho, size_t &lim, Generator &G){
	double s0 = 0.0;
	//Move from initial position and direction
	s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
	//diffuseLocal(p.u,G);
	//s0 = 2.0*R*p.u.u_z;
	p.s_tot += s0;
	p.N_w ++;

	while(p.N_w < lim){
		if(G.uniform() < rho){
			p.N_w ++;
			diffuseLocal(p.u, G);
			s0 = 2.0*R*p.u.u_z;
			p.s_tot += s0;
		}
		else{
			break;
		}
	}
}

void traceCollectionDiffuse(vector<Photon> &P, double &R, double &rho, size_t &lim, Generator &G){
	for(size_t i=0; i<P.size(); i++){
		traceSingleDiffuse(P[i], R, rho, lim, G);
	}
}



double meanNumberCollision(std::vector<Photon> &P){
	size_t col = 0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		col += P[i].N_w;
	}
	double collisions = (double) col;
	double N = (double) len;
	return (collisions/N);
}

double stddevNumberCollision(std::vector<Photon> &P, double &mean){
	double var = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		double col = (double) P[i].N_w;
		var += (col-mean)*(col-mean);
	}
	double N = (double) len;
	return sqrt(var/(N-1.0));
}

double meanPathLength(std::vector<Photon> &P){
	//using namespace arma; 
	double s = 0.0;
	//double s_test = 0.0;
	size_t len = P.size();
	//vec path(len);
	for(size_t i=0; i<len; i++){
		s += P[i].s_tot;
		//path(i) = P[i].s_tot;
	}
	/*
	path = sort(path);
	for(size_t i=0; i<len; i++){
		s_test += path(i);
	}
	*/
	double N = (double) len;
	return (s/N);
	//return (s_test/N);
}

double meanPathLengthCollected(std::vector<Photon> &photons){
	double s = 0.0;
	size_t N = photons.size();
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			s += photons[i].s_tot;
		}
	}
	double N_d = (double) N;
	return (s/N);
}

double stddevPathLength(std::vector<Photon> &P, double &mean){
	double var = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		var += (P[i].s_tot - mean)*(P[i].s_tot - mean);
	}
	double N = (double) len;
	return sqrt(var/(N-1.0));
}

double stddevPathLengthCollected(std::vector<Photon> &P, double &mean){
	double var = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		if(P[i].isCollected){
			var += (P[i].s_tot - mean)*(P[i].s_tot - mean);
		}
	}
	double N = (double) len;
	return sqrt(var/(N-1.0));
}

void histCollisions(vector<Photon> &photons, string fileName){

	using namespace arma;
	size_t N = photons.size();
	uvec collision_counter(N);
	for(size_t i=0; i<N; i++){
		collision_counter(i) = photons[i].N_w;
	}
	size_t min_N_w = collision_counter.min();
	size_t max_N_w = collision_counter.max();
	size_t range = max_N_w - min_N_w + 1;

	vec centers = linspace<vec>((double) min_N_w, (double) max_N_w, range);
	vec collisions = conv_to< vec >::from(collision_counter);
	uvec hist_int = hist(collisions, centers);
	vec hist_double = conv_to< vec >::from(hist_int);
	hist_double = hist_double/((double) N);

	mat collision_hist(range,2);	
	collision_hist.col(0) = centers;
	collision_hist.col(1) = hist_double;
	collision_hist.save(fileName, raw_ascii);
}

void histPathLength(vector<Photon> &photons, string fileName, double min_s, double max_s, size_t bins){
	using namespace arma;
	size_t N = photons.size();
	vec path(N);
	for(size_t i=0; i<N; i++){
		path(i) = photons[i].N_w;
	}
	double min_s0 = path.min();
	double max_s0 = path.max();

	vec centers = linspace<vec>(min_s0, max_s0, bins);
	vec paths = conv_to< vec >::from(path);
	uvec hist_int = hist(paths, centers);
	vec hist_double = conv_to< vec >::from(hist_int);
	hist_double = hist_double/((double) N);

	mat path_hist(bins,2);	
	path_hist.col(0) = centers;
	path_hist.col(1) = hist_double;
	path_hist.save(fileName, raw_ascii);
}

void savePhotons(std::vector<Photon> &photons, std::string fileName, double R, double rho, double N_photons, size_t seed){

	ofstream myFile;
	myFile.open(fileName.c_str());
	myFile << "#Radius: " << R << "\n" << "#Reflectivity: " << rho << "\n" << "#No. photons: "
	       << N_photons << "\n" << "#Seed: " << seed << "\n";
	myFile << "#N_w s_tot\n";
	for(size_t i=0; i<N_photons;i++){
		myFile << photons[i].N_w << " " << std::setprecision (17) << photons[i].s_tot << "\n";
	}
	myFile.close();
}


void initPhotons(std::vector<Photon> &photons, double R, double z_s, double cos_theta0, double a_s, double b_s, Generator &G){

	size_t N = photons.size();
	for(size_t i=0; i<N; i++){
		//random points inside unit circle
		double chi_x = G.uniform();
		double chi_y = G.uniform();
		while((chi_x*chi_x + chi_y*chi_y) > 1.0){
			chi_x = G.uniform();
			chi_y = G.uniform();
		}
		//scaling
		chi_x = a_s*chi_x;
		chi_y = a_s*chi_y;
		//set position
		photons[i].r.z = z_s;
		photons[i].r.x = chi_x - b_s;
		photons[i].r.y = chi_y;
		double chi_1 = G.uniform();
		double chi_2 = G.uniform();
		double z = (1.0 - cos_theta0)*chi_2 + cos_theta0;
		//set direction
		photons[i].u.u_z = -z;
		photons[i].u.u_x = cos(2.0*M_PI*chi_1)*sqrt(1.0 - z*z);
		photons[i].u.u_y = sin(2.0*M_PI*chi_1)*sqrt(1.0 - z*z);
		//cout << "u norm: " << photons[i].u.norm() << endl;
	}
}


void tracePhotonEmptyIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &lim, Generator &G){
	//Move from initial position and direction
	double dot_prod = p.r.x*p.u.u_x + p.r.y*p.u.u_y + p.r.z*p.u.u_z;
//	cout << "dot_prod: " << dot_prod << endl;
	double s0 = -dot_prod + sqrt(dot_prod*dot_prod + R*R - p.r.square());
//	cout << "s0 first: " << s0 << endl;


	p.s_tot += s0;
	p.N_w ++;

	while(p.N_w < lim){
		//Move to wall
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		//cout << "R: " << p.r.norm() << endl;
		//Reflect?
		if(G.uniform() < rho){
			//Yes
		//	p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			//Can photon propagate distance s0?
			if((p.r.z + s0*p.u.u_z) > z_s){
				if(p.u.u_z < cos_theta0){
					p.isAbsorbedEntrance = true;
					break;
				}
				else{
					double s = (z_s - p.r.z)/p.u.u_z;
					//cout << "s: " << s << endl;
					double x = p.r.x + s*p.u.u_x;
					double y = p.r.y + s*p.u.u_y;
					if(((x-b_p)*(x-b_p) + y*y) < (a_p*a_p)){
						p.isCollected = true;
						//cout << "collected!" << endl;
						p.s_tot += s;
						break;
					}
					else{
						p.isAbsorbedEntrance = true;
						break;
					}
				}
			}
			else{
				p.N_w ++;
				p.s_tot += s0;
			}
		}
		else{
			//No
			p.isAbsorbedWall = true;
			break;
		}
	}
}

void getStats(std::vector<Photon> &photons, double &eps_c, double &eps_e, double &eps_w, size_t &N){
	double N_d = (double) N;
	eps_c = 0.0;
	eps_e = 0.0;
	eps_w = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			eps_c += 1.0;
		}else if(photons[i].isAbsorbedEntrance){
			eps_e += 1.0;
		}else{ 
			eps_w += 1.0;
		}
	}
	eps_c = eps_c/N_d;
	eps_e = eps_e/N_d;
	eps_w = eps_w/N_d;
}
		
















