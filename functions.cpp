#include "functions.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

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
	double s = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		s += P[i].s_tot;
	}
	double N = (double) len;
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
