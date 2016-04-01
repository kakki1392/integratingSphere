#include "functions.h"
#include <cmath>

using namespace photon_IS;

void diffuseReflection(Dir &u, const Generator &G){
	double chi_1 = G.uniform();
	double chi_2 = G.uniform();
	double c = std::cos(2.0*M_PI*chi_1);
	double chi_2_sqrt = std::sqrt(chi_2);
	u.u_x = c*chi_2_sqrt;
	u.u_y = std::sqrt(1.0-c*c)*chi_2_sqrt;
	u.u_z = std::sqrt(1.0-chi_2);
	//CONTINUE
	//
}
	


