#include "Coord.h"
#include "Dir.h"
#include "Photon.h"
#include <iostream>

using namespace photon_IS;


int main(){
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

	std::cout << p2.u.u_z;






	
















	return 0;
}
