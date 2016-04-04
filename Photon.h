#ifndef PHOTON_H
#define PHOTON_H

#include "Coord.h"
#include "Dir.h"
#include <cstring> //for size_t

namespace photon_IS{

class Photon{
	public:
		Photon();
		Photon(const Coord &coord, const Dir &dir);

		Coord r;
		Dir u;
		double W;
		double s_tot; //Total distance travelled
		size_t N_w;
		size_t N_s;


	private:

};

} //End namespace photon_IS

#endif
