#include "Photon.h"
using namespace photon_IS;

Photon::Photon(): r(), u(), W(1.0), N_w(0), N_s(0) {

}


Photon::Photon(const Coord &coord, const Dir &dir): r(coord), u(dir), W(1.0), N_w(0), N_s() {

}


