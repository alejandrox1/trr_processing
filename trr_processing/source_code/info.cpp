#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "info.h"


info::info(std::string a, std::string r, unsigned int rn, double ax, double ay, double az) : atom(a), resname(r), resnum(rn), x(ax), y(ay), z(az)
{}

template <typename T>
void info::print(T & stream) const
{ 
	stream << atom <<'\t'<< resname <<'\t'<< resnum <<'\t'<< x<<'\t'<< y <<'\t'<< z << '\n'; 
}

double info::distance(info i,double box_size) 
{
        double delx = i.x - x - std::floor((i.x-x)/box_size + 0.5)*box_size;
        double dely = i.y - y - std::floor((i.y-y)/box_size + 0.5)*box_size;
        double delz = i.z - z - std::floor((i.z-z)/box_size + 0.5)*box_size;
	return delx*delx + dely*dely + delz*delz;
}

void info::read(std::istream& stream)
{
        stream >> atom >> resname >> resnum >> x >> y >> z;
}


void info::readAtom(std::string a)
{
	atom = a;
}
