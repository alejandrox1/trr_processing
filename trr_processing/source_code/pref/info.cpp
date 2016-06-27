#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <string>
#include "info.h"

void info::print(std::ostream& stream)
{ 
	stream << label <<'\t'<< type <<'\t'<< num <<'\t'<< x<<'\t'<<y<<'\t'<< z << '\n'; 
}
void info::read(std::istream& stream)
{ 
	stream >> label >> type >> num >> x >> y >> z; 
}
double info::distance(info *i,double L) 
{
        double delx = i->x - x - std::floor((i->x-x)/L + 0.5)*L;
        double dely = i->y - y - std::floor((i->y-y)/L + 0.5)*L;
        double delz = i->z - z - std::floor((i->z-z)/L + 0.5)*L;
	return delx*delx + dely*dely + delz*delz;
}
info::info(std::istream& stream)
{ 
	stream >> label >> type >> num >> x >> y >> z; 
}
info::info(const info& Ref)
{
    // copy the attributes
    this->label = Ref.label;
    this->type = Ref.type;
    this->num = Ref.num;
    this->x = Ref.x;
    this->y = Ref.y;
    this->z = Ref.z;
    this->velx = Ref.velx;
    this->vely = Ref.vely;
    this->velz = Ref.velz;
}
info& info::operator= (const info& Ref)
{
    // copy the attributes
        this->label = Ref.label;
        this->type = Ref.type;
        this->num = Ref.num;
        this->x = Ref.x;
        this->y = Ref.y;
        this->z = Ref.z;
        this->velx = Ref.velx;
        this->vely = Ref.vely;
        this->velz = Ref.velz;

    // return the existing object
    return *this;
}

/******************************************************************************/

void mod_info::mod_print(std::ostream& stream)
{
        stream << label <<'\t'<< type <<'\t'<< x<<'\t'<<y<<'\t'<< z << '\n';
}
void mod_info::mod_read(std::istream& stream)
{
        stream >> label >> type >> x >> y >> z;
}
double mod_info::mod_distance(info *i,double L)
{
        double delx = i->x - x - std::floor((i->x-x)/L + 0.5)*L;
        double dely = i->y - y - std::floor((i->y-y)/L + 0.5)*L;
        double delz = i->z - z - std::floor((i->z-z)/L + 0.5)*L;
        return delx*delx + dely*dely + delz*delz;
}
mod_info::mod_info(const mod_info& Ref)
{
    // copy the attributes
    this->label = Ref.label;
    this->type = Ref.type;
    this->x = Ref.x;
    this->y = Ref.y;
    this->z = Ref.z;
    this->velx = Ref.velx;
    this->vely = Ref.vely;
    this->velz = Ref.velz;
}
mod_info& mod_info::operator= (const mod_info& Ref)
{
    // copy the attributes
        this->label = Ref.label;
        this->type = Ref.type;
        this->x = Ref.x;
        this->y = Ref.y;
        this->z = Ref.z;
        this->velx = Ref.velx;
        this->vely = Ref.vely;
        this->velz = Ref.velz;

    // return the existing object
    return *this;
}

