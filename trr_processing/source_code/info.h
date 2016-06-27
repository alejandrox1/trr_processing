#ifndef INFO_H
#define INFO_H

class Info{
	public:
	
		Info(std::string atom="", std::string resname="", unsigned int resnum=0, double x=0, double y=0, double z=0);

		double distance(Info i,double box_size);
		template <typename T>
		void print(T & stream) const; 
		void read(std::istream& stream);
		
		/* MODIFIERS */
		void setAtom(std::string atm );
		void setResname(std::string residue );
		void setResnum(unsigned int res ); 
		void setCoords(double xx, double yy, double zz); 
		
		/* ACCESORIES */
		std::string getAtom() const;
		std::string getResname() const;
		unsigned int getResnum() const;
		std::vector<double> getCoords() const; 
	
	private:
		std::string atom; 
		std::string resname;
        	unsigned int resnum;
        	double x, y, z;

};
/********************************************************************************************************************************************************/

Info::Info(std::string a, std::string r, unsigned int rn, double ax, double ay, double az) : atom(a), resname(r), resnum(rn), x(ax), y(ay), z(az)
{}


template <typename T>
void Info::print(T & stream) const{
        stream << atom <<'\t'<< resname <<'\t'<< resnum <<'\t'<< x<<'\t'<< y <<'\t'<< z << '\n';
}


double Info::distance(Info i,double box_size){
        double delx = i.x - x - std::floor((i.x-x)/box_size + 0.5)*box_size;
        double dely = i.y - y - std::floor((i.y-y)/box_size + 0.5)*box_size;
        double delz = i.z - z - std::floor((i.z-z)/box_size + 0.5)*box_size;
        return delx*delx + dely*dely + delz*delz;
}


void Info::read(std::istream& stream){
        stream >> atom >> resname >> resnum >> x >> y >> z;
}


void Info::setAtom(std::string atm ){ 
	atom = atm; 
}


void Info::setResname(std::string residue ){ 
	resname = residue; 
}


void Info::setResnum(unsigned int res ){ 
	resnum = res; 
}


void Info::setCoords(double xx, double yy, double zz){ 
	x=xx; 
	y=yy; 
	z=zz; 
}


std::string Info::getAtom() const { 
	return atom; 
}


std::string Info::getResname() const { 
	return resname; 
}


unsigned int Info::getResnum() const { 
	return resnum; 
}


std::vector<double> Info::getCoords() const{
	std::vector<double> vec;
        vec.push_back(x); vec.push_back(y); vec.push_back(z);
        return vec;
}


#endif

