#ifndef info_H
#define info_H


class info{
public:
        std::string label, type;
        unsigned int num;
        double x, y, z, velx, vely, velz;
	void print(std::ostream& stream); 
        void read(std::istream& stream);
	double distance(info *i,double L); 
        info(std::istream& stream); 
        info(){};
	info(const info& Ref);
	info& operator= (const info& Ref);	
};

class mod_info{
public:
        std::string label, type;
        double x, y, z, velx, vely, velz;
        void mod_print(std::ostream& stream);
        void mod_read(std::istream& stream);
        double mod_distance(info *i,double L);
//        mod_info(std::istream& stream);
        mod_info(){};
        mod_info(const mod_info& Ref);
        mod_info& operator= (const mod_info& Ref);
};


#endif

/*
class info{
public:
        std::string label, type;
        unsigned int num;
        double x, y, z, velx, vely, velz;
        void print(std::ostream& stream){ stream << label <<'\t'<< type <<'\t'<< num <<'\t'<< x<<'\t'<<y<<'\t'<< z << '\n'; }
        void read(std::istream& stream){ stream >> label >> type >> num >> x >> y >> z; }
        double distance(info *i,double L) {
                double delx = i->x - x - std::floor((i->x-x)/L + 0.5)*L;
                double dely = i->y - y - std::floor((i->y-y)/L + 0.5)*L;
                double delz = i->z - z - std::floor((i->z-z)/L + 0.5)*L;
                return delx*delx + dely*dely + delz*delz;
        }
        info(std::istream& stream){ stream >> label >> type >> num >> x >> y >> z; }
        info(){}
};
*/
