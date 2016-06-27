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

int water_count(std::vector<info>& atoms){
	int counter=0;
	for (unsigned int n=0; n<atoms.size(); ++n){
		if (atoms[n].type=="OW"){ counter++; }
	}	
	return counter;
}

double molality(std::vector<info>& atoms, double moles_solute){	// m := moles of solute over the mass in kg of solvent.
	double mm_water = 0.018015;	// kg.
	double m = moles_solute/(water_count(atoms) * mm_water);
	return m;							/* PASS THE NUMBER OF MOLES OF SOLUTE BY ARGV*/
}

void make_protein(std::vector<info>& atoms, std::vector<info*>& protein, unsigned int m, unsigned int n){
	for (unsigned int i=m;i<n;++i){
		if (atoms[i].type!="HW1" && atoms[i].type!="HW2" && atoms[i].type!="CL"){	// strings are already pointers!
			protein.push_back(&atoms[i]);}
	}
}



double distance(info* in, info* out, double L){
	double x1 = in->x, y1 = in->y, z1 = in->z;
	double x2 = out->x, y2 = out->y, z2 = out->z;
	double delx = x2 - x1 - std::floor((x2-x1)/L + 0.5)*L;
	double dely = y2 - y1 - std::floor((y2-y1)/L + 0.5)*L;
	double delz = z2 - z1 - std::floor((z2-z1)/L + 0.5)*L;
        return std::sqrt( delx*delx + dely*dely + delz*delz );
}


void linspace(double start, double end, unsigned int intervals, std::vector<double> &vec){
	double delta = (end - start)/(intervals - 1);
	for (unsigned int i=0; i<intervals; ++i){ vec.push_back(start + delta*i); }
	return;
}


double pref_coef(std::vector<info> atoms, double N_co, double n3, double n1){
	double nwat = water_count(atoms);
	double c = (n3 -n1*(N_co/nwat));	  //(N_co - nwat*(n3/n1));
	return c;
} 


int main(int argc, char* argv[]){	
	unsigned int protein_size=304;	// CHANGE FOR PROTEINS OTHER THAN TRP CAGE
	unsigned int arg_size=27;	// CHANGE FOR COSOLVENTS OTHER THAN ARGININE
	double dimensions=4.70;
	

	std::vector<info> atoms;
	
	if (argc<2){ std::cerr<<"\nError: "<< argv[0]<<" needs a gro frile passed as argument to work.\nAlso the number of cosolvents present in the simulation must be specified.\n"; return 1; }
	
	/* All gro files in current directory */
	std::vector<std::string> directory;
	std::string file;
	
	std::system("ls *gro > currentdir.txt");
	std::ifstream dir("currentdir.txt");
	while (true){
		dir >> file;
		if (dir.fail()){ break; }
		directory.push_back(file);
		if (dir.eof()){ break; }
	} 
	dir.close();
	
	std::vector<double> positions(100), values(100);	// The number of entries is arbitrary, just change this and the number for vector<double> cutoff
        std::vector<double> cutoff;
	//linspace(0,dimensions,100,cutoff);
        std::vector<double> coefficients;
	for (unsigned int q=0; q<directory.size(); ++q){
		std::ifstream gro_file(directory[q].c_str());
		unsigned int N_co= std::atoi(argv[1]);	
		std::cout<< directory[q]<<'\n';		

		/* Body */
        	std::vector<info*> protein;
        	std::vector<info*> cosolvent;
        	std::vector<info*> solvent;

		// Read from file
		if (gro_file.is_open()){
			gro_file.ignore(1000,'\n');	// Ignores the header file until newline.
			
			int num_atoms;
			gro_file >> num_atoms;				
			
			for (int n=0; n<num_atoms; ++n){
				info list(gro_file);	//list.read(gro_file) same as, but nicer than  gro_file >> list.label >> list.type >> list.num >> list.x >> list.y >> list.z; 	
				atoms.push_back(list);
			}
			gro_file >> dimensions;       		
 
			if (gro_file.fail()){ std::cout<<"\nError: nothing on file.\n"; return 1; }
		/* Testing */
		//std::cout<<'\n'<<"Here it comes\t"<<dimensions<<'\n';
		//std::cout<<'\n'<<"number of atoms in file "<<num_atoms<<'\n';
		//std::cout<<'\n'<<"vector size =  "<<atoms.size()<<'\n';
		//std::cout<<'\n'<<"Number of waters in gro file = "<<water_count(atoms)<<'\n';			
	
			make_protein(atoms, protein, 0, protein_size); 	
			make_protein(atoms, cosolvent, protein_size, protein_size+(N_co*arg_size));	
			make_protein(atoms, solvent, protein_size+(N_co*arg_size), num_atoms);
		// test make_protein by writing to file	
		/*for (unsigned int n=0; n<protein.size(); ++n){ protein[n]->print(out); }		// for pointers (dereference then access)*/
	
			gro_file.close();			
		} else {
               		std::cerr << "\nError: Unable to open file.\n"<<'\n';
               		return 1;	
		}

		/* MAIN */
		//std::vector<double> cutoff;
		linspace(0,dimensions,100,cutoff);
		//std::vector<double> coefficients;
		for (unsigned int x=0; x<cutoff.size(); ++x){
			unsigned int n3=0, n1=0;
        	
			/* Local calculation for cosolvents */
        		for (unsigned int i=0; i<N_co; ++i){
				std::vector<double> local_co;
                		double minimum_calculated = dimensions*dimensions;
                		for (unsigned int n=(i*arg_size); n<((i+1)*arg_size); ++n){
                        		for (unsigned int j=0; j<protein.size(); ++j){
                                		double dist = cosolvent[n]->distance(protein[j],dimensions);
                                		if (dist < minimum_calculated ){ minimum_calculated = dist; }
                        		}
                        		local_co.push_back(minimum_calculated);
                		}
                		if (*min_element(local_co.begin(), local_co.end()) < (cutoff[x]*cutoff[x])){ ++n3; }
        		}

        		/* Local calculation for water */
			std::vector<double> local_wat;
        		for (unsigned int i=0; i<solvent.size(); ++i){
                		double minimum_calculated = dimensions*dimensions;
                		for (unsigned int j=0; j<protein.size(); ++j){
                        		double dist = solvent[i]->distance(protein[j],dimensions);
                        		if (dist < minimum_calculated){ minimum_calculated = dist; }
				}
				local_wat.push_back(minimum_calculated);
        		}
			for (unsigned int i=0; i<local_wat.size(); ++i){ if (local_wat[i] < (cutoff[x]*cutoff[x])){ ++n1;} }
			coefficients.push_back(pref_coef(atoms, N_co, n3, n1));
		}
		
		for (unsigned int i=0; i<cutoff.size(); ++i){
			values[i] += coefficients[i];
		}
		for (unsigned int i=0; i<cutoff.size(); ++i){ positions[i] += cutoff[i]; }
		coefficients.clear();
		cutoff.clear();
		atoms.clear();
		//for (unsigned int i=0; i<cutoff.size(); ++i){ pic << std::setw(10)<< (cutoff[i]*cutoff[i]) <<'\t'<< std::setw(10)<< coefficients[i] <<'\n'; }
	}
	for (unsigned int i=0; i<values.size(); ++i){ values[i] /= directory.size(); }
        for (unsigned int i=0; i<values.size(); ++i){ positions[i] /= directory.size(); }

	std::ofstream pic("pic_calculation.dat");	
	for (unsigned int i=0; i<values.size(); ++i){ pic << std::setw(10)<< positions[i] <<'\t'<< std::setw(10)<< values[i] <<'\n'; }
	pic.close();


	return 0;
}

