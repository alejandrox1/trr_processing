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
#include "info.h"               /* Class */
#include "pref_co.h"            /* water_count(), pref_coef() */
#include "file_pro.h"           /* make_protein() */
#include "spatial.h"            /* distance(), linspace() */
#include "concentration.h"      /* molality() */
#ifdef _OPENMP
#include <omp.h>
#endif

void make_protein(std::vector<info>& atoms, std::vector<info*>& protein, unsigned int m, unsigned int n)
{
        for (unsigned int i=m;i<n;++i)
	{
		if (atoms[i].type.find("HW1")==std::string::npos && atoms[i].type.find("HW2")==std::string::npos && atoms[i].type.find("CL")==std::string::npos)
		//  if (atoms[i].type!="HW1" && atoms[i].type!="HW2" && atoms[i].type!="CL")
		{       // strings are already pointers!
                        protein.push_back(&atoms[i]);
		}
        }
}
void mod_make_protein(std::vector<mod_info>& atoms, std::vector<mod_info*>& protein, unsigned int m, unsigned int n)
{
        for (unsigned int i=m;i<n;++i)
        {
                if (atoms[i].type.find("HW1")==std::string::npos && atoms[i].type.find("HW2")==std::string::npos && atoms[i].type.find("CL")==std::string::npos)
                //  if (atoms[i].type!="HW1" && atoms[i].type!="HW2" && atoms[i].type!="CL")
                {       // strings are already pointers!
                        protein.push_back(&atoms[i]);
                }
        }
}

void ls_lh_gro(std::vector<std::string> &directory, unsigned int skip)
{
	std::string file, ignore;

        std::system("ls *gro > currentdir.txt");
        std::ifstream dir("currentdir.txt");
        unsigned int counter=0;
	dir >> file; 
	directory.push_back(file);
	while (true){
		if (skip==counter)
		{
              		dir >> file;
			counter=0;

			if (dir.fail()){ break; }
        	        directory.push_back(file);
	                if (dir.eof()){ break; }

		} else {
			++counter;
			dir >> ignore;
		}
        }
        dir.close();
}

void parser(std::vector<info> &atoms, std::vector<mod_info> &mod_atoms, std::string file_handle, double &dimensions, int &num_atoms)
{ 
	std::ifstream gro_file(file_handle.c_str());

	if (gro_file.is_open())
        {
       		gro_file.ignore(1000,'\n');     // Ignores the header file until newline.

                gro_file >> num_atoms;
                for (int n=0; n<num_atoms; ++n)
                {	//list.read(gro_file) same as, but nicer than  gro_file >> list.label >> list.type >> list.num >> list.x >> list.y >> list.z;
                	if (n<9999)
                        {
                        	info list(gro_file);      
                        	atoms.push_back(list);
                        	//atoms[n]=list;
                        } else {
                        	mod_info mod_list;
                                mod_list.mod_read(gro_file);
                                //atoms[n]=list;
                                mod_atoms.push_back(mod_list);
                	}
		}
                gro_file >> dimensions;
                if (gro_file.fail()){ std::cout<<"\nError: nothing on file.\n"; abort(); }
                gro_file.close();
	} else {
        	std::cerr << "\nError: Unable to open file.\n"<<'\n';
                abort();
	}
}
