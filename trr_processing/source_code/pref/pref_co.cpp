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
#include "info.h"               /* class */

#include "pref_co.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int water_count(std::vector<info>& atoms)
{
        int counter=0;
        for (unsigned int n=0; n<atoms.size(); ++n)
	{
                //if (atoms[n].type=="OW"){ ++counter; }
        	if (atoms[n].type.find("OW") != std::string::npos){ ++counter; }
	}       
        return counter;
}
int mod_water_count(std::vector<mod_info>& atoms)
{
        int counter=0;
        for (unsigned int n=0; n<atoms.size(); ++n)
        {
                //if (atoms[n].type=="OW"){ ++counter; }
                if (atoms[n].type.find("OW") != std::string::npos){ ++counter; }
        }
        return counter;
}


double pref_coef(double N_co, double N_wat, double n3, double n1)
{
        double c = (n3 -n1*(N_co/N_wat));          //(N_co - nwat*(n3/n1));
        return c;
}

