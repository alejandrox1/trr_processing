#include <iostream>
#include <vector>
#include <cmath>
#include "functions_pref.h"


void linspace(double start, double end, unsigned int intervals, std::vector<double> &vec)
{
        double delta = (end - start)/(intervals - 1);
        for (unsigned int i=0; i<intervals; ++i){ vec.push_back((start + delta*i)); }
        return;
}

double molality(int N_solute, int N_wat)
{
        double w_mm = 0.018015;         /* molar mass of water in Kg. */
        double m = N_solute/(N_wat * w_mm);
        return m;
}

double molarity(int N_solute, double box_size)
{
        double Na=6.023e23;             /* Avogrado's number                    */
        double V=(box_size*box_size*box_size)*(1.0e-24);        /* Volume in nm/liter = 1.0e-3 m^3      */
        double M= N_solute/(Na*V);
        return M;
}

double pref_coef(double N_co, double N_wat, double n3, double n1)
{
        double c = (n3 -n1*(N_co/N_wat));          //(N_co - nwat*(n3/n1));
        return c;
}


