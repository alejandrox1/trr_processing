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
#include "info.h"       /* Class */
#include "pref_co.h"    /* water_count(), pref_coef() */
#include "file_pro.h"   /* make_protein() */
#include "spatial.h"    /* distance(), linspace() */

#include "concentration.h"

double molality_(std::vector<info>& atoms, double moles_solute)
{ // m := moles of solute over the mass in kg of solvent.
        double mm_water = 0.018015;     // kg.
        double m = moles_solute/(water_count(atoms) * mm_water);
        return m;                                                       /* PASS THE NUMBER OF MOLES OF SOLUTE BY ARGV*/
}


/* 	n_solute == Nmg		*/
double molality(int n_solute, int n_wat)
{
	double w_mm = 0.018015;		/* molar mass of water in Kg. */
	double m = n_solute/(n_wat * w_mm);
	return m;
}

double molarity(int n_solute, double L)
{
	double Na=6.023e23;		/* Avogrado's number 			*/
	double V=(L*L*L)*(1.0e-24);	/* Volume in nm/liter = 1.0e-3 m^3 	*/
	double M= n_solute/(Na*V);
	return M;
}
