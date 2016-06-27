#ifndef FUNCTIONS_PREF_H
#define FUNCTIONS_PREF_H

void linspace(double start, double end, unsigned int intervals, std::vector<double> &vec);

double molality(int N_solute, int N_wat);

double molarity(int N_solute, double box_size);

double pref_coef(double N_co, double N_wat, double n3, double n1);

#endif
