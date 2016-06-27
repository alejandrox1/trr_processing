#ifndef pref_co_H
#define pref_co_H

int water_count(std::vector<info>& atoms);
int mod_water_count(std::vector<mod_info>& atoms);

double pref_coef(double N_co, double N_wat, double n3, double n1);

#endif
