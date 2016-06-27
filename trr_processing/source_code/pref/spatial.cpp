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
#include "info.h"               /* Class                                                */
#include "pref_co.h"            /* preferential interaction coefficient subroutines     */
#include "file_pro.h"           /* Makes protein structure from pdb                     */

#include "spatial.h"
#ifdef _OPENMP
#include <omp.h>
#endif

double distance(info* in, info* out, double L)
{
        double x1 = in->x, y1 = in->y, z1 = in->z;
        double x2 = out->x, y2 = out->y, z2 = out->z;
        double delx = x2 - x1 - std::floor((x2-x1)/L + 0.5)*L;
        double dely = y2 - y1 - std::floor((y2-y1)/L + 0.5)*L;
        double delz = z2 - z1 - std::floor((z2-z1)/L + 0.5)*L;
        return std::sqrt( delx*delx + dely*dely + delz*delz );
}


void linspace(double start, double end, unsigned int intervals, std::vector<double> &vec)
{
        double delta = (end - start)/(intervals - 1);
        for (unsigned int i=0; i<intervals; ++i){ vec.push_back((start + delta*i)); }
        return;
}

