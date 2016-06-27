/*
 * $id$ 
 *
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
*/
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "info.h"
#include "functions_pref.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedefs.h"
#include "xdrf.h"
#include "gmxfio.h"
#include "xtcio.h"
#include "smalloc.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"

//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/typedefs.h"

int main(int argn, char* args[])
{
	static char *desc[] = {
		"This script will read the xtc trajectory from an MD simulation",
		"performed by GROMACS and output all data related to the analysis",
		"of preferential interactions.",
		"Thu Aug 20 16:49:21 EDT 2015."
	};

	unsigned int num_blocks=1;
	unsigned int protein_sequence=20; /* parse tpr.itp for info */
        unsigned int N_co=4; /* look at topology for number of co-solvents */
        unsigned int num_sol=7906; /* look at topology as well */
	unsigned int skip_nr=1;
	unsigned int granularity=2*100;
	unsigned int f0=0;

	t_pargs pa[] = {
		{"-f0", FALSE, etINT, {&f0}, "Read after frame f0."},
		{"-skip", FALSE, etINT, {&skip_nr}, "Only write every nr-th frame."},
		{"-block", FALSE, etINT, {&num_blocks}, "Number of blocks for the trajectory to be divided into."},
		{"-seq", FALSE, etINT,{&protein_sequence}, "Number of amino acids in the protein sequence."},
		{"-Nco", FALSE, etINT,{&N_co}, "Number of cosolvent molecules."},
		{"-Sol", FALSE, etINT,{&num_sol}, "Number of water molecules."}
		
	};	

	char 		title[STRLEN];	
	t_topology 	top;
	int		ePBC;
	rvec		*xtop;
	matrix		box;
	
	t_filenm fnm[] = {
		{ efTPS, NULL, NULL, ffREAD },	 /*  topology 	*/
		{ efTRX, "-f", NULL, ffREAD } 	 /*  trajcetory  */
	};

#define NFILE asize(fnm)

	/* Interface. Adds default options. */	
	CopyRight(stderr,args[0]);
	parse_common_args(&argn,args,PCA_CAN_TIME | PCA_CAN_VIEW, NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
	
	read_tps_conf( ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &xtop, NULL, box, TRUE);
	sfree(xtop);
/*******************************************************************************************************************/
/*	int iter;
	for (iter=0;iter<top.atoms.nr/100;++iter){
	printf("Atom name: %s\tAtom charge: %f\n", *(top.atoms.atomname[iter]), top.atoms.atom[iter].q );
	printf("Atom type: %d\tAtomic Residue NUmber: %d\n", top.atoms.atom[iter].type, top.atoms.atom[iter].resnr );
	printf("Chain Identifier: %u\tNr of residues names: %d\n", top.atoms.atom[iter].ptype, top.atoms.nres );
	printf("Residue name: %s\n\n", *(top.atoms.resname[iter]) );
	}
*/
/*******************************************************************************************************************/
int status;
t_trxframe fr;
int flags = TRX_READ_X;

	/* Frame Number */
        unsigned int bframes=0, frames=0;	
	int bwrite;
	read_first_frame(&status, ftp2fn(efTRX,NFILE,fnm), &fr, flags);
        do {
		bwrite = bframes % skip_nr; 
		++bframes; 
		if (bwrite==0){ ++frames; }
	} while ( read_next_frame(status,&fr) );
	
	std::vector<info> protein;
        std::vector<info> cosolvent;
	std::vector<info> solvent;

	std::vector<double> cutoff, coefficients;
	double molals=0, molars=0;

/* BLOCK AVERAGE */
typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;

	unsigned int partition = (frames)/num_blocks;
	unsigned int frame_counter=0;
	std::vector<double> positions(granularity), averages(granularity), averages_sqr(granularity), sigma(granularity);
	Matrix values(num_blocks,Row(granularity)), values_sqr(num_blocks,Row(granularity));

	unsigned int frame=0, twrite=0;
read_first_frame(&status, ftp2fn(efTRX,NFILE,fnm), &fr, flags);
do 
{
	twrite = frame % skip_nr;	/* if twrite==0 then read frame. */
	++frame;
	if ((frame > f0) && (twrite==0))
	{
	/* READING */
		unsigned int protein_co=protein_sequence+N_co;
		unsigned int solvated=protein_co+num_sol;	
		for (unsigned int n=0; n<top.atoms.nr; ++n)
		{
			/* Parsing */
			info atomo;
			atomo.atom = *(top.atoms.atomname[n]); 
			atomo.resnum = top.atoms.atom[n].resnr;
			atomo.resname = *(top.atoms.resname[atomo.resnum]); 
			atomo.x = fr.x[n][XX]; atomo.y = fr.x[n][YY]; atomo.z = fr.x[n][ZZ]; 
					
			/* Appending */
			if ( atomo.atom[0]!='H')
			{
				if ( atomo.resnum < protein_sequence )
				{ 
					protein.push_back(atomo);
				}
			
				if ( (atomo.resnum >= protein_sequence) && (atomo.resnum < protein_co) )
				{ 
					cosolvent.push_back(atomo); 
				}

				if ( (atomo.resnum >= protein_co) && (atomo.resnum < solvated) )
				{
					solvent.push_back(atomo);
				}
			}					
		}

	/* MAIN */      
		double box_size = fr.box[0][0];

	        /* Local calculation for cosolvents */
                std::vector<double> close_co;                                     /* Stores Distances of closest */
		for (unsigned int i=0; i<N_co; ++i)
                {
                        std::vector<double> local_co;
                        double minimum_calculated = box_size*box_size;
			unsigned int Nco_size=cosolvent.size()/N_co;
                        for (unsigned int n=(i*Nco_size); n<((i+1)*Nco_size); ++n)
                        {
                                for (unsigned int j=0; j<protein.size(); ++j)
                                {
                                        double dist = cosolvent[n].distance(protein[j],box_size);
                                        if (dist < minimum_calculated ){ minimum_calculated = dist; }
                                }
                                local_co.push_back(minimum_calculated);
                        }
                        close_co.push_back( *min_element(local_co.begin(), local_co.end()) );
                }	
		/* Local calculation for water */
                std::vector<double> local_wat;
                for (unsigned int i=0; i<solvent.size(); ++i)
                {
                        double minimum_calculated = box_size*box_size;
                        for (unsigned int j=0; j<protein.size(); ++j)
                        {
                                double dist = solvent[i].distance(protein[j],box_size);
                                if (dist < minimum_calculated){ minimum_calculated = dist; }
                        }
                        local_wat.push_back(minimum_calculated);
                } 

		/* Preferential Interaction Coefficient */
                double N_wat = solvent.size();
		linspace(0,box_size,granularity,cutoff);
                for (unsigned int x=0; x<granularity; ++x)
                {
                        unsigned int n3=0, n1=0;
                        for (unsigned int i=0; i<close_co.size(); ++i){ if (close_co[i] < (cutoff[x]*cutoff[x])) ++n3; }
                        for (unsigned int i=0; i<local_wat.size(); ++i){ if (local_wat[i] < (cutoff[x]*cutoff[x])) ++n1; }
                        coefficients.push_back(pref_coef(N_co, N_wat, n3, n1));
                }
		molals+= molality(N_co, (int) N_wat);
                molars+= molarity(N_co, box_size);

	        /* Storage */
	        for (unsigned int i=0; i<granularity; ++i){ positions[i] += cutoff[i]; }
                for (unsigned int j=0; j<granularity; ++j)
                {
                        values[(int) (frame_counter/partition)][j] += coefficients[j];
                        values_sqr[(int) (frame_counter/partition)][j] += coefficients[j]*coefficients[j];
                }


		cutoff.clear();
        	coefficients.clear();
		protein.clear();
		cosolvent.clear();
		solvent.clear();
	
		++frame_counter;
	}
} while ( read_next_frame(status,&fr) );
	thanx(stderr);

/* FINAL REUSLTS */
        double normalization = partition;
        molals /= (num_blocks*normalization);
        molars /= (num_blocks*normalization);
        for (unsigned int i=0; i<granularity; ++i){ positions[i] /= (num_blocks*normalization); }
        for (unsigned int i=0; i<num_blocks; ++i)
	{
        	for (unsigned int j=0; j<granularity; ++j)
                {
                	values[i][j] /= normalization;
                        values_sqr[i][j] /= normalization;
                }
        }

/* AVERAGES */
        for (unsigned int i=0; i<granularity; ++i)
        {
                for (unsigned int j=0; j<num_blocks; ++j)
                {
                        averages[i] += values[j][i];
                        averages_sqr[i] += values_sqr[j][i];
                }
                averages[i] /= (double) num_blocks;
                averages_sqr[i] /= (double) num_blocks;
        }

/* STANDARD DEVIATIONS */
        for (unsigned int i=0; i<granularity; ++i)
        {
                sigma[i] += sqrt( averages_sqr[i] - averages[i]*averages[i] );
        }

/* OUTPUT */
        std::ofstream distances("dimensions.dat");
        for (unsigned int i=0; i<granularity/2; ++i){ distances << std::setw(1) << positions[i] <<'\n'; }
        distances.close();

        std::string output1="preferential_val.dat";
        std::string output2="preferential_sig.dat";
        std::string result;
        std::ostringstream convert;
        convert << num_blocks; //skip;
        result = convert.str();
        output1.insert(16,result);
        output2.insert(16,result);
        std::ofstream pic1(output1.c_str());
        std::ofstream pic2(output2.c_str());
        for (unsigned int i=0; i<granularity/2; ++i)
        {
                pic1 << std::setw(10)<< averages[i] <<'\n';
                pic2 << std::setw(10)<< sigma[i] <<'\n';
        }
        pic1.close();
        pic2.close();

        std::ofstream concentrations("concentrations.dat");
        concentrations<< std::setw(10) << "Molality" <<'\t'<< std::setw(10) << "Molarity" <<'\n';
        concentrations<< std::setw(10) << molals <<'\t'<< std::setw(10) << molars << '\n';
        concentrations.close();

        std::cout << "ANALYSIS COMPLETED!" <<'\n';
	std::cout<<"Total number of frames: "<<frame_counter<<'\n';

return 0;
}


#ifdef __cplusplus
} 
#endif
