#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <unistd.h>             /* getopt */
#include "info.h"               /* Class */
#include "pref_co.h"            /* water_count(), pref_coef() */
#include "file_pro.h"           /* make_protein() */
#include "spatial.h"            /* distance(), linspace() */
#include "concentration.h"      /* molality() */

void test(int &check)
{
	++check;
	std::cout << "CHEKPOINT: " << check <<'\n';
}

int main(int argc, char* argv[]){

	//int check=0;
	
	unsigned int protein_size=0;	
	unsigned int co_size=0;
	unsigned int N_co=0;
	unsigned int granularity=2*100;
	unsigned int skip=0;
	unsigned int block_trans=0;
	unsigned int files_included=0;
/* Input */
	int c;
	while ( (c=getopt(argc,argv,"n: p: c: s: b:"))!=-1 )
	{
		switch (c)
		{
			case 'n':
				N_co=std::atoi(optarg);
				break;
			case 'p':
				protein_size=std::atoi(optarg);
				break;
			case 'c':
				co_size=std::atoi(optarg);
				break;
			case 's':
				skip=std::atoi(optarg);
				break;
			case 'b':
				block_trans=std::atoi(optarg);
				break;
			default:
				std::cerr<<"\nError:\n\t Need topology for number of co-solvent molecules and representative gro files for protein/co-solvent strucuture."<<'\n';
				break; 
		}
	}

/* All gro files in current directory */
	std::vector<std::string> directory;
	ls_lh_gro(directory, skip);
	std::vector<info> atoms;
	std::vector<mod_info> mod_atoms;
	
        std::vector<double> cutoff, coefficients;

/* BLOCK AVERAGE */
	unsigned int block_num=0;
	/*unsigned int block_partition=( directory.size()/ (int) std::pow(2, (double) block_trans) );*/
	unsigned int block_partition= ( directory.size()/ (block_trans+1) );
typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;

	double molals=0, molars=0;
	std::vector<double> positions(granularity), averages(granularity), averages_sqr(granularity), sigma(granularity);
	Matrix values( (block_trans+1), Row(granularity) ), values_sqr( (block_trans+1), Row(granularity) );

while ( block_num <= block_trans ) 
{
	for (unsigned int q= block_num*block_partition; q< (block_num+1)*block_partition; ++q) 
	{ 
		double box_size=0;
		int num_atoms=0;
        	std::vector<info*> protein;
        	std::vector<info*> cosolvent;
        	std::vector<info*> solvent;
		std::vector<mod_info*> solvent_2;
		
/* Read from file */
		parser(atoms, mod_atoms, directory[q], box_size, num_atoms);
	
		make_protein(atoms, protein, 0, protein_size); 
		make_protein(atoms, cosolvent, protein_size, protein_size+(N_co*co_size));	
		if ( num_atoms > 9999)
		{
			make_protein(atoms, solvent, protein_size+(N_co*co_size), 9999);
			mod_make_protein(mod_atoms, solvent_2, 0, num_atoms-9999); 
		} else {
			make_protein(atoms, solvent, protein_size+(N_co*co_size), num_atoms);
		}
/* MAIN */	
	/* Local calculation for cosolvents */
		std::vector<double> close_co;										/* Stores Distances of closest */
        	for (unsigned int i=0; i<N_co; ++i)
		{
			std::vector<double> local_co;
                	double minimum_calculated = box_size*box_size;
                	for (unsigned int n=(i*co_size); n<((i+1)*co_size); ++n)
			{
                       		for (unsigned int j=0; j<protein.size(); ++j)
				{
                               		double dist = cosolvent[n]->distance(protein[j],box_size);
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
                       		double dist = solvent[i]->distance(protein[j],box_size);
                       		if (dist < minimum_calculated){ minimum_calculated = dist; }
			}
			local_wat.push_back(minimum_calculated);
        	}
		if (num_atoms > 9999)
		{
			for (unsigned int i=0; i<solvent_2.size(); ++i)
                	{
                		double minimum_calculated = box_size*box_size;
                        	for (unsigned int j=0; j<protein.size(); ++j)
                        	{
                        		double dist = solvent_2[i]->mod_distance(protein[j],box_size);
                                	if (dist < minimum_calculated){ minimum_calculated = dist; }
                        	}
                        	local_wat.push_back(minimum_calculated);
                	}
		}
	/* Preferential Interaction Coefficient */
        	double N_wat; 
	        if (num_atoms > 9999)
		{ 
			N_wat = solvent.size() + solvent_2.size();
		} else { 
			N_wat = solvent.size(); 
		}
		linspace(0,box_size,granularity,cutoff);
		for (unsigned int x=0; x<cutoff.size(); ++x)
		{
			unsigned int n3=0, n1=0;
			for (unsigned int i=0; i<close_co.size(); ++i){ if (close_co[i] < (cutoff[x]*cutoff[x])){ ++n3;} }
			for (unsigned int i=0; i<local_wat.size(); ++i){ if (local_wat[i] < (cutoff[x]*cutoff[x])){ ++n1;} }
			coefficients.push_back(pref_coef(N_co, N_wat, n3, n1));
		}
		molals+= molality(N_co, (int) N_wat);
		molars+= molarity(N_co, box_size);  
		

	/* Storage */
		for (unsigned int i=0; i<granularity; ++i){ positions[i] += cutoff[i]; }	
		for (unsigned int j=0; j<granularity; ++j)
		{ 
			values[block_num][j] += coefficients[j]; 
			values_sqr[block_num][j] += coefficients[j]*coefficients[j];
		}
	/* Trash Collecting */	
		atoms.clear();
		if (num_atoms > 9999){ mod_atoms.clear(); }
		coefficients.clear();
		cutoff.clear();
	
		++files_included;	// WHASAAP WITH DAT!!!!!!!!!	
	}
	++block_num;
}
/* FINAL REUSLTS */
	double normalization = block_partition;
	std::cout << "norm: " << (block_trans+1)*normalization << ", actual: "<< files_included <<'\n';
	std::cout << "m: "<< molals <<", M: "<<molars <<'\n'; 
	molals /= ((block_trans+1)*normalization);
	molars /= ((block_trans+1)*normalization);
	for (unsigned int i=0; i<granularity; ++i){ positions[i] /= ((block_trans+1)*normalization); }
	for (unsigned int i=0; i<(block_trans+1); ++i)
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
        	for (unsigned int j=0; j<(block_trans+1); ++j)
		{ 
			averages[i] += values[j][i]; 
			averages_sqr[i] += values_sqr[j][i];
		}
                averages[i] /= (double) (block_trans+1);
		averages_sqr[i] /= (double) (block_trans+1);		
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
	convert << block_trans; //skip;
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
	return 0;
}

