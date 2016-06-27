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
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/xdrf.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/gmxfio.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/xtcio.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/smalloc.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/vec.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/futil.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/copyrite.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/statutil.h"
//#include "/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/tpxio.h"


int main(int argn, char* args[])
{
	static char *desc[] = {
		"This is a small test programn for using all GROMACS capabilities.",
		"May the force be with you in this long and winded journey!"
	};

/********************************************************************************************************/

	static int n=1;

	t_pargs pa[] = {
		{"-n", FALSE, etINT, {&n}, "Plot data for number n (starting on 1)"}
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

	int iter;
	for (iter=0;iter<top.atoms.nr;++iter){
	printf("Atom name: %s\tAtom charge: %f\n", *(top.atoms.atomname[iter]), top.atoms.atom[iter].q );
	printf("Atom type: %d\tAtomic Residue NUmber: %d\n", top.atoms.atom[iter].type, top.atoms.atom[iter].resnr );
	//printf("Chain Identifier: %u\tNr of residues names: %d\n", top.atoms.atom[iter].chain, top.atoms.nres );
	//printf("Residue name: %s\n\n", *(top.atoms.resname[iter]) );
	}


/*******************************************************************************************************************/
	int status;
	t_trxframe fr;
	int flags = TRX_READ_X;
	
	int nn=0;

	read_first_frame(&status, ftp2fn(efTRX,NFILE,fnm), &fr, flags);
	do {
		printf("Coordinates at step=%d; t=%f : %f %f %f\n", fr.step,fr.time,fr.x[nn][XX],fr.x[nn][YY],fr.x[nn][ZZ] );
	} while ( read_next_frame(status,&fr) );
	thanx(stderr);


return 0;
}


#ifdef __cplusplus
} 
#endif
