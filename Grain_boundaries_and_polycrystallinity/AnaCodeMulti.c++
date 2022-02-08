#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<fftw3.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_coupling.h>

#define MAX(a,b)	((a)>(b))?(a):(b)
#define MIN(a,b)	((a)<(b))?(a):(b)
#define cube(x)	((x)*(x)*(x))
#define sqr(x)		((x)*(x))
#define FactPrec 8
#define MaxKappa 30000

#include"functions/Array.h"
#include"Structures.h"
#include"reader/readLammps.h"

#include"functions/Gr.h"
#include"functions/Msd.h"
#include"functions/Isf.h"
#include"functions/OrPar.h"
#include"functions/OrParLocAvg.h"
#include"functions/Sq.h"
#include"functions/Coher.h"
#include"functions/G5r.h"
#include"functions/ClSize.h"
#include"functions/YgCoordNo.h"
#include"functions/Interface_v5.h"
#include<time.h>
main (int argc, char *argv[])
{
	struct SIMDAT s;
	struct DATASIM d;
	struct COORDNO cn; // check if this works
	char buffer[500];
	FILE *fp;
	int t;
	
	sprintf(s.output,"%s", argv[3]); // fourth name in command line is output file
	sprintf(s.buffer,"%s", argv[1]); // second is input file
	printf("%s - arguments come from here \n",s.buffer);
	printf("%s - output file \n",s.output);
	fp=fopen(s.buffer,"r");		
	if(fp==NULL)
	{
		printf("null pointer problems \n");
	}
	readPara(fp,&s,&d);
	printf("%d - test output\n",s.StartAt);
	fclose(fp);

/////////////////////////// Initialization ///////////////////////////////////////////

	//t=s.StartAt*s.Step;
	t=s.StartAt;
	sprintf(s.buffer,"%s%s%d.dat", argv[2],s.name,t);
	printf("%s - input file name %d time\n",s.buffer,t);
	fp=fopen(s.buffer,"r");	
	if(fp==NULL)
	{
		printf("null pointer problems\n");
	}	
	readLammps(fp,&s,&d,t);
	fclose(fp);
	printf("have we read coords at %d\n",t);
	if(s.Gr==1)
	Init_Gr(&s,&d,&densityK);
	if(s.Msd==1)
	Init_Msd(&s,&d,&MeanSquareDispl);
	if(s.Isf==1)
	Init_Isf(&s,&d,&IsfK);
	if(s.coordno_flag==1) // remember coordNo needs to be 1 for Orpar and clsize calcs
	{
	Init_CoordNo(&s,&d,&CoordNo);
	printf("calculating coord nos %d\n",s.coordno_flag);
	}

	if(s.OrPar==1)
	Init_OrPar(&s,&d,&order);
	if(s.Sq==1)
	Init_Sq(&s,&d,&esseq);
	if(s.Coher==1)
	Init_Coher(&s,&d,&coher);
	if(s.G5r==1)
	Init_G5r(&s,&d,&gCinqueR);
	if(s.ClSize==1)
	Init_ClSize(&s,&d,&clsize);
	if(s.Interface==1)
	Init_Interface(&s,&d,&grain);
    if(s.OrPar_la==1)
	Init_OrParla(&s,&d,&orderla);
	
    printf("done initialising\n");
/////////////////////////// Computation ///////////////////////////////////////////

	//for(t=s.StartAt*s.Step;t<s.StopAt*s.Step;t+=s.Step)
	for(t=s.StartAt;t<s.StopAt;t+=s.Step)
	{
		printf("Time: %d\n",(int)t/s.Step);
		sprintf(s.buffer,"%s%s%d.dat", argv[2],s.name,t);
		printf("%s - input file name - time %d\n",s.buffer,t);
		fp=fopen(s.buffer,"r");	
		if(fp==NULL)
		{
			printf("null pointer problems\n");
		}		
		readLammps(fp,&s,&d,t);
		fclose(fp);
		printf("have we read coords at %d\n",t);
		if(s.Gr==1)
		Compute_Gr(&s,&d,&densityK);
		if(s.Msd==1)
		{
			Compute_Msd(&s,&d,&MeanSquareDispl);
			printf("n files %d\n",MeanSquareDispl.counter);
			if(MeanSquareDispl.counter==s.Msd_limit) 
			{
				printf("printing msd msd limit flag**\n");
				printf("%s\n",s.output);
				printf("%d\n",MeanSquareDispl.fileCounter);
				//printf("%f %f\n",MeanSquareDispl.msd[0][0],MeanSquareDispl.msd[MeanSquareDispl.counter][0]);
				sprintf(buffer,"%sMsd_%d.dat",s.output,MeanSquareDispl.fileCounter);
				printf("the MSD file %s\n",buffer);
				print_Msd(&s,&d,&MeanSquareDispl);
				printf("printed msd\n");
			}
		}
		if(s.Isf==1)
		Compute_Isf(&s,&d,&IsfK);
		if(s.coordno_flag==1)
		{
		printf("computing coordno\n");
		Compute_CoordNo(&s,&d,&CoordNo);		
		}

		if(s.OrPar==1)
		Compute_OrPar(&s,&d,&order);
		if(s.Sq==1)
		Compute_Sq(&s,&d,&esseq);
		if(s.Coher==1)
		Compute_Coher(&s,&d,&coher);
		if(s.G5r==1)
		Compute_G5r(&s,&d,&gCinqueR);
		if(s.ClSize==1)
		{
		Compute_ClSize(&s,&d,&clsize);
		Compute_Cl_correl(&s,&d,&clsize);
		}
		if(s.Interface==1)
		compute_Interface(&s,&d,&grain,t);
        if(s.OrPar_la==1)
		Compute_OrParla(&s,&d,&orderla);
		// we need a separate if for the subroutine to calculate correlations.
	}
	
}


