#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_coupling.h>
#include "mpi.h"

#define MAX(a,b)	((a)>(b))?(a):(b)
#define MIN(a,b)	((a)<(b))?(a):(b)
#define cube(x)	((x)*(x)*(x))
#define sqr(x)		((x)*(x))

#define SKIN 		0.2
#define minGrow	0.95
#define N 			624
#define minDist2	0.8
#define LINESIZE 1024
#define alat 1.675


#include "functions/Structures.h"
#include "functions/FileInput.h"
#include "functions/Random.h"
#include "functions/NeighList.h"
#include "functions/Potential_Force_LJ_TS.h"
#include "functions/Create_Init.h"
#include "functions/Integrate_Gaussian.h"


int main (int argc, char *argv[])
{
	if(argc!=4)
	{
		printf("Error: Launch %s <Para File> <Source> <Output>\n", argv[0]);
		MPI_Finalize();
		exit(1);
	}

/*************************************************** MPI VARIABLES START **************************************************/

	int rank, size, domrank, received, notreceived;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	MPI_Request request;
	MPI_Group groupall, group;

	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_group( MPI_COMM_WORLD, &groupall);
	
/*************************************************** MPI VARIABLES END **************************************************/

/*************************************************** READING INPUT FILES ************************************************/

	struct SYSTEM structbox;
	struct SYSTEM tmpbox;
	struct APPOGGIO app;
	struct PARAINPUT p;
	
	FILE *fpLog,*fpTrack;
	int i;
	char filename[500],filepart[500],**fileList;
	char string[500];
	double **source;
    double baroLim=0.0;
	int ranspam=1;
	int spamint=1;
	double spamfloat=0.0;
	//printf("Start reading input file %s...\n",argv[1]);
	fpLog=fopen(argv[1],"r");					// Input file (paraUS.dat)
	if(fpLog==NULL) 
	{
		printf("Error: %s not found!\n",argv[1]);
		exit(1);
	}
	readPara(fpLog,&p);
	fclose(fpLog);
	printf("Done input file %s!\n",argv[1]);
    
	source=calloc(size,sizeof(double *));		
	fileList=calloc(size,sizeof(char *));
	for(i=0;i<size;i++)
	{
		source[i]=calloc(3,sizeof(double *));
		fileList[i]=calloc(500,sizeof(char *));
	}
    
    
    
    //printf("N %d ens %d eps %f sig %f rc %f md %f mv %f \n",)
    
	//printf("Start reading source file %s...\n",argv[2]);
	fpLog=fopen(argv[2],"r");		// Source file (source.dat)
	if(fpLog==NULL) 
	{
		printf("Error: %s not found!\n",argv[2]);
		exit(1);
	}
	readSource(fpLog,source,fileList,size);
	fclose(fpLog);
	printf("Done reading source file %s!\n",argv[2]);
    if(p.ens==0)
    {
        structbox.rho=source[rank][0];
        structbox.pressure=source[rank][1];
        structbox.tem=source[rank][2]; 
        
    }
    printf("%f %f %f %f\n",structbox.rho,structbox.pressure,structbox.tem,p.rc);
    if(p.ens==0) // NVT
    {
        sprintf(filepart,"%s/%s_N%d_rho%g_T%g",argv[3],p.label,p.nAto,structbox.rho,structbox.tem);
    }
    if(p.ens==1)
    {
        sprintf(filepart,"%s/%s_N%d_P%g_T%g",argv[3],p.label,p.nAto,structbox.pressure,structbox.tem);
    }
//     if(p.ens==2)
//     {
//         sprintf(filepart,"%s/%s_mu%g_rho%g_T%g",argv[3],p.label,structbox.mu,structbox.rho,structbox.tem);
//     }

	MPI_Barrier(MPI_COMM_WORLD);
	
	inizialized(&p,&structbox,&tmpbox,&app);
	//countList(&p,source,size);		// Conto quanti N0 e T ho in tutto. Mi serve per fare lo swap tra processori
	
	MPI_Barrier(MPI_COMM_WORLD);	
	
/***************************************************************** CREO E/O LEGGO SITUAZIONE INIZIALE **************************************************/

	int mdstart,lcount,mdstep;
    double time=0.0;
	
	if(p.initRead==0)					// initread 0 = Create atoms
	{
		init_genrand(p.seed,&app);
		creaAtom(&p,&structbox,&app,p.startCry); // create liquid
		init_vel(&p,&structbox,&app);
		mdstart=0.0;
	}

	if(p.initRead==1)						// initread 1 = start from STANDARD restart file
	{

		sprintf(filename,"%s-restart1.dat",filepart);
		
		structbox.fprestart=fopen(filename,"r");
		if(structbox.fprestart==NULL) 
		{
			printf("Error: %s not found!\n",filename);
			exit(1);
		}
		mdstart=readRestart(&p,&structbox,&app);
        init_vel(&p,&structbox,&app);
		fclose(structbox.fprestart);
			

    }

	lcount=vnlist(&p,&structbox);							// Build neighbor list
	structbox.potEne=energyForce(&structbox,&p);			// Compute potential energy
    kinEnergy(&structbox,&p);
    structbox.totEne=structbox.potEne+structbox.KE;
    printf("%d %f ready to start\n",lcount,structbox.potEne);
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("past the barrier\n");
/********************************************************************* INIZIO CODICE MONTE CARLO **************************************************/

	double random,pechng,deltaE;
	double deltaW,ran;
	int s,t,l,spam;
	spam = 0;
    double spam1,spam2=0.0;
    
    // rhok vars
    double spamx,spamy,spamz,rhok,denomrhok,fact,rsq;
    int nx,ny,nz;
    
	// for Long2
	sprintf(filename,"%s-thermo.dat",filepart);
	// for Long
	//sprintf(filename,"%s/%sHard%gK%g_T%g-thermo.dat",argv[3],p.label,structbox.kNsp,structbox.rho0,structbox.krho,structbox.tem);	
	if(p.initRead==1)						
	{
		structbox.fpthermo=fopen(filename,"a");			// Open thermo file and append
		if(structbox.fpthermo==NULL)
            printf("could not find thermo file\n");
	}
	else
	{
		structbox.fpthermo=fopen(filename,"w");			// Open new thermo file
	}
	fclose(structbox.fpthermo);
	//printf("start stepping\n");
    
    
    for(mdstep=mdstart;mdstep<=p.totRun;mdstep++)
    {
        if(mdstep%p.thermo==0)									// Thermo output
		{
            rhok=0.0;
            denomrhok=0.0;
            
//             printf("rho k %f %f\n",rhok,denomrhok); 
//             printf("entering thermo print %d\n",mdstep);
			sprintf(filename,"%s-thermo.dat",filepart);
			structbox.fpthermo=fopen(filename,"a");
			//deltaE=structbox.potEne-energyForce(&structbox,&p);		// Just an additional check. I need to compute virial
            structbox.pressure=structbox.rho*0.66*structbox.KE/p.nAto + structbox.virial*structbox.rho/p.nAto;
            //printf("step %d delta E? %f %f %f\n",mdstep,deltaE,spam1,spam2);
			//if(deltaE<0.0000001)
			//{
// 				printf("%d MC step at thermo %f %f\n",mdstep,structbox.potEne,structbox.delUdelLambda);
                if(structbox.acceptDenom>1)
				fprintf(structbox.fpthermo,"%d	%f	%f %f %f %f %f\n",mdstep,structbox.potEne,structbox.KE,structbox.totEne,2.0*structbox.KE/(3.0*p.nAto),structbox.pressure,structbox.rho);
                else
                fprintf(structbox.fpthermo,"%d	%f	%f %f %f %f %f\n",mdstep,structbox.potEne,structbox.KE,structbox.totEne,2.0*structbox.KE/(3.0*p.nAto),structbox.pressure,structbox.rho);    
				fflush(structbox.fpthermo);
			//}
// 			else
// 			{
// 				printf("Thermo Error! Difference in energy %le at MC:%d\n",deltaE,mdstep);
// 				exit(1);
// 			}
			fclose(structbox.fpthermo);
//             printf("done thermo printing\n");
			
		}

		if(mdstep%p.dump==0)									
		{
			sprintf(filename,"%s-restart%d.dat",filepart,ranspam);
			ranspam=abs(ranspam-1);
			structbox.fprestart=fopen(filename,"w");
			writeRestart(&p,&structbox,&app,mdstep);		// Sputo fuori il restart
			fclose(structbox.fprestart);
			
			sprintf(filename,"%s-LAMMPSdump_%d.dat",filepart,mdstep);
			structbox.fpLAMMPSdump=fopen(filename,"w");
			writeLAMMPSdump(&p,&structbox,&app,mdstep);
			fclose(structbox.fpLAMMPSdump);
			
			
			
		}
        
        
        
        
        
        
        
        integrateGaussian(&structbox,&p,0); // 0 means updated positions and stored old forces
        structbox.potEne=energyForce(&structbox,&p);
        integrateGaussian(&structbox,&p,1); // update velocity using old force and new force
        kinEnergy(&structbox,&p);
        structbox.totEne = structbox.potEne + structbox.KE;
        //printf("%f %f\n",time,structbox.totEne);
        time=(t+1)*p.dt;
//         if(mdstep%1000==0)
//             init_vel(&p,&structbox,&app);
        
		
		if(mdstep%p.dump==0)
		{
			printf("Rank %d: MD step %d completed!\n",rank,mdstep);
		}		
    }

	MPI_Finalize();

}
