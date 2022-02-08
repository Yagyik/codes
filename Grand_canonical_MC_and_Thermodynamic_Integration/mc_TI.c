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

#define SWaCon 	7.0495562770
#define SWbCon 	0.6022245584
#define SWlamb 	21.0
#define SWgamm 	1.20
#define SWcut 		1.80
#define SKIN 		0.3
#define minGrow	0.95
#define N 			624
#define minDist2	0.9025
#define LINESIZE 1024
#define alat 2.6

#define NspUppLim	50		// N1 max
#define rhoUppLim	200		// rho max
#define drho	((0.49-0.45)/(1.0*rhoUppLim))	// delta rho for histogram bins			

#include "functions/Structures.h"
#include "functions/FileInput.h"
#include "functions/Random.h"
#include "functions/NeighList.h"
#include "functions/PotentialSW_PAnnihilTI.h"
#include "functions/CreateMove_TI.h"
#include "functions/OrderParam.h"
#include "functions/Metropolis.h"
#include "functions/Cluster.h"
#include "functions/ParallelTempering.h"

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
	char filename[200],**fileList;
	char string[200];
	double **source;
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
		source[i]=calloc(2,sizeof(double *));
		fileList[i]=calloc(500,sizeof(char *));
	}

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

	if(rank==0)									// Writing info on log file
	{
		sprintf(filename,"%s/log.dat",argv[3]);
		fpLog=fopen(filename,"w");
		if(fpLog==NULL) 
		{
			printf("Error: %s not found!\n",filename);
			exit(1);
		}
		writePara(fpLog,&p);
		writesource(fpLog,source,fileList,size);
		fclose(fpLog);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	structbox.lambda=source[rank][0];		// lambda value at which to simulate
	structbox.tem=source[rank][1];			// Temperature

	printf("lambda - %f, temp - %f\n",structbox.lambda,structbox.tem);
	inizialized(&p,&structbox,&tmpbox,&app);
	//countList(&p,source,size);		// Conto quanti N0 e T ho in tutto. Mi serve per fare lo swap tra processori
	
	if(rank==0)
	{
		sprintf(filename,"%s/log.dat",argv[3]);
		fpLog=fopen(filename,"a");
		fprintf(fpLog,"Found %d values of n_0\n",p.nSp);
		fprintf(fpLog,"Found %d values of Temper\n",p.nT);
		fclose(fpLog);
	}
	printf("I got all needed input values\n");
	MPI_Barrier(MPI_COMM_WORLD);	
	
/***************************************************************** CREO E/O LEGGO SITUAZIONE INIZIALE **************************************************/

	int mcStart,lcount,MCstep;
	
	if(p.initRead==0)					// initread 0 = Create atoms
	{
		init_genrand(p.seed,&app);
		creaAtom(&p,&structbox,&app,p.startCry); // create liquid
		
		if(rank==0)
		{
			sprintf(filename,"%s/log.dat",argv[3]);
			fpLog=fopen(filename,"a");
			fprintf(fpLog,"Start from creating atoms\n");
			fprintf(fpLog,"Box: %f	%f	%f\n",structbox.lx,structbox.ly,structbox.lz);
			fprintf(fpLog,"Atom %d: %f	%f	%f\n",p.nAto,structbox.X[p.nAto],structbox.Y[p.nAto],structbox.Z[p.nAto]);
			fclose(fpLog);
		}
		mcStart=0;
	}

	if(p.initRead==1)						// initread 1 = start from STANDARD restart file
	{
	    // for Long_2
	    //printf("restarting sim from restart files\n");
	    // if restart1 gives trouble, read from restart0. ///////
		sprintf(filename,"%s/%sLamb%g_T%g-restart1.dat",argv[3],p.label,structbox.lambda,structbox.tem);
		// for Long
		//sprintf(filename,"%s/%sHard%gK%g_T%g-restart.dat",argv[3],p.label,structbox.kNsp,structbox.rho0,structbox.krho,structbox.tem);
		structbox.fprestart=fopen(filename,"r");
		if(structbox.fprestart==NULL) 
		{
			printf("Error: %s not found!\n",filename);
			exit(1);
		}
		mcStart=readRestart(&p,&structbox,&app);
		//printf("%d\n",mcStart);
		fclose(structbox.fprestart);

		if(mcStart>p.startHist+p.dump)		
		{   
		    // for Long2
			sprintf(filename,"%s/%sLamb%g_T%g-histo.dat",argv[3],p.label,structbox.lambda,structbox.tem);
			// for Long
			//sprintf(filename,"%s/%sHard%gK%g_T%g-histo.dat",argv[3],p.label,structbox.kNsp,structbox.rho0,structbox.krho,structbox.tem);
			structbox.fpHisto=fopen(filename,"r");
			if(structbox.fpHisto==NULL) 
			{
				printf("Error: %s not found at MCstep %d!\n",filename,mcStart);
				exit(1);
			}
			lcount=readHisto(&structbox);	// Reading previous histogram
			fclose(structbox.fpHisto);
			
			
		}
			
		if(rank==0)
		{
			sprintf(filename,"%s/log.dat",argv[3]);
			fpLog=fopen(filename,"a");
			fprintf(fpLog,"Start from STANDARD restart file. Next step: %d\n",mcStart);
			fprintf(fpLog,"Box: %f	%f	%f\n",structbox.lx,structbox.ly,structbox.lz);
			fprintf(fpLog,"Atom %d: %f	%f	%f\n",p.nAto,structbox.X[p.nAto],structbox.Y[p.nAto],structbox.Z[p.nAto]);
			fprintf(fpLog,"Read %d values from histogram\n",lcount);
			fclose(fpLog);
		}
		// read tracker files to find last value
    }

	if(p.initRead==2)						// initread 2 = start from INPUT restart file being in source
	{
		sprintf(filename,"%s",fileList[rank]);
		structbox.fprestart=fopen(filename,"r");
		if(structbox.fprestart==NULL) 
		{
			printf("Error: %s not found!\n",filename);
			exit(1);
		}
		mcStart=readRestart(&p,&structbox,&app);
		
		if(rank==0)
		{
			sprintf(filename,"%s/log.dat",argv[3]);
			fpLog=fopen(filename,"a");
			fprintf(fpLog,"Start from INPUT restart file: %s\n",filename);
			fprintf(fpLog,"MCstep zeroed: %d\n",mcStart);
			fprintf(fpLog,"Box: %f	%f	%f\n",structbox.lx,structbox.ly,structbox.lz);
			fprintf(fpLog,"Atom %d: %f	%f	%f\n",p.nAto,structbox.X[p.nAto],structbox.Y[p.nAto],structbox.Z[p.nAto]);
			fclose(fpLog);
		}
		mcStart=0;
		
	}

	lcount=vnlist(&p,&structbox);							// Build neighbor list
	double storeLambda=structbox.lambda;
    structbox.lambda=0.0;
	structbox.potEne=energy(1,&structbox,&p);			// Compute potential energy
    //structbox.potEne+=pow(structbox.lambda,p.n)*energyAnnihil(1,&structbox,&p);
	/*if(rank==17 || rank==32 || rank==47 || rank==116)
	{
	printf("pot ene %f for rank %d\n",structbox.potEne,rank);
	printf("%f box\n",structbox.lx);
	printf("%f %f %f\n",structbox.X[1],structbox.Y[1],structbox.Z[1]);
	printf("%f %f %f\n",structbox.X[p.nAto],structbox.Y[p.nAto],structbox.Z[p.nAto]);
	}*/
	//largestClus(&structbox,&p);	
	//printf("clust and PE calculated for ranke %d\n",rank);				// Compute solid- and liquid-like particles
	//structbox.W=0.5*structbox.kNsp*sqr(structbox.nSp_max-structbox.nSp_max0)+0.5*structbox.krho*sqr(structbox.rho-structbox.rho0);
	//savePre(&structbox,&p);									// Save configurations
	//printf("U = %f\n",structbox.potEne);
	//printf("Solid part = %f	%f\n",structbox.nSp,structbox.nSp_max);
	//printf("HDL part = %f	%f\n",structbox.nl5,structbox.nl5_max);
	//printf("LDL part = %f	%f\n",structbox.nl4,structbox.nl4_max);
	//printf("W = %f\n",structbox.W);
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("past the barrier\n");
/********************************************************************* INIZIO CODICE MONTE CARLO **************************************************/

	double random,pechng,deltaE;
	double deltaW,ran;
	int s,t,l,spam;
	spam = 0;
    double spam1,spam2=0.0;
	// for Long2
	sprintf(filename,"%s/%sLamb%g_T%g-thermo.dat",argv[3],p.label,storeLambda,structbox.tem);
	// for Long
	//sprintf(filename,"%s/%sHard%gK%g_T%g-thermo.dat",argv[3],p.label,structbox.kNsp,structbox.rho0,structbox.krho,structbox.tem);	
	if(p.initRead==1)						
	{
		structbox.fpthermo=fopen(filename,"a");			// Open thermo file and append
		if(structbox.fpthermo==NULL)
	    printf("%d rank is effed up with nSp_max %f, rho %f and temp %f\n",rank,structbox.nSp_max0,structbox.rho0,structbox.tem);
	}
	else
	{
		structbox.fpthermo=fopen(filename,"w");			// Open new thermo file
	}
	fclose(structbox.fpthermo);
	//printf("start stepping\n");
    
	for(MCstep=mcStart;MCstep<=p.totRun;MCstep++)
	{
	    //printf("%d %d rank step %f\n",rank,MCstep,structbox.lambda);
		random=genrand_real1(&app);
		if(random<0.95)						// 90% move atoms. Update potential energy inside routine
		{
//             printf("attempting MC sweep on %d\n",MCstep);
			pechng=moveAtom(&structbox,&p,MCstep,&app);
			structbox.potEne+=pechng;
//             printf("rank %d has successfully moved %d with %f total %f\n",rank,MCstep,pechng,structbox.potEne);

		}
		else
		{
// 		    printf("attempting vol move on %d\n",MCstep);
			moveVol(&structbox,&p,&app,MCstep);	// 10% move volume Update potential energy and volume inside routine
// 			printf("MC step %d %f %f %f\n",MCstep,structbox.potEne,structbox.pressure,structbox.rho);
		}
        if(MCstep==p.eqRun)
        {
            structbox.lambda=storeLambda;
            AddAtom(&p,&structbox,&app);
            printf("adding atom - %f %f %f\n",structbox.X[p.nAto+1],structbox.Y[p.nAto+1],structbox.Z[p.nAto+1]);
            spam1=energy(1,&structbox,&p);
            spam2=pow(structbox.lambda,p.n)*energyAnnihil(1,&structbox,&p);
            structbox.potEne=spam1+spam2;
            MPI_Barrier(MPI_COMM_WORLD);
        }
		
		//printf("rho swap branch exited\n");
		//printf("rank %d has got here\n",rank);
// 		if(MCstep>p.startHist && MCstep%(100*p.bias)==0)							// Update histogram
// 		{
// 			//printf("entering histograming %d %d\n",rank,MCstep);
// 			//if(structbox.nSp_max>=NspUppLim || structbox.rho>=rhoUppLim*drho+0.45)
// 			if(structbox.nSp_max>=NspUppLim)
// 			{
// 				printf("Error rank %d: nSp or rho too large!\n",rank);
// 				printf("Error rank %d: nSp-rho:%f\n",rank,structbox.nSp_max);
// 			}
// 			else
// 			{
// 				//printf(" %d %f %f %f %f %d\n",(int)((structbox.rho-0.45)/drho),structbox.rho-0.45,(structbox.rho-0.45)/drho,structbox.rho,drho,structbox.nSp_max);
// 				//printf(" %f %f %f %f %f\n",structbox.rho,structbox.rho-0.45,drho,1.0/drho,(structbox.rho-0.45)/drho);
// 				//if(structbox.rho>0.45)
// 				structbox.nSp_Histo[(int)structbox.nSp_max]+=1.0; // convert 0.4-0.5 to index between 0-100
// 				for(i=0;i<NspUppLim;i++)
// 				{
// 				//if(structbox.rho>0.45)
// 				structbox.n_Histo[i]+=structbox.currentHisto[i];//*exp((0.5*structbox.kNsp*sqr(structbox.nSp_max-structbox.nSp_max0))/structbox.tem);
// 				// rho bias will be removed later, but the nSp bias has to be removed now (Saika-Voivod, Poole, 2006, Eq 10)
// 
// 				}
// 			}
// 		}
		
		if(MCstep%p.thermo==0)									// Thermo output
		{
//             printf("entering thermo print %d\n",MCstep);
			sprintf(filename,"%s/%sLamb%g_T%g-thermo.dat",argv[3],p.label,storeLambda,structbox.tem);
			structbox.fpthermo=fopen(filename,"a");
			saveConf(&structbox,&tmpbox,&p);
//             printf("save conf!\n");
            spam1=spam2=0.0;
            spam1=energy(1,&tmpbox,&p);
            if(MCstep>p.eqRun)
                spam2=pow(structbox.lambda,p.n)*energyAnnihil(1,&tmpbox,&p);
            else
                spam2=0.0;
			deltaE=structbox.potEne-spam1-spam2;		// Just an additional check. I need to compute virial
            //printf("step %d delta E? %f %f %f\n",MCstep,deltaE,spam1,spam2);
			if(deltaE<0.0000001)
			{
				structbox.pressure=structbox.rho*structbox.tem+(tmpbox.virial+tmpbox.Avirial)/(structbox.lx*structbox.ly*structbox.lz);
                largestClus(&structbox,&p);
				structbox.q6=q6eval(&structbox,&p);
                if(MCstep>p.eqRun)
                    structbox.delUdelLambda=p.n*pow(structbox.lambda,p.n-1)*energyAnnihil(1,&structbox,&p);
// 				printf("%d MC step at thermo %f %f\n",MCstep,structbox.potEne,structbox.delUdelLambda);
				fprintf(structbox.fpthermo,"%d	%f	%f	%f	%f	%.1f	%.1f	%.1f	%.1f	%.1f	%.1f %f\n",MCstep,structbox.potEne,structbox.pressure,structbox.rho,structbox.q6,structbox.nSp,structbox.nSp_max,structbox.nl4,structbox.nl4_max,structbox.nl5,structbox.nl5_max,structbox.delUdelLambda);
				fflush(structbox.fpthermo);
			}
			else
			{
				sprintf(filename,"%s/log.dat",argv[3]);
				fpLog=fopen(filename,"a");
				fprintf(fpLog,"Thermo Error! Difference in energy %le at MC:%d\n",deltaE,MCstep);
				fclose(fpLog);
				exit(1);
			}
			fclose(structbox.fpthermo);
//             printf("done thermo printing\n");
			
		}

		if(MCstep%p.dump==0)									
		{
			sprintf(filename,"%s/%sLamb%g_T%g-restart%d.dat",argv[3],p.label,structbox.lambda,structbox.tem,ranspam);
			ranspam=abs(ranspam-1);
			structbox.fprestart=fopen(filename,"w");
			writeRestart(&p,&structbox,&app,MCstep);		// Sputo fuori il restart
			fclose(structbox.fprestart);

			if(MCstep>p.startHist)								// Sputo fuori l'istogramma
			{
				sprintf(filename,"%s/%sLamb%g_T%g-histo.dat",argv[3],p.label,structbox.lambda, structbox.tem);
				structbox.fpHisto=fopen(filename,"w");
				//writeHisto(&structbox);
				fclose(structbox.fpHisto);
				
			}
			
			sprintf(filename,"%s/%sLamb%g_T%g-traj.bin",argv[3],p.label,structbox.lambda,structbox.tem);
			structbox.fpTraj=fopen(filename,"a");
			writeDump(&p,&structbox,&app,MCstep);		// Sputo fuori il dump
			fclose(structbox.fpTraj);
			
			sprintf(filename,"%s/%sLamb%g_T%g-LAMMPSdump_%d.dat",argv[3],p.label,structbox.lambda,structbox.tem,MCstep);
			structbox.fpLAMMPSdump=fopen(filename,"w");
			writeLAMMPSdump(&p,&structbox,&app,MCstep);
			fclose(structbox.fpLAMMPSdump);
			
			
			
		}
		
		if(MCstep%p.dump==0)
		{
			printf("Rank %d: MC step %d completed!\n",rank,MCstep);
		}		
	}

	MPI_Finalize();

}
