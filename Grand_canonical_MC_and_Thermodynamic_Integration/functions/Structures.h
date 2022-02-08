struct SYSTEM
{
	// Thermo variables
	double tem,potEne,pressure,rho;
	double lx,ly,lz;

	//Neighbour list related//
	double maxdispsq,*dispX,*dispY,*dispZ;


	//Coordinates and Neighbourlist Variables//
	double *X,*Y,*Z;
    double *fx,*fy,*fz;
	//int *LIST,*POINT,*NEIGPAR;
	int **NEIGLIST,*NEIGPAR;
	int nlen;

	//Observable/Averages related Variables//
	double virial,Avirial;
    double acceptNum,acceptDenom;
    //crystal unit
    double *xunit, *yunit, *zunit;
	//FILES
	FILE *fpthermo,*fprestart,*fpLAMMPSdump,*fpHisto,*fpHisto2,*fpTraj;

    char ens[3];

};

struct PARAINPUT
{
	int nAto,eqRun,totRun,thermo,dump,seed,initRead;
	double press,dispMax,volMax;
	char label[100];
    
    int ens,startCry;
    double rc,sigma,eps;
    
    double n;
};

struct APPOGGIO
{
	//Random number related//
	unsigned long *mt;
	int mti;
	int **NEIGLISTold;
	float *buff_float;
	double *Xold,*Yold,*Zold,*dispXold,*dispYold,*dispZold;
	double lxold, lyold, lzold;
	double vold,vnew,lnvnew;
	double arg, resfactX, resfactY, resfactZ;
	double penew, pediff;
	double deltaBeta;
};


void inizialized(struct PARAINPUT *p,struct SYSTEM *structbox,struct SYSTEM *tmpbox,struct APPOGGIO *a)
{
	int i;
    

//     if(p->startCry==0)
//         structbox->rho0=0.48;
//     if(p->startCry==1)
//         structbox->rho0=0.455;
    
    structbox->acceptNum=0.0;
    structbox->acceptDenom=0.0;
	structbox->X=calloc(p->nAto,sizeof(double *));
	structbox->Y=calloc(p->nAto,sizeof(double *));
	structbox->Z=calloc(p->nAto,sizeof(double *));
    structbox->fx=calloc(p->nAto,sizeof(double *));
	structbox->fy=calloc(p->nAto,sizeof(double *));
	structbox->fz=calloc(p->nAto,sizeof(double *));
    // for bcc/fcc
    structbox->xunit=calloc(2,sizeof(double *));
	structbox->yunit=calloc(2,sizeof(double *));
	structbox->zunit=calloc(2,sizeof(double *));
    
    structbox->xunit[0]=0.0*alat;
	structbox->yunit[0]=0.0*alat;
	structbox->zunit[0]=0.0*alat;
    structbox->xunit[0]=0.5*alat;
	structbox->yunit[0]=0.5*alat;
	structbox->zunit[0]=0.5*alat;
    
// structbox->xunit[0]=0.0*alat;
// 	structbox->yunit[0]=0.0*alat;
// 	structbox->zunit[0]=0.0*alat;
// 	structbox->xunit[1]=0.5*alat;
// 	structbox->yunit[1]=0.0*alat;
// 	structbox->zunit[1]=0.5*alat;
// 	structbox->xunit[2]=0.25*alat;
// 	structbox->yunit[2]=0.25*alat;
// 	structbox->zunit[2]=0.25*alat;
// 	structbox->xunit[3]=0.75*alat;
// 	structbox->yunit[3]=0.25*alat;
// 	structbox->zunit[3]=0.75*alat;
// 	structbox->xunit[4]=0.0*alat;
// 	structbox->yunit[4]=0.5*alat;
// 	structbox->zunit[4]=0.5*alat;
// 	structbox->xunit[5]=0.5*alat;
// 	structbox->yunit[5]=0.5*alat;
// 	structbox->zunit[5]=0.0*alat;	
// 	structbox->xunit[6]=0.25*alat;
// 	structbox->yunit[6]=0.75*alat;
// 	structbox->zunit[6]=0.75*alat;
// 	structbox->xunit[7]=0.75*alat;
// 	structbox->yunit[7]=0.75*alat;
// 	structbox->zunit[7]=0.25*alat;	
    
	structbox->nlen=(int)(4*M_PI*cube(SKIN+p->rc)*structbox->rho*2.0/3.0);
//     structbox->nlen=p->nAto+1;
    printf("%d neig len\n",structbox->nlen);
	structbox->NEIGPAR=calloc(p->nAto+1,sizeof(int *));
	structbox->NEIGLIST=calloc(p->nAto+1,sizeof(int *));
	for(i=0;i<=p->nAto;i++)
	{
		structbox->NEIGLIST[i]=calloc(structbox->nlen,sizeof(int *));
	}

	structbox->dispX=calloc(p->nAto,sizeof(double *));
	structbox->dispY=calloc(p->nAto,sizeof(double *));
	structbox->dispZ=calloc(p->nAto,sizeof(double *));

	
	
	//Tempbox is used for checking in energy and pressure
	tmpbox->X=calloc(p->nAto,sizeof(double *));
	tmpbox->Y=calloc(p->nAto,sizeof(double *));
	tmpbox->Z=calloc(p->nAto,sizeof(double *));

	tmpbox->NEIGLIST=calloc(p->nAto+1,sizeof(int *));
	for(i=0;i<=p->nAto;i++)
	{
		tmpbox->NEIGLIST[i]=calloc(structbox->nlen,sizeof(int *));
	}

	 
	tmpbox->dispX=calloc(p->nAto,sizeof(double *));
	tmpbox->dispY=calloc(p->nAto,sizeof(double *));
	tmpbox->dispZ=calloc(p->nAto,sizeof(double *));

	a->mti=N+1;
	a->mt=calloc(N,sizeof(unsigned long));
	a->buff_float=calloc(p->nAto,sizeof(float *));

	a->Xold=calloc(p->nAto,sizeof(double *));
	a->Yold=calloc(p->nAto,sizeof(double *));
	a->Zold=calloc(p->nAto,sizeof(double *));

	a->NEIGLISTold=calloc(p->nAto+1,sizeof(int *));
	for(i=0;i<=p->nAto;i++)
	{
		a->NEIGLISTold[i]=calloc(structbox->nlen,sizeof(int *));
	}

	a->dispXold=calloc(p->nAto,sizeof(double *));
	a->dispYold=calloc(p->nAto,sizeof(double *));
	a->dispZold=calloc(p->nAto,sizeof(double *));
  
}
