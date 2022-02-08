void countList(struct PARAINPUT *p, double **source, int size)
{
	int i,j,k,n;
	double n0,rho,T;
	
	n0=source[0][0];
	rho=source[0][4];
	T=source[0][6];
	
	p->nT=1;
	p->nSp=1;
	p->nrho=1;
	for(i=1;i<size;i++)
	{
		if(n0==source[i][0] && rho==source[i][4])
		{
			p->nT++;
		}
	}
	
	for(i=1;i<size;i++)
	{
		if(T==source[i][6] && n0==source[i][0])
		{
			p->nrho++;
		}
	}	
	
	for(i=1;i<size;i++)
	{
		if(T==source[i][6] && rho==source[i][4])
		{
		    //printf("%f %f %f %f\n",T,source[i][6],rho,source[1][4]);
			p->nSp++;
		}
	}
	//printf("%d %d %d %d\n",size,p->nSp,p->nrho,p->nT);
	if(p->nrho*p->nSp*p->nT!=size)
	{
		printf("PT warning: grid of %d values, but using %d cores!!\n",p->nrho*p->nSp*p->nT,size);
	}
	
	p->procList=calloc(p->nSp,sizeof(int *));
	for(i=0;i<p->nSp;i++)
	{
		p->procList[i]=calloc(p->nrho,sizeof(int *));
		for(j=0;j<p->nrho;j++)
		{
			p->procList[i][j]=calloc(p->nT,sizeof(int *));
		}
	}
	
	n=0;
	for(k=0;k<p->nT;k++)
	{
	for(i=0;i<p->nSp;i++)
	{
	for(j=0;j<p->nrho;j++)
	{
		p->procList[i][j][k]=n;		// This order follows the makeSource.sh's one!!!
		//printf(" ijk %d %d %d - n %d\n",i,j,k,n);
		n++;
	}
	}
	}		
}

void TswapSend(struct SYSTEM *structbox, int proc,struct APPOGGIO *app)
{
	double Vswap;
	Vswap=structbox->lx*structbox->ly*structbox->lz;
	MPI_Send(&Vswap,1,MPI_DOUBLE,proc,0,MPI_COMM_WORLD);											// tag 0 = Volume
	MPI_Send(&structbox->potEne,1,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);							// tag 1 = U
	MPI_Send(&structbox->nSp_max,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD);									// tag 23 = n1
	MPI_Recv(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);			// tag 2 = esitoSwap
}
	
void TswapReceive(struct SYSTEM *structbox, int proc, double **source,struct APPOGGIO *app,struct PARAINPUT *p)
{	
	double Vswap,Uswap,ran;
	int flag;
	double n0Swap,n0uswap,n0lswap;
	MPI_Recv(&Vswap,1,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 0 = Volume
	MPI_Recv(&Uswap,1,MPI_DOUBLE,proc,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 1 = U
	MPI_Recv(&n0Swap,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 23 = n1
	
	app->deltaBeta=(1.0/structbox->tem-1.0/source[proc][6]); // find out what this does
	app->vold=structbox->lx*structbox->ly*structbox->lz-Vswap;
	ran=genrand_real1(&(*app));	// Attento: questa formula e' vera se le finestre sono esattamente uguali come k,N0!
	
		// for swap to happen, the bins should be overlapping
	// n0L and n0U are accessible via source[rank][2] and source[rank][3]
	// make a flag conditional on source rank and structbox.n0U, structbox.n0L
	flag=0;
	if(structbox->nSp_max >= source[proc][2] && structbox->nSp_max <= source[proc][3])
	{
	    if(n0Swap >= structbox->n0L && n0Swap <=structbox->n0U)
	    flag=1;
	}
	if(ran<exp(app->deltaBeta*(p->press*app->vold+(structbox->potEne-Uswap))) && flag==1)
	{
		app->esitoSwap=1;
		MPI_Send(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
	else
	{
		app->esitoSwap=0;
		MPI_Send(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
}


void N0SwapSend(struct SYSTEM *structbox, int proc,struct APPOGGIO *app)
{
	//printf("swap send to %d\n",proc);
	MPI_Send(&structbox->nSp_max,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD);									// tag 23 = n1
	MPI_Recv(&app->esitoSwap,1,MPI_DOUBLE,proc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);			// tag 2 = esitoSwap
}
	
void N0SwapRiceive(struct SYSTEM *structbox, int proc, double **source,struct APPOGGIO *app,struct PARAINPUT *p)
{	
	double ran,wN,w0;
	double n0Swap,flag;
	
	MPI_Recv(&n0Swap,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 23 = n1
	//printf("swap recv from %d\n",proc);
	//wN=0.5*source[proc][1]*sqr(structbox->nSp_max-source[proc][0])+0.5*structbox->kNsp*sqr(n0Swap-structbox->nSp_max0);
	//w0=0.5*source[proc][1]*sqr(n0Swap-source[proc][0])+0.5*structbox->kNsp*sqr(structbox->nSp_max-structbox->nSp_max0);
	ran=genrand_real1(&(*app));
	// for swap to happen, the bins should be overlapping
	// n0L and n0U are accessible via source[rank][2] and source[rank][3]
	// make a flag conditional on source rank and structbox.n0U, structbox.n0L
	flag=0;
	if(structbox->nSp_max >= source[proc][2] && structbox->nSp_max <= source[proc][3])
	{
	    printf("one way swap works %f %f %f %f\n",structbox->nSp_max,n0Swap,source[proc][2],structbox->n0L);
	    if(n0Swap >= structbox->n0L && n0Swap <=structbox->n0U)
	    {
	    flag=1;
	    printf("other way swap also works\n");
	    }
	}
	if(flag==1)
	//if(flag==1)
	{
		app->esitoSwap=1;
		MPI_Send(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
	else
	{
		app->esitoSwap=0;
		MPI_Send(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
}

void rhoSwapSend(struct SYSTEM *structbox, int proc,struct APPOGGIO *app)
{

	MPI_Send(&structbox->rho,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD);									// tag 23 = n1
	MPI_Recv(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);			// tag 2 = esitoSwap

}
	
void rhoSwapRiceive(struct SYSTEM *structbox, int proc, double **source,struct APPOGGIO *app,struct PARAINPUT *p)
{	
	double ran,wN,w0;
	double n1Swap;

	MPI_Recv(&n1Swap,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 23 = n1 (rho)
	
	wN=0.5*source[proc][3]*sqr(structbox->rho-source[proc][2])+0.5*structbox->krho*sqr(n1Swap-structbox->rho0);
	w0=0.5*source[proc][3]*sqr(n1Swap-source[proc][2])+0.5*structbox->krho*sqr(structbox->rho-structbox->rho0);
	ran=genrand_real1(&(*app));
	
	if(ran<exp(-(wN-w0)/structbox->tem))
	{
		app->esitoSwap=1;
		MPI_Send(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
	else
	{
		app->esitoSwap=0;
		MPI_Send(&app->esitoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}

}


void SendTmpBox(struct SYSTEM *structbox,struct SYSTEM *tmpbox,struct PARAINPUT *p,int proc)
{
	int i;
	
	tmpbox->lx=structbox->lx;
	tmpbox->ly=structbox->ly;
	tmpbox->lz=structbox->lz;
	tmpbox->nSp_maxTrack=structbox->nSp_maxTrack;
	tmpbox->rhoTrack=structbox->rhoTrack;
	tmpbox->temTrack=structbox->temTrack;
	for(i=1;i<=p->nAto;i++)
	{
		tmpbox->X[i]=structbox->X[i];
		tmpbox->Y[i]=structbox->Y[i];
		tmpbox->Z[i]=structbox->Z[i];
		tmpbox->dispX[i]=structbox->dispX[i];
		tmpbox->dispY[i]=structbox->dispY[i];
		tmpbox->dispZ[i]=structbox->dispZ[i];
		
	}

	MPI_Send(tmpbox->X,p->nAto+1,MPI_DOUBLE,proc,3,MPI_COMM_WORLD);		// tag 3 = tmpbox->X
	MPI_Send(tmpbox->Y,p->nAto+1,MPI_DOUBLE,proc,4,MPI_COMM_WORLD);		// tag 4 = tmpbox->Y
	MPI_Send(tmpbox->Z,p->nAto+1,MPI_DOUBLE,proc,5,MPI_COMM_WORLD);		// tag 5 = tmpbox->Z
	//printf("Proc 0; Send atom: %f	%f	%f\n",tmpbox->X[10],tmpbox->Y[1],tmpbox->Z[p->nAto]);
	
	MPI_Send(&tmpbox->lx,1,MPI_DOUBLE,proc,6,MPI_COMM_WORLD);									// tag 6 = tmpbox->lx
	MPI_Send(&tmpbox->ly,1,MPI_DOUBLE,proc,7,MPI_COMM_WORLD);									// tag 7 = tmpbox->ly
	MPI_Send(&tmpbox->lz,1,MPI_DOUBLE,proc,8,MPI_COMM_WORLD);									// tag 8 = tmpbox->lz
	//printf("Proc 0;Send Box: %f	%f	%f\n",tmpbox->lx,tmpbox->ly,tmpbox->lz);
		
	MPI_Send(tmpbox->dispX,p->nAto+1,MPI_DOUBLE,proc,10,MPI_COMM_WORLD);		// tag 10 = tmpbox->dispX
	MPI_Send(tmpbox->dispY,p->nAto+1,MPI_DOUBLE,proc,11,MPI_COMM_WORLD);		// tag 11 = tmpbox->dispY
	MPI_Send(tmpbox->dispZ,p->nAto+1,MPI_DOUBLE,proc,12,MPI_COMM_WORLD);		// tag 12 = tmpbox->dispZ
	//printf("Proc 0;Send Disp: %f	%f	%f\n",tmpbox->dispX[10],tmpbox->dispY[1],tmpbox->dispX[p->nAto-1]);
	
	MPI_Send(&tmpbox->nSp_maxTrack,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD);									// tag 23 = tmpbox->nSp_maxTrack
	MPI_Send(&tmpbox->rhoTrack,1,MPI_DOUBLE,proc,24,MPI_COMM_WORLD);									// tag 24 = tmpbox->rhoTrack
	MPI_Send(&tmpbox->temTrack,1,MPI_DOUBLE,proc,25,MPI_COMM_WORLD);									// tag 25 = tmpbox->temTrack
	// Sent all the info on tmpbox. Now I receive in structbox
	
	MPI_Recv(structbox->X,p->nAto+1,MPI_DOUBLE,proc,13,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 24 = structbox->X
	MPI_Recv(structbox->Y,p->nAto+1,MPI_DOUBLE,proc,14,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 25 = structbox->Y	
	MPI_Recv(structbox->Z,p->nAto+1,MPI_DOUBLE,proc,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 26 = structbox->Z
	//printf("Proc 0; Receive atom: %f	%f	%f\n",structbox->X[10],structbox->Y[1],structbox->Z[p->nAto]);		   		 

	MPI_Recv(&structbox->lx,1,MPI_DOUBLE,proc,16,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 27 = structbox->lx
	MPI_Recv(&structbox->ly,1,MPI_DOUBLE,proc,17,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 28 = structbox->ly
	MPI_Recv(&structbox->lz,1,MPI_DOUBLE,proc,18,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 29 = structbox->lz
	//printf("Proc 0; Receive Box: %f	%f	%f\n",structbox->lx,structbox->ly,structbox->lz);

	MPI_Recv(structbox->dispX,p->nAto+1,MPI_DOUBLE,proc,20,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 20 = structbox->dispX
	MPI_Recv(structbox->dispY,p->nAto+1,MPI_DOUBLE,proc,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 21 = structbox->dispY	
	MPI_Recv(structbox->dispZ,p->nAto+1,MPI_DOUBLE,proc,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 22 = structbox->dispZ
	//printf("Proc 0; Receive Disp: %f	%f	%f\n",structbox->dispX[10],structbox->dispY[1],structbox->dispX[p->nAto-1]);		    
	
	MPI_Recv(&structbox->nSp_maxTrack,1,MPI_DOUBLE,proc,26,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 26 = structbox->nSp_maxTrack
	MPI_Recv(&structbox->rhoTrack,1,MPI_DOUBLE,proc,27,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 27 = structbox->rhoTrack
	MPI_Recv(&structbox->temTrack,1,MPI_DOUBLE,proc,28,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 28 = structbox->temTrack
	//printf("Proc 0;Send Box: %f	%f	%f\n",tmpbox->lx,tmpbox->ly,tmpbox->lz);		    
}

void RiceiveTmpBox(struct SYSTEM *structbox,struct SYSTEM *tmpbox,struct PARAINPUT *p,int proc)
{
	int i;
	
	tmpbox->lx=structbox->lx;
	tmpbox->ly=structbox->ly;
	tmpbox->lz=structbox->lz;
	tmpbox->nSp_maxTrack=structbox->nSp_maxTrack;
	tmpbox->rhoTrack=structbox->rhoTrack;
	tmpbox->temTrack=structbox->temTrack;	
	for(i=1;i<=p->nAto;i++)
	{
		tmpbox->X[i]=structbox->X[i];
		tmpbox->Y[i]=structbox->Y[i];
		tmpbox->Z[i]=structbox->Z[i];
		tmpbox->dispX[i]=structbox->dispX[i];
		tmpbox->dispY[i]=structbox->dispY[i];
		tmpbox->dispZ[i]=structbox->dispZ[i];
	}
	
	MPI_Recv(structbox->X,p->nAto+1,MPI_DOUBLE,proc,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 3 = tmpbox->X
	MPI_Recv(structbox->Y,p->nAto+1,MPI_DOUBLE,proc,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 4 = tmpbox->Y	
	MPI_Recv(structbox->Z,p->nAto+1,MPI_DOUBLE,proc,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 5 = tmpbox->Z
	//printf("Proc 1; Receive atom: %f	%f	%f\n",structbox->X[10],structbox->Y[1],structbox->Z[p->nAto]);		   		 

	MPI_Recv(&structbox->lx,1,MPI_DOUBLE,proc,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 6 = tmpbox->lx
	MPI_Recv(&structbox->ly,1,MPI_DOUBLE,proc,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 7 = tmpbox->ly
	MPI_Recv(&structbox->lz,1,MPI_DOUBLE,proc,8,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 8 = tmpbox->lz
	//printf("Proc 1; Receive Box: %f	%f	%f\n",structbox->lx,structbox->ly,structbox->lz);	
	
	MPI_Recv(structbox->dispX,p->nAto+1,MPI_DOUBLE,proc,10,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 10 = structbox->dispX
	MPI_Recv(structbox->dispY,p->nAto+1,MPI_DOUBLE,proc,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 11 = structbox->dispY
	MPI_Recv(structbox->dispZ,p->nAto+1,MPI_DOUBLE,proc,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 12 = structbox->dispZ	
	//printf("Proc 1; Receive Disp: %f	%f	%f\n",structbox->dispX[10],structbox->dispY[1],structbox->dispX[p->nAto-1]);

	MPI_Recv(&structbox->nSp_maxTrack,1,MPI_DOUBLE,proc,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 23 = tmpbox->nSp_maxtrack
	MPI_Recv(&structbox->rhoTrack,1,MPI_DOUBLE,proc,24,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 24 = tmpbox->rhoTrack
	MPI_Recv(&structbox->temTrack,1,MPI_DOUBLE,proc,25,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 25 = tmpbox->temTrack
	
	//	I received all the info on structbox. Now I send tempbox

	MPI_Send(tmpbox->X,p->nAto+1,MPI_DOUBLE,proc,13,MPI_COMM_WORLD);		// tag 13 = structbox->X
	MPI_Send(tmpbox->Y,p->nAto+1,MPI_DOUBLE,proc,14,MPI_COMM_WORLD);		// tag 14 = structbox->Y
	MPI_Send(tmpbox->Z,p->nAto+1,MPI_DOUBLE,proc,15,MPI_COMM_WORLD);		// tag 15 = structbox->Z
	//printf("Proc 1; Send atom: %f	%f	%f\n",tmpbox->X[10],tmpbox->Y[1],tmpbox->Z[p->nAto]);
	
	MPI_Send(&tmpbox->lx,1,MPI_DOUBLE,proc,16,MPI_COMM_WORLD);				// tag 16 = structbox->lx
	MPI_Send(&tmpbox->ly,1,MPI_DOUBLE,proc,17,MPI_COMM_WORLD);				// tag 17 = structbox->ly
	MPI_Send(&tmpbox->lz,1,MPI_DOUBLE,proc,18,MPI_COMM_WORLD);				// tag 18 = structbox->lz
	//printf("Proc 1; Send Box: %f	%f	%f\n",tmpbox->lx,tmpbox->ly,tmpbox->lz);
	
	MPI_Send(tmpbox->dispX,p->nAto+1,MPI_DOUBLE,proc,20,MPI_COMM_WORLD);		// tag 20 = structbox->dispX
	MPI_Send(tmpbox->dispY,p->nAto+1,MPI_DOUBLE,proc,21,MPI_COMM_WORLD);		// tag 21 = structbox->dispY
	MPI_Send(tmpbox->dispZ,p->nAto+1,MPI_DOUBLE,proc,22,MPI_COMM_WORLD);		// tag 22 = structbox->dispZ

	MPI_Send(&tmpbox->nSp_maxTrack,1,MPI_DOUBLE,proc,26,MPI_COMM_WORLD);				// tag 26 = structbox->nSp_maxTrack
	MPI_Send(&tmpbox->rhoTrack,1,MPI_DOUBLE,proc,27,MPI_COMM_WORLD);				// tag 27 = structbox->rhoTrack
	MPI_Send(&tmpbox->temTrack,1,MPI_DOUBLE,proc,28,MPI_COMM_WORLD);				// tag 28 = structbox->temTrack

	//printf("Proc 1; Send Disp: %f	%f	%f\n",tmpbox->dispX[10],tmpbox->dispY[1],tmpbox->dispX[p->nAto-1]);				
}
