void creaAtom(struct PARAINPUT *p,struct SYSTEM *structbox,struct APPOGGIO *a)
{
	int i,j,k;
	int atomdone,atomnotdone;
	double dummyX,dummyY,dummyZ;
	double rX,rY,rZ,r;
	
	structbox->lx=structbox->ly=structbox->lz=pow((p->nAto/structbox->rho),(1.0/3.0));
    //printf("created box with lx %f\n",structbox->lx);
	for(i=1;i<=p->nAto;i++)
	{
		atomdone=0;
		while(atomdone==0)
		{
			dummyX=structbox->lx*genrand_real1(&(*a));
			dummyY=structbox->ly*genrand_real1(&(*a));
			dummyZ=structbox->lz*genrand_real1(&(*a));

			atomnotdone=0;
			for(j=1;j<=i-1;j++)
			{
				rX=(dummyX-structbox->X[j])-structbox->lx*lround((dummyX-structbox->X[j])/structbox->lx);
				rY=(dummyY-structbox->Y[j])-structbox->ly*lround((dummyY-structbox->Y[j])/structbox->ly);
				rZ=(dummyZ-structbox->Z[j])-structbox->lz*lround((dummyZ-structbox->Z[j])/structbox->lz);

				r=rX*rX+rY*rY+rZ*rZ;
				if(r<minDist2)
				{
					atomnotdone=1;
				}
			}
			atomdone=1-atomnotdone;

			structbox->X[i]=dummyX;
			structbox->Y[i]=dummyY;
			structbox->Z[i]=dummyZ;
		}
	}

	// Put origin in (0,0,0)
	for(i=1;i<=p->nAto;i++)
	{
		structbox->X[i]-=structbox->lx*0.5;
		structbox->Y[i]-=structbox->ly*0.5;
		structbox->Z[i]-=structbox->lz*0.5;
	}
}

void update(int par, struct SYSTEM *box,struct PARAINPUT *p)
{
	double displacesq;

	displacesq=sqr(box->dispX[par])+sqr(box->dispY[par])+sqr(box->dispZ[par]);
	box->maxdispsq=MAX(box->maxdispsq,displacesq);

	if(box->maxdispsq>sqr(0.5*SKIN))
	{
		box->maxdispsq=0.0;
		vnlist(&(*p),&(*box));
	}
}

double moveAtom(struct SYSTEM *box, struct PARAINPUT *p,int mc,struct APPOGGIO *a)
{
	int par,i,t;
	double delx,dely,delz;
	double Xold,Yold,Zold;
	double ene_before=0.0, ene_after=0.0;
	double ran,prob,pechng;
    int flag=0;
	pechng=0.0;

	for(t=1;t<=p->nAto;t++)
	{
		delx=(genrand_real1(&(*a))-0.5)*p->dispMax;
		dely=(genrand_real1(&(*a))-0.5)*p->dispMax;
		delz=(genrand_real1(&(*a))-0.5)*p->dispMax;

		par=(int)(genrand_real1(&(*a))*(p->nAto));
		if(par==(p->nAto))
		{
			par=p->nAto-1;
		}

		Xold=box->X[par];
		Yold=box->Y[par];
		Zold=box->Z[par];
		ene_before=oldparenergy(par,&(*box));	// Energy before displacement

		box->X[par]+=delx; 
		box->Y[par]+=dely;
		box->Z[par]+=delz;

		vnlistpar(par,&(*box),delx,dely,delz,&(*p));		// Build neighbor list
		ene_after=newparenergy(par,&(*box));	// Energy after displacement

		//Metropolis Technique
		ran=genrand_real1(&(*a));
		prob=exp(-(ene_after-ene_before)/box->tem);
       // printf("prob is %f for %d\n",prob,t);
		if(isnan(prob)!= 0 || isinf(prob)!= 0)
		{
			box->X[par]=Xold;		// Rejected
			box->Y[par]=Yold;
			box->Z[par]=Zold;
			printf("DISPLACE %d: Probability coming out to be nan or inf\t%f\t%f\n",mc,ene_before,ene_after);
			exit(0);
		}
		else
		{
		    //printf("move accept or reject here\n");
			if(ran>prob)
			{
				box->X[par]=Xold; 	// Rejected
				box->Y[par]=Yold;
				box->Z[par]=Zold;
			}	
			else
			{
			    flag+=1;
				box->dispX[par]+=delx;
				box->dispY[par]+=dely;
				box->dispZ[par]+=delz;

				pechng+=ene_after-ene_before;
				box->NEIGPAR[0]=0;
				update(par,&(*box),&(*p));
			}
		}
	}
	//printf("%d accepted moves\n",flag);
	return(pechng);
}

void moveVol(struct SYSTEM *structbox,struct PARAINPUT *p,struct APPOGGIO *a)
{
	int i,j,lcount;
	double arg,prob,ran;

	a->lxold=structbox->lx;
	a->lyold=structbox->ly;
	a->lzold=structbox->lz;
	a->vold=a->lxold*a->lyold*a->lzold;
	a->lnvnew=log(a->vold)+(genrand_real1(&(*a))-0.5)*p->volMax;
	a->vnew=exp(a->lnvnew);

	structbox->lx=structbox->ly=structbox->lz=pow(a->vnew,1.0/3.0);
	a->resfactX=(double)(structbox->lx/a->lxold);
	a->resfactY=(double)(structbox->ly/a->lyold);
	a->resfactZ=(double)(structbox->lz/a->lzold);

	for(i=1;i<=p->nAto;i++)				// Move the volume
	{
		a->Xold[i]=structbox->X[i];
		a->Yold[i]=structbox->Y[i];
		a->Zold[i]=structbox->Z[i];
		structbox->X[i]=a->Xold[i]*a->resfactX;
		structbox->Y[i]=a->Yold[i]*a->resfactY;
		structbox->Z[i]=a->Zold[i]*a->resfactZ;

		for(j=0;j<structbox->nlen;j++)
		{
			a->NEIGLISTold[i][j]=structbox->NEIGLIST[i][j];
		}

		a->dispXold[i]=structbox->dispX[i];
		a->dispYold[i]=structbox->dispY[i];
		a->dispZold[i]=structbox->dispZ[i];
	}	
	
	lcount=vnlist(&(*p),&(*structbox));
	a->penew=energy(1,&(*structbox),&(*p));
	//printf("%f delta E\n",a->penew-structbox->potEne);
	//printf("%f press and deltaV %f %f\n",p->press,a->vnew-a->vold,a->vold);
	arg=-((a->penew-structbox->potEne)+p->press*(a->vnew-a->vold))/structbox->tem+((p->nAto+1)*log(a->vnew/a->vold));
	prob=exp(arg);
	ran=genrand_real1(&(*a));

	if(ran>prob)
	{
		structbox->lx=a->lxold;				// Rejected
		structbox->ly=a->lyold;
		structbox->lz=a->lzold;
       // printf("Move vol rejected\n");
		for(i=1;i<=p->nAto;i++)				
		{
			structbox->X[i]=a->Xold[i];
			structbox->Y[i]=a->Yold[i];
			structbox->Z[i]=a->Zold[i];

			for(j=0;j<structbox->nlen;j++)
			{
				structbox->NEIGLIST[i][j]=a->NEIGLISTold[i][j];
			}

			structbox->dispX[i]=a->dispXold[i];
			structbox->dispY[i]=a->dispYold[i];
			structbox->dispZ[i]=a->dispZold[i];
			
		}
	}
	else						
	{
	   // printf("volume updated accepted %f %f\n",a->penew-structbox->potEne,a->vnew-a->vold);
		structbox->potEne=a->penew;				// Update potential energy
		structbox->rho=p->nAto/(structbox->lx*structbox->ly*structbox->lz);		// Update density
	}
}
