void creaAtom(struct PARAINPUT *p,struct SYSTEM *structbox,struct APPOGGIO *a, int flag)
{
	int i,j,k;
	int atomdone,atomnotdone;
	double dummyX,dummyY,dummyZ;
	double rX,rY,rZ,r;
	
	structbox->lx=structbox->ly=structbox->lz=pow((p->nAto/structbox->rho0),(1.0/3.0));
    //printf("created box with lx %f\n",structbox->lx);
    if(flag==0)
    {
        for(i=1;i<=p->nAto+1;i++)
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
    }
	
    
    // create crystal
    int spamcount,m;
    if(flag == 1)
	{	
        spamcount = 0;
		if(p->nAto != 512) 
		{
		printf("making crystal from scratch needs 512 atoms");
		exit(1);
		}
		for(i=0;i<4;i++)
		{
            for(j=0;j<4;j++)
            {
                for(k=0;k<4;k++)
                {
                    for(m=0;m<8;m++)
                    {
                    structbox->X[spamcount] = i*alat + structbox->xunit[m];
                    structbox->Y[spamcount] = j*alat + structbox->yunit[m]; 
                    structbox->Z[spamcount] = k*alat + structbox->zunit[m];            
                    spamcount++;
                    }
                }
            }
        }
	}
    
        
	// Put origin in (0,0,0)
	for(i=1;i<=p->nAto+1;i++)
	{
		structbox->X[i]-=structbox->lx*0.5;
		structbox->Y[i]-=structbox->ly*0.5;
		structbox->Z[i]-=structbox->lz*0.5;
	}
	structbox->rho=p->nAto/(structbox->lx*structbox->ly*structbox->lz);
}

void AddAtom(struct PARAINPUT *p,struct SYSTEM *structbox,struct APPOGGIO *a)
{
    int i,j,k;
	int atomdone,atomnotdone;
	double dummyX,dummyY,dummyZ;
	double rX,rY,rZ,r;
    
    atomdone=0;
    i=p->nAto+1;
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
    double spam1,spam2=0.0;
	for(t=1;t<=p->nAto+1;t++)
	{
		delx=(genrand_real1(&(*a))-0.5)*p->dispMax;
		dely=(genrand_real1(&(*a))-0.5)*p->dispMax;
		delz=(genrand_real1(&(*a))-0.5)*p->dispMax;
        ene_before=0.0;
        ene_after=0.0;
        spam1=0.0;
        spam2=0.0;
		par=(int)(genrand_real1(&(*a))*(p->nAto+2));
		if(par==0)
		{
			par=1;
		}
		if(par==(p->nAto+2))
		{
			par=p->nAto+1;
		}
        //printf("step %d moving atom %d\n",mc,par);
		Xold=box->X[par];
		Yold=box->Y[par];
		Zold=box->Z[par];
        //printf("obtained positions of %d\n",par);
        if(par<=p->nAto)
            spam1=oldparenergy(par,&(*box));	// Energy before displacement
        
        //spam2=pow(box->lambda,p->n)*oldparenergyAnnihil(par,&(*box),&(*p));
        if(mc>p->eqRun)
            spam2=pow(box->lambda,p->n)*energyAnnihil(0,&(*box),&(*p));
        ene_before=spam1 + spam2;
//         if(par==p->nAto+1)
         //printf("new ene before %f annihil ene %f but lambda %f\n",spam1,spam2,box->lambda);
		box->X[par]+=delx; 
		box->Y[par]+=dely;
		box->Z[par]+=delz;
        if(par<=p->nAto)
        {
            vnlistpar(par,&(*box),delx,dely,delz,&(*p));		// Build neighbor list
            spam1=newparenergy(par,&(*box));	// Energy after displacement
        }
        //spam2=pow(box->lambda,p->n)*newparenergyAnnihil(par,&(*box),&(*p));
        if(mc>p->eqRun)
            spam2=pow(box->lambda,p->n)*energyAnnihil(0,&(*box),&(*p));
        ene_after=spam1+spam2;
//         if(par==p->nAto+1)
         //printf("%d par new ene after %f annihil ene %f but lambda %f\n",par,spam1,spam2,box->lambda);
		//Metropolis Technique
		ran=genrand_real1(&(*a));
		prob=exp(-(ene_after-ene_before)/box->tem);
        //printf("prob is %f for %d\n",prob,t);
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
                if(par<=p->nAto)
                    updatehiSiTi(par,Xold,Yold,Zold,&(*box));
                
                //updatehiSiTiAnnihil(par,Xold,Yold,Zold,&(*box),&(*p));
				box->NEIGPAR[0]=0;
                if(par<=p->nAto)
                    update(par,&(*box),&(*p));
			}
		}
		//printf("managed to move with %f\n",pechng);
	}
	//printf("%d accepted moves\n",flag);
	return(pechng);
}

void moveVol(struct SYSTEM *structbox,struct PARAINPUT *p,struct APPOGGIO *a,int mc)
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

	for(i=1;i<=p->nAto+1;i++)				// Move the volume
	{
		a->Xold[i]=structbox->X[i];
		a->Yold[i]=structbox->Y[i];
		a->Zold[i]=structbox->Z[i];
		structbox->X[i]=a->Xold[i]*a->resfactX;
		structbox->Y[i]=a->Yold[i]*a->resfactY;
		structbox->Z[i]=a->Zold[i]*a->resfactZ;
		a->hiold[i]=structbox->hi[i];
		a->siXold[i]=structbox->siX[i];
		a->siYold[i]=structbox->siY[i];
		a->siZold[i]=structbox->siZ[i];
		a->TiXXold[i]=structbox->TiXX[i];
		a->TiXYold[i]=structbox->TiXY[i];
		a->TiXZold[i]=structbox->TiXZ[i];
		a->TiYYold[i]=structbox->TiYY[i];
		a->TiYZold[i]=structbox->TiYZ[i];
		a->TiZZold[i]=structbox->TiZZ[i];
        
        a->Ahiold[i]=structbox->Ahi[i];
		a->AsiXold[i]=structbox->AsiX[i];
		a->AsiYold[i]=structbox->AsiY[i];
		a->AsiZold[i]=structbox->AsiZ[i];
		a->ATiXXold[i]=structbox->ATiXX[i];
		a->ATiXYold[i]=structbox->ATiXY[i];
		a->ATiXZold[i]=structbox->ATiXZ[i];
		a->ATiYYold[i]=structbox->ATiYY[i];
		a->ATiYZold[i]=structbox->ATiYZ[i];
		a->ATiZZold[i]=structbox->ATiZZ[i];
        if(i<=p->nAto)
        {
            for(j=0;j<structbox->nlen;j++)
            {
                a->NEIGLISTold[i][j]=structbox->NEIGLIST[i][j];
            }
        }
		

		a->dispXold[i]=structbox->dispX[i];
		a->dispYold[i]=structbox->dispY[i];
		a->dispZold[i]=structbox->dispZ[i];
	}	
	
	lcount=vnlist(&(*p),&(*structbox));
	a->penew=energy(1,&(*structbox),&(*p));
    //printf("%f energy old\n",a->penew);
    if(mc>p->eqRun)
        a->penew+=pow(structbox->lambda,p->n)*energyAnnihil(1,&(*structbox),&(*p));
	//printf("%f delta E vol change %f\n",a->penew,a->penew-structbox->potEne);
	//printf("%f press and deltaV %f %f\n",p->press,a->vnew-a->vold,a->vold);
	arg=-((a->penew-structbox->potEne)+p->press*(a->vnew-a->vold))/structbox->tem+((p->nAto+1)*log(a->vnew/a->vold));
	prob=exp(arg);
    //printf("%f prob or prob not?\n",prob);
	ran=genrand_real1(&(*a));

	if(ran>prob)
	{
		structbox->lx=a->lxold;				// Rejected
		structbox->ly=a->lyold;
		structbox->lz=a->lzold;
       // printf("Move vol rejected\n");
		for(i=1;i<=p->nAto+1;i++)				
		{
            structbox->hi[i]=a->hiold[i];
			structbox->siX[i]=a->siXold[i];
			structbox->siY[i]=a->siYold[i];
			structbox->siZ[i]=a->siZold[i];
			structbox->TiXX[i]=a->TiXXold[i];
			structbox->TiXY[i]=a->TiXYold[i];
			structbox->TiXZ[i]=a->TiXZold[i];
			structbox->TiYY[i]=a->TiYYold[i];
			structbox->TiYZ[i]=a->TiYZold[i];
			structbox->TiZZ[i]=a->TiZZold[i];
            
			structbox->Ahi[i]=a->Ahiold[i];
			structbox->AsiX[i]=a->AsiXold[i];
			structbox->AsiY[i]=a->AsiYold[i];
			structbox->AsiZ[i]=a->AsiZold[i];
			structbox->ATiXX[i]=a->ATiXXold[i];
			structbox->ATiXY[i]=a->ATiXYold[i];
			structbox->ATiXZ[i]=a->ATiXZold[i];
			structbox->ATiYY[i]=a->ATiYYold[i];
			structbox->ATiYZ[i]=a->ATiYZold[i];
			structbox->ATiZZ[i]=a->ATiZZold[i];
			structbox->X[i]=a->Xold[i];
			structbox->Y[i]=a->Yold[i];
			structbox->Z[i]=a->Zold[i];

            if(i<=p->nAto)
            {
                for(j=0;j<structbox->nlen;j++)
                {
				structbox->NEIGLIST[i][j]=a->NEIGLISTold[i][j];
                }
                
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
