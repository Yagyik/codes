double gaussVar(double mu,double sig,struct APPOGGIO *a)
{
    double gauss_r,v1,v2,gauss_var;
    gauss_r=2.0;
    while(gauss_r>=1.0)
    {
        v1=2.0*genrand_real1(&(*a))-1.0;
        v2=2.0*genrand_real1(&(*a))-1.0;
//         printf("v1 v2 %f %f\n",v1,v2);
        gauss_r = v1*v1 + v2*v2;
    }
//     printf("gr log gr %f %f\n",gauss_r,log(gauss_r));
    gauss_var=v1*sqrt(-2.0*log(gauss_r)/gauss_r);
    gauss_var=mu+sig*gauss_var;
    
    return gauss_var;
}


void creaAtom(struct PARAINPUT *p,struct SYSTEM *structbox,struct APPOGGIO *a,int flag)
{
	int i,j,k;
	int atomdone,atomnotdone;
	double dummyX,dummyY,dummyZ;
	double rX,rY,rZ,r;
	
	structbox->lx=structbox->ly=structbox->lz=pow((p->nAto/structbox->rho),(1.0/3.0));
    // ad hoc
    
    double useLatt=structbox->lx/4.0;
    for(i=0;i<4;i++)
    {
        structbox->xunit[i]*=useLatt/alat;
        structbox->yunit[i]*=useLatt/alat;
        structbox->zunit[i]*=useLatt/alat;
    }
    printf("created box with lx %f\n",structbox->lx);
	if(flag==0)
    {
        for(i=0;i<p->nAto;i++)
        {
            atomdone=0;
            while(atomdone==0)
            {
                dummyX=structbox->lx*genrand_real1(&(*a));
                dummyY=structbox->ly*genrand_real1(&(*a));
                dummyZ=structbox->lz*genrand_real1(&(*a));

                atomnotdone=0;
                for(j=0;j<i;j++)
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
//                 printf("%d %f %f %f\n",i,dummyX,dummyY,dummyZ);
            }
        }  
    }
    
    // create crystal
    int spamcount,m;
    if(flag == 1)
	{	
        spamcount = 0;
		if(p->nAto != 256) 
		{
		printf("making crystal from scratch needs 256 atoms");
		exit(1);
		}
		for(i=0;i<4;i++)
		{
            for(j=0;j<4;j++)
            {
                for(k=0;k<4;k++)
                {
                    for(m=0;m<4;m++)
                    {
                    structbox->X[spamcount] = i*useLatt + structbox->xunit[m];
                    structbox->Y[spamcount] = j*useLatt + structbox->yunit[m]; 
                    structbox->Z[spamcount] = k*useLatt + structbox->zunit[m];   
//                     printf("%f %f %f\n",structbox->X[spamcount],structbox->Y[spamcount],structbox->Z[spamcount]);
                    spamcount++;
                    }
                }
            }
        }
	}
    
	// Put origin in (0,0,0)
	for(i=0;i<p->nAto;i++)
	{
		structbox->X[i]-=structbox->lx*0.5;
		structbox->Y[i]-=structbox->ly*0.5;
		structbox->Z[i]-=structbox->lz*0.5;
	}
}

int init_vel(struct PARAINPUT *p,struct SYSTEM *structbox,struct APPOGGIO *a)
{
    int i;
    double vxcom=0.0, vycom=0.0,vzcom=0.0;
    double sig=sqrt(structbox->tem);
    for(i=0;i<p->nAto;i++)
    {
        structbox->vx[i]=gaussVar(0,sig,&(*a));
        structbox->vy[i]=gaussVar(0,sig,&(*a));
        structbox->vz[i]=gaussVar(0,sig,&(*a));
        //printf("%f %f %f\n",structbox->vx[i],structbox->vy[i],structbox->vz[i]);
        vxcom+=structbox->vx[i];
        vycom+=structbox->vy[i];
        vzcom+=structbox->vz[i];
    }
    
    for(i=0;i<p->nAto;i++)
    {
        structbox->vx[i]-=vxcom/p->nAto;
        structbox->vy[i]-=vycom/p->nAto;
        structbox->vz[i]-=vzcom/p->nAto;
//         printf("velocities %f %f %f\n",structbox->vx[i],structbox->vy[i],structbox->vz[i]);
    }
}





/*
Real function gauss_var(l0,sig)
    Implicit none
    real(kind=8),intent(in):: l0,sig
    Real(kind=8)::gauss_r,v1,v2
    Gauss_r = 2.0d0

    Do while(gauss_r >=1.0d0)
       V1 = 2*rand()-1
       V2 = 2*rand()-1
       Gauss_r = v1**2 + v2**2
    Enddo
    Gauss_var = v1*sqrt(-2.0d0*log(gauss_r)/gauss_r)
    Gauss_var = l0+sig*gauss_var

  End function gauss_var*/



