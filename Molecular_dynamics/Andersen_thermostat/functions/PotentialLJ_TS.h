long int lround(double x);

double energy(int flag,struct SYSTEM *box,struct PARAINPUT *p,int flag)
{
	int i,j;
	int neigcount,neipar=0;
	double rij=0.0,rijsq=0.0,invrij=0.0;
	double diffX,diffY,diffZ;
	double *Virial;
    double en=0.0,ff=0.0;
	Virial=calloc(p->nAto,sizeof(double));
    
    
    double TS_fact=4.0*p->eps*(pow(p->sigma/p->rc,12.0)-pow(p->sigma/p->rc,6.0));

    
    //     printf("%f %f\n",p->eps,TS_fact);
    
    
    
    
    
    
    
    for(i=0;i<p->nAto;i++)
    {
        neigcount=box->NEIGLIST[i][0];
//         if(i==2 || i==213)
//         printf("full %d %d\n",i,neigcount);
		for(j=1;j<=neigcount;j++)
		{
            
			neipar=box->NEIGLIST[i][j];
            //printf("full %d %d %d\n",i,neipar,neigcount);
			if(neipar>i)
			{
				diffX=box->X[neipar]-box->X[i];
				diffX-=box->lx*lround(diffX/box->lx); // this line applies PBC
                
				diffY=box->Y[neipar]-box->Y[i];
				diffY-=box->ly*lround(diffY/box->ly);
                
				diffZ=box->Z[neipar]-box->Z[i];
				diffZ-=box->lz*lround(diffZ/box->lz);
			   
				rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
				if(rijsq<sqr(p->rc))
				{
					rij=sqrt(rijsq);
                    //printf("%f %f %f\n",rij,pow(p->sigma/rij,12.0),pow(p->sigma/rij,6.0));
                    en+=4.0*p->eps*(pow(p->sigma/rij,12.0) - pow(p->sigma/rij,6.0)) - TS_fact;
                    //printf("%f %f\n",rij,en);
                    
                }
            }

        }
       // printf("%d %f\n",i,en);
    }

	if(flag==1)
	{
        // virial[i] = sum_j fij.rij
        // box.virial= -1/6 (sum_i virial[i])
        // 
        box->virial=0.0;
        
        for(i=0;i<p->nAto;i++)
        {
            neigcount=box->NEIGLIST[i][0];
            for(j=1;j<=neigcount;j++)
            {
                ff=0.0;
                neipar=box->NEIGLIST[i][j];
                
                diffX=box->X[neipar]-box->X[i];
                diffX-=box->lx*lround(diffX/box->lx);
                diffY=box->Y[neipar]-box->Y[i];
                diffY-=box->ly*lround(diffY/box->ly);
                diffZ=box->Z[neipar]-box->Z[i];
                diffZ-=box->lz*lround(diffZ/box->lz);
            
                rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
                if(rijsq<sqr(p->rc))
                {
                    rij=sqrt(rijsq);
                    invrij=1.0/rij;
                    ff = 48.0*p->eps*(pow(p->sigma/rij,12)*invrij) - 24.0*p->eps*(pow(p->sigma/rij,6)*invrij);
                }

                box->fx[i]+=ff*diffX*invrij;
                box->fy[i]+=ff*diffY*invrij;
                box->fz[i]+=ff*diffZ*invrij;
                box->fx[neipar]-=ff*diffX*invrij;
                box->fy[neipar]-=ff*diffY*invrij;
                box->fz[neipar]-=ff*diffZ*invrij;
                
                
                Virial[i]+=ff*rij;
            }
            //Virial[i]+=box->fx[i]*box->X[i]+box->fy[i]*box->Y[i]+box->fz[i]*box->Z[i];
        }
        for(i=0;i<p->nAto;i++)
             box->virial+=Virial[i];
        //box->virial=box->virial/30.0;
        box->virial=-box->virial/3.0;
		
	}
// 	printf("en tot %f\n",en);
    free(Virial);
	return en;
}

double oldparenergy(int par, struct SYSTEM *box,struct PARAINPUT *p)
{
	int neigcount,i,j,neig,startlist,endlist;
	double rijsq,rij;
    double diffX,diffY,diffZ;
    double en=0.0;
    double TS_fact=4.0*p->eps*(pow(p->sigma/p->rc,12.0)-pow(p->sigma/p->rc,6.0));
    
    
    
    neigcount=box->NEIGLIST[par][0];
//     if(par==2 || par==213)
//     printf("\npar old %d %d \n",par,neigcount);
    for(i=1;i<=neigcount;i++)
    {
        neig=box->NEIGLIST[par][i];
        //printf("par neig count %d %d %d\n",par,neig,neigcount);
        diffX=box->X[neig]-box->X[par];
        diffX-=box->lx*lround(diffX/box->lx);
        diffY=box->Y[neig]-box->Y[par];
        diffY-=box->ly*lround(diffY/box->ly);
        diffZ=box->Z[neig]-box->Z[par];
        diffZ-=box->lz*lround(diffZ/box->lz);
        
        rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
        if(rijsq<sqr(p->rc))
        {
            rij=sqrt(rijsq);
            
            en+=4.0*p->eps*(pow(p->sigma/rij,12.0) - pow(p->sigma/rij,6.0)) - TS_fact;
            //printf("rij en %f %f\n",rij,en);
        }


    }
    
    return en;
    
}

double newparenergy(int par, struct SYSTEM *box, struct PARAINPUT *p)
{
	int i,j;
	int neig=0;
	double rij=0.0,rijsq=0.0;
    double diffX,diffY,diffZ;
	double TS_fact=4.0*p->eps*(pow(p->sigma/p->rc,12.0)-pow(p->sigma/p->rc,6.0));
    double en=0.0;
//     if(par==2 || par==213)
//     printf("par new %d %d\n",par,box->NEIGPAR[0]);
    for(i=1;i<=box->NEIGPAR[0];i++)
    {
        //printf("new %d %d %d\n",par,i,box->NEIGPAR[0]);
        neig=box->NEIGPAR[i];
        diffX=box->X[neig]-box->X[par];
        diffX-=box->lx*lround(diffX/box->lx);
        diffY=box->Y[neig]-box->Y[par];
        diffY-=box->ly*lround(diffY/box->ly);
        diffZ=box->Z[neig]-box->Z[par];
        diffZ-=box->lz*lround(diffZ/box->lz);
        
        rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
        if(rijsq<sqr(p->rc))
        {
            rij=sqrt(rijsq);
            en+=4.0*p->eps*(pow(p->sigma/rij,12.0) - pow(p->sigma/rij,6.0)) - TS_fact;
            //printf("new %f %f\n",rij,en);
        }
    }
    return en; 
}


