long int lround(double x);

double energy(int flag,struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,j;
	int neigcount,neipar=0;
	double rij=0.0,rijsq=0.0;
	double gij=0.0,delgij=0.0,delvij=0.0,delvijterm=0.0;
	double diffX,diffY,diffZ;
	double sisq=0.0,TrTisq=0.0;
	double TisqXX=0.0,TisqYY=0.0,TisqZZ=0.0;
	double U2,U3,U3Correct;
	double threebodycorrection=0.0;
	double cjicij=0.0,cjicijterm1=0.0,cjicijterm2=0.0,cjicijterm3=0.0;
	double dotproduct1=0.0,dotproduct2=0.0;
	double siminussjX=0.0,siminussjY=0.0,siminussjZ=0.0;
	double TiaddXX=0.0,TiaddXY=0.0,TiaddXZ=0.0,TiaddYY=0.0,TiaddYZ=0.0,TiaddZZ=0.0;
	double TirX=0.0,TirY=0.0,TirZ=0.0;
	double fij2X,fij2Y,fij2Z;
	double fij3X,fij3Y,fij3Z;
	double *Virial;

	Virial=calloc(p->nAto+2,sizeof(double));

	for(i=0;i<p->nAto+2;i++)
	{
		box->hi[i]=0.0;
		box->siX[i]=0.0;
		box->siY[i]=0.0;
		box->siZ[i]=0.0;
		box->TiXX[i]=0.0;
		box->TiXY[i]=0.0;
		box->TiXZ[i]=0.0;
		box->TiYY[i]=0.0;
		box->TiYZ[i]=0.0;
		box->TiZZ[i]=0.0;
		Virial[i]=0.0;
	}

	U3=U2=U3Correct=0.0;
	for(i=1;i<=p->nAto;i++)
	{
        
		neigcount=box->NEIGLIST[i][0];
        //printf(" %d neighbours for atom %d\n",neigcount,i);
		for(j=1;j<=neigcount;j++)
		{
			neipar=box->NEIGLIST[i][j];
            //printf("neig is %d\n",neipar);
			if(neipar>i)
			{
				diffX=box->X[neipar]-box->X[i];
				diffX-=box->lx*lround(diffX/box->lx);
				diffY=box->Y[neipar]-box->Y[i];
				diffY-=box->ly*lround(diffY/box->ly);
				diffZ=box->Z[neipar]-box->Z[i];
				diffZ-=box->lz*lround(diffZ/box->lz);
			   
				rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
				if(rijsq<sqr(SWcut))
				{
					rij=sqrt(rijsq);
					gij=exp(SWgamm/(rij-SWcut));

					U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
					U3Correct+=(16.0*SWlamb/9.0)*sqr(gij);

					box->hi[i]+=gij;
					box->hi[neipar]+=gij;
					box->siX[i]+=diffX*gij/rij;
					box->siY[i]+=diffY*gij/rij;
					box->siZ[i]+=diffZ*gij/rij;
					box->siX[neipar]-=diffX*gij/rij;
					box->siY[neipar]-=diffY*gij/rij;
					box->siZ[neipar]-=diffZ*gij/rij;
					box->TiXX[i]+=diffX*diffX*gij/sqr(rij);
					box->TiXY[i]+=diffX*diffY*gij/sqr(rij);
					box->TiXZ[i]+=diffX*diffZ*gij/sqr(rij);
					box->TiYY[i]+=diffY*diffY*gij/sqr(rij);
					box->TiYZ[i]+=diffY*diffZ*gij/sqr(rij);
					box->TiZZ[i]+=diffZ*diffZ*gij/sqr(rij);
					box->TiXX[neipar]+=diffX*diffX*gij/sqr(rij);
					box->TiXY[neipar]+=diffX*diffY*gij/sqr(rij);
					box->TiXZ[neipar]+=diffX*diffZ*gij/sqr(rij);
					box->TiYY[neipar]+=diffY*diffY*gij/sqr(rij);
					box->TiYZ[neipar]+=diffY*diffZ*gij/sqr(rij);
					box->TiZZ[neipar]+=diffZ*diffZ*gij/sqr(rij);
				}
			}
		}

		sisq=sqr(box->siX[i])+sqr(box->siY[i])+sqr(box->siZ[i]);
		TisqXX=sqr(box->TiXX[i])+sqr(box->TiXY[i])+sqr(box->TiXZ[i]);
		TisqYY=sqr(box->TiXY[i])+sqr(box->TiYY[i])+sqr(box->TiYZ[i]);
		TisqZZ=sqr(box->TiXZ[i])+sqr(box->TiYZ[i])+sqr(box->TiZZ[i]);
		
		U3+=SWlamb/18.0*sqr(box->hi[i])+SWlamb/3.0*sisq+SWlamb/2.0*(TisqXX+TisqYY+TisqZZ);
	}
	    //printf("done with energy\n");

	if(flag==1)
	{
		box->virial=0.0;
		for(i=1;i<=p->nAto;i++)
		{
			neigcount=box->NEIGLIST[i][0];
			for(j=1;j<=neigcount;j++)
			{
				neipar=box->NEIGLIST[i][j];
				if(neipar>i)
				{
					diffX=box->X[neipar]-box->X[i];
					diffX-=box->lx*lround(diffX/box->lx);
					diffY=box->Y[neipar]-box->Y[i];
					diffY-=box->ly*lround(diffY/box->ly);
					diffZ=box->Z[neipar]-box->Z[i];
					diffZ-=box->lz*lround(diffZ/box->lz);

					rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
					if(rijsq<sqr(SWcut))
					{
						rij=sqrt(rijsq);
				 		gij=exp(SWgamm/(rij-SWcut));
						delgij=-SWgamm*gij/sqr(rij-SWcut);

						delvijterm=-SWaCon*exp(1.0/(rij-SWcut))*((4.0*SWbCon/pow(rij,5.0))+((SWbCon/pow(rij,4.0)-1.0)/sqr(rij-SWcut)));
						fij2X=delvijterm/rij*diffX;
						fij2Y=delvijterm/rij*diffY;
						fij2Z=delvijterm/rij*diffZ;

						threebodycorrection=-(32.0/9.0)*SWlamb*gij*delgij/rij;
						siminussjX=box->siX[i]-box->siX[neipar];
						siminussjY=box->siY[i]-box->siY[neipar];
						siminussjZ=box->siZ[i]-box->siZ[neipar];
						TiaddXX=box->TiXX[i]+box->TiXX[neipar];
						TiaddXY=box->TiXY[i]+box->TiXY[neipar];
						TiaddXZ=box->TiXZ[i]+box->TiXZ[neipar];
						TiaddYY=box->TiYY[i]+box->TiYY[neipar];
						TiaddYZ=box->TiYZ[i]+box->TiYZ[neipar];
						TiaddZZ=box->TiZZ[i]+box->TiZZ[neipar];
						TirX=(TiaddXX*diffX)+(TiaddXY*diffY)+(TiaddXZ*diffZ);
						TirY=(TiaddXY*diffX)+(TiaddYY*diffY)+(TiaddYZ*diffZ);
						TirZ=(TiaddXZ*diffX)+(TiaddYZ*diffY)+(TiaddZZ*diffZ);

						dotproduct1=diffX*siminussjX+diffY*siminussjY+diffZ*siminussjZ;
						dotproduct2=diffX*TirX+diffY*TirY+diffZ*TirZ;
						cjicijterm1=(SWlamb/9.0)*delgij*(box->hi[i]+box->hi[neipar])/rij;
						cjicijterm2=(2.0*SWlamb/3.0)*(delgij-(gij/rij))*dotproduct1/sqr(rij);
						cjicijterm3=SWlamb*(delgij-2.0*(gij/rij))*dotproduct2/cube(rij);
						cjicij=cjicijterm1+cjicijterm2+cjicijterm3;

						fij3X=((cjicij+threebodycorrection)*diffX)+((2.0/3.0)*SWlamb*gij/rij*siminussjX)+(2.0*SWlamb*gij/sqr(rij)*TirX);
						fij3Y=((cjicij+threebodycorrection)*diffY)+((2.0/3.0)*SWlamb*gij/rij*siminussjY)+(2.0*SWlamb*gij/sqr(rij)*TirY);
						fij3Z=((cjicij+threebodycorrection)*diffZ)+((2.0/3.0)*SWlamb*gij/rij*siminussjZ)+(2.0*SWlamb*gij/sqr(rij)*TirZ);

						Virial[i]+=(fij2X+fij3X)*diffX+(fij2Y+fij3Y)*diffY+(fij2Z+fij3Z)*diffZ;
						Virial[neipar]+=(fij2X+fij3X)*diffX+(fij2Y+fij3Y)*diffY+(fij2Z+fij3Z)*diffZ;
					}
				}
			}
			box->virial+=Virial[i];
		}
		box->virial=-box->virial/6.0;
	}
	//printf("before free\n");
	free(Virial);
    //    printf("done with virial\n");
//     printf("energy %f %f %f %f\n",U2,U3,U3Correct,U2+U3-U3Correct);
	return(U2+U3-U3Correct);
}

double energyAnnihil(int flag,struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,j;
	int neigcount,neipar=0;
	double rij=0.0,rijsq=0.0;
	double gij=0.0,delgij=0.0,delvij=0.0,delvijterm=0.0;
	double diffX,diffY,diffZ;
	double sisq=0.0,TrTisq=0.0;
	double TisqXX=0.0,TisqYY=0.0,TisqZZ=0.0;
	double U2,U3,U3Correct;
	double threebodycorrection=0.0;
	double cjicij=0.0,cjicijterm1=0.0,cjicijterm2=0.0,cjicijterm3=0.0;
	double dotproduct1=0.0,dotproduct2=0.0;
	double siminussjX=0.0,siminussjY=0.0,siminussjZ=0.0;
	double TiaddXX=0.0,TiaddXY=0.0,TiaddXZ=0.0,TiaddYY=0.0,TiaddYZ=0.0,TiaddZZ=0.0;
	double TirX=0.0,TirY=0.0,TirZ=0.0;
	double fij2X,fij2Y,fij2Z;
	double fij3X,fij3Y,fij3Z;
	double *Virial;
	Virial=calloc(p->nAto+2,sizeof(double));
    double spam;
	for(i=0;i<p->nAto+2;i++)
	{
        spam=box->Ahi[i];
//         if(i==p->nAto+1)
//         printf("old vals for par %d %f %f %f %f %f %f %f %f %f %f\n",i,box->Ahi[i],box->AsiX[i],box->AsiY[i],box->AsiZ[i],box->ATiXX[i],box->ATiXY[i],box->ATiXZ[i],box->ATiYY[i],box->ATiYZ[i],box->ATiZZ[i]);
		box->Ahi[i]=0.0;
		box->AsiX[i]=0.0;
		box->AsiY[i]=0.0;
		box->AsiZ[i]=0.0;
		box->ATiXX[i]=0.0;
		box->ATiXY[i]=0.0;
		box->ATiXZ[i]=0.0;
		box->ATiYY[i]=0.0;
		box->ATiYZ[i]=0.0;
		box->ATiZZ[i]=0.0;
		Virial[i]=0.0;
        
	}

	U3=U2=U3Correct=0.0;
    i=p->nAto+1;
    for(j=1;j<=p->nAto;j++)
    {
        neipar=j;

        diffX=box->X[neipar]-box->X[i];
        diffX-=box->lx*lround(diffX/box->lx);
        diffY=box->Y[neipar]-box->Y[i];
        diffY-=box->ly*lround(diffY/box->ly);
        diffZ=box->Z[neipar]-box->Z[i];
        diffZ-=box->lz*lround(diffZ/box->lz);
			   
        rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
        if(rijsq<sqr(SWcut))
        {
            rij=sqrt(rijsq);
            gij=exp(SWgamm/(rij-SWcut));
		
            U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
            U3Correct+=(16.0*SWlamb/9.0)*sqr(gij);

            box->Ahi[i]+=gij;
            box->Ahi[neipar]+=gij;
            box->AsiX[i]+=diffX*gij/rij;
            box->AsiY[i]+=diffY*gij/rij;
            box->AsiZ[i]+=diffZ*gij/rij;
            box->AsiX[neipar]-=diffX*gij/rij;
            box->AsiY[neipar]-=diffY*gij/rij;
            box->AsiZ[neipar]-=diffZ*gij/rij;
            box->ATiXX[i]+=diffX*diffX*gij/sqr(rij);
            box->ATiXY[i]+=diffX*diffY*gij/sqr(rij);
            box->ATiXZ[i]+=diffX*diffZ*gij/sqr(rij);
            box->ATiYY[i]+=diffY*diffY*gij/sqr(rij);
            box->ATiYZ[i]+=diffY*diffZ*gij/sqr(rij);
            box->ATiZZ[i]+=diffZ*diffZ*gij/sqr(rij);
            box->ATiXX[neipar]+=diffX*diffX*gij/sqr(rij);
            box->ATiXY[neipar]+=diffX*diffY*gij/sqr(rij);
            box->ATiXZ[neipar]+=diffX*diffZ*gij/sqr(rij);
            box->ATiYY[neipar]+=diffY*diffY*gij/sqr(rij);
            box->ATiYZ[neipar]+=diffY*diffZ*gij/sqr(rij);
            box->ATiZZ[neipar]+=diffZ*diffZ*gij/sqr(rij);
        }
			
    }
    for(i=p->nAto+1;i<=p->nAto+1;i++) // avoid overcounting
    {
//         printf("new vals for par %d %f %f %f %f %f %f %f %f %f %f\n",i,box->Ahi[i],box->AsiX[i],box->AsiY[i],box->AsiZ[i],box->ATiXX[i],box->ATiXY[i],box->ATiXZ[i],box->ATiYY[i],box->ATiYZ[i],box->ATiZZ[i]);
        sisq=sqr(box->AsiX[i])+sqr(box->AsiY[i])+sqr(box->AsiZ[i]);
        TisqXX=sqr(box->ATiXX[i])+sqr(box->ATiXY[i])+sqr(box->ATiXZ[i]);
        TisqYY=sqr(box->ATiXY[i])+sqr(box->ATiYY[i])+sqr(box->ATiYZ[i]);
        TisqZZ=sqr(box->ATiXZ[i])+sqr(box->ATiYZ[i])+sqr(box->ATiZZ[i]);
        U3+=SWlamb/18.0*sqr(box->Ahi[i])+SWlamb/3.0*sisq+SWlamb/2.0*(TisqXX+TisqYY+TisqZZ);
    }
    

	if(flag==1)
	{
		box->Avirial=0.0;
        i=p->nAto+1;
        for(j=1;j<=p->nAto;j++)
        {
            neipar=j;
            diffX=box->X[neipar]-box->X[i];
            diffX-=box->lx*lround(diffX/box->lx);
            diffY=box->Y[neipar]-box->Y[i];
            diffY-=box->ly*lround(diffY/box->ly);
            diffZ=box->Z[neipar]-box->Z[i];
            diffZ-=box->lz*lround(diffZ/box->lz);

            rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
            if(rijsq<sqr(SWcut))
            {
                rij=sqrt(rijsq);
                gij=exp(SWgamm/(rij-SWcut));
                delgij=-SWgamm*gij/sqr(rij-SWcut);

                delvijterm=-SWaCon*exp(1.0/(rij-SWcut))*((4.0*SWbCon/pow(rij,5.0))+((SWbCon/pow(rij,4.0)-1.0)/sqr(rij-SWcut)));
                fij2X=delvijterm/rij*diffX;
                fij2Y=delvijterm/rij*diffY;
                fij2Z=delvijterm/rij*diffZ;

                threebodycorrection=-(32.0/9.0)*SWlamb*gij*delgij/rij;
                siminussjX=box->AsiX[i]-box->AsiX[neipar];
                siminussjY=box->AsiY[i]-box->AsiY[neipar];
                siminussjZ=box->AsiZ[i]-box->AsiZ[neipar];
                TiaddXX=box->ATiXX[i]+box->ATiXX[neipar];
                TiaddXY=box->ATiXY[i]+box->ATiXY[neipar];
                TiaddXZ=box->ATiXZ[i]+box->ATiXZ[neipar];
                TiaddYY=box->ATiYY[i]+box->ATiYY[neipar];
                TiaddYZ=box->ATiYZ[i]+box->ATiYZ[neipar];
                TiaddZZ=box->ATiZZ[i]+box->ATiZZ[neipar];
                TirX=(TiaddXX*diffX)+(TiaddXY*diffY)+(TiaddXZ*diffZ);
                TirY=(TiaddXY*diffX)+(TiaddYY*diffY)+(TiaddYZ*diffZ);
                TirZ=(TiaddXZ*diffX)+(TiaddYZ*diffY)+(TiaddZZ*diffZ);

                dotproduct1=diffX*siminussjX+diffY*siminussjY+diffZ*siminussjZ;
                dotproduct2=diffX*TirX+diffY*TirY+diffZ*TirZ;
                cjicijterm1=(SWlamb/9.0)*delgij*(box->Ahi[i]+box->Ahi[neipar])/rij;
                cjicijterm2=(2.0*SWlamb/3.0)*(delgij-(gij/rij))*dotproduct1/sqr(rij);
                cjicijterm3=SWlamb*(delgij-2.0*(gij/rij))*dotproduct2/cube(rij);
                cjicij=cjicijterm1+cjicijterm2+cjicijterm3;

                fij3X=((cjicij+threebodycorrection)*diffX)+((2.0/3.0)*SWlamb*gij/rij*siminussjX)+(2.0*SWlamb*gij/sqr(rij)*TirX);
                fij3Y=((cjicij+threebodycorrection)*diffY)+((2.0/3.0)*SWlamb*gij/rij*siminussjY)+(2.0*SWlamb*gij/sqr(rij)*TirY);
                fij3Z=((cjicij+threebodycorrection)*diffZ)+((2.0/3.0)*SWlamb*gij/rij*siminussjZ)+(2.0*SWlamb*gij/sqr(rij)*TirZ);

                Virial[i]+=(fij2X+fij3X)*diffX+(fij2Y+fij3Y)*diffY+(fij2Z+fij3Z)*diffZ;
                Virial[neipar]+=(fij2X+fij3X)*diffX+(fij2Y+fij3Y)*diffY+(fij2Z+fij3Z)*diffZ;
            }

			box->Avirial+=pow(box->lambda,p->n)*Virial[j];
		}
		box->Avirial+=pow(box->lambda,p->n)*Virial[p->nAto+1];
		box->Avirial=-box->Avirial/6.0;
	}
	free(Virial);
//     printf("energy annihil %f %f %f %f\n",U2,U3,U3Correct,U2+U3-U3Correct);
	return(U2+U3-U3Correct);
}

double oldparenergy(int par, struct SYSTEM *box)
{
	int neigcount,j,neig,startlist,endlist;
	double rijsq,rij;
	double diffX,diffY,diffZ;
	double gij=0.0;
	double hisq=0.0,diffhisq=0.0;
	double sisq=0.0,diffsisq=0.0;
	double TisqXX=0.0,TisqYY=0.0,TisqZZ=0.0;
	double diffTisqXX=0.0,diffTisqYY=0.0,diffTisqZZ=0.0;
	double TrTisq=0.0,diffTrTisq=0.0;
	double U3=0.0,U2=0.0,U3Correct=0.0;
	double dUhi=0.0,dUSi=0.0,dUTi=0.0;

	sisq=sqr(box->siX[par])+sqr(box->siY[par])+sqr(box->siZ[par]);
	TisqXX=sqr(box->TiXX[par])+sqr(box->TiXY[par])+sqr(box->TiXZ[par]);
	TisqYY=sqr(box->TiXY[par])+sqr(box->TiYY[par])+sqr(box->TiYZ[par]);
	TisqZZ=sqr(box->TiXZ[par])+sqr(box->TiYZ[par])+sqr(box->TiZZ[par]);

	dUhi=(double)(SWlamb/18.0)*sqr(box->hi[par]);
	dUSi=(double)(SWlamb/3.0)*sisq;
	dUTi=(double)(SWlamb/2.0)*(TisqXX+TisqYY+TisqZZ);

	neigcount=box->NEIGLIST[par][0];
	for(j=1;j<=neigcount;j++)
	{
		neig=box->NEIGLIST[par][j];

		diffX=box->X[neig]-box->X[par];
		diffX-=box->lx*lround(diffX/box->lx);
		diffY=box->Y[neig]-box->Y[par];
		diffY-=box->ly*lround(diffY/box->ly);
		diffZ=box->Z[neig]-box->Z[par];
		diffZ-=box->lz*lround(diffZ/box->lz);

		rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
		if(rijsq<sqr(SWcut))
		{
			rij=sqrt(rijsq);
			gij=exp(SWgamm/(rij-SWcut));			 
		   
			box->dhi[neig]=gij;
			box->dSiX[neig]=-diffX*gij/rij;
			box->dSiY[neig]=-diffY*gij/rij;
			box->dSiZ[neig]=-diffZ*gij/rij;
			box->dTiXX[neig]=diffX*diffX*gij/sqr(rij);
			box->dTiXY[neig]=diffX*diffY*gij/sqr(rij);
			box->dTiXZ[neig]=diffX*diffZ*gij/sqr(rij);
			box->dTiYY[neig]=diffY*diffY*gij/sqr(rij);
			box->dTiYZ[neig]=diffY*diffZ*gij/sqr(rij);
			box->dTiZZ[neig]=diffZ*diffZ*gij/sqr(rij);

			U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
			U3Correct+=(16.0*SWlamb/9.0)*sqr(gij);

			diffhisq=(2.0*box->hi[neig]-box->dhi[neig])*box->dhi[neig];
			diffsisq=(2.0*box->siX[neig]-box->dSiX[neig])*box->dSiX[neig]+(2.0*box->siY[neig]-box->dSiY[neig])*box->dSiY[neig]+(2.0*box->siZ[neig]-box->dSiZ[neig])*box->dSiZ[neig];
			diffTisqXX=(2.0*box->TiXX[neig]-box->dTiXX[neig])*box->dTiXX[neig]+(2.0*box->TiXY[neig]-box->dTiXY[neig])*box->dTiXY[neig]+(2.0*box->TiXZ[neig]-box->dTiXZ[neig])*box->dTiXZ[neig];
			diffTisqYY=(2.0*box->TiXY[neig]-box->dTiXY[neig])*box->dTiXY[neig]+(2.0*box->TiYY[neig]-box->dTiYY[neig])*box->dTiYY[neig]+(2.0*box->TiYZ[neig]-box->dTiYZ[neig])*box->dTiYZ[neig];
			diffTisqZZ=(2.0*box->TiXZ[neig]-box->dTiXZ[neig])*box->dTiXZ[neig]+(2.0*box->TiYZ[neig]-box->dTiYZ[neig])*box->dTiYZ[neig]+(2.0*box->TiZZ[neig]-box->dTiZZ[neig])*box->dTiZZ[neig];
			diffTrTisq = diffTisqXX + diffTisqYY + diffTisqZZ;

			dUhi+=(double)(SWlamb/18.0)*diffhisq;
			dUSi+=(double)(SWlamb/3.0)*diffsisq;
			dUTi+=(double)(SWlamb/2.0)*diffTrTisq;
		}
		else
		{
			box->dhi[neig]=0.0;
			box->dSiX[neig]=0.0;
			box->dSiY[neig]=0.0;
			box->dSiZ[neig]=0.0;
			box->dTiXX[neig]=0.0;
			box->dTiXY[neig]=0.0;
			box->dTiXZ[neig]=0.0;
			box->dTiYY[neig]=0.0;
			box->dTiYZ[neig]=0.0;
			box->dTiZZ[neig]=0.0;
		}
	}
	U3=dUhi+dUSi+dUTi;

	return(U2+U3-U3Correct);
}

double oldparenergyAnnihil(int par, struct SYSTEM *box,struct PARAINPUT *p)
{
	int neigcount,j,neig,startlist,endlist;
	double rijsq,rij;
	double diffX,diffY,diffZ;
	double gij=0.0;
	double hisq=0.0,diffhisq=0.0;
	double sisq=0.0,diffsisq=0.0;
	double TisqXX=0.0,TisqYY=0.0,TisqZZ=0.0;
	double diffTisqXX=0.0,diffTisqYY=0.0,diffTisqZZ=0.0;
	double TrTisq=0.0,diffTrTisq=0.0;
	double U3=0.0,U2=0.0,U3Correct=0.0;
	double dUhi=0.0,dUSi=0.0,dUTi=0.0;

	sisq=sqr(box->AsiX[par])+sqr(box->AsiY[par])+sqr(box->AsiZ[par]);
	TisqXX=sqr(box->ATiXX[par])+sqr(box->ATiXY[par])+sqr(box->ATiXZ[par]);
	TisqYY=sqr(box->ATiXY[par])+sqr(box->ATiYY[par])+sqr(box->ATiYZ[par]);
	TisqZZ=sqr(box->ATiXZ[par])+sqr(box->ATiYZ[par])+sqr(box->ATiZZ[par]);
    
	dUhi=(double)(SWlamb/18.0)*sqr(box->Ahi[par]);
	dUSi=(double)(SWlamb/3.0)*sisq;
	dUTi=(double)(SWlamb/2.0)*(TisqXX+TisqYY+TisqZZ);
    if(par==p->nAto+1)
    {
        for(j=1;j<=p->nAto;j++)
        {
            neig=j;

            diffX=box->X[neig]-box->X[par];
            diffX-=box->lx*lround(diffX/box->lx);
            diffY=box->Y[neig]-box->Y[par];
            diffY-=box->ly*lround(diffY/box->ly);
            diffZ=box->Z[neig]-box->Z[par];
            diffZ-=box->lz*lround(diffZ/box->lz);

            rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
            if(rijsq<sqr(SWcut))
            {
                rij=sqrt(rijsq);
                gij=exp(SWgamm/(rij-SWcut));			 
		   
                box->Adhi[neig]=gij;
                box->AdSiX[neig]=-diffX*gij/rij;
                box->AdSiY[neig]=-diffY*gij/rij;
                box->AdSiZ[neig]=-diffZ*gij/rij;
                box->AdTiXX[neig]=diffX*diffX*gij/sqr(rij);
                box->AdTiXY[neig]=diffX*diffY*gij/sqr(rij);
                box->AdTiXZ[neig]=diffX*diffZ*gij/sqr(rij);
                box->AdTiYY[neig]=diffY*diffY*gij/sqr(rij);
                box->AdTiYZ[neig]=diffY*diffZ*gij/sqr(rij);
                box->AdTiZZ[neig]=diffZ*diffZ*gij/sqr(rij);

                U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
                U3Correct+=(16.0*SWlamb/9.0)*sqr(gij);

                diffhisq=(2.0*box->Ahi[neig]-box->Adhi[neig])*box->Adhi[neig];
                diffsisq=(2.0*box->AsiX[neig]-box->AdSiX[neig])*box->AdSiX[neig]+(2.0*box->AsiY[neig]-box->AdSiY[neig])*box->AdSiY[neig]+(2.0*box->AsiZ[neig]-box->AdSiZ[neig])*box->AdSiZ[neig];
                diffTisqXX=(2.0*box->ATiXX[neig]-box->AdTiXX[neig])*box->AdTiXX[neig]+(2.0*box->ATiXY[neig]-box->AdTiXY[neig])*box->AdTiXY[neig]+(2.0*box->ATiXZ[neig]-box->AdTiXZ[neig])*box->AdTiXZ[neig];
                diffTisqYY=(2.0*box->ATiXY[neig]-box->AdTiXY[neig])*box->AdTiXY[neig]+(2.0*box->ATiYY[neig]-box->AdTiYY[neig])*box->AdTiYY[neig]+(2.0*box->ATiYZ[neig]-box->AdTiYZ[neig])*box->AdTiYZ[neig];
                diffTisqZZ=(2.0*box->ATiXZ[neig]-box->AdTiXZ[neig])*box->AdTiXZ[neig]+(2.0*box->ATiYZ[neig]-box->AdTiYZ[neig])*box->AdTiYZ[neig]+(2.0*box->ATiZZ[neig]-box->AdTiZZ[neig])*box->AdTiZZ[neig];
                diffTrTisq = diffTisqXX + diffTisqYY + diffTisqZZ;

                dUhi+=(double)(SWlamb/18.0)*diffhisq;
                dUSi+=(double)(SWlamb/3.0)*diffsisq;
                dUTi+=(double)(SWlamb/2.0)*diffTrTisq;
            }
            else
            {
                box->Adhi[neig]=0.0;
                box->AdSiX[neig]=0.0;
                box->AdSiY[neig]=0.0;
                box->AdSiZ[neig]=0.0;
                box->AdTiXX[neig]=0.0;
                box->AdTiXY[neig]=0.0;
                box->AdTiXZ[neig]=0.0;
                box->AdTiYY[neig]=0.0;
                box->AdTiYZ[neig]=0.0;
                box->AdTiZZ[neig]=0.0;
            }
        }
        //printf("%d par with %f %f %f\n",par,U2,U3,U3Correct);
    }
	if(par<=p->nAto)
    {
        neig=p->nAto+1;
        diffX=box->X[neig]-box->X[par];
        diffX-=box->lx*lround(diffX/box->lx);
        diffY=box->Y[neig]-box->Y[par];
        diffY-=box->ly*lround(diffY/box->ly);
        diffZ=box->Z[neig]-box->Z[par];
        diffZ-=box->lz*lround(diffZ/box->lz);

        rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
        if(rijsq<sqr(SWcut))
        {
            rij=sqrt(rijsq);
            gij=exp(SWgamm/(rij-SWcut));			 
            //printf("%d par %f rij\n",par,rij);
            box->Adhi[neig]=gij;
            box->AdSiX[neig]=-diffX*gij/rij;
            box->AdSiY[neig]=-diffY*gij/rij;
            box->AdSiZ[neig]=-diffZ*gij/rij;
            box->AdTiXX[neig]=diffX*diffX*gij/sqr(rij);
            box->AdTiXY[neig]=diffX*diffY*gij/sqr(rij);
            box->AdTiXZ[neig]=diffX*diffZ*gij/sqr(rij);
            box->AdTiYY[neig]=diffY*diffY*gij/sqr(rij);
            box->AdTiYZ[neig]=diffY*diffZ*gij/sqr(rij);
            box->AdTiZZ[neig]=diffZ*diffZ*gij/sqr(rij);

            U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
            U3Correct+=(16.0*SWlamb/9.0)*sqr(gij);

            diffhisq=(2.0*box->Ahi[neig]-box->Adhi[neig])*box->Adhi[neig];
            diffsisq=(2.0*box->AsiX[neig]-box->AdSiX[neig])*box->AdSiX[neig]+(2.0*box->AsiY[neig]-box->AdSiY[neig])*box->AdSiY[neig]+(2.0*box->AsiZ[neig]-box->AdSiZ[neig])*box->AdSiZ[neig];
            diffTisqXX=(2.0*box->ATiXX[neig]-box->AdTiXX[neig])*box->AdTiXX[neig]+(2.0*box->ATiXY[neig]-box->AdTiXY[neig])*box->AdTiXY[neig]+(2.0*box->ATiXZ[neig]-box->AdTiXZ[neig])*box->AdTiXZ[neig];
            diffTisqYY=(2.0*box->ATiXY[neig]-box->AdTiXY[neig])*box->AdTiXY[neig]+(2.0*box->ATiYY[neig]-box->AdTiYY[neig])*box->AdTiYY[neig]+(2.0*box->ATiYZ[neig]-box->AdTiYZ[neig])*box->AdTiYZ[neig];
            diffTisqZZ=(2.0*box->ATiXZ[neig]-box->AdTiXZ[neig])*box->AdTiXZ[neig]+(2.0*box->ATiYZ[neig]-box->AdTiYZ[neig])*box->AdTiYZ[neig]+(2.0*box->ATiZZ[neig]-box->AdTiZZ[neig])*box->AdTiZZ[neig];
            diffTrTisq = diffTisqXX + diffTisqYY + diffTisqZZ;

            dUhi+=(double)(SWlamb/18.0)*diffhisq;
            dUSi+=(double)(SWlamb/3.0)*diffsisq;
            dUTi+=(double)(SWlamb/2.0)*diffTrTisq;
        }
        else
        {
            box->Adhi[neig]=0.0;
            box->AdSiX[neig]=0.0;
            box->AdSiY[neig]=0.0;
            box->AdSiZ[neig]=0.0;
            box->AdTiXX[neig]=0.0;
            box->AdTiXY[neig]=0.0;
            box->AdTiXZ[neig]=0.0;
            box->AdTiYY[neig]=0.0;
            box->AdTiYZ[neig]=0.0;
            box->AdTiZZ[neig]=0.0;
        }
    }
	U3=dUhi+dUSi+dUTi;
    //if(par==p->nAto+1)
    //printf("%d par with %f %f %f\n",par,U2,U3,U3Correct);
	return(U2+U3-U3Correct);
}

double newparenergy(int par, struct SYSTEM *box)
{
	int j;
	int neig=0;
	double rij=0.0,rijsq=0.0;
	double gij=0.0;
	double diffX,diffY,diffZ;
	double U3=0.0,U2=0.0,U3Correct=0.0;
	double hisq=0.0,diffhisq=0.0;
	double sisq=0.0,diffsisq=0.0;
	double TisqXX=0.0,TisqYY=0.0,TisqZZ=0.0;
	double diffTisqXX=0.0,diffTisqYY=0.0,diffTisqZZ=0.0;
	double TrTisq=0.0,diffTrTisq=0.0;
	double dUhi=0.0,dUSi=0.0,dUTi=0.0;

	box->dhinew[par]=0.0;
	box->dSiXnew[par]=0.0;
	box->dSiYnew[par]=0.0;
	box->dSiZnew[par]=0.0;
	box->dTiXXnew[par]=0.0;
	box->dTiXYnew[par]=0.0;
	box->dTiXZnew[par]=0.0;
	box->dTiYYnew[par]=0.0;
	box->dTiYZnew[par]=0.0;
	box->dTiZZnew[par]=0.0;

	for(j=1;j<=box->NEIGPAR[0];j++)
	{
		neig=box->NEIGPAR[j];

		diffX=box->X[neig]-box->X[par];
		diffX-=box->lx*lround(diffX/box->lx);
		diffY=box->Y[neig]-box->Y[par];
		diffY-=box->ly*lround(diffY/box->ly);
		diffZ=box->Z[neig]-box->Z[par];
		diffZ-=box->lz*lround(diffZ/box->lz);

		rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
		if(rijsq<sqr(SWcut))
		{
			rij=sqrt(rijsq);
			gij=exp(SWgamm/(rij-SWcut));			

			box->dhinew[par]+=gij;
			box->dSiXnew[par]+=diffX*gij/rij;
			box->dSiYnew[par]+=diffY*gij/rij;
			box->dSiZnew[par]+=diffZ*gij/rij;
			box->dTiXXnew[par]+=diffX*diffX*gij/sqr(rij);
			box->dTiXYnew[par]+=diffX*diffY*gij/sqr(rij);
			box->dTiXZnew[par]+=diffX*diffZ*gij/sqr(rij);
			box->dTiYYnew[par]+=diffY*diffY*gij/sqr(rij);
			box->dTiYZnew[par]+=diffY*diffZ*gij/sqr(rij);
			box->dTiZZnew[par]+=diffZ*diffZ*gij/sqr(rij);
			
			box->dhinew[neig]=gij;
			box->dSiXnew[neig]=-diffX*gij/rij;
			box->dSiYnew[neig]=-diffY*gij/rij;
			box->dSiZnew[neig]=-diffZ*gij/rij;
			box->dTiXXnew[neig]=diffX*diffX*gij/sqr(rij);
			box->dTiXYnew[neig]=diffX*diffY*gij/sqr(rij);
			box->dTiXZnew[neig]=diffX*diffZ*gij/sqr(rij);
			box->dTiYYnew[neig]=diffY*diffY*gij/sqr(rij);
			box->dTiYZnew[neig]=diffY*diffZ*gij/sqr(rij);
			box->dTiZZnew[neig]=diffZ*diffZ*gij/sqr(rij);

			U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
			U3Correct+=16.0*SWlamb/9.0*sqr(gij);

			diffhisq=(2.0*(box->hi[neig]-box->dhi[neig])+box->dhinew[neig])*box->dhinew[neig];
			diffsisq=((2.0*(box->siX[neig]-box->dSiX[neig])+box->dSiXnew[neig])*box->dSiXnew[neig])+((2.0*(box->siY[neig]-box->dSiY[neig])+box->dSiYnew[neig])*box->dSiYnew[neig])+((2.0*(box->siZ[neig]-box->dSiZ[neig])+box->dSiZnew[neig])*box->dSiZnew[neig]);
			diffTisqXX=((2.0*(box->TiXX[neig]-box->dTiXX[neig])+box->dTiXXnew[neig])*box->dTiXXnew[neig])+((2.0*(box->TiXY[neig]-box->dTiXY[neig])+box->dTiXYnew[neig])*box->dTiXYnew[neig])+((2.0*(box->TiXZ[neig]-box->dTiXZ[neig])+box->dTiXZnew[neig])*box->dTiXZnew[neig]);
			diffTisqYY=((2.0*(box->TiXY[neig]-box->dTiXY[neig])+box->dTiXYnew[neig])*box->dTiXYnew[neig])+((2.0*(box->TiYY[neig]-box->dTiYY[neig])+box->dTiYYnew[neig])*box->dTiYYnew[neig])+((2.0*(box->TiYZ[neig]-box->dTiYZ[neig])+box->dTiYZnew[neig])*box->dTiYZnew[neig]);
			diffTisqZZ=((2.0*(box->TiXZ[neig]-box->dTiXZ[neig])+box->dTiXZnew[neig])*box->dTiXZnew[neig])+((2.0*(box->TiYZ[neig]-box->dTiYZ[neig])+box->dTiYZnew[neig])*box->dTiYZnew[neig])+((2.0*(box->TiZZ[neig]-box->dTiZZ[neig])+box->dTiZZnew[neig])*box->dTiZZnew[neig]);
			diffTrTisq=diffTisqXX+diffTisqYY+diffTisqZZ;

			dUhi+=(double)(SWlamb/18.0)*(diffhisq);
			dUSi+=(double)(SWlamb/3.0)*(diffsisq);
			dUTi+=(double)(SWlamb/2.0)*(diffTrTisq);
		}
		else
		{
			box->dhinew[neig]=0.0;
			box->dSiXnew[neig]=0.0;
			box->dSiYnew[neig]=0.0;
			box->dSiZnew[neig]=0.0;
			box->dTiXXnew[neig]=0.0;
			box->dTiXYnew[neig]=0.0;
			box->dTiXZnew[neig]=0.0;
			box->dTiYYnew[neig]=0.0;
			box->dTiYZnew[neig]=0.0;
			box->dTiZZnew[neig]=0.0;
		}
	}

	sisq=sqr(box->dSiXnew[par])+sqr(box->dSiYnew[par])+sqr(box->dSiZnew[par]);
	TisqXX=sqr(box->dTiXXnew[par])+sqr(box->dTiXYnew[par])+sqr(box->dTiXZnew[par]);
	TisqYY=sqr(box->dTiXYnew[par])+sqr(box->dTiYYnew[par])+sqr(box->dTiYZnew[par]);
	TisqZZ=sqr(box->dTiXZnew[par])+sqr(box->dTiYZnew[par])+sqr(box->dTiZZnew[par]);
	TrTisq=TisqXX+TisqYY+TisqZZ;

	dUhi+=(double)(SWlamb/18.0)*sqr(box->dhinew[par]);
	dUSi+=(double)(SWlamb/3.0)*sisq;
	dUTi+=(double)(SWlamb/2.0)*TrTisq;

	U3=dUhi+dUSi+dUTi;
	return(U2+U3-U3Correct);
}

double newparenergyAnnihil(int par, struct SYSTEM *box,struct PARAINPUT *p)
{
	int j;
	int neig=0;
	double rij=0.0,rijsq=0.0;
	double gij=0.0;
	double diffX,diffY,diffZ;
	double U3=0.0,U2=0.0,U3Correct=0.0;
	double hisq=0.0,diffhisq=0.0;
	double sisq=0.0,diffsisq=0.0;
	double TisqXX=0.0,TisqYY=0.0,TisqZZ=0.0;
	double diffTisqXX=0.0,diffTisqYY=0.0,diffTisqZZ=0.0;
	double TrTisq=0.0,diffTrTisq=0.0;
	double dUhi=0.0,dUSi=0.0,dUTi=0.0;

	box->Adhinew[par]=0.0;
	box->AdSiXnew[par]=0.0;
	box->AdSiYnew[par]=0.0;
	box->AdSiZnew[par]=0.0;
	box->AdTiXXnew[par]=0.0;
	box->AdTiXYnew[par]=0.0;
	box->AdTiXZnew[par]=0.0;
	box->AdTiYYnew[par]=0.0;
	box->AdTiYZnew[par]=0.0;
	box->AdTiZZnew[par]=0.0;
    if(par==p->nAto+1)
    {
        for(j=1;j<=p->nAto;j++)
        {
            neig=j;

            diffX=box->X[neig]-box->X[par];
            diffX-=box->lx*lround(diffX/box->lx);
            diffY=box->Y[neig]-box->Y[par];
            diffY-=box->ly*lround(diffY/box->ly);
            diffZ=box->Z[neig]-box->Z[par];
            diffZ-=box->lz*lround(diffZ/box->lz);

            rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
            if(rijsq<sqr(SWcut))
            {
                rij=sqrt(rijsq);
                gij=exp(SWgamm/(rij-SWcut));			

                box->Adhinew[par]+=gij;
                box->AdSiXnew[par]+=diffX*gij/rij;
                box->AdSiYnew[par]+=diffY*gij/rij;
                box->AdSiZnew[par]+=diffZ*gij/rij;
                box->AdTiXXnew[par]+=diffX*diffX*gij/sqr(rij);
                box->AdTiXYnew[par]+=diffX*diffY*gij/sqr(rij);
                box->AdTiXZnew[par]+=diffX*diffZ*gij/sqr(rij);
                box->AdTiYYnew[par]+=diffY*diffY*gij/sqr(rij);
                box->AdTiYZnew[par]+=diffY*diffZ*gij/sqr(rij);
                box->AdTiZZnew[par]+=diffZ*diffZ*gij/sqr(rij);
			
                box->Adhinew[neig]=gij;
                box->AdSiXnew[neig]=-diffX*gij/rij;
                box->AdSiYnew[neig]=-diffY*gij/rij;
                box->AdSiZnew[neig]=-diffZ*gij/rij;
                box->AdTiXXnew[neig]=diffX*diffX*gij/sqr(rij);
                box->AdTiXYnew[neig]=diffX*diffY*gij/sqr(rij);
                box->AdTiXZnew[neig]=diffX*diffZ*gij/sqr(rij);
                box->AdTiYYnew[neig]=diffY*diffY*gij/sqr(rij);
                box->AdTiYZnew[neig]=diffY*diffZ*gij/sqr(rij);
                box->AdTiZZnew[neig]=diffZ*diffZ*gij/sqr(rij);

                U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
                U3Correct+=16.0*SWlamb/9.0*sqr(gij);

                diffhisq=(2.0*(box->Ahi[neig]-box->Adhi[neig])+box->Adhinew[neig])*box->Adhinew[neig];
                
                diffsisq=((2.0*(box->AsiX[neig]-box->AdSiX[neig])+box->AdSiXnew[neig])*box->AdSiXnew[neig])+((2.0*(box->AsiY[neig]-box->AdSiY[neig])+box->AdSiYnew[neig])*box->AdSiYnew[neig])+((2.0*(box->AsiZ[neig]-box->AdSiZ[neig])+box->AdSiZnew[neig])*box->AdSiZnew[neig]);
                
                diffTisqXX=((2.0*(box->ATiXX[neig]-box->AdTiXX[neig])+box->AdTiXXnew[neig])*box->AdTiXXnew[neig])+((2.0*(box->ATiXY[neig]-box->AdTiXY[neig])+box->AdTiXYnew[neig])*box->AdTiXYnew[neig])+((2.0*(box->ATiXZ[neig]-box->AdTiXZ[neig])+box->AdTiXZnew[neig])*box->AdTiXZnew[neig]);
                
                diffTisqYY=((2.0*(box->ATiXY[neig]-box->AdTiXY[neig])+box->AdTiXYnew[neig])*box->AdTiXYnew[neig])+((2.0*(box->ATiYY[neig]-box->AdTiYY[neig])+box->AdTiYYnew[neig])*box->AdTiYYnew[neig])+((2.0*(box->ATiYZ[neig]-box->AdTiYZ[neig])+box->AdTiYZnew[neig])*box->AdTiYZnew[neig]);
                
                diffTisqZZ=((2.0*(box->ATiXZ[neig]-box->AdTiXZ[neig])+box->AdTiXZnew[neig])*box->AdTiXZnew[neig])+((2.0*(box->ATiYZ[neig]-box->AdTiYZ[neig])+box->AdTiYZnew[neig])*box->AdTiYZnew[neig])+((2.0*(box->ATiZZ[neig]-box->AdTiZZ[neig])+box->AdTiZZnew[neig])*box->AdTiZZnew[neig]);
                
                diffTrTisq=diffTisqXX+diffTisqYY+diffTisqZZ;

                dUhi+=(double)(SWlamb/18.0)*(diffhisq);
                dUSi+=(double)(SWlamb/3.0)*(diffsisq);
                dUTi+=(double)(SWlamb/2.0)*(diffTrTisq);
            }
            else
            {
                box->Adhinew[neig]=0.0;
                box->AdSiXnew[neig]=0.0;
                box->AdSiYnew[neig]=0.0;
                box->AdSiZnew[neig]=0.0;
                box->AdTiXXnew[neig]=0.0;
                box->AdTiXYnew[neig]=0.0;
                box->AdTiXZnew[neig]=0.0;
                box->AdTiYYnew[neig]=0.0;
                box->AdTiYZnew[neig]=0.0;
                box->AdTiZZnew[neig]=0.0;
            }
        }  
    }
    if(par<=p->nAto)
    {
        neig=p->nAto+1;

		diffX=box->X[neig]-box->X[par];
		diffX-=box->lx*lround(diffX/box->lx);
		diffY=box->Y[neig]-box->Y[par];
		diffY-=box->ly*lround(diffY/box->ly);
		diffZ=box->Z[neig]-box->Z[par];
		diffZ-=box->lz*lround(diffZ/box->lz);

		rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
		if(rijsq<sqr(SWcut))
		{
			rij=sqrt(rijsq);
			gij=exp(SWgamm/(rij-SWcut));			

			box->Adhinew[par]+=gij;
			box->AdSiXnew[par]+=diffX*gij/rij;
			box->AdSiYnew[par]+=diffY*gij/rij;
			box->AdSiZnew[par]+=diffZ*gij/rij;
			box->AdTiXXnew[par]+=diffX*diffX*gij/sqr(rij);
			box->AdTiXYnew[par]+=diffX*diffY*gij/sqr(rij);
			box->AdTiXZnew[par]+=diffX*diffZ*gij/sqr(rij);
			box->AdTiYYnew[par]+=diffY*diffY*gij/sqr(rij);
			box->AdTiYZnew[par]+=diffY*diffZ*gij/sqr(rij);
			box->AdTiZZnew[par]+=diffZ*diffZ*gij/sqr(rij);
			
			box->Adhinew[neig]=gij;
			box->AdSiXnew[neig]=-diffX*gij/rij;
			box->AdSiYnew[neig]=-diffY*gij/rij;
			box->AdSiZnew[neig]=-diffZ*gij/rij;
			box->AdTiXXnew[neig]=diffX*diffX*gij/sqr(rij);
			box->AdTiXYnew[neig]=diffX*diffY*gij/sqr(rij);
			box->AdTiXZnew[neig]=diffX*diffZ*gij/sqr(rij);
			box->AdTiYYnew[neig]=diffY*diffY*gij/sqr(rij);
			box->AdTiYZnew[neig]=diffY*diffZ*gij/sqr(rij);
			box->AdTiZZnew[neig]=diffZ*diffZ*gij/sqr(rij);

			U2+=(SWaCon*((SWbCon*pow(rij,-4.0))-1.0)*exp(1.0/(rij-SWcut)));
			U3Correct+=16.0*SWlamb/9.0*sqr(gij);

			diffhisq=(2.0*(box->Ahi[neig]-box->Adhi[neig])+box->Adhinew[neig])*box->Adhinew[neig];
            diffsisq=((2.0*(box->AsiX[neig]-box->AdSiX[neig])+box->AdSiXnew[neig])*box->AdSiXnew[neig])+((2.0*(box->AsiY[neig]-box->AdSiY[neig])+box->AdSiYnew[neig])*box->AdSiYnew[neig])+((2.0*(box->AsiZ[neig]-box->AdSiZ[neig])+box->AdSiZnew[neig])*box->AdSiZnew[neig]);
            
            diffTisqXX=((2.0*(box->ATiXX[neig]-box->AdTiXX[neig])+box->AdTiXXnew[neig])*box->AdTiXXnew[neig])+((2.0*(box->ATiXY[neig]-box->AdTiXY[neig])+box->AdTiXYnew[neig])*box->AdTiXYnew[neig])+((2.0*(box->ATiXZ[neig]-box->AdTiXZ[neig])+box->AdTiXZnew[neig])*box->AdTiXZnew[neig]);
            
            diffTisqYY=((2.0*(box->ATiXY[neig]-box->AdTiXY[neig])+box->AdTiXYnew[neig])*box->AdTiXYnew[neig])+((2.0*(box->ATiYY[neig]-box->AdTiYY[neig])+box->AdTiYYnew[neig])*box->AdTiYYnew[neig])+((2.0*(box->ATiYZ[neig]-box->AdTiYZ[neig])+box->AdTiYZnew[neig])*box->AdTiYZnew[neig]);
            
            diffTisqZZ=((2.0*(box->ATiXZ[neig]-box->AdTiXZ[neig])+box->AdTiXZnew[neig])*box->AdTiXZnew[neig])+((2.0*(box->ATiYZ[neig]-box->AdTiYZ[neig])+box->AdTiYZnew[neig])*box->AdTiYZnew[neig])+((2.0*(box->ATiZZ[neig]-box->AdTiZZ[neig])+box->AdTiZZnew[neig])*box->AdTiZZnew[neig]);
            
            diffTrTisq=diffTisqXX+diffTisqYY+diffTisqZZ;

			dUhi+=(double)(SWlamb/18.0)*(diffhisq);
			dUSi+=(double)(SWlamb/3.0)*(diffsisq);
			dUTi+=(double)(SWlamb/2.0)*(diffTrTisq);
		}
		else
		{
			box->Adhinew[neig]=0.0;
			box->AdSiXnew[neig]=0.0;
			box->AdSiYnew[neig]=0.0;
			box->AdSiZnew[neig]=0.0;
			box->AdTiXXnew[neig]=0.0;
			box->AdTiXYnew[neig]=0.0;
			box->AdTiXZnew[neig]=0.0;
			box->AdTiYYnew[neig]=0.0;
			box->AdTiYZnew[neig]=0.0;
			box->AdTiZZnew[neig]=0.0;
		}
    }
	

	sisq=sqr(box->AdSiXnew[par])+sqr(box->AdSiYnew[par])+sqr(box->AdSiZnew[par]);
	TisqXX=sqr(box->AdTiXXnew[par])+sqr(box->AdTiXYnew[par])+sqr(box->AdTiXZnew[par]);
	TisqYY=sqr(box->AdTiXYnew[par])+sqr(box->AdTiYYnew[par])+sqr(box->AdTiYZnew[par]);
	TisqZZ=sqr(box->AdTiXZnew[par])+sqr(box->AdTiYZnew[par])+sqr(box->AdTiZZnew[par]);
	TrTisq=TisqXX+TisqYY+TisqZZ;

	dUhi+=(double)(SWlamb/18.0)*sqr(box->Adhinew[par]);
	dUSi+=(double)(SWlamb/3.0)*sisq;
	dUTi+=(double)(SWlamb/2.0)*TrTisq;

	U3=dUhi+dUSi+dUTi;
	return(U2+U3-U3Correct);
}

void updatehiSiTi( int par, double Xold, double Yold, double Zold, struct SYSTEM *box )
{
	int neigcount,nei,j,k;
	int startlist,endlist;
	double diffX,diffY,diffZ;
	double rij,rijsq,gij=0.0;
	double hitmp,siXtmp,siYtmp,siZtmp;
	double TiXXtmp,TiXYtmp,TiXZtmp,TiYYtmp,TiYZtmp,TiZZtmp;

	neigcount=box->NEIGLIST[par][0];
	for(j=1;j<=neigcount;j++)
	{
		nei=box->NEIGLIST[par][j];

		diffX=Xold-box->X[nei];
		diffX-=box->lx*lround(diffX/box->lx);
		diffY=Yold-box->Y[nei];
		diffY-=box->ly*lround(diffY/box->ly);
		diffZ=Zold-box->Z[nei];
		diffZ-=box->lz*lround(diffZ/box->lz);
		
		rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
		if(rijsq<sqr(SWcut))
		{
			box->hi[nei]-=box->dhi[nei];
			box->siX[nei]-=box->dSiX[nei];
			box->siY[nei]-=box->dSiY[nei];
			box->siZ[nei]-=box->dSiZ[nei];
			box->TiXX[nei]-=box->dTiXX[nei];
			box->TiXY[nei]-=box->dTiXY[nei];
			box->TiXZ[nei]-=box->dTiXZ[nei];
			box->TiYY[nei]-=box->dTiYY[nei];
			box->TiYZ[nei]-=box->dTiYZ[nei];
			box->TiZZ[nei]-=box->dTiZZ[nei];
		}
	}

	for(j=1;j<=box->NEIGPAR[0];j++)
	{
		nei=box->NEIGPAR[j];

		diffX=box->X[par]-box->X[nei];
		diffX-=box->lx*lround(diffX/box->lx);
		diffY=box->Y[par]-box->Y[nei];
		diffY-=box->ly*lround(diffY/box->ly);
		diffZ=box->Z[par]-box->Z[nei];
		diffZ-=box->lz*lround(diffZ/box->lz);
	
		rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
		if(rijsq<sqr(SWcut))
		{
			box->hi[nei]+=box->dhinew[nei];
			box->siX[nei]+=box->dSiXnew[nei];
			box->siY[nei]+=box->dSiYnew[nei];
			box->siZ[nei]+=box->dSiZnew[nei];
			box->TiXX[nei]+=box->dTiXXnew[nei];
			box->TiXY[nei]+=box->dTiXYnew[nei];
			box->TiXZ[nei]+=box->dTiXZnew[nei];
			box->TiYY[nei]+=box->dTiYYnew[nei];
			box->TiYZ[nei]+=box->dTiYZnew[nei];
			box->TiZZ[nei]+=box->dTiZZnew[nei];
		}
	}

	box->hi[par]=box->dhinew[par];
	box->siX[par]=box->dSiXnew[par];
	box->siY[par]=box->dSiYnew[par];
	box->siZ[par]=box->dSiZnew[par];
	box->TiXX[par]=box->dTiXXnew[par];
	box->TiXY[par]=box->dTiXYnew[par];
	box->TiXZ[par]=box->dTiXZnew[par];
	box->TiYY[par]=box->dTiYYnew[par];
	box->TiYZ[par]=box->dTiYZnew[par];
	box->TiZZ[par]=box->dTiZZnew[par];
}

void updatehiSiTiAnnihil( int par, double Xold, double Yold, double Zold, struct SYSTEM *box, struct PARAINPUT *p)
{
	int neigcount,nei,j,k;
	int startlist,endlist;
	double diffX,diffY,diffZ;
	double rij,rijsq,gij=0.0;
	double hitmp,siXtmp,siYtmp,siZtmp;
	double TiXXtmp,TiXYtmp,TiXZtmp,TiYYtmp,TiYZtmp,TiZZtmp;

    if(par==p->nAto+1)
    {
        for(j=1;j<=p->nAto;j++)
        {
            nei=j;

            diffX=Xold-box->X[nei];
            diffX-=box->lx*lround(diffX/box->lx);
            diffY=Yold-box->Y[nei];
            diffY-=box->ly*lround(diffY/box->ly);
            diffZ=Zold-box->Z[nei];
            diffZ-=box->lz*lround(diffZ/box->lz);
            
            rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
            if(rijsq<sqr(SWcut))
            {
                box->Ahi[nei]-=box->Adhi[nei];
                box->AsiX[nei]-=box->AdSiX[nei];
                box->AsiY[nei]-=box->AdSiY[nei];
                box->AsiZ[nei]-=box->AdSiZ[nei];
                box->ATiXX[nei]-=box->AdTiXX[nei];
                box->ATiXY[nei]-=box->AdTiXY[nei];
                box->ATiXZ[nei]-=box->AdTiXZ[nei];
                box->ATiYY[nei]-=box->AdTiYY[nei];
                box->ATiYZ[nei]-=box->AdTiYZ[nei];
                box->ATiZZ[nei]-=box->AdTiZZ[nei];
            }
        }

        for(j=1;j<=p->nAto;j++)
        {
            nei=j;

            diffX=box->X[par]-box->X[nei];
            diffX-=box->lx*lround(diffX/box->lx);
            diffY=box->Y[par]-box->Y[nei];
            diffY-=box->ly*lround(diffY/box->ly);
            diffZ=box->Z[par]-box->Z[nei];
            diffZ-=box->lz*lround(diffZ/box->lz);
        
            rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
            if(rijsq<sqr(SWcut))
            {
                box->Ahi[nei]+=box->Adhinew[nei];
                box->AsiX[nei]+=box->AdSiXnew[nei];
                box->AsiY[nei]+=box->AdSiYnew[nei];
                box->AsiZ[nei]+=box->AdSiZnew[nei];
                box->ATiXX[nei]+=box->AdTiXXnew[nei];
                box->ATiXY[nei]+=box->AdTiXYnew[nei];
                box->ATiXZ[nei]+=box->AdTiXZnew[nei];
                box->ATiYY[nei]+=box->AdTiYYnew[nei];
                box->ATiYZ[nei]+=box->AdTiYZnew[nei];
                box->ATiZZ[nei]+=box->AdTiZZnew[nei];
            }
        }

        box->Ahi[par]=box->Adhinew[par];
        box->AsiX[par]=box->AdSiXnew[par];
        box->AsiY[par]=box->AdSiYnew[par];
        box->AsiZ[par]=box->AdSiZnew[par];
        box->ATiXX[par]=box->AdTiXXnew[par];
        box->ATiXY[par]=box->AdTiXYnew[par];
        box->ATiXZ[par]=box->AdTiXZnew[par];
        box->ATiYY[par]=box->AdTiYYnew[par];
        box->ATiYZ[par]=box->AdTiYZnew[par];
        box->ATiZZ[par]=box->AdTiZZnew[par];
    }
    if(par<=p->nAto)
    {
        nei=p->nAto+1;

        diffX=Xold-box->X[nei];
        diffX-=box->lx*lround(diffX/box->lx);
        diffY=Yold-box->Y[nei];
        diffY-=box->ly*lround(diffY/box->ly);
        diffZ=Zold-box->Z[nei];
        diffZ-=box->lz*lround(diffZ/box->lz);
        
        rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
        if(rijsq<sqr(SWcut))
        {
            box->Ahi[nei]-=box->Adhi[nei];
            box->AsiX[nei]-=box->AdSiX[nei];
            box->AsiY[nei]-=box->AdSiY[nei];
            box->AsiZ[nei]-=box->AdSiZ[nei];
            box->ATiXX[nei]-=box->AdTiXX[nei];
            box->ATiXY[nei]-=box->AdTiXY[nei];
            box->ATiXZ[nei]-=box->AdTiXZ[nei];
            box->ATiYY[nei]-=box->AdTiYY[nei];
            box->ATiYZ[nei]-=box->AdTiYZ[nei];
            box->ATiZZ[nei]-=box->AdTiZZ[nei];
            
            //
            
            box->Ahi[nei]+=box->Adhinew[nei];
            box->AsiX[nei]+=box->AdSiXnew[nei];
            box->AsiY[nei]+=box->AdSiYnew[nei];
            box->AsiZ[nei]+=box->AdSiZnew[nei];
            box->ATiXX[nei]+=box->AdTiXXnew[nei];
            box->ATiXY[nei]+=box->AdTiXYnew[nei];
            box->ATiXZ[nei]+=box->AdTiXZnew[nei];
            box->ATiYY[nei]+=box->AdTiYYnew[nei];
            box->ATiYZ[nei]+=box->AdTiYZnew[nei];
            box->ATiZZ[nei]+=box->AdTiZZnew[nei];
        }
        box->Ahi[par]=box->Adhinew[par];
        box->AsiX[par]=box->AdSiXnew[par];
        box->AsiY[par]=box->AdSiYnew[par];
        box->AsiZ[par]=box->AdSiZnew[par];
        box->ATiXX[par]=box->AdTiXXnew[par];
        box->ATiXY[par]=box->AdTiXYnew[par];
        box->ATiXZ[par]=box->AdTiXZnew[par];
        box->ATiYY[par]=box->AdTiYYnew[par];
        box->ATiYZ[par]=box->AdTiYZnew[par];
        box->ATiZZ[par]=box->AdTiZZnew[par]; 
    }
	
}

