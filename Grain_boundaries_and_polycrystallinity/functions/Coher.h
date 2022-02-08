struct COHER
{
	double ***cohe;		// Final cohe
	int *numK;			// Number of Wavevectors for each time
	int nWaveVec;		// Number of Wavevectors
	int fileCounter;
	int counter; 
} coher;

void Init_Coher(SIMDAT *s, DATASIM *d, COHER *coh)
{
	int t,k;

	coh->cohe=AllocTensR(s->Coher_limit,MaxKappa,2);
	coh->numK=AllocVecI(s->Coher_limit);
	for(t=0;t<s->Coher_limit;t++)
	{
	for(k=0;k<MaxKappa;k++)
	{
		coh->cohe[t][k][0]=0.0;		// Sum of cos
		coh->cohe[t][k][1]=0.0;		// Sum of sin
	}
	}

	coh->counter=0;
	coh->fileCounter=0;
}

void print_Coher(SIMDAT *s, DATASIM *d, COHER *coh)
{	
	FILE *fp;
	int t,k;
	char buffer[500];
	double trash,trash0;

	sprintf(buffer,"%sCoher%g_%d.dat",s->output,s->Coher_value,coh->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=0;t<s->Coher_limit;t++)
	{
		trash=0.0;
		for(k=0;k<coh->numK[t];k++)
		{
			trash+=coh->cohe[t][k][0]*coh->cohe[0][k][0]+coh->cohe[t][k][1]*coh->cohe[0][k][1];
		}
		trash/=(coh->numK[t]*s->Ato);

		if(t==0)
		{
			trash0=trash;
		}
		else
		{
			fprintf(fp,"%le	%le	%le\n",s->DeltaT*s->Step*t,trash/trash0,trash);
		}
	}
	fclose(fp);

	coh->fileCounter++;
	coh->counter=0;
	
}

void Compute_Coher(SIMDAT *s, DATASIM *d, COHER *coh)
{
	int a,x;
	int qx,qy,qz,qthr;
	double q;

	coh->nWaveVec=0;												// Reset number of K
	qthr=rint(d->lBox[0]*(s->Coher_value+0.05)/(2.0*M_PI));			// Tolerance of 0.05

	for(qx=-qthr;qx<=qthr;qx++)
	{
	for(qy=-qthr;qy<=qthr;qy++)
	{
	for(qz=-qthr;qz<=qthr;qz++)
	{
		q=2.0*M_PI/d->lBox[0]*sqrt(qx*qx+qy*qy+qz*qz); 				// Modulus of K
		if(q<(s->Coher_value+0.05) && q>(s->Coher_value-0.05))		
		{
			if(coh->nWaveVec==MaxKappa)
			{
				printf("Error in Coher: MaxKappa must be greater than %d\n",MaxKappa);
				exit(0);
			}
			
			for(a=0;a<s->Ato;a++)
			{
				coh->cohe[coh->counter][coh->nWaveVec][0]+=cos(qx*2.0*M_PI/d->lBox[0]*d->r[a][0]+qy*2.0*M_PI/d->lBox[1]*d->r[a][1]+qz*2.0*M_PI/d->lBox[2]*d->r[a][2]);
				coh->cohe[coh->counter][coh->nWaveVec][1]+=sin(qx*2.0*M_PI/d->lBox[0]*d->r[a][0]+qy*2.0*M_PI/d->lBox[1]*d->r[a][1]+qz*2.0*M_PI/d->lBox[2]*d->r[a][2]);  
			}
			coh->nWaveVec++;										// Counter number of K
		}
	}
	}
	}

	coh->numK[coh->counter]=coh->nWaveVec;
	coh->counter++; 
  
	if(coh->counter==s->Coher_limit) 
	{
		print_Coher(s,d,coh);
	}
}

