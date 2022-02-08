struct MSD
{
	double **pos0;
	double *msd;
	int counter;            
	int fileCounter;        
} MeanSquareDispl;
// refer the isf file for this.
// the data structure msd is 1d here, different from seba version
void Init_Msd(SIMDAT *s, DATASIM *d, MSD *ms)
{
	int t;

	ms->pos0=AllocMatR(s->Ato,3);
	ms->msd=AllocVecR(s->Msd_limit);			// r^2 , r^4
	for(t=0;t<s->Msd_limit;t++)
	{
		ms->msd[t]=0.0;		// r^2
	}

	ms->counter=0;
	ms->fileCounter=0;
}
// above, ms->msd is allocated as a vector

void Compute_Msd(SIMDAT *s, DATASIM *d, MSD *ms)
{
	int a,k;
	double r2;

	if(ms->counter==0)
	{	
		for(a=0;a<s->Ato;a++)
		{
		for(k=0;k<3;k++)
		{
			ms->pos0[a][k]=d->r[a][k];
		}
		}
	}
	// initial is stored above
	for(a=0;a<s->Ato;a++)
	{	
		// for each atom
		// reset r2 which is the squared distance
		r2=0.0;
		for(k=0;k<3;k++)
		{
			//for each direction
			r2+=sqr(d->r[a][k]-ms->pos0[a][k]);
			// ath atom, kth coord and time counter and time 0
			// add the square of the difference in coordinate
			// msd = sig (x(t) - x(0))^2 + y.... 
		}
		ms->msd[ms->counter]+=r2;		
		//msd of counter has an added component because of ath atom
	}
	
	ms->msd[ms->counter]/=s->Ato;
	// divide by number of atoms

	
	ms->counter++; 
	// next time for when the function is called next
	if(ms->counter==s->Msd_limit) 
	{
		print_Msd(s,d,ms);
	}
	// if above for when we're done and we want to print
}


void print_Msd(SIMDAT *s, DATASIM *d, MSD *ms)
{
	FILE *fp;
	int t;
	char buffer[150];

	sprintf(buffer,"%sMsd_%d.dat",s->output,ms->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=1;t<s->Msd_limit;t++)
	{
		fprintf(fp,"%le	%le	\n",s->DeltaT*s->Step*t,ms->msd[t]);
	}
	fclose(fp);
	// at this stage, Seba's code has another quantity that is being calculated that we need to look up.
	ms->fileCounter++;
	ms->counter=0;
}


