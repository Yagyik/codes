struct MSD
{
	double **pos0;
	double **msd;
	int counter;            
	int fileCounter;        
} MeanSquareDispl;

void Init_Msd(SIMDAT *s, DATASIM *d, MSD *ms)
{
	int t;
	//printf("init msd\n");
	ms->pos0=AllocMatR(s->Ato,3);
	//printf("init - %f %f\n",ms->pos0[0][0],ms->pos0[s->Msd_limit][0]);
	ms->msd=AllocMatR(s->Msd_limit,2);			// r^2 , r^4
	//AllocMatR(s->Msd_limit,2,ms->msd);
	//printf("init - %f %f\n",ms->msd[0][0],ms->msd[s->Msd_limit][0]);
	printf("msd limit %d\n",s->Msd_limit);
	for(t=0;t<s->Msd_limit;t++)
	{
		//printf("in this loop??\n");
		ms->msd[t][0]=0.0;		// r^2
		ms->msd[t][1]=0.0;		// r^4
		//printf("but then?? %f %f\n",ms->msd[t][1],ms->msd[t][1]);
	}
	//printf("init - %f %f\n",ms->msd[0][0],ms->msd[s->Msd_limit][0]);
	ms->counter=0;
	ms->fileCounter=0;
}

void print_Msd(SIMDAT *s, DATASIM *d, MSD *ms)
{
	printf("for real, not even here??\n");
	FILE *fp;
	int t;
	char buffer[150];

	printf("%sMsd_%d.dat\n",s->output,ms->fileCounter);
	sprintf(buffer,"%sMsd_%d.dat",s->output,ms->fileCounter);
	printf("the MSD file %s\n",buffer);
	fp=fopen(buffer,"w");
	if(fp==NULL)
	printf("bad file open\n");
	for(t=1;t<s->Msd_limit;t++)
	{
		fprintf(fp,"%le	%le	%le\n",s->DeltaT*s->Step*t,ms->msd[t][0],ms->msd[t][1]);
	}
	fclose(fp);

	sprintf(buffer,"%sAlfa2_%d.dat",s->output,ms->fileCounter);
	fp=fopen(buffer,"w");

	for(t=1;t<s->Msd_limit;t++)
	{
		fprintf(fp,"%le	%le\n",s->DeltaT*s->Step*t,3.0*ms->msd[t][1]/(5.0*sqr(ms->msd[t][0]))-1.0);
	}

	fclose(fp);

	ms->fileCounter++;
	ms->counter=0;
}

void Compute_Msd(SIMDAT *s, DATASIM *d, MSD *ms)
{
	int a,k;
	double r2;
	//printf("computing msd\n");
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
	//printf("atoms atoms %d\n",s->Ato);
	//printf("last of the atomcans - %f %f %f\n",d->r[s->Ato-1][0],d->r[s->Ato-1][1],d->r[s->Ato-1][2]);
	for(a=0;a<s->Ato;a++)
	{
		r2=0.0;
		for(k=0;k<3;k++)
		{
			r2+=sqr(d->r[a][k]-ms->pos0[a][k]);
		}
		ms->msd[ms->counter][0]+=r2;		
		ms->msd[ms->counter][1]+=r2*r2;
	}
	//printf("did we calc all atoms??\n");
	ms->msd[ms->counter][0]/=s->Ato;
	ms->msd[ms->counter][1]/=s->Ato;
	printf("done assigning??\n");
	printf("computing properly or not??? %f %f\n",ms->msd[ms->counter][0],ms->msd[ms->counter][1]);
	
	ms->counter++; 
	printf("how many files?? %d\n",s->Msd_limit);
	printf("we there??%d\n",ms->counter);
//	if(ms->counter==s->Msd_limit) 
//	{
//		printf("printing msd\n");
//		print_Msd(s,d,ms);
//		printf("printed msd\n");
//	}
}

