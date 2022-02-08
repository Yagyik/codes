struct ISF
{
	double *isf;		// Final ISF
	double **pos0;		// Initial position
	int nWaveVec;		// Number of Wavevectors
	int fileCounter;
	int counter; 
} IsfK;
// isf stores info v time
// pos0 is initial position
// nWaveVec is for normalisation
//fileCounter unknown

// we need a subroutine to initialise all the values of type structure

void Init_Isf(SIMDAT *s, DATASIM *d, ISF *isf)
{
	int t;

	isf->pos0=AllocMatR(s->Ato,3); 
	isf->isf=AllocVecR(s->Isf_limit);
	for(t=0;t<s->Isf_limit;t++)
	{
		isf->isf[t]=0.0;
	}

	isf->counter=0;
	isf->fileCounter=0;
}


// computing subroutine
void Compute_Isf(SIMDAT *s, DATASIM *d, ISF *isf)
{
	int a,x;
	int qx,qy,qz,qthr;
	double q;
	// if the first entry of isf, we store the pos0 as whatever the first file value is
	// the structures s contains the time step, the structure d the positions, the structure isf, the place to keep the output
	if(isf->counter==0)
	{
		for(a=0;a<s->Ato;a++)
		{
		for(x=0;x<3;x++)
		{
			isf->pos0[a][x]=d->r[a][x];
		}
		}
	}
	
	isf->nWaveVec=0;												// Reset number of K
	qthr=rint(d->lBox[0]*(s->IsfK_value+0.05)/(2.0*M_PI));			// Tolerance of 0.05
	// define bounds above - tolerance picked arbitrarily
	for(qx=-qthr;qx<=qthr;qx++)
	{
	for(qy=-qthr;qy<=qthr;qy++)
	{
	for(qz=-qthr;qz<=qthr;qz++)
	{	
		// nested loop above for all q vectors
		// we pick the combinations such that they are within some tolerance of isfk_value
		// tolerance value arbitrary
		q=2.0*M_PI/d->lBox[0]*sqrt(qx*qx+qy*qy+qz*qz); 				// Modulus of K
		if(q<(s->IsfK_value+0.05) && q>(s->IsfK_value-0.05))		
		{
			isf->nWaveVec++;										// Counter number of K
			// values that satisfy add to the count of the wave vectors
			for(a=0;a<s->Ato;a++)
			{	
				// for each of those wave vectors, we calculate e^(q*2p/l.(r(t)-r(0)))
				// make note of the fact that during taking the dot product, the components are split for q and l as well
				// we take the cos + isin form
				// store it in isf of that time step.
				// we need to store time of counter as well - or we need isf to be 2d 
				//with the values in one column and time in other
				isf->isf[isf->counter]+=cos(qx*2.0*M_PI/d->lBox[0]*(d->r[a][0]-isf->pos0[a][0])+qy*2.0*M_PI/d->lBox[1]*(d->r[a][1]-isf->pos0[a][1])+qz*2.0*M_PI/d->lBox[2]*(d->r[a][2]-isf->pos0[a][2]));		// q dot DeltaR
			}
		}
	}
	}
	}
	// normalisation below
	isf->isf[isf->counter]/=s->Ato*isf->nWaveVec;
	
	isf->counter++; 
    // increment counter
	if(isf->counter==s->Isf_limit) 
	{
		print_Isf(s,d,isf);
	}
	// print after some value
	// remember that looped function call is external so we 
	//need not worry about how many times  code above will run.
	
}

void print_Isf(SIMDAT *s, DATASIM *d, ISF *isf)
{	
	FILE *fp;
	int t;
	char buffer[150];

	sprintf(buffer,"%sIsfK%g_%d.dat",s->output,s->IsfK_value,isf->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=1;t<s->Isf_limit;t++)
	{
		// t stores the counter info
		// DeltaT tells us the LAMMPS time step.
		// s->Step tells us the no. of time steps skipped
		// we print time, isf(t)
		fprintf(fp,"%le	%le\n",s->DeltaT*s->Step*t,isf->isf[t]);
	}
	fclose(fp);
	// reset counter once we have printed
	// file counter tells us how many files we have made. 
	isf->fileCounter++;
	isf->counter=0;
	
}

