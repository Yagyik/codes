struct COORDNO
{
	int counter;          
	int fileCounter;
	double *coordDist; // histogram shows frequency of each coordination number
	int **coordList; // stores coord number of each atom and their neighbours
	double maxCoord; // max value for binning
	double *avgCoordNo; // stores the average coordination no. of the sample at each time step
} CoordNo;

void Init_CoordNo(SIMDAT *s, DATASIM *d, COORDNO *cn)
{
	int i,j;
	char buffer[500];
	cn->maxCoord = 30;
	cn->coordDist=AllocVecR(cn->maxCoord); // init real becaus normalising for histograms
	cn->coordList=AllocMatI(s->Ato,cn->maxCoord); // list of particles, no. neigh and nerigh ids.
	cn->avgCoordNo=AllocVecR(s->cn_limit);
	for(i=0;i<cn->maxCoord;i++)
	{
		cn->coordDist[i]=0;
	}
	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<cn->maxCoord;j++)
		{
		cn->coordList[i][j]=0;
		}

	}
	for(i=0;i<s->cn_limit;i++)
	{
	cn->avgCoordNo[i]=0.0;
	}
	cn->counter=0;
	cn->fileCounter=0;
	printf("coordination number, all allocated??\n");
}

void print_CoordNo(SIMDAT *s, DATASIM *d, COORDNO *cn)
{
	FILE *fp;
	FILE *fp2;	
	int t;
	char buffer[500];
	
	sprintf(buffer,"%sCoordNoVTime%g-%d.dat",s->output,s->cn_cutoff,cn->fileCounter);
  	fp=fopen(buffer,"w");
  	// here we shall print the average vs time
  	// we shall also print the histogram!
	for(t=1;t<s->cn_limit;t++)
	{
		fprintf(fp,"%le	%f\n",s->DeltaT*s->Step*t,cn->avgCoordNo[t]);
	}
	fclose(fp);
	sprintf(buffer,"%sCoordNoDist%g-%d.dat",s->output,s->cn_cutoff,cn->fileCounter);
	fp2=fopen(buffer,"w");
	for(t=1;t<cn->maxCoord;t++)
	{	
		cn->coordDist[t] /= s->cn_limit; // normalise against
		fprintf(fp2,"%d %f\n",t,cn->coordDist[t]);
	
	}
	fclose(fp2);
	cn->fileCounter++;

}

void Compute_CoordNo(SIMDAT *s, DATASIM *d, COORDNO *cn)
{
	int i,j,k,NeiCount,temp;
	printf("coord no, calculation\n");
	double posR[4];
	// reset the coord/neighbour list
	for(i=0;i<s->Ato;i++)
	{	
		for(j=0;j<cn->maxCoord;j++)
		{
			cn->coordList[i][j]=0;
		}
	}
	// for each atom
	NeiCount =0;
	for(i=0;i<s->Ato;i++)
	{
		
		// for each j other atoms
		NeiCount = 1;
		for(j=0;j<s->Ato;j++)
		{
			if(i!=j)
			{
				// if dist betw i and j < cutoff
				for(k=0;k<3;k++)
				{
					posR[k]=d->r[i][k]-d->r[j][k];
					posR[k]-=d->lBox[k]*rint(posR[k]/d->lBox[k]);
				}
				posR[3]=sqrt(sqr(posR[0])+sqr(posR[1])+sqr(posR[2]));
				if(posR[3]<s->cn_cutoff) // cn_cutoff is cut-off in anafile, read in structures.h
				{
					printf("entered if for %d and %d\n",i,j);
					// store j as neigh of i, list[i][0] will be incremented for one more neighbour.
					cn->coordList[i][0]++; // one more neighbour for i
					cn->coordList[i][NeiCount]=j; // id of NeiCountth neighbour
					NeiCount++; // for the next neighbour
				
				}
			}
		}
//		printf("finished atom no %d\n",i);
		// now we know how many neighbours i has. We can update the histogram now
		// first find out which bin it will go in
		temp = cn->coordList[i][0]; // <-is an int
		// then increment the value at that bin index.
		cn->coordDist[temp]++; 
		cn->avgCoordNo[cn->counter] = cn->avgCoordNo[cn->counter] + cn->coordList[i][0]/s->Ato;
		// the average is over particles
	}	
	cn->counter++;
	printf("coord no exited\n");
	if(cn->counter==s->cn_limit)
	{
		print_CoordNo(s,d,cn);
	}

}
