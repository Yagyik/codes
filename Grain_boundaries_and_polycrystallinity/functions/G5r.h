struct G5R
{
	double *g5r;			// Fifth neighbor distribution
	double *quinto;
	int nNeigh;
	int numField;			// Per ora metto 50
	int counter;            
	int fileCounter;        	
} gCinqueR;

void Init_G5r(SIMDAT *s, DATASIM *d, G5R *r)
{
	int n;
	
	r->numField=50;
	r->g5r=AllocVecR(r->numField);
	r->quinto=AllocVecR(30);

	for(n=0;n<r->numField;n++)
	{ 
		r->g5r[n]=0.0;
	}
	
	r->counter=0;
	r->fileCounter=0;
}

void sort(G5R *r)
{
	int i,j,min;
	double temp;

	for(i=0;i<r->nNeigh-1;i++)
	{
		min=i;
		for(j=i+1;j<r->nNeigh;j++)
		{
			if(r->quinto[j] < r->quinto[min]) //cambiare questa condizione per invertire l'ordine
   			min = j;
		}
		temp=r->quinto[min];
		r->quinto[min]=r->quinto[i];
		r->quinto[i]=temp;
	}
}

void print_G5r(SIMDAT *s, DATASIM *d, G5R *r)
{
	FILE *fp;
	char buffer[500];
	int i;
	double deltaR;

	sprintf(buffer,"%sG5r_%d.dat",s->output,r->fileCounter);
	fp=fopen(buffer,"w");		

	deltaR=s->G5r_max/r->numField;

	for(i=0;i<r->numField;i++)
	{
		fprintf(fp,"%le	%le	%le\n",((double)i+0.5)*deltaR,r->g5r[i]/s->Ato/r->counter,r->g5r[i]);
	}
	fclose(fp);
	
	r->fileCounter++;
	r->counter=0;
}

void Compute_G5r(SIMDAT *s, DATASIM *d, G5R *r)
{
	int a1,a2,k,n;
	double posR[4],temp;

	for(a1=0;a1<s->Ato;a1++)
	{
		r->nNeigh=0;
		for(a2=0;a2<s->Ato;a2++)
		{
			if(a1!=a2) 
			{
				for(k=0;k<3;k++)
				{
					posR[k]=d->r[a2][k]-d->r[a1][k];
					posR[k]-=d->lBox[k]*rint(posR[k]/d->lBox[k]);
				}
				posR[3]=sqrt(sqr(posR[0])+sqr(posR[1])+sqr(posR[2]));
				if(posR[3]<s->G5r_max)
				{
					if(r->nNeigh>29)
					{
						printf("Errore G5R: cambiare max numero di vicini (30). Trovati %d\n",r->nNeigh);
						exit(0);
					}

					r->quinto[r->nNeigh]=posR[3];
					r->nNeigh++;
				}		
			}
		}
		if(r->nNeigh<5)
		{
			printf("Errore G5R: cambiare max R. Trovati solo %d vicini per atomo %d\n",r->nNeigh,a1);
			exit(0);
		}
		sort(r);

		temp=r->quinto[4]/s->G5r_max*r->numField;
		n=(int)(temp);
		r->g5r[n]++;
	}
		
	r->counter++;
			
	if(r->counter==s->G5r_limit)
	{
		print_G5r(s,d,r);
	}
}


