struct ESSEQ
{
	double **Sq;
	int counter;            
	int fileCounter;        	
} esseq;

void Init_Sq(SIMDAT *s, DATASIM *d, ESSEQ *sq)
{
	int i;

	sq->Sq=AllocMatR(s->Sq_bin,2);
	
	for(i=0;i<s->Sq_bin;i++)
	{
		sq->Sq[i][0]=sq->Sq[i][1]=0.0;
	}
	
	sq->counter=0;
	sq->fileCounter=0;
}

void print_Sq(SIMDAT *s, DATASIM *d, ESSEQ *sq)
{
	FILE *fp;
	char buffer[400];
	int q;
	double deltaQ;

	sprintf(buffer,"%sSq_%d.dat",s->output,sq->fileCounter);
	fp=fopen(buffer,"w");
	printf("%s - output file -> print Sq \n",buffer);
	deltaQ=s->Sq_max/s->Sq_bin;
	for(q=0;q<s->Sq_bin;q++)
	{
		if(sq->Sq[q][1]!=0)
		{
			sq->Sq[q][0]/=sq->Sq[q][1];
		}
		fprintf(fp,"%f	%f	%f\n",deltaQ*(q+0.5),sq->Sq[q][0],sq->Sq[q][1]);
	}

	fclose(fp);

	sq->fileCounter++;
	sq->counter=0;
}

void Compute_Sq(SIMDAT *s, DATASIM *d, ESSEQ *sq)
{
	int qx,qy,qz,qthr,x,a;
	double q,sumCos,sumSin;

	qthr=(int)(s->Sq_max*d->lBox[0]/(2.0*M_PI));

	for(qx=-qthr;qx<=qthr;qx++)
	{
	for(qy=-qthr;qy<=qthr;qy++)
	{
	for(qz=0;qz<=qthr;qz++)
	{
		q=2.0*M_PI/d->lBox[0]*sqrt(qx*qx+qy*qy+qz*qz); 				// Modulus of K
		if(q<s->Sq_max && !(qx==0 && qy==0 && qz==0))				// Avoid qx=qy=qz=0.0
		{
			sumCos=sumSin=0.0;
			for(a=0;a<s->Ato;a++)
			{
				sumCos+=cos(qx*2.0*M_PI/d->lBox[0]*d->r[a][0]+qy*2.0*M_PI/d->lBox[1]*d->r[a][1]+qz*2.0*M_PI/d->lBox[2]*d->r[a][2]);		// cos(q r_i)
				sumSin+=sin(qx*2.0*M_PI/d->lBox[0]*d->r[a][0]+qy*2.0*M_PI/d->lBox[1]*d->r[a][1]+qz*2.0*M_PI/d->lBox[2]*d->r[a][2]);		// sin(q r_i)
			}
			
			sumCos=(sqr(sumCos)+sqr(sumSin))/s->Ato;			// Sq= 1/N [ cos()^2 + sin()^2 ]
			
			sumSin=q*s->Sq_bin/s->Sq_max;						// For histogram
			x=(int)sumSin;
			
			if(qz==0)
			{
				sq->Sq[x][0]+=sumCos;
				sq->Sq[x][1]++;
			}
			else
			{
				sq->Sq[x][0]+=2.0*sumCos;
				sq->Sq[x][1]+=2;
			}
			
		}
	}
	}
	}

	sq->counter++;

	if(sq->counter==s->Sq_limit)
	{
		print_Sq(s,d,sq);
	}
}


