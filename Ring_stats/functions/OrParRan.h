struct ORDER
{
	int counter;          
	int fileCounter;
	int neigh;
	int bond;
	double **qBarG;
	double **qBarL;  
	double **timeG;
   	double **timeL;
	double **locDist;
	double **qloc;
} order;

void Init_OrPar(SIMDAT *s, DATASIM *d, ORDER *op)
{
	int l,m;

	srand(time(NULL));

	op->counter=0;
	op->fileCounter=0;
	op->qBarG=AllocMatR(10,21);
	op->qBarL=AllocMatR(10,21);
	op->timeG=AllocMatR(s->OrPar_limit,10);
	op->timeL=AllocMatR(s->OrPar_limit,10);
	op->locDist=AllocMatR(10,s->OrPar_bin);
	op->qloc=AllocMatR(s->Ato,10);

	for(l=0;l<10;l++)
	{
		for(m=0;m<21;m++)
		{
			op->qBarG[l][m]=0.0;
			op->qBarL[l][m]=0.0;
		}
		for(m=0;m<s->OrPar_limit;m++)
		{
			op->timeG[m][l]=0.0;
			op->timeL[m][l]=0.0;
		}
		for(m=0;m<s->OrPar_bin;m++)
		{
			op->locDist[l][m]=0.0;
		}
		for(m=0;m<s->Ato;m++)
		{
			op->qloc[m][l]=0.0;
		}
	}
}

int cutOff(SIMDAT *s, DATASIM *d, int a1, int a2)
{	
	int k;
	double r=0.0,trash;

	for(k=0;k<3;k++)
	{			
		trash=d->r[a2][k]-d->r[a1][k];
		trash-=d->lBox[k]*rint(trash/d->lBox[k]);	
		r+=trash*trash;
	}
	r=sqrt(r);
	if(r<=s->OrPar_r)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void qlm(DATASIM *d, ORDER *op, int a1, int a2)
{
	int l,m;
	double c1,c2,c3,trash;
	double r[4],cosTheta,phi;
	double random;

	for(l=0;l<3;l++)
	{
		r[l]=d->r[a2][l]-d->r[a1][l];
		r[l]-=d->lBox[l]*rint(r[l]/d->lBox[l]);
	}
	r[3]=sqrt(sqr(r[0])+sqr(r[1])+sqr(r[2]));

	cosTheta=r[2]/r[3];														// Spherical coordinate
	phi=atan2(r[1],r[0]);
	//phi=asin(r[1]/(sqrt(r[0]*r[0]+r[1]*r[1])));

	random=(double)rand() / (double)RAND_MAX;
	cosTheta=2.0*random-1.0;
	random=(double)rand() / (double)RAND_MAX;
	phi=2.0*M_PI*random;

	for(l=1;l<11;l++)
	{
		for(m=0;m<=l;m++)
		{
			op->qBarL[l-1][m]+=gsl_sf_legendre_sphPlm(l,m,cosTheta)*cos(m*phi);
		}
		for(m=1;m<=l;m++)
		{
			op->qBarL[l-1][m+l]+=gsl_sf_legendre_sphPlm(l,m,cosTheta)*cos(m*phi)*pow(-1,m);		
		}
	}
}

void local(ORDER *op, int a1)
{
	int l,m;

	for(l=0;l<10;l++)
	{
		for(m=0;m<21;m++)
		{
			op->qBarG[l][m]+=op->qBarL[l][m];
			op->qBarL[l][m]/=op->neigh;
			op->qBarL[l][m]*=op->qBarL[l][m];
			op->qloc[a1][l]+=op->qBarL[l][m];
			op->qBarL[l][m]=0.0;
		}
		op->qloc[a1][l]*=4*M_PI/(2*l+3);
		op->qloc[a1][l]=sqrt(op->qloc[a1][l]);		
	}
}

void global(ORDER *op)
{
	int l,m;

	for(l=0;l<10;l++)
	{
		for(m=0;m<21;m++)
		{
			op->qBarG[l][m]/=op->bond;
			op->qBarG[l][m]*=op->qBarG[l][m];
			op->timeG[op->counter][l]+=op->qBarG[l][m];
			op->qBarG[l][m]=0.0;
		}
		op->timeG[op->counter][l]*=4*M_PI/(2*l+3);
		op->timeG[op->counter][l]=sqrt(op->timeG[op->counter][l]);		
	}
}


void print_OrPar(SIMDAT *s, DATASIM *d, ORDER *op)
{
	FILE *fp;	
	int t,l;
	char buffer[150];
	double averL[10],averG[10];
	
	for(l=0;l<10;l++)
	{
		averL[l]=0.0;
		averG[l]=0.0;
	}

	sprintf(buffer,"%sQGtime%g-%d.dat",s->output,s->OrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w");
	
	for(t=1;t<s->OrPar_limit;t++)
	{
		fprintf(fp,"%f",s->DeltaT*s->Step*t);
		for(l=0;l<10;l++)
		{
			fprintf(fp,"	%f",op->timeG[t][l]);
			averG[l]+=op->timeG[t][l];
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(buffer,"%sQLtime%g-%d.dat",s->output,s->OrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w");
	
	for(t=1;t<s->OrPar_limit;t++)
	{
		fprintf(fp,"%f",s->DeltaT*s->Step*t);
		for(l=0;l<10;l++)
		{
			fprintf(fp,"	%f",op->timeL[t][l]/s->Ato);
			averL[l]+=op->timeL[t][l];
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(buffer,"%sQglob%g-%d.dat",s->output,s->OrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	 
	for(l=0;l<10;l++)
	{
		fprintf(fp,"%d	%f\n",l+1,averG[l]/s->OrPar_limit);
	}
	
	fclose(fp);

	sprintf(buffer,"%sQloc%g-%d.dat",s->output,s->OrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	 
	for(l=0;l<10;l++)
	{
		fprintf(fp,"%d	%f\n",l+1,averL[l]/(s->Ato*s->OrPar_limit));
	}
	
	fclose(fp);

	sprintf(buffer,"%sQlocDist%g-%d.dat",s->output,s->OrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	
	for(t=0;t<s->OrPar_bin;t++)
	{
		fprintf(fp,"%f",1.0/s->OrPar_bin*t);
		for(l=0;l<10;l++)
		{
			op->locDist[l][t]/= s->Ato*s->OrPar_limit;
			fprintf(fp,"	%f",op->locDist[l][t]*s->OrPar_bin);
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);

	op->fileCounter++;

}

void Compute_OrPar(SIMDAT *s, DATASIM *d, ORDER *op)
{
	int a1,a2,l,m;
	double r;

	op->bond=0;
	
	for(a1=0;a1<s->Ato;a1++)
	{
		op->neigh=0;
		for(a2=0;a2<s->Ato;a2++)
		{
			if(cutOff(s,d,a1,a2)==1 && a1!=a2)
			{
				op->neigh++;
				op->bond++;
				qlm(d,op,a1,a2);
			}
			
		}
		local(op,a1);
	}
	global(op);

	for(l=0;l<10;l++)
	{
		for(a1=0;a1<s->Ato;a1++)
		{
			op->timeL[op->counter][l]+=op->qloc[a1][l];
			op->qloc[a1][l]*=s->OrPar_bin;
			m=(int)op->qloc[a1][l];
			op->locDist[l][m]++;
			op->qloc[a1][l]=0.0;
		}
	}

	op->counter++;

	if(op->counter==s->OrPar_limit)
	{
		print_OrPar(s,d,op);
	}
}




