struct DENSITY
{
	double *gdr;			// RDF
	double *sdk;			// SSF
	int numField;			// Number of fields g(r)
	int counter;            
	int fileCounter;        	
} densityK;

void Init_Gr(SIMDAT *s, DATASIM *d, DENSITY *r)
{
	int n;

	r->numField=1;    
	while(r->numField<20*sqrt(s->Ato))
	{
		r->numField*=2;
	} 
	
	r->numField/=2;
	r->gdr=AllocVecR(r->numField);
	r->sdk=AllocVecR(r->numField*FactPrec);

	for(n=0;n<r->numField;n++)
	{ 
		r->gdr[n]=0.0;
	}
	for(n=r->numField;n<r->numField*FactPrec;n++)
	{
		r->sdk[n]=0.0;
	}

	r->counter=0;
	r->fileCounter=0;
}

void prepare_Fact(SIMDAT *s, DATASIM *d, DENSITY *r)
{
	int n;
	double deltaR,a,*input,*output;
	fftw_plan p;

	input=(double *) malloc(sizeof(double)*FactPrec*r->numField);
	output=(double *) malloc(sizeof(double)*FactPrec*r->numField);
	deltaR=0.5*d->lBox[0]/r->numField;
	
	for(n=0;n<r->numField;n++)
	{	 
		output[n]=0.0;
		input[n]=(r->gdr[n]-1.0)*(n+1)*deltaR;
	}  
	for(n=r->numField;n<FactPrec*r->numField;n++)
	{ 
		output[n]=0.0;
		input[n]=0.0;
	}

	p=fftw_plan_r2r_1d(FactPrec*r->numField,input,output,FFTW_RODFT11,FFTW_MEASURE);
	fftw_execute(p);
	
	a=2.0*M_PI*s->Ato*FactPrec/(r->numField*sqr(d->lBox[0]));

	for(n=0;n<FactPrec*r->numField;n++)
	{
		r->sdk[n]=output[n]*a/(n+1);
	}
	
	fftw_destroy_plan(p);
	fftw_free(input); fftw_free(output);
}

void print_Gr(SIMDAT *s, DATASIM *d, DENSITY *r)
{
	FILE *fp;
	char buffer[400];
	int i;
	double deltaR,norma,rho;
	double *I;
	
	I=AllocVecR(r->numField);

	sprintf(buffer,"%sGr_%d.dat",s->output,r->fileCounter);
	fp=fopen(buffer,"w");		
	printf("%s - output file ->print gr\n",buffer);
	deltaR=0.5*d->lBox[0]/r->numField;
	norma=2.0*pow(r->numField,3.0)/(M_PI*sqr(s->Ato)*s->Gr_limit);

	for(i=0;i<r->numField;i++)
	{
		r->gdr[i]*=norma/sqr(i+1);
		fprintf(fp,"%le	%le\n",((double)i+1)*deltaR,r->gdr[i]);
		I[i]=0.0;
	}
	fclose(fp);
	
	// here we include a subroutine to perform simpson's integration on g(r)
	sprintf(buffer,"%sNr_%d.dat",s->output,r->fileCounter);
	fp=fopen(buffer,"w");
	rho = s->Ato/(d->lBox[0]*d->lBox[1]*d->lBox[2]);	
	printf("density %f\n",rho);
	for(i=2;i<r->numField;i=i+1) // i+1 is very important! leads to errors if i+2
	{
		I[i]=I[i-1] + deltaR/6*(r->gdr[i-2]*sqr((i-1)*deltaR) + 4*r->gdr[i-1]*sqr((i)*deltaR) + r->gdr[i]*sqr((i+1)*deltaR));
		// when i=i+1 then increment for each i... (i+1)*deltaR corresponds to gdr(i) - need to remember this
		fprintf(fp,"%f %f\n",((double)i+1)*deltaR,4*M_PI*rho*I[i]);
	}
	fclose(fp);
	
	
	
	prepare_Fact(s,d,r);

	norma=2.0*M_PI/(d->lBox[0]*FactPrec);   
	sprintf(buffer,"%sFourGr_%d.dat",s->output,r->fileCounter);
	fp=fopen(buffer,"w");
	for(i=0;i<r->numField*FactPrec/6.0;i++)
	{
		fprintf(fp,"%le	%le\n",((double)i+1)*norma,r->sdk[i]+1.0);
	}
	fclose(fp);

	for(i=0;i<r->numField;i++)
	{ 
		r->gdr[i]=0.0;
	}
	for(i=r->numField;i<r->numField*FactPrec;i++)
	{
		r->sdk[i]=0.0;
	}

	r->fileCounter++;
	r->counter=0;
}

void Compute_Gr(SIMDAT *s, DATASIM *d, DENSITY *r)
{
	int a1,a2,k,n;
	double posR[4];
	printf("starting g(r) %d\n",r->counter);
	for(a1=0;a1<s->Ato;a1++)
	{
	for(a2=0;a2<s->Ato;a2++)
	{
		if(a1!=a2) 
		{
			for(k=0;k<3;k++)
			{
				posR[k]=d->r[a2][k]-d->r[a1][k];
				posR[k]-=d->lBox[k]*rint(posR[k]/d->lBox[k]);
			}
			posR[3]=sqrt(sqr(posR[0])+sqr(posR[1])+sqr(posR[2]))*2.0*r->numField/d->lBox[0];
			n=rint(posR[3]);
			if(n<=r->numField)
			{
				r->gdr[n-1]++;
			}
		}
	}
	}
	
	r->counter++;
			
	if(r->counter==s->Gr_limit)
	{
		printf("Gr_limit!!\n");
		print_Gr(s,d,r);
	}
}


