struct ENTROPY
{
	double **gm;
	int **neighlist;
	double *s, *sbar, *f, *I;
	int neighmax, counter, fileCounter;
	double dr, rho, Iaccum;
} entropy;

void Init_entropy(SIMDAT *s, DATASIM *d, ENTROPY *ep)
{
	double spam;
	int i,j;
	spam = s->eprm*100;
	ep->neighmax=(int) spam;
	ep->dr = 1.0*s->eprm/(s->epnbin-1); // the first bin corresponding to 0 is ignored to avoid NaN
	// therefore total number of bins is one fewer and total size of bins consequently larger
	ep->gm = AllocMatR(s->Ato,s->epnbin);
	ep->neighlist = AllocMatI(s->Ato,ep->neighmax);
	ep->s=AllocVecR(s->Ato);
	ep->sbar=AllocVecR(s->Ato);
	ep->f=AllocVecR(s->epnbin); // we need values either side of s->ra for continuity.
	ep->I=AllocVecR(s->epnbin);
	// therefore, and for other reasons, s->rm > s->ra always always
	printf("entered Init_entropy\n");
	for(i=0;i<s->Ato;i++)
	{
		ep->s[i]=0.0;
		ep->sbar[i]=0.0;
		for(j=0;j<s->epnbin;j++)
		{
			ep->gm[i][j]=0.0;
		}
		for(j=0;j<ep->neighmax;j++)
		{	
			ep->neighlist[i][j]=0;	
		}
	}
	for(i=0;i<s->epnbin;i++)
	{
		ep->f[i]=0.0;
		ep->I[i]=0.0;
	}
	ep->rho = s->Ato/(d->lBox[0]*d->lBox[1]*d->lBox[2]);
	printf("%f rho\n",ep->rho);
	printf("%f log(10)\n",gsl_sf_log(10));
	ep->Iaccum = 0.0;
	ep->counter=0;
	ep->fileCounter=0;
}

void print_entropy(SIMDAT *s, DATASIM *d, ENTROPY *ep)
{
	//print s and sbar atom-wise
	FILE *fp;	
	int t,i, spam;
	char buffer[300];
	
	sprintf(buffer,"%slocal_avg_entropy-%f-sig-%d.dat",s->output,s->epsigma,ep->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%d %f\n",t,ep->sbar[t]);
	}	
	fclose(fp);
	
	ep->fileCounter++;

}

void compute_f(SIMDAT *s, DATASIM *d, ENTROPY *ep)
{
	// we populate the f vector
	int i;
	double sr;
	FILE *fp;
	char buffer[300];
	sprintf(buffer,"%sfunction-f-%f.dat",s->output,s->epsigma);
	fp=fopen(buffer,"w");
	for(i=0;i<s->epnbin;i++)
	{
		sr = (i+1)*ep->dr; // i+1 to avoid nan elsewhere
		ep->f[i] = 1 - pow(sr/s->epra,6.0);
		ep->f[i] = ep->f[i]/(1 - pow(sr/s->epra,12.0));
		fprintf(fp,"%f %f\n",sr,ep->f[i]);
	}
	fclose(fp);
}

double feval( ENTROPY *ep,int i, int k) // evaluates the integrand
{
	double spr, It;
	spr = (k+1)*ep->dr; // k+1 to avoid nan while dividing by spr
	//printf("%f %f\n",spr, sqr(spr));
	if(ep->gm[i][k]!=0)
	{
	It = (ep->gm[i][k]*gsl_sf_log(ep->gm[i][k]) - ep->gm[i][k] + 1)*sqr(spr);
	}
	//if(k==2)
	//{	
	//printf("%f %f gm It for 2\n",ep->gm[i][k],It);
	//}
	return It;
}


void compute_entropy(SIMDAT *s, DATASIM *d, ENTROPY *ep)
{
	// first neighlist populate
	// then calculate the f make the lookup table
	// populate gm
	// calc s
	// calc sbar
	
	int i,j,k,si,sj,spam;
	double six,siy,siz,sjx,sjy,sjz,sdx,sdy,sdz,dist;
	double sumsf,sumf, spr,sig,spamsum,norma;
	FILE *fp, *fp2, *fp3;
	char buffer[300], buffer2[300], buffer3[300];
	sumsf = 0.0;
	sumf = 0.0;
	
	//populate neighlist
	s->epgm_lim = s->eprm; // vet this
	for(i=0;i<s->Ato-1;i++)
	{
		for(j=i+1;j<s->Ato;j++)
		{
			si = i; // atom id of i is just i
			sj = j; // atom id of j is just j
		    six = d->r[si][0];
		    siy = d->r[si][1];
		    sjx = d->r[sj][0];
		    sjy = d->r[sj][1];
		    siz = d->r[si][2];
		    sjz = d->r[sj][2];
		    sdx = six - sjx;
		    sdx-= d->lBox[0]*rint(sdx/d->lBox[0]);
		    sdy = siy - sjy;
		    sdy-= d->lBox[1]*rint(sdy/d->lBox[1]);
		    sdz = siz - sjz;
		    sdz-= d->lBox[2]*rint(sdz/d->lBox[2]);
		    dist = sqrt(sqr(sdx) + sqr(sdy) + sqr(sdz));
		    if(dist <=s->epgm_lim) // we use the limit to which gm_i(r) is calculated
		    {			 
   				ep->neighlist[si][0]++; // one more neigh for i
				ep->neighlist[sj][0]++; // one more neigh for j
				ep->neighlist[si][ep->neighlist[si][0]]=sj; // one more neighbour for i
				ep->neighlist[sj][ep->neighlist[sj][0]]=si; // one more neighbour for j
				//printf("%d %d neighed\n",ep->neighlist[si][0],ep->neighlist[sj][0]);
			}			 	
		}
		//printf("%d %d n neigh\n",i,ep->neighlist[i][0]);
	}

	// populate the f vector
	compute_f(s,d,ep); 
	sig = 2.0*s->epsigma*s->epsigma; 
	//calculate gm
	sprintf(buffer,"%sdiagnostic-gm-%f-sig.dat",s->output,s->epsigma);
	fp = fopen(buffer,"w");
	sprintf(buffer2,"%sdist_histo-gm-%f-sig.dat",s->output,s->epsigma);
	fp2 = fopen(buffer2,"w");
	sprintf(buffer3,"%sone-gaussian-test-%f-sig.dat",s->output,s->epsigma);
    fp3=fopen(buffer3,"w");
	for(i=0;i<s->Ato;i++)
	{
		//printf("%d %d \n",s->Ato,i);
		for(j=1;j<=ep->neighlist[i][0];j++)
		{
			si = i;
			sj = ep->neighlist[i][j];
			//printf("%d %d si sj\n",si,sj);
			six = d->r[si][0];
		    siy = d->r[si][1];
		    sjx = d->r[sj][0];
		    sjy = d->r[sj][1];
		    siz = d->r[si][2];
		    sjz = d->r[sj][2];
		    sdx = six - sjx;
		    sdx-= d->lBox[0]*rint(sdx/d->lBox[0]);
		    sdy = siy - sjy;
		    sdy-= d->lBox[1]*rint(sdy/d->lBox[1]);
		    sdz = siz - sjz;
		    sdz-= d->lBox[2]*rint(sdz/d->lBox[2]);
		    dist = sqrt(sqr(sdx) + sqr(sdy) + sqr(sdz));
		    fprintf(fp2,"%f\n",dist);
		    
		    for(k=0;k<s->epnbin;k++)
		    {
		    	spr = ep->dr*(k+1); // r position from bin position
		    	//printf("%d %f %f %f %f k dist, spr, spr-dist and exp\n",k,dist,spr,spr-dist,exp(-1.0*sqr(spr-dist)/sig));
		    	ep->gm[i][k] += exp((-1.0*sqr(spr-dist))/sig)/(sqrt(M_PI*sig)); // gaussian mean at dist
		    	//printf("%d %d %f epgm\n",i,k,ep->gm[i][k]);
		    	if(i==0 && j ==1)
		    	{		
		    		//printf("over here now!! %d %f %f\n",k,spr,ep->gm[i][k]);
		    		fprintf(fp3,"%f %f\n",spr,ep->gm[i][k]);	
		    	}
		    	// we add the contribution of the jth particle and divide by normalising factor;
		    }		
		}
		
		for(j=0;j<s->epnbin;j++)
		{
			spr = ep->dr*(j+1); // this j+1 is essential since we are dividing by spr
			ep->gm[i][j]=ep->gm[i][j]/(4.0*M_PI*ep->rho*spr*spr);
			if(i==0)
			fprintf(fp,"%f %f\n",spr,ep->gm[i][j]);
			// since this external normalisation is a function of r (spr) we need to do this in a loop
		}
		
	}
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	
	// print sum of all gaussians
	
	sprintf(buffer,"%ssum-gaussian-%f-sig.dat",s->output,s->epsigma);
	fp = fopen(buffer,"w");
	spamsum = 0.0;
	for(k=0;k<s->epnbin;k++)
	{
		spamsum = 0.0;
		spr = (k+1)*ep->dr;

		for(i=0;i<s->Ato;i++)
		{
			spamsum = spamsum + ep->gm[i][k];
		}
		spamsum = spamsum/s->Ato; // divide to take avg 
		fprintf(fp,"%f %f\n",spr,spamsum);
	}
	fclose(fp);
	// calculate s[i] by doing simpsons integration for each i
	//printf("%d %d %d epnbin int rint\n",s->epnbin,(int)(0.5*(3+s->epnbin)),(int)rint(0.5*(3+s->epnbin)));
	
	sprintf(buffer,"%sentropy-non-average-%f-sig.dat",s->output,s->epsigma);
	fp = fopen(buffer,"w");
	for(i=0;i<s->Ato;i++)
	{

		ep->Iaccum = 0.0;

		
		// formal Composite Simpsons 1-3 rule
		
		for(k=1;k<=(int)(0.5*s->epnbin);k++)
		{
		
			ep->Iaccum = ep->Iaccum + (feval(ep,i,2*k-2) + 4*feval(ep,i,2*k-1) + feval(ep,i,2*k));
		}
		ep->Iaccum =ep->Iaccum*(ep->dr/3.0);
		ep->s[i] = -1.0*ep->Iaccum*2*M_PI*ep->rho;
		fprintf(fp,"%d %f\n",i,ep->s[i]);
		

		////////// IMPORTANT TO CHECK THE VALUE OF KB IN REDUCED ENERGY UNITS ///////////
		////////// or make the distinction that this entropy is calculated in kb units///
	}
	fclose(fp);
	
	
	

	// now calculate sbar[i] taking summations etc etc
	for(i=0;i<s->Ato;i++)
	{
		sumsf = 0.0;
		sumf = 0.0;
		for(j=1;j<=ep->neighlist[i][0];j++)
		{
			si = i;
			sj = ep->neighlist[i][j];
			six = d->r[si][0];
		    siy = d->r[si][1];
		    sjx = d->r[sj][0];
		    sjy = d->r[sj][1];
		    siz = d->r[si][2];
		    sjz = d->r[sj][2];
		    sdx = six - sjx;
		    sdx-= d->lBox[0]*rint(sdx/d->lBox[0]);
		    sdy = siy - sjy;
		    sdy-= d->lBox[1]*rint(sdy/d->lBox[1]);
		    sdz = siz - sjz;
		    sdz-= d->lBox[2]*rint(sdz/d->lBox[2]);
		    dist = sqrt(sqr(sdx) + sqr(sdy) + sqr(sdz));
		    // find which bin dist belongs in
		    spam = (int)(dist/ep->dr); // check if this binning works ok
		    //printf("%f %f %d dist ep->dr spam\n",dist,ep->dr,spam);
		    //printf("%f %f epra f\n",s->epra,ep->f[spam]);
		    sumsf += ep->s[sj]*ep->f[spam];
		    sumf +=ep->f[spam];
		}
		//printf("%f %f\n",sumsf,sumf);
		ep->sbar[i] = (sumsf + ep->s[i])/(sumf + 1);	
	}
	ep->counter++;
	printf("entropy counter%d\n",ep->counter);
	if(ep->counter==s->entropy_limit)
	{
		printf("entering print ENTROPY!!\n");
		print_entropy(s,d,ep);
	}
	
}
