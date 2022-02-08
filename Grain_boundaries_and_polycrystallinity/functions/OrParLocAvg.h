//#include"YgCoordNo.h"
struct ORDERLA
{
	int counter;          
	int fileCounter;
	int neigh;
	int bond;
	double dr;
	double *f;
	double **qBarG_R;
	double **qBarL_R;
	double **qBarG_I;
	double **qBarL_I;  
	double **timeG;
   	double **timeL;
	double **locDist, **qlocavgDist;
	double **qloc, **qlocavg;
	double **q3CoordHisto;
	int **coordList;
	int maxCoord;
} orderla;

void Init_OrParla(SIMDAT *s, DATASIM *d, ORDERLA *op)
{
	int i,j,l,m;
	op->maxCoord = 30;
	op->counter=0;
	op->fileCounter=0;
	op->dr = 0.5*d->lBox[0]/s->laqnbin;
	printf("%f\n",op->dr);
	op->qBarG_R=AllocMatR(10,21);
	op->qBarL_R=AllocMatR(10,21);
	op->qBarG_I=AllocMatR(10,21);
	op->qBarL_I=AllocMatR(10,21);
	op->timeG=AllocMatR(s->laOrPar_limit,10);
	op->timeL=AllocMatR(s->laOrPar_limit,10);
	op->locDist=AllocMatR(10,s->laOrPar_bin);
	op->qlocavgDist=AllocMatR(10,s->laOrPar_bin);
	op->q3CoordHisto=AllocMatR(7,s->laOrPar_bin); //for coordination numbers!
	op->qloc=AllocMatR(s->Ato,10);
	op->qlocavg=AllocMatR(s->Ato,10);
	op->f=AllocVecR(s->laqnbin);
	op->coordList=AllocMatI(s->Ato,op->maxCoord); // list of particles, no. neigh and nerigh ids.
	for(l=0;l<10;l++)
	{
		for(m=0;m<21;m++)
		{
			op->qBarG_R[l][m]=0.0;
			op->qBarL_R[l][m]=0.0;
			op->qBarG_I[l][m]=0.0;
			op->qBarL_I[l][m]=0.0;
		}
		for(m=0;m<s->laOrPar_limit;m++)
		{
			op->timeG[m][l]=0.0;
			op->timeL[m][l]=0.0;
		}
		for(m=0;m<s->laOrPar_bin;m++)
		{
			op->locDist[l][m]=0.0; // check what this does
			op->qlocavgDist[l][m]=0.0;
		}
		for(m=0;m<s->Ato;m++)
		{
			op->qloc[m][l]=0.0; // q_local for each particle
			op->qlocavg[m][l]=0.0;
		}
	}
	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<op->maxCoord;j++)
		{
		op->coordList[i][j]=0;
		}

	}
	for(i=0;i<7;i++)
	{
		for(j=0;j<s->laOrPar_bin;j++)
		{
		op->q3CoordHisto[i][j]=0.0;		
		}

	}
	for(i=0;i<s->laqnbin;i++)
	{
		op->f[i]=0.0;
	}
}

void compute_fla(SIMDAT *s, DATASIM *d, ORDERLA *op)
{
	// we populate the f vector
	int i;
	double sr;
	FILE *fp;
	char buffer[500];
	sprintf(buffer,"%sfunction-f-%f.dat",s->output,s->laOrPar_ra);
	fp=fopen(buffer,"w");
	printf("remember that we're using s->laqnbin here!!!!!\n");
	for(i=0;i<s->laqnbin;i++)
	{
		sr = (i+1)*op->dr; // i+1 to avoid nan elsewhere
		op->f[i] = 1 - pow(sr/s->laOrPar_ra,6.0);
		op->f[i] = op->f[i]/(1 - pow(sr/s->laOrPar_ra,12.0));
		fprintf(fp,"%f %f\n",sr,op->f[i]);
	}
	fclose(fp);
}


int cutOffla(SIMDAT *s, DATASIM *d, int a1, int a2)
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
	if(r<=s->laOrPar_r)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	// flags whether j is within cut-off distance of i
}

void qlmla(DATASIM *d, ORDERLA *op, int a1, int a2)
{
	int l,m;
	double c1,c2,c3,trash;
	double r[4],cosTheta,phi;

	for(l=0;l<3;l++)
	{
		r[l]=d->r[a2][l]-d->r[a1][l];
		r[l]-=d->lBox[l]*rint(r[l]/d->lBox[l]);
	}
	r[3]=sqrt(sqr(r[0])+sqr(r[1])+sqr(r[2]));

	cosTheta=r[2]/r[3];														// Spherical coordinate
	phi=atan2(r[1],r[0]);
	
	for(l=1;l<11;l++)
	{
		for(m=0;m<=l;m++)
		{
			op->qBarL_R[l-1][m]+=gsl_sf_legendre_sphPlm(l,m,cosTheta)*cos(m*phi);
			op->qBarL_I[l-1][m]+=gsl_sf_legendre_sphPlm(l,m,cosTheta)*sin(m*phi);
		}
		for(m=1;m<=l;m++)
		{
			op->qBarL_R[l-1][m+l]+=gsl_sf_legendre_sphPlm(l,m,cosTheta)*cos(-m*phi)*pow(-1,m);
			op->qBarL_I[l-1][m+l]+=gsl_sf_legendre_sphPlm(l,m,cosTheta)*sin(-m*phi)*pow(-1,m);		
		}
	}
	// calculates and updates local order parameter
	// must pay attention to how many times this subroutine is called.
	// every time it's called, it makes a set of increments to the different l and m components
}

void localla(ORDERLA *op, int a1)
{
	int l,m;
	// this subroutine is specific to each atom
	// so once we calculate qBar etc, 
	// we update the global by adding the contribution for a th atom
	// local is simply whatever it is for that atom.
	// normalise local (global normalise later)
	// and qloc(a)(l) is simply the real part
	// reset the qBarL for next iteration
	for(l=0;l<10;l++)
	{
		for(m=0;m<21;m++)
		{
			op->qBarG_R[l][m]+=op->qBarL_R[l][m];
			op->qBarG_I[l][m]+=op->qBarL_I[l][m];
			if(op->neigh!=0)
			{
			op->qBarL_R[l][m]/=op->neigh;
			op->qBarL_I[l][m]/=op->neigh;
			}
			// if op->neigh is 0 then we don't divide, we have a 0/0 form.
			// as long as we have at least one neighbour
			op->qBarL_R[l][m]=sqr(op->qBarL_R[l][m])+sqr(op->qBarL_I[l][m]); // square of mod??
			op->qloc[a1][l]+=op->qBarL_R[l][m];
			op->qBarL_R[l][m]=0.0;
			op->qBarL_I[l][m]=0.0;
		}
		op->qloc[a1][l]*=4*M_PI/(2*l+1);
		op->qloc[a1][l]=sqrt(op->qloc[a1][l]);		
	}
	// this subroutine does the normalising that would be appropriate.
}

void globalla(ORDERLA *op)
{
	int l,m;
    printf("entering global la\n");
	for(l=0;l<10;l++)
	{
		for(m=0;m<21;m++)
		{
			if(op->bond!=0)
			{
			op->qBarG_R[l][m]/=op->bond;
			op->qBarG_I[l][m]/=op->bond;
			}
			else
			printf("ooops!! had to worry about Qglob\n");
			op->qBarG_R[l][m]=sqr(op->qBarG_R[l][m])+sqr(op->qBarG_I[l][m]);
			op->timeG[op->counter][l]+=op->qBarG_R[l][m];
			op->qBarG_R[l][m]=0.0;
			op->qBarG_I[l][m]=0.0;
		}
		op->timeG[op->counter][l]*=4*M_PI/(2*l+1);
		op->timeG[op->counter][l]=sqrt(op->timeG[op->counter][l]);
		printf("%d %f how's it??\n",op->counter,op->timeG[op->counter][l]);		
	}
	// this subroutine is for normalising as well... but for global
}


void print_OrParla(SIMDAT *s, DATASIM *d, ORDERLA *op)
{
	FILE *fp;	
	int t,l;
	double sx,sy,sz;
	char buffer[500];
	double averL[10],averG[10];
	
	for(l=0;l<10;l++)
	{
		averL[l]=0.0;
		averG[l]=0.0;
	}

	sprintf(buffer,"%sQGtime%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w");
	
	for(t=0;t<s->laOrPar_limit;t++)
	{
		printf("we need to bother ourselves with this??\n");
		fprintf(fp,"%f",s->DeltaT*s->Step*t);
		for(l=0;l<10;l++)
		{
			fprintf(fp,"	%f",op->timeG[t][l]);
			averG[l]+=op->timeG[t][l];
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(buffer,"%sQLtime%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w");
	
	for(t=0;t<s->laOrPar_limit;t++)
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

	sprintf(buffer,"%sQglob%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	 
	for(l=0;l<10;l++)
	{
		fprintf(fp,"%d	%f\n",l+1,averG[l]/(s->laOrPar_limit));
	}
	
	fclose(fp);

	sprintf(buffer,"%sQloc%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	 
	for(l=0;l<10;l++)
	{
		fprintf(fp,"%d	%f\n",l+1,averL[l]/(s->Ato*(s->laOrPar_limit)));
	}
	
	fclose(fp);

	sprintf(buffer,"%sQlocDist%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	
	for(t=0;t<s->laOrPar_bin;t++)
	{
		fprintf(fp,"%f",1.0/s->laOrPar_bin*t);
		for(l=0;l<10;l++)
		{
			op->locDist[l][t]/= s->Ato*s->laOrPar_limit;
			fprintf(fp,"	%f",op->locDist[l][t]*s->laOrPar_bin);
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);
	
	sprintf(buffer,"%sQlocAvgDist%g-ra%g-%d.dat",s->output,s->laOrPar_r,s->laOrPar_ra,op->fileCounter);
  	fp=fopen(buffer,"w"); 
	
	for(t=0;t<s->laOrPar_bin;t++)
	{
		fprintf(fp,"%f",1.0/s->laOrPar_bin*t);
		for(l=0;l<10;l++)
		{
			op->qlocavgDist[l][t]/= s->Ato*s->laOrPar_limit;
			fprintf(fp,"	%f",op->qlocavgDist[l][t]*s->laOrPar_bin);
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);
	
	sprintf(buffer,"%sQlocSepDistOvito%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
  	fp=fopen(buffer,"w"); 
  	fprintf(fp,"%d\n",s->Ato);
  	fprintf(fp,"random text\n");
  	for(t=0;t<s->Ato;t++)
  	{
  		// print ovito file where we separate out the atoms with high q3
  		sx = d->r[t][0] + 0.5*d->lBox[0]; // put all in bounds of 0 - lBox
  		sy = d->r[t][1] + 0.5*d->lBox[1];
  		sz = d->r[t][2] + 0.5*d->lBox[2];
  		if(op->qlocavg[t][2]>0.8)
  		{
  			fprintf(fp,"SiHigh %f %f %f\n",sx,sy,sz);
  		}
  		if(op->qlocavg[t][2]<0.8 && op->qloc[t][2]>0.7)
  		{
  			fprintf(fp,"Sic %f %f %f\n",sx,sy,sz);
  		}
  		else
  		{
  			fprintf(fp,"SiH %f %f %f\n",sx,sy,sz);
  		}
  	}
  	fclose(fp);
	
	
	
	sprintf(buffer,"%sQloc3CoordHisto%g-%d.dat",s->output,s->laOrPar_r,op->fileCounter);
	fp = fopen(buffer,"w");
	for(t=0;t<s->laOrPar_bin;t++)
	{
		fprintf(fp,"%f",1.0*t/s->laOrPar_bin); // you needd the 1.0 for float conversion
		for(l=0;l<7;l++) // this is for the 7 coordination numbers
		{
			op->q3CoordHisto[l][t]/= s->Ato*s->laOrPar_limit; //normalised wrt no atoms and no files
			fprintf(fp," %f",op->q3CoordHisto[l][t]*s->laOrPar_bin); // size of bin needed for op
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	op->fileCounter++;
    printf("done printing or par local average\n");
}

void Compute_OrParla(SIMDAT *s, DATASIM *d, ORDERLA *op)
{
	int a1,a2,l,m,mspam;
	double r, trash,spam,sumof,sumf;
	int i,j,k,NeiCount,temp,spamint;	
	double posR[4];
	op->bond=0;
	printf("entering loc avg order param cals\n");
	for(a1=0;a1<s->Ato;a1++)
	{
		op->neigh=0;
        //printf("atom %d is getting done\n",a1);
		for(a2=0;a2<s->Ato;a2++)
		{
            //printf("calculating cutoffla!!\n");
			if(cutOffla(s,d,a1,a2)==1 && a1!=a2)
			{
				op->neigh++;
				op->bond++;
				qlmla(d,op,a1,a2);
			}
			
		}
		localla(op,a1);
	}
    printf("calculating global la!\n ");

	globalla(op);
	printf("done with compute orpar la\n");
	
	// insert code to calc coordination number of all atoms here

	// reset the coord/neighbour list
	for(i=0;i<s->Ato;i++)
	{	
		for(j=0;j<op->maxCoord;j++)
		{
			op->coordList[i][j]=0;
		}
	}
	// for each atom
	printf("or par coord no entered\n");
	NeiCount =0;
	for(i=0;i<s->Ato;i++)
	{
		//printf("really now, what's happening!!!!\n");
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
					//printf("entered if for %d and %d with cut-off %f %f\n",i,j,posR[3],s->cn_cutoff);
					// store j as neigh of i, list[i][0] will be incremented for one more neighbour.
					//printf("entered if for %d and %d with cut-off %f %f\n",i,j,posR[3],s->cn_cutoff);
					op->coordList[i][0]++; // one more neighbour for i
					op->coordList[i][NeiCount]=j; // id of NeiCountth neighbour
					NeiCount++; // for the next neighbour
				
				}
			}
		}
		//printf("finished atom no %d\n",i);
		// now we know how many neighbours i has. We can update the histogram now
		// first find out which bin it will go in
		temp = op->coordList[i][0]; // <-is an int
		// then increment the value at that bin index.

	}	
	printf("or par coord no exited\n");
	
	
	
	
	// the code below updates 2 quantities
	// 1) the time counter adds the qloc of each particle for each l (subsequent division before print)
	// 2) finds the bin to which the qloc of the given particle belongs
	// 3) updates the distribution for the given l at a given bin location
	// 4) resets value for next run?
	
	// changes by yagyik
	// 1) we have 7 vs OrPar_bin array
	// there we shall store separate histograms for 7 coordination numbers for q3 only.
	// we shall see what qllocal values are taken by atoms wih different coord nos.
	// we are actually only interested in q3 but let's see. 
	
	compute_fla(s,d,op);
	printf("compute OPLA f\n");
	sumof=0.0;
	sumf=0.0;
	for(l=0;l<10;l++)
	{
		for(a1=0;a1<s->Ato;a1++)
		{
			sumof = 0.0;
			sumf = 0.0;
			for(a2=0;a2<s->Ato;a2++)
			{
				for(k=0;k<3;k++)
				{			
					trash=d->r[a2][k]-d->r[a1][k];
					trash-=d->lBox[k]*rint(trash/d->lBox[k]);	
					r+=trash*trash;
				}
				r=sqrt(r);
				spamint = (int)(r/op->dr); // check if this binning works ok
		    	//printf("%f %f %d dist ep->dr spam\n",dist,ep->dr,spam);
		    	//printf("%f %f epra f\n",s->epra,ep->f[spam]);
		    	sumof += op->qloc[a2][l]*op->f[spamint];
		    	//if(a1==0 && l==2)
		    	//printf("%d %f %f %f\n",a2,op->f[spamint],op->qloc[a2][l],sumof);
		    	sumf +=op->f[spamint];
			}
			op->qlocavg[a1][l] = (sumof + op->qloc[a1][l])/(sumf + 1);
			op->qlocavg[a1][l]*=s->laOrPar_bin;
			m = (int)op->qlocavg[a1][l];
			op->qlocavgDist[l][m]++;
			op->qlocavg[a1][l]=0.0;
		}
	}
	
	printf("done with opla calcs\n");
	
	
	for(l=0;l<10;l++)
	{
		for(a1=0;a1<s->Ato;a1++)
		{
			op->timeL[op->counter][l]+=op->qloc[a1][l];
			op->qloc[a1][l]*=s->laOrPar_bin;
			//printf("what is this a1 thing?? %d %f\n",a1,op->qloc[a1][l]);
			m=(int)op->qloc[a1][2];
			// now we have the bin to which that qloc belongs
			// we pull up the coordination no of that atom
			mspam = (int)op->coordList[a1][0]; // this is the row index
			//printf("%d %d how are we doing here?\n",mspam,m); 
			if(mspam<7)
			{
			op->q3CoordHisto[mspam][m]++; // one more entry for that coord no for that q3loc bin			
			}
			m=(int)op->qloc[a1][l];
			op->locDist[l][m]++;
			op->qloc[a1][l]=0.0;
			
		}
	}
	printf("done with the standard qloc dist\n" );
	op->counter++;

	if(op->counter==s->laOrPar_limit)
	{
		printf("printing orpar bobo la!!");
		print_OrParla(s,d,op);
	}
}





