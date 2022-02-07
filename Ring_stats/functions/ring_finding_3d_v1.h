struct SPRING
{
	// visited, inPath, complete,pathStack, hitlist
	// AllNeigh,ringStats,pairDist,pathStack,nucleo
	int *visited,*inpath,*complete,*pathstack,*hitlist,*ringcount, *neigh;
	int **AllNeigh, **ringstats,**nucleo, **pairDist,**neighList;
	int stacklength,counter,fileCounter;
	
	double *r, *q3;
	double **q3mR,**q3mI;
	//inRings, 
	
	
	// for actual, q3 related stuff comes in here. 

} spring;

void Init_spring_3d(SIMDAT *s, DATASIM *d, SPRING *sp)
{
	int i,j;
	// pairDist init to +infty and diagonals to 0 
	sp->visited=AllocVecI(s->Ato);
	sp->inpath=AllocVecI(s->Ato);
	sp->complete=AllocVecI(s->Ato);
	sp->pathstack=AllocVecI(20); // this number 20 is max ring length of interest
	// 20 should probably be replaced with some parameter but we'll leave it be for now
	sp->hitlist=AllocVecI(s->Ato);
	sp->ringcount=AllocVecI(20);
	sp->AllNeigh=AllocMatI(s->Ato,15);
	sp->ringstats=AllocMatI(s->Ato,20);
	sp->pairDist=AllocMatI(s->Ato,s->Ato);
	sp->nucleo=AllocMatI(s->Ato+1,2); // only 2 kinds
	sp->neighList=AllocMatI(s->Ato,15);
	sp->neigh=AllocVecI(s->Ato);
	sp->r=AllocVecR(4);
	sp->q3=AllocVecR(s->Ato);
	sp->q3mR=AllocMatR(s->Ato,7);
	sp->q3mI=AllocMatR(s->Ato,7);
	sp->stacklength=0;
	sp->counter=0;
	sp->fileCounter=0;
	for(i=0;i<s->Ato;i++)
	{
		sp->visited[i]=0;
		sp->inpath[i]=0;
		sp->complete[i]=0;
		sp->hitlist[i]=0;
		sp->nucleo[i][0]=0;
		sp->nucleo[i][1]=0;
		sp->neigh[i]=0;
		sp->q3[i]=0.0;
		for(j=0;j<20;j++) // need to remember to change this 20 as well
		{
			sp->pathstack[j]=0;
			sp->ringcount[j]=0;
			sp->ringstats[i][j]=0;
			if(j<10)
			{
				sp->AllNeigh[i][j]=0;
				sp->neighList[i][j]=0;
			}
		}
		for(j=0;j<s->Ato;j++)
		{
			if(i==j)
			sp->pairDist[i][j]=0;
			else
			sp->pairDist[i][j]=999;
		}
		for(j=0;j<7;j++)
		{
			sp->q3mR[i][j]=0.0;
			sp->q3mI[i][j]=0.0;
		}
	}
	for(i=0;i<4;i++)
	{	
		sp->r[i]=0.0;
	}
	sp->nucleo[s->Ato][0]=0;
	sp->nucleo[s->Ato][1]=0;
}

void print_spring_3d(SIMDAT *s, DATASIM *d, SPRING *sp)
{
	// here we collect and print shit.
	// for now, we print ringstats sum and ring count.
	// later we'll want spatial distribution of rings and some such.
	FILE *fp;	
	int t,i, Ato_accum,spam;
	double Ato_avg;
	char buffer[300];
	
	sprintf(buffer,"%satomavg_ring_3d-%d.dat",s->output,sp->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=1;t<=sp->nucleo[0][0];t++)
	{
		// for each atom, what is the average ring size?
		Ato_avg = 0.0;
		Ato_accum=0;
		spam=sp->nucleo[t][0];
		for(i=0;i<20;i++) /// this 20 may have to be replaced
		{	
			//printf("%d \n",sp->ringstats[spam][i]);
			printf("%d %d %d spam\n",spam,sp->nucleo[0][0],s->Ato);
			Ato_avg+=i*sp->ringstats[spam][i]; // weighted averaging
			Ato_accum+=sp->ringstats[spam][i]; // accumulate total number of rings 
		}
		//printf("%d %d %f %d\n",t,spam,Ato_avg,Ato_accum);
		if(Ato_accum>0)
		Ato_avg/=1.0*Ato_accum;
		fprintf(fp,"%d %f\n",t,Ato_avg); // make histogram of this series
	}
	
	
	fclose(fp);
	
	// histogram of ring lengths (total)
	
	sprintf(buffer,"%sringstats_3d-%d.dat",s->output,sp->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=0;t<20;t++)
	{
		fprintf(fp,"%d %d\n",t,sp->ringcount[t]);
	}
	fclose(fp);
	

	
	printf("done printing SP RINGGGGG!!!!\n");
	sp->fileCounter++;
}

void complete_ring(SIMDAT *s, DATASIM *d, SPRING *sp, int m)
{
	int i;
	// this subroutine updates stats, both vertex wise and for the spatial density thing
	printf("%d %d\n",sp->stacklength,m);
	printf("ENTERED RING COMPLETION!!!!!!!\n");
	for(i=0;i<sp->stacklength;i++)
	{
		sp->ringstats[sp->pathstack[i]][m]++; // one more ring of size m for corresponding atom
	}	
	sp->ringcount[m]++;
 	// this will do for now.
}

void compute_spring_3d(SIMDAT *s, DATASIM *d, SPRING *sp)
{
	int i,j,k,neighFlag,pdf,pdl,m,spTet,spN,rflag,rows,cols;
	int vert,si,sj,start,current,sp1,nBond;
	double six,siy,siz,sjx,sjy,sjz,sdx,sdy,sdz,dist,diam,q3q3R,q3q3I,cosTheta,phi;
	
	int *pdflag, *prev, *pdqueue;
   //	need to allocate and initialise here
	// find AllNeigh
	nBond=0;
	for(i=0;i<s->Ato-1;i++)
	{
		for(j=i+1;j<s->Ato;j++)
		{
			neighFlag = 0; // some legacy shit - ignore and don't get triggered
			if(neighFlag==0)
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
			    //if(i>10 && i<20)
			    //printf("wtf distances are we seeing?!!%d %d %f\n",i,j,dist);
			    if(dist <=s->Cl_cut) // we use the cluster formation cut-off here
			    {
			    //printf("distances and participants?? %f %d %d\n",dist,si,sj);
			    //printf("found a neighbour?? %d %d %f\n",si,sj,dist);
   				sp->AllNeigh[si][0]++; // one more neigh for i
				sp->AllNeigh[sj][0]++; // one more neigh for j
				sp->AllNeigh[si][sp->AllNeigh[si][0]]=sj; // one more neighbour for i
				sp->AllNeigh[sj][sp->AllNeigh[sj][0]]=si; // one more neighbour for j
				//printf("how many neigh?? %d %d %d %d\n",si,gr->AllNeigh[si][0],sj,gr->AllNeigh[sj][0]);
			    }
			}			 	
		}
	}
	/////////////////////////////////////////////////////////////
	// first section, find the types of atoms and populate nucleo
	/////////////////////////////////////////////////////////////
	
	
	for(i=0;i<s->Ato;i++)
	{	
	//printf("sanity 1 cluster size %d\n",i);
		nBond=0;
		for(j=0;j<s->Ato;j++)
		{
			
			if(i!=j)
			{
				for(k=0;k<3;k++)
				{
					sp->r[k]=d->r[j][k]-d->r[i][k];
					sp->r[k]-=d->lBox[k]*rint(sp->r[k]/d->lBox[k]);
				}
				sp->r[3]=sqrt(sqr(sp->r[0])+sqr(sp->r[1])+sqr(sp->r[2]));
				if(sp->r[3]<s->Cl_r)			// Atomo dentro la prima shell
				{
					nBond++;
				
					cosTheta=sp->r[2]/sp->r[3];
					phi=atan2(sp->r[1],sp->r[0]);
					
					for(m=0;m<=3;m++)
					{
						sp->q3mR[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(m*phi);
						sp->q3mI[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(m*phi);
					}
					for(m=1;m<=3;m++)
					{
						sp->q3mR[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(-m*phi)*pow(-1,m);
						sp->q3mI[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(-m*phi)*pow(-1,m);
					}
				}
			}
		}
			
		for(m=0;m<7;m++)
		{
			sp->q3mR[i][m]/=nBond;
			sp->q3mI[i][m]/=nBond;
			sp->q3[i]+=sp->q3mR[i][m]*sp->q3mR[i][m]+sp->q3mI[i][m]*sp->q3mI[i][m];
		}
		sp->q3[i]=sqrt(sp->q3[i]*4.0*M_PI/7.0);
		sp->neigh[i]=nBond;								// Adesso ho il q3(i) ed il numero di bond
	}
	
	printf("Finiti q3 e bond\n");
	
	spTet=0.0;
	spN=0.0;
	sp1=0;
	for(i=0;i<s->Ato;i++)
	{				
		sp1=0;			
		for(j=0;j<s->Ato;j++)
		{
			if(i!=j)
			{
				for(k=0;k<3;k++)
				{
					sp->r[k]=d->r[j][k]-d->r[i][k];
					sp->r[k]-=d->lBox[k]*rint(sp->r[k]/d->lBox[k]);
				}
				sp->r[3]=sqrt(sqr(sp->r[0])+sqr(sp->r[1])+sqr(sp->r[2]));
				if(sp->r[3]<s->Cl_r)			// Atomo dentro la prima shell
				{			
					q3q3R=0.0;
					q3q3I=0.0;
					
					for(m=0;m<7;m++)
					{
						q3q3R+=(sp->q3mR[i][m]*sp->q3mR[j][m])+(sp->q3mI[i][m]*sp->q3mI[j][m]);
						q3q3I+=(sp->q3mI[i][m]*sp->q3mR[j][m])-(sp->q3mR[i][m]*sp->q3mI[j][m]);
					}
	
					if(q3q3R<-0.23) // bonded atoms should be in neighList.
					{
						sp1++;
						
					}
				}
			}
		}
		
		if(sp1>=1)
		{
			spTet++;
			sp->nucleo[spTet][0]=i;
		}
		else
		{
			spN++;
			sp->nucleo[spN][1]=i;
		}
	}
	
	sp->nucleo[0][0]=spTet;
	sp->nucleo[0][1]=spN;


	printf("number of atomssss %d %d\n",spTet,spN);
	
	
	
	////////////////////////////
	// populate neighList
	///////////////////////////
	rows = s->Ato; // we still need this over-representing matrix to avoid reference problems
	cols = 10;

	prev=AllocVecI(rows);
	pdflag=AllocVecI(rows);
	pdqueue=AllocVecI(2*rows);
	for(i=0;i<s->Ato;i++)
	{
		prev[i]=0;
		pdflag[i]=0;
		pdqueue[i]=0;
		pdqueue[rows+i]=0;

	}	
	//printf("Checking allocation %d %d %d %d\n",prev[s->Ato-1],pdflag[s->Ato-1],pdqueue[s->Ato],pdqueue[2*s->Ato-1]);
	for(i=1;i<=sp->nucleo[0][0]-1;i++) // for each crystalline atom
	{
		for(j=i+1;j<=sp->nucleo[0][0];j++)
		{
			neighFlag = 0; // same legacy shit
			if(neighFlag==0)
			{
				si = sp->nucleo[i][0]; // atom id of i
				sj = sp->nucleo[j][0]; // atom id of j
			    six = d->r[si][0];
			    siy = d->r[si][1];
			    siz = d->r[si][2];
			    sjx = d->r[sj][0];
			    sjy = d->r[sj][1];
			    sjz = d->r[sj][2];
			    sdx = six - sjx;
			    sdx-= d->lBox[0]*rint(sdx/d->lBox[0]);
			    sdy = siy - sjy;
			    sdy-= d->lBox[1]*rint(sdy/d->lBox[1]);	
			    sdz = siz - sjz;
			    sdz-= d->lBox[2]*rint(sdz/d->lBox[2]);	 		  
			    dist = sqrt(sqr(sdx) + sqr(sdy) + sqr(sdz));		
			    if(dist <=s->Cl_cut) // we use the cluster formation cut-off here
			    {
			    //printf("distance check neighlist %d %d %f\n",si,sj,dist);
   				sp->neighList[si][0]++; // one more neigh for i
				sp->neighList[sj][0]++; // one more neigh for j
				sp->neighList[si][sp->neighList[si][0]]=sj; // one more neighbour for i
				sp->neighList[sj][sp->neighList[sj][0]]=si; // one more neighbour for j
				}
			}	
		}
	}
	//printf("%d\n",sp->neighList[466][0]);
	
	
	////////////////////////////
	// find pairDist matrix using breadth-dirst search with FIFO stack
	////////////////////////////
	

	for(i=1;i<=sp->nucleo[0][0];i++)
	{
		start=sp->nucleo[i][0];
		for(j=1;j<=sp->nucleo[0][0];j++)
		{
		pdflag[j]=0;
		prev[j]=-1;
		}
		pdf=0;
		pdl=0;
		
		// put the start in the queue
		//printf("start %d\n",start);
		pdflag[start]=1;
		pdqueue[pdl]=start;
		pdl++;
		while(pdf<=pdl)
		{
			vert = pdqueue[pdf]; // last item popped
			pdf++;
			for(j=1;j<=sp->neighList[vert][0];j++)
			{
				// for each neighbour
				sj = sp->neighList[vert][j];
				//printf("BFS start vert neighvert pdl j sj %d %d %d %d %d %d\n",start,vert,sp->neighList[vert][0],pdl,j,sj);
				if(pdflag[sj]==0)
				{
					//printf("for the sake of sanity BFS\n");
					pdflag[sj]=1;
					prev[sj]=vert;
					sp->pairDist[start][sj]=sp->pairDist[start][vert]+1;
					sp->pairDist[sj][start]=sp->pairDist[start][sj]; // redundant line
					//printf("did we update any distances at all?? %d %d %d %d \n",vert,sj,sp->pairDist[start][sj],sp->pairDist[start][vert]);
					pdqueue[pdl]=sj;
					pdl++;
				}
			}
			pdflag[vert]=2; // this vertex complete - no use for this just yet
		}		
	}
	//printf("distance pairs %d %d\n",sp->pairDist[3][7],sp->pairDist[7][3]);
	//printf("distance pairs %d %d\n",sp->pairDist[0][2],sp->pairDist[0][4]);
	/////////////////////////////////////////////////
	// now we have pairDist we start our ring-finding
	/////////////////////////////////////////////////
	
	// ordering of action - pick m, pick start (those not in hit list) and go!
	
	for(m=3;m<8;m++)
	{
		// for this ring size, pick start locations from among those that aren't done yet
		// serial will do
		//printf("\n\nthis m %d\n\n",m);
		for(i=0;i<s->Ato;i++)
		{
		sp->hitlist[i]=0; // reset hitlist for next ring size to consider.
		}
		if(m%2==0)
		{
			diam =(int) m/2;
		}
		else
		{
			diam = (int) (m+1)/2;
		}
		for(i=1;i<=sp->nucleo[0][0];i++)
		{
			start = sp->nucleo[i][0];
			current = start;
			sp->stacklength=0;
			sp->pathstack[sp->stacklength]=start;
			sp->stacklength++;
			sp->inpath[start]=1;
			while(sp->stacklength>0)
			{
				// 3 cases, ring completion, less than diameter, 
				if(sp->stacklength == m)
				{
					printf("\nwe here at almost complete rings??\n");
					if(sp->pairDist[current][start]==1)
					{
						complete_ring(s,d,sp,m); // updates stats for all vertices in the ring
						sp->complete[current]=1;
						current=sp->pathstack[sp->stacklength-1]; // current moves one step back
						sp->stacklength--; // finished moving stack one step back
					}					
					else
					{
						// check if we ever really get here
						sp->complete[current]=1;
						current=sp->pathstack[sp->stacklength-1]; // current moves one step back
						sp->stacklength--; // finished moving stack one step back
					}
				}
				else if(sp->stacklength < diam)
				{
					// find a neighbour and check if it's appropriate
					neighFlag=0;
					//printf("1 - increasing ring length %d\n",sp->stacklength);
					for(j=1;j<=sp->neighList[current][0];j++)
					{
						sj = sp->neighList[current][j];
						//printf("2 - current and neigh %d %d\n",current,sj);
						if(sp->hitlist[sj]==1)
						{
							continue;
						}
						//printf("3 - current and neigh innards %d %d\n",sp->pairDist[current][start],sp->pairDist[sj][start]);
						if(sp->pairDist[sj][start] == sp->pairDist[current][start]+1 && sp->complete[sj]==0 && sp->inpath[sj]==0)
						{
							neighFlag=1;
							break;
						}
					}
					if(neighFlag==1)
					{

						sp->pathstack[sp->stacklength]=sj; // add to stack
						sp->stacklength++; 
						//printf("4a - if neighFlag %d\n",sp->stacklength);
						sp->inpath[sj]=1;
						current = sj;
					}
					else
					{

						// no viable neighbours from current
						sp->complete[current]=1;
						sp->inpath[current]=0; // no longer in the path stack
						if(sp->stacklength>0)
						{
						sp->stacklength--;
						current = sp->pathstack[sp->stacklength-1];
						//printf("friendless dude so this line ends before it begins\n")					;
						}
						

						//printf("4b  - else neighflag %d\n",sp->stacklength);						
					}
				}
				else if(sp->stacklength==diam)
				{
					if(m%2==0) // if even, find a neighbour further away or go back
					{
					
					neighFlag=0;
					//printf("*1 - increasing ring length %d\n",sp->stacklength);
					for(j=1;j<=sp->neighList[current][0];j++)
					{
						sj = sp->neighList[current][j];
						//printf("*2 - current and neigh %d %d\n",current,sj);
						if(sp->hitlist[sj]==1)
						{
							continue;
						}
						//printf("*3 - current and neigh innards %d %d\n",sp->pairDist[current][start],sp->pairDist[sj][start]);
						if(sp->pairDist[sj][start] == sp->pairDist[current][start]+1 && sp->complete[sj]==0 && sp->inpath[sj]==0)
						{
							neighFlag=1;
							break;
						}
					}
					if(neighFlag==1)
					{
						sp->pathstack[sp->stacklength]=sj; // add to stack
						sp->stacklength++; 
						//printf("*4a - if neighFlag %d\n",sp->stacklength);
						sp->inpath[sj]=1;
						current = sj;
					}
					else
					{
						// no viable neighbours from current
						sp->complete[current]=1;
						sp->inpath[current]=0; // no longer in the path stack
						sp->stacklength--;
						current = sp->pathstack[sp->stacklength-1];

						//printf("*4b  - else neighflag %d\n",sp->stacklength);												
					}
					}
					
					
					else // if odd, find a neighbour equidistant from start
					{
					
					neighFlag=0;
					for(j=1;j<=sp->neighList[current][0];j++)
					{
						sj = sp->neighList[current][j];
						if(sp->hitlist[sj]==1 || (sp->pairDist[sp->pathstack[sp->stacklength-2]][sj]==1 &&m>3))
						{
							// special little condition for triangle formation in odd-m rings at diam distance from start
							continue;
						}
						if(sp->pairDist[sj][start] == sp->pairDist[current][start] && sp->complete[sj]==0 && sp->inpath[sj]==0)
						{
							//printf("\n(*)(*)we get here??\n");
							neighFlag=1;
							break;
						}
					}
					if(neighFlag==1)
					{
						//printf("\n ODD FOUND NEIGH!\n");
						sp->pathstack[sp->stacklength]=sj; // add to stack
						sp->stacklength++; 
						sp->inpath[sj]=1;
						current = sj;
					}
					else
					{
						//printf("BAD entry\n");
						// no viable neighbours from current
						sp->complete[current]=1;
						sp->inpath[current]=0; // no longer in the path stack
						sp->stacklength--;
						current = sp->pathstack[sp->stacklength-1];

					}
					}
				}
				else if(sp->stacklength >diam)
				{
					//printf("worse entry!!\n");
					neighFlag=0;
					for(j=1;j<=sp->neighList[current][0];j++)
					{
						sj = sp->neighList[current][j];
						if(sp->hitlist[sj]==1 || (sp->pairDist[sp->pathstack[sp->stacklength-2]][sj]==1 &&m>3))
						{
							// special little condition for triangle formation in odd-m rings at diam distance from start
							continue;
						}
						if(sp->pairDist[sj][start] == sp->pairDist[current][start]-1 && sp->complete[sj]==0 && sp->inpath[sj]==0)
						{
							//printf("\n(*)(*)we get here??\n");
							neighFlag=1;
							break;
						}
					}
					if(neighFlag==1)
					{
						sp->pathstack[sp->stacklength]=sj; // add to stack
						sp->stacklength++; 
						sp->inpath[sj]=1;
						current = sj;
					}
					else
					{
						//printf("worse but not worse\n");
						// no viable neighbours from current
						sp->complete[current]=1;
						sp->inpath[current]=0; // no longer in the path stack
						sp->stacklength--;
						current = sp->pathstack[sp->stacklength-1];

					}
				}
			}
			// end while - to complete one start location - all rings of size m with this vertex are enumerated.
			// we need to delete start from the set
			sp->hitlist[start]=1; // start will now be skipped from among available neighbours			
			// we then need to reset complete, inpath, stacklength and pathstack (check if last 3 auto fixed first)
			for(j=0;j<s->Ato;j++)
			{
				sp->complete[j]=0;
				sp->inpath[j]=0;
			}
			printf("%d the stack length \n",sp->stacklength);
			sp->stacklength=0;
		}
		// all rings of size m enumerated
	}
	sp->counter++;
	printf("%d\n",sp->counter);
	if(sp->counter==s->Ring_limit)
	{
		printf("entering print sp RINGGGG!!\n");
		print_spring_3d(s,d,sp);
	}
	free(pdflag);
	free(prev);
	free(pdqueue);
	//for(i=0;i<rows;i++)
	//{
	//	free(neighList[i]);
	//}
	
}


