// This code finds interfaces between crystallites that are close but do not necessarily fit on the same grid.
struct GRAIN
{
	int stack_last, cl_counter, in_clusters,total_step,fileCounter,dev_count,counter;
	int *neigh, *path_stack, *cl_index, *visited, *complete, *boundary;
	int **AllNeigh, **nucleo;
	double *q3, *dev_histo, *r;
	double **q3mR, **q3mI,  **q3q3M, **frame;
	double q3q3R,q3q3I;
} grain;

void Init_Interface(SIMDAT *s, DATASIM *d, GRAIN *gr)
{
	int i,j;
 	//int	**neighlist
 	// double **q3q3, **template
 	//int *visited, *complete, *clIndexlist, *border, *path
 	// double *q3
 	gr->stack_last=0;
 	gr->cl_counter=0;
 	gr->in_clusters=0;
 	gr->total_step=0;
 	gr->fileCounter=0;
 	gr->dev_count=0;
 	gr->q3q3R=0.0;
 	gr->q3q3I=0.0;
 	gr->q3mR=AllocMatR(s->Ato,7);
	gr->q3mI=AllocMatR(s->Ato,7);
	gr->q3=AllocVecR(s->Ato);
	gr->r=AllocVecR(4);
	gr->neigh=AllocVecI(s->Ato); // number of neighbours for each atom
	gr->path_stack=AllocVecI(s->Ato*s->Ato*s->Ato); // max length of the path stack
	gr->cl_index=AllocVecI(s->Ato); // number of crystals can be max number of atoms
	gr->visited=AllocVecI(s->Ato); // flag for each atom
	gr->complete=AllocVecI(s->Ato); // flag for each atom
	gr->boundary=AllocVecI(s->Ato); // flag for each atom
	gr->nucleo=AllocMatI(s->Ato+1,3);	//+1 because the first entry for each of the 3 types is how many of that type
	gr->AllNeigh=AllocMatI(s->Ato,10); // here we have assumed that each atom has max 10 neighbours
 	gr->q3q3M=AllocMatR(s->Ato,s->Ato); // pairwise q3q3 for each pair of particles
 	gr->frame=AllocMatR(4,3);
 	gr->dev_histo=AllocVecR(s->Ato*s->Ato*s->Ato); // some absurd long length
 	// some conflict of interest scenes here since we've assumed max no of neighbours for an atom is 10
 	// but the number of neighbours used to make and compare templates is 4. Check effects later.
 	
 	
 	
 	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<7;j++)
		{
			gr->q3mR[i][j]=0.0;
			gr->q3mI[i][j]=0.0;
		}
		gr->q3[i]=0.0;
		gr->neigh[i]=0;
		gr->cl_index[i]=0;
		gr->visited[i]=0;
		gr->complete[i]=0;
		gr->boundary[i]=0;
		for(j=0;j<s->Ato;j++)
		{
			gr->q3q3M[i][j]=0.0;
		}
		for(j=0;j<10;j++)
		{
			gr->AllNeigh[i][j]=0;
		}
	}
 	for(i=0;i<4;i++)
 	{
 		for(j=0;j<3;j++)
 		{
 			gr->frame[i][j]=0.0;
 		}
 		gr->r[i]=0.0;
 	}
 	for(i=0;i<s->Ato+1;i++)
 	{
 		for(j=0;j<3;j++)
 		{
 			gr->nucleo[i][j]=0;
 		}
 	}
 	for(i=0;i<s->Ato*s->Ato*s->Ato;i++)
 	{
 		gr->dev_histo[i]=0;
 		gr->path_stack[i]=0;
 	}

}

void print_Interface(SIMDAT *s, DATASIM *d, GRAIN *gr)
{

	// first we print out the dev file
	// then we print out the total number of steps
	FILE *fp;	
	int t,si;
	char buffer[300];
	
	
	printf("got the basics - printing interface! %d %d\n",gr->nucleo[0][0],gr->in_clusters); 
	printf("%s s->output\n",s->output);
	sprintf(buffer,"%sinterface_histo-%d.dat",s->output,gr->fileCounter);
  	fp=fopen(buffer,"w");
  	printf("%s\n",buffer);
  	if(fp==NULL)
  	printf("null pointer problems cluster\n");
  	
  	printf("%d dev-count thing %f \n",gr->dev_count,d->lBox[0]);
  	for(t=0;t<gr->dev_count;t++)
  	{
  		//printf("%d - here we go, t\n",t);
  		fprintf(fp,"%d %f\n",t,gr->dev_histo[t]);
  	}
  	fclose(fp);
	sprintf(buffer,"%sq3_histo-%d.dat",s->output,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%f\n",gr->q3[t]);
	}
	fclose(fp);
	sprintf(buffer,"%sneigh_count_histo-%d.dat",s->output,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%d\n",gr->AllNeigh[t][0]);
	}
	fclose(fp);
	
	printf("FINALLY!!! %d %d\n",gr->cl_counter,gr->in_clusters);
	for(t=1;t<=gr->nucleo[0][0];t++)
	{
		si=gr->nucleo[t][0];
		//printf("VnC %d %d %d %d\n",t,si,gr->visited[si],gr->complete[si]);
	}
	printf("\n");
	gr->fileCounter++;
}

double CompareTemplate(SIMDAT *s, DATASIM *d, GRAIN *gr,int current)
{
	// calculate deviation of positions with template positions
	// remember that gr->nucleo[current][0]	is the actual atom index
	// remember that gr->neighList[current][i] is spami
	// further remember that gr->nucleo[spami][0] is the actual neighbour
	// Here we calculate distance of each neighbour from all template locations from gr->frame
	// from current pull neighbours, calculate distance to each [x,y,z] in frame
	// find min distance for each neighbour and find the total deviation
	// return total deviation.
	int cid,natom,aid;
	int k,n_count,neinei;
	double cx,cy,cz,ax,ay,az,dx,dy,dz,min,dist,sq_dev;
	cid = current; // actual atom id is what is passed as argument
	cx = d->r[cid][0];
	cy = d->r[cid][1];
	cz = d->r[cid][2];
	n_count = gr->AllNeigh[current][0]; // number of neighbours
	while(n_count>0)
	{
		neinei = gr->AllNeigh[current][0] - n_count; // which neighbour
		natom = gr->AllNeigh[current][neinei];
		aid = natom; // since neighList stores actual atom ids we don't need to transform again
		ax = d->r[aid][0];
		ay = d->r[aid][1];
		az = d->r[aid][2];
		min=d->lBox[0]*d->lBox[0];
		for(k=0;k<4;k++) // 4 atoms form the frame/template
		{
			// first the position of neighbour relative to current with PBC
			// 25/05 edit -> No PBC cus this is an altogether different reference frame.
			dx = ax - cx;
			dx-=d->lBox[k]*rint(dx/d->lBox[k]);
			dy = ay - cy;
			dy-=d->lBox[k]*rint(dy/d->lBox[k]);
			dz = az - cz;
			dz-=d->lBox[k]*rint(dz/d->lBox[k]);
			dist = sqr(dx-gr->frame[k][0]) + sqr(dy-gr->frame[k][1]) + sqr(dz-gr->frame[k][2]);
			if(dist<min)
			{	
				min = dist;	
			}
			// found min deviation of neighbour from frame/template positions
		}
		sq_dev +=min; // add the min deviation to total deviation of neighbours
		n_count--;
	}
	if(gr->AllNeigh[current][0]!=0)
	{
		sq_dev=sq_dev/gr->AllNeigh[current][0]; // normalised against actual number of neighbours
		return sqrt(sq_dev); // return the square root of deviation
	}
	else
	{
		return 0;
	}

}

int GoTo(SIMDAT *s, DATASIM *d, GRAIN *gr,int target,int current)
{
	// this sub does very little. but it updates the stack and new current location
	int newCurrent,i,j,k;
	newCurrent = target;
	printf("track line track current target %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[newCurrent]);
	gr->visited[newCurrent]=1;
	gr->stack_last++; // one more item in stack
	gr->path_stack[gr->stack_last] = current; // source gets added to the stack
	printf("track line TCT GoTo 2 %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[newCurrent]);
	return newCurrent;
}

int PickNeigh(SIMDAT *s, DATASIM *d, GRAIN *gr,int current,int **neighList,int rows, int cols)
{
 	// this sub finds a viable neighbour to go to, updates location and path stack
 	// it calls a subrouting GoTo that does a bit of the work.
	int flag_gone,i,j,k;
	int spami,candidate;
	flag_gone = 0;
	for(i=1;i<=neighList[current][0];i++) // for each neighbour of i
	{
		spami = neighList[current][i];
		printf("neighList check %d %d %d\n",i,neighList[current][0],neighList[current][i]);
		printf("track line TCT 3 %d %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[spami],spami);
		if(gr->visited[spami]==0) // if unvisited, go there 
		{	
			printf("visited of current %d %d \n",current,gr->visited[current]);
			current = GoTo(s,d,gr,spami,current); // subroutine GoTo updates the stack
			printf("visited of new current %d %d \n",current,gr->visited[current]);
			flag_gone=1;
			break;
		}
			printf("track line TCT 4 %d %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[spami],spami);		
	}

	// case in "if" below is contingent - we don't expect to find ourselves here
	if(flag_gone==0) // if all neighbours already visited
	{	
		for(i=1;i<=neighList[current][0];i++)
		{
			spami = neighList[current][i];
			printf("track line TCT 5 %d %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[spami],spami);
			if(gr->complete[spami]==0 && gr->path_stack[gr->stack_last]!=spami) 
			{
				// not complete and didn't come from here
				if(gr->visited[spami]>1)
				{
				printf("some kind of sick pickneigh joke second if? %d %d %d %d\n",current,spami,gr->visited[current],gr->visited[spami]);
				printf("more laughs, complete, stacked %d %d \n",gr->complete[spami],gr->path_stack[gr->stack_last-1]);
				printf("black magic woman %d\n",gr->stack_last);
				}				

				current = GoTo(s,d,gr,spami,current); //
				flag_gone=1;
				break;
			}
		printf("track line TCT 6 %d %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[spami],spami);		
		}
	}
	if(flag_gone==0) // nowhere to go but back
	{
		gr->complete[current]=1;
		gr->cl_index[current]=gr->cl_counter; // current belongs to the "cl_counter" cluster now
		gr->in_clusters++; // one more atom belongs to a cluster
		//printf("old visited and complete %d %d \n",gr->visited[current],gr->complete[current]);
		candidate = current;
		current = gr->path_stack[gr->stack_last];
		printf("track line TCT 7 %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[candidate]);
		//printf("new and old v and c %d %d %d %d \n",gr->visited[current],gr->complete[current],gr->visited[candidate],gr->complete[candidate]);
		gr->stack_last = gr->stack_last -1; // last item popped out of stack
	}
	printf("track line TCT 8 %d %d %d %d\n",gr->visited[1],gr->visited[current],gr->visited[spami],spami);	
 return current;
}

void compute_Interface(SIMDAT *s, DATASIM *d, GRAIN *gr)
{
 FILE *fp;	
 int t;
 char buffer[300];
 // vars
 // current, flag_gone, flag_complete, neigh_flag, seed_flag, inClust, ClustIndex, next
 int i,j,k,m,neighFlag,histo_spam;
 int sp1, nBond,spLiq4,spLiq5,spSolid;
 int si,sj,seed_flag,neid,current,seed,lchind,spflag;	
 int *list_checked;
 int **neighList;
 int rows,cols;
 double cosTheta,phi,spam;
 double six,siy,siz,sjx,sjy,sjz,sdx,sdy,sdz,dist,sx,sy,sz,nex,ney,nez,dev;
 // first calculate q3 for all particles.
 // also figure out each particle's type
 // take from ClSize.h modify the bond part and segregation part
 
 
 for(i=0;i<s->Ato;i++)
	{	
		nBond = 0;
		gr->q3[i]=0.0;
		for(m=0;m<7;m++)
		{
			gr->q3mR[i][m]=0.0;
			gr->q3mI[i][m]=0.0;
			
		}
		for(j=0;j<s->Ato;j++)
		{			
			if(i!=j)
			{
				for(k=0;k<3;k++)
				{
					gr->r[k]=d->r[j][k]-d->r[i][k];
					gr->r[k]-=d->lBox[k]*rint(gr->r[k]/d->lBox[k]);
				}
				gr->r[3]=sqrt(sqr(gr->r[0])+sqr(gr->r[1])+sqr(gr->r[2]));
				if(gr->r[3]<s->Cl_r)			// Atomo dentro la prima shell
				{	
					nBond++;
					cosTheta=gr->r[2]/gr->r[3];
					phi=atan2(gr->r[1],gr->r[0]);
					
					for(m=0;m<=3;m++)
					{
						gr->q3mR[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(m*phi);
						gr->q3mI[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(m*phi);
					}
					for(m=1;m<=3;m++)
					{
						gr->q3mR[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(-m*phi)*pow(-1,m);
						gr->q3mI[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(-m*phi)*pow(-1,m);
					}
				}
			}
		}
			
		for(m=0;m<7;m++)
		{
			gr->q3mR[i][m]/=nBond;
			gr->q3mI[i][m]/=nBond;
			gr->q3[i]+=gr->q3mR[i][m]*gr->q3mR[i][m]+gr->q3mI[i][m]*gr->q3mI[i][m];
		}
		gr->q3[i]=sqrt(gr->q3[i]*4.0*M_PI/7.0);
		gr->neigh[i]=nBond;		
	}
 	// now update q3q3 and figure out all atoms types;
 	
 	printf("found q3 and number of bonds total\n");
	
	spLiq5=spLiq4=spSolid=0.0;
	for(i=0;i<s->Ato;i++)
	{				
		sp1=0;			
		for(j=0;j<s->Ato;j++)
		{
			gr->q3q3R=0.0;
			gr->q3q3I=0.0;
			if(i!=j)
			{
				for(k=0;k<3;k++)
				{
					gr->r[k]=d->r[j][k]-d->r[i][k];
					gr->r[k]-=d->lBox[k]*rint(gr->r[k]/d->lBox[k]);
				}
				gr->r[3]=sqrt(sqr(gr->r[0])+sqr(gr->r[1])+sqr(gr->r[2]));
				if(gr->r[3]<s->Cl_r)			// Atom in first shell
				{			
					for(m=0;m<7;m++)
					{
						gr->q3q3R+=(gr->q3mR[i][m]*gr->q3mR[j][m])+(gr->q3mI[i][m]*gr->q3mI[j][m]);
						gr->q3q3I+=(gr->q3mI[i][m]*gr->q3mR[j][m])-(gr->q3mR[i][m]*gr->q3mI[j][m]);
					}
					gr->q3q3M[i][j]=gr->q3q3R;
					if(gr->q3q3R<-0.23)
					{
						sp1++;
					}
				}
			}
		}
		
		if(sp1>=s->Cl_solid) // 3 or more bonds
		{	
			spSolid++;
			gr->nucleo[spSolid][0]=i;
			// set current type of i here
		}
		else
		{
			if(gr->q3[i]>=0.6 && gr->neigh[i]==4)
			{
				spLiq4++;
				gr->nucleo[spLiq4][2]=i; 
				// set current type of i here
			}
			else
			{
				spLiq5++;
				gr->nucleo[spLiq5][1]=i; 
				// set current type of i here
			}
		}
	}
	// we have all the types and the q3q3[i,j]
		
 	
 	gr->nucleo[0][0]=spSolid;
 	gr->nucleo[0][1]=spLiq4;
 	gr->nucleo[0][2]=spLiq5;
 	// here we add to neighbourlists of i and j 
	// then we calculate q3[i]
	// we have another loop after these two to calculate q3q3[i,j]
	rows = s->Ato; // we still need this over-representing matrix to avoid reference problems
	cols = 10;
	neighList=AllocMatI(rows,cols);

	list_checked=AllocVecI(gr->nucleo[0][0]);
	for(i=0;i<=gr->nucleo[0][0];i++)
	{
		list_checked[i]=0;
	}
	for(i=0;i<s->Ato;i++)
	{	
		for(j=0;j<10;j++)
		{
			neighList[i][j]=0;
		}
	}
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
   				gr->AllNeigh[si][0]++; // one more neigh for i
				gr->AllNeigh[sj][0]++; // one more neigh for j
				gr->AllNeigh[si][gr->AllNeigh[si][0]]=sj; // one more neighbour for i
				gr->AllNeigh[sj][gr->AllNeigh[sj][0]]=si; // one more neighbour for j
			    }
			}			 	
		}
	}
	for(i=1;i<=gr->nucleo[0][0]-1;i++) // for each crystalline atom
	{
		for(j=i+1;j<=gr->nucleo[0][0];j++)
		{
			neighFlag = 0; // same legacy shit
			if(neighFlag==0)
			{
				si = gr->nucleo[i][0]; // atom id of i
				sj = gr->nucleo[j][0]; // atom id of j
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
   				neighList[si][0]++; // one more neigh for i
				neighList[sj][0]++; // one more neigh for j
				neighList[si][neighList[si][0]]=sj; // one more neighbour for i
				neighList[sj][neighList[sj][0]]=si; // one more neighbour for j
			    }	
			}
			
		}
	}
	// remember now that when we loop over neighbours we will have to
	// -> note that the indices on the neighbourlist give in turn indices for nucleo
	// -> with indices from nucleo we can print final locations from actual atom ids.
 	// now start master loop -> re-check in case any prep is needed
 	histo_spam = 0;
 	printf("how many crystalline?? %d\n",gr->nucleo[0][0]);
 	sprintf(buffer,"%sstack_length_track-%d.dat",s->output,gr->fileCounter);
	fp=fopen(buffer,"w");
 	while(gr->in_clusters<=gr->nucleo[0][0])
 	{
 		seed_flag = 0;
 		lchind = 0;
 		while(seed_flag==0 && lchind<=gr->nucleo[0][0])
 		{
 			seed = (int) gr->nucleo[0][0]*drand48();
 			seed = gr->nucleo[seed][0];
 			spflag = 0;
 			gr->stack_last=0;
 			gr->path_stack[gr->stack_last]=seed;
 			for(i=0;i<=lchind;i++)
 			{
 				if(seed==list_checked[i]) // we've already tried this seed
 				{
 					spflag=1;
 					break;
 				}
 			}
 			if(gr->visited[seed]==0 && gr->q3[seed]>=0.6 && spflag==0)// unvisited, unchecked crystal atom
 			{
 				gr->visited[seed]=1;
 				printf("updated visited seed %d %d\n",seed,gr->visited[seed]);
 				seed_flag=1;
 				break;
 			}
 			else
 			{
 				list_checked[lchind]=seed; // this one's been checked
 				lchind++; // one more rejected
 			}
 		}
 		// if seed_flag still 0 here then all clusters visited
 		if(seed_flag==0)
 		{
 			printf("Houston we have a problem finding seeds\n");
 			break;
 		}
 		// assuming we found a seed. 
 		// first we find the template taking 4 nearest neighbours from neighbourlist

 		for(i=1;i<5;i++) // zeroeth of neighList is number of neighbours/
 		{
 			// find xyz of seed
 			neid = gr->AllNeigh[seed][i]; // neid is the ith neighbour of the seed
 			sx = d->r[seed][0];
 			sy = d->r[seed][1];
 			sz = d->r[seed][2];
 			nex = d->r[neid][0];
 			ney = d->r[neid][1];
 			nez = d->r[neid][2];
 			gr->frame[i-1][0] = nex - sx;
 			gr->frame[i-1][0]-= d->lBox[0]*rint(gr->frame[i-1][0]/d->lBox[0]);
 			gr->frame[i-1][1] = ney - sy;
 			gr->frame[i-1][1]-= d->lBox[1]*rint(gr->frame[i-1][1]/d->lBox[1]);
 			gr->frame[i-1][2] = nez - sz;
 			gr->frame[i-1][2]-= d->lBox[2]*rint(gr->frame[i-1][2]/d->lBox[2]);			
 		}
 		// found frame of seed. 
 		// now we propagate from seed asking if each successive neighbour and neighbour of neighbour fit with the frame.
 		current = seed;

 		// now we begin!
 		while(gr->complete[seed]==0)
 		{
 			gr->total_step++; // one more visitation!
 			dev = CompareTemplate(s,d,gr,current); // returns square root of sum of squared individual deviations
 			gr->dev_histo[gr->dev_count]=dev; // append the deviation to the histogram
 			gr->dev_count++;
 			if(dev<s->interface_tolerance)
 			{
 				printf("current %d\n",current);
 				for(i=1;i<=neighList[current][0];i++)
 				{
 					printf("neighbour %d ",neighList[current][i]);
 					
 				}
 				printf("\n");
 				for(i=1;i<=neighList[current][0];i++)
 				{
 					printf("neighbour v %d ",gr->visited[neighList[current][i]]);
 				}
 				printf("\n");
 				current=PickNeigh(s,d,gr,current,neighList,rows,cols); // updates position
 				printf("current second set %d\n",current);
 				for(i=1;i<=neighList[current][0];i++)
 				{
 					printf("neighbour second set %d ",neighList[current][i]);
 					
 				}
 				printf("\n");
 				for(i=1;i<=neighList[current][0];i++)
 				{
 					printf("neighbour v second set %d ",gr->visited[neighList[current][i]]);
 				}
 				printf("\n");
 			}
 			if(dev>=s->interface_tolerance)// nowhere to go but back
			{
				printf("here? %d %d\n",current,gr->visited[current]);
				
				gr->complete[current]=1;
				gr->boundary[current]=1; // this is a boundary atom
				gr->cl_index[current]=gr->cl_counter; // current belongs to the "cl_counter" cluster now
				gr->in_clusters++; // one more atom belongs to a cluster
				current = gr->path_stack[gr->stack_last];
				gr->stack_last = gr->stack_last -1; // last item popped out of stack
			}
			for(i=1;i<=neighList[current][0];i++)
 			{
 				if(gr->visited[neighList[current][i]]>1)
 				{
 				printf("exiting");
 				exit(0);
 				}

 			}	
 		}
 			
 	}
 	fclose(fp);
 	
	gr->counter++;
	printf("%d\n",gr->counter);
	if(gr->counter==s->Interface_limit)
	{
		printf("entering print gr\n");
		print_Interface(s,d,gr);
	}
//	free(list_checked);
//	for(i=0;i<s->Ato;i++)
//	{
//		free(neighList[i]);
//	}
//	free(neighList);

}
