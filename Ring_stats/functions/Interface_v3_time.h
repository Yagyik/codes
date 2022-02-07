// This code finds interfaces between crystallites that are close but do not necessarily fit on the same grid.
#include<time.h>
struct GRAIN
{
	int stack_last, cl_counter, in_clusters,total_step,fileCounter,dev_count,counter,frameCount;
	int *neigh, *path_stack, *cl_index, *visited, *complete, *boundary, *inPath, *cl_size;
	int **AllNeigh, **nucleo, **bulkbound;
	double *q3, *dev_histo, *r;
	double **q3mR, **q3mI,  **q3q3M, **frame;
	double q3q3R,q3q3I,boundary_sum;
} grain;

void Init_Interface(SIMDAT *s, DATASIM *d, GRAIN *gr)
{
	srand(time(0));
	int i,j;
 	//int	**neighlist
 	// double **q3q3, **template
 	//int *visited, *complete, *clIndexlist, *border, *path
 	// double *q3
 	gr->boundary_sum;
 	gr->stack_last=0;
 	gr->cl_counter=0;
 	gr->in_clusters=0;
 	gr->total_step=0;
 	gr->fileCounter=0;
 	gr->dev_count=0;
 	gr->frameCount=0;
 	gr->q3q3R=0.0;
 	gr->q3q3I=0.0;
 	gr->q3mR=AllocMatR(s->Ato,7);
	gr->q3mI=AllocMatR(s->Ato,7);
	gr->q3=AllocVecR(s->Ato);
	gr->r=AllocVecR(4);
	gr->neigh=AllocVecI(s->Ato); // number of neighbours for each atom
	gr->path_stack=AllocVecI(s->Ato*s->Ato); // max length of the path stack
	gr->cl_index=AllocVecI(s->Ato); // number of crystals can be max number of atoms
	gr->visited=AllocVecI(s->Ato); // flag for each atom
	gr->complete=AllocVecI(s->Ato); // flag for each atom
	gr->boundary=AllocVecI(s->Ato); // flag for each atom
	gr->inPath=AllocVecI(s->Ato);
	gr->cl_size=AllocVecI(s->Ato); // size of each cluster
	gr->nucleo=AllocMatI(s->Ato+1,3);	//+1 because the first entry for each of the 3 types is how many of that type
	gr->AllNeigh=AllocMatI(s->Ato,15); // here we have assumed that each atom has max 10 neighbours
	gr->bulkbound=AllocMatI(s->Interface_limit,2);
 	gr->q3q3M=AllocMatR(s->Ato,s->Ato); // pairwise q3q3 for each pair of particles
 	gr->frame=AllocMatR(4,3);
 	gr->dev_histo=AllocVecR(s->Ato*s->Ato); // some absurd long length
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
		gr->inPath[i]=0;
		gr->cl_size[i]=0;
		gr->bulkbound[i][0]=0;
		gr->bulkbound[i][1]=0;
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
 	for(i=0;i<s->Ato*s->Ato;i++)
 	{
 		gr->dev_histo[i]=0;
 		gr->path_stack[i]=0;
 	}

}

void print_Interface(SIMDAT *s, DATASIM *d, GRAIN *gr)
{

	// first we print out the dev file
	// then we print out the total number of steps
	FILE *fp, *fp2;	
	int t,si,boundary_sum,m;
	char buffer[300],buffer2[300];
	
	
	printf("got the basics - printing interface! %d %d\n",gr->nucleo[0][0],gr->in_clusters); 
	printf("%s s->output\n",s->output);
	sprintf(buffer,"%sinterface_histo-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
  	fp=fopen(buffer,"w");
  	printf("%s\n",buffer);
  	if(fp==NULL)
  	printf("null pointer problems cluster\n");
  	fclose(fp);
	sprintf(buffer,"%sq3_histo-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%f\n",gr->q3[t]);
	}
	fclose(fp);
	sprintf(buffer,"%sneigh_count_histo-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%d\n",gr->AllNeigh[t][0]);
	}
	fclose(fp);
	
	sprintf(buffer,"%sin_Path-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%d\n",gr->inPath[t]);
	}
	fclose(fp);
	
	sprintf(buffer,"%sClSize_interface-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		fprintf(fp,"%d\n",gr->cl_size[t]);
	}
	fclose(fp);
	
	
	sprintf(buffer,"%sDev_histo-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<gr->dev_count;t++)
	{
		fprintf(fp,"%f\n",gr->dev_histo[t]);
	}
	fclose(fp);
	sprintf(buffer,"%sBoundary-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{
		if(gr->boundary[t]==1)
		{
			fprintf(fp,"1.1 %f %f %f\n",d->r[t][0],d->r[t][1],d->r[t][2]);
		}
	}
	fclose(fp);
	
	sprintf(buffer,"%scry-%f-clr-%f-tol-%d.dat",s->output,s->Cl_r,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=1;t<=gr->nucleo[0][0];t++)
	{	
		si=gr->nucleo[t][0];
		fprintf(fp,"cry %f %f %f\n",d->r[si][0],d->r[si][1],d->r[si][2]);
	}
	fclose(fp);
	
	sprintf(buffer,"%shdl-%f-clr-%f-tol-%d.dat",s->output,s->Cl_r,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=1;t<=gr->nucleo[0][1];t++)
	{	
		si=gr->nucleo[t][1];
		fprintf(fp,"hdl %f %f %f\n",d->r[si][0],d->r[si][1],d->r[si][2]);
	}
	fclose(fp);
	
	sprintf(buffer,"%sldl-%f-clr-%f-tol-%d.dat",s->output,s->Cl_r,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=1;t<=gr->nucleo[0][2];t++)
	{	
		si=gr->nucleo[t][2];
		fprintf(fp,"ldl %f %f %f\n",d->r[si][0],d->r[si][1],d->r[si][2]);
	}
	fclose(fp);
	
	sprintf(buffer,"%stype_cry-%f-clr-%f-tol-%d.dat",s->output,s->Cl_r,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Ato;t++)
	{	
		if(gr->boundary[t]==1)
		{
			fprintf(fp,"bound %f %f %f\n",d->r[si][0],d->r[si][1],d->r[si][2]);
		}
		
		else
		{
			fprintf(fp,"bulk %f %f %f\n",d->r[si][0],d->r[si][1],d->r[si][2]);
		}
		
	}
	fclose(fp);

	
	boundary_sum=0;
	for(t=0;t<s->Ato;t++)
	{
		//printf("b sum %d\n",boundary_sum);
		boundary_sum +=gr->boundary[t];
	}
	printf("Number of boundary atoms %f %d\n",s->interface_tolerance,boundary_sum);
	printf("FINALLY!!! %d %d\n",gr->cl_counter,gr->in_clusters);
	for(t=1;t<=gr->nucleo[0][0];t++)
	{
		si=gr->nucleo[t][0];
		//printf("VnC %d %d %d %d\n",t,si,gr->visited[si],gr->complete[si]);
	}
	printf("\n");
	
	for(m=0;m<3;m++)
	{
		sprintf(buffer,"%stypelist_%d-%f-tol-%d.dat",s->output,m,s->interface_tolerance,gr->fileCounter);
		fp=fopen(buffer,"w");
		fprintf(fp,"type%d\n",m);
		fprintf(fp,"%d\n",gr->nucleo[0][m]);
		for(t=1;t<=gr->nucleo[0][m];t++)
		{
			si = gr->nucleo[t][m];
			//printf("%d %d\n",m,si);
			fprintf(fp,"%d\n",si);
		}
		fclose(fp);
	}
	sprintf(buffer,"%scry_boundlist_-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	sprintf(buffer2,"%scry_bulklist_-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp2=fopen(buffer2,"w");
	fprintf(fp2,"type_bulk\n");
	fprintf(fp2,"%d\n",gr->nucleo[0][0]-boundary_sum);
	fprintf(fp,"type_boundary\n");
	fprintf(fp,"%d\n",boundary_sum);
	for(t=1;t<=gr->nucleo[0][0];t++)
	{
		si = gr->nucleo[t][0];
		if(gr->boundary[si]==1)
		fprintf(fp,"%d\n",si);
		else
		fprintf(fp2,"%d\n",si);
	}
	fclose(fp);
	fclose(fp2);
	
	sprintf(buffer,"%sbulkboundvTime_-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<s->Interface_limit;t++)
	{
		fprintf(fp,"%f %d %d\n",s->DeltaT*s->Step*t,gr->bulkbound[t][0],gr->bulkbound[t][1]);
	}
	fclose(fp);
	gr->fileCounter++;
}

double CompareTemplate(SIMDAT *s, DATASIM *d, GRAIN *gr,int current, int spamfl)
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
	sq_dev=0;
	while(n_count>0)
	{
		neinei = gr->AllNeigh[current][0] +1 - n_count; // which neighbour
		natom = gr->AllNeigh[current][neinei];
		aid = natom; // since neighList stores actual atom ids we don't need to transform again
		//printf("%d aid\n",aid);
		ax = d->r[aid][0];
		ay = d->r[aid][1];
		az = d->r[aid][2];
		min=d->lBox[0]*d->lBox[0];
		dx = ax - cx;
		dx-=d->lBox[0]*rint(dx/d->lBox[0]);
		dy = ay - cy;
		dy-=d->lBox[1]*rint(dy/d->lBox[1]);
		dz = az - cz;
		dz-=d->lBox[2]*rint(dz/d->lBox[2]);
		
		
		/// for loop definition below has issues tied to the definition of frame size 
		// refer line {//for(i=1;i<=gr->AllNeigh[seed][0];i++) }
		//for(k=0;k<gr->frameCount;k++) // 4 atoms form the frame/template
		//printf("%d\n",gr->frameCount);
		for(k=0;k<4;k++)
		{
			// first the position of neighbour relative to current with PBC
			// 25/05 edit -> No PBC cus this is an altogether different reference frame.

			dist = sqr(dx-gr->frame[k][0]) + sqr(dy-gr->frame[k][1]) + sqr(dz-gr->frame[k][2]);
			if(dist<=min)
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
	gr->visited[newCurrent]=1;
	gr->stack_last++; // one more item in stack
	gr->path_stack[gr->stack_last] = current; // source gets added to the stack
	gr->inPath[current]=gr->inPath[current] + 1; // current is now in the path stack
	return newCurrent;
}

int PickNeigh(SIMDAT *s, DATASIM *d, GRAIN *gr,int current,int **neighList,int rows, int cols)
{
 	// this sub finds a viable neighbour to go to, updates location and path stack
 	// it calls a subrouting GoTo that does a bit of the work.
	int flag_gone,i,j,k,found_flag;
	int spami,candidate;
	flag_gone = 0;
	found_flag = 0;
	for(i=1;i<=neighList[current][0];i++) // for each neighbour of i
	{
		spami = neighList[current][i];
		if(gr->visited[spami]==0) // if unvisited, go there 
		{	
			//printf("are and going %d %d %d\n",i, current, spami);
			current = GoTo(s,d,gr,spami,current); // subroutine GoTo updates the stack
			flag_gone=1;
			break;
		}
	}
	if(flag_gone==0) // nowhere to go but back
	{
		//if(neighList[current][0]<gr->frameCount) // BOGEY CONDITION #################################
		//gr->boundary[current]=1; 				 // CONSIDER CAREFULLY
		gr->complete[current]=1;
		gr->cl_index[current]=gr->cl_counter; // current belongs to the "cl_counter" cluster now
		gr->in_clusters++; // one more atom belongs to a cluster
		gr->cl_size[gr->cl_counter]++; // one more atom belongs to this particular cluster
		candidate = current;
		current = gr->path_stack[gr->stack_last];
		gr->stack_last = gr->stack_last -1; // last item popped out of stack
		gr->inPath[current]=gr->inPath[current]-1;
		if(gr->inPath[current]<0)
		{
			printf("EXCEPTION EXCEPTION!!!! %d %d\n",current,gr->stack_last);
		}
	}
	//printf(" gone?? %d\n",current);
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
 int si,sj,seed_flag,neid,current,seed,lchind,spflag,found_flag;	
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
				//printf("HDL %d \n",i);
				// set current type of i here
			}
		}
	}
	// we have all the types and the q3q3[i,j]
		

 	gr->nucleo[0][0]=spSolid;
 	gr->nucleo[0][1]=spLiq5;
 	gr->nucleo[0][2]=spLiq4;
 	 	printf("found types?? %d %d %d\n",gr->nucleo[0][0],gr->nucleo[0][1],gr->nucleo[0][2]);
 	// here we add to neighbourlists of i and j 
	// then we calculate q3[i]
	// we have another loop after these two to calculate q3q3[i,j]
	rows = s->Ato; // we still need this over-representing matrix to avoid reference problems
	cols = 15;
	found_flag=0;
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
			    if(dist <=s->Cl_r) // we use the first min g(r) cut-off here
			    {
   				gr->AllNeigh[si][0]++; // one more neigh for i
				gr->AllNeigh[sj][0]++; // one more neigh for j
				gr->AllNeigh[si][gr->AllNeigh[si][0]]=sj; // one more neighbour for i
				gr->AllNeigh[sj][gr->AllNeigh[sj][0]]=si; // one more neighbour for j
			    }
			}			 	
		}
		//printf("how many neigh %d %d\n",si,gr->AllNeigh[si][0]);
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
			    // the cut-off we use here has to be the long cutoff because this is the traversal list
			    // crystalline particles may not be in each other's coordination shells but should still be evaluate
			    if(dist <=s->Cl_cut) // we have been using the cluster formation cut-off here 
			    {
			    	//if(neighList[si][0]<=10 && neighList[sj][0] <=10)
			    	//{
			    	neighList[si][0]++; // one more neigh for i
					neighList[sj][0]++; // one more neigh for j
					if(neighList[si][0] ==0 || neighList[sj][0] ==0)
					printf("SPAM SPAM SPAM!!\n");
					neighList[si][neighList[si][0]]=sj; // one more neighbour for i
					neighList[sj][neighList[sj][0]]=si; // one more neighbour for j
			    	//}

			    }	
			}
			
		}
		//printf("how many neigh %d %d\n",si,neighList[si][0]);
	}
	printf("neigh lists made\n");
	sprintf(buffer,"%straverse_histo-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
	for(i=1;i<=gr->nucleo[0][0];i++)
	{
		si=gr->nucleo[i][0];
		printf("%d\n",si);
		fprintf(fp,"%d %d\n",si,neighList[si][0]);
	}
	fclose(fp);
	// remember now that when we loop over neighbours we will have to
	// -> note that the indices on the neighbourlist give in turn indices for nucleo
	// -> with indices from nucleo we can print final locations from actual atom ids.
 	// now start master loop -> re-check in case any prep is needed
 	histo_spam = 0;
 	printf("how many crystalline?? %d\n",gr->nucleo[0][0]);
 	sprintf(buffer,"%sstack_length_track-%f-tol-%d.dat",s->output,s->interface_tolerance,gr->fileCounter);
	fp=fopen(buffer,"w");
 	while(gr->in_clusters<=gr->nucleo[0][0])
 	{
 		printf("\n \n \nhow far along %d %d\n",gr->in_clusters,gr->nucleo[0][0]);
 		seed_flag = 0;
 		lchind = 0;
 		while(seed_flag==0 && lchind<=gr->nucleo[0][0])
 		{
 			seed = 5 + (int) gr->nucleo[0][0]*drand48();
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
 			if(gr->visited[seed]==0 && gr->q3[seed]>=0.7 && spflag==0)// unvisited, unchecked crystal atom
 			{
 				gr->visited[seed]=1;
 				gr->cl_counter++; // one more cluster!
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
 		printf("found seed %d\n\n",seed);
 		// assuming we found a seed. 
 		// first we find the template taking 4 nearest neighbours from neighbourlist
		for(i=0;i<s->Ato*s->Ato;i++)
 		{
 		gr->path_stack[i]=0; // reset for new seed
 		}
 		for(i=0;i<s->Ato;i++)
 		{
 			gr->inPath[i]=0; // reset for new seed
 		}
 		gr->frameCount=gr->AllNeigh[seed][0];
 		// first version of for loop here is an issue bec the space allocated to frame is for 4 neighbours 
 		// and some crystalline particles, for reasons hard to understand have 5
 		// this is particularly an issue if we choose a 5-coordinated particle as crystal seed
 		//for(i=1;i<=gr->AllNeigh[seed][0];i++) // zeroeth of neighList is number of neighbours/
 		for(i=1;i<=4;i++) // we store only the first 4 particles in cry neighbourhood - could give rise to issues since we don't 
 		// know which are the ordered ones. 
 		{
 			// find xyz of seed
 			neid = gr->AllNeigh[seed][i]; // neid is the ith neighbour of the seed
 			//printf("the neighbour in frame %d %d\n",seed,neid);
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
 			//printf("frame position for it %f %f %f\n",gr->frame[i-1][0],gr->frame[i-1][1],gr->frame[i-1][2]);			
 		}
 		// found frame of seed. 
 		// now we propagate from seed asking if each successive neighbour and neighbour of neighbour fit with the frame.
 		current = seed;
//		for(i=1;i<=gr->AllNeigh[seed][0];i++)
//		{
//			printf("  all neigh seed %d",gr->AllNeigh[seed][i]);
//		}
//		printf("\n");
 		// now we begin!
 		while(gr->complete[seed]==0)
 		{
 			gr->total_step++; // one more visitation!
 			if(current==seed)
 			dev = CompareTemplate(s,d,gr,current,1); // returns square root of sum of squared individual deviations
 			else
 			dev = CompareTemplate(s,d,gr,current,0); // returns square root of sum of squared individual deviations
 			gr->dev_histo[gr->dev_count]=dev; // append the deviation to the histogram
 			gr->dev_count++;
 			if(dev<s->interface_tolerance)
 			{
 				//printf("moving along! %d %d\n",seed,current);
 				current=PickNeigh(s,d,gr,current,neighList,rows,cols); // updates position
 				//printf("moved along %d %d \n",seed,current);
 			}
 			if(dev>=s->interface_tolerance)// nowhere to go but back
			{				
				gr->complete[current]=1;
				//if(neighList[current][0] < gr->frameCount)
				gr->boundary[current]=1; // t// BOGEY CONDITION REPEATED        #################################
				gr->cl_index[current]=gr->cl_counter; // current belongs to the "cl_counter" cluster now
				gr->in_clusters++; // one more atom belongs to a cluster
				gr->cl_size[gr->cl_counter]++; // one more atom belongs to this particular cluster
				current = gr->path_stack[gr->stack_last];
				gr->stack_last = gr->stack_last - 1; // last item popped out of stack
				gr->inPath[current]=gr->inPath[current]-1;
				if(gr->inPath[current]<0)
				{
					printf("prep exception %d %d %d\n",seed,current,gr->inPath[current]);
					printf("EXCEPTION EXCEPTION compute!!!!\n");
				}
			}	
 		}
 			
 	}
 	fclose(fp);
 	gr->boundary_sum=0;
	for(i=0;i<s->Ato;i++)
	{
		//printf("b sum %d\n",boundary_sum);
		gr->boundary_sum +=gr->boundary[t];
	}
	gr->bulkbound[gr->counter][0]=gr->boundary_sum;
	gr->bulkbound[gr->counter][1]=gr->nucleo[0][0]-gr->boundary_sum;
	gr->counter++;
	printf("%d\n",gr->counter);
	if(gr->counter==s->Interface_limit)
	{
		printf("entering print gr\n");
		print_Interface(s,d,gr);
	}
	free(list_checked);
	for(i=0;i<s->Ato;i++)
	{
		free(neighList[i]);
	}
	free(neighList);

}
