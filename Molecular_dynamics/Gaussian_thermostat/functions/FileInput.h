void readPara(FILE *fp, struct PARAINPUT *p)
{
	char string[100];

	if(fscanf(fp, "%s",string)>0); 				// nAto
	if(fscanf(fp, "%s",string)>0);
	p->nAto=atoi(string);
    if(fscanf(fp, "%s",string)>0);				// Press
	if(fscanf(fp, "%s",string)>0); 	
	p->ens=atoi(string);
	if(fscanf(fp, "%s",string)>0);				// Press
	if(fscanf(fp, "%s",string)>0); 	
	p->eps=atof(string);
    if(fscanf(fp, "%s",string)>0);				// Press
	if(fscanf(fp, "%s",string)>0); 	
	p->sigma=atof(string);
    if(fscanf(fp, "%s",string)>0);				// Press
	if(fscanf(fp, "%s",string)>0); 	
	p->rc=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// DispMax
	if(fscanf(fp, "%s",string)>0);
	p->dt=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// EqRun
	if(fscanf(fp, "%s",string)>0);
	p->eqRun=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// TotRun
	if(fscanf(fp, "%s",string)>0);
	p->totRun=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// Dump
	if(fscanf(fp, "%s",string)>0);
	p->dump=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// Thermo
	if(fscanf(fp, "%s",string)>0);
	p->thermo=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// Seed
	if(fscanf(fp, "%s",string)>0);
	p->seed=atoi(string);
	if(fscanf(fp, "%s",string)>0);				// InitRead
	if(fscanf(fp, "%s",string)>0);
	p->initRead=atoi(string);
	if(fscanf(fp, "%s",string)>0);				// Label
	if(fscanf(fp, "%s",p->label)>0);
    if(fscanf(fp, "%s",string)>0);				// InitRead
	if(fscanf(fp, "%s",string)>0);
	p->startCry=atoi(string);
    if(fscanf(fp, "%s",string)>0);				// InitRead
	if(fscanf(fp, "%s",string)>0);
	p->nu=atof(string);
}


void readSource(FILE *fp, double **source, char **fileList, int nSyst)
{
	char *line;
	int i,length,not_space;
	int count,vals;
	
	line=(char *) malloc(sizeof(char) * LINESIZE);
	if(!line)
	{
		printf("Read Source: Couldn't allocate space for line\n");
		exit(-1);
	}

	rewind(fp);
	line=fgets(line,LINESIZE,fp);
	count=0;		// Count how many lines are there
	
	while(line!=NULL && count<nSyst)
	{
		if(line[0]=='#')	// It is a comment
		{
			line=fgets(line,LINESIZE,fp);
		}
		else 	
		{
			// verify it's not a blank line
			length=strlen(line);		// string length
			not_space=0;
			i=0;
			while(!not_space && i<length)
			{	
				if(!isspace(line[i]))	// Checks whether line[i] is a white-space character.
				{
					not_space=1;
				}
				i++;
			}
			
			if(not_space)
			{
				vals=sscanf(line,"%lf  %lf %lf",&source[count][0],&source[count][1],&source[count][2]);
				if(vals!=3)
				{
					printf("Source file: reading %d values, but expected 3 values!\n",vals);
					exit(-1);
				}
				count++;
				line=fgets(line,LINESIZE,fp);
			}
			else
			{
				line=fgets(line,LINESIZE,fp);
			}			
		}		
	}
	
	if(count!=nSyst)
	{
		printf("Source Error: expected %d lines but found %d lines\n",nSyst,count);
		exit(-1);
	}
	
	free(line);
}


void writeLAMMPSdump(struct PARAINPUT *p, struct SYSTEM *structbox, struct APPOGGIO *a, int MCstep)
{
	// subroutine to write output in format that mimics LAMMPS
	// this is so that AnaCode can be used on configs printed out of this code as well
	int i;
	fprintf(structbox->fpLAMMPSdump,"ITEM: TIMESTEP\n");
	fprintf(structbox->fpLAMMPSdump,"%d\n",MCstep);		// TIMESTEP	
	fprintf(structbox->fpLAMMPSdump,"ITEM: NUMBER OF ATOMS\n");
	fprintf(structbox->fpLAMMPSdump,"%d\n",p->nAto);
	fprintf(structbox->fpLAMMPSdump,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(structbox->fpLAMMPSdump,"%f %f\n",-0.5*structbox->lx,0.5*structbox->lx);
	fprintf(structbox->fpLAMMPSdump,"%f %f\n",-0.5*structbox->ly,0.5*structbox->ly);
	fprintf(structbox->fpLAMMPSdump,"%f %f\n",-0.5*structbox->lz,0.5*structbox->lz);
	fprintf(structbox->fpLAMMPSdump,"ITEM: ATOMS id xu yu zu vx vy vz\n");
	for(i=0;i<p->nAto;i++)
	{
		fprintf(structbox->fpLAMMPSdump,"%d %lf %lf %lf %lf %lf %lf\n",i+1,structbox->X[i],structbox->Y[i],structbox->Z[i],structbox->vx[i],structbox->vy[i],structbox->vz[i]);
	}
}

void writeRestart(struct PARAINPUT *p, struct SYSTEM *structbox, struct APPOGGIO *a, int MCstep)
{
	int i;

	fprintf(structbox->fprestart,"%d\n",MCstep);
	fprintf(structbox->fprestart,"%d\n",a->mti);

	for(i=0;i<N;i++)
	{
		fprintf(structbox->fprestart,"%lu\n",a->mt[i]);
	}
	fprintf(structbox->fprestart,"%2.14f\n",structbox->lx);
	for(i=1;i<=p->nAto+1;i++)
	{
		fprintf(structbox->fprestart,"%2.14f %2.14f	%2.14f\n",structbox->X[i],structbox->Y[i],structbox->Z[i]);
	}
}

int readRestart(struct PARAINPUT *p, struct SYSTEM *structbox, struct APPOGGIO *a)
{
	int mc,i;
	char string[100];

	if(fscanf(structbox->fprestart,"%s",string)>0);				// MCstep
	mc=atoi(string);
	mc++;
	if(fscanf(structbox->fprestart,"%s",string)>0);				// mti
	a->mti=atoi(string);
	for(i=0;i<N;i++)
	{
		if(fscanf(structbox->fprestart,"%s",string)>0);
		a->mt[i]=strtoul(string,NULL,0);
	}
	if(fscanf(structbox->fprestart,"%s",string)>0);				// lx
	structbox->lx=atof(string);			
	structbox->ly=structbox->lx;								// ly
	structbox->lz=structbox->lx;								// lx
	structbox->rho=p->nAto/(structbox->lx*structbox->ly*structbox->lz);	// rho

	for(i=1;i<=p->nAto+1;i++)					// X,Y,Z
	{
		if(fscanf(structbox->fprestart,"%s",string)>0);
		structbox->X[i]=atof(string);	
		if(fscanf(structbox->fprestart,"%s",string)>0);
		structbox->Y[i]=atof(string);
		if(fscanf(structbox->fprestart,"%s",string)>0);
		structbox->Z[i]=atof(string);
	}
    //printf("%f %f %f\n",structbox->X[p->nAto],structbox->Y[p->nAto],structbox->Z[p->nAto]);
	return(mc);
}



