void readLammps(FILE *fp,struct DATASIM *d,int t)
{
	int n,i,j;
	char datoSTR[100];
	printf("readLammps - %s\n",d->name);
	if(fp == NULL) 
	{
		printf("%s not found!\n",d->name);
		exit(0);
	}
	
	fscanf(fp, "%s", datoSTR); 		// ITEM:
	//printf("%s, wtf\n",datoSTR);
	fscanf(fp, "%s", datoSTR);		// TIMESTEP
	//printf("%s, wtf\n",datoSTR);
	fscanf(fp, "%s", datoSTR);		
	n=atoi(datoSTR);
	//printf("time actually, not number of atoms %d\n",n);
	if(n!=t)
	{
		printf("Error!\nI am reading file %d instead of %d\n",n,t);
		exit(0);
	}

	for(n=0;n<5;n++)				// ITEM: NUMBER OF ATOMS
	{
		fscanf(fp, "%s", datoSTR); 	
	}
	n=atoi(datoSTR);
	printf("is this number of atoms?? %d\n",n);
	
	if(n!=d->nAto)
	{
		printf("Error!\nI am reading %d atoms instead of %d\n",n,d->nAto);
		exit(0);
	}
	for(n=0;n<6;n++)				// ITEM: BOX BOUNDS pp pp pp
	{
		fscanf(fp, "%s", datoSTR); 
		printf("%s, wtf\n",datoSTR);	
	}
	for(n=0;n<3;n++)
	{
		fscanf(fp, "%s", datoSTR);
		d->lBox[n]=atof(datoSTR);
		fscanf(fp, "%s", datoSTR);
		d->lBox[n]-=atof(datoSTR);
		d->lBox[n]*=-1;
		printf("%f\n",d->lBox[n]);
		//printf("the box dimension of %d is %f\n",n,d->lBox[n]);
	}
	
//	for(n=0;n<10;n++)	// ITEM: ATOMS id type xu yu zu vx vy vz
//	{
//		fscanf(fp, "%s", datoSTR); 	
//	}
//	for(n=0;n<s->Ato;n++)
//	{
//		
//		fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf\n",&i,&j,&d->r[n][0],&d->r[n][1],&d->r[n][2],&d->rv[n][0],&d->rv[n][1],&d->rv[n][2]);
//	}

//******************************************************************************//
// below we find august code written by the august chinna-mari

	for(n=0;n<9;n++) // ITEM: ATOMS id xu yu zu vx vy vz
	{
		fscanf(fp, "%s",datoSTR);
	}
	//printf("nAto %d",s->Ato);
	for(n=0;n<d->nAto;n++) // xu yu and zu
	{	
		fscanf(fp,"%d",&i); // reads the integer id and then the "cursor" moves to the space and 
		//printf("%d - the id\n",i);
		// we get the line below!!! where it reads the rest of the line acc. format below
		
		/////////////////////////
		// format with type info also
		fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&d->x[i-1],&d->y[i-1],&d->z[i-1],&d->vx[i-1],&d->vy[i-1],&d->vz[i-1]);
		/////
		//fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&d->r[i-1][0],&d->r[i-1][1],&d->r[i-1][2],&d->rv[i-1][0],&d->rv[i-1][1],&d->rv[i-1][2]);
		// above is own-code format
		// alternate format - Vishwas binary reconstruct format.
		
		//fscanf(fp,"%lf %lf %lf\n",&d->r[i-1][0],&d->r[i-1][1],&d->r[i-1][2]);
		//printf("%lf %lf %lf\n",d->r[i-1][0],d->r[i-1][1],d->r[i-1][2]);
		//////////////////////
		 
		//printf("the position - x,y,z %lf %lf %lf of %d\n ",d->r[i-1][0],d->r[i-1][1],d->r[i-1][2],i);
		// we need the acrobatics above to make sure that the data pertaining to "i" atom is stored in the "i" line.
		// since we need to pull out the id before we read the rest, we need to split the line read into two statements
		// important for calling MSD etc.
		// make sure for extra security that "dump modify sort id" is used in the LAMMPS script.
		// CHECK THIS LINE!
	}
	//printf("seg fault here?? %d\n",s->Ato);
	//printf("last of the atomcans - %f %f %f\n",d->r[s->Ato-1][0],d->r[s->Ato-1][1],d->r[s->Ato-1][2]);
	//printf("first of the atomcans - %f %f %f\n",d->r[0][0],d->r[0][1],d->r[0][2]);
	printf("\nread atom positions?? - %d\n",t);
}
