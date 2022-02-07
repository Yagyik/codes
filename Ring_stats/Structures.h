struct SIMDAT  // Variables never changing during analysis
{
	int Step;			// Every stepAvg data are stored
	int Ato;     		// Number of Atoms
	double DeltaT;		// Integration time
	char name[300];		// Name of data files or path
	char buffer[300];	// Useful for file name
	char output[300];	// Output path
	int StartAt;		// First file for analysis
	int StopAt;			// Last file for analysis
	int Gr;				// Radial Distrib Funct
	int Gr_limit;
	int Msd;			// Mean Square Displ
	int Msd_limit;
	int Isf;			// Inter Scatt Funct Self
	int Isf_limit;	
	double IsfK_value;	// K of ISF-self
	int OrPar;			// Steinhardt parameters
	int OrPar_limit;	
	double OrPar_r;		// Cutoff for neighbors
	int OrPar_bin;		// Number of bin for histogram
	int Sq;				// Static structure factor
	int Sq_limit;
	double Sq_max;		// Max value of Q
	int Sq_bin;
	int Coher;			// Inter Scatt Funct Dist
	int Coher_limit;	
	double Coher_value;	// K of ISF-dist
	int G5r;					// Fifth neighbor distribution 
	int G5r_limit;
	double G5r_max;		// Max value of r
	int ClSize;				// Max Cluster size 
	int Cl_limit;
	double Cl_r;			// Cutoff for neighbors
	int Cl_solid;			// Number of solid particle
	double Cl_cut;			// Cutoff for solid particle
	int coordno_flag;
	int cn_limit;
	double cn_cutoff; 
	int Interface; 		// flag for whether we id interfaces.
	int Interface_limit; 
	double interface_tolerance;
	int Ring;
	int Ring_limit;
	double ring_cut;
	int entropy;
	int entropy_limit;
	int epnbin;
	double epra;
	double eprm;
	double epgm_lim;
	double epsigma;
	int OrPar_la;
	int laOrPar_limit;
	int laOrPar_bin;
	int laqnbin;
	double laOrPar_r;
	double laOrPar_ra;
	
	
	
				
} simdat; 

struct DATASIM  // Variables updating each step
{
	double **r;     // Position
	double **rv;    // Velocity
	double lBox[3];	// Lx - Ly - Lz
} data;

void readPara(FILE *fp, SIMDAT *s, DATASIM *d)		// You must modify when add a function
{
	char string[100];

	if(fp == NULL) 
	{
		printf("%s not found!\n",s->buffer);
		exit(0);
	}
	
	fscanf(fp, "%s",string); 			// Step
	fscanf(fp, "%s",string);
	s->Step=atoi(string);
	fscanf(fp, "%s",string);			// Ato
	fscanf(fp, "%s",string); 	
	s->Ato=atoi(string);
	printf("n atoms %d\n",s->Ato);
	d->r=AllocMatR(s->Ato,3);	
	d->rv=AllocMatR(s->Ato,3);	
	fscanf(fp, "%s",string); 			// DeltaT
	fscanf(fp, "%s",string);
	s->DeltaT=atof(string);
	fscanf(fp, "%s",string);					
	fscanf(fp, "%s", s->name); 			// FileName
	fscanf(fp, "%s",string);
	fscanf(fp, "%s",string); 			// StartAt
	s->StartAt=atoi(string);					
	fscanf(fp, "%s",string); 			// StopAt
	fscanf(fp, "%s",string);
	s->StopAt=atoi(string);
	fscanf(fp, "%s",string); 			// Gr
	fscanf(fp, "%s",string);
	s->Gr=atoi(string);					
	fscanf(fp, "%s",string); 			// Gr_limit
	fscanf(fp, "%s",string);
	s->Gr_limit=atoi(string);
	fscanf(fp, "%s",string); 			// Msd
	fscanf(fp, "%s",string);
	s->Msd=atoi(string);					
	fscanf(fp, "%s",string); 			// Msd_limit
	fscanf(fp, "%s",string);
	s->Msd_limit=atoi(string);
	fscanf(fp, "%s",string); 			// Isf
	fscanf(fp, "%s",string);
	s->Isf=atoi(string);					
	fscanf(fp, "%s",string); 			// Isf_limit
	fscanf(fp, "%s",string);
	s->Isf_limit=atoi(string);
	fscanf(fp, "%s",string); 			// IsfK_value
	fscanf(fp, "%s",string);										
	s->IsfK_value=atof(string);
	fscanf(fp, "%s",string); 			// OrPar
	fscanf(fp, "%s",string);
	s->OrPar=atoi(string);					
	fscanf(fp, "%s",string); 			// OrPar_limit
	fscanf(fp, "%s",string);
	s->OrPar_limit=atoi(string);
	fscanf(fp, "%s",string); 			// OrPar_r
	fscanf(fp, "%s",string);										
	s->OrPar_r=atof(string);
	fscanf(fp, "%s",string); 			// OrPar_bin
	fscanf(fp, "%s",string);										
	s->OrPar_bin=atoi(string);
	fscanf(fp, "%s",string); 			// Sq
	fscanf(fp, "%s",string);
	s->Sq=atoi(string);					
	fscanf(fp, "%s",string); 			// Sq_limit
	fscanf(fp, "%s",string);
	s->Sq_limit=atoi(string);
	fscanf(fp, "%s",string); 			// Sq_max
	fscanf(fp, "%s",string);
	s->Sq_max=atof(string);					
	fscanf(fp, "%s",string); 			// Sq_bin
	fscanf(fp, "%s",string);
	s->Sq_bin=atoi(string);
	fscanf(fp, "%s",string); 			// Coher
	fscanf(fp, "%s",string);
	s->Coher=atoi(string);					
	fscanf(fp, "%s",string); 			// Coher_limit
	fscanf(fp, "%s",string);
	s->Coher_limit=atoi(string);
	fscanf(fp, "%s",string); 			// Coher_value
	fscanf(fp, "%s",string);										
	s->Coher_value=atof(string);
	fscanf(fp, "%s",string); 			// G5r
	fscanf(fp, "%s",string);
	s->G5r=atoi(string);					
	fscanf(fp, "%s",string); 			// G5r_limit
	fscanf(fp, "%s",string);
	s->G5r_limit=atoi(string);
	fscanf(fp, "%s",string); 			// G5r_max
	fscanf(fp, "%s",string);
	s->G5r_max=atof(string);
	fscanf(fp, "%s",string); 			// ClSize
	fscanf(fp, "%s",string);
	s->ClSize=atoi(string);					
	fscanf(fp, "%s",string); 			// Cl_limit
	fscanf(fp, "%s",string);
	s->Cl_limit=atoi(string);
	fscanf(fp, "%s",string); 			// Cl_r
	fscanf(fp, "%s",string);
	s->Cl_r=atof(string);
	fscanf(fp, "%s",string); 			// Cl_solid
	fscanf(fp, "%s",string);
	s->Cl_solid=atoi(string);
	fscanf(fp, "%s",string); 			// Cl_cut
	fscanf(fp, "%s",string);
	s->Cl_cut=atof(string);
	fscanf(fp, "%s",string); 			// CoordNo yes/no
	fscanf(fp, "%s",string);
	s->coordno_flag=atoi(string);	
	printf("%d the flag \n",s->coordno_flag);	
	fscanf(fp, "%s",string); 			// cn_limit REMEMBER that this must be equal to OrPar lim and clsize lim
	fscanf(fp, "%s",string);
	//s->cn_limit=atoi(string);
	s->cn_limit = s->Cl_limit; 			// this needs to be true always
	fscanf(fp, "%s",string); 			// cn_cutoff
	fscanf(fp, "%s",string);
	s->cn_cutoff=atof(string);
	fscanf(fp, "%s",string); 			// Interface
	fscanf(fp, "%s",string);
	s->Interface = atoi(string);
	fscanf(fp, "%s",string); 			// Interface limit
	fscanf(fp, "%s",string);
	s->Interface_limit = atoi(string);
	fscanf(fp, "%s",string); 			// Interface tolerance
	fscanf(fp, "%s",string);
	s->interface_tolerance = atof(string);
	fscanf(fp, "%s",string); 			// Ring 3d
	fscanf(fp, "%s",string);
	s->Ring = atoi(string);
	fscanf(fp, "%s",string); 			// Ring 3d limit
	fscanf(fp, "%s",string);
	s->Ring_limit = atoi(string);
	fscanf(fp, "%s",string); 			// Ring 3d cutoff
	fscanf(fp, "%s",string);
	s->ring_cut = atof(string);
	fscanf(fp, "%s",string); 			// entropy
	fscanf(fp, "%s",string);
	s->entropy = atoi(string);
	fscanf(fp, "%s",string); 			// entropy limit
	fscanf(fp, "%s",string);
	s->entropy_limit = atoi(string);
	fscanf(fp, "%s",string); 			// epra
	fscanf(fp, "%s",string);
	s->epra = atof(string);
	fscanf(fp, "%s",string); 			// eprm
	fscanf(fp, "%s",string);
	s->eprm = atof(string);
	fscanf(fp, "%s",string); 			// first coord shell
	fscanf(fp, "%s",string);
	s->epgm_lim = atof(string);
	fscanf(fp, "%s",string); 			// epnbin
	fscanf(fp, "%s",string);
	s->epnbin = atoi(string);
	fscanf(fp, "%s",string); 			// epsigma
	fscanf(fp, "%s",string);
	s->epsigma = atof(string);
	fscanf(fp, "%s",string); 			// Order_la
	fscanf(fp, "%s",string);
	s->OrPar_la = atoi(string);
	fscanf(fp, "%s",string); 			// Order_la_limit
	fscanf(fp, "%s",string);
	s->laOrPar_limit = atoi(string);
	fscanf(fp, "%s",string); 			// Order_la_bin
	fscanf(fp, "%s",string);
	s->laOrPar_bin = atoi(string);
	fscanf(fp, "%s",string); 			// OrPar la distance bin
	fscanf(fp, "%s",string);
	s->laqnbin = atoi(string);
	fscanf(fp, "%s",string); 			// Order_la_r
	fscanf(fp, "%s",string);
	s->laOrPar_r = atof(string);
	fscanf(fp, "%s",string); 			// Order_la_ra
	fscanf(fp, "%s",string);
	s->laOrPar_ra = atof(string);
	
}




