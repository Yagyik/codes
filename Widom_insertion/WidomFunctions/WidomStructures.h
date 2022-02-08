struct DATASIM  // Variables updating each step
{
     int StartAt,StopAt,Step,nAto,nInsert;
	double *x,*y,*z;     // Position
	double *vx,*vy,*vz;
	
    double eps,sigma,rc,T;
    
	double lBox[3];	// Lx - Ly - Lz
	char name[500];
} data;
