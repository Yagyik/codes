#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"WidomFunctions/WidomStructures.h"
#include"WidomFunctions/readLammps_Widom.h"


int main (int argc, char *argv[])
{
    struct DATASIM d;
    
    FILE *fpPara,*fpConf,*fpOut;
    char buffer[500],string[500],buffer2[500];
    int i,j,t;
    fpPara=fopen(argv[1],"r");
    if(fpPara==NULL)
    {
        printf("didn't fine para file\n");
        
    }
   
    double spamx,spamy,spamz,rijsq,rij;
    double dx,dy,dz,en;
    
    fscanf(fpPara, "%s",string); 			// Step
	fscanf(fpPara, "%s",string);
	d.Step=atoi(string);
	fscanf(fpPara, "%s",string);			// Ato
	fscanf(fpPara, "%s",string); 	
	d.nAto=atoi(string);
	printf("n atoms %d\n",d.nAto);
	fscanf(fpPara, "%s",string);					
	fscanf(fpPara, "%s", d.name); 			// FileName
	fscanf(fpPara, "%s",string);
	fscanf(fpPara, "%s",string); 			// StartAt
	d.StartAt=atoi(string);					
	fscanf(fpPara, "%s",string); 			// StopAt
	fscanf(fpPara, "%s",string);
	d.StopAt=atoi(string);
	fscanf(fpPara, "%s",string); 			// nInsert
	fscanf(fpPara, "%s",string);
    d.nInsert=atoi(string);
    fscanf(fpPara, "%s",string); 			// eps
	fscanf(fpPara, "%s",string);
    d.eps=atof(string);
    fscanf(fpPara, "%s",string); 			// sigma
	fscanf(fpPara, "%s",string);
    d.sigma=atof(string);
    fscanf(fpPara, "%s",string); 			// rc
	fscanf(fpPara, "%s",string);
    d.rc=atof(string);
    fscanf(fpPara, "%s",string); 			// rc
	fscanf(fpPara, "%s",string);
    d.T=atof(string);
    
    d.x=calloc(d.nAto,sizeof(double));
    d.y=calloc(d.nAto,sizeof(double));
    d.z=calloc(d.nAto,sizeof(double));
    d.vx=calloc(d.nAto,sizeof(double));
    d.vy=calloc(d.nAto,sizeof(double));
    d.vz=calloc(d.nAto,sizeof(double));
    
    
    double TS_fact=4.0*d.eps*(pow(d.sigma/d.rc,12.0)-pow(d.sigma/d.rc,6.0));
    double avgDeltaU=0.0;
    sprintf(buffer2,"%s",argv[3]);
    fpOut=fopen(buffer2,"w");
    for(int t=d.StartAt;t<=d.StopAt;t+=d.Step)
    {
        printf("attempting config %d\n",t);
        sprintf(buffer,"%s/%s%d.dat", argv[2],d.name,t);
        fpConf=fopen(buffer,"r");	
		if(fpConf==NULL)
		{
			printf("null pointer problems\n");
		}		
		readLammps(fpConf,&d,t);
        fclose(fpConf);
        avgDeltaU=0.0;
        for(i=0;i<d.nInsert;i++)
        {
            spamx=d.lBox[0]*(drand48()-0.5);
            spamy=d.lBox[1]*(drand48()-0.5);
            spamz=d.lBox[2]*(drand48()-0.5);
            
            en=0.0;
            
            for(j=0;j<d.nAto;j++)
            {
                //printf("%d %f\n",j,d.x[j]);
                dx=d.x[j]-spamx;
                dx-=d.lBox[0]*lround(dx/d.lBox[0]);
                dy=d.y[j]-spamy;
                dy-=d.lBox[1]*lround(dy/d.lBox[1]);
                dz=d.z[j]-spamz;
                dz-=d.lBox[2]*lround(dz/d.lBox[2]);
                
                rijsq=dx*dx+dy*dy+dz*dz;
                if(rijsq<d.rc)
                {
                    rij=sqrt(rijsq);
                    en+=4.0*d.eps*(pow(d.sigma/rij,12.0)-pow(d.sigma/rij,6.0)) - TS_fact;
                    
                }
  
            }
            avgDeltaU+=exp(-en/d.T);
            //printf("did one insert round %d\n",i);
        }
        //printf("done inserting printing now\n");
        fprintf(fpOut,"%d %f\n",t,avgDeltaU/d.nInsert);
        
    }
    printf("done all configs\n");
    fclose(fpOut);
    

}
