void savePre(struct SYSTEM *structbox,struct PARAINPUT *p)
{
	int i,j;

	structbox->nSp_pre=structbox->nSp;
	structbox->nSp_maxPre=structbox->nSp_max;
	structbox->nl4_pre=structbox->nl4;
	structbox->nl4_maxPre=structbox->nl4_max;
	structbox->nl5_pre=structbox->nl5;
	structbox->nl5_maxPre=structbox->nl5_max;
	
	structbox->Wpre=structbox->W;
	structbox->potEne_pre=structbox->potEne;
	structbox->lxpre=structbox->lx;
	structbox->lypre=structbox->ly;
	structbox->lzpre=structbox->lz;
	structbox->rho_pre=structbox->rho;

	for(i=1;i<=p->nAto+1;i++)
	{
		structbox->Xpre[i]=structbox->X[i];
		structbox->Ypre[i]=structbox->Y[i];
		structbox->Zpre[i]=structbox->Z[i];
        if(i<=p->nAto)
        {
            for(j=0;j<structbox->nlen;j++)
            {
                structbox->NEIGLISTpre[i][j]=structbox->NEIGLIST[i][j];
            }
        }
		

		structbox->hipre[i]=structbox->hi[i];
		structbox->siXpre[i]=structbox->siX[i];
		structbox->siYpre[i]=structbox->siY[i];
		structbox->siZpre[i]=structbox->siZ[i];
		structbox->TiXXpre[i]=structbox->TiXX[i];
		structbox->TiXYpre[i]=structbox->TiXY[i];
		structbox->TiXZpre[i]=structbox->TiXZ[i];
		structbox->TiYYpre[i]=structbox->TiYY[i];
		structbox->TiYZpre[i]=structbox->TiYZ[i];
		structbox->TiZZpre[i]=structbox->TiZZ[i];
		structbox->dispXpre[i]=structbox->dispX[i];
		structbox->dispYpre[i]=structbox->dispY[i];
		structbox->dispZpre[i]=structbox->dispZ[i];
        
//         structbox->Ahipre[i]=structbox->Ahi[i];
// 		structbox->AsiXpre[i]=structbox->AsiX[i];
// 		structbox->AsiYpre[i]=structbox->AsiY[i];
// 		structbox->AsiZpre[i]=structbox->AsiZ[i];
// 		structbox->ATiXXpre[i]=structbox->ATiXX[i];
// 		structbox->ATiXYpre[i]=structbox->ATiXY[i];
// 		structbox->ATiXZpre[i]=structbox->ATiXZ[i];
// 		structbox->ATiYYpre[i]=structbox->ATiYY[i];
// 		structbox->ATiYZpre[i]=structbox->ATiYZ[i];
// 		structbox->ATiZZpre[i]=structbox->ATiZZ[i];
// 		structbox->AdispXpre[i]=structbox->AdispX[i];
// 		structbox->AdispYpre[i]=structbox->AdispY[i];
// 		structbox->AdispZpre[i]=structbox->AdispZ[i];
	}
	
	for(i=0;i<NspUppLim;i++)
	structbox->currentHistopre[i]=structbox->currentHisto[i];
}

void replacePre(struct SYSTEM *structbox,struct PARAINPUT *p)
{
	int i,j;

	
	structbox->nSp=structbox->nSp_pre;
	structbox->nSp_max=structbox->nSp_maxPre;
	structbox->nl4=structbox->nl4_pre;
	structbox->nl4_max=structbox->nl4_maxPre;
	structbox->nl5=structbox->nl5_pre;
	structbox->nl5_max=structbox->nl5_maxPre;	
	
	structbox->W=structbox->Wpre;
	structbox->potEne=structbox->potEne_pre;
	structbox->lx=structbox->lxpre;
	structbox->ly=structbox->lypre;
	structbox->lz=structbox->lzpre;
	structbox->rho=structbox->rho_pre;

	for(i=1;i<=p->nAto+1;i++)
	{
		structbox->X[i]=structbox->Xpre[i];
		structbox->Y[i]=structbox->Ypre[i];
		structbox->Z[i]=structbox->Zpre[i];
        if(i<=p->nAto)
        {
            for(j=0;j<structbox->nlen;j++)
            {
                structbox->NEIGLIST[i][j]=structbox->NEIGLISTpre[i][j];
            }
        }
		

		structbox->hi[i]=structbox->hipre[i];
		structbox->siX[i]=structbox->siXpre[i];
		structbox->siY[i]=structbox->siYpre[i];
		structbox->siZ[i]=structbox->siZpre[i];
		structbox->TiXX[i]=structbox->TiXXpre[i];
		structbox->TiXY[i]=structbox->TiXYpre[i];
		structbox->TiXZ[i]=structbox->TiXZpre[i];
		structbox->TiYY[i]=structbox->TiYYpre[i];
		structbox->TiYZ[i]=structbox->TiYZpre[i];
		structbox->TiZZ[i]=structbox->TiZZpre[i];
		structbox->dispX[i]=structbox->dispXpre[i];
		structbox->dispY[i]=structbox->dispYpre[i];
		structbox->dispZ[i]=structbox->dispZpre[i];
        
//         structbox->Ahi[i]=structbox->Ahipre[i];
// 		structbox->AsiX[i]=structbox->AsiXpre[i];
// 		structbox->AsiY[i]=structbox->AsiYpre[i];
// 		structbox->AsiZ[i]=structbox->AsiZpre[i];
// 		structbox->ATiXX[i]=structbox->ATiXXpre[i];
// 		structbox->ATiXY[i]=structbox->ATiXYpre[i];
// 		structbox->ATiXZ[i]=structbox->ATiXZpre[i];
// 		structbox->ATiYY[i]=structbox->ATiYYpre[i];
// 		structbox->ATiYZ[i]=structbox->ATiYZpre[i];
// 		structbox->ATiZZ[i]=structbox->ATiZZpre[i];
// 		structbox->AdispX[i]=structbox->AdispXpre[i];
// 		structbox->AdispY[i]=structbox->AdispYpre[i];
// 		structbox->AdispZ[i]=structbox->AdispZpre[i];
	}
	
	for(i=0;i<NspUppLim;i++)
	structbox->currentHisto[i]=structbox->currentHistopre[i];
}

void saveConf(struct SYSTEM *structbox,struct SYSTEM *tmpbox,struct PARAINPUT *p)
{
	int i,j;
	
	tmpbox->lx=structbox->lx;
	tmpbox->ly=structbox->ly;
	tmpbox->lz=structbox->lz;

	for(i=1;i<=p->nAto+1;i++)
	{
		tmpbox->X[i]=structbox->X[i];
		tmpbox->Y[i]=structbox->Y[i];
		tmpbox->Z[i]=structbox->Z[i];
        if(i<=p->nAto)
        {
            for(j=0;j<structbox->nlen;j++)
            {
                tmpbox->NEIGLIST[i][j]=structbox->NEIGLIST[i][j];
            } 
        }
		
	}
}


