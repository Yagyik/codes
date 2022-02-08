int *AllocVecI (int size) 
{
	int *vec;

	vec = (int *) malloc (size * sizeof (int));
	// allocates memory to a pointer thing
	if (vec == NULL)
	{
		printf("Allocation failed in AllocVecI\n");
	}

	return (vec);
} 


double *AllocVecR (int size) 
{
	double *vec;

	vec=(double *) malloc(size * sizeof (double));
	if(vec==NULL) 
	{
		printf("Allocation failed in AllocVecR\n");
	}

	return (vec);
}

int **AllocMatI (int size1, int size2)
{
	int **v;
	unsigned k;
	
	v = (int **) malloc (size1 * sizeof (int *));

	if (v == NULL)
	{
		printf ("Allocation failed in AllocMatR\n");
	}

	v[0] = (int *) malloc (size1 * size2 * sizeof (int));
	if (v[0] == NULL)
	{
		printf ("Allocation failed in AllocMatR\n");
	}

	for (k=0;k < (size1-1); k++)
	{
		v[k+1] = v[k] + size2;
	}

	return (v);
} 

double **AllocMatR (int size1, int size2)
{
	double **v;
	unsigned k;
	
	v = (double **) malloc (size1 * sizeof (double *));

	if (v == NULL)
	{
		printf ("Allocation failed in AllocMatR\n");
	}

	v[0] = (double *) malloc (size1 * size2 * sizeof (double));
	if (v[0] == NULL)
	{
		printf ("Allocation failed in AllocMatR\n");
	}

	for (k=0;k < (size1-1); k++)
	{
		v[k+1] = v[k] + size2;
	}

	return (v);
}

double ***AllocTensR (int size1, int size2, int size3) 
{
	double ***v;
	unsigned k;

	v = (double ***) malloc (size1 * sizeof (double **));
	
	if (v == NULL)
	{
		printf ("Allocation failed in AllocTensR\n");
	}

	v[0] = (double **) malloc (size1 * size2 * sizeof (double *));
	if (v[0] == NULL)
	{
		printf ("Allocation failed in AllocTensR\n");
	}

	v[0][0]=(double *) malloc (size1 * size2 * size3*sizeof (double));
	if (v[0] == NULL)
	{
		printf ("Allocation failed in AllocTensR\n");
	}

	for (k=0; k< (size2*size1-1) ;k++)
	{
		v[0][k+1]=v[0][k] +size3;
	}
	for (k=0;k < (size1-1);k++)
	{
		v[k+1]=v[k] + size2 ; 
	}

	return (v);
}

double ****AllocTensRang4R (int size1, int size2, int size3, int size4) 
{
	double ****v;
	unsigned k;

	v = (double ****) malloc (size1 * sizeof (double ***));
	if (v == NULL)
	{
		printf ("Allocation failed in AllocTens4R\n");
	}

	v[0] = (double ***) malloc (size1 * size2 * sizeof (double **));
	if (v[0] == NULL)
	{
		printf ("Allocation failed in AllocTens4R\n");
	}

	v[0][0]=(double **) malloc (size1 * size2 * size3*sizeof (double *));
	if (v[0][0] == NULL)
	{
		printf ("Allocation failed in AllocTens4R\n");
	}

	v[0][0][0] = (double *) malloc (size1 * size2 * size3 * size4*sizeof (double));
	if (v[0][0][0] == NULL) 
	{
		printf ("Allocation failed in AllocTens4R\n");
	}

	for (k=0; k<(size3*size2*size1-1); k++)
	{
		v[0][0][k+1]=v[0][0][k]+size4;
	}
	for (k=0; k<(size2*size1-1);k++)
	{
		v[0][k+1]=v[0][k]+size3;
	}
	for (k=0; k<(size1-1);k++)
	{
		v[k+1]=v[k]+size2;
	}

	return (v);
}

double *****AllocTensRang5R (int size1, int size2, int size3, int size4, int size5) 
{
	double *****v;
	unsigned k;

	v = (double *****) malloc (size1 * sizeof (double ****));
	if (v == NULL)
	{
		printf ("Allocation failed in AllocTens5R\n");
	}

	v[0] = (double ****) malloc (size1 * size2 * sizeof (double ***));
	if (v[0] == NULL)
	{
		printf ("Allocation failed in AllocTens5R\n");
	}

	v[0][0]=(double ***) malloc (size1 * size2 * size3*sizeof (double **));
	if (v[0][0] == NULL)
	{
		printf ("Allocation failed in AllocTens5R\n");
	}

	v[0][0][0] = (double **) malloc (size1 * size2 * size3 * size4*sizeof (double *));
	if (v[0][0][0] == NULL) 
	{
		printf ("Allocation failed in AllocTens5R\n");
	}

	v[0][0][0][0] = (double *) malloc (size1 * size2 * size3 * size4 * size5*sizeof (double));
	if (v[0][0][0][0] == NULL)
	{
		printf ("Allocation failed in AllocTens5R\n");
	}

	for (k=0; k<(size4*size3*size2*size1-1); k++)
	{
		v[0][0][0][k+1]=v[0][0][0][k] + size5;
	}
	for (k=0; k<(size3*size2*size1-1);k++)
	{
		v[0][0][k+1]=v[0][0][k]+size4;
	}
	for (k=0; k<(size2*size1-1);k++)
	{
		v[0][k+1]=v[0][k]+size3;
	}
	for (k=0; k<(size1-1);k++)
	{
		v[k+1]=v[k]+size2;
	}

	return (v);
}

void FreeVecR (double *v)
{
  free (v);
}

void FreeMatR (double **v)
{
  free (v[0]);
  free (v);
}

void FreeTensR (double ***v)
{
  free (v[0][0]);
  free (v[0]);
  free (v);
}

void FreeTens4R (double ****v)
{
  free (v[0][0][0]);
  free (v[0][0]);
  free (v[0]);
  free (v);
}

void FreeTens5R (double *****v)
{
  free (v[0][0][0][0]);
  free (v[0][0][0]);
  free (v[0][0]);
  free (v[0]);
  free (v);
}


