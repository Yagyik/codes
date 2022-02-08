double gauss_var_AT(double mu,double sig,struct APPOGGIO *a)
{
    double gauss_r,v1,v2,gauss_var;
    gauss_r=2.0;
    while(gauss_r>=1.0)
    {
        v1=2.0*genrand_real1(&(*a))-1.0;
        v2=2.0*genrand_real1(&(*a))-1.0;
//         printf("v1 v2 %f %f\n",v1,v2);
        gauss_r = v1*v1 + v2*v2;
    }
//     printf("gr log gr %f %f\n",gauss_r,log(gauss_r));
    gauss_var=v1*sqrt(-2.0*log(gauss_r)/gauss_r);
    gauss_var=mu+sig*gauss_var;
    
    return gauss_var;
}

void integrateVelVerlet(struct SYSTEM *box,struct PARAINPUT *p, int step)
{
    // this subroutine performs a Velocity Verlet integration time step for each of our "particles"
    
    int i;
    double displacesq;
    double rand=0.0;
    double sig=sqrt(structbox->tem);
    if(step==0)
    {
       for(i=0;i<p->nAto;i++)
        {
            box->X[i] = box->X[i] + p->dt*box->vx[i] + p->dt*p->dt*0.5*p->m*box->fx[i];
            box->Y[i] = box->Y[i] + p->dt*box->vy[i] + p->dt*p->dt*0.5*p->m*box->fy[i];
            box->Z[i] = box->Z[i] + p->dt*box->vz[i] + p->dt*p->dt*0.5*p->m*box->fz[i];
            box->oldfx[i]=box->fx[i];
            box->oldfy[i]=box->fy[i];
            box->oldfz[i]=box->fz[i];
            
            // check neighbourlist update condition
//             box->dispX[i]+=p->dt*box->vx[i] + p->dt*p->dt*0.5*p->m*box->fx[i];
//             box->dispY[i]+=p->dt*box->vy[i] + p->dt*p->dt*0.5*p->m*box->fy[i];
//             box->dispZ[i]+=p->dt*box->vz[i] + p->dt*p->dt*0.5*p->m*box->fz[i];
//             displacesq=sqr(box->dispX[i])+sqr(box->dispY[i])+sqr(box->dispZ[i]);
//             box->maxdispsq=MAX(box->maxdispsq,displacesq);
//             //printf("%d %f %f\n",par,box->maxdispsq,displacesq);
//             if(box->maxdispsq>sqr(0.5*SKIN))
//             {
//                 box->maxdispsq=0.0;
//                 vnlist(&(*p),&(*box));
//             }
            
        }
        
        
    }
    if(step==1)
    {
        for(i=0;i<p->nAto;i++)
        {
            box->vx[i] = box->vx[i] + p->dt*0.5*p->m*(box->oldfx[i] + box->fx[i]);
            box->vy[i] = box->vy[i] + p->dt*0.5*p->m*(box->oldfy[i] + box->fy[i]);
            box->vz[i] = box->vz[i] + p->dt*0.5*p->m*(box->oldfz[i] + box->fz[i]);
            
            
        }
        // sample velocities from heat bath with probability
        for(i=0;i<p->nAto;i++)
        {
            rand=genrand_real1(&(*a));
            if(rand<p->nu*p->dt)
            {
                box->vx[i]=gauss_var_AT(0,sig);
            }
        }
        
    }
    
    
}

