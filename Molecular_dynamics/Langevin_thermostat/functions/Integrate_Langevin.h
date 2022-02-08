void integrateLangevin(struct SYSTEM *box,struct PARAINPUT *p,struct APPOGGIO *a)
{
    // this subroutine performs a Velocity Verlet integration time step for each of our "particles"
    
    int i;
    double displacesq;
    double rgx,rgy,rgz,fRx,fRy,fRz;
    
    for(i=0;i<p->nAto;i++)
    {

        
        box->vx[i] = box->vx[i]*p->fricCoeff + 0.5*(box->oldfx[i]+box->fx[i])*p->dt;
        box->vy[i] = box->vy[i]*p->fricCoeff + 0.5*(box->oldfy[i]+box->fy[i])*p->dt;
        box->vz[i] = box->vz[i]*p->fricCoeff + 0.5*(box->oldfz[i]+box->fz[i])*p->dt;
        
//      mol[i].px = mol[i].px + mol[i].vx * MDdeltaTd + 0.5 * mol[i].ax * MDdeltaTd *MDdeltaTd;
        
        box->X[i] = box->X[i] + box->vx[i]*p->dt + 0.5*box->fx[i]*p->dt*p->dt;
        box->Y[i] = box->Y[i] + box->vy[i]*p->dt + 0.5*box->fy[i]*p->dt*p->dt;
        box->Z[i] = box->Z[i] + box->vz[i]*p->dt + 0.5*box->fz[i]*p->dt*p->dt;
        
        // include random component to forces
        
//         rgx = 2.0*curand_uniform (&localState)-1.0;
//         fRx = (sqrt(6.0* gamma* temperature / (MDdeltaT))) * rgx;
        
        box->oldfx[i]=box->fx[i];
        box->oldfy[i]=box->fy[i];
        box->oldfz[i]=box->fz[i];
        
        rgx = 2.0*genrand_real1(&(*a))-1.0;
        rgy = 2.0*genrand_real1(&(*a))-1.0;
        rgz = 2.0*genrand_real1(&(*a))-1.0;
        
        fRx = sqrt(6.0*p->gamma*box->tem/p->dt)*rgx;
        fRy = sqrt(6.0*p->gamma*box->tem/p->dt)*rgy;
        fRz = sqrt(6.0*p->gamma*box->tem/p->dt)*rgz;
        
        box->fx[i]=fRx;
        box->fy[i]=fRy;
        box->fz[i]=fRz;
        
    }
    
    
}

