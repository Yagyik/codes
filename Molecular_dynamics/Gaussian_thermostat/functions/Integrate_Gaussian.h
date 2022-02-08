void integrateGaussian(struct SYSTEM *box,struct PARAINPUT *p, int step)
{
    // this subroutine performs a Velocity Verlet integration time step for each of our "particles"
    
    int i;
    double displacesq;
    double alpha,denom,betaSq,h,f;
    
    
    
    
    if(step==0)
    {
        alpha=0.0;
        denom=0.0;
        betaSq=0.0;
        h=0.0;
        f=0.0;
        // CARE constant mass assumed
        for(i=0;i<p->nAto;i++)
        {
            alpha += box->fx[i]*box->vx[i] + box->fy[i]*box->vy[i] + box->fz[i]*box->vz[i];
            
            betaSq += box->fx[i]*box->fx[i] + box->fy[i]*box->fy[i] + box->fz[i]*box->fz[i];
            
            denom += box->vx[i]*box->vx[i] + box->vy[i]*box->vy[i] + box->vz[i]*box->vz[i];
                        
        }
        betaSq=betaSq/denom;
        alpha=alpha/denom;
        f=exp(-sqrt(betaSq)*0.5*p->dt);
        h = (alpha+sqrt(betaSq))/(alpha - sqrt(betaSq));
        
        for(i=0;i<p->nAto;i++)
        {
            box->vx[i] = (1-h)/(f-(h/f))*(box->vx[i] + box->fx[i]*(1+h-f-(h/f))/((1-h)*sqrt(betaSq)));
            box->vy[i] = (1-h)/(f-(h/f))*(box->vy[i] + box->fy[i]*(1+h-f-(h/f))/((1-h)*sqrt(betaSq)));
            box->vz[i] = (1-h)/(f-(h/f))*(box->vz[i] + box->fz[i]*(1+h-f-(h/f))/((1-h)*sqrt(betaSq)));
        }
        
        for(i=0;i<p->nAto;i++)
        {
            box->X[i] = box->X[i] + p->dt*box->vx[i];
            box->Y[i] = box->Y[i] + p->dt*box->vy[i];
            box->Z[i] = box->Z[i] + p->dt*box->vz[i];
            box->oldfx[i]=box->fx[i];
            box->oldfy[i]=box->fy[i];
            box->oldfz[i]=box->fz[i];
            
        }
    }
    
    if(step==1)
    {
        alpha=0.0;
        denom=0.0;
        betaSq=0.0;
        h=0.0;
        f=0.0;
        // CARE constant mass assumed
        for(i=0;i<p->nAto;i++)
        {
            alpha += box->fx[i]*box->vx[i] + box->fy[i]*box->vy[i] + box->fz[i]*box->vz[i];
            
            betaSq += box->fx[i]*box->fx[i] + box->fy[i]*box->fy[i] + box->fz[i]*box->fz[i];
            
            denom += box->vx[i]*box->vx[i] + box->vy[i]*box->vy[i] + box->vz[i]*box->vz[i];
                        
        }
        betaSq=betaSq/denom;
        alpha=alpha/denom;
        f=exp(-sqrt(betaSq)*0.5*p->dt);
        h = (alpha+sqrt(betaSq))/(alpha - sqrt(betaSq));
        
        for(i=0;i<p->nAto;i++)
        {
            box->vx[i] = (1-h)/(f-(h/f))*(box->vx[i] + box->fx[i]*(1+h-f-(h/f))/((1-h)*sqrt(betaSq)));
            box->vy[i] = (1-h)/(f-(h/f))*(box->vy[i] + box->fy[i]*(1+h-f-(h/f))/((1-h)*sqrt(betaSq)));
            box->vz[i] = (1-h)/(f-(h/f))*(box->vz[i] + box->fz[i]*(1+h-f-(h/f))/((1-h)*sqrt(betaSq)));
        }
        
    }
        
    
    
}

