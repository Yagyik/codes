!!!spamming keeps moba allliiiiiivvvvveee!!!
!!!allliiiivvvvvvvveeeee!!
!!!! psam to stay alive!!!!
!!!! Oh-whoa I heard, I'm still aliveeee! 



Module define_a2_q2
!!parameters
integer,parameter::N = 1000
real, parameter:: T=1.32
real,parameter:: beta = 1/1.32
real,parameter::kb=1
integer,parameter::maxsteps=50000
real,parameter::eps=1
real,parameter::sigma=1
real,parameter::rv = 3.0*sigma
real,parameter::rc = 2.5*sigma
real,parameter::L = 6.21
real,parameter::dr = 0.02*L
integer,parameter::nexc = 400 !!! ??? what do we put here? the less it is the more we move stuff 
real,parameter::h = 0.186 !!! some funky argon units from Shiladitya's thesis
!! but nexc should not change what our final density is since whether we add or not depends on arggggghhhh!


!!loop indices and accumulators and flags and SPAMMMMMMM!!!
integer::i,i1,ii,j1,j,k,k1,rhoi,dri,step,flag,o,fileind,spam,mui,npav,npart,moveorexc
real :: dist,dx,enn,eno,en,r,rsq,V,arg,rho,xen,mu,mu_id,mu_ex,mu_tot,dmu,rhoe,zz,vol,lambda,coin,rho_avg

!!arrays
real,dimension(N,3)::x
real,dimension(3)::xn,xr,xi,xj,xtest


!!!dummy dum dum
 character(99)::filename
End module define_a2_q2

Program gcmc
Use define_a2_q2
Implicit none
!open(unit=55,file='density_relax.dat',action='write')
open(unit=66,file='density_mu_50k.dat',action='write')
dmu = 0.15
Do mui = 1,40 !!! need to know how many data points
Mu = -1.8 + mui*dmu !! -> which mu is this, total or excess??
! calc npav for each mu from expected density according to the widom thing
Rhoe = 0.40 + 0.015*mui !! rho predicted from Widom curve
! calc npav for each mu from expected density according to the widom thing
Npav = nint(rhoe*(L**3))
npart = 108 !! check the implications of this initialisation
vol = l**3
!print*,npav,"what's helping us decide??"

!!!!section calculating zz and lambda
!lambda = h/sqrt(2*3.14*T)
lambda = 1.0 !! we do this because h seems to be dominating everything
zz = exp(beta*mu)/(lambda**3)


!!!! the value of h is a bit of an issue
!!!! further, whether the mu there is mu_tot or mu_ex is another issue
!!!! if it were mu_total then it would depend on rho - clearly rhoe -> the expected value
!!!! one possible conclusion is that the mu we are looping over is simply the excess part as is the mu in zz.
rho_avg = 0.0
Call init()
print*,mui,rho_avg,mu,npav,"next sample!!"
  Do step = 1,maxsteps !! -> max steps not mc steps
  Moveorexc = int(rand()*(npav+nexc)) + 1 !! -> npav is what?? Update :fixed

    If (moveorexc < npart)then ! move stuff
      !print*,step,"moving??"
      !print*,x(1,:),x(npart,:),"before"
      Call mcmove(dr)
      !print*,x(1,:),x(npart,:),"after"
    Else
      Call mcexc
    Endif
!!! ⇒ ⇒ do we sample here or what?
!!! ⇒ first keep track of density.
  Rho = real(npart)/(L**3)
!!! ⇒ keep track of ideal part of mu
!!! maybe print it to a file also
  Mu_id = -T*3*log(lambda) + T*log(rho) !!! lambda_factor is -kb*T*ln(lambda**3)
  mu_ex = mu - mu_id ! first possibility
  mu_tot = mu + mu_id ! if we walk a different path where the mu we vary is mu_ex
  !if(mod(step,50)==0)then
    !print*,npart,rho,L
 ! endif
!! and maybe, just maybe we only care about mu and some post-equilibration average value of rho or the final value of rho
  !write(55,*)step, rho
!!! here we will accumulate after some 15000 MC steps till 20000
  if(step>20000)then
    rho_avg = rho_avg + rho ! div by 5000 for weighting
  endif

  End do !! step
!! after a long number of steps, if rho has settled down, 
!!We plot it against mu - mu_id or simply mu (in which case mu is just mu_ex) 
!! store pairs in file
!print*,npav,"npav"
!print*,lambda,zz,"first up"
write(66,*)rho_avg/(maxsteps-20000),mu

Enddo !! mui
 !close(55)
 !close(66)

Contains

subroutine init
implicit none
x(:,:) = 0.0 !!-> X(N,3) array
!print*,x,"x zeros"
x(1,:) = [rand()*L,rand()*L,rand()*L]

Do i=2,npav !! we put npav here because we are worried about time of init
  spam=0
  do while(.true.)
  spam = spam+1
  x(i,:) = [rand()*L,rand()*L,rand()*L]
  flag = 0
  Do j=1,i-1
  dist = 0.0
    do k=1,3
    dx = x(i,k) - x(j,k)
    dx = dx - L*anint(dx/L) !!!!!!!! put in PBC
    dist = dist + (dx)**2
    enddo 
  if(spam > npav)then
  print*,"spam exit    ----->>> this is an ERROR"
  exit
  end if
  if(dist < 0.64*sigma**2) then
  flag = 1
  exit
  endif
  End do !!– j
  If(flag/=1) then !!!! check this syntax
  exit ! means x(i,:) has been fixed
  endif
  end do ! at the end of each while, we loop over one more "i"
end do !!“I”
end subroutine init


real function phi(i,j,xi,xj)
implicit none 
integer,intent(in)::i,j
real,dimension(3),intent(in)::xi,xj
integer::d1
real::dx1,rsq1,r1
dx1  =0.0 ; r1   =0.0 ; rsq1 =0.0

do d1=1,3
dx1 = xi(d1) - xj(d1)
dx1 = dx1 - L*anint(dx1/L)
rsq1 = rsq1 + dx1**2
enddo

r1 = sqrt(rsq1)
if(r1<=rc)then
phi = 4*eps*((sigma/r1)**12 - (sigma/r1)**6) - 4*eps*((sigma/rc)**12 - (sigma/rc)**6)
elseif(r1>rc)then
phi = 0.0
endif
end function phi

!!!calc energy!!!

subroutine energy(i1,xi,en)
implicit none
real,intent(inout)::en
integer,intent(in)::i1
real,dimension(3),intent(in)::xi
en = 0.0
do j1 = 1,npart
xen=0.0
if(j1 /= i1) then
en  = en + phi(i1,j1,xi,x(j1,:))
end if
enddo
end subroutine

!!!!Monte Carlo movin'

Subroutine mcmove(dr)
implicit none
real,intent(in)::dr
do ii = 1,npart
o = int(rand()*npart)+1
xn(:)=0
do k1=1,3
  xn(k1) = x(o,k1) + (rand()-0.5)*dr !!- -> put in PBC after
  xn(k1) = xn(k1) - nint(xn(k1)/L)*L
enddo
  Call energy(o,x(o,:),eno)
  Call energy(o,xn(:),enn)
arg = exp(-beta*(enn-eno))
!print*,arg,"always always always print arrrrrrrgggggghhhh!!"
If(rand() < arg) then
x(o,:) = xn(:)
Endif
enddo !"ii"
end subroutine mcmove

Subroutine mcexc
Implicit none
coin = rand()
If(coin<0.5) then !!!remove a particle
  !print*,step,"gonna remove!!"
  if (npart /=0)then
    o = int(npart*rand()) + 1
    call energy(o,x(o,:),eno)
    Arg = real(npart)*exp(beta*eno)/(zz*vol)
    !print*,eno,zz,vol,"eno, zz and vol"
    !print*,arg,"alwayss print arggghhh - removing"
    if (rand()<arg) then
      x(o,:) = x(npart,:) !!! exchange oth with last
      Npart = npart - 1 !!! stop keeping track of the last one
      !print*,"removed"
    Endif
  End if
Else
  !print*,step,coin,"gonna add!!"
  xn = [L*rand(),L*rand(),L*rand()] ! new particle position
  Call energy(npart+1,xn,enn)
  Arg = exp(-beta*enn)*zz*vol/real(npart+1)
  !print*,arg,"alllllwaaaayyyyysssss prinntttt arggggghhhhh - adding!"
  if(rand()<arg) then !!!add the particle
      x(npart+1,:) = xn
      Npart = npart+1
      !print*,"removed"
  Endif
endif
End subroutine mcexc

end program gcmc

