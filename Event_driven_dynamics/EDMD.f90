Module define_q4
implicit none

!!parameters
integer,parameter::N=10
integer,parameter::NenBin=50
real(kind=8),parameter::a=2
integer,parameter::coll_max = 20000
integer,parameter::coll_block=1000
real(kind=8),parameter::vyy=1
real(kind=8),parameter::sigma = 0.2


!!!sim vars
real(kind=8),dimension(N*N,2)::x,v,vold
real(kind=8),dimension(N*N)::ke
real(kind=8),dimension(NenBin)::pke
real(kind=8),dimension(coll_max/1000,coll_max*N*N)::vx,vy
real(kind=8),dimension(N*N)::tmin
integer,dimension(N*N)::indices
!!!!loop indices, counters and dummies
integer::i,j,j11,k,p,i1,block_count,coll_count,spamcount,first,second,fileind,fileind2,spam
real(kind=8)::vijsq,rijsq,bij,spam_r,t_tot,tij,rij,tgmin,spam_x,energy,enBin,enMax,hfunc,denom
 character(99)::filename,filename2


End module define_q4



Program eventful
Use define_q4
Implicit none

vold(:,:)=0
coll_count=0
Block_count = 1
call init()
fileind = 10+block_count
fileind2 = 30+block_count
write(filename,'("velocity_profile_x_",I1,".dat")') block_count
write(filename2,'("velocity_profile_y_",I1,".dat")') block_count
print*,fileind,fileind2,"indices first"
open(unit=fileind,file=filename,action='write')
open(unit=fileind2,file=filename2,action='write')
open(unit=99,file='coll_times.dat',action='write')
open(unit=88,file='energy.dat',action='write')
open(unit=79,file='hFunc.dat',action='write')
Do while(coll_count<=coll_max)
  !print*,coll_count,"colliding?"
  Do i=1,N*N
    call find_min(i)
  Enddo
  tgmin = 100000
  first = 0
  second = 0
  Do i=1,N*N
    If(tmin(i) < tgmin)then
        tgmin = tmin(i)
        first = i
        second = indices(i)
    Endif
  Enddo
  write(99,*)coll_count,tgmin
  !print*,tgmin,first,second,"colliders"
!! found global minimum
  !print*,x(:,1),v(:,1),"before evolution"
  Call evolve(tgmin)
  !print*,x(:,1),v(:,1),"after evolution"
  !print*,first,second,tgmin,"colliding!!!!"
  Call collide(first, second)
  !print*,x(:,1),v(:,1),"after collision"
  energy = 0.0
  do i=1,N*N
!     if((abs(v(i,1)-vold(i,1)) > 1).or.(abs(v(i,2)-vold(i,2))>100))then
!     !print*,i,vold(i,1),vold(i,2),"oldddd!"
!     !print*,i,v(i,1),v(i,2),"tag!! you're it!!"
!     endif
    vold(i,1) = v(i,1)
    vold(i,2) = v(i,2)
    energy = energy + 0.5*(v(i,1)**2 + v(i,2)**2)
  enddo
  write(88,*)coll_count,energy
  coll_count = coll_count+1  
  t_tot = t_tot + tgmin
!!!!include stats
  If((coll_count< block_count*coll_block).and.(coll_count>(block_count-1)*coll_block))then
    Do i=1,N*N
        vx(block_count,spamcount) = v(i,1)
        vy(block_count,spamcount) = v(i,2)
        !! debug checks here... spamcount index going to N*N*1000??
        !print*,fileind,fileind2,"indices files"
        write(fileind,*) spamcount,vx(block_count,spamcount)
        write(fileind2,*) spamcount,vy(block_count,spamcount)
        spamcount = spamcount+1
        !print*,spamcount,"on purpose"
    Enddo
  !print*,spamcount,"how many prints"
  endif
  !! check the following condition for switching blocks
  If(mod(coll_count,coll_block)==0)then
    !!here we change the block_count => means next block and next file => means close old file and open new one... 
    !! debug checks - for 10 blocks are we opening 10 files? do they have the right identifier? 
    close(fileind)
    close(fileind2)
    Block_count = block_count + 1
    print*,coll_count,block_count,fileind,"closing??"
    !! every time we change blocks, we reset spamcount to 0
    spamcount=1
    ! debug check, see if we are resetting spamcount and updating block_count at correct times wrt coll_count
    write(filename,'("velocity_profile_x_",I1,".dat")') block_count
    write(filename2,'("velocity_profile_y_",I1,".dat")') block_count
    fileind = 10+block_count
    fileind2 = 30+block_count
    print*,fileind,fileind2,"indices files in new file"
    open(unit=fileind,file=filename,action='write')
    open(unit=fileind2,file=filename2,action='write')
    call calc_hfunc()
    write(79,*)t_tot,hfunc
    
  endif
Enddo !! for the while loop
 close(99)




contains

Subroutine find_min(i1)
implicit none
integer,intent(in)::i1
!!!! remember to check nearest image dist
tmin(i1)=100000
indices(i1)=0
!print*,i1,"called, anon!!!"
Do j11=1,N*N
bij = 0.0
spam_r=0.0
If (i1/=j11)then
  !print*,i1,j,coll_count,"who and when?"
  spam_r =0.0
  Do k=1,2
   spam_r = x(i1,k)-x(j11,k)
   !print*,spam_r,"distances"
   spam_r = spam_r - anint(spam_r/(N*a))*N*a !!! --->>> Put in PBC
   !print*,spam_r,"fixed distances"
   bij = bij + (spam_r)*(v(i1,k)-v(j11,k))
   !print*,bij,"realllyyyy??????"
  Enddo
  
  vijsq= 0.0
  rijsq=0.0
  spam_r = 0.0
  Do k=1,2
    Vijsq = vijsq + (v(i1,k)-v(j11,k))**2
    spam_r = x(i1,k) - x(j11,k)
    spam_r = spam_r - anint(spam_r/(N*a))*N*a   !!! ->>> minimum image distance
    rijsq = rijsq + (spam_r)**2
  Enddo 
  

!print*,i1,j,vijsq,rijsq,bij,"vijsq,rijsq,bij"
if(bij > 0)then
    !print*,i1,j11,"cycling,not colliding"
    Cycle  !!----->>> loop continue here to skip calc the rest for current j and go to next j
  elseif(bij**2 - vijsq*(rijsq - (sigma)**2)<0)then
    !print*,i1,j11,"cycling"
    Cycle
  else 
  tij = (-bij - sqrt(bij**2 - vijsq*(rijsq - (sigma)**2)))/vijsq
  if(tij<0)then
    print*,i1,j11,"the suspects"
    print*,bij,bij**2 - vijsq*(rijsq - (sigma)**2),"the weapons"
    print*,vijsq,rijsq,"parts of a machete"
    print*,tij,"the crime"
  endif
  endif
 
!print*,i1,j,bij**2 - vijsq*(rijsq - sigma**2),"are we supposed to be here??"

  If(tij <tmin(i1)) then
  tmin(i1) = tij
  indices(i1) = j11
  Endif

Endif !!- (j/=i

End do !!! j11
end subroutine find_min

Subroutine evolve(tgmin)
implicit none
real(kind=8),intent(in)::tgmin
!!!!remember to put in PBC
Do i=1,N*N
Do k=1,2
!print*,i,x(i,k),v(i,k),tgmin,"old"
x(i,k) = x(i,k) + v(i,k)*tgmin
!print*,i,x(i,k),"new"
!print*,x(i,k),N,a,floor(x(i,k)/(N*a)),"the innards"
x(i,k) = x(i,k) - int(x(i,k)/(N*a))*N*a !! ---> everything back in the box
!print*,x(i,k),"in box"
Enddo
Enddo
end subroutine evolve

Subroutine collide(i,j)
implicit none
integer,intent(in)::i,j
!print*,i,j,"impostors??"
!print*,x(i,1),x(j,1),"where is capt. picard??"
bij=0.0
rijsq=0.0
rij=0.0
spam_r = 0.0
Do k=1,2
  spam_r = x(i,k)-x(j,k)
  !print*,spam_r,floor(spam_r/(N*a)),int(spam_r/(N*a)),anint(spam_r/(N*a)),"base, floor, int and anint"
  !print*,spam_r,"old spam_r"
  spam_r = spam_r - anint(spam_r/(N*a))*N*a !!! --->>> Put in PBC
  !print*,spam_r,spam_r**2,"should they be big??"
  bij = bij + spam_r*(v(i,k)-v(j,k))
  rijsq = rijsq + (spam_r)**2
Enddo
rij = sqrt(rijsq)
!print*,i,j,bij,rij,"are these biggg?"
!print*,bij*rij/(sigma**2),"speechless"
spam_x = 0
Do k=1,2
!print*,v(i,k),v(j,k),"old vels"
spam_x = x(i,k)-x(j,k)
spam_x = spam_x - anint(spam_x/(N*a))*N*a
v(i,k) = v(i,k) - bij*spam_x/(sigma**2)
v(j,k) = v(j,k) + bij*spam_x/(sigma**2)
!print*,v(i,k),v(j,k),"new vels"
Enddo
End subroutine

Subroutine init()
x(:,:)=0
v(:,:)=0
Do i=1,N
Do j = 1,N
P = j + (i-1)*N
x(p,1) = (j-1)*a
!print*,p,x(p,1),j*a,"vei vei wong"
x(p,2) = i*a
!print*,p,x(p,2),i*a,"so vei wong"
v(p,1) = 0.001*(rand()-0.5)
v(p,2) = ((-1)**i)*vyy + 0.001*(rand()-0.5)
Enddo
Enddo
!print*,x(:,1),x(:,2),"init posn"
!print*,v(:,1),v(:,2),"init velocities"
End subroutine

subroutine calc_hFunc()
hfunc=0.0
energy=0.0
enMax=-1
do i=1,N*N
ke(i)=0.5*(v(i,1)**2 + v(i,2)**2)
if(ke(i)>enMax)then
enMax=ke(i)
endif
enddo

! reinit
enBin=enMax/NenBin
do i=1,NenBin
pke(i)=0.0
enddo
denom=0.0

! KE histogram
do i=1,N*N
spam=anint(ke(i)/enBin)+1
if((spam>0).and.(spam<NenBin+1))then
pke(spam)=pke(spam)+1.0
denom=denom+1.0
endif
enddo

if(denom/=N*N)then
print*,"denom coming",denom,N*N
endif

do i=1,NenBin
pke(i)=pke(i)/denom
if(pke(i)>0.0)then
hfunc=hfunc+pke(i)*(log(pke(i)/sqrt((i-0.5)*enBin)) -1)
endif
enddo



end subroutine


End program eventful






