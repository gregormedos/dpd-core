module parameters
  !ensemble parameters
  integer :: N                  !number of particles
  real*8 :: dens                !density
  real*8 :: temp                !temperature
  !simulation parameters
  real*8 :: lbox                !length of box
  real*8 :: a1,a2,b1,b2,c1      !force field parameters
  integer :: cutoff             !flag to switch on or off potential cutoff
  real*8 :: rcut                !cutoff distance for potential
  integer :: ran                !random seed
  integer :: nseries            !number of series of samples
  integer :: nsampl             !number of time steps for sampling
  integer :: f_snap,f_sam,f_frame     !flags to switch on or off snapshots, energies, frames
  integer :: nsam               !number of time steps per energy saved
  integer :: nframes            !number of frames saved
  integer :: nframe             !number of time steps per frame saved
  !DPD parameters
  integer :: nesteps            !number of times steps for equilibration
  integer :: ntsteps            !number of time steps for sampling
  real*8 :: tstep,dt2,d2t2      !time step, half step, half squared step
  !number Pi
  real*8 :: Pi
  PARAMETER (Pi=4.0d0*datan(1.0d0))
endmodule parameters

module particles
  real*8,dimension(:,:),allocatable :: pos               !matrix of position vectors
  real*8,dimension(:,:),allocatable :: vel               !matrix of velocity vectors
  real*8,dimension(:,:),allocatable :: acc               !matrix of acceleration vectors
endmodule particles

module potential
  real*8,dimension(:,:),allocatable :: pot               !matrix of potential energies
  real*8,dimension(:,:),allocatable :: tla               !matrix of virial pressures
  real*8,dimension(:,:,:),allocatable :: sil             !matrix of forces
  real*8,dimension(:,:),allocatable :: dis               !matrix of distances
  real*8,dimension(:,:),allocatable :: dis2              !matrix of distances squared
  real*8,dimension(:,:,:),allocatable :: rrr             !matrix of distance vectors
endmodule potential

module microstate
  real*8 :: ekin                    !kinetic energy
  real*8 :: tem                     !equipartion theorem temperature
  real*8 :: u                       !excess internal energy
  real*8 :: etot                    !total internal energy
  real*8 :: p                       !virial theorem pressure
endmodule microstate

module macrostate
  real*8 :: eekin                   !average kinetic energy
  real*8 :: ttem                    !average temperature
  real*8 :: uu,uu2                  !average excess internal energy
  real*8 :: eetot                   !average internal energy
  real*8 :: pp                      !average pressure
  integer :: nsampls                !number of samples
endmodule macrostate

module correlation
  real*8 :: grinterv                            !g(r) radial interval
  real*8,dimension(:,:),allocatable :: gr       !pair correlation function g(r) samples
  real*8,dimension(:),allocatable :: vt,shear,flux
  real*8,dimension(:,:),allocatable :: vel0
  real*8,dimension(:),allocatable :: shear0
  real*8,dimension(:),allocatable :: flux0
endmodule correlation


program core
  use parameters
  use particles
  use potential
  use macrostate
  use correlation
  implicit none
  character :: tekst
  character(8) :: fmt
  integer :: i,j
  real*8 :: cv                                  !excess heat capacity
  real*8 :: c                                   !numerical constant
  real*8 :: stru                                !translational order parameter
  real*8,dimension(:),allocatable :: vvt,vvt2,ssh,ssh2,ffl,ffl2
  real*8 :: dif                                 !diffusion coefficient
  real*8 :: vis                                 !viscosity
  real*8 :: thec                                !thermal conductivity
  real*8,dimension(:),allocatable :: ggr        !average g(r)
  real*8,dimension(:),allocatable :: dgr        !estimated error of g(r)

  open(1001,file='dpd_param')
  read(1001,*)tekst
  read(1001,*)tekst,N
  read(1001,*)tekst,dens
  read(1001,*)tekst,temp
  read(1001,*)tekst
  read(1001,*)tekst,a1,a2,b1,b2
  read(1001,*)tekst,cutoff
  read(1001,*)tekst,rcut
  read(1001,*)tekst,ran
  read(1001,*)tekst,nseries
  read(1001,*)tekst,nsampl
  read(1001,*)tekst,f_snap
  read(1001,*)tekst,f_sam
  read(1001,*)tekst,nsam
  read(1001,*)tekst,f_frame
  read(1001,*)tekst,nframes
  read(1001,*)tekst,nframe
  read(1001,*)tekst
  read(1001,*)tekst,nesteps
  read(1001,*)tekst,ntsteps
  read(1001,*)tekst,tstep
  read(1001,*)tekst
  read(1001,*)tekst,grinterv
  close(1001)

  allocate( pos(N,2) )
  allocate( vel(N,2) )
  allocate( acc(N,2) )
  allocate( pot(N,N) )
  allocate( tla(N,N) )
  allocate( sil(N,N,2) )
  allocate( dis(N,N) )
  allocate( dis2(N,N) )
  allocate( rrr(N,N,2) )
  allocate( gr(nseries,10000) )
  allocate( ggr(10000) )
  allocate( dgr(10000) )
  allocate( vt(0:ntsteps) )
  allocate( shear(0:ntsteps) )
  allocate( flux(0:ntsteps) )
  allocate( vel0(N,2) )
  allocate( shear0(N) )
  allocate( flux0(N) )
  allocate( vvt(0:ntsteps) )
  allocate( vvt2(0:ntsteps) )
  allocate( ssh(0:ntsteps) )
  allocate( ssh2(0:ntsteps) )
  allocate( ffl(0:ntsteps) )
  allocate( ffl2(0:ntsteps) )
  lbox=dsqrt(dble(N)/dens)
  dt2=tstep/2.0d0
  d2t2=tstep**2/2.0d0
  c1=b1**2/2.0d0/temp !Fluctutaion-dissipation theorem (thermostat is built into the Langevin stochastic dynamics)

  open(1002,file='dpd.log')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'N=',N
  write(1002,*)'V=',lbox**2
  write(1002,*)'T=',temp
  write(1002,*)'density=',dens
  write(1002,*)'length of box=',lbox
  write(1002,*)'random seed=',ran
  write(1002,*)'time step=',tstep
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'Dissipative Particle Dynamics (DPD) simulation of 2D Lennard-Jones (LJ) clusters'
  write(1002,*)'-------------------------------------------------------------------------'


  !----------------------------------------------------------
  !DPD
  !----------------------------------------------------------
  call random_pos()
  if (f_snap.eq.1) call snapshot('sta')
  write(1002,*)'clusters randomly put in box'
  write(1002,*)'-------------------------------------------------------------------------'

  call init_pos()
  if (f_snap.eq.1) call snapshot('int')
  write(1002,*)'positions initialized'
  write(1002,*)'-------------------------------------------------------------------------'

  call init_vel()
  write(1002,*)'velocities initialized'
  write(1002,*)'-------------------------------------------------------------------------'

  call equil()
  if (f_snap.eq.1) call snapshot('ekv')
  write(1002,*)'DPD equilibration successful'
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'uu/N=',uu/dble(N)
  write(1002,*)'pp=',pp
  write(1002,*)'eekin=',eekin
  write(1002,*)'ttem=',ttem
  write(1002,*)'eetot=',eetot
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'DPD simulation initialized'
  write(1002,*)'-------------------------------------------------------------------------'

  open(223,file='dpd_td.dat')
  write(223,'(10a16)')'temperature','density','ex energy/N','ex energy**2/N','ex heat cap','pressure', &
  & 'trans order','diffusion','viscosity','therm cond'
  !correlation module
  gr=0.0d0
  vvt=0.0d0
  vvt2=0.0d0
  ssh=0.0d0
  ssh2=0.0d0
  ffl=0.0d0
  ffl2=0.0d0
  do j=1,nseries
    call dpdseries(j)
    if (f_snap.eq.1) then
      write(fmt,'(i3.3)')j
      call snapshot(trim(fmt))
    endif
    cv=(uu2-uu**2)/temp**2/dble(N)
    do i=1,10000
      c=Pi*dble(2*i-1)*grinterv**2
      gr(j,i)=gr(j,i)/dble(N-1)/dble(nsampls)/dens/c
    enddo
    stru=0.0d0
    do i=1,10000
      if ( (dble(i+1)-0.5d0)*grinterv .lt. lbox/2.0d0 ) stru=stru+(dabs(gr(j,i)-1.0d0)+dabs(gr(j,i+1)-1.0d0))*grinterv/2.0d0
    enddo
    stru=stru*dsqrt(dens)
    do i=0,ntsteps
      vvt(i)=vvt(i)+vt(i)
      vvt2(i)=vvt2(i)+vt(i)**2
      ssh(i)=ssh(i)+shear(i)
      ssh2(i)=ssh2(i)+shear(i)**2
      ffl(i)=ffl(i)+flux(i)
      ffl2(i)=ffl2(i)+flux(i)**2
    enddo
    dif=0.0d0
    vis=0.0d0
    thec=0.0d0
    do i=1,ntsteps
      dif=dif+(vt(i)+vt(i-1))*dt2
      vis=vis+(shear(i)+shear(i-1))*dt2
      thec=thec+(flux(i)+flux(i-1))*dt2
    enddo
    dif=dif/4.0d0
    vis=vis*dens/temp
    thec=thec*dens/temp**2
    write(223,'(10e16.7)')temp,dens,uu/dble(N),uu2/dble(N),cv,pp,stru,dif,vis,thec
    write(1002,*) j,'-th DPD run successful'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'uu2/N=',uu2/dble(N)
    write(1002,*)'cv=',cv
    write(1002,*)'pp=',pp
    write(1002,*)'stru=',stru
    write(1002,*)'dif=',dif
    write(1002,*)'vis=',vis
    write(1002,*)'thec=',thec
    write(1002,*)'eekin=',eekin
    write(1002,*)'ttem=',ttem
    write(1002,*)'eetot=',eetot
    write(1002,*)'-------------------------------------------------------------------------'
  enddo
  close(223)

  open(202,file='dpd_corr.dat')
  ggr=0.0d0
  dgr=0.0d0
  do i=1,10000
    do j=1,nseries
      ggr(i)=ggr(i)+gr(j,i)
    enddo
    ggr(i)=ggr(i)/dble(nseries)
    do j=1,nseries
      dgr(i)=dgr(i)+(gr(j,i)-ggr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dble(nseries))
    if ( (dble(i)-0.5d0)*grinterv .lt. lbox ) write(202,'(3e16.7)')(dble(i)-0.5d0)*grinterv,ggr(i),dgr(i)
  enddo
  close(202)

  open(202,file='dpd_acorr.dat')
  vvt=vvt/dble(nseries)
  vvt2=vvt2/dble(nseries)
  ssh=ssh/dble(nseries)
  ssh2=ssh2/dble(nseries)
  ffl=ffl/dble(nseries)
  ffl2=ffl2/dble(nseries)
  do i=1,ntsteps
    write(202,'(13e16.7)')dble(i)*tstep, &
    & vvt(i),dsqrt(vvt2(i)-vvt(i)**2),vvt(i)/vvt(0), &
    & ssh(i),dsqrt(ssh2(i)-ssh(i)**2),ssh(i)/ssh(0), &
    & ffl(i),dsqrt(ffl2(i)-ffl(i)**2),ffl(i)/ffl(0)
  enddo
  close(202)

  write(1002,*)'DONE'
  write(1002,*)'-------------------------------------------------------------------------'
  close(1002)

  deallocate( pos )
  deallocate( vel )
  deallocate( acc )
  deallocate( pot )
  deallocate( tla )
  deallocate( sil )
  deallocate( dis )
  deallocate( dis2 )
  deallocate( rrr )
  deallocate( gr )
  deallocate( ggr )
  deallocate( dgr )

  open(1002,file='dpd_DONE')
  close(1002)
endprogram core

!----------------------------------------------------------------------------------------------------
!Generates random initial positions for clusters inside the box without overlapping
!----------------------------------------------------------------------------------------------------
subroutine random_pos()
  use parameters
  use particles
  implicit none
  real*8 :: ran3,image !functions
  real*8 :: a
  real*8 :: x,y,r2
  integer :: i,j,k

  !first particle
  a=ran3(ran)
  pos(1,1)=(a-0.5d0)*lbox
  a=ran3(ran)
  pos(1,2)=(a-0.5d0)*lbox
  
  !prevent overlapping
  do i=2,N
    j=0
    do while (j.lt.0.5)
      a=ran3(ran)
      pos(i,1)=(a-0.5d0)*lbox
      a=ran3(ran)
      pos(i,2)=(a-0.5d0)*lbox
      j=1
      do k=1,i-1
        x=image(pos(i,1),pos(k,1),lbox)
        y=image(pos(i,2),pos(k,2),lbox)
        r2=x**2+y**2
        if (r2.lt.0.4d0) j=0
      enddo
    enddo
  enddo
endsubroutine random_pos

!----------------------------------------------------------------------------------------------------
!Ran3 random number generator generates real numbers [0,1)
!----------------------------------------------------------------------------------------------------
function ran3(idum)
  implicit none
  integer :: idum
  integer :: MBIG,MSEED,MZ
  real*8 :: ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.0d0/dble(MBIG))
  integer :: i,iff,ii,inext,inextp,k
  integer :: mj,mk,ma(55)
  SAVE iff,inext,inextp,ma
  DATA iff /0/
  
  if (idum.lt.0.or.iff.eq.0) then
    iff=1
    mj=MSEED-iabs(idum)
    mj=mod(mj,MBIG)
    ma(55)=mj
    mk=1
    do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if (mk.lt.MZ) mk=mk+MBIG
      mj=ma(ii)
    enddo
    do k=1,4
      do i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if (ma(i).lt.MZ) ma(i)=ma(i)+MBIG
      enddo
    enddo
    inext=0
    inextp=31
    idum=1
  endif
  inext=inext+1
  if (inext.eq.56) inext=1
  inextp=inextp+1
  if (inextp.eq.56) inextp=1
  mj=ma(inext)-ma(inextp)
  if (mj.lt.MZ) mj=mj+MBIG
  ma(inext)=mj
  ran3=mj*FAC
  return
endfunction ran3

!----------------------------------------------------------------------------------------------------
!Minimum image convention for periodic boundary conditions
!----------------------------------------------------------------------------------------------------
function image(xa,xb,ll)
  implicit none
  real*8 :: xa,xb,ll
  real*8 :: xc
  real*8 :: image

  xc=xb-xa
  xc=xc-ll*dnint(xc/ll)      !rounds double to nearest double
  image=xc
  return
endfunction image

!----------------------------------------------------------------------------------------------------
!Calculates interactions between all clusters
!----------------------------------------------------------------------------------------------------
subroutine interactions()
  use parameters
  use particles
  use potential
  use microstate
  implicit none
  real*8 :: image,ran3,force,potent !functions
  real*8 :: x,y,r,r2
  real*8 :: vx,vy,v
  real*8 :: f
  real*8 :: x1,x2,s
  integer :: i,j

  !potential module
  pot=0.0d0
  tla=0.0d0
  sil=0.0d0
  dis=0.0d0
  dis2=0.0d0
  rrr=0.0d0

  !microstate module
  u=0.0d0
  p=0.0d0

  !interactions
  do i=1,N-1
    do j=i+1,N
      x=image(pos(j,1),pos(i,1),lbox)
      y=image(pos(j,2),pos(i,2),lbox)
      r2=x**2+y**2
      r=dsqrt(r2)
      if ( (cutoff.ne.1).or.(r.lt.rcut) ) then
        vx=vel(i,1)-vel(j,1)
        vy=vel(i,2)-vel(j,2)
        v=(vx*x+vy*y)/r
        !Box-Muller transform with Marsaglia polar method
        x1=0.0d0
        x2=0.0d0
        s=2.0d0
        do while ((s.eq.0.0d0).or.(s.ge.1.0d0))
          x1=2.0d0*ran3(ran)-1.0d0
          x2=2.0d0*ran3(ran)-1.0d0
          s=x1**2+x2**2
        enddo
        s=x1*dsqrt(-2.0d0*dlog(s)/s)
        !potential module
        f=force(r2,r,v,a1,a2,b1,b2,c1,s,potent)
        pot(i,j)=potent
        tla(i,j)=f*r
        sil(i,j,1)=f*x/r
        sil(i,j,2)=f*y/r
        dis(i,j)=r
        dis2(i,j)=r2
        rrr(i,j,1)=x
        rrr(i,j,2)=y
        pot(j,i)=pot(i,j)
        tla(j,i)=tla(i,j)
        sil(j,i,1)=-sil(i,j,1)
        sil(j,i,2)=-sil(i,j,2)
        dis(j,i)=r
        dis2(j,i)=r2
        rrr(j,i,1)=-x
        rrr(j,i,2)=-y
        !microstate module
        u=u+pot(i,j)
        p=p+tla(i,j)
      else
        !potential moduleclusters
        dis(i,j)=r
        dis2(i,j)=r2
        rrr(i,j,1)=x
        rrr(i,j,2)=y
        dis(j,i)=r
        dis2(j,i)=r2
        rrr(j,i,1)=-x
        rrr(j,i,2)=-y
      endif
    enddo
  enddo

  !particles module
  do i=1,N
    acc(i,1)=sum(sil(i,:,1))
    acc(i,2)=sum(sil(i,:,2))
  enddo
endsubroutine interactions

!----------------------------------------------------------------------------------------------------
!Calculates the force between two clusters
!----------------------------------------------------------------------------------------------------
function force(r2,r,v,a1,a2,b1,b2,c1,s,potent)
  implicit none
  !using reduced units
  real*8 :: r2,r,v !distance squared, distance, projection of the relative velocity onto the direction between the clusters
  real*8 :: a1,a2,b1,b2,c1
  real*8 :: s !random Gaussian number with unit variance
  real*8 :: w !weight function
  real*8 :: x
  real*8 :: force,potent

  if (r.lt.b2) then
    w=1.0d0-r/b2
  else
    w=0.0d0
  endif
  x=a2**2/r2
  x=x**3
  force=48.0d0*a1*x*(x-0.5d0)/r+b1*w*s-c1*w**2*v
  potent=4.0d0*a1*x*(x-1.0d0)
  return
endfunction force

!---------------------------------------------------------------------------------------------
!Initializes positions for a homogeneous fluid
!---------------------------------------------------------------------------------------------
subroutine init_pos()
  use parameters
  use particles
  implicit none
  integer :: i
  real*8 :: image !function
  real*8 :: sumpos(2)

  !setting the sum of positions to zero
  sumpos=1.0d0
  do while (sum(sumpos*sumpos).gt.1.0d-5)
    sumpos(1)=sum(pos(:,1))/dble(N)
    sumpos(2)=sum(pos(:,2))/dble(N)
    do i=1,N
      pos(i,1)=image(sumpos(1),pos(i,1),lbox)
      pos(i,2)=image(sumpos(2),pos(i,2),lbox)
    enddo
  enddo
endsubroutine init_pos

!-------------------------------------------------------------------------------------------------------------------------------
!Initializes velocities for an isotropic fluid
!Box-Muller transform with Marsaglia polar method for Maxwell-Boltzmann distribution
!-------------------------------------------------------------------------------------------------------------------------------
subroutine init_vel()
  use parameters
  use particles
  use microstate
  implicit none
  integer :: i
  real*8 :: ran3 !function
  real*8 :: sigma,x1,x2,s,z1,z2
  real*8 :: sumvel(2)
  real*8 :: rescale

  !Box-Muller transform with Marsaglia polar method
  x1=0.0d0
  x2=0.0d0
  sigma=dsqrt(temp)
  do i=1,N
    s=2.0d0
    do while ((s.eq.0.0d0).or.(s.ge.1.0d0))
      x1=2.0d0*ran3(ran)-1.0d0
      x2=2.0d0*ran3(ran)-1.0d0
      s=x1**2+x2**2
    enddo
    z1=x1*dsqrt(-2.0d0*dlog(s)/s)
    z2=x2*dsqrt(-2.0d0*dlog(s)/s)
    vel(i,1)=sigma*z1
    vel(i,2)=sigma*z2
  enddo

  !setting the sum of velocities to zero
  sumvel=1.0d0
  do while (sum(sumvel*sumvel).gt.1.0d-5)
    sumvel(1)=sum(vel(:,1))/dble(N)
    sumvel(2)=sum(vel(:,2))/dble(N)
    vel(:,1)=vel(:,1)-sumvel(1)
    vel(:,2)=vel(:,2)-sumvel(2)
  enddo

  !calculate the temperature
  call temperature()

  !rescale the particle velocities to the desired temperature
  rescale=dsqrt(temp/tem)
  vel=vel*rescale
endsubroutine init_vel

!----------------------------------------------------------------------------------------------------
!Calculates current temperature using the equipartition theorem
!----------------------------------------------------------------------------------------------------
subroutine temperature()
  use parameters
  use particles
  use microstate
  implicit none
  integer :: i

  ekin=0.0d0
  do i=1,N
    ekin=ekin+(vel(i,1)**2+vel(i,2)**2)/2.0d0
  enddo
  !(x,y) -> 2 degrees of freedom
  !1/2 T for each degree of freedom
  !2N degrees of freedom (N particles with (x,y))
  tem=ekin/dble(N)
endsubroutine temperature

!-------------------------------------------------------------------------
!Velocity Verlet algorithm for calculating velocities
!Newton's equations of motion
!-------------------------------------------------------------------------
subroutine move()
  use parameters
  use particles
  use microstate
  implicit none
  integer :: i
  real*8 :: image !function

  do i=1,N
    !half-update the particle velocities by tstep/2
    vel(i,1)=vel(i,1)+acc(i,1)*dt2
    vel(i,2)=vel(i,2)+acc(i,2)*dt2
    !update the particle positions by tstep
    pos(i,1)=pos(i,1)+vel(i,1)*tstep
    pos(i,2)=pos(i,2)+vel(i,2)*tstep
    pos(i,1)=image(0.0d0,pos(i,1),lbox)
    pos(i,2)=image(0.0d0,pos(i,2),lbox)
  enddo

  !update the interactions
  call interactions()

  do i=1,N
    !update the particle velocities by tstep/2
    vel(i,1)=vel(i,1)+acc(i,1)*dt2
    vel(i,2)=vel(i,2)+acc(i,2)*dt2
  enddo
  
  !update the temperature
  call temperature()
endsubroutine move

!----------------------------------------------------------------------------------------------------
!DPD equilibration
!----------------------------------------------------------------------------------------------------
subroutine equil()
  use parameters
  use particles
  use microstate
  use macrostate
  implicit none
  integer :: i
  real*8 :: t

  call interactions()
  call temperature()
  t=0.0d0

  !macrostate module
  eekin=0.0d0
  ttem=0.0d0
  uu=0.0d0
  eetot=0.0d0
  pp=0.0d0
  nsampls=0

  !macrostate = average over microstates
  if (f_sam.eq.1) open(202,file='dpd_sam_ekv')
  do i=1,nesteps
    call move()
!rescale the particle velocities to the desired temperature
vel=vel*dsqrt(temp/tem)
    t=t+tstep
    if (mod(i,nsampl).eq.0) then
      etot=ekin+u
      if ( (f_sam.eq.1).and.(mod(i,nsam).eq.0) ) write(202,'(6e16.7)')t,ekin,tem,u,etot,p
      uu=uu+u
      pp=pp+p
      nsampls=nsampls+1
      eekin=eekin+ekin
      ttem=ttem+tem
      eetot=eetot+etot
    endif
  enddo
  if (f_sam.eq.1) close(202)

  !macrostate module
  uu=uu/dble(nsampls)
  pp=pp/dble(nsampls)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
  eekin=eekin/dble(nsampls)
  ttem=ttem/dble(nsampls)
  eetot=eetot/dble(nsampls)
endsubroutine equil

!----------------------------------------------------------------------------------------------------
!DPD sampling
!----------------------------------------------------------------------------------------------------
subroutine dpdseries(idum)
  use parameters
  use particles
  use potential
  use microstate
  use macrostate
  use correlation
  implicit none
  integer :: idum
  character(8) :: fmt
  integer :: i,j
  real*8 :: t

  call interactions()
  call temperature()
  t=0.0d0

  !macrostate module
  uu=0.0d0
  uu2=0.0d0
  pp=0.0d0
  eekin=0.0d0
  ttem=0.0d0
  eetot=0.0d0
  nsampls=0

  !correlation module
  vel0=vel
  do i=1,N
    shear0(i)=vel(i,1)*vel(i,2)+pos(i,1)*acc(i,2)
    flux0(i)=0.5d0*vel(i,1)*(vel(i,1)**2+vel(i,2)**2+sum(pot(:,i)))
    do j=1,N
      flux0(i)=flux0(i)+rrr(i,j,1)*(sil(i,j,1)*vel(i,1)+sil(i,j,2)*vel(i,2))
    enddo
  enddo
  call acorr(0)

  !macrostate = average over microstates
  if (f_sam.eq.1) then
    write(fmt,'(i3.3)')idum
    open(202,file='dpd_sam_'//trim(fmt))
  endif
  do i=1,ntsteps
    call move()
    t=t+tstep
    if (mod(i,nsampl).eq.0) then
      etot=ekin+u
      if ( (f_sam.eq.1).and.(mod(i,nsam).eq.0) ) write(202,'(6e16.7)')t,ekin,tem,u,etot,p
      uu=uu+u
      uu2=uu2+u**2
      pp=pp+p
      nsampls=nsampls+1
      eekin=eekin+ekin
      ttem=ttem+tem
      eetot=eetot+etot
      call corr(idum)
      call acorr(i)
    endif
    if ( (f_frame.eq.1).and.(idum.eq.1).and.(mod(i,nframe).eq.0).and.(i/nframe.lt.nframes+1) ) then
      write(fmt,'(i3.3,i5.5)')idum,i/nframe
      call movie(fmt)
    endif
  enddo
  if (f_sam.eq.1) close(202)

  !macrostate module
  uu=uu/dble(nsampls)
  pp=pp/dble(nsampls)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
  eekin=eekin/dble(nsampls)
  ttem=ttem/dble(nsampls)
  eetot=eetot/dble(nsampls)
endsubroutine dpdseries

!----------------------------------------------------------------------------------------------------
!Pair correlation function g(r) -- Radial distribution function
!----------------------------------------------------------------------------------------------------
subroutine corr(idum)
  use parameters
  use potential
  use correlation
  implicit none
  integer :: idum
  integer :: i,j
  integer :: k
  real*8 :: r

  do i=1,N-1
    do j=i+1,N
      r=dis(i,j)
      k=int(r/grinterv)        !floors double to integer
      if (k.lt.10000) then
        gr(idum,k+1)=gr(idum,k+1)+2.0d0
      endif
    enddo
  enddo
endsubroutine corr

!----------------------------------------------------------------------------------------------------
!Autocorrelation function v(t)v(0)
!----------------------------------------------------------------------------------------------------
subroutine acorr(i)
  use parameters
  use particles
  use potential
  use correlation
  implicit none
  integer :: i
  integer :: j,k

  do j=1,N
    vt(i)=vt(i)+vel(j,1)*vel0(j,1)+vel(j,2)*vel0(j,2)
    shear(i)=shear(i)+(vel(j,1)*vel(j,2)+pos(j,1)*acc(j,2))*shear0(j)
    flux(i)=flux(i)+0.50d0*vel(j,1)*(vel(j,1)**2+vel(j,2)**2+sum(pot(j,:)))*flux0(j)
    do k=1,N
      flux(i)=flux(i)+(rrr(j,k,1)*(sil(j,k,1)*vel(j,1)+sil(j,k,2)*vel(j,2)))*flux0(j)
    enddo
  enddo

  !correlation module
  vt(i)=vt(i)/dble(N)
  shear(i)=shear(i)/dble(N)
  flux(i)=flux(i)/dble(N)
  end subroutine acorr

!-------------------------------------------------------
!Snapshot
!-------------------------------------------------------
subroutine snapshot(fmt)
  use parameters
  use particles
  implicit none
  character(3) :: fmt
  integer :: i

  open(111,file='dpd_snap_'//fmt)
  do i=1,N
    write(111,'(3e16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine snapshot

!--------------------------------------------------------
!Movie
!--------------------------------------------------------
subroutine movie(fmt)
  use parameters
  use particles
  implicit none
  character(8) :: fmt
  integer :: i

  open(111,file='dpd_frame_'//fmt)
  do i=1,N
    write(111,'(3e16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine movie
