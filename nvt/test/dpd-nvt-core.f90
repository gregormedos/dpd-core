module parameters
  !ensemble parameters
  integer :: N                  !number of particles (N)
  real*8 :: dens                !density
  real*8 :: temp                !temperature
  !simulation parameters
  real*8 :: lbox                !length of box
  real*8 :: a1,a2               !LJ potential parameters
  integer :: cutoff             !flag to switch on or off potential cutoff
  real*8 :: rcut                !cutoff distance for potential
  integer :: ran                !random seed
  integer :: nseries            !number of series of sampling (MAX 20)
  integer :: snap,sam,frame     !flags to switch on or off snapshots, energies, frames
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
  real*8, dimension(:,:), allocatable :: pos               !matrix of position vectors
  real*8, dimension(:,:), allocatable :: vel               !matrix of velocity vectors
  real*8, dimension(:,:), allocatable :: acc               !matrix of acceleration vectors
endmodule particles

module potential
  real*8, dimension(:,:), allocatable :: pot               !matrix of potential energies
  real*8, dimension(:,:), allocatable :: tla               !matrix of virial pressures
  real*8, dimension(:,:,:), allocatable :: sil             !matrix of forces
  real*8, dimension(:,:), allocatable :: dis               !matrix of distances
  real*8, dimension(:,:), allocatable :: dis2              !matrix of distances squared
  real*8, dimension(:,:,:), allocatable :: rrr             !matrix of distance vectors
endmodule potential

module energy
  real*8 :: ekin                    !kinetic energy
  real*8 :: tem                     !equipartion theorem temperature
  real*8 :: u,uu,uu2                !potential energy, excess internal energy, excess internal energy squared
  real*8 :: etot                    !total internal energy
  real*8 :: p,pp                    !virial pressure
  integer :: m                      !number of samples
endmodule energy

module correlation
  real*8, dimension(10000) :: gr                !correlation function g(r)
  real*8 :: grinterv                            !g(r) radial interval
  integer :: ncorr                              !cycles/steps per sampling
  integer :: mcorr                              !number of samples
endmodule correlation


program core
  use parameters
  use particles
  use potential
  use energy
  use correlation
  implicit none
  real*8 :: time1,time2                         !CPU_time
  integer, dimension(8) :: values1,values2      !value(1)=year, value(2)=month, value(3)=day, value(4)=time difference with UTC in minutes, value(5)=hour, value(6)=minute, value(7)=second, value(8)=milisecond
  character :: tekst
  character(8) :: fmt
  integer :: i,j
  real*8 :: cv                                  !heat capacity
  real*8 :: dgr(10000),ggr(10000,20)            !correlation function

  call cpu_time(time1)
  call date_and_time(VALUES=values1)
  open(1002,file='dpd.log')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'start date and time YYYY, MM, DD, UTC, HH, MIN, SS'
  write(1002,'(8i5)')values1

  open(1001,file='param')
  read(1001,*)tekst
  read(1001,*)tekst,N
  read(1001,*)tekst,dens
  read(1001,*)tekst,temp
  read(1001,*)tekst
  read(1001,*)tekst,a1,a2
  read(1001,*)tekst,cutoff
  read(1001,*)tekst,rcut
  read(1001,*)tekst,ran
  read(1001,*)tekst,nseries
  read(1001,*)tekst,snap
  read(1001,*)tekst,sam
  read(1001,*)tekst,nsam
  read(1001,*)tekst,frame
  read(1001,*)tekst,nframes
  read(1001,*)tekst,nframe
  read(1001,*)tekst
  read(1001,*)tekst,nesteps
  read(1001,*)tekst,ntsteps
  read(1001,*)tekst,tstep
  read(1001,*)tekst
  read(1001,*)tekst,grinterv
  read(1001,*)tekst,ncorr
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
  lbox=dsqrt(dble(N)/dens)
  dt2=tstep/2.0d0
  d2t2=tstep**2/2.0d0

  write(1002,*)'N=',N
  write(1002,*)'V=',lbox**2
  write(1002,*)'T=',temp
  write(1002,*)'density=',dens
  write(1002,*)'length of box=',lbox
  write(1002,*)'random seed=',ran
  write(1002,*)'time step=',tstep
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'Dissipative Particle Dynamics (DPD) simulation of 2D Lennard-Jones (LJ) clusters'

  !----------------------------------------------------------
  !DPD
  !----------------------------------------------------------
  call random_pos()
  if (snap.eq.1) call snapshot('sta')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'particles randomly put in box'

  call init_pos()
  if (snap.eq.1) call snapshot('int')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'positions initialized'

  call init_vel()
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'velocities initialized'

  call dpdequil()
  if (snap.eq.1) call snapshot('ekv')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'DPD equilibration successful'
  write(1002,*)'uu/N=',uu/dble(N)
  write(1002,*)'pp=',pp
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'DPD simulation initialized'

  open(223,file='dpd-td.dat')
  write(223,'(4a16)')'temperature','density','ex energy/N','ex heat cap','pressure'
  do j=1,nseries
    call dpdseries(j)
    write(fmt,'(i3.3)')j
    if (snap.eq.1) call snapshot(trim(fmt))
    cv=(uu2-uu**2)/temp**2
    ggr(:,j)=gr(:)
    write(223,'(5e16.7)')temp,dens,uu/dble(N),cv,pp
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*) j,'-th DPD run successful'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'cv=',cv
    write(1002,*)'pp=',pp
  enddo
  close(223)

  open(202,file='dpd-gr.dat')
  gr=0.0d0
  dgr=0.0d0
  do i=1,10000
    do j=1,nseries
      gr(i)=gr(i)+ggr(i,j)
    enddo
    gr(i)=gr(i)/dble(nseries)
    do j=1,nseries
      dgr(i)=dgr(i)+(ggr(i,j)-gr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dble(nseries))
    write(202,'(3f16.7)')(dble(i)-0.5d0)*grinterv,gr(i),dgr(i)
  enddo
  close(202)

  call cpu_time(time2)
  call date_and_time(VALUES=values2)
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'CPU simulation time=',time2-time1
  write(1002,*)'start and finish date and time YYYY, MM, DD, UTC, HH, MIN, SS'
  write(1002,'(8i5)')values1
  write(1002,'(8i5)')values2
  time2=(values2(8)-values1(8))/1000.0d0+values2(7)-values1(7)
  time2=time2+(values2(6)-values1(6))*60.0d0
  time2=time2+(values2(5)-values1(5))*60.0d0*60.0d0
  time2=time2+(values2(3)-values1(3))*60.0d0*60.0d0*24.0d0
  write(1002,*)'real simulation time=',time2
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

  open(1002,file='dpd-DONE')
  close(1002)
endprogram core

!----------------------------------------------------------------------------------------------------
!Generates random initial positions for particles inside the box without overlapping
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
!Calculates interactions between all particles
!----------------------------------------------------------------------------------------------------
subroutine interactions()
  use parameters
  use particles
  use potential
  use energy
  implicit none
  real*8 :: image,ljpot,tl !functions
  real*8 :: x,y,r,r2
  integer :: i,j

  !potential module
  pot=0.0d0
  tla=0.0d0
  sil=0.0d0
  dis=0.0d0
  dis2=0.0d0
  rrr=0.0d0

  !energy module
  u=0.0d0
  p=0.0d0

  !interactions
  do i=1,N-1
    do j=i+1,N
      x=image(pos(i,1),pos(j,1),lbox)
      y=image(pos(i,2),pos(j,2),lbox)
      r2=x**2+y**2
      r=dsqrt(r2)
      if ( (cutoff.eq.0).or.(r.lt.rcut) ) then
        !potential module
        pot(i,j)=ljpot(r2,a1,a2,tl)
        tla(i,j)=tl
        sil(i,j,1)=tl*x/r2
        sil(i,j,2)=tl*y/r2
        sil(j,i,1)=-sil(i,j,1)
        sil(j,i,2)=-sil(i,j,2)
        dis(i,j)=r
        dis2(i,j)=r2
        rrr(i,j,1)=x
        rrr(i,j,2)=y
        pot(j,i)=pot(i,j)
        tla(j,i)=tl
        dis(j,i)=r
        dis2(j,i)=r2
        rrr(j,i,1)=-x
        rrr(j,i,2)=-y
        !energy module
        u=u+pot(i,j)
        p=p+tla(i,j)
      else
        !potential module
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
  do j=1,N
    acc(j,1)=sum(sil(:,j,1))
    acc(j,2)=sum(sil(:,j,2))
  enddo
endsubroutine interactions

!----------------------------------------------------------------------------------------------------
!Calculates the CS potential between two particles
!----------------------------------------------------------------------------------------------------
function ljpot(r2,a1,a2,p)
  implicit none
  !using reduced units
  real*8 :: r2
  real*8 :: a1,a2
  !a1=epsilon = absolute value of the minimum value of the LJ potential (depth of the potential well)
  !a2=sigma = distance at which the potential becomes positive
  real*8 :: x
  real*8 :: ljpot,p !CS potential, virial pressure

  x=a2**2/r2
  x=x**3
  ljpot=4.0d0*a1*x*(x-1.0d0)
  p=48.0d0*a1*x*(x-0.5d0)        ! -dU/dr * r = F * r
  return
endfunction ljpot

!----------------------------------------------------------------------------------------------------
!Pair correlation function g(r) -- Radial distribution function
!----------------------------------------------------------------------------------------------------
subroutine corr()
  use parameters
  use potential
  use correlation
  implicit none
  integer :: i,j
  integer :: k
  real*8 :: r

  do i=1,N-1
    do j=i+1,N
      r=dis(i,j)
      k=int(r/grinterv)        !floors double to integer
      if (k.lt.10000) then
        gr(k+1)=gr(k+1)+2.0d0
      endif
    enddo
  enddo
endsubroutine corr

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
!Box-Muller transform for Maxwell-Boltzmann distribution
!-------------------------------------------------------------------------------------------------------------------------------
subroutine init_vel()
  use parameters
  use particles
  use energy
  implicit none
  integer :: i
  real*8 :: ran3 !function
  real*8 :: sigma,x1,x2,s,z1,z2
  real*8 :: sumvel(2)
  real*8 :: rescale

  !Box-Muller transform
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
  use energy
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
!Rescaling thermostat
!Newton's equations of motion
!-------------------------------------------------------------------------
subroutine move()
  use parameters
  use particles
  use energy
  implicit none
  integer :: i
  real*8 :: image !function
  real*8 :: rescale

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

  !rescale the particle velocities to the desired temperature
  rescale=dsqrt(temp/tem)
  vel=vel*rescale
endsubroutine move

!----------------------------------------------------------------------------------------------------
!DPD equilibration
!----------------------------------------------------------------------------------------------------
subroutine dpdequil()
  use parameters
  use particles
  use energy
  implicit none
  integer :: i
  real*8 :: t

  t=0.0d0
  uu=0.0d0
  pp=0.0d0
  m=0

  open(202,file='dpd-ekv.sam')
  call interactions()
  call temperature()
  do i=1,nesteps
    call move()
    t=t+tstep
    etot=ekin+u
    if ( (sam.eq.1).and.(mod(i,nsam).eq.0) ) write(202,'(6f16.7)')t,u,ekin,etot,p,tem
    uu=uu+u
    pp=pp+p
    m=m+1
  enddo
  close(202)

  uu=uu/dble(m)
  pp=pp/dble(m)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
endsubroutine dpdequil

!----------------------------------------------------------------------------------------------------
!DPD sampling
!----------------------------------------------------------------------------------------------------
subroutine dpdseries(idum)
  use parameters
  use particles
  use energy
  use correlation
  implicit none
  integer :: idum
  character(8) :: fmt
  real*8 :: c                                   !numerical constant
  integer :: i
  real*8 :: t

  t=0.0d0
  uu=0.0d0
  uu2=0.0d0
  pp=0.0d0
  m=0

  gr=0.0d0
  mcorr=0

  write(fmt,'(i3.3)')idum
  if (sam.eq.1) open(202,file='dpd-results'//trim(fmt)//'.sam')
  call interactions()
  call temperature()
  do i=1,ntsteps
    call move()
    t=t+tstep
    etot=ekin+u
    if ( (sam.eq.1).and.(mod(i,nsam).eq.0) ) write(202,'(6f16.7)')t,u,ekin,etot,p,tem
    uu=uu+u
    uu2=uu2+u**2
    pp=pp+p
    m=m+1
    if (mod(i,ncorr).eq.0) then
      call corr()
      mcorr=mcorr+1
    endif
    if ( (frame.eq.1).and.(idum.eq.1).and.(mod(i,nframe).eq.0).and.(i/nframe.lt.(nframes+1)) ) then
      write(fmt,'(i3.3,i5.5)')idum,i/nframe
      call movie(fmt)
    endif
  enddo
  if (sam.eq.1) close(202)

  uu=uu/dble(m)
  uu2=uu2/dble(m)
  pp=pp/dble(m)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp

  do i=1,10000
    c=Pi*dble(2*i-1)*grinterv**2
    gr(i)=gr(i)/dble(N-1)/dble(mcorr)/dens/c
  enddo
endsubroutine dpdseries

!-------------------------------------------------------
!Snapshot
!-------------------------------------------------------
subroutine snapshot(fmt)
  use parameters
  use particles
  implicit none
  character(3) :: fmt
  integer :: i

  open(111,file='dpd-s'//fmt//'.snap')
  do i=1,N
    write(111,'(3f16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
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

  open(111,file='dpd-f'//fmt//'.frame')
  do i=1,N
    write(111,'(3f16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine movie