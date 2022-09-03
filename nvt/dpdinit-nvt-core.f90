module parameters
    !ensemble parameters
    integer :: N                  !number of particles
    integer :: NC                 !number of clusters
    real*8 :: dens                !particle density
    real*8 :: densc               !cluster density
    real*8 :: temp                !temperature
    !simulation parameters
    real*8 :: lbox                !length of box
    real*8 :: a1,a2               !potential parameters
    integer :: cutoff             !flag to switch on or off potential cutoff
    real*8 :: rcut                !cutoff distance for potential
    integer :: ran                !random seed
    integer :: nseries            !number of series of samples
    integer :: nsampl             !number of time steps for sampling
    integer :: f_snap,f_sam,f_frame     !flags to switch on or off snapshots, energies, frames
    integer :: nsam               !number of time steps per energy saved
    integer :: nframes            !number of frames saved
    integer :: nframe             !number of time steps per frame saved
    !MD parameters
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
    integer,dimension(:),allocatable :: clust              !each particle is assigned to a cluster index
  endmodule particles

  module clusters
    real*8,dimension(:,:),allocatable :: posc
    real*8,dimension(:,:),allocatable :: velc
    real*8,dimension(:,:),allocatable :: accc
    integer,dimension(:),allocatable :: part               !number of particles assigned to each cluster index
  endmodule clusters
  
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
    real*8 :: ekinc
    real*8 :: temc
    real*8 :: u                       !excess internal energy
    real*8 :: etot                    !total internal energy
    real*8 :: p                       !virial theorem pressure
  endmodule microstate
  
  module macrostate
    real*8 :: eekin                   !average kinetic energy
    real*8 :: ttem                    !average temperature
    real*8 :: uu                      !average excess internal energy
    real*8 :: eetot                   !average internal energy
    real*8 :: pp                      !average pressure
    integer :: nsampls                !number of samples
  endmodule macrostate
  
  module correlation
    real*8 :: grinterv                             !g(r) radial interval
    real*8,dimension(:,:),allocatable :: grc       !pair correlation function g(r) samples
  endmodule correlation
  
  
  program core
    use parameters
    use particles
    use clusters
    use potential
    use macrostate
    use correlation
    implicit none
    character :: tekst
    character(8) :: fmt
    integer :: i,j
    real*8 :: c                                    !numerical constant
    real*8,dimension(:),allocatable :: ggrc        !average g(r)
    real*8,dimension(:),allocatable :: dgrc        !estimated error of g(r)
  
    open(1001,file='dpdinit_param')
    read(1001,*)tekst
    read(1001,*)tekst,N
    read(1001,*)tekst,NC
    read(1001,*)tekst,dens
    read(1001,*)tekst,temp
    read(1001,*)tekst
    read(1001,*)tekst,a1,a2
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
    allocate( clust(N) )
    allocate( posc(NC,2) )
    allocate( velc(NC,2) )
    allocate( accc(NC,2) )
    allocate( part(NC) )
    allocate( pot(N,N) )
    allocate( tla(N,N) )
    allocate( sil(N,N,2) )
    allocate( dis(N,N) )
    allocate( dis2(N,N) )
    allocate( rrr(N,N,2) )
    allocate( grc(nseries,10000) )
    allocate( ggrc(10000) )
    allocate( dgrc(10000) )
    lbox=dsqrt(dble(N)/dens)
    densc=dens*dble(NC)/dble(N)
    dt2=tstep/2.0d0
    d2t2=tstep**2/2.0d0
  
    open(1002,file='dpdinit.log')
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*)'N=',N
    write(1002,*)'NC=',NC
    write(1002,*)'V=',lbox**2
    write(1002,*)'T=',temp
    write(1002,*)'density=',dens
    write(1002,*)'length of box=',lbox
    write(1002,*)'random seed=',ran
    write(1002,*)'time step=',tstep
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*)'Molecular Dynamics (MD) simulation of 2D Lennard-Jones (LJ) discs'
    write(1002,*)'-------------------------------------------------------------------------'
  
  
    !----------------------------------------------------------
    !MD
    !----------------------------------------------------------
    call random_pos()
    if (f_snap.eq.1) call snapshot(0,'sta')
    write(1002,*)'particles randomly put in box'
    write(1002,*)'-------------------------------------------------------------------------'
  
    call init_pos()
    if (f_snap.eq.1) call snapshot(0,'int')
    write(1002,*)'positions initialized'
    write(1002,*)'-------------------------------------------------------------------------'
  
    call init_vel()
    write(1002,*)'velocities initialized'
    write(1002,*)'-------------------------------------------------------------------------'
  
    call equil()
    if (f_snap.eq.1) call snapshot(0,'ekv')
    write(1002,*)'MD equilibration successful'
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'pp=',pp
    write(1002,*)'eekin=',eekin
    write(1002,*)'ttem=',ttem
    write(1002,*)'eetot=',eetot
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*)'MD simulation initialized'
    write(1002,*)'-------------------------------------------------------------------------'
  
    open(223,file='dpdinit_td.dat')
    write(223,'(4a16)')'temperature','density','ex energy/N','pressure'
    !correlation module
    grc=0.0d0
    do j=1,nseries
      call mdseries(j)
      if (f_snap.eq.1) then
        write(fmt,'(i3.3)')j
        call snapshot(1,trim(fmt))
      endif
      do i=1,10000
        c=Pi*dble(2*i-1)*grinterv**2
        grc(j,i)=grc(j,i)/dble(NC-1)/dble(nsampls)/densc/c
      enddo
      write(223,'(4e16.7)')temp,dens,uu/dble(N),pp
      write(1002,*) j,'-th MD run successful'
      write(1002,*)'uu/N=',uu/dble(N)
      write(1002,*)'pp=',pp
      write(1002,*)'eekin=',eekin
      write(1002,*)'ttem=',ttem
      write(1002,*)'eetot=',eetot
      write(1002,*)'-------------------------------------------------------------------------'
    enddo
    close(223)
  
    open(202,file='dpdinit_corr.dat')
    ggrc=0.0d0
    dgrc=0.0d0
    do i=1,10000
      do j=1,nseries
        ggrc(i)=ggrc(i)+grc(j,i)
      enddo
      ggrc(i)=ggrc(i)/dble(nseries)
      do j=1,nseries
        dgrc(i)=dgrc(i)+(grc(j,i)-ggrc(i))**2
      enddo
      dgrc(i)=dsqrt(dgrc(i)/dble(nseries))
      if ( (dble(i)-0.5d0)*grinterv .lt. lbox ) write(202,'(3e16.7)')(dble(i)-0.5d0)*grinterv,ggrc(i),dgrc(i)
    enddo
    close(202)
  
    write(1002,*)'DONE'
    write(1002,*)'-------------------------------------------------------------------------'
    close(1002)
  
    deallocate( pos )
    deallocate( vel )
    deallocate( acc )
    deallocate( clust )
    deallocate( posc )
    deallocate( velc )
    deallocate( accc )
    deallocate( part )
    deallocate( pot )
    deallocate( tla )
    deallocate( sil )
    deallocate( dis )
    deallocate( dis2 )
    deallocate( rrr )
    deallocate( grc )
    deallocate( ggrc )
    deallocate( dgrc )
  
    open(1002,file='dpdinit_DONE')
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
  subroutine interactions(idum)
    use parameters
    use particles
    use clusters
    use potential
    use microstate
    implicit none
    integer :: idum
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
          !potential module
          pot(i,j)=ljpot(r2,a1,a2,tl)
          tla(i,j)=tl
          sil(i,j,1)=tl*x/r2
          sil(i,j,2)=tl*y/r2
          dis(i,j)=r
          dis2(i,j)=r2
          rrr(i,j,1)=x
          rrr(i,j,2)=y
          pot(j,i)=pot(i,j)
          tla(j,i)=tl
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
    do i=1,N
      acc(i,1)=sum(sil(i,:,1))
      acc(i,2)=sum(sil(i,:,2))
    enddo

    if (idum.eq.1) then
      !clusters module
      accc=0.0d0
      do i=1,N
        accc(clust(i),1)=accc(clust(i),1)+acc(i,1)
        accc(clust(i),2)=accc(clust(i),2)+acc(i,2)
      enddo
      do i=1,NC
        accc(i,:)=accc(i,:)/dble(part(i))
      enddo
    endif
  endsubroutine interactions
  
  !----------------------------------------------------------------------------------------------------
  !Calculates the Lennard-Jones potential between two particles
  !----------------------------------------------------------------------------------------------------
  function ljpot(r2,a1,a2,p)
    implicit none
    !using reduced units
    real*8 :: r2
    real*8 :: a1,a2
    !a1=epsilon = absolute value of the minimum value of the LJ potential (depth of the potential well)
    !a2=sigma = distance at which the potential becomes positive
    real*8 :: x
    real*8 :: ljpot,p !LJ potential, virial pressure
  
    x=a2**2/r2
    x=x**3
    ljpot=4.0d0*a1*x*(x-1.0d0)
    p=48.0d0*a1*x*(x-0.5d0)        ! -dU/dr * r = F * r
    return
  endfunction ljpot
  
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
    call temperature(0)
  
    !rescale the particle velocities to the desired temperature
    rescale=dsqrt(temp/tem)
    vel=vel*rescale
  endsubroutine init_vel
  
  !----------------------------------------------------------------------------------------------------
  !Calculates current temperature using the equipartition theorem
  !----------------------------------------------------------------------------------------------------
  subroutine temperature(idum)
    use parameters
    use particles
    use clusters
    use microstate
    implicit none
    integer :: idum
    integer :: i
  
    ekin=0.0d0
    do i=1,N
      ekin=ekin+(vel(i,1)**2+vel(i,2)**2)/2.0d0
    enddo
    !(x,y) -> 2 degrees of freedom
    !1/2 T for each degree of freedom
    !2N degrees of freedom (N particles with (x,y))
    tem=ekin/dble(N)

    if (idum.eq.1) then
      ekinc=0.0d0
      do i=1,NC
        ekinc=ekinc+part(i)*(velc(i,1)**2+velc(i,2)**2)/2.0d0
      enddo
      temc=ekinc/dble(NC)
    endif
  endsubroutine temperature
  
  !-------------------------------------------------------------------------
  !Velocity Verlet algorithm for calculating velocities
  !Rescaling thermostat
  !Newton's equations of motion
  !-------------------------------------------------------------------------
  subroutine move(idum)
    use parameters
    use particles
    use clusters
    use microstate
    implicit none
    integer :: idum
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
    if (idum.eq.1) then
      do i=1,NC
        !half-update the cluster velocities by tstep/2
        velc(i,1)=velc(i,1)+accc(i,1)*dt2
        velc(i,2)=velc(i,2)+accc(i,2)*dt2
        !update the cluster positions by tstep
        posc(i,1)=posc(i,1)+velc(i,1)*tstep
        posc(i,2)=posc(i,2)+velc(i,2)*tstep
        posc(i,1)=image(0.0d0,posc(i,1),lbox)
        posc(i,2)=image(0.0d0,posc(i,2),lbox)
      enddo
      !update the clusters
      call assign_clusters()
    endif
  
    !update the interactions
    call interactions(idum)
  
    do i=1,N
      !update the particle velocities by tstep/2
      vel(i,1)=vel(i,1)+acc(i,1)*dt2
      vel(i,2)=vel(i,2)+acc(i,2)*dt2
    enddo
    if (idum.eq.1) then
      do i=1,NC
        !update the particle velocities by tstep/2
        velc(i,1)=velc(i,1)+accc(i,1)*dt2
        velc(i,2)=velc(i,2)+accc(i,2)*dt2
      enddo
    endif
    
    !update the temperature
    call temperature(idum)
  
    !rescale the particle velocities to the desired temperature
    rescale=dsqrt(temp/tem)
    vel=vel*rescale
    if (idum.eq.1) then
      rescale=dsqrt(temp/temc)
      velc=velc*rescale
    endif
  endsubroutine move
  
  !----------------------------------------------------------------------------------------------------
  !MD equilibration
  !----------------------------------------------------------------------------------------------------
  subroutine equil()
    use parameters
    use particles
    use microstate
    use macrostate
    implicit none
    integer :: i
    real*8 :: t
  
    call interactions(0)
    call temperature(0)
    t=0.0d0
  
    !macrostate module
    eekin=0.0d0
    ttem=0.0d0
    uu=0.0d0
    eetot=0.0d0
    pp=0.0d0
    nsampls=0
  
    !macrostate = average over microstates
    if (f_sam.eq.1) open(202,file='dpdinit_sam_ekv')
    do i=1,nesteps
      call move(0)
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
  !MD sampling
  !----------------------------------------------------------------------------------------------------
  subroutine mdseries(idum)
    use parameters
    use particles
    use clusters
    use potential
    use microstate
    use macrostate
    use correlation
    implicit none
    integer :: idum
    character(8) :: fmt
    integer :: i
    real*8 :: t
  
    !clusters module
    do i=1,NC
      posc(i,:)=pos(i*N/NC,:)
    enddo
    velc=0.0d0
    call assign_clusters()

    call interactions(1)
    call temperature(1)
    t=0.0d0
  
    !macrostate module
    uu=0.0d0
    pp=0.0d0
    eekin=0.0d0
    ttem=0.0d0
    eetot=0.0d0
    nsampls=0
    
    !macrostate = average over microstates
    if (f_sam.eq.1) then
      write(fmt,'(i3.3)')idum
      open(202,file='dpdinit_sam_'//trim(fmt))
    endif
    do i=1,ntsteps
      call move(1)
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
        call corr(idum)
      endif
      if ( (f_frame.eq.1).and.(idum.eq.1).and.(mod(i,nframe).eq.0).and.(i/nframe.lt.nframes+1) ) then
        write(fmt,'(i3.3,i5.5)')idum,i/nframe
        call movie(1,fmt)
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
  endsubroutine mdseries

  !----------------------------------------------------------------------------------------------------
  !Assigns each particle to a cluster index
  !----------------------------------------------------------------------------------------------------
  subroutine assign_clusters()
    use parameters
    use particles 
    use clusters
    implicit none
    real*8 :: image !functions
    integer :: i,j
    real*8 :: x,y,r2,r2min

    part=0
    do i=1,N
      x=image(pos(i,1),posc(1,1),lbox)
      y=image(pos(i,2),posc(1,2),lbox)
      r2min=x**2+y**2
      clust(i)=1
      do j=2,NC
        x=image(pos(i,1),posc(j,1),lbox)
        y=image(pos(i,2),posc(j,2),lbox)
        r2=x**2+y**2
        if (r2.lt.r2min) then
          r2min=r2
          clust(i)=j
        endif
      enddo
      part(clust(i))=part(clust(i))+1
    enddo
  endsubroutine assign_clusters
  
  !----------------------------------------------------------------------------------------------------
  !Pair correlation function g(r) -- Radial distribution function
  !----------------------------------------------------------------------------------------------------
  subroutine corr(idum)
    use parameters
    use clusters
    use correlation
    implicit none
    integer :: idum
    real*8 :: image !function
    integer :: i,j
    integer :: k
    real*8 :: x,y,r
  
    do i=1,NC-1
      do j=i+1,NC
        x=image(posc(i,1),posc(j,1),lbox)
        y=image(posc(i,2),posc(j,2),lbox)
        r=dsqrt(x**2+y**2)
        k=int(r/grinterv)        !floors double to integer
        if (k.lt.10000) then
          grc(idum,k+1)=grc(idum,k+1)+2.0d0
        endif
      enddo
    enddo
  endsubroutine corr
  
  !-------------------------------------------------------
  !Snapshot
  !-------------------------------------------------------
  subroutine snapshot(idum,fmt)
    use parameters
    use particles
    use clusters
    implicit none
    integer :: idum
    character(3) :: fmt
    integer :: i
  
    open(111,file='dpdinit_snap_'//fmt)
    do i=1,N
      write(111,'(3e16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
    enddo
    close(111)
    if (idum.eq.1) then
      open(111,file='dpdinit_snapc_'//fmt)
      do i=1,NC
        write(111,'(3e16.7)')posc(i,1)/lbox,posc(i,2)/lbox,(4.0d0*Pi/3.0d0*(a2/2.0d0)**3)**(1.0d0/3.0d0)/lbox
      enddo
      close(111)
    endif
  endsubroutine snapshot
  
  !--------------------------------------------------------
  !Movie
  !--------------------------------------------------------
  subroutine movie(idum,fmt)
    use parameters
    use clusters
    use particles
    implicit none
    integer :: idum
    character(8) :: fmt
    integer :: i
  
    open(111,file='dpdinit_frame_'//fmt)
    do i=1,N
      write(111,'(3e16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
    enddo
    close(111)
    if (idum.eq.1) then
      open(111,file='dpdinit_framec_'//fmt)
      do i=1,NC
        write(111,'(3e16.7)')posc(i,1)/lbox,posc(i,2)/lbox,(4.0d0*Pi/3.0d0*(a2/2.0d0)**3)**(1.0d0/3.0d0)/lbox
      enddo
      close(111)
    endif
  endsubroutine movie
  