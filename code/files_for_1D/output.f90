!>
!! \brief This module contains routines for file output
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2012-10-13 (but older)
!!
!! \b Version: 1D, radial (spherically symmetric) grid.

module output_module
  
  ! This file contains routines having to do with the output
  ! of the C2-Ray program.
  
  ! setup_out : open files
  ! close_down : close files
  ! output : write output

  use precision, only: dp
  use my_mpi, only: rank
  use file_admin, only: stdinput, results_dir, file_input, logf
  use tped

  implicit none
  
  private

  ! To controle what kind of output is produced we use an array of flags.
  ! There can be at most max_output_streams types of output.
  integer,parameter :: max_output_streams=5 !< maximum number of output streams
  integer,dimension(max_output_streams) :: streams !< flag array for output streams

  public :: setup_output, output, close_down

contains
  !----------------------------------------------------------------------------

  !> Initializes global output (photonstatistics, analytical solution)
  !!
  !! Analytical: ASCII table
  !! \li time (s)
  !! \li Numerical position (cm)
  !! \li Analytical position (cm)
  !! \li Relative error between analytical and numerical
  !! \li Uncertainty due to front width
  !! \li Error due to finite cell size.
  !!
  !! PhotonCounts: ASCII table
  !! \li time
  !! \li Number of (ionizations + recombinations) / photons during time step
  !! \li Number of ionizations /(ionizations + recombinations)
  !! \li Number of recombinations /(ionizations + recombinations)
  !! \li Number of (ionizations + recombinations) / photons 
  !! \li Number of (ionizations + recombinations) / photons since t=0

  subroutine setup_output ()
    
    ! Sets up output
    
    ! photon statistics
    use photonstatistics, only: do_photonstatistics, &
         initialize_photonstatistics

    if (rank == 0) then
       ! Open files
       if (do_photonstatistics) then
          open(unit=90,file=trim(adjustl(results_dir))//"PhotonCounts.out", &
               form="formatted",status="unknown",position="append")
          ! Write header with description of what is in the file
          write(90,*) "Columns: Rate of change for H0, H+, He0, He+, He2+, " &
               "total photon losses, rate of ionizations + recombinations".

       !GM/121013: Commented out as it is not used
       !open(unit=91,file=trim(adjustl(results_dir))//'Analytical.out', &
       !     form='formatted',status='unknown')
       !write(91,*) 'Time - Numerical Solution - Analytical solution -', &
       !     'Fractional Error'

       ! GM/121013: Extra diagnostic files
       ! nit: number of iterations?
       ! Outputtimes.out: simulation time corresponding to outputs
       ! Inputnumbers.out: ?
       open(unit=48,file='nit',status='unknown')
       open(unit=50,file='Outputtimes.out',status='unknown')
       open(unit=52,file='Inputnumbers.out',status='unknown')

    endif
    if (do_photonstatistics) call initialize_photonstatistics ()
    
  end subroutine setup_output
  
  !-----------------------------------------------------------------------------

  !> Closes global output files which have been open the entire run

  subroutine close_down ()
    
    ! Closes down
    
    if (rank == 0) then
       close(unit=90)
       !close(unit=91)
       close(unit=48)
       close(unit=50)
       close(unit=52)
    endif

  end subroutine close_down
  
  !----------------------------------------------------------------------------

  !> Produce output for a time frame
  subroutine output (step,time,dt,end_time) !* (step time dt end_time)

    ! Simple output routine.

    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    real(kind=dp),intent(in) :: end_time !< end simulation at this time
    integer, intent(in)      :: step    

    call write_stream1 (step,time,dt,end_time)

    call write_photonstatistics(step,time,dt,end_time) 

    call write_analytical (step,time,dt,end_time) 

  end subroutine output

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  !!
  !! Output format: ASCII table
  !! \li r-position
  !! \li neutral hydrogen fraction
  !! \li ionized hydrogen fraction
  !! \li temperature
  !! \li density
  !!
  !! No time information is available in the output, but the
  !! file names are Ifront_xx.xxx.dat, where xx.xxx is the
  !! fraction of the simulation time passed (initial condition 0.000,
  !! last output 1.000).
 
  subroutine write_stream1 (step,time,dt,end_time) !* (step time dt end_time)

    ! Simple output routine.

    !      Output format:
    !      The output stream (see setup_output)
    !      stream 1: r-position
    !                neutral hydrogen fraction
    !                ionized hydrogen fraction
    !                temperature
    !                density
    !                neutral helium fraction   
    !                single ionized helium fraction
    !                double ionized helium fraction  
    !      Data from subsequent time steps are stacked without
    !      any separators. No time information is available in
    !      the output.
    
    use sizes, only: mesh
    use grid, only: x,dr
    use material
    use photonstatistics
    use radiation, only: teff,rstar,lstar,S_star
    
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    real(kind=dp),intent(in) :: end_time !< end simulation at this time
    integer, intent(in)      :: step    

    integer :: i,j,k,ns
    character(len=40) :: file1

    if (rank == 0) then
       ! Stream 1
!MMF real(time)/end_time
       write(file1,*) step
       file1=trim(adjustl(results_dir))//'Ifront1_'//trim(adjustl(file1))//'.dat'
       open(unit=51,file=file1,form='formatted',status='unknown')
   
       do i=1,mesh(1)
          write(51,'(8(1pe11.4,1x))') x(i),xh(i,0),xh(i,1), temper(i),ndens(i,1,1), &
       xhe(i,0), xhe(i,1),xhe(i,2)
       enddo
    endif
  end subroutine write_stream1

  !----------------------------------------------------------------------------

  subroutine write_photonstatistics (step,time,dt,end_time) !* (step time dt end_time)

    ! Simple output routine.
    
    !      PhotonCounts: time
    !                    Number of (ionizations + recombinations) / photons 
    !                     during time step
    !                    Number of ionizations /(ionizations + recombinations)
    !                    Number of recombinations /(ionizations + recombinations)
    !                    Number of (ionizations + recombinations) / photons 
    !                    Number of (ionizations + recombinations) / photons 
    !                     since t=0

    use sizes, only: mesh
    use grid, only: x,dr
    use material
    use photonstatistics
    use radiation, only: teff,rstar,lstar,S_star
    
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    real(kind=dp),intent(in) :: end_time !< end simulation at this time
    integer, intent(in)      :: step    

    real(kind=dp) :: totalsrc
    logical crossing,recording_photonstats

    if (rank == 0) then
       if (do_photonstatistics .and. time > 0.0) then
          ! Photon Statistics
          total_ion=total_ion!+photon_loss
          grtotal_ion=grtotal_ion+total_ion-totcollisions
          if (time > 0.0) then
             write(90,'(7(1pe13.5))') &
 !                  real(time)/end_time, &
 !                 (total_ion-totcollisions)/(s_star*dt), &
 !                  total_ion, &
                   dh0/dt, &
                   dh1/dt, &
                   dhe0/dt, &
                   dhe1/dt, &
                   dhe2/dt, &
                   photon_loss, &
                   (2.0*dhe0+dhe1+dh0+totrec)/dt
 !                  (dh0+dh1+dhe0+dhe1+dhe2)/dt
 !                  totcollisions
 !                 grtotal_ion/(s_star*time)
		   
          endif
       endif
    endif

  end subroutine write_photonstatistics

  !----------------------------------------------------------------------------

  subroutine write_analytical (step,time,dt,end_time) !* (step time dt end_time)

    ! Simple output routine.
    
    !      Analytical: time (s)
    !                  Numerical position (cm)
    !                  Analytical position (cm)
    !                  Relative error between analytical and numerical
    !                  Uncertainty due to front width
    !                  Error due to finite cell size.
    
    use sizes, only: mesh
    use grid, only: x,dr
    use material
    use photonstatistics
    use radiation, only: teff,rstar,lstar,S_star
    
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    real(kind=dp),intent(in) :: end_time !< end simulation at this time
    integer, intent(in)      :: step    

    integer :: i,j,k,ns
    real(kind=dp) :: ana_front,num_front,num_front1,num_front2

    if (rank == 0) then

       ! Compare to analytical solution (ejr, 04072004)
!       if (time > 0.0) then     
!          call calc_ana_front(ana_front,time)
!          call calc_num_front(num_front1,0.1d0)
!          call calc_num_front(num_front,0.5d0)
!          call calc_num_front(num_front2,0.9d0)
!          if (num_front == 0.0) num_front=1.0
!          if (ana_front == 0.0) ana_front=1.0
!          write(91,'(6(1pe10.3,1x))') time,num_front, &
!               ana_front, &
!               (num_front-ana_front)/ana_front, &
!               (num_front1-num_front2)/num_front, &
!               dr/num_front
!       endif

    endif

  end subroutine write_analytical

  !---------------------------------------------------------------------------

  !> Calculate the location of the ionization front analytically\n
  !! \b Authors: Erik-Jan Rijkhorst, Ilian Iliev, Garrelt Mellema\n
  !! \b Date: 20-Oct-2004 (5-Aug-2004)\n
  !! This routine needs several external functions:
  !!  - LambertW(x): the LambertW function.(real argument and ouput)
  !!  - expint(n,x): the exponentional integral of order n

  subroutine calc_ana_front(ana_front,time)
    !     
    !     Calculate the location of the ionization front.
    !     
    !     Authors: Erik-Jan Rijkhorst, Ilian Iliev, Garrelt Mellema
    !     Date: 20-Oct-2004 (5-Aug-2004)
    
    !     This routine needs several external functions:
    !     LambertW(x): the LambertW function.(real argument and ouput)
    !     expint(n,x): the exponentional integral of order n

    use cosmology
    use material
    use radiation
    use grid
    
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(out) :: ana_front !< front position
    
    real(kind=dp) :: StromgrenRadius
    real(kind=dp) :: typical_dens
    real(kind=dp) :: L,K,a,S1,t_reccore
    real(kind=dp) :: tratio
    !real(kind=dp) :: LambertW
    !real(kind=dp) :: expint     

    ! The analytical solution depends on the input:
    ! test 1: uniform density
    ! test 2: 1/r density
    ! test 3: constant density for r<r_core, 1/r^2 for r>r_core,
    !          and L close to 0 (see below).
    ! test 4: uniform density in expanding universe


    select case (testnum)
    case(1) ! constant density solution
       typical_dens=ndens(1,1,1)
       StromgrenRadius = (3.0*S_star/ &
            (4.0*pi*typical_dens*typical_dens* &
            clumping*bh00))**(1.0/3.0)
       ana_front = StromgrenRadius*(1.0- &
            exp(-typical_dens*clumping*bh00*time))**(1.0/3.0)
    case(2)
       L=S_star/(4.0*pi*dens_core*r_core)
       K=dens_core*r_core*clumping*bh00
       ana_front=L/K*(1.0+LambertW(-exp(-K*K*time/L-1.0)))
    case(3)
       !write(*,*) S_star/(4.0*pi*dens_core*r_core*r_core), &
       !     4./3.*dens_core*r_core*clumping*bh00
       L=S_star/(4.0*pi*dens_core*r_core*r_core) &
            -4./3.*dens_core*r_core*clumping*bh00
       K=dens_core*r_core*r_core*clumping*bh00
       t_reccore=1./(dens_core*clumping*bh00)
       if(abs(L)/(4./3.*dens_core*r_core*clumping*bh00) < 1e-3) then
          ana_front=r_core*sqrt(1.+2.*time/t_reccore)
       else
          write(*,*) 'No analytical solution implemented for', &
               ' these parameters.'
          write(*,*) abs(L)/(4./3.*dens_core*r_core*clumping*bh00)
          ana_front=0.0
          ! The following analytical solution seems to be wrong??
          ! ana_front=K/L*(1.0+ &
          ! LambertW(-(L*r_core/K+1.)/exp(L*L*time/K+L*r_core/K+1.0d0)))
              
       endif
       typical_dens=ndens(1,1,1)
       StromgrenRadius = (3.0*S_star/ &
            (4.0*pi*typical_dens*typical_dens* &
            clumping*bh00))**(1.0/3.0)
       if (time < -t_reccore*log(1.-(r_core/StromgrenRadius)**3)) then
          ana_front = StromgrenRadius*(1.0- &
               exp(-typical_dens*clumping*bh00*time))**(1.0/3.0)
       endif

    case(4)
       StromgrenRadius = (3.0*S_star/ &
            (4.0*pi*dens_core*dens_core* &
            clumping*bh00))**(1.0/3.0) !comoving r_S 
       ! (Note: not necessarily reached!!)
       tratio=t0/(t0+time)
       ana_front = StromgrenRadius*(eta/(1.+zred_t0)**3* & !exp(eta*tratio)* 
            (expint(2,eta*tratio,tratio*eta)/tratio &
            -expint(2,eta,tratio*eta)))**(1./3.) &
            /(1.+zred)
    case default
       write(*,*) 'Unknown test problem'
       ana_front=0.0
    end select

  end subroutine calc_ana_front
      
  !---------------------------------------------------------------------------
      
  !> finds the the location of the numerical ionization front.\n
  !! \b Author: Erik-Jan Rijkhorst\n
  !! \b Date: 5-Aug-2004

  subroutine calc_num_front(num_front,xlimit)
    !     
    !     Calculate the location of the ionization front.
    !     
    !     Author: Erik-Jan Rijkhorst
    !     Date: 5-Aug-2004

    use grid, only: x,dr
    use material, only: xh

    !> ionization fraction that defines front location
    real(kind=dp),intent(in) :: xlimit 
    real(kind=dp),intent(out) :: num_front !< front position
      
    integer :: i, i1, i2

    num_front = x(1)-0.5*dr(1)
    i=1
    i1=i
    i2=i
    do while(xh(i,1).ge.xlimit)
       i1 = i
       i2 = i+1
       i=i+1
    enddo

    if (xh(i1,1) == 0.0 .and. xh(i2,1) == 0.0) then
       num_front = x(1)-0.5*dr(1)
    else
       !remember: problems if too less ionizing photons
       num_front = (xlimit-xh(i1,1)) * (x(i1)-x(i2)) / &
            (xh(i1,1)-xh(i2,1)) + x(i1)
    endif

  end subroutine calc_num_front
  
  !---------------------------------------------------------------------------

  !> Calculates LambertW function for z
  !
  !      /* Lambert W function. 
  !      Was ~/C/LambertW.c written K M Briggs Keith dot Briggs at 
  !      bt dot com 97 May 21.  
  !      Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 
  !      01 Apr 09
  !      Translated to Fortran 77 by Garrelt Mellema, 04 Sep 20.
  
  !      Computes Lambert W function, principal branch.
  !      See LambertW1.c for -1 branch.
  
  !      Returned value W(z) satisfies W(z)*exp(W(z))=z
  !      test data...
  !      W(1)= 0.5671432904097838730
  !      W(2)= 0.8526055020137254914
  !      W(20)=2.2050032780240599705
  !      To solve (a+b*R)*exp(-c*R)-d=0 for R, use
  !      R=-(b*W(-exp(-a*c/b)/b*d*c)+a*c)/b/c
  
  !      Test: 
  !      gcc -DTESTW LambertW.c -o LambertW -lm && LambertW
  !      Library:
  !      gcc -O3 -c LambertW.c */
  
  !      double LambertW(const double z);
  !      const int dbgW=0;
  
  function LambertW(z)
    
    real(kind=dp) :: LambertW
    real(kind=dp),intent(in) :: z !< input value, > -0.367879
    integer i
    real(kind=dp) :: eps,em1
    parameter(eps=4.0e-16)
    parameter(em1=0.3678794411714423215955237701614608)
    real(kind=dp) :: p,e,t,w
    real(kind=dp) :: q,r,q2,q3
    if (z < -em1) then
       write(*,*) "LambertW: bad argument ",z," exiting"
       return
    endif
    if (z == 0.0) then
       LambertW=0.0
    elseif (z < -em1+1e-4) then !series near -em1 in sqrt(q)
       q=z+em1
       r=sqrt(q)
       q2=q*q
       q3=q2*q
       LambertW=-1.0 &
            +2.331643981597124203363536062168*r &
            -1.812187885639363490240191647568*q &
            +1.936631114492359755363277457668*r*q &
            -2.353551201881614516821543561516*q2 &
            +3.066858901050631912893148922704*r*q2 &
            -4.175335600258177138854984177460*q3 &
            +5.858023729874774148815053846119*r*q3 &
            -8.401032217523977370984161688514*q3*q !! error approx 1e-16
       
       ! initial approx for iteration... */
    else
       if (z < 1.0) then    !! /* series near 0 */
          p=sqrt(2.0* &
               (2.7182818284590452353602874713526625*z+1.0))
          w=-1.0+p*(1.0+p*(-0.333333333333333333333+ &
               p*0.152777777777777777777777))
       else 
          w=log(z)          !! /* asymptotic */
       endif
       if (z > 3.0) w=w-log(w) !! /* useful? */
       do i=0,19             !! /* Halley iteration */
          e=exp(w)
          t=w*e-z
          p=w+1.0
          t=t/(e*p-0.5*(p+1.0)*t/p)
          w=w-t
          if (abs(t) < eps*(1.0+abs(w))) then
             LambertW=w    !! /* rel-abs error */
             goto 100
          endif
          !     /* should never get here */
       enddo
       write(*,*) "LambertW: No convergence at z= ",z," exiting."
       write(*,*) abs(t),eps*(1.0+abs(w))
    endif
      
100 continue
  end function LambertW

  !---------------------------------------------------------------------------

  !> calculates exponential integral of order n for x

  FUNCTION expint(n,x,etatratio)
    

    REAL(kind=dp) :: expint

    INTEGER,intent(in) :: n !< order
    REAL(kind=dp),intent(in) :: x !< argument
    REAL(kind=dp),intent(in) ::etatratio !< asymptotic value

    integer,parameter :: MAXIT=100
    REAL(kind=dp),parameter :: EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649

    INTEGER :: i,ii,nm1
    REAL(kind=dp) :: a,b,c,d,del,fact,h,psi

    nm1=n-1
    ! print*,'expint called',n,x,etatratio
    ! print*,'args',n,x
    if(n<0 .or. x<0. .or. (x == 0. .and. (n==0 .or. n==1)))then
       pause 'bad arguments in expint'
    else if(n == 0)then
       expint=exp(-x)/x
    else if(x == 0.)then
       expint=1./nm1
    else if(x > 1.)then
       b=x+n
       c=1./FPMIN
       d=1./b
       h=d
       do i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.) < EPS)then
             expint=h*exp(-x+etatratio)
             ! print*,'expint done',expint,x,etatratio
             return
          endif
       enddo
       pause 'continued fraction failed in expint'
    else
       if(nm1 /= 0)then
          expint=1./nm1
       else
          expint=-log(x)-EULER
       endif
       fact=1.
       do i=1,MAXIT
          fact=-fact*x/i
          if(i /= nm1)then
             del=-fact/(i-nm1)
          else
             psi=-EULER
             do ii=1,nm1
                psi=psi+1./ii
             enddo
             del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          ! if(abs(del) < abs(expint)*EPS) return
          if (abs(del) < abs(expint)*EPS) then
             expint=expint*exp(etatratio)
             return
          end if
       enddo
       pause 'series failed in expint'
    endif
    return
  END FUNCTION expint
  !  (C) Copr. 1986-92 Numerical Recipes Software 

end module output_module
