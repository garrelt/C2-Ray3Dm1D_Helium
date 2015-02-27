module cooling

  ! This module file contains subroutines for radiative cooling
  ! - coolin     - calculate cooling rate
  ! - cooling_initialization - setup cooling table

  ! Version: H-only cooling

  use precision, only: dp
  use abundances, only: abu_he
  use my_mpi
  use file_admin, only: logf

  implicit none

  integer,parameter,private :: temppoints=801 !< number of points in cooling tables
  real(kind=dp),dimension(temppoints),private :: h0_cool !< H0 cooling table
  real(kind=dp),dimension(temppoints),private :: h1_cool !< H+ cooling table
  real(kind=dp),dimension(temppoints),private :: he0_cool !< He0 cooling table 
  real(kind=dp),dimension(temppoints),private :: he1_cool !< He1 cooling table 
  real(kind=dp),dimension(temppoints),private :: he2_cool !< He2 cooling table 
  real(kind=dp),private :: mintemp !< lowest log10(temperature) in table
  real(kind=dp),private :: dtemp !< steps in log10(temperature) in table
  real(kind=dp),private :: dumi

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the cooling rate
  function coolin(nucldens,eldens,xHI,xHII,xHeI,xHeII,xHeIII,temp0)
    
    implicit none

    real(kind=dp) :: coolin
    real(kind=dp),intent(in) :: nucldens !< number density
    real(kind=dp),intent(in) :: eldens !< electron density
    real(kind=dp),intent(in) :: xHI
    real(kind=dp),intent(in) :: xHII
    real(kind=dp),intent(in) :: xHeI
    real(kind=dp),intent(in) :: xHeII
    real(kind=dp),intent(in) :: xHeIII
    real(kind=dp),intent(in) :: temp0 !< temperature

    real(kind=dp) :: tpos, dtpos
    integer :: itpos,itpos1

    tpos=(log10(temp0)-mintemp)/dtemp+1.0d0
    itpos=min(temppoints-1,max(1,int(tpos)))
    dtpos=tpos-real(itpos)
    itpos1=min(temppoints,itpos+1)

    ! Cooling curve
    coolin=nucldens*eldens*(( &
         xHI*(h0_cool(itpos)+ &
         (h0_cool(itpos1)-h0_cool(itpos))*dtpos)+ &
         xHII*(h1_cool(itpos)+ &
         (h1_cool(itpos1)-h1_cool(itpos))*dtpos))*(1.0_dp-abu_he)+ &
         (xHeI*(he0_cool(itpos)+ &
         (he0_cool(itpos1)-he0_cool(itpos))*dtpos)+&
         xHeII*(he1_cool(itpos)+ &
         (he1_cool(itpos1)-he1_cool(itpos))*dtpos)+&
         xHeIII*(he2_cool(itpos)+ &
         (he2_cool(itpos1)-he2_cool(itpos))*dtpos))*abu_he)

  end function coolin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read in and set up the cooling table(s)
  subroutine cooling_initialization()

    implicit none

    real(kind=dp),dimension(temppoints) :: temp
    integer :: itemp
    integer :: element,ion,nchck
    character(len=512) :: coolingfile

    if (rank .eq. 0) then
      write(logf,*) "Beginning of cooling initialization" 
    endif

    coolingfile = "tables/H0-cool.tab"
    if (rank .eq. 0) then
      write(logf,*) "Reading cooling file from ",trim(adjustl(coolingfile))
    endif

    ! Open cooling table (H0)
    open(unit=22,file=coolingfile,status='old')
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),h0_cool(itemp)
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    close(22)
    
    mintemp=temp(1)
    dtemp=temp(2)-temp(1)
    
    coolingfile = "tables/H1-cool-B.tab"
    if (rank .eq. 0) then
      write(logf,*) "Reading cooling file from ",trim(adjustl(coolingfile))
    endif

    ! Open cooling table (H1)
    open(unit=22,file=coolingfile,status='old')
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),h1_cool(itemp)
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    close(22)

    ! Open cooling table (He0)
    !open(unit=22,file='He0-cool.tab',status='old')
    ! Since it is not clear what is in the old table, we use now:
    ! The data from the new table is compiled from:
    ! collisional ionization from Hui&Gnedin 1997
    coolingfile = "tables/He0-cool_new.tab"
    if (rank .eq. 0) then
      write(logf,*) "Reading cooling file from ",trim(adjustl(coolingfile))
    endif

    open(unit=22,file=coolingfile,status='old')    
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),he0_cool(itemp) 
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    close(22)


    ! Open cooling table (He1)
    !open(unit=22,file='He1-cool.tab',status='old')   
    ! Since it is not clear what is in the old table, we use now:     
    ! The data from the new table is compiled from:
    ! free-free and recombination B from Hummer&Storey 1998
    ! collisional excitaion and ionization + dielectronic recomb from Hui&Gnedin 1997
    coolingfile = "tables/He1-cool_new_nocollion.tab"
    if (rank .eq. 0) then
      write(logf,*) "Reading cooling file from ",trim(adjustl(coolingfile))
    endif

    open(unit=22,file=coolingfile,status='old')  
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),he1_cool(itemp)
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    close(22)

    coolingfile = "tables/He2-cool.tab"
    if (rank .eq. 0) then
      write(logf,*) "Reading cooling file from ",trim(adjustl(coolingfile))
    endif

    ! Open cooling table (He2)
    open(unit=22,file=coolingfile,status='old')  
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),he2_cool(itemp)
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    
    close(22)
    
    ! Convert cooling to from log to linear 
    do itemp=1,temppoints
       h0_cool(itemp)=10.0d0**h0_cool(itemp)
       h1_cool(itemp)=10.0d0**h1_cool(itemp)
       he0_cool(itemp)=10.0d0**he0_cool(itemp)
       he1_cool(itemp)=10.0d0**he1_cool(itemp)
       he2_cool(itemp)=10.0d0**he2_cool(itemp)
    enddo

    if (rank .eq. 0) then
      write(logf,*) "End of cooling initialization" 
      write(logf,*)      
      flush(logf) 
    endif

  end subroutine cooling_initialization
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cooling
