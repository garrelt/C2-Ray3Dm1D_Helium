
!>
!! \brief This module contains routines for calculating the ionization and temperature evolution of the entire grid (1D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 1D version similar to the 3D version.

module evolve

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of the entire grid (1D).

  ! Version:
  ! 1D version similar to the 3D version.

  use precision, only: dp
  use file_admin, only: logf
  use my_mpi ! supplies all the MPI definitions
  use sizes, only: Ndim, mesh
  use grid, only: x,vol,dr
  use material, only: ndens, xh, xhe, temper
  use photonstatistics, only: state_before, calculate_photon_statistics, &
       photon_loss
  use c2ray_parameters, only: minimum_fractional_change
  use c2ray_parameters, only: minimum_fraction_of_atoms
  use c2ray_parameters, only: minium_fraction_of_photons
  use c2ray_parameters, only: convergence_fraction
  use c2ray_parameters, only: subboxsize, max_subbox
  use c2ray_parameters, only: epsilon
  use abundances, only: abu_he

  implicit none

  !private

  public :: evolve1D !> evolve 1D grid

  !> H column density at the back end of the cell
  real(kind=dp),dimension(mesh(1)) :: coldensh_out
  !> He0 column density at the back end of the cell
  real(kind=dp),dimension(mesh(1)) :: coldenshe0_out
  !> He1 column density at the back end of the cell
  real(kind=dp),dimension(mesh(1)) :: coldenshe1_out
  integer :: smoothlim=0

contains

  ! =======================================================================

  !> Calculates the evolution of the hydrogen ionization state
  !! and temperature for the whole 1D grid\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 21-Aug-2006 (f77/OMP: 13-Jun-2005)\n
  !! \b Version: Simple 1D

  subroutine evolve1D (dt)

    ! Calculates the evolution of the hydrogen ionization state

    ! Author: Garrelt Mellema

    ! Date: 21-Aug-2006 (f77/OMP: 13-Jun-2005)

    ! Version: Simple 1D

    ! History:
    ! 11-Jun-2004 (GM) : grid arrays now passed via common (in grid.h)
    !    and material arrays also (in material.h).
    ! 11-Jun-2004 (GM) : adapted for multiple sources.
    !  3-Jan-2005 (GM) : reintegrated with updated Ifront3D
    ! 20-May-2005 (GM) : split original eveolve0D into two routines
    ! 13-Jun-2005 (HM) : OpenMP version : Hugh Merz
    ! 31-Mar-2008 (GM) : clean up.

    !> The time step
    real(kind=dp),intent(in) :: dt

    !> Will contain the integer position of the cell being treated
    integer,dimension(Ndim) :: pos

    ! Loop variables
    integer :: i

    ! Flag variable (passed back from evolve0D_global)
    integer :: conv_flag

    ! End of declarations

    ! Initial state (for photon statistics)
    call state_before ()

    ! reset photon loss counter
    photon_loss=0.0

    ! Put column densities to zero
    coldensh_out(:)=0.0
    coldenshe0_out(:)=0.0
    coldenshe1_out(:)=0.0

    ! Loop through grid
    do i=1,mesh(1)
       pos(1)=i
       ! if (pos(1).ge.7) pause
       call evolve0D(dt,pos)
       ! write(*,*) '--'

    end do
    !    pause 
    ! Calculate photon statistics
    call calculate_photon_statistics (dt)

    return
  end subroutine evolve1D


  !=======================================================================

  !> Calculates the evolution of the hydrogen ionization state in a
  !! single cell.\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 21-Aug-2006 (20-May-2005, 5-Jan-2005, 02 Jun 2004)\n
  !! \b Version: single source

  subroutine evolve0D(dt,rtpos)

    ! Calculates the evolution of the hydrogen ionization state in a
    ! single cell.

    ! Author: Garrelt Mellema, Martina M. Friedrich (helium)

    ! Date: 21-Aug-2006 (20-May-2005, 5-Jan-2005, 02 Jun 2004)

    ! Version: single sources, fixed temperature

    ! Multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated. 

    use mathconstants, only: pi
    use tped, only: electrondens
    use doric_module, only: doric, coldens, coldens_bndry_HI, &
         coldens_bndry_HeI,coldens_bndry_HeII
    use radiation, only: photoion_rates, photrates
    use material, only: isothermal, ionstates, gamma_uvb  
    use thermalevolution, only: thermal
    use cgsconstants,only : ini_rec_colion_factors

    implicit none

    real(kind=dp),parameter :: max_coldensh=2.0e26 !< column density for stopping chemisty

    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    !> The time step
    real(kind=dp),intent(in) :: dt 
    !> mesh position of cell being done
    integer,dimension(Ndim),intent(in) :: rtpos 

    integer :: nx,nit ! loop counters
    integer,dimension(Ndim) :: pos
    real(kind=dp) :: coldensh_in
    real(kind=dp) :: coldenshe0_in
    real(kind=dp) :: coldenshe1_in
    real(kind=dp) :: coldensh_cell
    real(kind=dp) :: coldenshe0_cell
    real(kind=dp) :: coldenshe1_cell
    real(kind=dp),dimension(0:1) :: hecolum_in,hecolum_out
    real(kind=dp) :: path, hcolum_out
    real(kind=dp) :: de,oldhav,oldhe0av,oldhe1av
    real(kind=dp) :: ndens_p
    real(kind=dp) :: avg_temper
    real(kind=dp) :: yfrac,zfrac
    real(kind=dp) :: y2afrac,y2bfrac
    real(kind=dp) :: dist,vol_ph
    real(kind=dp) :: temper0,temper1,temper2,ionh0old,ionh1old
    real(kind=dp) :: yh0_av_old,yh1_av_old,yhe0_av_old,yhe1_av_old,yhe2_av_old
    real(kind=dp) :: ionhe0old,ionhe1old,ionhe2old,ionheavold
    type(ionstates)  :: ion
    type(photrates) :: phi
    integer :: testpos

    pos(1)=rtpos(1)
    ! Initialize local ionization states to the global ones
    do nx=0,1
       testpos=max(pos(1),1)
       ion%h(nx)=xh(testpos,nx)!ion%h(nx)=xh(pos(1),nx)
       ion%h_old(nx)=xh(pos(1),nx)
       ion%h_av(nx)=xh(testpos,nx)!ion%h_av(nx)=xh(pos(1),nx) ! use calculated xh_av
    enddo

    do nx=0,2
       ion%he(nx)=xhe(testpos,nx)!ion%he(nx)=xhe(pos(1),nx)
       ion%he_old(nx)=xhe(pos(1),nx)
       ion%he_av(nx)=xhe(testpos,nx)!ion%he_av(nx)=xhe(pos(1),nx) 
    enddo

    ! Initialize local temperature and density
    avg_temper=temper(pos(1))
    ! Initialize scalars temper0 and temper1
    temper0=temper(pos(1))    ! original temperature as scalar
    temper1=temper0           ! will contain the new temperature

    ndens_p=ndens(pos(1),1,1)

    ! Find the column density at the entrance point of the cell (short
    ! characteristics)
    if (pos(1).eq.1) then
       coldensh_in = coldens_bndry_HI()
       coldenshe0_in = coldens_bndry_HeI()
       coldenshe1_in = coldens_bndry_HeII()
    else
       coldensh_in = coldensh_out(pos(1)-1)
       coldenshe0_in = coldenshe0_out(pos(1)-1)
       coldenshe1_in = coldenshe1_out(pos(1)-1)
    endif

    path=dr(1)

    ! Find the distance to the source
    dist=x(pos(1))

    ! Find the volume of the shell this cell is part of (dilution factor)
    vol_ph=vol(pos(1))

    ! Iterate to get mean ionization state (column density / optical depth) 
    ! in cell
    nit=0


    if (coldensh_in.le.max_coldensh) then

       do 

          nit=nit+1


          ! Save the value of yh_av found in the previous iteration
          yh0_av_old=ion%h_av(0)
          yh1_av_old=ion%h_av(1)
          yhe0_av_old=ion%he_av(0)
          yhe1_av_old=ion%he_av(1)
          yhe2_av_old=ion%he_av(2)
          temper2=temper1

          !-----------PHOTOION BLOCK-----------------------------------------------
          coldensh_cell  =coldens(path,ion%h_av(0),ndens_p,(1.0_dp-abu_he))
          coldenshe0_cell=coldens(path,ion%he_av(0),ndens_p,abu_he)
          coldenshe1_cell=coldens(path,ion%he_av(1),ndens_p,abu_he)
          hcolum_out=coldensh_in+coldensh_cell
          hecolum_in= (/coldenshe0_in,coldenshe1_in/)
          hecolum_out=(/coldenshe0_in+coldenshe0_cell ,coldenshe1_in+coldenshe1_cell/)       
          ! Calculate (photon-conserving) photo-ionization rate
          phi=photoion_rates(coldensh_in,hcolum_out, &
               hecolum_in(0), hecolum_out(0), &
               hecolum_in(1), hecolum_out(1), &
               vol_ph,1,ion%h_av(1))

          phi%photo_cell_HI=phi%photo_cell_HI/(ion%h_av(0)*ndens_p* &
               (1.0_dp-abu_he))
          phi%photo_cell_HeI=phi%photo_cell_HeI/(ion%he_av(nx)*ndens_p*abu_he)
          phi%photo_cell_HeII=phi%photo_cell_HeII/(ion%he_av(nx)*ndens_p*abu_he)

          !write(logf,*) "phiheat= ",phi%heat, pos
          ! Add the UV background
          phi%photo_cell_HI=phi%photo_cell_HI + gamma_uvb(1)
          phi%photo_cell_HeI=phi%photo_cell_HeI + gamma_uvb(2)
          phi%photo_cell_HeII=phi%photo_cell_HeII + gamma_uvb(3)

          !   phi%h=phi%h/(ndens_p*abu_he*ion%he_av(1)+ndens_p*abu_he*ion%he_av(0)+ndens_p*(1.0_dp-abu_he)*ion%h_av(0))
          !   phi%he(0)=phi%he(0)/(ndens_p*(1.0_dp-abu_he)*ion%h_av(0)+ndens_p*abu_he*ion%he_av(0))
          !phi%he(1)=phi%he(1)/(ndens_p*(1.0_dp-abu_he)*ion%h_av(0)+ndens_p*abu_he*ion%he_av(0)+ndens_p*abu_he*ion%he_av(1))

          !       phi%h=phi%h/ndens_p
          !       phi%he=phi%he/ndens_p
          !--------------------------------------------------------------------

          de=electrondens(ndens_p,ion%h_av,ion%he_av)

          ! initialize the collisional ion and recomb rates for new T
          if (.not.isothermal) call ini_rec_colion_factors(avg_temper) 

          !----------------DORIC BLOCK-----------------------------------------
          coldensh_cell  =coldens(path,ion%h(0),ndens_p,(1.0_dp-abu_he))   ! average or not average?
          coldenshe0_cell=coldens(path,ion%he(0),ndens_p,abu_he)           ! average or not average?
          coldenshe1_cell=coldens(path,ion%he(1),ndens_p,abu_he)           ! average or not average?

          ! Prepare factors needed by doric
          call prepare_doric_factors(coldensh_cell, coldenshe0_cell, &
               coldenshe1_cell, yfrac, zfrac, y2afrac,y2bfrac)

          call doric(dt,de,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac)! calculate de new?  (2)
          !--------------------------------------------------------------------

          de=electrondens(ndens_p,ion%h_av,ion%he_av)

          !----------------DORIC BLOCK------------------------------------------

          coldensh_cell  =coldens(path,ion%h(0),    ndens_p,(1.0_dp-abu_he))  ! average or not average?
          coldenshe0_cell=coldens(path,ion%he(0),ndens_p,         abu_he)  ! average or not average?
          coldenshe1_cell=coldens(path,ion%he_av(1),ndens_p,abu_he)

          ! Prepare factors needed by doric
          call prepare_doric_factors(coldensh_cell, coldenshe0_cell, &
               coldenshe1_cell, yfrac, zfrac, y2afrac,y2bfrac)

          ionh0old=ion%h(0)
          ionh1old=ion%h(1)
          ionhe0old=ion%he(0)
          ionhe1old=ion%he(1)
          ionhe2old=ion%he(2)
          ionheavold=de 
          oldhav=ion%h_av(0)
          oldhe0av=ion%he_av(0)
          oldhe1av=ion%he_av(1)

          call doric(dt,de,ndens_p,ion,phi, yfrac, zfrac,y2afrac,y2bfrac) 

          ion%h(0)=0.5*(ion%h(0)+ionh0old)
          ion%h(1)=0.5*(ion%h(1)+ionh1old)
          ion%he(0)=0.5*(ion%he(0)+ionhe0old)
          ion%he(1)=0.5*(ion%he(1)+ionhe1old)
          ion%he(2)=0.5*(ion%he(2)+ionhe2old)
          ion%h_av(0)=0.5*(ion%h_av(0)+oldhav)
          ion%he_av(0)=0.5*(ion%he_av(0)+oldhe0av)
          ion%he_av(1)=0.5*(ion%he_av(1)+oldhe1av)
          de=electrondens(ndens_p,ion%h_av,ion%he_av)


          temper1=temper0       ! set temper1 to the original temperature

!          if (nit.eq.1.and.abs(ion%h_old(0)-ion%h(0)).gt.0.1 .and. &
!               abs(ion%he_old(0)-ion%he(0)).gt.0.1) then
!             continue
!          else
             if (.not.isothermal) then
                if (pos(1) == mesh(1)/2) write(logf,*) "phiheat= ",phi%heat
                call thermal(dt,temper1,avg_temper,de,ndens_p, &
                     ion,phi)               
             endif
!          endif

          if(  ( (abs(ion%h_av(0)-yh0_av_old)/ion%h_av(0).lt.minimum_fractional_change) &
               .or.                       &
               (ion%h_av(0).lt.minimum_fraction_of_atoms)   &
               )   &
               .and.                             &
               ( (abs(ion%he_av(1)-yhe1_av_old)/ion%he_av(1).lt.minimum_fractional_change) &
               .or.                        &
               (ion%he_av(1).lt.minimum_fraction_of_atoms)   &
               )    &
               .and.                              &
               ( (abs(ion%he_av(2)-yhe2_av_old)/ion%he_av(2).lt.minimum_fractional_change) &
               .or.                        &
               (ion%he_av(2).lt.minimum_fraction_of_atoms)   &
               )    &
               .and.                              &
               (  (abs(ion%he_av(0)-yhe0_av_old)/ion%he_av(0).lt.minimum_fractional_change) &
               .or.                        &
               (ion%he_av(0).lt.minimum_fraction_of_atoms)   &
               )    & 
               .and. &
               abs(temper1-temper2)/temper1.lt.minimum_fractional_change &
               ) then 
             exit
          else

             ! Warn about non-convergence
             if (nit.gt.4000) then
                write(logf,*) 'Convergence failing'
                write(logf,*) 'nit=',nit
                write(logf,*) 'xh: ',ion%h_av(0),yh0_av_old!, ion%h(0)
                write(logf,*) 'xhe0: ',ion%he_av(0), yhe0_av_old!, ion%he(0)
                write(logf,*) 'xhe0: ',ion%he(0) 
                write(logf,*) 'xhe1: ',ion%he_av(1), yhe1_av_old!, ion%he(1)
                write(logf,*) 'xhe2: ',ion%he_av(2), yhe2_av_old!, ion%he(2)
                write(logf,*) 'temper: ',temper1,temper2
                write(logf,*) 'pos(1):', pos(1)
                write(logf,*) 'convcrits',abs(ion%h_av(0)-yh0_av_old)/ion%h_av(0), &
                     abs(ion%he_av(1)-yhe1_av_old)/ion%he_av(1)
                write(logf,*) 'convcrits',abs(ion%he_av(2)-yhe2_av_old)/ion%he_av(2), &
                     abs(ion%he_av(0)-yhe0_av_old)/ion%he_av(0)
                !  pause
                exit
             endif
          endif   !if convergence exit loop      

       enddo ! end of iteration
    else 
       phi%photo_cell_HI=0.0_dp
       phi%photo_out_HI=0.0_dp
       phi%photo_cell_HeI=0.0_dp
       phi%photo_cell_HeII=0.0_dp
       phi%photo_out_HeI=0.0_dp
       phi%photo_out_HeII=0.0_dp
       !GM/130801: not sure what this is useful for
       !phi%int_out(:)=0.0_dp
    endif


    write(48,*) nit

    ! Copy ionic abundances back

    xh(pos(1),:)=ion%h(:)
    xhe(pos(1),:)=ion%he(:)

    ! Copy temperature back
    temper(pos(1))=temper1

    ! Add the (time averaged) column density of this cell
    ! to the total column density (for this source)
    coldensh_out(pos(1))=coldensh_in+ &
         coldens(path,ion%h_av(0),ndens_p,(1.0_dp-abu_he))
    coldenshe0_out(pos(1))=coldenshe0_in+ &
         coldens(path,ion%he_av(0),ndens_p,abu_he)
    coldenshe1_out(pos(1))=coldenshe1_in+ &
         coldens(path,ion%he_av(1),ndens_p,abu_he)
    ! Photon statistics: register number of photons leaving the grid
    if (pos(1).eq.mesh(1)) &
         photon_loss=0.0_dp!photon_loss+ &
    !(phi%int1_out+phi%int2_out+phi%int3_out)*vol(pos(1))/vol_ph
    !(sum(phi%int_out(:)))*vol(pos(1))/vol_ph
  end subroutine evolve0D

  ! ===========================================================================

  subroutine prepare_doric_factors(NHI,NHeI,NHeII,yfrac,zfrac,y2afrac,y2bfrac)

    use cgsphotoconstants, only: sigma_H_heth, sigma_H_heLya,sigma_He_heLya
    use cgsphotoconstants, only: sigma_HeI_at_ion_freq, sigma_He_he2
    use cgsphotoconstants, only: sigma_HeII_at_ion_freq,sigma_He_he2,sigma_H_he2

    real(kind=dp),intent(in) :: NHI ! H0 column density
    real(kind=dp),intent(in) :: NHeI ! He column densities
    real(kind=dp),intent(in) :: NHeII ! He column densities

    real(kind=dp),intent(out) :: yfrac,zfrac
    real(kind=dp),intent(out) :: y2afrac,y2bfrac
 
    real(kind=dp) :: tau_H_heth
    real(kind=dp) :: tau_He_heth
    real(kind=dp) :: tau_H_heLya
    real(kind=dp) :: tau_He_heLya 
    real(kind=dp) :: tau_He2_he2th
    real(kind=dp) :: tau_He_he2th
    real(kind=dp) :: tau_H_he2th

    tau_H_heth  = NHI*sigma_H_heth ! opt depth of HI at HeI ion threshold
    tau_He_heth = NHeI*sigma_HeI_at_ion_freq ! opt depth of HeI at HeI ion threshold
    tau_H_heLya = NHI*sigma_H_heLya ! opt depth of H  at he+Lya (40.817eV)
    tau_He_heLya= NHeI*sigma_He_heLya ! opt depth of He at he+Lya (40.817eV) 
    tau_H_he2th = NHI*sigma_H_he2 ! opt depth of H at HeII ion threshold
    tau_He_he2th = NHeI*sigma_He_he2 ! opt depth of HeI at HeII ion threshold
    tau_He2_he2th = NHeII*sigma_HeII_at_ion_freq ! opt depth of HeII at HeII ion threshold
    
    ! Ratios of these optical depths needed in doric
    yfrac= tau_H_heth /(tau_H_heth +tau_He_heth)
    zfrac= tau_H_heLya/(tau_H_heLya+tau_He_heLya)
    y2afrac=  tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
    y2bfrac=  tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
    
  end subroutine prepare_doric_factors

end module evolve
