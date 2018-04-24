!>
!! \brief This module contains routines for calculating the ionization and 
!! temperature evolution of a single point on the grid (3D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2013-09-05
!!
!! \b Version: 3D, MPI & OpenMP

module evolve_point

  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  ! - evolve0D : calculate the rates for one grid point, possibly
  !              apply them
  ! - evolve0D_global: apply the rates to a single point
  ! - do_chemistry: apply the rates

  ! Needs:
  ! doric : ionization calculation for one point + photo-ionization rates
  ! tped : temperature,pressure,electron density calculation

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,timefile,iterdump, results_dir, dump_dir
  use mathconstants, only: pi
  use cgsconstants, only: ini_rec_colion_factors
  use abundances, only: abu_he
  use c2ray_parameters, only: minimum_fractional_change
  use c2ray_parameters, only: minimum_fraction_of_atoms
  use c2ray_parameters, only: epsilon
  use c2ray_parameters, only: type_of_clumping, use_LLS,type_of_LLS
  use c2ray_parameters, only: add_photon_losses
  use sizes, only: Ndim, mesh
  use grid, only: vol,dr
  use material, only: ndens
  use material, only: xh,xhe
  use material, only: xh_hot
  use material, only: special
  use material, only: get_temperature_point, set_temperature_point
  use material, only: set_final_temperature_point, isothermal
  use material, only: ionstates
  use material, only: clumping_point
  use material, only: coldensh_LLS, LLS_point
  use sourceprops, only: srcpos
  use radiation_photoionrates, only: photrates, photoion_rates
  use radiation_photoionrates, only: prepare_column_densities
  use radiation_photoionrates, only: set_photrates_to_zero
  use thermalevolution, only: thermal
  use photonstatistics, only: photon_loss, total_LLS_loss
  use tped, only: electrondens
  use doric_module, only: doric, prepare_doric_factors, coldens
  use evolve_data, only: bb_phih_grid, bb_phihe_grid, bb_phiheat_grid
#ifdef QUASARS
  use evolve_data, only: qpl_phih_grid, qpl_phihe_grid, qpl_phiheat_grid
#endif
#ifdef PL
  use evolve_data, only: pl_phih_grid, pl_phihe_grid, pl_phiheat_grid
#endif
  use evolve_data, only: xh_av, xhe_av, xh_intermed, xhe_intermed
  use evolve_data, only: xh_hot_av, xh_hot_intermed
  use evolve_data, only: coldensh_out, coldenshe_out
  use evolve_data, only: photon_loss_src_thread
  use evolve_data, only: last_l,last_r
  use evolve_data, only: tn
  use evolve_data, only: epsilon_dx, limit_partial_cells
  use column_density, only: cinterp

  implicit none

  save

  private

  ! Flag to know whether the do_chemistry routine was called with local option
  logical,public ::local_chemistry=.false.

  public evolve0d, evolve0d_global

contains

  !=======================================================================

  !> Calculates the photo-ionization rate for one cell due to one source
  !! and adds this contribution to the collective rate.
  subroutine evolve0D(dt,rtpos,ns,niter,phase)
    
    ! Note for multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated.
    ! For the first pass (niter = 1) it makes sense to DO update the
    ! ionization fractions since this will increase convergence speed
    ! in the case of isolated sources.

    ! column density for stopping chemistry !***how should this criterion be for including he and more than one freq bands?
    ! for the moment, leave it as it is, it's probably ok. 
    real(kind=dp),parameter :: max_coldensh=2e29!2.0e22_dp!2e19_dp 
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    ! subroutine arguments
    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: rtpos ! cell position (for RT)
    integer,intent(in)      :: ns ! source number 
    integer,intent(in)      :: niter ! global iteration number
    character(len=1),intent(in) :: phase
    
    integer :: nx,nd,idim ! loop counters
    integer,dimension(Ndim) :: pos
    real(kind=dp) :: dist2,path,vol_ph
    real(kind=dp) :: xs,ys,zs
    real(kind=dp) :: coldensh_in
    real(kind=dp),dimension(0:1) :: coldenshe_in,coldenshe_out_temp
    real(kind=dp) :: ndens_p
    
    type(photrates) :: bb_phi
#ifdef PL
    type(photrates) :: pl_phi
#endif
#ifdef QUASARS
    type(photrates) :: qpl_phi
#endif
    type(ionstates) :: ion
    logical :: simple_flag

    ! Map pos to mesh pos, assuming a periodic mesh
    do idim=1,Ndim
       pos(idim)=modulo(rtpos(idim)-1,mesh(idim))+1
    enddo

    ! If coldensh_out is zero, we have not done this point
    ! yet, so do it. Otherwise do nothing. (grid is set to 0 for every source)
    if (coldensh_out(pos(1),pos(2),pos(3)) == 0.0) then

       ! Initialize local ionization states to the global ones
       do nx=0,1
          ion%h_av(nx)=max(xh_av(pos(1),pos(2),pos(3),nx),epsilon)
       enddo

       do nx=0,2
          ion%he_av(nx)=max(xhe_av(pos(1),pos(2),pos(3),nx),epsilon)
          ion%he(nx)=max(xhe_intermed(pos(1),pos(2),pos(3),nx),epsilon)
          ion%he_old(nx)=max(xhe(pos(1),pos(2),pos(3),nx),epsilon)         
       enddo

       ! Modify in case of partially ionized cells
       if (xh_hot(pos(1),pos(2),pos(3)) <= limit_partial_cells) then
          ion%h_av(0)=ion%h_av(0)*(1.0-xh_hot_av(pos(1),pos(2),pos(3)))
          ion%h_av(1)=ion%h_av(1)*(1.0-xh_hot_av(pos(1),pos(2),pos(3))) + &
               xh_hot_av(pos(1),pos(2),pos(3))
          ion%he_av(0)=ion%he_av(0)*(1.0-xh_hot_av(pos(1),pos(2),pos(3)))
          ion%he_av(1)=ion%he_av(1)*(1.0-xh_hot_av(pos(1),pos(2),pos(3))) + &
               xh_hot_av(pos(1),pos(2),pos(3))*(1.0d0-ion%he_av(2))
       endif
      
       ! Initialize local density
       ndens_p=ndens(pos(1),pos(2),pos(3))

       ! Find the column density at the entrance point of the cell (short
       ! characteristics)
       
       if ( all( rtpos(:) == srcpos(:,ns) ) ) then
          ! Do not call cinterp for the source point.
          ! Set coldensh and path by hand
          coldensh_in=0.0
          coldenshe_in(:)=0.0
          path=0.5*dr(1)

          ! Find the distance to the source (average?)
          !dist2=0.5*dr(1) !NOT NEEDED         ! this makes vol=dx*dy*dz
          !vol_ph=4.0/3.0*pi*dist2**3
          vol_ph=dr(1)*dr(2)*dr(3)
          
       else
          
          ! For all other points call cinterp to find the column density
          call cinterp(rtpos,srcpos(:,ns),coldensh_in,coldenshe_in(0), &
               coldenshe_in(1),path)

          path=path*dr(1)
          
          ! Find the distance to the source
          xs=dr(1)*real(rtpos(1)-srcpos(1,ns))
          ys=dr(2)*real(rtpos(2)-srcpos(2,ns))
          zs=dr(3)*real(rtpos(3)-srcpos(3,ns))
          dist2=xs*xs+ys*ys+zs*zs
          
          ! Find the volume of the shell this cell is part of 
          ! (dilution factor).
          vol_ph=4.0*pi*dist2*path

          ! Add LLS opacity
          ! GM/110224: previously we added this to coldensh_out,
          ! but this gives funny results for the photo-ionization
          ! rate which is based on the difference between the
          ! in and out column density. To mimick a LLS fog, it
          ! may be better to add it here.
          ! Initialize local LLS (if type of LLS is appropriate)
          if (use_LLS) then
             if (type_of_LLS == 2) call LLS_point (pos(1),pos(2),pos(3))
             coldensh_in = coldensh_in + coldensh_LLS * path/dr(1)
          endif          
       endif

       ! Add the (time averaged) column density of this cell
       ! to the total column density (for this source)
       coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
            coldens(path,ion%h_av(0),ndens_p,(1.0_dp-abu_he))

       coldenshe_out(pos(1),pos(2),pos(3),0)=coldenshe_in(0) + &
            coldens(path,ion%he_av(0),ndens_p,abu_he)

       coldenshe_out(pos(1),pos(2),pos(3),1)=coldenshe_in(1) + &
            coldens(path,ion%he_av(1),ndens_p,abu_he)
       
       ! Calculate (photon-conserving) photo-ionization rate from the
       ! column densities.

       if (phase == "H") then
          ! For the H or hot phase we only need the H ionizing photoionization
          ! rate from "soft" stellar or "B" sources.
          ! This is used elsewehere to calculate the extent of the ionized
          ! regions and helps to establish which cells are only partly inside an
          ! ionized (hot) region. These cells need to be treated differently.

          ! Limit the calculation to a certain maximum column density (hydrogen)
          if (coldensh_in < max_coldensh) then 

             coldenshe_out_temp=coldenshe_out(pos(1),pos(2),pos(3),:)
             ! photoion_rates the structure of rates (photo and heating)
             simple_flag=.true. ! only get the H photoion rate
             ! Calculate the optical depth quantities needed
             call prepare_column_densities(coldensh_in, &
                  coldensh_out(pos(1),pos(2),pos(3)), &
                  coldenshe_in(0),coldenshe_out_temp(0), &
                  coldenshe_in(1),coldenshe_out_temp(1),simple_flag)
             ! Find the photoionization rates
             bb_phi=photoion_rates(vol_ph,ns,ion%h_av(1),"B",simple_flag)

             if (pos(1) == 49 .and. pos(2) == 50 .and. pos(3) == 50) then
                write(*,*) phase, bb_phi%photo_cell_HI, coldensh_in, &
                     coldensh_out(pos(1),pos(2),pos(3))
             endif
             ! Divide the photo-ionization rates by the appropriate neutral
             ! density (part of the photon-conserving rate prescription)
             call divide_rate_by_density (bb_phi,ion,ndens_p)

             ! Calculate the losses due to LLSs.
             ! GM/110224: Add the factor vol/vol_ph to phi, just as we do for
             ! the photon losses.
             ! GM/110302: Use phi%h_in (not out) here since this is where we
             ! draw the photons from, above.
             if (use_LLS) call total_LLS_loss(bb_phi%photo_in_HI*vol/vol_ph, &
                  coldensh_LLS * path/dr(1))

          else
             
             ! If the H0 column density is above the maximum, set rates to zero
             call set_photrates_to_zero (bb_phi)

          endif

          ! Add value to grid of values
          bb_phih_grid(pos(1),pos(2),pos(3)) = &
            bb_phih_grid(pos(1),pos(2),pos(3))+bb_phi%photo_cell_HI

          ! Photon statistics: register number of photons leaving the grid
          ! Note: This is only the H0 photo-ionization rate
          if ( (any(rtpos(:) == last_l(:))) .or. &
               (any(rtpos(:) == last_r(:))) ) then
             photon_loss_src_thread(tn)=photon_loss_src_thread(tn) + &
                  bb_phi%photo_out*vol/vol_ph
             !photon_loss_src(1,tn)=photon_loss_src(1,tn) + phi%h_out*vol/vol_ph
          endif
          
       else
          ! If not for the "hot" phase we need to find all the photoionization
          ! and heating rates, for all types of sources.

          coldenshe_out_temp=coldenshe_out(pos(1),pos(2),pos(3),:)
          simple_flag=.false. ! get all the rates

          call prepare_column_densities(coldensh_in, &
                  coldensh_out(pos(1),pos(2),pos(3)), &
                  coldenshe_in(0),coldenshe_out_temp(0), &
                  coldenshe_in(1),coldenshe_out_temp(1),simple_flag)

          if (coldensh_in < max_coldensh) then 

             ! photoion_rates the structure of rates (photo and heating)
             bb_phi=photoion_rates(vol_ph,ns,ion%h_av(1),"B",simple_flag)

             if (pos(1) == 49 .and. pos(2) == 50 .and. pos(3) == 50) then
                write(*,*) phase, bb_phi%photo_cell_HI, coldensh_in, &
                     coldensh_out(pos(1),pos(2),pos(3))
             endif
          else
             
             ! If the H0 column density is above the maximum, set rates to zero
             call set_photrates_to_zero (bb_phi)

          endif
       

#ifdef QUASARS
          ! photoion_rates the structure of rates (photo and heating)
          qpl_phi=photoion_rates(vol_ph,ns,ion%h_av(1),"Q",simple_flag)
#endif
#ifdef PL
          ! photoion_rates the structure of rates (photo and heating)
          pl_phi=photoion_rates(vol_ph,ns,ion%h_av(1),"P",simple_flag)
#endif
          
          ! Divide the photo-ionization rates by the appropriate neutral density
          ! (part of the photon-conserving rate prescription)
          call divide_rate_by_density (bb_phi,ion,ndens_p)

#ifdef QUASARS
          call divide_rate_by_density (qpl_phi,ion,ndens_p)
#endif
#ifdef PL
          call divide_rate_by_density (pl_phi,ion,ndens_p)
#endif
          ! Calculate the losses due to LLSs.
          ! GM/110224: Add the factor vol/vol_ph to phi, just as we do for
          ! the photon losses.
          ! GM/110302: Use phi%h_in (not out) here since this is where we draw
          ! the photons from, above.
          if (use_LLS) call total_LLS_loss(bb_phi%photo_in_HI*vol/vol_ph, &
               coldensh_LLS * path/dr(1))          

          ! Add values to grid of values
          bb_phih_grid(pos(1),pos(2),pos(3))= &
               bb_phih_grid(pos(1),pos(2),pos(3))+bb_phi%photo_cell_HI
          bb_phihe_grid(pos(1),pos(2),pos(3),0)=&
             bb_phihe_grid(pos(1),pos(2),pos(3),0)+bb_phi%photo_cell_HeI
          bb_phihe_grid(pos(1),pos(2),pos(3),1)=&
               bb_phihe_grid(pos(1),pos(2),pos(3),1)+bb_phi%photo_cell_HeII
          bb_phiheat_grid(pos(1),pos(2),pos(3))=&
               bb_phiheat_grid(pos(1),pos(2),pos(3))+bb_phi%heat
#ifdef QUASARS
          qpl_phih_grid(pos(1),pos(2),pos(3))= &
               qpl_phih_grid(pos(1),pos(2),pos(3))+qpl_phi%photo_cell_HI
          qpl_phihe_grid(pos(1),pos(2),pos(3),0)=&
             qpl_phihe_grid(pos(1),pos(2),pos(3),0)+qpl_phi%photo_cell_HeI
          qpl_phihe_grid(pos(1),pos(2),pos(3),1)=&
               qpl_phihe_grid(pos(1),pos(2),pos(3),1)+qpl_phi%photo_cell_HeII
          qpl_phiheat_grid(pos(1),pos(2),pos(3))=&
               qpl_phiheat_grid(pos(1),pos(2),pos(3))+qpl_phi%heat
#endif
#ifdef PL
          pl_phih_grid(pos(1),pos(2),pos(3))= &
               pl_phih_grid(pos(1),pos(2),pos(3))+pl_phi%photo_cell_HI
          pl_phihe_grid(pos(1),pos(2),pos(3),0)=&
             pl_phihe_grid(pos(1),pos(2),pos(3),0)+pl_phi%photo_cell_HeI
          pl_phihe_grid(pos(1),pos(2),pos(3),1)=&
               pl_phihe_grid(pos(1),pos(2),pos(3),1)+pl_phi%photo_cell_HeII
          pl_phiheat_grid(pos(1),pos(2),pos(3))=&
               pl_phiheat_grid(pos(1),pos(2),pos(3))+pl_phi%heat
#endif          
          ! Photon statistics: register number of photons leaving the grid
          ! Note: This is only the H0 photo-ionization rate
          if ( (any(rtpos(:) == last_l(:))) .or. &
               (any(rtpos(:) == last_r(:))) ) then
             photon_loss_src_thread(tn)=photon_loss_src_thread(tn) + &
                  vol/vol_ph* &
                  (bb_phi%photo_out &
#ifdef QUASARS
                  + qpl_phi%photo_out &
#endif
#ifdef PL
                  + pl_phi%photo_out &
#endif
                  )
          endif
       endif

       if (pos(1) == 49 .and. pos(2) == 50 .and. pos(3) == 50) then
          write(*,*) phase, bb_phi%photo_cell_HI, ion%h_av(0)
       endif

    endif ! end of coldens test
    
  end subroutine evolve0D

  ! =======================================================================

  subroutine divide_rate_by_density (this_phi,ion,numdens)

    type(photrates) :: this_phi
    type(ionstates) :: ion
    real(kind=dp) :: numdens
    
    this_phi%photo_cell_HI=this_phi%photo_cell_HI/ &
         (ion%h_av(0)*numdens*(1.0_dp-abu_he))
    this_phi%photo_cell_HeI=this_phi%photo_cell_HeI/ &
         (ion%he_av(0)*numdens*abu_he)
    this_phi%photo_cell_HeII=this_phi%photo_cell_HeII/&
         (ion%he_av(1)*numdens*abu_he)
    
  end subroutine divide_rate_by_density

  ! =======================================================================

  !> Calculates the evolution of the ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolve0D_global(dt,pos,conv_flag,phase)

    ! Calculates the evolution of the hydrogen + helium ionization state for
    ! one cell (pos) and multiple sources.

    ! Author: Garrelt Mellema

    ! Date: 11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence.

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: conv_flag ! convergence counter
    character(len=1),intent(in) :: phase

    integer :: conv_flag_check=0

    ! Initialize convergence check variable to zero
    conv_flag_check=0
    
    if (phase == "H") then
       call evolve0D_global_partial(dt,pos,conv_flag,phase)
    else
       ! The other phase is always processed with the evolve0D_global_full
       ! routine.
       call evolve0D_global_full(dt,pos,conv_flag,phase)
    endif

  end subroutine evolve0D_global
  
  ! =======================================================================

  !> Calculates the evolution of ionized partial volume of a cell containing
  !! an ionization front of a stellar (soft spectrum / BB) sources
  subroutine evolve0D_global_partial(dt,pos,conv_flag,phase)

    ! Author: Garrelt Mellema

    ! Date: 24-Nov-2017 (11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionized
    ! partial volume is calculated, as well as its time averaged value.
    ! We assume that within the ionized partial volume H is fully ionized
    ! and He is singly ionized. If needed we assume that T=10^4 K in this
    ! region.

    use cgsconstants, only: bh00, colh0, temph0

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: conv_flag ! convergence counter

    integer :: conv_flag_check=0
    real(kind=dp) :: allrates, deltht, ee, avg_factor, ndens_p, eqxfh1
    real(kind=dp) :: recombinations, ionizations, collisions=0.0
    real(kind=dp), pointer :: xh_hot_av_point
    real(kind=dp) :: xh_hot_av_old
    character(len=1),intent(in) :: phase
    real(kind=dp) :: avg_temper, temper0, temper1,temper2,temper_inter
    
    ! Initialize conv_flag_check to zero
    conv_flag_check=0
    
    ! Store previously calculated value of xh_hot_av
    xh_hot_av_old = xh_hot_av(pos(1),pos(2),pos(3))
    xh_hot_av_point => xh_hot_av(pos(1),pos(2),pos(3))
    ndens_p=ndens(pos(1),pos(2),pos(3))
    ! Initialize local temperature
    call get_temperature_point (pos(1),pos(2),pos(3), &
         temper_inter,avg_temper,temper1)

    ! Grow the ionized region (partial volume) as if there were no
    ! recombinations
    !deltht=bb_phih_grid(pos(1),pos(2),pos(3))*dt
    ! We do use recombinations
    recombinations=xh_hot_av_point*ndens_p*bh00
    !collisions=colh0*sqrt(avg_temper)*exp(-temph0/avg_temper)
    ionizations=bb_phih_grid(pos(1),pos(2),pos(3))+collisions
    allrates=(ionizations+recombinations)
    if (allrates <= 0.0d0) then
       !write(*,*) pos(1),pos(2),pos(3), bb_phih_grid(pos(1),pos(2),pos(3)), &
       !     collisions,recombinations,avg_temper
       eqxfh1=0.0
    else
       eqxfh1=ionizations/allrates
    endif
    deltht=allrates*dt
    ee=exp(-deltht)
    xh_hot_intermed(pos(1),pos(2),pos(3)) = eqxfh1 + &
         (xh_hot(pos(1),pos(2),pos(3))-eqxfh1) * ee

    !if (pos(1) == 49 .and. pos(2) == 50 .and. pos(3) == 50) then
    !   write(*,*) ionizations, dt
    !   write(*,*) xh_hot_intermed(pos(1),pos(2),pos(3)),xh_hot(pos(1),pos(2),pos(3))
    !   write(*,*) eqxfh1, ee, recombinations
    !endif
    ! Determine average ionized region size (partial volume) over the time step
    ! Mind fp fluctuations. (1.0-ee)/deltht should go to 1.0 for
    ! small deltht, but finite precision leads to values slightly
    ! above 1.0 and for very small values even to 0.0.
    if (deltht < 1.0e-8) then
       avg_factor=1.0
    else
       avg_factor=(1.0-ee)/deltht
    endif
    xh_hot_av(pos(1),pos(2),pos(3)) = eqxfh1 + &
         (xh_hot(pos(1),pos(2),pos(3))-eqxfh1) * avg_factor
    !xh_hot_av(pos(1),pos(2),pos(3)) = 1.0 + &
    !     (xh_hot(pos(1),pos(2),pos(3))-1.0) * avg_factor

    ! Check hydrogen ionized fraction for convergence
    if ( (abs((xh_hot_av_point-xh_hot_av_old)) > minimum_fractional_change &
         .and. &
         abs((xh_hot_av_point - xh_hot_av_old)/xh_hot_av_point) > minimum_fractional_change &
         .and. &
         (1.0-xh_hot_av_point > minimum_fraction_of_atoms)  )) conv_flag_check=1

    ! Add local check to convergence flag. If any of the above checks failed
    ! the conv_flag will be increased by 1
    conv_flag=conv_flag+conv_flag_check
    
  end subroutine evolve0D_global_partial

  ! =======================================================================

  !> Calculates the evolution of the ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolve0D_global_full(dt,pos,conv_flag,phase)

    ! Calculates the evolution of the hydrogen + helium ionization state for
    ! one cell (pos) and multiple sources.

    ! Author: Garrelt Mellema

    ! Date: 11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence.

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: conv_flag ! convergence counter
    character(len=1),intent(in) :: phase

    integer :: nx,nit ! loop counters
    integer :: conv_flag_check
    real(kind=dp) :: de ! electron density
!    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0 ! ionization fractions
    real(kind=dp) :: yh0_av_old, yhe0_av_old, yhe1_av_old, yhe2_av_old,yh1_av_old 
    real(kind=dp) :: avg_temper, temper ! temperature
    real(kind=dp) :: ndens_p ! local number density
    real(kind=dp) :: temper_old   
    real(kind=dp) :: temper1   
    real(kind=dp) :: temp_av_old,temp_av_new,temper_inter
    

    real(kind=dp) :: phih ! local H photo-ionization rate (only non-zero when local=.false.!)
    real(kind=dp) :: phih_total ! local total photo-ionization rate (including
                                ! photon loss term)
    real(kind=dp),dimension(0:1) :: phihe ! local He photo-ionization rate (only non-zero when local=.false.!)
    real(kind=dp),dimension(0:1) :: phihe_total ! local total photo-ionization rate (including
                                ! photon loss term)
    real(kind=dp),dimension(0:1) :: dummy=(/0.0_dp,0.0_dp/) ! dummy array
    real(kind=dp) :: convergence
    type(ionstates) :: ion    
    type(photrates) :: phi 
    type(photrates) :: bb_phi 
#ifdef PL
    type(photrates) :: pl_phi
#endif
#ifdef QUASARS
    type(photrates) :: qpl_phi
#endif

    ! Initialize conv_flag_check to zero
    conv_flag_check=0

    ! Initialize ionization rate phi to zero
    call set_photrates_to_zero (phi)

    ! Initialize local ionization states to global ones
    do nx=0,1
       ion%h(nx)=max(epsilon,xh_intermed(pos(1),pos(2),pos(3),nx))
       ion%h_old(nx)=max(epsilon,xh(pos(1),pos(2),pos(3),nx))
       ion%h_av(nx)=max(epsilon,xh_av(pos(1),pos(2),pos(3),nx)) ! use calculated xh_av
    enddo
       
    do nx=0,2
       ion%he(nx)=max(epsilon,xhe_intermed(pos(1),pos(2),pos(3),nx))
       ion%he_old(nx)=max(epsilon,xhe(pos(1),pos(2),pos(3),nx))
       ion%he_av(nx)=max(epsilon,xhe_av(pos(1),pos(2),pos(3),nx)) ! use calculated xhe_av
    enddo

    ! Initialize local scalars for density and temperature
    ndens_p=ndens(pos(1),pos(2),pos(3))
    call get_temperature_point (pos(1),pos(2),pos(3),temper_inter, &
         temp_av_old,temper_old)
    
    ! Use the collected photo-ionization rates
    bb_phi%photo_cell_HI=bb_phih_grid(pos(1),pos(2),pos(3))
    bb_phi%photo_cell_HeI=bb_phihe_grid(pos(1),pos(2),pos(3),0)
    bb_phi%photo_cell_HeII=bb_phihe_grid(pos(1),pos(2),pos(3),1)
    if (.not.isothermal) bb_phi%heat=bb_phiheat_grid(pos(1),pos(2),pos(3))

#ifdef QUASARS
    qpl_phi%photo_cell_HI=qpl_phih_grid(pos(1),pos(2),pos(3))
    qpl_phi%photo_cell_HeI=qpl_phihe_grid(pos(1),pos(2),pos(3),0)
    qpl_phi%photo_cell_HeII=qpl_phihe_grid(pos(1),pos(2),pos(3),1)
    if (.not.isothermal) qpl_phi%heat=qpl_phiheat_grid(pos(1),pos(2),pos(3))
#endif

#ifdef PL
    pl_phi%photo_cell_HI=pl_phih_grid(pos(1),pos(2),pos(3))
    pl_phi%photo_cell_HeI=pl_phihe_grid(pos(1),pos(2),pos(3),0)
    pl_phi%photo_cell_HeII=pl_phihe_grid(pos(1),pos(2),pos(3),1)
    if (.not.isothermal) pl_phi%heat=pl_phiheat_grid(pos(1),pos(2),pos(3))
#endif    

    ! Construct phi taking into account the hot phase

    ! All cells will be exposed to the quasar and power law rates
#ifdef QUASARS
    !phi=phi+qpl_phi
    phi%photo_cell_HI=phi%photo_cell_HI+qpl_phi%photo_cell_HI
    phi%photo_cell_HeI=phi%photo_cell_HeI+qpl_phi%photo_cell_HeI
    phi%photo_cell_HeII=phi%photo_cell_HeII+qpl_phi%photo_cell_HeII
    phi%heat=phi%heat+qpl_phi%heat
#endif
#ifdef PL
    !phi=phi+pl_phi
    phi%photo_cell_HI=phi%photo_cell_HI+pl_phi%photo_cell_HI
    phi%photo_cell_HeI=phi%photo_cell_HeI+pl_phi%photo_cell_HeI
    phi%photo_cell_HeII=phi%photo_cell_HeII+pl_phi%photo_cell_HeII
    phi%heat=phi%heat+pl_phi%heat
#endif

    ! Find out which cells are special and if to change their status
    if (special(pos(1),pos(2),pos(3)) == 0 .and. &
         bb_phi%photo_cell_HI > epsilon_dx/dt) &
         special(pos(1),pos(2),pos(3)) = 1

    if (special(pos(1),pos(2),pos(3)) /= 1) then
       ! Other cells are either not received by stellar photons, or are
       ! fully ionized, or are recombining... here we add bb rates
       !phi=phi+bb_phi
       phi%photo_cell_HI=phi%photo_cell_HI+bb_phi%photo_cell_HI
       phi%photo_cell_HeI=phi%photo_cell_HeI+bb_phi%photo_cell_HeI
       phi%photo_cell_HeII=phi%photo_cell_HeII+bb_phi%photo_cell_HeII
       phi%heat=phi%heat+bb_phi%heat
    endif

    ! Test if there are stellar photons and also check whether
    ! the fraction of the cell which is a stellar HII region ("hot") is small
    !if (bb_phi%photo_cell_HI > epsilon_dx/dt .and. &
    !     xh_hot(pos(1),pos(2),pos(3)) <= limit_partial_cells) then
       ! This cell is partly "hot" and needs to be processed separately
       ! The bb photons have already been used to make the hot HII region,
       ! here we only process the "cold" phase which is affected by x-rays
       !special(pos(1),pos(2),pos(3))=.true.
    !elseif (bb_phi%photo_cell_HI < epsilon_dx/dt .and. special(pos(1),pos(2),pos(3))=.true.)
    !else
       ! Other cells are either not received by stellar photons, or are
       ! fully ionized, or are recombining... here we add bb rates
       !phi=phi+bb_phi
    !   phi%photo_cell_HI=phi%photo_cell_HI+bb_phi%photo_cell_HI
    !   phi%photo_cell_HeI=phi%photo_cell_HeI+bb_phi%photo_cell_HeI
    !   phi%photo_cell_HeII=phi%photo_cell_HeII+bb_phi%photo_cell_HeII
    !   phi%heat=phi%heat+bb_phi%heat
    !   special(pos(1),pos(2),pos(3))=.false.
    !endif
    
    ! I think instead of calling here twice get_temp, it is perhaps better to pass t_new
    ! and t_old as arguments from/to do_chemistry. (?)
    call do_chemistry (dt, ndens_p, ion, phi, 0.0_dp, &
         dummy,  1.0_dp, 0.0_dp, pos, 0 , .false., phase)

    ! Test for global convergence using the time-averaged neutral fraction.
    ! For low values of this number assume convergence
    yh0_av_old=xh_av(pos(1),pos(2),pos(3),0) ! use previously calculated xh_av
    yh1_av_old=xh_av(pos(1),pos(2),pos(3),1)
    yhe0_av_old=xhe_av(pos(1),pos(2),pos(3),0)
    yhe1_av_old=xhe_av(pos(1),pos(2),pos(3),1)
    yhe2_av_old=xhe_av(pos(1),pos(2),pos(3),2)

    ! Check hydrogen neutral fraction for convergence
    if ( (abs((ion%h_av(0)-yh0_av_old)) > minimum_fractional_change &
         .and. &
         abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) > minimum_fractional_change &
         .and. &
         (ion%h_av(0) > minimum_fraction_of_atoms)  )) conv_flag_check=1
    
    ! Check helium neutral fractions for convergence
    if ( (abs((ion%he_av(0)-yhe0_av_old)) > minimum_fractional_change &
         .and. &
         abs((ion%he_av(0)-yhe0_av_old)/ion%he_av(0)) > minimum_fractional_change &
         .and. &
         (ion%he_av(0) > minimum_fraction_of_atoms)  ) ) conv_flag_check=1

    ! Check helium 2+ fractions for convergence
    if ((abs((ion%he_av(2)-yhe2_av_old)) > minimum_fractional_change &
         .and. &
         abs((ion%he_av(2)-yhe2_av_old)/ion%he_av(2)) > minimum_fractional_change &
         .and. &
         (ion%he_av(2) > minimum_fraction_of_atoms)  ) ) conv_flag_check=1
    
    ! Check temperature for convergence
    call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,temp_av_new,temper_old)
    if (abs((temp_av_old-temp_av_new)/temp_av_new) > 1.0e-1_dp &
         .and. &
         (abs(temp_av_new-temp_av_old) > 100.0_dp)) conv_flag_check=1
    
    ! Add local check to convergence flag. If any of the above checks failed
    ! the conv_flag will be increased by 1
    conv_flag=conv_flag+conv_flag_check
    

    ! Copy ion fractions to the global arrays.
    do nx=0,1
       xh_intermed(pos(1),pos(2),pos(3),nx)=ion%h(nx)
       xh_av(pos(1),pos(2),pos(3),nx)=ion%h_av(nx)
    enddo

    do nx=0,2
       xhe_intermed(pos(1),pos(2),pos(3),nx)=ion%he(nx)
       xhe_av(pos(1),pos(2),pos(3),nx)=ion%he_av(nx)
    enddo
    
  end subroutine evolve0D_global_full

  ! ===========================================================================

  subroutine do_chemistry (dt, ndens_p, ion, &
       phi, coldensh_in, coldenshe_in, path, vol_ph, pos, ns, local, &
       phase)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: ndens_p
    real(kind=dp),intent(in) :: coldensh_in
    real(kind=dp),dimension(0:1),intent(in) :: coldenshe_in
    real(kind=dp),intent(in) :: path
    real(kind=dp),intent(in) :: vol_ph
    integer,dimension(Ndim),intent(in) :: pos !< position on mesh
    integer,intent(in)      :: ns !< source number 
    logical,intent(in) :: local !< true if doing a non-global calculation.
    character(len=1),intent(in) :: phase

    real(kind=dp) :: avg_temper, temper0, temper1,temper2,temper_inter
    real(kind=dp) :: yh0_av_old,oldhe1av,oldhe0av,oldhav
    real(kind=dp) :: yh1_av_old
    real(kind=dp) :: yhe0_av_old
    real(kind=dp) :: yhe1_av_old
    real(kind=dp) :: yhe2_av_old
    real(kind=dp) :: de
    real(kind=dp) :: coldensh_cell
    real(kind=dp), dimension(0:1) :: coldenshe_cell, coldensheout
    real(kind=dp) :: yfrac 
    real(kind=dp) :: zfrac
    real(kind=dp) :: y2afrac 
    real(kind=dp) :: y2bfrac 

    real(kind=dp) :: ionh0old,ionh1old,ionhe0old,ionhe1old,ionhe2old
    integer :: nx
    integer :: nit

    type(photrates) :: phi
    type(ionstates) :: ion

    ! Initialize local temperature
    call get_temperature_point (pos(1),pos(2),pos(3), &
         temper_inter,avg_temper,temper1)
    temper0 = temper1
  
    ! Initialize local clumping (if type of clumping is appropriate)
    if (type_of_clumping == 5) call clumping_point (pos(1),pos(2),pos(3))
    
    nit=0
    do 
       nit=nit+1
       temper2=temper1  ! This is the temperature from last iteration

       ! Save the values of yh_av found in the previous iteration
       yh0_av_old=ion%h_av(0)
       yh1_av_old=ion%h_av(1)
       yhe0_av_old=ion%he_av(0)
       yhe1_av_old=ion%he_av(1)
       yhe2_av_old=ion%he_av(2)

       ! Calculate (mean) electron density
       de=electrondens(ndens_p,ion%h_av,ion%he_av)
       
       ! initialize the collisional ionization and recombinations rates 
       ! (temperature dependent)
       if (.not.isothermal) call ini_rec_colion_factors(avg_temper)
       
       ! Add photon losses to the photo-ionization rates
       if (add_photon_losses) then
          call distribute_photon_losses(ion,phi,ndens_p,vol)
       endif
       
       ! Calculate the new and mean ionization states
       coldensh_cell    =coldens(path,ion%h(0),ndens_p,(1.0_dp-abu_he))
       coldenshe_cell(0)=coldens(path,ion%he(0),ndens_p,abu_he)
       coldenshe_cell(1)=coldens(path,ion%he(1),ndens_p,abu_he)
       
       call prepare_doric_factors(coldensh_cell,coldenshe_cell,yfrac,zfrac, &
            y2afrac,y2bfrac)
       
       call doric(dt,de,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac)

       de=electrondens(ndens_p,ion%h_av,ion%he_av)
       
       ! Update column density for 2nd call to doric
       coldensh_cell  =coldens(path,ion%h(0),ndens_p,(1.0_dp-abu_he))
       coldenshe_cell(0)=coldens(path,ion%he(0),ndens_p,abu_he)
       coldenshe_cell(1)=coldens(path,ion%he(1),ndens_p,abu_he)
       
       ! Prepare factors needed by doric
       call prepare_doric_factors(coldensh_cell,coldenshe_cell,yfrac,zfrac, &
            y2afrac,y2bfrac)
       
       ! Save results of first pass over doric
       ionh0old=ion%h(0)
       ionh1old=ion%h(1)
       ionhe0old=ion%he(0)
       ionhe1old=ion%he(1)
       ionhe2old=ion%he(2)
       oldhav=ion%h_av(0)
       oldhe0av=ion%he_av(0)
       oldhe1av=ion%he_av(1)        
       
       call doric(dt,de,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac)
       
       ! Average the answers from the two passes over doric
       ion%h(0)=(ion%h(0)+ionh0old)/2.0_dp
       ion%h(1)=(ion%h(1)+ionh1old)/2.0_dp
       ion%he(0)=(ion%he(0)+ionhe0old)/2.0_dp
       ion%he(1)=(ion%he(1)+ionhe1old)/2.0_dp
       ion%he(2)=(ion%he(2)+ionhe2old)/2.0_dp
       ion%h_av(0)=(ion%h_av(0)+oldhav)/2.0_dp
       ion%he_av(0)=(ion%he_av(0)+oldhe0av)/2.0_dp
       ion%he_av(1)=(ion%he_av(1)+oldhe1av)/2.0_dp
       
       de=electrondens(ndens_p,ion%h_av,ion%he_av)
       !  endif
       
       temper1=temper0 
       if (.not.isothermal) &
          call thermal(dt,temper1,avg_temper,de,ndens_p, &
               ion,phi,phase)
       
       ! Test for convergence on time-averaged neutral fraction
       ! For low values of this number assume convergence
       if ((abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) < &
            minimum_fractional_change .or. &
            (ion%h_av(0) < minimum_fraction_of_atoms)).and. &
            
            !(abs((ion%h_av(1)-yh1_av_old)/ion%h_av(1)) < convergence2         &
            !.or. (ion%h_av(1) < convergence_frac)).and.                       &
            
            (abs((ion%he_av(0)-yhe0_av_old)/ion%he_av(0)) < &
            minimum_fractional_change .or. &
            (ion%he_av(0) < minimum_fraction_of_atoms)).and. &
            
	    ! (abs((ion%he_av(1)-yhe1_av_old)/ion%he_av(1)) < convergence2      &
            ! .or. (ion%he_av(1) < convergence_frac)).and.                      &
            
            (abs((ion%he_av(2)-yhe2_av_old)/ion%he_av(2)) < & 
            minimum_fractional_change .or. &
            (ion%he_av(2) < minimum_fraction_of_atoms)) .and. &
            
            (abs((temper1-temper2)/temper1) < minimum_fractional_change) & 
            ) then  
          exit
       endif
       
       ! Warn about non-convergence and terminate iteration
       if (nit > 400) then
          if (rank == 0) then   
             write(logf,*) 'Convergence failing (global) nit=', nit
             write(logf,*) 'x',ion%h_av(0), ion%he_av(0:2)
             write(logf,*) 'h',yh0_av_old, yhe0_av_old,yhe1_av_old, yhe2_av_old
             write(logf,*) abs(ion%h_av(0)-yh0_av_old),abs(ion%he_av(0)-yhe0_av_old),abs(ion%he_av(1))
          endif
          exit
       endif
    enddo

    ! Update temperature
    ! GM/130815: Why is this done here?
    if (.not. isothermal .and. phase .ne. "H") &
         call set_temperature_point (pos(1),pos(2),pos(3),temper1,avg_temper)
    
  end subroutine do_chemistry
  
  ! ===========================================================================

  ! This subroutine is supposed to add the photon losses to the photo-ionization
  ! rates. It is however not up to date as it only uses 7 frequency bands.
  ! It should not be used until this is solved.
  ! Refer to the H-only version to see another way to process photon_losses
  subroutine distribute_photon_losses(ion,phi,ndens_p,vol)

    !use radiation, only: scale_int2,scale_int3
    use radiation_photoionrates, only: scale_int2,scale_int3

    type(ionstates),intent(in) :: ion
    type(photrates),intent(inout) :: phi
    real(kind=dp),intent(in) :: ndens_p
    real(kind=dp),intent(in) :: vol


    real(kind=dp) :: Nhe0,Nhe1,Nh0, ovtotneutral
    real(kind=dp) :: N_He0_ov_H0,N_He1_ov_H0,N_He0_ov_He1,N_H0_ov_He1
    real(kind=dp) :: N_H0_ov_He0,N_He1_ov_He0
    real(kind=dp) :: fH_2a,fH_2b,fHe0_2a,fHe0_2b
    real(kind=dp) :: fH_3a,fH_3b,fH_3c,fH_3d
    real(kind=dp) :: fHe0_3a,fHe0_3b,fHe0_3c,fHe0_3d
    real(kind=dp) :: fHe1_3a,fHe1_3b,fHe1_3c,fHe1_3d

    if (sum(photon_loss) /= 0.0) then
       Nhe0= abu_he*ion%he(0)              ! *** probably use average fractions
       Nhe1= abu_he*ion%he(1)              ! *** probably use average fractions
       Nh0 = (1.0-abu_he)*ion%h(0)
       N_He0_ov_H0         = max(epsilon,Nhe0)/max(epsilon,Nh0)
       N_He0_ov_He1        = max(epsilon,Nhe0)/max(epsilon,Nhe1)
       N_H0_ov_He0         = max(epsilon,Nh0)  /max(epsilon,Nhe0)
       N_H0_ov_He1         = max(epsilon,Nh0)  /max(epsilon,Nhe1)
       N_He1_ov_He0        = max(epsilon,Nhe1)/max(epsilon,Nhe0)
       N_He1_ov_H0         = max(epsilon,Nhe1)/max(epsilon,Nh0)
       !call scale_int3(fH_3a,fHe0_3a,fHe1_3a, &
       !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
       !1,2,3)
       !call scale_int3(fH_3b,fHe0_3b,fHe1_3b, &
       !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
       !4,5,6)
       !call scale_int3(fH_3c,fHe0_3c,fHe1_3c, &
       !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
       !7,8,9)
       !call scale_int3(fH_3d,fHe0_3d,fHe1_3d, &
       !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
       !10,11,12)
       !call scale_int2(fH_2a,fHe0_2a,N_He0_ov_H0,1,2)
       !call scale_int2(fH_2b,fHe0_2b,N_He0_ov_H0,3,4)
       !*** it is not so clear to me if I should first scale the photon loss in every interval with the amount of 
       !*** "neutral" ions those photons could possibly ionize, or if I should first sclae them and ascribe them to 
       !*** either H, He0 or He+ and then scale them with only the species in question. 
       !*** for now, I'll just take the first option (scale and then devide)
       ovtotneutral= 1.0_dp/(vol*ion%h_av(0)*ndens_p*(1.0-abu_he))
       photon_loss(1)=photon_loss(1)*ovtotneutral
       ovtotneutral=1.0_dp/(vol*ndens_p*((1.0-abu_he)*ion%h_av(0)+abu_he*ion%he_av(0)))
       photon_loss(2)=photon_loss(2)*ovtotneutral
       photon_loss(3)=photon_loss(3)*ovtotneutral
       ovtotneutral= 1.0_dp/(vol*ndens_p*((1.0-abu_he)*ion%h_av(0)+abu_he*(ion%he_av(0)+ion%he_av(1))))
       photon_loss(4)=photon_loss(4)*ovtotneutral
       photon_loss(5)=photon_loss(5)*ovtotneutral
       photon_loss(6)=photon_loss(6)*ovtotneutral
       photon_loss(7)=photon_loss(7)*ovtotneutral
       phi%photo_cell_HI  = phi%photo_cell_HI   + &
            (photon_loss(1) +photon_loss(2)* fH_2a  + &
            photon_loss(4)*fH_3a+photon_loss(5)*fH_3b + & 
            photon_loss(3)* fH_2b  + &
            photon_loss(6)*fH_3c + &
            photon_loss(7)*fH_3d) 
       phi%photo_cell_HeI= phi%photo_cell_HeI + &
            (photon_loss(2)*fHe0_2a  + &
            photon_loss(4)* fHe0_3a + &
            photon_loss(5)* fHe0_3b + &
            photon_loss(3)*fHe0_2b + &
            photon_loss(6)* fHe0_3c + &
            photon_loss(7)* fHe0_3d)
       phi%photo_cell_HeII= phi%photo_cell_HeII + &
            (photon_loss(4)*fHe1_3a + &
            photon_loss(5)*fHe1_3b + &
            photon_loss(6)*fHe1_3c + &
            photon_loss(7)*fHe1_3d)
    endif
    
  end subroutine distribute_photon_losses
  
end module evolve_point
