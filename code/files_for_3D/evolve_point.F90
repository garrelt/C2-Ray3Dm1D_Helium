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
  use material, only: ndens, xh,xhe
  use material, only: get_temperature_point, set_temperature_point
  use material, only: set_final_temperature_point, isothermal
  use material, only: ionstates
  use material, only: clumping_point
  use material, only: coldensh_LLS, LLS_point
  use sourceprops, only: srcpos
  use radiation_photoionrates, only:photrates,photoion_rates
    !individual_photoion_rates,&
  use thermalevolution, only: thermal
  use photonstatistics, only: photon_loss, total_LLS_loss
  use tped, only: electrondens
  use doric_module, only: doric, prepare_doric_factors, coldens
#if defined(QUASARS) && defined(PL)
  use evolve_data, only: phih_grid, phihe_grid, phiheat, pl_phih_grid,&
   pl_phiheat, qpl_phih_grid, qpl_phiheat, bb_phih_grid, bb_phiheat
#elif defined(QUASARS)
  use evolve_data, only: phih_grid, phihe_grid, phiheat,qpl_phih_grid,&
   qpl_phiheat, bb_phih_grid, bb_phiheat
#elif defined(PL)
  use evolve_data, only: phih_grid, phihe_grid, phiheat, pl_phih_grid,&
    pl_phiheat, bb_phih_grid, bb_phiheat
#else
  use evolve_data, only: phih_grid, phihe_grid, phiheat
#endif
  use evolve_data, only: xh_av, xhe_av, xh_intermed, xhe_intermed
  use evolve_data, only: coldensh_out, coldenshe_out
  use evolve_data, only: photon_loss_src_thread
  use evolve_data, only: last_l,last_r
  use evolve_data, only: tn

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
  subroutine evolve0D(dt,rtpos,ns,niter,sourcetype)
    
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
    character(len=1),intent(in) :: sourcetype
    
    integer :: nx,nd,idim ! loop counters
    integer,dimension(Ndim) :: pos
    real(kind=dp) :: dist2,path,vol_ph
    real(kind=dp) :: xs,ys,zs
    real(kind=dp) :: coldensh_in
    real(kind=dp),dimension(0:1) :: coldenshe_in,coldenshe_out_temp
    real(kind=dp) :: ndens_p
    
    type(photrates) :: phi, dummiphi
#if defined(PL) || defined(QUASARS)
    real(kind=dp),dimension(2) :: bb_phi
#endif
#ifdef PL
    real(kind=dp),dimension(2) :: pl_phi
#endif
#ifdef QUASARS
    real(kind=dp),dimension(2) :: qpl_phi
#endif
    type(ionstates) :: ion

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
          ion%h(nx)=max(xh_intermed(pos(1),pos(2),pos(3),nx),epsilon)
          ion%h_old(nx)=max(xh(pos(1),pos(2),pos(3),nx),epsilon)         
       enddo

       do nx=0,2
          ion%he_av(nx)=max(xhe_av(pos(1),pos(2),pos(3),nx),epsilon)
          ion%he(nx)=max(xhe_intermed(pos(1),pos(2),pos(3),nx),epsilon)
          ion%he_old(nx)=max(xhe(pos(1),pos(2),pos(3),nx),epsilon)         
       enddo
      
       ! Initialize local density and temperature
       ndens_p=ndens(pos(1),pos(2),pos(3))
!       call get_temperature_point(pos(1),pos(2),pos(3),temper)
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

       ! Only ray trace and exit. Do not touch the ionization
       ! fractions. They are updated using phih_grid in evolve0d_global.
       ! Only do chemistry if this is the first pass over the sources,
       ! and if column density is below the maximum.
       ! On the first global iteration pass it may be beneficial to assume 
       ! isolated sources, but on later passes the effects of multiple sources 
       ! has to be taken into account. 
       ! Therefore no changes to xh, xh_av, etc. should happen on later passes!
       ! This option is temporarily disabled by testing niter == -1 which
       ! is always false.
       if (niter == -1 .and. coldensh_in < max_coldensh) then
          local_chemistry=.true.
          dummiphi%photo_cell_HI=0.0_dp

          call do_chemistry (dt, ndens_p, ion, dummiphi, &
               coldensh_in,coldenshe_in, path, vol_ph, pos, ns, local=.true.)

          ! Copy ion fractions to global arrays.
          ! This will speed up convergence if
          ! the sources are isolated and only ionizing up.
          ! In other cases it does not make a difference.
          xh_intermed(pos(1),pos(2),pos(3),1)=max(ion%h(1), &
               xh_intermed(pos(1),pos(2),pos(3),1))
          xh_intermed(pos(1),pos(2),pos(3),0)=max(epsilon,1.0- &
               xh_intermed(pos(1),pos(2),pos(3),1))
          xh_av(pos(1),pos(2),pos(3),1)=max(ion%h_av(1), &
               xh_av(pos(1),pos(2),pos(3),1))
          xh_av(pos(1),pos(2),pos(3),0)=max(epsilon,1.0_dp- &
               xh_av(pos(1),pos(2),pos(3),1))

          xhe_intermed(pos(1),pos(2),pos(3),0)=min(ion%he(0), &
               xhe_intermed(pos(1),pos(2),pos(3),0))
          xhe_intermed(pos(1),pos(2),pos(3),2)=max(ion%he(2), &
               xhe_intermed(pos(1),pos(2),pos(3),2))
          xhe_intermed(pos(1),pos(2),pos(3),0)=max(epsilon,1.0_dp- &
               xhe_intermed(pos(1),pos(2),pos(3),0) - &
               xhe_intermed(pos(1),pos(2),pos(3),2))

          xhe_av(pos(1),pos(2),pos(3),0)=min(ion%he_av(0), &
               xhe_av(pos(1),pos(2),pos(3),0))
          xhe_av(pos(1),pos(2),pos(3),2)=max(ion%he_av(2), &
               xhe_av(pos(1),pos(2),pos(3),2))

          xhe_av(pos(1),pos(2),pos(3),1)=max(epsilon,1.0_dp- &
               xhe_av(pos(1),pos(2),pos(3),0) - &
               xhe_av(pos(1),pos(2),pos(3),2))          
       endif  !if niter==!

       ! Add the (time averaged) column density of this cell
       ! to the total column density (for this source)
       ! and add the LLS column density to this.
       ! GM/110224: No! This messes up phi since phi is based
       !  upon the difference between the in and out column density.
       !  Instead add the LLS to coldensh_in, see above
       coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
            coldens(path,ion%h_av(0),ndens_p,(1.0_dp-abu_he))

       coldenshe_out(pos(1),pos(2),pos(3),0)=coldenshe_in(0) + &
            coldens(path,ion%he_av(0),ndens_p,abu_he)

       coldenshe_out(pos(1),pos(2),pos(3),1)=coldenshe_in(1) + &
            coldens(path,ion%he_av(1),ndens_p,abu_he)
       
       ! Calculate (photon-conserving) photo-ionization rate from the
       ! column densities.

       ! Limit the calculation to a certain maximum column density (hydrogen)
       if (coldensh_in < max_coldensh) then 

          coldenshe_out_temp=coldenshe_out(pos(1),pos(2),pos(3),:)
          ! photoion_rates the structure of rates (photo and heating)
          phi=photoion_rates(coldensh_in,coldensh_out(pos(1),pos(2),pos(3)), &
			 coldenshe_in(0),coldenshe_out_temp(0), &
			 coldenshe_in(1),coldenshe_out_temp(1), &
			 vol_ph,ns,ion%h_av(1),sourcetype)

#if defined(PLs) || defined(QUASARS)
          if (sourcetype == "B") then
             bb_phi(1)=phi%bb_photo_cell_HI
             bb_phi(2)=phi%bb_heat
          endif
#endif

#ifdef PL
          if (sourcetype == "P") then
             pl_phi(1)=phi%pl_photo_cell_HI
             pl_phi(2)=phi%pl_heat
          endif
#endif

#ifdef QUASARS
          if (sourcetype = "Q") then
             qpl_phi(1)=phi%qpl_photo_cell_HI
             qpl_phi(2)=phi%qpl_heat
          endif
#endif

          ! Divide the photo-ionization rates by the appropriate neutral density
          ! (part of the photon-conserving rate prescription)
          phi%photo_cell_HI=phi%photo_cell_HI/(ion%h_av(0)*ndens_p*(1.0_dp-abu_he))
          phi%photo_cell_HeI=phi%photo_cell_HeI/(ion%he_av(0)*ndens_p*abu_he)
          phi%photo_cell_HeII=phi%photo_cell_HeII/(ion%he_av(1)*ndens_p*abu_he)

#if defined(PLs) || defined(QUASARS)
          if (source_type == "B") &
             bb_phi(1)=bb_phi(1)/(ion%h_av(0)*ndens_p*(1.0_dp-abu_he))
             
#endif

#ifdef PL
          if (source_type == "P") &
               pl_phi(1)=pl_phi(1)/(ion%h_av(0)*ndens_p*(1.0_dp-abu_he))
#endif

#ifdef QUASARS
          if (source_type == "Q") &
               qpl_phi(1)=qpl_phi(1)/(ion%h_av(0)*ndens_p*(1.0_dp-abu_he))
#endif
          
          ! Calculate the losses due to LLSs.
          ! GM/110224: Add the factor vol/vol_ph to phi, just as we do for
          ! the photon losses.
          ! GM/110302: Use phi%h_in (not out) here since this is where we draw
          ! the photons from, above.
          if (use_LLS) call total_LLS_loss(phi%photo_in_HI*vol/vol_ph, &
               coldensh_LLS * path/dr(1))          

       else

          ! If the H0 column density is above the maximum, set rates to zero
          phi%photo_cell_HI = 0.0_dp
          phi%photo_cell_HeI = 0.0_dp
          phi%photo_cell_HeII = 0.0_dp
          phi%photo_out_HI = 0.0_dp
          phi%photo_out_HeI = 0.0_dp
          phi%photo_out_HeII = 0.0_dp
          phi%heat = 0.0_dp
          phi%photo_in = 0.0_dp
          phi%photo_out = 0.0_dp

#if defined(PLs) || defined(QUASARS)
          if (source_type == "B") then
             bb_phi(1) = 0.0_dp
             bb_phi(2) = 0.0_dp
          endif
#endif

#ifdef PL
          if (source_type == "P") then
             pl_phi(1) = 0.0_dp
             pl_phi(2) = 0.0_dp
          endif
#endif

#ifdef QUASARS
          if (source_type == "Q") then
             qpl_phi(1) = 0.0_dp
             qpl_phi(2) = 0.0_dp
          endif
#endif
       endif
       
       phih_grid(pos(1),pos(2),pos(3))= &
            phih_grid(pos(1),pos(2),pos(3))+phi%photo_cell_HI
       phihe_grid(pos(1),pos(2),pos(3),0)=&
             phihe_grid(pos(1),pos(2),pos(3),0)+phi%photo_cell_HeI
       phihe_grid(pos(1),pos(2),pos(3),1)=&
             phihe_grid(pos(1),pos(2),pos(3),1)+phi%photo_cell_HeII

#if defined(PLs) || defined(QUASARS)
       if (source_type == "B") &
            bb_phih_grid(pos(1),pos(2),pos(3))= &
            bb_phih_grid(pos(1),pos(2),pos(3))+bb_phi(1)
#endif

#ifdef PL
       if (source_type == "P") &
            pl_phih_grid(pos(1),pos(2),pos(3))= &
            pl_phih_grid(pos(1),pos(2),pos(3))+pl_phi(1)
#endif

#ifdef QUASARS
       if (source_type == "Q") &
            qpl_phih_grid(pos(1),pos(2),pos(3))= &
            qpl_phih_grid(pos(1),pos(2),pos(3))+qpl_phi(1)
#endif   
       if (.not. isothermal) then
          phiheat(pos(1),pos(2),pos(3))=phiheat(pos(1),pos(2),pos(3))+phi%heat

#if defined(PLs) || defined(QUASARS)
          if (source_type == "B") &
               bb_phiheat(pos(1),pos(2),pos(3))=bb_phiheat(pos(1),pos(2),pos(3))+bb_phi(2)
#endif

#ifdef PL
          if (source_type == "P") &
               pl_phiheat(pos(1),pos(2),pos(3))=pl_phiheat(pos(1),pos(2),pos(3))+pl_phi(2)
#endif

#ifdef QUASARS
          if (source_type == "Q") & 
               qpl_phiheat(pos(1),pos(2),pos(3))=qpl_phiheat(pos(1),pos(2),pos(3))+qpl_phi(2)
#endif
       endif
       
       ! Photon statistics: register number of photons leaving the grid
       ! Note: This is only the H0 photo-ionization rate
       if ( (any(rtpos(:) == last_l(:))) .or. &
            (any(rtpos(:) == last_r(:))) ) then
          photon_loss_src_thread(tn)=photon_loss_src_thread(tn) + &
               phi%photo_out*vol/vol_ph
          !photon_loss_src(1,tn)=photon_loss_src(1,tn) + phi%h_out*vol/vol_ph
       endif

    endif ! end of coldens test
    
  end subroutine evolve0D

  ! =======================================================================

  !> Calculates the evolution of the ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolve0D_global(dt,pos,conv_flag,sourcetype)

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
    character(len=1),intent(in) :: sourcetype

    integer :: nx,nit ! loop counters
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
#if defined(QUASARS) || defined(PL)
    real(kind=dp),dimension(2) :: bb_phi
#endif
#ifdef PL
    real(kind=dp),dimension(2) :: pl_phi 
#endif
#ifdef QUASARS
    real(kind=dp),dimension(2) :: qpl_phi 
#endif

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
    !avg_temper=temper

    ! Use the collected photo-ionization rates
    if (source_type == "A") then
       phi%photo_cell_HI=phih_grid(pos(1),pos(2),pos(3))
       phi%photo_cell_HeI=phihe_grid(pos(1),pos(2),pos(3),0)
       phi%photo_cell_HeII=phihe_grid(pos(1),pos(2),pos(3),1)
    endif
#if defined(QUASARS) || defined(PL)
    if (source_type == "B") then
       phi%photo_cell_HI=bb_phih_grid(pos(1),pos(2),pos(3))
       phi%photo_cell_HeI=bb_phihe_grid(pos(1),pos(2),pos(3),0)
       phi%photo_cell_HeII=bb_phihe_grid(pos(1),pos(2),pos(3),1)
    endif
#endif
#ifdef PL
    if (source_type == "P") then
       phi%photo_cell_HI=pl_phih_grid(pos(1),pos(2),pos(3))
       phi%photo_cell_HeI=pl_phihe_grid(pos(1),pos(2),pos(3),0)
       phi%photo_cell_HeII=pl_phihe_grid(pos(1),pos(2),pos(3),1)
    endif
#endif
#ifdef QUASARS
    if (source_type == "Q") then
       phi%photo_cell_HI=qpl_phih_grid(pos(1),pos(2),pos(3))
       phi%photo_cell_HeI=qpl_phihe_grid(pos(1),pos(2),pos(3),0)
       phi%photo_cell_HeII=qpl_phihe_grid(pos(1),pos(2),pos(3),1)
    endif
#endif
             
    if(.not.isothermal) then 
        phi%heat=phiheat(pos(1),pos(2),pos(3))
#if defined(QUASARS) || defined(PL)
        if (source_type == "B") phi%heat=bb_phiheat(pos(1),pos(2),pos(3))
#endif
#ifdef PL
        if (source_type == "P") phi%heat=pl_phiheat(pos(1),pos(2),pos(3))
#endif
#ifdef QUASARS
        if (source_type == "Q") qpl_phi(2)=qpl_phiheat(pos(1),pos(2),pos(3))
#endif
    endif
    ! I think instead of calling here twice get_temp, it is perhaps better to pass t_new
    ! and t_old as arguments from/to do_chemistry. (?)
    call do_chemistry (dt, ndens_p, ion, phi, 0.0_dp, &
         dummy,  1.0_dp, 0.0_dp, pos, 0 , local=.false., source_type)
    
    ! Test for global convergence using the time-averaged neutral fraction.
    ! For low values of this number assume convergence
    yh0_av_old=xh_av(pos(1),pos(2),pos(3),0) ! use previously calculated xh_av
    yh1_av_old=xh_av(pos(1),pos(2),pos(3),1)
    yhe0_av_old=xhe_av(pos(1),pos(2),pos(3),0)
    yhe1_av_old=xhe_av(pos(1),pos(2),pos(3),1)
    yhe2_av_old=xhe_av(pos(1),pos(2),pos(3),2)
    call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,temp_av_new,temper_old)

    if ( (abs((ion%h_av(0)-yh0_av_old)) > minimum_fractional_change                .and. &
          abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) > minimum_fractional_change   .and. &
              (ion%h_av(0) > minimum_fraction_of_atoms)  ).or.                       &
         (abs((ion%he_av(0)-yhe0_av_old)) > minimum_fractional_change .and. &
          abs((ion%he_av(0)-yhe0_av_old)/ion%he_av(0)) > minimum_fractional_change .and. &
              (ion%he_av(0) > minimum_fraction_of_atoms)  ).or.                      &
    !     (abs((ion%he_av(1)-yhe1_av_old)) > minimum_fractional_change .and. &
    !      abs((ion%he_av(1)-yhe1_av_old)/ion%he_av(1)) > minimum_fractional_change .and. &
    !          (ion%he_av(1) >  minimum_fraction_of_atoms)  ).or.                      & 
         (abs((ion%he_av(2)-yhe2_av_old)) > minimum_fractional_change              .and. &
          abs((ion%he_av(2)-yhe2_av_old)/ion%he_av(2)) > minimum_fractional_change .and. &
              (ion%he_av(2) > minimum_fraction_of_atoms)  ).or.                & 
         !(abs((temper1-temper_old)/temper1) > 1.0e-1_dp).and.              &
         !(abs(temper1-temper_old) >     100.0_dp)                          &
         (abs((temp_av_old-temp_av_new)/temp_av_new) > 1.0e-1_dp).and.              &
         (abs(temp_av_new-temp_av_old) >     100.0_dp)                          &                  
                                                                      ) then
       conv_flag=conv_flag+1
    endif

    ! Copy ion fractions to the global arrays.
    do nx=0,1
       xh_intermed(pos(1),pos(2),pos(3),nx)=ion%h(nx)
       xh_av(pos(1),pos(2),pos(3),nx)=ion%h_av(nx)
    enddo

    do nx=0,2
       xhe_intermed(pos(1),pos(2),pos(3),nx)=ion%he(nx)
       xhe_av(pos(1),pos(2),pos(3),nx)=ion%he_av(nx)
    enddo
    
    ! This was already done in do_chemistry
    !if (.not.isothermal) call set_temperature_point (pos(1),pos(2),pos(3),temper1,av)
    
  end subroutine evolve0D_global

  ! ===========================================================================

  subroutine do_chemistry (dt, ndens_p, ion, &
       phi, coldensh_in, coldenshe_in, path, vol_ph, pos, ns, local, &
       source_type)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: ndens_p
    real(kind=dp),intent(in) :: coldensh_in
    real(kind=dp),dimension(0:1),intent(in) :: coldenshe_in
    real(kind=dp),intent(in) :: path
    real(kind=dp),intent(in) :: vol_ph
    integer,dimension(Ndim),intent(in) :: pos !< position on mesh
    integer,intent(in)      :: ns !< source number 
    logical,intent(in) :: local !< true if doing a non-global calculation.
    character(len=1),intent(in) :: source_type

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
    call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,avg_temper,temper1)
    !avg_temper=temper1
    temper0   =temper1
  
    ! Initialize local clumping (if type of clumping is appropriate)
    if (type_of_clumping == 5) call clumping_point (pos(1),pos(2),pos(3))
    
    nit=0
    do 
       nit=nit+1
       temper2   =temper1  ! This is the temperature from last iteration

       ! Save the values of yh_av found in the previous iteration
       yh0_av_old=ion%h_av(0)
       yh1_av_old=ion%h_av(1)
       yhe0_av_old=ion%he_av(0)
       yhe1_av_old=ion%he_av(1)
       yhe2_av_old=ion%he_av(2)

       ! Copy ionic abundances back to initial values (doric assumes
       ! that it contains this) !*** this is not longer needed because I have the
       ! *** old values in ion%
       !yh(:)=yh0(:)
              
       ! Calculate (mean) electron density
       de=electrondens(ndens_p,ion%h_av,ion%he_av)

       ! Find total photo-ionization rate
       if (local) then
          
          ! Calculate (time averaged) column density of cell
          coldensh_cell=coldens(path,ion%h_av(0),ndens_p,1.0_dp-abu_he)
          do nx=0,1
             coldenshe_cell(nx)=coldens(path,ion%he_av(nx),ndens_p,abu_he)
          enddo
          coldensheout(:)=coldenshe_in(:)+coldenshe_cell(:)

          ! Calculate (photon-conserving) photo-ionization rate
          phi=photoion_rates(coldensh_in,coldensh_in+coldensh_cell, &
               coldenshe_in(0), coldensheout(0), &
               coldenshe_in(1), coldensheout(1), &
               vol_ph,ns,ion%h_av(1))
          
          phi%photo_cell_HI=phi%photo_cell_HI/(ion%h_av(0)*ndens_p* &
               (1.0_dp-abu_he))
          phi%photo_cell_HeI=phi%photo_cell_HeI/(ion%he_av(nx)*ndens_p*abu_he)
          phi%photo_cell_HeII=phi%photo_cell_HeII/(ion%he_av(nx)*ndens_p*abu_he)
          !do nx=0,1
          !   phi%he(nx)=phi%he(nx)/(ion%he_av(nx)*ndens_p*abu_he)
          !enddo

       else

          ! (direct plus photon losses)
          ! DO THIS HERE, yh_av is changing
          ! (if the cell is ionized, add a fraction of the lost photons)
          !if (xh_intermed(pos(1),pos(2),pos(3),1) > 0.5)
          !phi%photo_cell_HI=phih_cell
          !phi%photo_cell_HeI=phihe_cell(0)
          !phi%photo_cell_HeII=phihe_cell(1)
          !phi%heat=phihv_cell

          ! initialize the collisional ionization and recombinations rates 
          ! (temperature dependent)
          if (.not.isothermal) call ini_rec_colion_factors(avg_temper) 
          
          ! Add photon losses to the photo-ionization rates
          if (add_photon_losses) then
             call distribute_photon_losses(ion,phi,ndens_p,vol)
          endif

       endif     ! local if/else

       !      Calculate the new and mean ionization states
       !*** DO THIS ONLY IF PHI IS NOT ==0
       ! if (phi%h.ne.0.0_dp) then
       
       coldensh_cell    =coldens(path,ion%h(0),ndens_p,(1.0_dp-abu_he))
       coldenshe_cell(0)=coldens(path,ion%he(0),ndens_p,abu_he)
       coldenshe_cell(1)=coldens(path,ion%he(1),ndens_p,abu_he)
       
       call prepare_doric_factors(coldensh_cell,coldenshe_cell,yfrac,zfrac, &
            y2afrac,y2bfrac)
       
       call doric(dt,de,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac, &
            source_type)!,local)! 
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

       call doric(dt,de,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac, &
            source_type)!,local)!

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
            ion,phi,source_type)    
       
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
    if (.not. isothermal) call set_temperature_point (pos(1),pos(2),pos(3),temper1,avg_temper)
    
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
