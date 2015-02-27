module AT_photoionization

  use my_mpi
  use precision, only: dp
  use array, only: global_photoionization_HI_array, &
                   global_photoionization_HeI_array, global_photoionization_HeII_array, &
                   number_density_array, xHI_array, xHII_array, xHeI_array, xHeII_array, temperature_array, &
                   non_screened_Ifront_array
  use parameter, only: mesh, half, do_periodic,&
                       length_pos_x, length_neg_x, &
                       length_pos_y, length_neg_y, &
                       length_pos_z, length_neg_z
  use grid, only: cell_size
  use photo_rate, only: photoion_rates
  use type_definition, only: photrates
  use doric_module, only: coldens
  use mathconstant, only: pi
  use sourceprops, only: srcpos
  use abundances, only: abu_he

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_get_photoionization_rate(ns)

    implicit none

    integer, intent(in) :: ns
    integer :: i,j,k
    integer :: i_global_frame
    integer :: j_global_frame
    integer :: k_global_frame
    ! Incoming column density
    real(kind=dp) :: colum_in_HI
    ! Incoming column density
    real(kind=dp) :: colum_in_HeI
    ! Incoming column density
    real(kind=dp) :: colum_in_HeII
    ! outgoing column density
    real(kind=dp) :: colum_out_HI
    ! outgoing column density
    real(kind=dp) :: colum_out_HeI
    ! outgoing column density
    real(kind=dp) :: colum_out_HeII
    ! Size of a cell
    real(kind=dp) :: dr
    ! Electron density of the cell
    real(kind=dp) :: ndens_electron_cell
    ! Ionized fractions of I-front
    real(kind=dp) :: xHI
    real(kind=dp) :: xHII
    real(kind=dp) :: xHeI
    real(kind=dp) :: xHeII
    ! Temperature of I-front
    real(kind=dp) :: temperature
    ! Atom density of I-front
    real(kind=dp) :: ndens_atom_cell
    ! Volume of the cell
    real(kind=dp) :: vol_shell
    ! distance from the source
    real(kind=dp) :: distance_from_source
    ! distance traversed in the cell
    real(kind=dp) :: distance_traverse
    ! Photoionization rate and heating rate
    type(photrates) :: phi

    if (non_screened_Ifront_array(0,0,0).eq.1) then

      i_global_frame = srcpos(1,ns)
      j_global_frame = srcpos(2,ns)
      k_global_frame = srcpos(3,ns)

      ndens_atom_cell = number_density_array(i_global_frame,j_global_frame,k_global_frame)
      xHI = xHI_array(i_global_frame,j_global_frame,k_global_frame)
      xHII = xHII_array(i_global_frame,j_global_frame,k_global_frame)
      xHeI = xHeI_array(i_global_frame,j_global_frame,k_global_frame)
      xHeII = xHeII_array(i_global_frame,j_global_frame,k_global_frame)
      temperature = temperature_array(i_global_frame,j_global_frame,k_global_frame)
      dr = cell_size
      distance_traverse = half*dr

      colum_in_HI = 0
      colum_out_HI = colum_in_HI+coldens(distance_traverse,xHI,ndens_atom_cell,1.0_dp-abu_he)
 
      colum_in_HeI = 0
      colum_out_HeI = colum_in_HeI+coldens(distance_traverse,xHeI,ndens_atom_cell,abu_he)

      colum_in_HeII = 0
      colum_out_HeII = colum_in_HeII+coldens(distance_traverse,xHeII,ndens_atom_cell,abu_he)

      vol_shell = cell_size*cell_size*cell_size

      call photoion_rates (phi,colum_in_HI,colum_out_HI,colum_in_HeI,colum_out_HeI, &
	       	       colum_in_HeII,colum_out_HeII,vol_shell,ns,xHII)


      phi%photo_cell_HI = phi%photo_cell_HI/(xHI*ndens_atom_cell*(1.0_dp-abu_he))
      phi%photo_cell_HeI = phi%photo_cell_HeI/(xHeI*ndens_atom_cell*abu_he)
      phi%photo_cell_HeII = phi%photo_cell_HeII/(xHeII*ndens_atom_cell*abu_he)

      global_photoionization_HI_array(i_global_frame,j_global_frame,k_global_frame) = &
          global_photoionization_HI_array(i_global_frame,j_global_frame,k_global_frame) + phi%photo_cell_HI
      global_photoionization_HeI_array(i_global_frame,j_global_frame,k_global_frame) = &
          global_photoionization_HeI_array(i_global_frame,j_global_frame,k_global_frame) + phi%photo_cell_HeI
      global_photoionization_HeII_array(i_global_frame,j_global_frame,k_global_frame) = &
          global_photoionization_HeII_array(i_global_frame,j_global_frame,k_global_frame) + phi%photo_cell_HeII
   
    else

!$OMP PARALLEL PRIVATE(i,j,k,i_global_frame,j_global_frame,k_global_frame,ndens_atom_cell,xHI,xHII,xHeI,xHeII,temperature,&
!$OMP &                dr,colum_in_HI,colum_in_HeI,colum_in_HeII,distance_from_source,distance_traverse,vol_shell,&
!$OMP &                colum_out_HI,colum_out_HeI,colum_out_HeII,phi)
!$OMP DO SCHEDULE(DYNAMIC,1)

      ! true for periodic case only
      do i = length_neg_x, length_pos_x
        do j = length_neg_y, length_pos_y
          do k = length_neg_z, length_pos_z
            if (non_screened_Ifront_array(i,j,k).eq.1) then
              ! transformation from source frame to global frame
              if (do_periodic) then
                i_global_frame = modulo(i+srcpos(1,ns)-1,mesh)+1
                j_global_frame = modulo(j+srcpos(2,ns)-1,mesh)+1
                k_global_frame = modulo(k+srcpos(3,ns)-1,mesh)+1
              else
                i_global_frame = i + srcpos(1,ns)
                j_global_frame = j + srcpos(2,ns)
                k_global_frame = k + srcpos(3,ns)
              endif

              ndens_atom_cell = number_density_array(i_global_frame,j_global_frame,k_global_frame)
              xHI = xHI_array(i_global_frame,j_global_frame,k_global_frame)
              xHII = xHII_array(i_global_frame,j_global_frame,k_global_frame)
              xHeI = xHeI_array(i_global_frame,j_global_frame,k_global_frame)
              xHeII = xHeII_array(i_global_frame,j_global_frame,k_global_frame)
              temperature = temperature_array(i_global_frame,j_global_frame,k_global_frame)
              dr = cell_size
              colum_in_HI = 0
              colum_in_HeI = 0
              colum_in_HeII = 0

              !distance measurement
              distance_from_source = sqrt(real(i)*real(i)+real(j)*real(j)+real(k)*real(k))*dr
              !distance traversed in cell
              if ( (abs(k).ge.abs(i)) .and. (abs(k).ge.abs(j)) ) then
                distance_traverse = sqrt(1.0+(real(i)*real(i)+real(j)*real(j))/(real(k)*real(k)))*dr
              else if ( (abs(i).ge.abs(j)) .and. (abs(i).ge.abs(k)) ) then
                distance_traverse = sqrt(1.0+(real(j)*real(j)+real(k)*real(k))/(real(i)*real(i)))*dr
              else
                distance_traverse = sqrt(1.0+(real(k)*real(k)+real(i)*real(i))/(real(j)*real(j)))*dr
              endif
              vol_shell = 4.0*pi*distance_from_source*distance_from_source*distance_traverse

              colum_out_HI = colum_in_HI+coldens(distance_traverse,xHI,ndens_atom_cell,1.0_dp-abu_he)
              colum_out_HeI = colum_in_HeI+coldens(distance_traverse,xHeI,ndens_atom_cell,abu_he)
              colum_out_HeII = colum_in_HeII+coldens(distance_traverse,xHeII,ndens_atom_cell,abu_he)

              call photoion_rates (phi,colum_in_HI,colum_out_HI,colum_in_HeI,colum_out_HeI, &
                                   colum_in_HeII,colum_out_HeII,vol_shell,ns,xHII)

              phi%photo_cell_HI = phi%photo_cell_HI/(xHI*ndens_atom_cell*(1.0_dp-abu_he))
              phi%photo_cell_HeI = phi%photo_cell_HeI/(xHeI*ndens_atom_cell*abu_he)
              phi%photo_cell_HeII = phi%photo_cell_HeII/(xHeII*ndens_atom_cell*abu_he)

              global_photoionization_HI_array(i_global_frame,j_global_frame,k_global_frame) = &
                        global_photoionization_HI_array(i_global_frame,j_global_frame,k_global_frame) + phi%photo_cell_HI
              global_photoionization_HeI_array(i_global_frame,j_global_frame,k_global_frame) = &
                        global_photoionization_HeI_array(i_global_frame,j_global_frame,k_global_frame) + phi%photo_cell_HeI
              global_photoionization_HeII_array(i_global_frame,j_global_frame,k_global_frame) = &
                        global_photoionization_HeII_array(i_global_frame,j_global_frame,k_global_frame) + phi%photo_cell_HeII

            endif 
          enddo
        enddo
      enddo

!$OMP END DO
!$OMP END PARALLEL

    endif

  end subroutine AT_get_photoionization_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_photoionization
