module evolve_source

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,timefile,iterdump, results_dir, dump_dir
  use abundances, only: abu_he
  use c2ray_parameters, only: subboxsize, max_subbox
  use parameter, only: mesh
  use parameter, only: number_of_source
  use sourceprops, only: srcpos
  use array, only: coldens_in_HI,coldens_out_HI
  use array, only: coldens_in_HeI,coldens_out_HeI
  use array, only: coldens_in_HeII,coldens_out_HeII
  use array, only: ring_vol
  use evolve_point, only: photo_photo
  use column_density_0D, only: column_density_0D_point
  use column_density_1D, only: column_density_1D_axis
  use column_density_2D, only: column_density_2D_plane
  use column_density_3D, only: column_density_3D_octant

  implicit none

  ! mesh positions of end points for RT


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine do_source_photo(source_index)

    integer, intent(in) :: source_index !< number of the source being done

    integer :: naxis,nplane,nquadrant
    integer :: k
    integer :: nnt
    integer :: logf1
    real(kind=dp) :: total_source_flux
    integer,dimension(1:3) :: rtpos
    integer :: tn



    !!!! 0D !!!!
    rtpos(:)=srcpos(:,source_index)
    call photo_photo(rtpos,source_index)

    !!!! 1D !!!!
    do naxis=1,6
      call evolve1D_photo(source_index,naxis)
    enddo

    !!!! 2D !!!!
    do nplane=1,12
      call evolve2D_photo(source_index,nplane)
    end do

    !!!! 3D !!!!
    do nquadrant=1,8
      call evolve3D_photo(source_index,nquadrant)
    end do

  end subroutine do_source_photo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Ray tracing for the axes going through the source point
  ! should be called after having done the source point
  subroutine evolve1D_photo(source_index,naxis)

    integer,intent(in) :: source_index           ! current source
    integer,intent(in) :: naxis        ! axis to do

  integer,dimension(1:3) :: last_l !< mesh position of left end point for RT
  integer,dimension(1:3) :: last_r !< mesh position of right end point for RT

    integer :: i,j,k
    integer,dimension(1:3) :: rtpos ! mesh position

    last_r(:)=srcpos(:,source_index)+min(max_subbox,mesh/2-1+mod(mesh,2))
    last_l(:)=srcpos(:,source_index)-min(max_subbox,mesh/2)

    select case (naxis)
    case(1)
       ! sweep in +i direction
       rtpos(2:3)=srcpos(2:3,source_index)
       do i=srcpos(1,source_index)+1,last_r(1)
          rtpos(1)=i
          call photo_photo(rtpos,source_index) !# `positive' i
       enddo
    case(2)
       ! sweep in -i direction
       rtpos(2:3)=srcpos(2:3,source_index)
       do i=srcpos(1,source_index)-1,last_l(1),-1
          rtpos(1)=i
          call photo_photo(rtpos,source_index) !# `negative' i
       end do
    case(3)
       ! sweep in +j direction
       rtpos(1)=srcpos(1,source_index)
       rtpos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)+1,last_r(2)
          rtpos(2)=j
          call photo_photo(rtpos,source_index) !# `positive' j
       end do
    case(4)
       ! sweep in -j direction
       rtpos(1)=srcpos(1,source_index)
       rtpos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)-1,last_l(2),-1
          rtpos(2)=j
          call photo_photo(rtpos,source_index) !# `negative' j
       end do
    case(5)
       ! sweep in +k direction
       rtpos(1:2)=srcpos(1:2,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          call photo_photo(rtpos,source_index) !# `positive' k
       end do
    case(6)
       ! sweep in -k direction
       rtpos(1:2)=srcpos(1:2,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          call photo_photo(rtpos,source_index) !# `negative' k
       end do
    end select
    
  end subroutine evolve1D_photo

  ! ===========================================================================

  !> Ray tracing for planes containing the source point
  !! should be called after evolve1D_column
  subroutine evolve2D_photo(source_index,nplane)

    integer,intent(in) :: source_index           ! current source
    integer,intent(in) :: nplane        ! plane to do

  integer,dimension(1:3) :: last_l !< mesh position of left end point for RT
  integer,dimension(1:3) :: last_r !< mesh position of right end point for RT

    integer :: i,j,k
    integer,dimension(1:3) :: rtpos ! mesh position

    last_r(:)=srcpos(:,source_index)+min(max_subbox,mesh/2-1+mod(mesh,2))
    last_l(:)=srcpos(:,source_index)-min(max_subbox,mesh/2)

    select case (nplane)
    case(1)
       ! sweep in +i,+j direction
       rtpos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,source_index)+1,last_r(1)
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(2)
       ! sweep in +i,-j direction
       rtpos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,source_index)+1,last_r(1)
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(3)
       ! sweep in -i,+j direction
       rtpos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,source_index)-1,last_l(1),-1
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(4)
       ! sweep in -i,-j direction
       rtpos(3)=srcpos(3,source_index)
       do j=srcpos(2,source_index)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,source_index)-1,last_l(1),-1
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(5)
       ! sweep in +i,+k direction
       rtpos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,source_index)+1,last_r(1)
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(6)
       ! sweep in -i,+k direction
       rtpos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,source_index)-1,last_l(1),-1
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(7)
       ! sweep in -i,-k direction
       rtpos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,source_index)-1,last_l(1),-1
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(8)
       ! sweep in +i,-k direction
       rtpos(2)=srcpos(2,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,source_index)+1,last_r(1)
             rtpos(1)=i
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(9) 
       ! sweep in +j,+k direction
       rtpos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             rtpos(2)=j
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(10) 
       ! sweep in -j,+k direction
       rtpos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             rtpos(2)=j
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(11) 
       ! sweep in +j,-k direction
       rtpos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             rtpos(2)=j
             call photo_photo(rtpos,source_index)
          enddo
       enddo
    case(12) 
       ! sweep in -j,-k direction
       rtpos(1)=srcpos(1,source_index)
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             rtpos(2)=j
             call photo_photo(rtpos,source_index)
          enddo
       enddo
       
    end select
    
  end subroutine evolve2D_photo

  ! ===========================================================================

  !> Ray tracing for the 8 octants 
  !! should be called after evolve2D_column
  subroutine evolve3D_photo(source_index,nquadrant)

    ! find column density for a z-plane srcpos(3) by sweeping in x and y
    ! directions
    integer,intent(in) :: source_index           ! current source
    integer,intent(in) :: nquadrant    ! which quadrant to do    

  integer,dimension(1:3) :: last_l !< mesh position of left end point for RT
  integer,dimension(1:3) :: last_r !< mesh position of right end point for RT

    integer :: i,j,k
    integer,dimension(1:3) :: rtpos ! mesh position

    last_r(:)=srcpos(:,source_index)+min(max_subbox,mesh/2-1+mod(mesh,2))
    last_l(:)=srcpos(:,source_index)-min(max_subbox,mesh/2)

    select case (nquadrant)
    case (1)
       ! sweep in +i,+j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,source_index)+1,last_r(1)
                rtpos(1)=i
                call photo_photo(rtpos,source_index)
             end do
          enddo
       enddo
    case (2)
       ! sweep in -i,+j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `negative' i
             end do
          end do
       enddo
    case (3)
       ! sweep in +i,-j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,source_index)+1,last_r(1)
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `negative' i
             end do
          end do
       enddo
    case(4)
       ! sweep in -i,-j,+k direction
       do k=srcpos(3,source_index)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `negative' i
             end do
          end do
       enddo
    case (5)
       ! sweep in +i,+j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,source_index)+1,last_r(1)
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `positive' i
             end do
          enddo
       enddo
    case (6)
       ! sweep in -i,+j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,source_index)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `negative' i
             end do
          end do
       enddo
    case (7)
       ! sweep in +i,-j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,source_index)+1,last_r(1)
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `negative' i
             end do
          end do
       enddo
    case(8)
       ! sweep in -i,-j,-k direction
       do k=srcpos(3,source_index)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,source_index)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,source_index)-1,last_l(1),-1
                rtpos(1)=i
                call photo_photo(rtpos,source_index) !# `negative' i
             end do
          end do
       enddo
    end select

  end subroutine evolve3D_photo

end module evolve_source
