module total_photo_rate

  use precision, only: dp 
  use my_mpi
  use array, only: photoionization_HI_array, &
                   photoionization_HeI_array, &
                   photoionization_HeII_array, &
                   heating_array, &
                   buffer
  use parameter, only: mesh

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine total_photo_rate_from_all_sources()

    implicit none

    integer :: mympierror

#ifdef MPI

     call MPI_ALLREDUCE(photoionization_HI_array, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_NEW, mympierror)        
     photoionization_HI_array(:,:,:) = buffer(:,:,:)

     call MPI_ALLREDUCE(photoionization_HeI_array, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_NEW, mympierror)        
     photoionization_HeI_array(:,:,:) = buffer(:,:,:)

     call MPI_ALLREDUCE(photoionization_HeII_array, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_NEW, mympierror)        
     photoionization_HeII_array(:,:,:) = buffer(:,:,:)

     call MPI_ALLREDUCE(heating_array, buffer, mesh*mesh*mesh,MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_NEW, mympierror)
     heating_array(:,:,:) = buffer(:,:,:)

#endif

!write(*,*) 'sum PHI  ',sum(photoionization_HI_array)
!write(*,*) 'sum PHeI ',sum(photoionization_HeI_array)
!write(*,*) 'sum PHeII',sum(photoionization_HeII_array)
!write(*,*) 'sum Heat ',sum(heating_array)
!write(*,*) 'some PHI  ',photoionization_HI_array(15,15,14)
!write(*,*) 'some PHeI ',photoionization_HeI_array(15,15,14)
!write(*,*) 'some PHeII',photoionization_HeII_array(15,15,14)


   end subroutine total_photo_rate_from_all_sources

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module total_photo_rate
