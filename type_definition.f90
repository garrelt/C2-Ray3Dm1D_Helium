module type_definition

  use precision, only: dp
  use parameter, only: NumFreqBnd

  implicit none

  ! photrates contains all the photo-ionization rates and heating rates
  type photrates    
     real(kind=dp) :: photo_cell_HI          ! HI photoionization rate of the cell    
     real(kind=dp) :: photo_cell_HeI         ! HeI photoionization rate of the cell    
     real(kind=dp) :: photo_cell_HeII        ! HeII photoionization rate of the cell    
     real(kind=dp) :: heat_cell_HI           ! HI heating rate of the cell       
     real(kind=dp) :: heat_cell_HeI          ! HeI heating rate of the cell    
     real(kind=dp) :: heat_cell_HeII         ! HeII heating rate of the cell          
     real(kind=dp) :: photo_in_HI            ! HI photoionization rate incoming to the cell    
     real(kind=dp) :: photo_in_HeI           ! HeI photoionization rate incoming to the cell
     real(kind=dp) :: photo_in_HeII          ! HeII photoionization rate incoming to the cell
     real(kind=dp) :: heat_in_HI             ! HI heating rate incoming to the cell
     real(kind=dp) :: heat_in_HeI            ! HeI heating rate incoming to the cell
     real(kind=dp) :: heat_in_HeII           ! HeII heating rate incoming to the cell 
     real(kind=dp) :: photo_out_HI           ! HI photoionization rate outgoing from the cell
     real(kind=dp) :: photo_out_HeI          ! HeI photoionization rate outgoing from the cell
     real(kind=dp) :: photo_out_HeII         ! HeII photoionization rate outgoing from the cell 
     real(kind=dp) :: heat_out_HI            ! HI heating rate outgoing from the cell
     real(kind=dp) :: heat_out_HeI           ! HeI heating rate outgoing from the cell
     real(kind=dp) :: heat_out_HeII          ! HeII heating rate outgoing from the cell
     real(kind=dp) :: heat                   ! Total heating rate of the cell
     real(kind=dp) :: photo_in               ! Total photoionization rate incoming to the cell
     real(kind=dp) :: photo_out               ! Total photoionization rate incoming to the cell
  end type photrates

  ! tablepos helps to locate correct position of the photoionization and heating tables
  type tablepos
    real(kind=dp), dimension(NumFreqBnd) :: tau            
    real(kind=dp), dimension(NumFreqBnd) :: odpos          
    real(kind=dp), dimension(NumFreqBnd) :: residual       
    integer, dimension(NumFreqBnd)       :: ipos           
    integer, dimension(NumFreqBnd)       :: ipos_p1        
  end type tablepos 

   type ionstates    
     real(kind=dp) :: end_HI 
     real(kind=dp) :: end_HII
     real(kind=dp) :: end_HeI
     real(kind=dp) :: end_HeII
     real(kind=dp) :: end_HeIII
     real(kind=dp) :: avg_HI
     real(kind=dp) :: avg_HII
     real(kind=dp) :: avg_HeI
     real(kind=dp) :: avg_HeII
     real(kind=dp) :: avg_HeIII
     real(kind=dp) :: begin_HI
     real(kind=dp) :: begin_HII
     real(kind=dp) :: begin_HeI
     real(kind=dp) :: begin_HeII
     real(kind=dp) :: begin_HeIII
  end type ionstates

end module type_definition
