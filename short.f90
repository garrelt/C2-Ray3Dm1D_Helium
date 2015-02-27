module short

  use precision, only: dp
  use mathconstant, only: four_over_three_pi
  use parameter, only: mesh
  use array, only: coldens_in_HI
  use array, only: coldens_out_HI
  use array, only: coldens_in_HeI
  use array, only: coldens_out_HeI
  use array, only: coldens_in_HeII
  use array, only: coldens_out_HeII
  use array, only: number_density_array
  use array, only: ring_vol
  use parameter, only : cell_size
  use abundances, only: abu_he
  use doric_module, only: coldens

  implicit none

contains

  ! ===========================================================================

    !subroutine short_characteristic (ray_traced_pos,source_pos,global_pos,cdensi,cdensihe0,cdensihe1,path)
    subroutine short_characteristic (ray_traced_pos,source_pos,global_pos,xHI,xHeI,xHeII)
  
    ! Author: Garrelt Mellema
    
    ! Date: 21-Mar-2006 (06-Aug-2004)
    
    ! History:
    ! Original routine written by Alex Raga, Garrelt Mellema, Jane Arthur
    ! and Wolfgang Steffen in 1999.
    ! This version: Modified for use with a grid based approach.
    ! Better handling of the diagonals.
    ! Fortran90
    
    ! does the interpolation to find the column density at pos
    ! as seen from the source point source_pos. the interpolation
    ! depends on the orientation of the ray. The ray crosses either
    ! a z-plane, a y-plane or an x-plane.
    
    implicit none

    integer,dimension(1:3),intent(in) :: ray_traced_pos 
    integer,dimension(1:3),intent(in) :: source_pos 
    integer, dimension(1:3), intent(in) :: global_pos
    real(kind=dp) :: cdensi !< column density to cell
    real(kind=dp) :: cdensihe0 !< column density to cell
    real(kind=dp) :: cdensihe1 !< column density to cell
    !real(kind=dp),intent(out) :: cdensihe1 !< column density to cell
    !real(kind=dp),intent(out) :: path !< path length over cell
    real(kind=dp),intent(in) :: xHI
    real(kind=dp),intent(in) :: xHeI
    real(kind=dp),intent(in) :: xHeII


    real(kind=dp),parameter :: sqrt3=sqrt(3.0)
    real(kind=dp),parameter :: sqrt2=sqrt(2.0)

    integer :: i,j,k,i0,j0,k0
    real(kind=dp) :: number_density_atom
    real(kind=dp) :: path
    integer :: idel,jdel,kdel
    integer :: idela,jdela,kdela
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgni,sgnj,sgnk
    real(kind=dp) :: alam,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1,c2,c3,c4,c1He0,c2He0,c3He0,c4He0,c1He1,c2He1,c3He1,c4He1
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1,w2,w3,w4,w1He0,w2He0,w3He0,w4He0,w1He1,w2He1,w3He1,w4He1
    real(kind=dp) :: di,dj,dk
    real(kind=dp) :: in_x, in_y, in_z
    real(kind=dp) :: out_x, out_y, out_z
    real(kind=dp) :: in_distance, out_distance


    ! map to local variables (should be pointers ;)
    i=ray_traced_pos(1)
    j=ray_traced_pos(2)
    k=ray_traced_pos(3)
    i0=source_pos(1)
    j0=source_pos(2)
    k0=source_pos(3)
 
    ! calculate the distance between the source point (i0,j0,k0) and 
    ! the destination point (i,j,k)
    idel=i-i0
    jdel=j-j0
    kdel=k-k0
    idela=abs(idel)
    jdela=abs(jdel)
    kdela=abs(kdel)
     
if (kdela.ge.idela .and. kdela.ge.jdela) then
in_x = (idela)*(kdela-0.5)/(kdela)
in_y = (jdela)*(kdela-0.5)/(kdela)
in_z = kdela-0.5	
out_x = (idela)*(kdela+0.5)/(kdela)
out_y = (jdela)*(kdela+0.5)/(kdela)				
out_z = kdela+0.5
elseif (jdela.ge.idela .and. jdela.ge.kdela) then
in_x = (idela)*(jdela-0.5)/(jdela)
in_y = jdela-0.5	
in_z = (kdela)*(jdela-0.5)/(jdela)				
out_x = (idela)*(jdela+0.5)/(jdela)
out_y = jdela+0.5
out_z = (kdela)*(jdela+0.5)/(jdela)
elseif (idela.ge.jdela .and. idela.ge.kdela) then
in_x = idela-0.5
in_y = (jdela)*(idela-0.5)/(idela)
in_z = (kdela)*(idela-0.5)/(idela)	
out_x = idela+0.5
out_y = (jdela)*(idela+0.5)/(idela)
out_z = (kdela)*(idela+0.5)/(idela)				  			  
endif
in_distance = in_x*in_x + in_y*in_y + in_z*in_z
in_distance = in_distance**0.5
in_distance = in_distance*cell_size
out_distance = out_x*out_x + out_y*out_y + out_z*out_z
out_distance = out_distance**0.5
out_distance = out_distance*cell_size
ring_vol(global_pos(1),global_pos(2),global_pos(3)) = four_over_three_pi*(out_distance*out_distance*out_distance-&
in_distance*in_distance*in_distance)



    ! Find coordinates of points closer to source
    sgni=sign(1,idel)
    sgnj=sign(1,jdel)
    sgnk=sign(1,kdel)
    im=i-sgni
    jm=j-sgnj
    km=k-sgnk
    di=real(idel)
    dj=real(jdel)
    dk=real(kdel)
 
    ! Z plane (bottom and top face) crossing
    ! we find the central (c) point (xc,xy) where the ray crosses 
    ! the z-plane below or above the destination (d) point, find the 
    ! column density there through interpolation, and add the contribution
    ! of the neutral material between the c-point and the destination
    ! point.
    
    if (kdela >= jdela.and.kdela >= idela) then
       
       ! alam is the parameter which expresses distance along the line s to d
       ! add 0.5 to get to the interface of the d cell.
       alam=(real(km-k0)+sgnk*0.5)/dk
              
       xc=alam*di+real(i0) ! x of crossing point on z-plane 
       yc=alam*dj+real(j0) ! y of crossing point on z-plane
       
       dx=2.0*abs(xc-(real(im)+0.5*sgni)) ! distances from c-point to
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj)) ! the corners.
       
       s1=(1.-dx)*(1.-dy)    ! interpolation weights of
       s2=(1.-dy)*dx         ! corner points to c-point
       s3=(1.-dx)*dy
       s4=dx*dy
        
       ip =modulo(i-1, mesh)+1
       imp=modulo(im-1,mesh)+1
       jp =modulo(j-1, mesh)+1
       jmp=modulo(jm-1,mesh)+1
       kmp=modulo(km-1,mesh)+1

       c1=     coldens_out_HI(imp,jmp,kmp)    !# column densities at the
       c2=     coldens_out_HI(ip,jmp,kmp)     !# four corners
       c3=     coldens_out_HI(imp,jp,kmp)
       c4=     coldens_out_HI(ip,jp,kmp)

       c1he0=coldens_out_HeI(imp,jmp,kmp)    !# column densities at the
       c2he0=coldens_out_HeI(ip,jmp,kmp)     !# four corners
       c3he0=coldens_out_HeI(imp,jp,kmp)
       c4he0=coldens_out_HeI(ip,jp,kmp)

       c1he1=coldens_out_HeII(imp,jmp,kmp)    !# column densities at the
       c2he1=coldens_out_HeII(ip,jmp,kmp)     !# four corners
       c3he1=coldens_out_HeII(imp,jp,kmp)
       c4he1=coldens_out_HeII(ip,jp,kmp)
       
       ! extra weights for better fit to analytical solution
       w1=   s1*weightf(c1,0)
       w2=   s2*weightf(c2,0)
       w3=   s3*weightf(c3,0)
       w4=   s4*weightf(c4,0)

       w1he0=s1*weightf(c1he0,1)
       w2he0=s2*weightf(c2he0,1)
       w3he0=s3*weightf(c3he0,1)
       w4he0=s4*weightf(c4he0,1)

       w1he1=s1*weightf(c1he1,2)
       w2he1=s2*weightf(c2he1,2)
       w3he1=s3*weightf(c3he1,2)
       w4he1=s4*weightf(c4he1,2)

       ! column density at the crossing point
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4) 
       cdensihe0=(c1he0*w1he0+c2he0*w2he0+c3he0*w3he0+c4he0*w4he0)/(w1he0+w2he0+w3he0+w4he0) 
       cdensihe1=(c1he1*w1he1+c2he1*w2he1+c3he1*w3he1+c4he1*w4he1)/(w1he1+w2he1+w3he1+w4he1) 
 
       if (kdela .eq. 1 .and. (idela .eq. 1 .or. jdela .eq. 1)) then
          if (idela .eq. 1 .and. jdela .eq. 1) then
             cdensi=sqrt3*cdensi
             cdensihe0=sqrt3*cdensihe0
             cdensihe1=sqrt3*cdensihe1
          else
             cdensi=sqrt2*cdensi
             cdensihe0=sqrt2*cdensihe0
             cdensihe1=sqrt2*cdensihe1
          endif
       endif
       ! if (kdela == 1) then
       ! if ((w3 == 1.0).or.(w2 == 1.0)) cdensi=sqrt(2.0)*cdensi
       ! if (w1 == 1.0) cdensi=sqrt(3.0)*cdensi
       ! write(logf,*) idela,jdela,kdela
       !endif

       ! Path length from c through d to other side cell.
       !dxp=di/dk
       !dyp=dj/dk
       
       !path=sqrt((di*di+dj*dj)/(dk*dk)+1.0) ! pathlength from c to d point  
 
coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) = cdensi
coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) = cdensihe0
coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) = cdensihe1

path = cell_size*sqrt((di*di+dj*dj)/(dk*dk)+1.0)
number_density_atom = number_density_array(global_pos(1),global_pos(2),global_pos(3))

coldens_out_HI(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHI,number_density_atom,(1.0_dp-abu_he))

coldens_out_HeI(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHeI,number_density_atom,abu_he)

coldens_out_HeII(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHeII,number_density_atom,abu_he)


       ! y plane (left and right face) crossing
       ! (similar approach as for the z plane, see comments there)
    elseif (jdela >= idela.and.jdela >= kdela) then
           
       alam=(real(jm-j0)+sgnj*0.5)/dj
       zc=alam*dk+real(k0)
       xc=alam*di+real(i0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dx=2.0*abs(xc-(real(im)+0.5*sgni))
       s1=(1.-dx)*(1.-dz)
       s2=(1.-dz)*dx
       s3=(1.-dx)*dz
       s4=dx*dz
       ip=modulo(i-1,mesh)+1
       imp=modulo(im-1,mesh)+1
       jmp=modulo(jm-1,mesh)+1
       kp=modulo(k-1,mesh)+1
       kmp=modulo(km-1,mesh)+1

       c1=  coldens_out_HI(imp,jmp,kmp)
       c1he0=coldens_out_HeI(imp,jmp,kmp)
       c1he1=coldens_out_HeII(imp,jmp,kmp)

       c2=  coldens_out_HI(ip,jmp,kmp)
       c2he0=coldens_out_HeI(ip,jmp,kmp)
       c2he1=coldens_out_HeII(ip,jmp,kmp)

       c3=  coldens_out_HI(imp,jmp,kp)
       c3he0=coldens_out_HeI(imp,jmp,kp)
       c3he1=coldens_out_HeII(imp,jmp,kp)

       c4=  coldens_out_HI(ip,jmp,kp)
       c4he0=coldens_out_HeI(ip,jmp,kp)
       c4he1=coldens_out_HeII(ip,jmp,kp)

       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1,0)
       w2=s2*weightf(c2,0)
       w3=s3*weightf(c3,0)
       w4=s4*weightf(c4,0)

       w1he0=s1*weightf(c1he0,1)
       w2he0=s2*weightf(c2he0,1)
       w3he0=s3*weightf(c3he0,1)
       w4he0=s4*weightf(c4he0,1)

       w1he1=s1*weightf(c1he1,2)
       w2he1=s2*weightf(c2he1,2)
       w3he1=s3*weightf(c3he1,2)
       w4he1=s4*weightf(c4he1,2)  
     
       cdensi=   (c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4)
       cdensihe0=(c1he0*w1he0+c2he0*w2he0+c3he0*w3he0+c4he0*w4he0)/(w1he0+w2he0+w3he0+w4he0) 
       cdensihe1=(c1he1*w1he1+c2he1*w2he1+c3he1*w3he1+c4he1*w4he1)/(w1he1+w2he1+w3he1+w4he1)        
       ! Take care of diagonals
       if (jdela .eq. 1 .and. (idela .eq. 1 .or. kdela .eq. 1)) then
          if (idela .eq. 1 .and. kdela .eq. 1) then
             !write(logf,*) 'error',i,j,k
             cdensi=   sqrt3*cdensi
             cdensihe0=sqrt3*cdensihe0
             cdensihe1=sqrt3*cdensihe1
          else
             !write(logf,*) 'diagonal',i,j,k
             cdensi=   sqrt2*cdensi
             cdensihe0=sqrt2*cdensihe0
             cdensihe1=sqrt2*cdensihe1
          endif
       endif

       !path=sqrt((di*di+dk*dk)/(dj*dj)+1.0)
       
coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) = cdensi
coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) = cdensihe0
coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) = cdensihe1

path = cell_size*sqrt((di*di+dk*dk)/(dj*dj)+1.0)
number_density_atom = number_density_array(global_pos(1),global_pos(2),global_pos(3))

coldens_out_HI(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHI,number_density_atom,(1.0_dp-abu_he))

coldens_out_HeI(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHeI,number_density_atom,abu_he)

coldens_out_HeII(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHeII,number_density_atom,abu_he)

       ! x plane (front and back face) crossing
       ! (similar approach as with z plane, see comments there)

    elseif(idela >= jdela.and.idela >= kdela) then

       alam=(real(im-i0)+sgni*0.5)/di
       zc=alam*dk+real(k0)
       yc=alam*dj+real(j0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj))
       s1=(1.-dz)*(1.-dy)
       s2=(1.-dz)*dy
       s3=(1.-dy)*dz
       s4=dy*dz
  
       imp=modulo(im-1,mesh)+1
       jp= modulo(j-1,mesh)+1
       jmp=modulo(jm-1,mesh)+1
       kp= modulo(k-1,mesh)+1
       kmp=modulo(km-1,mesh)+1
       c1=  coldens_out_HI(imp,jmp,kmp)
       c2=  coldens_out_HI(imp,jp,kmp)
       c3=  coldens_out_HI(imp,jmp,kp)
       c4=  coldens_out_HI(imp,jp,kp)

       c1he0=coldens_out_HeI(imp,jmp,kmp)
       c2he0=coldens_out_HeI(imp,jp,kmp)
       c3he0=coldens_out_HeI(imp,jmp,kp)
       c4he0=coldens_out_HeI(imp,jp,kp)

       c1he1=coldens_out_HeII(imp,jmp,kmp)
       c2he1=coldens_out_HeII(imp,jp,kmp)
       c3he1=coldens_out_HeII(imp,jmp,kp)
       c4he1=coldens_out_HeII(imp,jp,kp)

       ! extra weights for better fit to analytical solution
       w1   =s1*weightf(c1,0)
       w2   =s2*weightf(c2,0)
       w3   =s3*weightf(c3,0)
       w4   =s4*weightf(c4,0)

       w1he0=s1*weightf(c1he0,1)
       w2he0=s2*weightf(c2he0,1)
       w3he0=s3*weightf(c3he0,1)
       w4he0=s4*weightf(c4he0,1)

       w1he1=s1*weightf(c1he1,2)
       w2he1=s2*weightf(c2he1,2)
       w3he1=s3*weightf(c3he1,2)
       w4he1=s4*weightf(c4he1,2)      
 
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       cdensihe0=(c1he0*w1he0+c2he0*w2he0+c3he0*w3he0+c4he0*w4he0)/(w1he0+w2he0+w3he0+w4he0) 
       cdensihe1=(c1he1*w1he1+c2he1*w2he1+c3he1*w3he1+c4he1*w4he1)/(w1he1+w2he1+w3he1+w4he1)       
       if ( idela .eq. 1 .and. ( jdela .eq. 1 .or. kdela .eq. 1 ) ) then
          if ( jdela .eq. 1 .and. kdela .eq. 1 ) then
             cdensi=   sqrt3*cdensi
             cdensihe0=sqrt3*cdensihe0            
             cdensihe1=sqrt3*cdensihe1
          else
             cdensi   =sqrt2*cdensi
             cdensihe0=sqrt2*cdensihe0
             cdensihe1=sqrt2*cdensihe1
          endif
       endif
        
       !dyp=dj/di
       !dzp=dk/di
       !path=sqrt(1.0+dyp*dyp+dzp*dzp)
       !path=sqrt(1.0+(dj*dj+dk*dk)/(di*di))
       
coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) = cdensi
coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) = cdensihe0
coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) = cdensihe1

path = cell_size*sqrt((dj*dj+dk*dk)/(di*di)+1.0)
number_density_atom = number_density_array(global_pos(1),global_pos(2),global_pos(3))

coldens_out_HI(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHI,number_density_atom,(1.0_dp-abu_he))

coldens_out_HeI(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHeI,number_density_atom,abu_he)

coldens_out_HeII(global_pos(1),global_pos(2),global_pos(3))=coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) + &
coldens(path,xHeII,number_density_atom,abu_he)

    end if
  
  end subroutine short_characteristic


  ! =========================================================================

  !> Weight function for interpolation in short_characteristic
  real(kind=dp) function weightf (cd,id)

    use cgsphotoconstants, only: sigma_HI_at_ion_freq, sigma_HeI_at_ion_freq, &
         sigma_HeII_at_ion_freq
    real(kind=dp):: sig
    real(kind=dp),intent(in) :: cd
    integer,intent(in) :: id
    real(kind=dp),parameter :: minweight=1.0_dp/0.6_dp

    !weightf=1.0
    ! weightf=1.0/max(1.0d0,cd**0.54)
    ! weightf=exp(-min(700.0,cd*0.15*6.3d-18))
    select case (id)
    case(0)
       sig=sigma_HI_at_ion_freq
    case(1)
       sig=sigma_HeI_at_ion_freq
    case(2)
       sig=sigma_HeII_at_ion_freq
    end select

    weightf=1.0/max(0.6_dp,cd*sig)

    ! weightf=1.0/log(max(e_ln,cd))

  end function weightf

end module short
