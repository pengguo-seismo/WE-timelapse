
   subroutine free_surface_denise(nx_pml,nz_pml,iz,dt,dx,fd_order,lamda,mu,xx,zz,xz,u,w)

     use global

     implicit none
     integer:: nx_pml,nz_pml
     integer:: ii, ix, iz, lv, fd_order
     
     real:: dt, dx
     real:: tmp_u, tmp_w

     real,dimension(1:nz_pml,1:nx_pml)::lamda,mu,xx,zz,xz,u,w

     ! The free surface is located exactly in y=(ndepth-1/2)*dh
     ! consistent with SOFI and DENISE
     do ix=1,nx_pml
        zz(iz, ix)=0.0
     enddo

     lv = fd_order/2

     !2nd order accuracy
     do ix = 2, nx_pml-1
     
      do ii = 1, lv   
        zz(iz-ii, ix) = -zz(iz+ii,ix)
        xz(iz-ii,ix) = -xz(iz+ii-1,ix)
      enddo 

      tmp_u = 0.0
      tmp_w = 0.0

      if (fd_order == 2) then 

       tmp_u = c1_elastic_2th*(u(iz,ix) -u(iz,ix-1))/dx
       tmp_w = c1_elastic_2th*(w(iz,ix) -w(iz-1,ix))/dx

      elseif (fd_order == 4) then

       tmp_u = c1_elastic_4th*(u(iz,ix) -u(iz,ix-1))/dx + c2_elastic_4th*(u(iz,ix+1) -u(iz,ix-2))/dx
       tmp_w = c1_elastic_4th*(w(iz,ix) -w(iz-1,ix))/dx + c2_elastic_4th*(w(iz+1,ix) -w(iz-2,ix))/dx

      endif 

       xx(iz,ix) = xx(iz,ix) - dt*((tmp_u*(lamda(iz,ix)*lamda(iz,ix)))/(lamda(iz,ix)+2.0*mu(iz,ix))+lamda(iz,ix)*tmp_w)
     enddo

    end subroutine free_surface_denise

   subroutine free_surface_legacy(nx_pml,nz_pml,iz,dt,dx,zz)

    implicit none
    integer:: nx_pml,nz_pml
    integer:: ix, iz     
    real:: dt, dx

    real,dimension(1:nz_pml,1:nx_pml)::zz

     do ix=1,nx_pml
        zz(iz, ix)=0.0
     enddo   

   end subroutine free_surface_legacy

   subroutine free_surface_JIMU_ss(nx_pml,nz_pml,iz,dt,dx,fd_order,zz,xz)

    implicit none
    integer:: nx_pml,nz_pml
    integer:: ii, ix, iz, fd_order, lv
    real:: dt, dx
    real,dimension(1:nz_pml,1:nx_pml)::zz,xz

     do ix=1,nx_pml
        zz(iz, ix)=0.0
     enddo

     lv = fd_order/2

     do ix = 2, nx_pml-1

       do ii = 1, lv
         zz(iz-ii, ix) = -zz(iz+ii,ix)
         xz(iz-ii,ix) = -xz(iz+ii-1,ix)
       enddo

     enddo 

  !   xz(iz,:)=-xz(iz+1,:)

   end subroutine free_surface_JIMU_ss

!   subroutine free_surface_JIMU_uw(nx_pml,nz_pml,iz,dt,dx,lamda,mu,xx,zz,xz) 
!
!    implicit none
!    integer:: nx_pml,nz_pml
!    integer:: ix, iz
!    real:: dt, dx
!
!    f%vz(iz,ix)= f%vz(iz+1,ix) + lamda(1,ix)*(f%vx(1,ix)-f%vx(1,ix+1))*dz_dx/lm%ldap2mu(1,ix)
!
!   end subroutine free_surface_JIMU_uw

