!! module do 2-8 finite-difference computations for attenuation
module a2d_staggered_elastic_cpml_kernel

 use global
 use datatype
 use math
 use string
 use io

 implicit none

  real,  allocatable :: a_x(:,:), b_x(:,:), k_x(:), a_z(:,:), b_z(:,:), k_z(:)  !! PML damping along X and Z direction, local matrix
  real,  allocatable :: a_x_half(:,:), b_x_half(:,:), k_x_half(:), &  !! for half grid value
                          a_z_half(:,:), b_z_half(:,:), k_z_half(:)

  real,  allocatable ::   memory_du_dx(:,:), memory_dv_dz(:,:), &
                          memory_du_dz(:,:), memory_dv_dx(:,:), &
                          memory_dtxx_dx(:,:), memory_dtzz_dz(:,:), &
                          memory_dtxz_dx(:,:), memory_dtxz_dz(:,:)

 integer, private                    :: isx,isz,igx,igz,iz,ix,it,ig

contains

subroutine cpml_abc_alloc(nx_pml,nz_pml)

  allocate(a_x(nz_pml,nx_pml), b_x(nz_pml,nx_pml), k_x(nx_pml),&
           a_z(nz_pml,nx_pml), b_z(nz_pml,nx_pml), k_z(nz_pml))

  allocate(a_x_half(nz_pml,nx_pml), b_x_half(nz_pml,nx_pml), k_x_half(nx_pml),&
           a_z_half(nz_pml,nx_pml), b_z_half(nz_pml,nx_pml), k_z_half(nz_pml))

end subroutine cpml_abc_alloc

subroutine cpml_abc_coef(nx_pml,nz_pml,npml,dx,dz,deltat,f0,c0)

  implicit none

  integer::nx_pml,nz_pml,npml

  real::c0(nz_pml,nx_pml)

  real::thickness_PML_x,thickness_PML_z,xoriginleft,xoriginright,&
             zoriginbottom,zorigintop,zero

  real::Rcoef,xval,zval,abscissa_in_PML,&
             abscissa_normalized,alpha_max_pml,quasi_cp_max

  real::deltat,dx,dz,pi
  
  real::d_x(1:nx_pml),d_z(1:nz_pml),alpha_x(1:nx_pml),alpha_z(1:nz_pml)
  real::d_x_half(1:nx_pml),d_z_half(1:nz_pml),alpha_x_half(1:nx_pml),alpha_z_half(1:nz_pml)

  real::d0_x(1:nz_pml,1:nx_pml),d0_z(1:nz_pml,1:nx_pml)

  real::npower,k_max_pml,f0

  integer::i,j,nxpmls,nzpmls

  nxpmls=npml
  nzpmls=npml

  a_x(:,:)=0.0
  b_x(:,:)=0.0
  k_x(:)=1.0
  a_z(:,:)=0.0
  b_z(:,:)=0.0
  k_z(:)=1.0

  a_x_half(:,:)=0.0
  b_x_half(:,:)=0.0
  k_x_half(:)=1.0
  a_z_half(:,:)=0.0
  b_z_half(:,:)=0.0
  k_z_half(:)=1.0

  d_x(:)=0.0
  d_z(:)=0.0
  alpha_x(:)=0.0
  alpha_z(:)=0.0
  d_x_half(:)=0.0
  d_z_half(:)=0.0
  alpha_x_half(:)=0.0
  alpha_z_half(:)=0.0  

  zero=0.0

  pi = 4.0 * atan(1.0)

  npower=2.0
  k_max_pml=1.0

  quasi_cp_max=c0(1,1)

  thickness_PML_x = nxpmls * DX
  thickness_PML_z = nzpmls * Dz

   Rcoef = 0.0001

   maxvp=maxval(winvp)

   if (maxvp > 10000.0) stop

   do ix=1,nx_pml
      do iz=1,nz_pml
         d0_x(iz,ix) = - (NPOWER + 1) * c0(iz,ix) * log(Rcoef) / (2.0 * thickness_PML_x)
         d0_z(iz,ix) = - (NPOWER + 1) * c0(iz,ix) * log(Rcoef) / (2.0 * thickness_PML_z)
      enddo
   enddo

   xoriginleft = thickness_PML_x
   xoriginright = (nx_pml-1)*dx - thickness_PML_x

   do j=1, nz_pml
        do i=1, nx_pml
        !left
             xval = dlen * dble(i-1)
             abscissa_in_PML = xoriginleft - xval

             if(abscissa_in_PML >= zero) then
               abscissa_normalized = abscissa_in_PML / thickness_PML_x
               d_x(i) = d0_x(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
               k_x(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
               alpha_x(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)
             endif

             abscissa_in_PML = xoriginleft - (xval+dlen/2.e0)
             if(abscissa_in_PML >= zero) then
               abscissa_normalized = abscissa_in_PML / thickness_PML_x
               d_x_half(i) = d0_x(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
               k_x_half(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
               alpha_x_half(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)
             endif

         ! right
! define damping profile at the grid points
             abscissa_in_PML = xval - xoriginright
             if(abscissa_in_PML >= zero) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x(i) = d0_x(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                K_x(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                alpha_x(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)
             endif

             abscissa_in_PML = xval + dlen/2.e0 - xoriginright
             if(abscissa_in_PML >= zero) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x_half(i) = d0_x(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                K_x_half(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                alpha_x_half(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)
             endif

         !! just in case

             if(alpha_x(i) < ZERO) alpha_x(i) = ZERO
             if(alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

             b_x(j,i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
             b_x_half(j,i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

! this to avoid division by zero outside the PML
             if(abs(d_x(i)) > 1.d-6) a_x(j,i) = d_x(i) * (b_x(j,i) - 1.0) / (K_x(i) * (d_x(i) &
                                           + K_x(i) * alpha_x(i)))
             if(abs(d_x_half(i)) > 1.d-6) a_x_half(j,i) = d_x_half(i) * &
                      (b_x_half(j,i) - 1.0) / (K_x_half(i) * (d_x_half(i) &
                      + K_x_half(i) * alpha_x_half(i)))

        enddo
   enddo

! origin of the PML layer (position of right edge minus thickness, in meters)
   zoriginbottom = thickness_PML_z
   zorigintop = (nz_pml-1)*dz - thickness_PML_z

   do j = 1,nz_pml
      do i=1,nx_pml

! abscissa of current grid point along the damping profile
             zval = dlen * real(j-1)

!---------- top edge
! define damping profile at the grid points
             abscissa_in_PML = zoriginbottom - zval
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z(j) = d0_z(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)

             endif

             abscissa_in_PML = zoriginbottom - (zval+dlen/2.e0)
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z_half(j) = d0_z(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z_half(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z_half(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)

             endif

!---------- bottom edge
! define damping profile at the grid points
             abscissa_in_PML = zval - zorigintop
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z(j) = d0_z(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)

             endif

             abscissa_in_PML = zval+dlen/2.e0 - zorigintop
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z_half(j) = d0_z(j,i) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z_half(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z_half(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)

             endif
!           endif

             b_z(j,i) = exp(- (d_z(j) / K_z(j) + alpha_z(j)) * DELTAT)
             b_z_half(j,i) = exp(- (d_z_half(j) / K_z_half(j) &
                           + alpha_z_half(j)) * DELTAT)

! this to avoid division by zero outside the PML
             if(abs(d_z(j)) > 1.d-6) a_z(j,i) = d_z(j) * (b_z(j,i) - 1.0) &
                                     / (K_z(j) * (d_z(j) + K_z(j) * alpha_z(j)))
             if(abs(d_z_half(j)) > 1.d-6) a_z_half(j,i) = d_z_half(j) * &
                                   (b_z_half(j,i) - 1.0) / (K_z_half(j) &
                               * (d_z_half(j) + K_z_half(j) * alpha_z_half(j)))

      enddo
    enddo

!bottom and up
 endif

end subroutine cpml_abc_coef 

subroutine iso_els_step_uw_cpml(nx_pml,nz_pml,npml,is,par,dtdx,den2,damp,u,w,xx,zz,xz,&
                                memory_dxx_dx,memory_dxz_dz,)

  implicit none

  type(param),       intent(in)  :: par 
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: u,w,xx,zz,xz
  real,dimension(nz_pml,nx_pml) :: den2
  real,              intent(in) :: damp(:,:)

  !! Update the particle velocity
  !! -------------------------
  !$omp parallel do private(iz,ix,alpha,value_dxx_dx,value_dxz_dz)
  do ix=2,nx_pml-2
  do iz=2,nz_pml-2

     alpha=dtdx/den2(iz,ix)
!     kappa1=damp(iz,ix)*par%dt

     value_dxx_dx = c1_elastic_2th*(xx(iz,ix+1)-xx(iz,ix))
     value_dxz_dz = c1_elastic_2th*(xz(iz,ix)-xz(iz-1,ix))

     memory_dxx_dx(iz,ix) = b_x_half(iz,ix) * memory_dxx_dx(iz,ix) &
                          + a_x_half(iz,ix) * value_dxx_dx

     memory_dxz_dz(iz,ix) = b_z_half(iz,ix) * memory_dxz_dz(iz,ix) &
                          + a_z_half(iz,ix) * value_dxz_dz

     value_dxx_dx = value_dxx_dx / k_x_half(ix) &
                  + memory_dxx_dx(iz,ix)

     value_dxz_dz = value_dxz_dz / k_z_half(iz) &
                  + memory_dxz_dz(iz,ix)

     u(iz,ix)=u(iz,ix)+alpha*(c1_elastic_2th*(xx(iz,ix+1)-xx(iz,ix)+xz(iz,ix)-xz(iz-1,ix)))! +   &
                                           !c2_elastic_8th*(xx(iz,ix+2)-xx(iz,ix-1)+xz(iz+1,ix)-xz(iz-2,ix)) + &
                                           !c3_elastic_8th*(xx(iz,ix+3)-xx(iz,ix-2)+xz(iz+2,ix)-xz(iz-3,ix)) + &
                                           !c4_elastic_8th*(xx(iz,ix+4)-xx(iz,ix-3)+xz(iz+3,ix)-xz(iz-4,ix)))
  enddo
  enddo
  !$omp end parallel do

!  if(fs==1)then
!     do ix=1,nx_pml
!        zz(pad_top,ix)=0.0
!     enddo
!  endif

  !$omp parallel do private(iz,ix,alpha,value_dxz_dx,value_dzz_dz)
  do ix=2,nx_pml-2
  do iz=2,nz_pml-2
     alpha=dtdx/den2(iz,ix)
!     kappa1=damp(iz,ix)*par%dt

     value_dxz_dx = c1_elastic_2th*(xz(iz,ix)-xz(iz,ix-1))
     value_dzz_dz = c1_elastic_2th*(zz(iz+1,ix)-zz(iz,ix))

     memory_dxz_dx(iz,ix) = b_x(iz,ix) * memory_dxz_dx(iz,ix) &
                          + a_x(iz,ix) * value_dxz_dx

     memory_dzz_dz(iz,ix) = b_z(iz,ix) * memory_dzz_dz(iz,ix) &
                          + a_z(iz,ix) * value_dzz_dz

     w(iz,ix)=w(iz,ix)+alpha*(c1_elastic_2th*(zz(iz+1,ix)-zz(iz,ix)+xz(iz,ix)-xz(iz,ix-1)))! +   &
                                           !c2_elastic_8th*(zz(iz+2,ix)-zz(iz-1,ix)+xz(iz,ix+1)-xz(iz,ix-2)) + &
                                           !c3_elastic_8th*(zz(iz+3,ix)-zz(iz-2,ix)+xz(iz,ix+2)-xz(iz,ix-3)) + &
                                           !c4_elastic_8th*(zz(iz+4,ix)-zz(iz-3,ix)+xz(iz,ix+3)-xz(iz,ix-4)))

  enddo
  enddo
  !$omp end parallel do        

end subroutine iso_els_step_uw_cpml


subroutine iso_els_step_sigma_cpml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)

  implicit none

  type(param),       intent(in)  :: par
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: u,w,xx,zz,xz
  real,dimension(nz_pml,nx_pml) :: lamda, mu
  real,              intent(in) :: damp(:,:)

  !! Update stress
  !! --------------------------
  !$omp parallel do private(iz,ix,alpha,alpha2,value_dvx_dx,value_dvz_dz)
  do ix=2,nx_pml-2
  do iz=2,nz_pml-2
     alpha=dtdx*(lamda(iz,ix)+2.0*mu(iz,ix))
     alpha2=dtdx*lamda(iz,ix)
!     kappa1=damp(iz,ix)*par%dt

     value_dvx_dx = c1_elastic_2th*(u(iz,ix)-u(iz,ix-1))     
     value_dvz_dz = c1_elastic_2th*(w(iz,ix)-w(iz-1,ix))

     memory_dvx_dx(iz,ix) = b_x(iz,ix) * memory_dvx_dx(iz,ix) + a_x(iz,ix) * value_dvx_dx
     memory_dvz_dz(iz,ix) = b_z(iz,ix) * memory_dvz_dz(iz,ix) + a_z(iz,ix) * value_dvz_dz

     value_dvx_dx = value_dvx_dx / K_x(ix) + memory_dvx_dx(iz,ix)
     value_dvz_dz = value_dvz_dz / K_z(iz) + memory_dvz_dz(iz,ix)

     xx(iz,ix)=xx(iz,ix)+alpha*(c1_elastic_2th*(u(iz,ix)-u(iz,ix-1))) +    &
                                            ! c2_elastic_8th*(u(iz,ix+1)-u(iz,ix-2))+   &
                                            ! c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3))+   &
                                            ! c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4))) + &
                                    + alpha2*(c1_elastic_2th*(w(iz,ix)-w(iz-1,ix)))! +    &
                                            ! c2_elastic_8th*(w(iz+1,ix)-w(iz-2,ix))+   &
                                            ! c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+   &
                                            ! c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix)))    

     zz(iz,ix)=zz(iz,ix)+alpha2*(c1_elastic_2th*(u(iz,ix)-u(iz,ix-1))) +   &
                                             !c2_elastic_8th*(u(iz,ix+1)-u(iz,ix-2)) + &
                                             !c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3)) + &
                                             !c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4)))+ &
                                    + alpha*(c1_elastic_2th*(w(iz,ix)-w(iz-1,ix)))! +   &
                                            ! c2_elastic_8th*(w(iz+1,ix)-w(iz-2,ix))+  &
                                            ! c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+  &
                                            ! c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix)))

  enddo
  enddo
  !$omp end parallel do


  !$omp parallel do private(iz,ix,alpha,alpha2,value_dvz_dx,value_dvx_dz)
  do ix=2,nx_pml-2
  do iz=2,nz_pml-2
     alpha=dtdx*mu(iz,ix)
!     kappa1=damp(iz,ix)*par%dt
 
     value_dvz_dx = c1_elastic_2th*(w(iz,ix+1)-w(iz,ix))
     value_dvx_dz = c1_elastic_2th*(u(iz+1,ix)-u(iz,ix))

     memory_dvz_dx(iz,ix) = b_x_half(iz,ix) * memory_dvz_dx(iz,ix) &
                        + a_x_half(iz,ix) * value_dvz_dx

     memory_dvx_dz(iz,ix) = b_z_half(iz,ix) * memory_dvx_dz(iz,ix) &
                                 + a_z_half(iz,ix) * value_dvx_dz

     xz(iz,ix)=xz(iz,ix)+alpha*(c1_elastic_2th*(u(iz+1,ix)-u(iz,ix)+w(iz,ix+1)-w(iz,ix)))!+  &
                                             !c2_elastic_8th*(u(iz+2,ix)-u(iz-1,ix)+w(iz,ix+2)-w(iz,ix-1))+  &
                                             !c3_elastic_8th*(u(iz+3,ix)-u(iz-2,ix)+w(iz,ix+3)-w(iz,ix-2))+  &
                                             !c4_elastic_8th*(u(iz+4,ix)-u(iz-3,ix)+w(iz,ix+4)-w(iz,ix-3)))
  enddo
  enddo
  !$omp end parallel do        

end subroutine iso_els_step_sigma_cpml

!subroutine iso_els_step_uw_cpml()
!
!end subroutine iso_els_step_uw_cpml


!subroutine iso_els_step_sigma_cpml()
!
!end subroutine iso_els_step_sigma_cpml

subroutine iso_els_step_uw_adj_lagran_cpml(nx_pml,nz_pml,npml,is,par,dtdx,den2,lamda,mu,damp,ub,wb,xxb,zzb,xzb)

  implicit none

  type(param),       intent(in)  :: par
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: ub,wb,xxb,zzb,xzb
  real,dimension(nz_pml,nx_pml) :: den2,lamda,mu
  real,              intent(in) :: damp(:,:)

 !! Update the particle velocity
 !! -------------------------
 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=2,nx_pml-2
 do iz=2,nz_pml-2
    alpha=dtdx/den2(iz,ix)
    kappa1=damp(iz,ix)*par%dt
    ub(iz,ix)=(1.0-kappa1)*ub(iz,ix)+alpha*(c1_elastic_2th*((lamda(iz,ix+1)+2.0*mu(iz,ix+1))*xxb(iz,ix+1)-(lamda(iz,ix)+2.0*mu(iz,ix))*xxb(iz,ix)+lamda(iz,ix+1)*zzb(iz,ix+1)-lamda(iz,ix)*zzb(iz,ix)+mu(iz,ix)*xzb(iz,ix)-mu(iz-1,ix)*xzb(iz-1,ix)))! +     &
                                            !c2_elastic_8th*((lamda(iz,ix+2)+2.0*mu(iz,ix+2))*xxb(iz,ix+2)-(lamda(iz,ix-1)+2.0*mu(iz,ix-1))*xxb(iz,ix-1)+lamda(iz,ix+2)*zzb(iz,ix+2)-lamda(iz,ix-1)*zzb(iz,ix-1)+mu(iz+1,ix)*xzb(iz+1,ix)-mu(iz-2,ix)*xzb(iz-2,ix)) + &
                                            !c3_elastic_8th*((lamda(iz,ix+3)+2.0*mu(iz,ix+3))*xxb(iz,ix+3)-(lamda(iz,ix-2)+2.0*mu(iz,ix-2))*xxb(iz,ix-2)+lamda(iz,ix+3)*zzb(iz,ix+3)-lamda(iz,ix-2)*zzb(iz,ix-2)+mu(iz+2,ix)*xzb(iz+2,ix)-mu(iz-3,ix)*xzb(iz-3,ix)) + &
                                            !c4_elastic_8th*((lamda(iz,ix+4)+2.0*mu(iz,ix+4))*xxb(iz,ix+4)-(lamda(iz,ix-3)+2.0*mu(iz,ix-3))*xxb(iz,ix-3)+lamda(iz,ix+4)*zzb(iz,ix+4)-lamda(iz,ix-3)*zzb(iz,ix-3)+mu(iz+3,ix)*xzb(iz+3,ix)-mu(iz-4,ix)*xzb(iz-4,ix)))
 enddo
 enddo
 !$omp end parallel do

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=2,nx_pml-2
 do iz=2,nz_pml-2
    alpha=dtdx/den2(iz,ix)
    kappa1=damp(iz,ix)*par%dt
    wb(iz,ix)=(1.0-kappa1)*wb(iz,ix)+alpha*(c1_elastic_2th*((lamda(iz+1,ix)+2.0*mu(iz+1,ix))*zzb(iz+1,ix)-(lamda(iz,ix)+2.0*mu(iz,ix))*zzb(iz,ix)+lamda(iz+1,ix)*xxb(iz+1,ix)-lamda(iz,ix)*xxb(iz,ix)+mu(iz,ix)*xzb(iz,ix)-mu(iz,ix-1)*xzb(iz,ix-1)))! +     &
                                           ! c2_elastic_8th*((lamda(iz+2,ix)+2.0*mu(iz+2,ix))*zzb(iz+2,ix)-(lamda(iz-1,ix)+2.0*mu(iz-1,ix))*zzb(iz-1,ix)+lamda(iz+2,ix)*xxb(iz+2,ix)-lamda(iz-1,ix)*xxb(iz-1,ix)+mu(iz,ix+1)*xzb(iz,ix+1)-mu(iz,ix-2)*xzb(iz,ix-2)) + &
                                           ! c3_elastic_8th*((lamda(iz+3,ix)+2.0*mu(iz+3,ix))*zzb(iz+3,ix)-(lamda(iz-2,ix)+2.0*mu(iz-2,ix))*zzb(iz-2,ix)+lamda(iz+3,ix)*xxb(iz+3,ix)-lamda(iz-2,ix)*xxb(iz-2,ix)+mu(iz,ix+2)*xzb(iz,ix+2)-mu(iz,ix-3)*xzb(iz,ix-3)) + &
                                           ! c4_elastic_8th*((lamda(iz+4,ix)+2.0*mu(iz+4,ix))*zzb(iz+4,ix)-(lamda(iz-3,ix)+2.0*mu(iz-3,ix))*zzb(iz-3,ix)+lamda(iz+4,ix)*xxb(iz+4,ix)-lamda(iz-3,ix)*xxb(iz-3,ix)+mu(iz,ix+3)*xzb(iz,ix+3)-mu(iz,ix-4)*xzb(iz,ix-4)))
 enddo
 enddo
 !$omp end parallel do        

end subroutine iso_els_step_uw_adj_lagran_cpml

subroutine iso_els_step_sigma_adj_lagran_cpml(nx_pml,nz_pml,npml,is,par,dtdx,damp,ub,wb,xxb,zzb,xzb)

  implicit none

  type(param),       intent(in)  :: par
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: ub,wb,xxb,zzb,xzb
  real,              intent(in) :: damp(:,:)

 !! Update stress
 !! --------------------------
 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=2,nx_pml-2
 do iz=2,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    xxb(iz,ix)=(1.0-kappa1)*xxb(iz,ix)+alpha*(c1_elastic_2th*(ub(iz,ix)-ub(iz,ix-1)))! +  &
                                              !c2_elastic_8th*(ub(iz,ix+1)-ub(iz,ix-2))+ &
                                              !c3_elastic_8th*(ub(iz,ix+2)-ub(iz,ix-3))+ &
                                              !c4_elastic_8th*(ub(iz,ix+3)-ub(iz,ix-4)))
 enddo
 enddo
 !$omp end parallel do

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=2,nx_pml-2
 do iz=2,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    zzb(iz,ix)=(1.0-kappa1)*zzb(iz,ix)+alpha*(c1_elastic_2th*(wb(iz,ix)-wb(iz-1,ix)))! +   &
                                              !c2_elastic_8th*(wb(iz+1,ix)-wb(iz-2,ix)) + &
                                              !c3_elastic_8th*(wb(iz+2,ix)-wb(iz-2,ix)) + &
                                              !c4_elastic_8th*(wb(iz+3,ix)-wb(iz-3,ix)))
 enddo
 enddo
 !$omp end parallel do

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=2,nx_pml-2
 do iz=2,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    xzb(iz,ix)=(1.0-kappa1)*xzb(iz,ix)+alpha*(c1_elastic_2th*(ub(iz+1,ix)-ub(iz,ix)+wb(iz,ix+1)-wb(iz,ix)))! +   &
                                              !c2_elastic_8th*(ub(iz+2,ix)-ub(iz-1,ix)+wb(iz,ix+2)-wb(iz,ix-1)) + &
                                              !c3_elastic_8th*(ub(iz+3,ix)-ub(iz-2,ix)+wb(iz,ix+3)-wb(iz,ix-2)) + &
                                              !c4_elastic_8th*(ub(iz+4,ix)-ub(iz-3,ix)+wb(iz,ix+4)-wb(iz,ix-3)))
 enddo
 enddo
 !$omp end parallel do

end subroutine iso_els_step_sigma_adj_lagran_cpml


end module a2d_staggered_elastic_cpml_kernel
