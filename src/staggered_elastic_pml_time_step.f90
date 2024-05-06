module staggered_elastic_pml_time_step

!            
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        sigma_xz   |
!       v_z  +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!         sigma_xx    v_x
!         sigma_zz
!

 use global
 use datatype
 use math
 use string
 use io

 implicit none

 integer, private                    :: isx,isz,igx,igz,iz,ix,it,ig

contains

subroutine iso_els_step_uw_pml(nx_pml,nz_pml,npml,is,par,dtdx,den2,damp,u,w,xx,zz,xz)

  implicit none

  type(param),       intent(in)  :: par 
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: u,w,xx,zz,xz
  real,dimension(nz_pml,nx_pml) :: den2
  real,              intent(in) :: damp(:,:)

  integer :: ix1, iz1

  if (par%fd_order==2) then
      ix1 = 2
      iz1 = 2
  elseif (par%fd_order==4) then
      ix1 = 3
      iz1 = 3
  endif

  if (par%free_surface == 1) then
     ix1=3
     iz1=par%npml+1
  endif

  !! Update the particle velocity
  !! -------------------------
  if (par%fd_order==2) then 

  !$omp parallel do private(iz,ix,alpha,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx/den2(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     u(iz,ix)=(1.0-kappa1)*u(iz,ix)+alpha*(c1_elastic_2th*(xx(iz,ix+1)-xx(iz,ix)+xz(iz,ix)-xz(iz-1,ix)))! +   &
                                           !c2_elastic_8th*(xx(iz,ix+2)-xx(iz,ix-1)+xz(iz+1,ix)-xz(iz-2,ix)) + &
                                           !c3_elastic_8th*(xx(iz,ix+3)-xx(iz,ix-2)+xz(iz+2,ix)-xz(iz-3,ix)) + &
                                           !c4_elastic_8th*(xx(iz,ix+4)-xx(iz,ix-3)+xz(iz+3,ix)-xz(iz-4,ix)))
  enddo
  enddo
  !$omp end parallel do

  elseif (par%fd_order==4) then

  !$omp parallel do private(iz,ix,alpha,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx/den2(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     u(iz,ix)=(1.0-kappa1)*u(iz,ix)+alpha*(c1_elastic_4th*(xx(iz,ix+1)-xx(iz,ix)+xz(iz,ix)-xz(iz-1,ix))    &
                                         + c2_elastic_4th*(xx(iz,ix+2)-xx(iz,ix-1)+xz(iz+1,ix)-xz(iz-2,ix))) 
                                           !c3_elastic_8th*(xx(iz,ix+3)-xx(iz,ix-2)+xz(iz+2,ix)-xz(iz-3,ix)) + &
                                           !c4_elastic_8th*(xx(iz,ix+4)-xx(iz,ix-3)+xz(iz+3,ix)-xz(iz-4,ix)))
  enddo
  enddo
  !$omp end parallel do

  endif 

!  if(fs==1)then
!     do ix=1,nx_pml
!        zz(pad_top,ix)=0.0
!     enddo
!  endif

  if (par%fd_order==2) then

  !$omp parallel do private(iz,ix,alpha,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx/den2(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     w(iz,ix)=(1.0-kappa1)*w(iz,ix)+alpha*(c1_elastic_2th*(zz(iz+1,ix)-zz(iz,ix)+xz(iz,ix)-xz(iz,ix-1)))! +   &
                                           !c2_elastic_8th*(zz(iz+2,ix)-zz(iz-1,ix)+xz(iz,ix+1)-xz(iz,ix-2)) + &
                                           !c3_elastic_8th*(zz(iz+3,ix)-zz(iz-2,ix)+xz(iz,ix+2)-xz(iz,ix-3)) + &
                                           !c4_elastic_8th*(zz(iz+4,ix)-zz(iz-3,ix)+xz(iz,ix+3)-xz(iz,ix-4)))

  enddo
  enddo
  !$omp end parallel do        

  elseif (par%fd_order==4) then

  !$omp parallel do private(iz,ix,alpha,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx/den2(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     w(iz,ix)=(1.0-kappa1)*w(iz,ix)+alpha*(c1_elastic_4th*(zz(iz+1,ix)-zz(iz,ix)+xz(iz,ix)-xz(iz,ix-1))   &
                                         + c2_elastic_4th*(zz(iz+2,ix)-zz(iz-1,ix)+xz(iz,ix+1)-xz(iz,ix-2))) 
                                           !c3_elastic_8th*(zz(iz+3,ix)-zz(iz-2,ix)+xz(iz,ix+2)-xz(iz,ix-3)) + &
                                           !c4_elastic_8th*(zz(iz+4,ix)-zz(iz-3,ix)+xz(iz,ix+3)-xz(iz,ix-4)))

  enddo
  enddo
  !$omp end parallel do

  endif

end subroutine iso_els_step_uw_pml


subroutine iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)

  implicit none

  type(param),       intent(in)  :: par
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  integer :: ix1, iz1

  real,dimension(nz_pml,nx_pml) :: u,w,xx,zz,xz
  real,dimension(nz_pml,nx_pml) :: lamda, mu
  real,              intent(in) :: damp(:,:)

  if (par%fd_order==2) then
      ix1 = 2
      iz1 = 2
  elseif (par%fd_order==4) then
      ix1 = 3
      iz1 = 3
  endif

  if (par%free_surface == 1) then
     ix1=3
     iz1=par%npml+1
  endif

  !! Update stress
  !! --------------------------
  if (par%fd_order==2) then 

  !$omp parallel do private(iz,ix,alpha,alpha2,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx*(lamda(iz,ix)+2.0*mu(iz,ix))
     alpha2=dtdx*lamda(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     xx(iz,ix)=(1.0-kappa1)*xx(iz,ix)+alpha*(c1_elastic_2th*(u(iz,ix)-u(iz,ix-1)))     &
                                            ! c2_elastic_8th*(u(iz,ix+1)-u(iz,ix-2))+   &
                                            ! c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3))+   &
                                            ! c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4))) + &
                                     +alpha2*(c1_elastic_2th*(w(iz,ix)-w(iz-1,ix)))! +    &
                                            ! c2_elastic_8th*(w(iz+1,ix)-w(iz-2,ix))+   &
                                            ! c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+   &
                                            ! c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix)))     
  enddo
  enddo
  !$omp end parallel do

  elseif (par%fd_order==4) then

  !$omp parallel do private(iz,ix,alpha,alpha2,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx*(lamda(iz,ix)+2.0*mu(iz,ix))
     alpha2=dtdx*lamda(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     xx(iz,ix)=(1.0-kappa1)*xx(iz,ix)+alpha*(c1_elastic_4th*(u(iz,ix)-u(iz,ix-1))     &
                                           + c2_elastic_4th*(u(iz,ix+1)-u(iz,ix-2)))  &
                                            ! c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3))+   &
                                            ! c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4))) + &
                                    +alpha2*(c1_elastic_4th*(w(iz,ix)-w(iz-1,ix))    &
                                          +  c2_elastic_4th*(w(iz+1,ix)-w(iz-2,ix)))   
                                            ! c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+   &
                                            ! c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix)))     
  enddo
  enddo
  !$omp end parallel do

  endif


  if (par%fd_order==2) then

  !$omp parallel do private(iz,ix,alpha,alpha2,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx*(lamda(iz,ix)+2.0*mu(iz,ix))
     alpha2=dtdx*lamda(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     zz(iz,ix)=(1.0-kappa1)*zz(iz,ix)+alpha2*(c1_elastic_2th*(u(iz,ix)-u(iz,ix-1)))    &
                                             !c2_elastic_8th*(u(iz,ix+1)-u(iz,ix-2)) + &
                                             !c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3)) + &
                                             !c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4)))+ &
                                     +alpha*(c1_elastic_2th*(w(iz,ix)-w(iz-1,ix)))! +   &
                                            ! c2_elastic_8th*(w(iz+1,ix)-w(iz-2,ix))+  &
                                            ! c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+  &
                                            ! c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix)))
  enddo
  enddo
  !$omp end parallel do

  elseif (par%fd_order==4) then

  !$omp parallel do private(iz,ix,alpha,alpha2,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx*(lamda(iz,ix)+2.0*mu(iz,ix))
     alpha2=dtdx*lamda(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     zz(iz,ix)=(1.0-kappa1)*zz(iz,ix)+alpha2*(c1_elastic_4th*(u(iz,ix)-u(iz,ix-1))    &
                                            + c2_elastic_4th*(u(iz,ix+1)-u(iz,ix-2)))  &
                                             !c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3)) + &
                                             !c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4)))+ &
                                     +alpha* (c1_elastic_4th*(w(iz,ix)-w(iz-1,ix)) & ! +   &
                                            + c2_elastic_4th*(w(iz+1,ix)-w(iz-2,ix)))
                                            ! c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+  &
                                            ! c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix)))
  enddo
  enddo
  !$omp end parallel do

  endif

  if (par%fd_order==2) then

  !$omp parallel do private(iz,ix,alpha,alpha2,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx*mu(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     xz(iz,ix)=(1.0-kappa1)*xz(iz,ix)+alpha*(c1_elastic_2th*(u(iz+1,ix)-u(iz,ix)+w(iz,ix+1)-w(iz,ix)))!+  &
                                             !c2_elastic_8th*(u(iz+2,ix)-u(iz-1,ix)+w(iz,ix+2)-w(iz,ix-1))+  &
                                             !c3_elastic_8th*(u(iz+3,ix)-u(iz-2,ix)+w(iz,ix+3)-w(iz,ix-2))+  &
                                             !c4_elastic_8th*(u(iz+4,ix)-u(iz-3,ix)+w(iz,ix+4)-w(iz,ix-3)))
  enddo
  enddo
  !$omp end parallel do        

  elseif (par%fd_order==4) then

  !$omp parallel do private(iz,ix,alpha,alpha2,kappa1)
  do ix=ix1,nx_pml-2
  do iz=iz1,nz_pml-2
     alpha=dtdx*mu(iz,ix)
     kappa1=damp(iz,ix)*par%dt
     xz(iz,ix)=(1.0-kappa1)*xz(iz,ix)+alpha*(c1_elastic_4th*(u(iz+1,ix)-u(iz,ix)+w(iz,ix+1)-w(iz,ix))  &
                                           + c2_elastic_4th*(u(iz+2,ix)-u(iz-1,ix)+w(iz,ix+2)-w(iz,ix-1)))
                                             !c3_elastic_8th*(u(iz+3,ix)-u(iz-2,ix)+w(iz,ix+3)-w(iz,ix-2))+  &
                                             !c4_elastic_8th*(u(iz+4,ix)-u(iz-3,ix)+w(iz,ix+4)-w(iz,ix-3)))
  enddo
  enddo
  !$omp end parallel do

  endif

end subroutine iso_els_step_sigma_pml

!subroutine iso_els_step_uw_cpml()
!
!end subroutine iso_els_step_uw_cpml


!subroutine iso_els_step_sigma_cpml()
!
!end subroutine iso_els_step_sigma_cpml

subroutine iso_els_step_uw_adj_lagran_pml(nx_pml,nz_pml,npml,is,par,dtdx,den2,lamda,mu,damp,ub,wb,xxb,zzb,xzb)

  implicit none

  type(param),       intent(in)  :: par
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: ub,wb,xxb,zzb,xzb
  real,dimension(nz_pml,nx_pml) :: den2,lamda,mu
  real,              intent(in) :: damp(:,:)

  integer :: ix1, iz1

  if (par%fd_order==2) then
      ix1 = 2
      iz1 = 2
  elseif (par%fd_order==4) then
      ix1 = 3
      iz1 = 3
  endif

  if (par%free_surface == 1) then
     ix1=3
     iz1=par%npml+1
  endif

!  if (par%free_surface == 1) then
!     ix1=3
!     iz1=par%npml+1
!  endif

 !! Update the particle velocity
 !! -------------------------
 if (par%fd_order==2) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx/den2(iz,ix)
    kappa1=damp(iz,ix)*par%dt
    ub(iz,ix)=(1.0-kappa1)*ub(iz,ix)+alpha*(c1_elastic_2th*((lamda(iz,ix+1)+2.0*mu(iz,ix+1))*xxb(iz,ix+1)-(lamda(iz,ix)+2.0*mu(iz,ix))*xxb(iz,ix) &   
                                            +lamda(iz,ix+1)*zzb(iz,ix+1)-lamda(iz,ix)*zzb(iz,ix)+mu(iz,ix)*xzb(iz,ix)-mu(iz-1,ix)*xzb(iz-1,ix)))! +     &
                                            !c2_elastic_8th*((lamda(iz,ix+2)+2.0*mu(iz,ix+2))*xxb(iz,ix+2)-(lamda(iz,ix-1)+2.0*mu(iz,ix-1))*xxb(iz,ix-1)+lamda(iz,ix+2)*zzb(iz,ix+2)-lamda(iz,ix-1)*zzb(iz,ix-1)+mu(iz+1,ix)*xzb(iz+1,ix)-mu(iz-2,ix)*xzb(iz-2,ix)) + &
                                            !c3_elastic_8th*((lamda(iz,ix+3)+2.0*mu(iz,ix+3))*xxb(iz,ix+3)-(lamda(iz,ix-2)+2.0*mu(iz,ix-2))*xxb(iz,ix-2)+lamda(iz,ix+3)*zzb(iz,ix+3)-lamda(iz,ix-2)*zzb(iz,ix-2)+mu(iz+2,ix)*xzb(iz+2,ix)-mu(iz-3,ix)*xzb(iz-3,ix)) + &
                                            !c4_elastic_8th*((lamda(iz,ix+4)+2.0*mu(iz,ix+4))*xxb(iz,ix+4)-(lamda(iz,ix-3)+2.0*mu(iz,ix-3))*xxb(iz,ix-3)+lamda(iz,ix+4)*zzb(iz,ix+4)-lamda(iz,ix-3)*zzb(iz,ix-3)+mu(iz+3,ix)*xzb(iz+3,ix)-mu(iz-4,ix)*xzb(iz-4,ix)))
 enddo
 enddo
 !$omp end parallel do

 elseif (par%fd_order==4) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx/den2(iz,ix)
    kappa1=damp(iz,ix)*par%dt
    ub(iz,ix)=(1.0-kappa1)*ub(iz,ix)+alpha*(c1_elastic_4th*((lamda(iz,ix+1)+2.0*mu(iz,ix+1))*xxb(iz,ix+1)-(lamda(iz,ix)+2.0*mu(iz,ix))*xxb(iz,ix) & 
                                                           + lamda(iz,ix+1)*zzb(iz,ix+1)-lamda(iz,ix)*zzb(iz,ix) + mu(iz,ix)*xzb(iz,ix)-mu(iz-1,ix)*xzb(iz-1,ix))  &
                                          + c2_elastic_4th*((lamda(iz,ix+2)+2.0*mu(iz,ix+2))*xxb(iz,ix+2)-(lamda(iz,ix-1)+2.0*mu(iz,ix-1))*xxb(iz,ix-1) &
                                                           + lamda(iz,ix+2)*zzb(iz,ix+2)-lamda(iz,ix-1)*zzb(iz,ix-1) + mu(iz+1,ix)*xzb(iz+1,ix)-mu(iz-2,ix)*xzb(iz-2,ix)))
                                            !c3_elastic_8th*((lamda(iz,ix+3)+2.0*mu(iz,ix+3))*xxb(iz,ix+3)-(lamda(iz,ix-2)+2.0*mu(iz,ix-2))*xxb(iz,ix-2)+lamda(iz,ix+3)*zzb(iz,ix+3)-lamda(iz,ix-2)*zzb(iz,ix-2)+mu(iz+2,ix)*xzb(iz+2,ix)-mu(iz-3,ix)*xzb(iz-3,ix)) + &
                                            !c4_elastic_8th*((lamda(iz,ix+4)+2.0*mu(iz,ix+4))*xxb(iz,ix+4)-(lamda(iz,ix-3)+2.0*mu(iz,ix-3))*xxb(iz,ix-3)+lamda(iz,ix+4)*zzb(iz,ix+4)-lamda(iz,ix-3)*zzb(iz,ix-3)+mu(iz+3,ix)*xzb(iz+3,ix)-mu(iz-4,ix)*xzb(iz-4,ix)))
 enddo
 enddo
 !$omp end parallel do

 endif

 if (par%fd_order==2) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx/den2(iz,ix)
    kappa1=damp(iz,ix)*par%dt
    wb(iz,ix)=(1.0-kappa1)*wb(iz,ix)+alpha*(c1_elastic_2th*((lamda(iz+1,ix)+2.0*mu(iz+1,ix))*zzb(iz+1,ix)-(lamda(iz,ix)+2.0*mu(iz,ix))*zzb(iz,ix) &           
                                                            +lamda(iz+1,ix)*xxb(iz+1,ix)-lamda(iz,ix)*xxb(iz,ix)+mu(iz,ix)*xzb(iz,ix)-mu(iz,ix-1)*xzb(iz,ix-1)))! +     &
                                           ! c2_elastic_8th*((lamda(iz+2,ix)+2.0*mu(iz+2,ix))*zzb(iz+2,ix)-(lamda(iz-1,ix)+2.0*mu(iz-1,ix))*zzb(iz-1,ix)+lamda(iz+2,ix)*xxb(iz+2,ix)-lamda(iz-1,ix)*xxb(iz-1,ix)+mu(iz,ix+1)*xzb(iz,ix+1)-mu(iz,ix-2)*xzb(iz,ix-2)) + &
                                           ! c3_elastic_8th*((lamda(iz+3,ix)+2.0*mu(iz+3,ix))*zzb(iz+3,ix)-(lamda(iz-2,ix)+2.0*mu(iz-2,ix))*zzb(iz-2,ix)+lamda(iz+3,ix)*xxb(iz+3,ix)-lamda(iz-2,ix)*xxb(iz-2,ix)+mu(iz,ix+2)*xzb(iz,ix+2)-mu(iz,ix-3)*xzb(iz,ix-3)) + &
                                           ! c4_elastic_8th*((lamda(iz+4,ix)+2.0*mu(iz+4,ix))*zzb(iz+4,ix)-(lamda(iz-3,ix)+2.0*mu(iz-3,ix))*zzb(iz-3,ix)+lamda(iz+4,ix)*xxb(iz+4,ix)-lamda(iz-3,ix)*xxb(iz-3,ix)+mu(iz,ix+3)*xzb(iz,ix+3)-mu(iz,ix-4)*xzb(iz,ix-4)))
 enddo
 enddo
 !$omp end parallel do        

 elseif (par%fd_order==4) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx/den2(iz,ix)
    kappa1=damp(iz,ix)*par%dt
    wb(iz,ix)=(1.0-kappa1)*wb(iz,ix)+alpha*(c1_elastic_4th*((lamda(iz+1,ix)+2.0*mu(iz+1,ix))*zzb(iz+1,ix)-(lamda(iz,ix)+2.0*mu(iz,ix))*zzb(iz,ix) &
                                                           + lamda(iz+1,ix)*xxb(iz+1,ix)-lamda(iz,ix)*xxb(iz,ix) + mu(iz,ix)*xzb(iz,ix)-mu(iz,ix-1)*xzb(iz,ix-1)) &
                                          + c2_elastic_4th*((lamda(iz+2,ix)+2.0*mu(iz+2,ix))*zzb(iz+2,ix)-(lamda(iz-1,ix)+2.0*mu(iz-1,ix))*zzb(iz-1,ix) &                                
                                                           + lamda(iz+2,ix)*xxb(iz+2,ix)-lamda(iz-1,ix)*xxb(iz-1,ix) + mu(iz,ix+1)*xzb(iz,ix+1)-mu(iz,ix-2)*xzb(iz,ix-2))) 
                                           ! c3_elastic_8th*((lamda(iz+3,ix)+2.0*mu(iz+3,ix))*zzb(iz+3,ix)-(lamda(iz-2,ix)+2.0*mu(iz-2,ix))*zzb(iz-2,ix)+lamda(iz+3,ix)*xxb(iz+3,ix)-lamda(iz-2,ix)*xxb(iz-2,ix)+mu(iz,ix+2)*xzb(iz,ix+2)-mu(iz,ix-3)*xzb(iz,ix-3)) + &
                                           ! c4_elastic_8th*((lamda(iz+4,ix)+2.0*mu(iz+4,ix))*zzb(iz+4,ix)-(lamda(iz-3,ix)+2.0*mu(iz-3,ix))*zzb(iz-3,ix)+lamda(iz+4,ix)*xxb(iz+4,ix)-lamda(iz-3,ix)*xxb(iz-3,ix)+mu(iz,ix+3)*xzb(iz,ix+3)-mu(iz,ix-4)*xzb(iz,ix-4)))
 enddo
 enddo
 !$omp end parallel do         

 endif

end subroutine iso_els_step_uw_adj_lagran_pml

subroutine iso_els_step_sigma_adj_lagran_pml(nx_pml,nz_pml,npml,is,par,dtdx,damp,ub,wb,xxb,zzb,xzb)

  implicit none

  type(param),       intent(in)  :: par
  integer,           intent(in)  :: is,nx_pml,nz_pml,npml

  real :: alpha,alpha2,kappa1,dtdx

  real,dimension(nz_pml,nx_pml) :: ub,wb,xxb,zzb,xzb
  real,              intent(in) :: damp(:,:)

  integer :: ix1, iz1

  if (par%fd_order==2) then
      ix1 = 2
      iz1 = 2
  elseif (par%fd_order==4) then
      ix1 = 3
      iz1 = 3
  endif

  if (par%free_surface == 1) then
     ix1=3
     iz1=par%npml+1
  endif

 !! Update stress
 !! --------------------------
 if (par%fd_order==2) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    xxb(iz,ix)=(1.0-kappa1)*xxb(iz,ix)+alpha*(c1_elastic_2th*(ub(iz,ix)-ub(iz,ix-1)))! +  &
                                              !c2_elastic_8th*(ub(iz,ix+1)-ub(iz,ix-2))+ &
                                              !c3_elastic_8th*(ub(iz,ix+2)-ub(iz,ix-3))+ &
                                              !c4_elastic_8th*(ub(iz,ix+3)-ub(iz,ix-4)))
 enddo
 enddo
 !$omp end parallel do

 elseif (par%fd_order==4) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    xxb(iz,ix)=(1.0-kappa1)*xxb(iz,ix)+alpha*(c1_elastic_4th*(ub(iz,ix)-ub(iz,ix-1)) &
                                            + c2_elastic_4th*(ub(iz,ix+1)-ub(iz,ix-2)))
                                              !c3_elastic_8th*(ub(iz,ix+2)-ub(iz,ix-3))+ &
                                              !c4_elastic_8th*(ub(iz,ix+3)-ub(iz,ix-4)))
 enddo
 enddo
 !$omp end parallel do         

 endif

 if (par%fd_order==2) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    zzb(iz,ix)=(1.0-kappa1)*zzb(iz,ix)+alpha*(c1_elastic_2th*(wb(iz,ix)-wb(iz-1,ix)))! +   &
                                              !c2_elastic_8th*(wb(iz+1,ix)-wb(iz-2,ix)) + &
                                              !c3_elastic_8th*(wb(iz+2,ix)-wb(iz-2,ix)) + &
                                              !c4_elastic_8th*(wb(iz+3,ix)-wb(iz-3,ix)))
 enddo
 enddo
 !$omp end parallel do

 elseif (par%fd_order==4) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    zzb(iz,ix)=(1.0-kappa1)*zzb(iz,ix)+alpha*(c1_elastic_4th*(wb(iz,ix)-wb(iz-1,ix)) &
                                            + c2_elastic_4th*(wb(iz+1,ix)-wb(iz-2,ix))) 
                                              !c3_elastic_8th*(wb(iz+2,ix)-wb(iz-2,ix)) + &
                                              !c4_elastic_8th*(wb(iz+3,ix)-wb(iz-3,ix)))
 enddo
 enddo
 !$omp end parallel do

 endif

 if (par%fd_order==2) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    xzb(iz,ix)=(1.0-kappa1)*xzb(iz,ix)+alpha*(c1_elastic_2th*(ub(iz+1,ix)-ub(iz,ix)+wb(iz,ix+1)-wb(iz,ix)))! +   &
                                              !c2_elastic_8th*(ub(iz+2,ix)-ub(iz-1,ix)+wb(iz,ix+2)-wb(iz,ix-1)) + &
                                              !c3_elastic_8th*(ub(iz+3,ix)-ub(iz-2,ix)+wb(iz,ix+3)-wb(iz,ix-2)) + &
                                              !c4_elastic_8th*(ub(iz+4,ix)-ub(iz-3,ix)+wb(iz,ix+4)-wb(iz,ix-3)))
 enddo
 enddo
 !$omp end parallel do

 elseif (par%fd_order==4) then

 !$omp parallel do private(iz,ix,alpha,kappa1)
 do ix=ix1,nx_pml-2
 do iz=iz1,nz_pml-2
    alpha=dtdx
    kappa1=damp(iz,ix)*par%dt
    xzb(iz,ix)=(1.0-kappa1)*xzb(iz,ix)+alpha*(c1_elastic_4th*(ub(iz+1,ix)-ub(iz,ix)+wb(iz,ix+1)-wb(iz,ix)) &
                                            + c2_elastic_4th*(ub(iz+2,ix)-ub(iz-1,ix)+wb(iz,ix+2)-wb(iz,ix-1)))
                                              !c3_elastic_8th*(ub(iz+3,ix)-ub(iz-2,ix)+wb(iz,ix+3)-wb(iz,ix-2)) + &
                                              !c4_elastic_8th*(ub(iz+4,ix)-ub(iz-3,ix)+wb(iz,ix+4)-wb(iz,ix-3)))
 enddo
 enddo
 !$omp end parallel do

 endif

end subroutine iso_els_step_sigma_adj_lagran_pml


end module staggered_elastic_pml_time_step
