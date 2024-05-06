
module staggered_elastic_kernel

 use global
 use datatype
 use math
 use string
 use io

 use staggered_elastic_pml_time_step

 use wave_bound

 implicit none

 ! real, private,         allocatable  :: u(:,:),w(:,:),p(:,:),u1(:,:),w1(:,:),p1(:,:)
 integer, private                    :: isx,isz,igx,igz,iz,ix,it,ig

contains


          !elastic_modeling(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
!          wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_true_u,seis_true_w,store_wave_flag)
subroutine modeling_elastic_kernel(is,par,coord,s,vp,vs,den,fs,fd_order,nx_pml,nz_pml,npml,damp,&
           wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_u,seis_w,store_wave_flag)

 implicit none

 type(param),       intent(in)  :: par
 type(acquisition), intent(in)  :: coord
 type(wave_boundary_param), intent(inout)  :: wave_boundary
 integer,           intent(in)  :: is,fs,nx_pml,nz_pml,npml,fd_order
 real,              intent(in)  :: s(:),vp(:,:),vs(:,:),den(:,:),damp(:,:)
! integer,           intent(in)  :: save_wavefield, save_wavefield_boundary
 real,              intent(out) :: seis_u(:,:),seis_w(:,:),wave_pux(:,:,:),wave_puz(:,:,:),wave_pwx(:,:,:),wave_pwz(:,:,:)

 real,allocatable :: lamda(:,:),mu(:,:),den2(:,:)
 real,allocatable :: u(:,:),w(:,:),xx(:,:),zz(:,:),xz(:,:)
 real :: alpha,alpha2,kappa1,dtdx

 real :: tmp_u, tmp_w

 integer  ::  pad_top

 character snap_str*5

 logical :: store_wave_flag


 pad_top = npml

 allocate(lamda(nz_pml,nx_pml),mu(nz_pml,nx_pml),den2(nz_pml,nx_pml))
 lamda=0.0;mu=0.0;den2=0.0
 do ix=1,nx_pml
 do iz=1,nz_pml
    den2(iz,ix)=den(iz,ix)
    mu(iz,ix)=den(iz,ix)*vs(iz,ix)*vs(iz,ix)
    lamda(iz,ix)=den(iz,ix)*vp(iz,ix)*vp(iz,ix)-2.0*mu(iz,ix)
 enddo
 enddo

 print *, 'elastic_modeling test0', 'par%lv', par%lv

 print *, 'elastic_modeling save_wavefield', par%save_wavefield, is
 print *, 'par%nx, par%nz, par%nt in elastic_modeling ', par%nx, par%nz, par%nt

 !! Initialize the outputs
 !! ----------------------
 seis_u=0.0;seis_w=0.0

 if (store_wave_flag) then  

 if(par%save_wavefield) then
    !$omp parallel do private(iz,ix)
    do ix=1,par%nx
    do iz=1,par%nz
       wave_pux(iz,ix,1:par%nt)=0.0; wave_puz(iz,ix,1:par%nt)=0.0
       wave_pwx(iz,ix,1:par%nt)=0.0; wave_pwz(iz,ix,1:par%nt)=0.0
    enddo
    enddo
    !$omp end parallel do
 endif

 endif

 print *, 'finish initialization in elastic_modeling'

 !! Memory allocations
 !! ------------------
 allocate(u(nz_pml,nx_pml),w(nz_pml,nx_pml),xx(nz_pml,nx_pml),zz(nz_pml,nx_pml),xz(nz_pml,nx_pml))

 !! Initialize arrays
 !! -----------------
 !$omp parallel do private(iz,ix)
 do ix=1,nx_pml
 do iz=1,nz_pml
    u(iz,ix)=0.0;w(iz,ix)=0.0
    xx(iz,ix)=0.0;zz(iz,ix)=0.0;xz(iz,ix)=0.0
 enddo
 enddo
 !$omp end parallel do

 print *, 'elastic_modeling test1'

 !! Set up source position
 !! ----------------------
 if(fs==0)then
    isx=npml+int(coord%xs(is)/par%dx)
    isz=npml+int(coord%zs(is)/par%dx)
 else
    isx=npml+int(coord%xs(is)/par%dx)
    isz=pad_top+int(coord%zs(is)/par%dx)
 endif    
 
 dtdx=par%dt/par%dx !! Local variable for reducing flops

 !! free surface condition
! if(fs==1)then
!   do ix=1,nx_pml
!      mu(pad_top+1,ix)=1.0*mu(pad_top+1,ix)
!      lamda(pad_top+1,ix)=0.0
!      den2(pad_top+1,ix)=den2(pad_top+1,ix)/2.0     
!   enddo
! endif        

 !! Time marching loop
 !! ------------------
 do it=1,par%nt + par%nt_shift

  !! Add source
  !! ----------
  if (it < par%nt + 1) then 
    if (par%source_mec == 1) then

       u(isz,isx)=u(isz,isx)+s(it)
       w(isz,isx)=w(isz,isx)+s(it)

    elseif (par%source_mec == 2) then
   
       w(isz,isx)=w(isz,isx)+s(it)
  
    elseif (par%source_mec == 3) then
            
       xx(isz,isx) = xx(isz,isx) + s(it)
       zz(isz,isx) = zz(isz,isx) + s(it)

    endif        
  endif

  !! Update the particle velocity
  !! -------------------------
  call iso_els_step_uw_pml(nx_pml,nz_pml,npml,is,par,dtdx,den2,damp,u,w,xx,zz,xz) 

  if (fs==1) then

    if (par%fs_zero == 1) then
     u(:pad_top-1,:) = 0.
     w(:pad_top-1,:) = 0.      
    endif

  endif

  !! Update stress
  !! --------------------------
      !iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)
  call iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)

!  print *, 'fsfs', fs

  !free surface bc
  if(fs==1) then

      if (mod((it-par%nt_shift),1000).eq.0) then 
        print *, 'fs', fs, 'par%fs_type', par%fs_type, 'par%fs_zero',par%fs_zero
      endif

     iz = pad_top+1

     if (par%fs_type == 1) then

            !! currently use this one
            !free_surface_denise(nx_pml,nz_pml,iz,dt,dx,fd_order,lamda,mu,xx,zz,xz,u,w)
        call free_surface_denise(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,lamda,mu,xx,zz,xz,u,w)

     elseif (par%fs_type == 2) then
            !free_surface_legacy(nx_pml,nz_pml,iz,dt,dx,zz)
        call free_surface_legacy(nx_pml,nz_pml,iz,par%dt,par%dx,zz) 

     elseif (par%fs_type == 3) then
            !free_surface_JIMU_ss(nx_pml,nz_pml,iz,dt,dx,zz,xz)
        call free_surface_JIMU_ss(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,zz,xz) 

     endif

     if (par%fs_zero == 1) then
      xx(:pad_top-1,:) = 0.
      zz(:pad_top-1,:) = 0.
      xz(:pad_top-1,:) = 0.
     endif

  endif  !! free surface bc 
  
  !! Output pressure seismogram 
  !! --------------------------
  if (it > par%nt_shift) then 
  
    do ig=1,coord%ng(is)
       if(fs==0)then
          igx=npml+int(coord%xg(is,ig)/par%dx)
          igz=npml+int(coord%zg(is,ig)/par%dx)
       else
          igx=npml+int(coord%xg(is,ig)/par%dx)
          igz=pad_top+int(coord%zg(is,ig)/par%dx)
       endif
       seis_u(it-par%nt_shift,ig)=u(igz,igx)
       seis_w(it-par%nt_shift,ig)=w(igz,igx)
    enddo

    !! only wave wavefield (either the entire wavefield or only boundaries) 
    !! if not calculating step length
    
    if (store_wave_flag) then 

      if(par%save_wavefield) then
        if (mod((it-par%nt_shift), 5000).eq.0) print *, 'fuck', par%save_wavefield, is
        !$omp parallel do private(iz,ix)
        do ix=npml+1,nx_pml-npml
           do iz=pad_top+1,nz_pml-npml
           !do iz=npml+1,nz_pml-npml
            wave_pux(iz-pad_top,ix-npml,it-par%nt_shift)=c1_elastic_2th*(u(iz,ix)-u(iz,ix-1))!+c2_elastic_8th*(u(iz,ix+1)-u(iz,ix-2))+  &
                                                          !c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3))+c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4))
            wave_puz(iz-pad_top,ix-npml,it-par%nt_shift)=c1_elastic_2th*(u(iz+1,ix)-u(iz,ix))!+c2_elastic_8th*(u(iz+2,ix)-u(iz-1,ix))+  &
                                                          !c3_elastic_8th*(u(iz+3,ix)-u(iz-2,ix))+c4_elastic_8th*(u(iz+4,ix)-u(iz-3,ix))
            wave_pwx(iz-pad_top,ix-npml,it-par%nt_shift)=c1_elastic_2th*(w(iz,ix+1)-w(iz,ix))!+c2_elastic_8th*(w(iz,ix+2)-w(iz,ix-1))+  &
                                                          !c3_elastic_8th*(w(iz,ix+3)-w(iz,ix-2))+c4_elastic_8th*(w(iz,ix+4)-w(iz,ix-3))
            wave_pwz(iz-pad_top,ix-npml,it-par%nt_shift)=c1_elastic_2th*(w(iz,ix)-w(iz-1,ix))!+c2_elastic_8th*(w(iz+1,ix)-w(iz-2,ix))+  &
                                                          !c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix))                        
           enddo
        enddo
        !$omp end parallel do

!      write(snap_str,'(i5.5)')  it
!
!      open(99999, file='snap_0_zz'//snap_str, access='direct', form='unformatted', recl=4*(par%nz))
!      do ix=1,par%nx
!         write(99999,rec=ix) ((wave_pux(iz,ix,it)+wave_pwz(iz,ix,it)),iz=1,par%nz)
!      enddo
!      close(99999)

      endif !! end if par%save_wavefield

      if(par%save_wavefield_boundary) then

         call record_wave_boundary_uw(it-par%nt_shift,par%nt,u,w,&
                                      wave_boundary%u_bl,wave_boundary%u_br,wave_boundary%u_bt,wave_boundary%u_bb,&
                                      wave_boundary%w_bl,wave_boundary%w_br,wave_boundary%w_bt,wave_boundary%w_bb,&
                                      nx_pml,nz_pml, par%nx, par%nz, par%lv,npml)

         if((it-par%nt_shift).eq.par%nt) then 
 
            wave_boundary%last_u(:,:) = u(npml+1:nz_pml-npml, npml+1:nx_pml-npml)
            wave_boundary%last_w(:,:) = w(npml+1:nz_pml-npml, npml+1:nx_pml-npml)

         endif !!! if((it-par%nt_shift).eq.par%nt)

         !!! for test
      !! it works but not necessary to save stress wavefield
      !   call record_wave_boundary_stress(it-par%nt_shift,par%nt,xx,zz,xz,&
      !                                wave_boundary%xx_bl,wave_boundary%xx_br,wave_boundary%xx_bt,wave_boundary%xx_bb,&
      !                                wave_boundary%zz_bl,wave_boundary%zz_br,wave_boundary%zz_bt,wave_boundary%zz_bb,&
      !                                wave_boundary%xz_bl,wave_boundary%xz_br,wave_boundary%xz_bt,wave_boundary%xz_bb,&
      !                                nx_pml,nz_pml, par%nx, par%nz, par%lv,npml)

         if((it-par%nt_shift).eq.par%nt) then
        
            wave_boundary%last_xx(:,:) = xx(npml+1:nz_pml-npml, npml+1:nx_pml-npml)
            wave_boundary%last_zz(:,:) = zz(npml+1:nz_pml-npml, npml+1:nx_pml-npml)
            wave_boundary%last_xz(:,:) = xz(npml+1:nz_pml-npml, npml+1:nx_pml-npml)

         endif !!! if((it-par%nt_shift).eq.par%nt)
         !!! end for test

      endif !! endif save_wavefield_boundary

    endif !! (store_wave_flag)

  endif !! endif (it > par%nt_shift)

  if (.false.) then 
!  if(par%save_wavefield_boundary) then

    if(is==10 .and. mod((it-par%nt_shift),500).eq.0.and.it >par%nt_shift) then
      write(snap_str,'(i5.5)')  (it-par%nt_shift)

      open(99999, file='snap_u'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (u(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo
      close(99999)

      open(99999, file='snap_w'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (w(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo
      close(99999)

      open(99999, file='snap_xx'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (xx(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo
      close(99999)

      open(99999, file='snap_zz'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (zz(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo
      close(99999)

      open(99999, file='snap_xz'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (xz(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo
      close(99999)

    endif

!  endif
  endif !! endif false

 enddo  !! End of time-marching loop

 print *, 'elastic_modeling test3', 'par%lv', par%lv

 print *, 'elastic_modeling save_wavefield', par%save_wavefield, is


 !deallocate(u,w,p)
 deallocate(lamda,mu,den2)
 deallocate(u,w,xx,zz,xz)

end subroutine modeling_elastic_kernel


!! compute the gradient of elastic FWI for vs
!          staggered_elastic_fwi_kernel(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
!          wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,shot_virtual,energy,gk)
subroutine staggered_elastic_fwi_kernel(is,par,coord,s,vp,vs,den,fs,fd_order,nx_pml,nz_pml,npml,damp,&
           wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_v,energy,mig)

 implicit none

 type(param),       intent(in)    :: par
 type(acquisition), intent(in)    :: coord
 integer,           intent(in)    :: is,fs,nx_pml,nz_pml,npml,fd_order
 real,              intent(in)    :: s(:),vp(:,:),vs(:,:),den(:,:),damp(:,:), &
                                     seis_v(:,:),wave_pux(:,:,:),wave_puz(:,:,:),wave_pwx(:,:,:),wave_pwz(:,:,:)
 real,              intent(inout) :: mig(:,:),energy(:,:)

 real,allocatable                 :: lamda(:,:),mu(:,:),den2(:,:)
 real,allocatable                 :: ub(:,:),wb(:,:),xxb(:,:),zzb(:,:),xzb(:,:)
 real,allocatable                 :: u(:,:),w(:,:),xx(:,:),zz(:,:),xz(:,:)
 real                             :: alpha,alpha2,kappa1,dtdx
 integer                          :: pad_top 

! integer,           intent(in)  :: save_wavefield,save_wavefield_boundary

  character is_str*4

 type(wave_boundary_param), intent(in)  :: wave_boundary

 real,allocatable ::wave_pux_fly(:,:),wave_puz_fly(:,:),wave_pwx_fly(:,:),wave_pwz_fly(:,:)

 !!
 real,allocatable ::grdla(:,:),grdmu(:,:),grdvp(:,:),grdvs(:,:),grdmu_a(:,:),grdmu_b(:,:),grdmu_c(:,:)

 !! 
 real,allocatable ::souresp(:)

 character snap_str*5


 pad_top=npml

 !! Initialize the outputs
 !! ----------------------
 !$omp parallel do private(ix,iz)
 do ix=1,nx_pml
 do iz=1,nz_pml
    mig(iz,ix)=0.0; energy(iz,ix)=0.0
 enddo
 enddo
 !$omp end parallel do

 print *, 'staggered_elastic_fwi test0'

 allocate(lamda(nz_pml,nx_pml),mu(nz_pml,nx_pml),den2(nz_pml,nx_pml))
 lamda=0.0;mu=0.0;den2=0.0
 do ix=1,nx_pml
 do iz=1,nz_pml
    den2(iz,ix)=den(iz,ix)
    mu(iz,ix)=den(iz,ix)*vs(iz,ix)*vs(iz,ix)
    lamda(iz,ix)=den(iz,ix)*vp(iz,ix)*vp(iz,ix)-2.0*mu(iz,ix)
 enddo
 enddo

 !! Memory allocations
 !! ------------------ 
 allocate(ub(nz_pml,nx_pml),wb(nz_pml,nx_pml),xxb(nz_pml,nx_pml),zzb(nz_pml,nx_pml),xzb(nz_pml,nx_pml))

 if(par%save_wavefield_boundary) then 
    allocate(u(nz_pml,nx_pml),w(nz_pml,nx_pml),xx(nz_pml,nx_pml),zz(nz_pml,nx_pml),xz(nz_pml,nx_pml))
    u = 0.; w = 0.
    xx = 0.; zz = 0.; xz = 0.

    allocate(wave_pux_fly(1:par%nz,1:par%nx),wave_puz_fly(1:par%nz,1:par%nx),&
             wave_pwx_fly(1:par%nz,1:par%nx),wave_pwz_fly(1:par%nz,1:par%nx))
    wave_pux_fly = 0.
    wave_puz_fly = 0.
    wave_pwx_fly = 0.
    wave_pwz_fly = 0.
 endif

! if (par%adj_lagran==0) then
! endif

 allocate(grdla(nz_pml,nx_pml),grdmu(nz_pml,nx_pml),&
          grdmu_a(nz_pml,nx_pml),&
          grdmu_b(nz_pml,nx_pml),grdmu_c(nz_pml,nx_pml))

 allocate(grdvp(nz_pml,nx_pml),grdvs(nz_pml,nx_pml))

 grdla = 0.0
 grdmu = 0.0
 grdmu_a = 0.0
 grdmu_b = 0.0
 grdmu_c = 0.0

 grdvp = 0.0
 grdvs = 0.0

 allocate(souresp(1:coord%ng(is)))
 souresp = 0.

 !$omp parallel do private(ix,iz)
 do ix=1,nx_pml
 do iz=1,nz_pml
    ub(iz,ix)=0.0; wb(iz,ix)=0.0
    xxb(iz,ix)=0.0; zzb(iz,ix)=0.0; xzb(iz,ix)=0.0
 enddo
 enddo
 !$omp end parallel do

 dtdx=par%dt/par%dx !! Local variable for reducing flops

 !! free surface condition
! if(fs==1)then
!   do ix=1,nx_pml
!      mu(pad_top+1,ix)=1.0*mu(pad_top+1,ix)
!      lamda(pad_top+1,ix)=0.0
!      den2(pad_top+1,ix)=den2(pad_top+1,ix)/2.0
!   enddo
! endif

 print *, 'staggered_elastic_fwi test1'

 dtdx = -dtdx

 print *, 'par%dt',par%dt, 'fs', fs, 'par%fs_adj', par%fs_adj, 'par%fs_type', par%fs_type, &
          'par%adj_lagran', par%adj_lagran  

 if (is.eq.40) then 
 if(.True.) then
    write(is_str,'(i4.4)')  is
    open(88888,file='shot_w_adj.bin'//is_str,access='direct', form='unformatted', recl=4*par%nt)
    do ig=1,coord%ng(is)
       write(88888, rec=ig) (seis_v(it,ig), it=1,par%nt)
    enddo
    close(88888)
 endif
 endif

 !! Time marching loop 
 !! ------------------
 do it=par%nt+par%nt_shift,par%nt_shift+1,-1

 !free surface bc
!  if(fs==1)then
!     do ix=1,nx_pml
!        zzb(pad_top,ix)=0.0
!     enddo
!  endif

 !! Update the particle velocity
 !! -------------------------
 if (par%adj_lagran==1) then 

    call iso_els_step_uw_adj_lagran_pml(nx_pml,nz_pml,npml,is,par,dtdx,den2,lamda,mu,damp,ub,wb,xxb,zzb,xzb)

 else

    call iso_els_step_uw_pml(nx_pml,nz_pml,npml,is,par,dtdx,den2,damp,ub,wb,xxb,zzb,xzb) 

 endif

 !! Add seismogram
 !! ---------------
 if (par%adj_lagran==1) then

   do ig=1,coord%ng(is)
      igx=npml+int(coord%xg(is,ig)/par%dx)
      igz=pad_top+int(coord%zg(is,ig)/par%dx)
      wb(igz,igx)=wb(igz,igx)+seis_v(it-par%nt_shift,ig)*par%dt/den(igz,igx)
   enddo

 else

   souresp = souresp + seis_v(it-par%nt_shift,:) 

   do ig=1,coord%ng(is)
      igx=npml+int(coord%xg(is,ig)/par%dx)
      igz=pad_top+int(coord%zg(is,ig)/par%dx)
      wb(igz,igx)=wb(igz,igx) + souresp(ig)
   enddo

 endif
!!
 !! Update stress
 !! --------------------------
 if (par%adj_lagran==1) then 
 
    call iso_els_step_sigma_adj_lagran_pml(nx_pml,nz_pml,npml,is,par,dtdx,damp,ub,wb,xxb,zzb,xzb)

 else
        !iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)
    call iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,ub,wb,xxb,zzb,xzb) 

 endif

 !free surface bc
!  if(fs==1)then
!     do ix=1,nx_pml
!        zzb(pad_top,ix)=0.0
!     enddo
!  endif

  if(fs==1) then

     iz = pad_top + 1

     if (par%fs_adj == 1) then 

        if (par%fs_type == 1) then
                                 !(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,lamda,mu,xx,zz,xz,u,w)
          call free_surface_denise(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,lamda,mu,xxb,zzb,xzb,ub,wb)

        elseif (par%fs_type == 2) then
             !free_surface_legacy(nx_pml,nz_pml,iz,dt,dx,zz)
          call free_surface_legacy(nx_pml,nz_pml,iz,par%dt,par%dx,zzb)

        elseif (par%fs_type == 3) then
            !free_surface_JIMU_ss(nx_pml,nz_pml,iz,dt,dx,zz,xz)
          call free_surface_JIMU_ss(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,zzb,xzb)

        endif

     endif   !! endif adjoint free surface

  endif


  !! save the entire wavefield is only applicable when the model is very small
  if(par%save_wavefield) then 

   if (par%adj_lagran==1) then
    !! Dot product imaging condition
    !! -----------------------------
    !$omp parallel do private(iz,ix)
    do ix=npml+1,nx_pml-npml
    do iz=pad_top+1,nz_pml-npml
  !   mig(iz,ix)=mig(iz,ix)+4.0*den(iz,ix)*vs(iz,ix)*(wave_pwz(iz-pad_top,ix-npml,it)*xxb(iz,ix)+wave_pux(iz-pad_top,ix-npml,it)*zzb(iz,ix))-  & 
  !                         2.0*den(iz,ix)*vs(iz,ix)*(wave_puz(iz-pad_top,ix-npml,it)+wave_pwx(iz-pad_top,ix-npml,it))*xzb(iz,ix)
       mig(iz,ix)=mig(iz,ix)+par%dt*4.0*vs(iz,ix)*(wave_pux(iz-pad_top,ix-npml,it-par%nt_shift)+wave_pwz(iz-pad_top,ix-npml,it-par%nt_shift))*(xxb(iz,ix)+zzb(iz,ix)) & 
                          -par%dt*2.0*vs(iz,ix)*2.0*(xxb(iz,ix)*wave_pux(iz-pad_top,ix-npml,it-par%nt_shift)+zzb(iz,ix)*wave_pwz(iz-pad_top,ix-npml,it-par%nt_shift)+0.5*xzb(iz,ix)*(wave_puz(iz-pad_top,ix-npml,it-par%nt_shift)+wave_pwx(iz-pad_top,ix-npml,it-par%nt_shift)))
  !     mig(iz,ix)=mig(iz,ix)+2.0*den(iz,ix)*c(iz,ix)*wave(iz-npml,ix-npml,it)*p1(iz,ix)
    enddo
    enddo
    !$omp end parallel do

    !! Save illumination
    !! -----------------
    !$omp parallel do private(iz,ix)
    do ix=npml+1,nx_pml-npml
    do iz=pad_top+1,nz_pml-npml
       energy(iz,ix)=energy(iz,ix)+(wave_pux(iz-pad_top,ix-npml,it-par%nt_shift)+wave_pwz(iz-pad_top,ix-npml,it-par%nt_shift))**2
    enddo
    enddo
    !$omp end parallel do

   else

!      if (par%grad_twist==1) then

        call crosco(par%nx,par%nz,npml,nx_pml,nz_pml,pad_top,xxb,zzb,xzb,xx,zz,xz,grdla,grdmu_a,grdmu_b,grdmu_c)

!      else

!        call crosco_denise(par%nx,par%nz,npml,nx_pml,nz_pml,pad_top,xxb,zzb,xzb,xx,zz,xz,grdla,grdmu_a,grdmu_b,grdmu_c)

!      endif !! endif par%grad_twist

     !! Save illumination
     !! -----------------

     !$omp parallel do private(iz,ix)
     do ix=npml+1,nx_pml-npml
        do iz=npml+1,nz_pml-npml
           energy(iz,ix)=energy(iz,ix)+(u(iz,ix)*u(iz,ix)+w(iz,ix)*w(iz,ix))
        enddo
     enddo
     !$omp end parallel do

   endif !! endif par%adj_lagran==1 

  elseif(par%save_wavefield_boundary) then 


  if(.true.) then 

  !! free surface boundary condition 
  !! use with care
  if(fs==1) then

     iz = pad_top + 1

     if (par%fs_type == 1) then
                               !(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,lamda,mu,xx,zz,xz,u,w)
        call free_surface_denise(nx_pml,nz_pml,iz,-par%dt,par%dx,fd_order,lamda,mu,xx,zz,xz,u,w)

      elseif (par%fs_type == 2) then
            !free_surface_legacy(nx_pml,nz_pml,iz,dt,dx,zz)
        call free_surface_legacy(nx_pml,nz_pml,iz,par%dt,par%dx,zz)

      elseif (par%fs_type == 3) then
            !free_surface_JIMU_ss(nx_pml,nz_pml,iz,dt,dx,zz,xz)
        call free_surface_JIMU_ss(nx_pml,nz_pml,iz,par%dt,par%dx,fd_order,zz,xz)

      endif

!      xx(:pad_top-1,:) = 0.
!      zz(:pad_top-1,:) = 0.
!      xz(:pad_top-1,:) = 0.

  endif 

  endif !! endif false 

  !! Update stress
  !! --------------------------
      !iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)
  call iso_els_step_sigma_pml(nx_pml,nz_pml,npml,is,par,dtdx,lamda,mu,damp,u,w,xx,zz,xz)

  !! Update the particle velocity
  !! -------------------------
  call iso_els_step_uw_pml(nx_pml,nz_pml,npml,is,par,dtdx,den2,damp,u,w,xx,zz,xz) 

  if (fs==1) then

!     u(:pad_top-1,:) = 0.
!     w(:pad_top-1,:) = 0.      

  endif

!!!! bvr
!  if (.false.) then 
!  do iz=1, par%nz
!     do ix=1, par%lv
!!!      left
!        u(iz+npml,npml+ix+1)=wave_boundary%u_bl(iz,ix,it-par%nt_shift)
!        w(iz+npml,npml+ix+1)=wave_boundary%w_bl(iz,ix,it-par%nt_shift)
!!!      right
!        u(iz+npml,nx_pml-npml-par%lv+ix-1)=wave_boundary%u_br(iz,ix,it-par%nt_shift)
!        w(iz+npml,nx_pml-npml-par%lv+ix-1)=wave_boundary%w_br(iz,ix,it-par%nt_shift)
!
!!!      left
!        xx(iz+npml,npml+ix+1)=wave_boundary%xx_bl(iz,ix,it-par%nt_shift)
!        zz(iz+npml,npml+ix+1)=wave_boundary%zz_bl(iz,ix,it-par%nt_shift)
!        xz(iz+npml,npml+ix+1)=wave_boundary%xz_bl(iz,ix,it-par%nt_shift)
!!!      right
!        xx(iz+npml,nx_pml-npml-par%lv+ix-1)=wave_boundary%xx_br(iz,ix,it-par%nt_shift)
!        zz(iz+npml,nx_pml-npml-par%lv+ix-1)=wave_boundary%zz_br(iz,ix,it-par%nt_shift)
!        xz(iz+npml,nx_pml-npml-par%lv+ix-1)=wave_boundary%xz_br(iz,ix,it-par%nt_shift)
!
!     enddo
!  enddo
!i,j
!  do ix=1, par%nx
!     do iz=1, par%lv
!!!       top
!        u(npml+iz+1,ix+npml)=wave_boundary%u_bt(ix,iz,it-par%nt_shift)
!        w(npml+iz+1,ix+npml)=wave_boundary%w_bt(ix,iz,it-par%nt_shift)
!!!       bottom
!        u(nz_pml-npml-par%lv+iz-1,ix+npml)=wave_boundary%u_bb(ix,iz,it-par%nt_shift)
!        w(nz_pml-npml-par%lv+iz-1,ix+npml)=wave_boundary%w_bb(ix,iz,it-par%nt_shift)
!
!!!       top
!        xx(npml+iz+1,ix+npml)=wave_boundary%xx_bt(ix,iz,it-par%nt_shift)
!        zz(npml+iz+1,ix+npml)=wave_boundary%zz_bt(ix,iz,it-par%nt_shift)
!        xz(npml+iz+1,ix+npml)=wave_boundary%xz_bt(ix,iz,it-par%nt_shift)
!!!       bottom
!        xx(nz_pml-npml-par%lv+iz-1,ix+npml)=wave_boundary%xx_bb(ix,iz,it-par%nt_shift)
!        zz(nz_pml-npml-par%lv+iz-1,ix+npml)=wave_boundary%zz_bb(ix,iz,it-par%nt_shift)
!        xz(nz_pml-npml-par%lv+iz-1,ix+npml)=wave_boundary%xz_bb(ix,iz,it-par%nt_shift)
!
!     enddo
!  enddo
!  endif

  if (.true.) then
!!!! bvr
  do iz=1, par%nz
     do ix=1, par%lv
!!      left
        u(iz+npml,npml+ix)=wave_boundary%u_bl(iz,ix,it-par%nt_shift)
        w(iz+npml,npml+ix)=wave_boundary%w_bl(iz,ix,it-par%nt_shift)
!!      right
        u(iz+npml,nx_pml-npml-par%lv+ix)=wave_boundary%u_br(iz,ix,it-par%nt_shift)
        w(iz+npml,nx_pml-npml-par%lv+ix)=wave_boundary%w_br(iz,ix,it-par%nt_shift)

!!      left
!        xx(iz+npml,npml+ix)=wave_boundary%xx_bl(iz,ix,it-par%nt_shift)
!        zz(iz+npml,npml+ix)=wave_boundary%zz_bl(iz,ix,it-par%nt_shift)
!        xz(iz+npml,npml+ix)=wave_boundary%xz_bl(iz,ix,it-par%nt_shift)
!!      right
!        xx(iz+npml,nx_pml-npml-par%lv+ix)=wave_boundary%xx_br(iz,ix,it-par%nt_shift)
!        zz(iz+npml,nx_pml-npml-par%lv+ix)=wave_boundary%zz_br(iz,ix,it-par%nt_shift)
!        xz(iz+npml,nx_pml-npml-par%lv+ix)=wave_boundary%xz_br(iz,ix,it-par%nt_shift)

     enddo
  enddo
!i,j
  do ix=1, par%nx
     do iz=1, par%lv
!!       top
        u(npml+iz,ix+npml)=wave_boundary%u_bt(ix,iz,it-par%nt_shift)
        w(npml+iz,ix+npml)=wave_boundary%w_bt(ix,iz,it-par%nt_shift)
!!       bottom
        u(nz_pml-npml-par%lv+iz,ix+npml)=wave_boundary%u_bb(ix,iz,it-par%nt_shift)
        w(nz_pml-npml-par%lv+iz,ix+npml)=wave_boundary%w_bb(ix,iz,it-par%nt_shift)

!!       top
! it works but not necessary to save stress wavefield
!        xx(npml+iz,ix+npml)=wave_boundary%xx_bt(ix,iz,it-par%nt_shift)
!        zz(npml+iz,ix+npml)=wave_boundary%zz_bt(ix,iz,it-par%nt_shift)
!        xz(npml+iz,ix+npml)=wave_boundary%xz_bt(ix,iz,it-par%nt_shift)
!!       bottom
!        xx(nz_pml-npml-par%lv+iz,ix+npml)=wave_boundary%xx_bb(ix,iz,it-par%nt_shift)
!        zz(nz_pml-npml-par%lv+iz,ix+npml)=wave_boundary%zz_bb(ix,iz,it-par%nt_shift)
!        xz(nz_pml-npml-par%lv+iz,ix+npml)=wave_boundary%xz_bb(ix,iz,it-par%nt_shift)

     enddo
  enddo
  endif

  if(it.eq.par%nt+par%nt_shift) then
      u(npml+1:nz_pml-npml,npml+1:nx_pml-npml)=wave_boundary%last_u(:,:)
      w(npml+1:nz_pml-npml,npml+1:nx_pml-npml)=wave_boundary%last_w(:,:)

      xx(npml+1:nz_pml-npml,npml+1:nx_pml-npml)=wave_boundary%last_xx(:,:)
      zz(npml+1:nz_pml-npml,npml+1:nx_pml-npml)=wave_boundary%last_zz(:,:)
      xz(npml+1:nz_pml-npml,npml+1:nx_pml-npml)=wave_boundary%last_xz(:,:)
  endif

  if (par%adj_lagran==1) then

     !! Dot product imaging condition
     !! -----------------------------
     !$omp parallel do private(iz,ix)
     do ix=npml+1,nx_pml-npml
        do iz=npml+1,nz_pml-npml
           wave_pux_fly(iz-npml,ix-npml)=c1_elastic_2th*(u(iz,ix)-u(iz,ix-1))!+c2_elastic_8th*(u(iz,ix+1)-u(iz,ix-2))+  &
                                                          !c3_elastic_8th*(u(iz,ix+2)-u(iz,ix-3))+c4_elastic_8th*(u(iz,ix+3)-u(iz,ix-4))
           wave_puz_fly(iz-npml,ix-npml)=c1_elastic_2th*(u(iz+1,ix)-u(iz,ix))!+c2_elastic_8th*(u(iz+2,ix)-u(iz-1,ix))+  &
                                                          !c3_elastic_8th*(u(iz+3,ix)-u(iz-2,ix))+c4_elastic_8th*(u(iz+4,ix)-u(iz-3,ix))
           wave_pwx_fly(iz-npml,ix-npml)=c1_elastic_2th*(w(iz,ix+1)-w(iz,ix))!+c2_elastic_8th*(w(iz,ix+2)-w(iz,ix-1))+  &
                                                          !c3_elastic_8th*(w(iz,ix+3)-w(iz,ix-2))+c4_elastic_8th*(w(iz,ix+4)-w(iz,ix-3))
           wave_pwz_fly(iz-npml,ix-npml)=c1_elastic_2th*(w(iz,ix)-w(iz-1,ix))!+c2_elastic_8th*(w(iz+1,ix)-w(iz-2,ix))+  &
                                                          !c3_elastic_8th*(w(iz+2,ix)-w(iz-3,ix))+c4_elastic_8th*(w(iz+3,ix)-w(iz-4,ix))
  !   mig(iz,ix)=mig(iz,ix)+4.0*den(iz,ix)*vs(iz,ix)*(wave_pwz(iz-pad_top,ix-npml,it)*xxb(iz,ix)+wave_pux(iz-pad_top,ix-npml,it)*zzb(iz,ix))-  &
  !                         2.0*den(iz,ix)*vs(iz,ix)*(wave_puz(iz-pad_top,ix-npml,it)+wave_pwx(iz-pad_top,ix-npml,it))*xzb(iz,ix)
           mig(iz,ix)=mig(iz,ix)+par%dt*4.0*vs(iz,ix)*(wave_pux_fly(iz-npml,ix-npml)+wave_pwz_fly(iz-npml,ix-npml))*(xxb(iz,ix)+zzb(iz,ix)) &
                     -par%dt*2.0*vs(iz,ix)*2.0*(xxb(iz,ix)*wave_pux_fly(iz-npml,ix-npml)+zzb(iz,ix)*wave_pwz_fly(iz-npml,ix-npml)+&
                      0.5*xzb(iz,ix)*(wave_puz_fly(iz-npml,ix-npml)+wave_pwx_fly(iz-npml,ix-npml)))
  !     mig(iz,ix)=mig(iz,ix)+2.0*den(iz,ix)*c(iz,ix)*wave(iz-npml,ix-npml,it)*p1(iz,ix)
        enddo
     enddo
     !$omp end parallel do

  !! Save illumination
  !! -----------------

     !$omp parallel do private(iz,ix)
     do ix=npml+1,nx_pml-npml
        do iz=npml+1,nz_pml-npml
           energy(iz,ix)=energy(iz,ix)+(wave_pux_fly(iz-npml,ix-npml)+wave_pwz_fly(iz-npml,ix-npml))**2
        enddo
     enddo
     !$omp end parallel do

  else

!      if (par%grad_twist==1) then

         !   crosco(nx,nz,npml,nx_pml,nz_pml,pad_top,txx,tzz,txz,rectxx,rectzz,rectxz,grdla,grdmu_a,grdmu_b,grdmu_c)
        call crosco(par%nx,par%nz,npml,nx_pml,nz_pml,pad_top,xxb,zzb,xzb,xx,zz,xz,grdla,grdmu_a,grdmu_b,grdmu_c)

!      endif !! endif par%grad_twist

     !! Save illumination
     !! -----------------

     !$omp parallel do private(iz,ix)
     do ix=npml+1,nx_pml-npml
        do iz=npml+1,nz_pml-npml
           energy(iz,ix)=energy(iz,ix)+(u(iz,ix)*u(iz,ix)+w(iz,ix)*w(iz,ix))
        enddo
     enddo
     !$omp end parallel do

  endif

  if(.false.) then 
 
    if(is==1 .and. mod((it-par%nt_shift),500).eq.0) then 
      write(snap_str,'(i5.5)')  (it-par%nt_shift)

      open(99999, file='snap_rec_u'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (u(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo 
      close(99999)

      open(99999, file='snap_rec_w'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (w(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo 
      close(99999)

      open(99999, file='snap_ub'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (ub(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo 
      close(99999)

      open(99999, file='snap_wb'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
      do ix=npml+1,nx_pml-npml
         write(99999,rec=ix-npml) (wb(iz,ix),iz=pad_top+1,nz_pml-npml)
      enddo 
      close(99999)

!      open(99999, file='snap_rec_xx'//snap_str, access='direct', form='unformatted', recl=4*(nz_pml-npml-pad_top))
!      do ix=npml+1,nx_pml-npml
!         write(99999,rec=ix-npml) (xx(iz,ix),iz=pad_top+1,nz_pml-npml)
!      enddo 
!      close(99999)

    endif !! if(is==10)

  endif !! if(.true.)

 endif !! save wavefield

 enddo  ! End of time-marching loop

 if (.not.(par%adj_lagran==1)) then
! if ((par%adj_lagran==1)) then

!    if (par%grad_twist==1) then

          !grconv_crosco(nx,nz,npml,nx_pml,nz_pml,pad_top,lambda,mu,grdla,grdmu,grdmu_a,grdmu_b,grdmu_c)
      call grconv_crosco(par%nx,par%nz,npml,nx_pml,nz_pml,pad_top,lamda,mu,grdla,grdmu,grdmu_a,grdmu_b,grdmu_c)
          !grconv(nx,nz,npml,nx_pml,nz_pml,pad_top,dt,rho,winvp,winvs,grdvp,grdvs,grdla,grdmu)
      call grconv(par%nx,par%nz,npml,nx_pml,nz_pml,pad_top,par%dt,den2,vp,vs,grdvp,grdvs,grdla,grdmu)

!    endif

    mig = grdvs

 endif

 !! removing potential edge effects from wavefield reconstruction 
 !! left
 ix = npml+1
 do iz=npml+1,nz_pml-npml
    mig(iz,ix) = mig(iz,ix+1)
 enddo 

 !! right 
 ix = nx_pml-npml
 do iz=npml+1,nz_pml-npml
    mig(iz,ix) = mig(iz,ix-1)
 enddo  

 !! top 
 iz = npml+1
 do ix=npml+1,nx_pml-npml
    mig(iz,ix) = mig(iz+1,ix)
 enddo 

 !! bot
 iz = nz_pml-npml
 do ix=npml+1,nx_pml-npml
    mig(iz,ix) = mig(iz-1,ix)
 enddo

 if (is.eq.1) then
 if(.False.) then
    write(is_str,'(i4.4)')  is
 
    call filename(output,trim(par%gradfile),par%iter,'.shot_mig.H@')
    call write_binfile(output,mig(npml+1:npml+par%nz,npml+1:npml+par%nx),par%nz,par%nx)

 endif
 endif

 print *, 'staggered_elastic_fwi test3', 'par%lv', par%lv


! 999 continue
 deallocate(lamda,mu,den2)
 deallocate(ub,wb,xxb,zzb,xzb)

 if(par%save_wavefield_boundary) then 
    deallocate(u,w,xx,zz,xz)
    !wave_pux_fly(:,:),wave_puz_fly(:,:),wave_pwx_fly(:,:),wave_pwz_fly(:,:)
    deallocate(wave_pux_fly,wave_puz_fly,wave_pwx_fly,wave_pwz_fly)
 endif

 deallocate(grdla,grdmu,grdmu_a,grdmu_b,grdmu_c)

 deallocate(grdvp,grdvs)

 deallocate(souresp)

end subroutine staggered_elastic_fwi_kernel


end module staggered_elastic_kernel

