
module fwi_elastic_subroutine

 use staggered_elastic_kernel
 use global
 use string
 use math
 use string 
 use io
 use datatype
 use pml
 use mmi_mpi
 use gradient_precon

 implicit none

! include 'fftw3.f'

 real, private,          allocatable  :: seis_true(:,:),seis(:,:),seis_temp1(:,:),taper(:)
 integer, private                     :: is_pre,is_pre2,it,ig

 integer, dimension(1)           :: indx_obs,indx_calc
 real                            :: distance,t_shift

 contains

subroutine compute_gradient_elastic(coord,par,is,fd_order,max_offset,s,vp,vs,den,nx_pml,nz_pml,npml,damp,energy,gk,res)

 implicit none

 type(acquisition), intent(in)   :: coord
 type(param),       intent(in)   :: par
 integer,           intent(in)   :: is,nx_pml,nz_pml,npml,fd_order,max_offset
 real,              intent(in)   :: s(:),vp(:,:),vs(:,:),den(:,:),damp(:,:)
                                    
 real,              intent(out)  :: energy(nz_pml,nx_pml),gk(nz_pml,nx_pml)
 double precision,  intent(out)  :: res

 real,allocatable :: seis_true(:,:),seis_w(:,:),seis_u(:,:)
 real,allocatable :: wave_pux(:,:,:),wave_puz(:,:,:),wave_pwx(:,:,:),wave_pwz(:,:,:)

 integer          :: kk,np,nf
 real,allocatable :: shot_virtual(:,:)

 real :: data_eps
 real :: norm_mod, norm_obs
 real :: norm_w_syn_t0, norm_true_t1
 real,allocatable :: modary_tmp(:), obsary_tmp(:), resary_tmp(:)

 real,allocatable :: seis_w_syn_t0_tmp(:), seis_true_t1_tmp(:), seis_true_t0_tmp(:)

 real, allocatable :: seis_true_t1(:,:), seis_true_t0(:,:), &
                      seis_w_syn_t0(:,:)

 type(wave_boundary_param)::wave_boundary

 real::deltat

 integer::i_start,i_end,ishift
 real::cc_max

 real,allocatable::obs_s_vel(:),mod_s_vel(:),obs_s_vel_shift(:)

 integer:: num_trace_l,num_trace_r

 real::pi,df

 real :: max_obs,max_syn

 ! logical :: save_wavefield 

 character is_str*4

 real,allocatable::obsary_tmp_t1(:),obsary_tmp_t0(:),&
                   modary_tmp_t1(:),modary_tmp_t0(:),&
                   modary_tmp_t1_cc(:),modary_tmp_t0_cc(:),&
                   modary_tmp_t1_cc_vel(:),modary_tmp_t0_cc_vel(:)

 integer::ishift_obs,ishift_mod
 real::cc_max_obs,cc_max_mod
 real::tshift_mod,tshift_obs,ddtshift_cc

 logical :: store_wave_flag

! par%lv = par%fd_order

 !! Initialize the output
 !! ---------------------
 res=0.D0;energy=0.0;gk=0.0

 ! save_wavefield = .True.

 !! Allocate and initialize the variables
 !! -------------------------------------
 allocate(seis_true(par%nt,coord%ng(is)),seis_w(par%nt,coord%ng(is)),seis_u(par%nt,coord%ng(is)))

 if (par%save_wavefield)  then 
     allocate(wave_pux(par%nz,par%nx,par%nt),wave_puz(par%nz,par%nx,par%nt),wave_pwx(par%nz,par%nx,par%nt),wave_pwz(par%nz,par%nx,par%nt))
 
     wave_pux = 0.0
     wave_puz = 0.0
     wave_pwx = 0.0
     wave_pwz = 0.0

 elseif (par%save_wavefield_boundary) then

     allocate(wave_pux(par%nz,par%nx,1),wave_puz(par%nz,par%nx,1),wave_pwx(par%nz,par%nx,1),wave_pwz(par%nz,par%nx,1))
     wave_pux = 0.0
     wave_puz = 0.0
     wave_pwx = 0.0
     wave_pwz = 0.0

     allocate(wave_boundary%u_bl(par%nz,par%lv,par%nt), wave_boundary%u_br(par%nz,par%lv,par%nt), &
              wave_boundary%u_bt(par%nx,par%lv,par%nt), wave_boundary%u_bb(par%nx,par%lv,par%nt))

     allocate(wave_boundary%w_bl(par%nz,par%lv,par%nt), wave_boundary%w_br(par%nz,par%lv,par%nt), &
              wave_boundary%w_bt(par%nx,par%lv,par%nt), wave_boundary%w_bb(par%nx,par%lv,par%nt))

     wave_boundary%u_bl = 0.0; wave_boundary%u_br = 0.0; wave_boundary%u_bt = 0.0; wave_boundary%u_bb = 0.0
     wave_boundary%w_bl = 0.0; wave_boundary%w_br = 0.0; wave_boundary%w_bt = 0.0; wave_boundary%w_bb = 0.0
     
     allocate(wave_boundary%last_u(par%nz,par%nx), wave_boundary%last_w(par%nz,par%nx))

     wave_boundary%last_u = 0.
     wave_boundary%last_w = 0.

   !!! this is for test
   !  allocate(wave_boundary%xx_bl(par%nz,par%lv,par%nt), wave_boundary%xx_br(par%nz,par%lv,par%nt), &
   !           wave_boundary%xx_bt(par%nx,par%lv,par%nt), wave_boundary%xx_bb(par%nx,par%lv,par%nt))      

   !  allocate(wave_boundary%zz_bl(par%nz,par%lv,par%nt), wave_boundary%zz_br(par%nz,par%lv,par%nt), &
   !           wave_boundary%zz_bt(par%nx,par%lv,par%nt), wave_boundary%zz_bb(par%nx,par%lv,par%nt))      

   !  allocate(wave_boundary%xz_bl(par%nz,par%lv,par%nt), wave_boundary%xz_br(par%nz,par%lv,par%nt), &
   !           wave_boundary%xz_bt(par%nx,par%lv,par%nt), wave_boundary%xz_bb(par%nx,par%lv,par%nt))          

   !  wave_boundary%xx_bl = 0.0; wave_boundary%xx_br = 0.0; wave_boundary%xx_bt = 0.0; wave_boundary%xx_bb = 0.0
   !  wave_boundary%zz_bl = 0.0; wave_boundary%zz_br = 0.0; wave_boundary%zz_bt = 0.0; wave_boundary%zz_bb = 0.0
   !  wave_boundary%xz_bl = 0.0; wave_boundary%xz_br = 0.0; wave_boundary%xz_bt = 0.0; wave_boundary%xz_bb = 0.0

     allocate(wave_boundary%last_xx(par%nz,par%nx), &
              wave_boundary%last_zz(par%nz,par%nx), &
              wave_boundary%last_xz(par%nz,par%nx))

     wave_boundary%last_xx = 0.0
     wave_boundary%last_zz = 0.0
     wave_boundary%last_xz = 0.0

 endif

 if (.not.allocated(wave_pux)) print *, 'wave_pux not allocated'
 if (.not.allocated(wave_puz)) print *, 'wave_puz not allocated'
 if (.not.allocated(wave_pwx)) print *, 'wave_pwx not allocated'
 if (.not.allocated(wave_pwz)) print *, 'wave_pwz not allocated'

 seis_true=0.0;  seis=0.0
 seis_w = 0.0
 seis_u = 0.0

 print *, 'compute_gradient_elastic ', 'test0'

 is_pre=999
 !! Read observed data and generate the synthetic data
 !! ------------------------------------------------- 

!  if (save_wavefield)  then
 store_wave_flag = .true.
 call modeling_elastic_kernel(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
      wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_u,seis_w,store_wave_flag)

! elseif (save_wavefield_boundary) then
!
!     call elastic_modeling_bvr(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,save_wavefield_boundary,wave_boundary,seis_u,seis_w)

! endif

 print *, 'compute_gradient_elastic ', 'test1'

 call filename(output,par%csg_in,is,'.wH@')

 print *, 'output in compute_gradient_elastic ', output
 call read_binfile(output,seis_true,par%nt,coord%ng(is))

 print *, 'time_lapse', par%time_lapse

 if (par%time_lapse) then 

   print *, 'now begin time lapse'

   allocate(seis_true_t1(par%nt,coord%ng(is)), &
            seis_true_t0(par%nt,coord%ng(is)), &
            seis_w_syn_t0(par%nt,coord%ng(is)))

   seis_true_t1 = 0.0
   seis_true_t0 = 0.0
   seis_w_syn_t0 = 0.0

   allocate(seis_w_syn_t0_tmp(1:par%nt), seis_true_t1_tmp(1:par%nt))
   allocate(seis_true_t0_tmp(1:par%nt))

   seis_w_syn_t0_tmp = 0.0
   seis_true_t1_tmp = 0.0
   seis_true_t0_tmp = 0.0

   call filename(output,par%csg_t1_in,is,'.wH@')
   call read_binfile(output,seis_true_t1,par%nt,coord%ng(is))

   call filename(output,par%csg_syn_t0_in,is,'.wH@')
   call read_binfile(output,seis_w_syn_t0,par%nt,coord%ng(is))

 endif 

 allocate(shot_virtual(par%nt,coord%ng(is)))
 shot_virtual=0.0

 max_obs=maxval(seis_true)
 !! here I only use the vertical component
 max_syn=maxval(seis_w)

 print *, 'max_obs', max_obs, 'max_syn', max_syn

 print *, 'obtype', par%obtype 

 if(.True.) then
    call filename(output,par%csg_out,is,'.wH@')
    call write_binfile(output,seis_w,par%nt,coord%ng(is))
 endif

 if (par%time_lapse) then 
!   seis_true_t1 = seis_true_t1/ maxval(seis_true_t1) 
!   seis_w_syn_t0 = seis_w_syn_t0/ maxval(seis_w_syn_t0) 
 endif

 do ig=1,coord%ng(is)
 do it=1,par%nt
!    seis_true(it,ig)=seis_true(it,ig)/max_obs
!    seis_w(it,ig)=seis_w(it,ig)/max_syn
 enddo
 enddo

!  ## begin section 1:
! if (par%time_lapse) then
!    seis_true_t0 = seis_true
!    seis_true = 0.0
!    seis_true = seis_w_syn_t0 + (seis_true_t1 - seis_true_t0)
! endif
! ## end section 1

 data_eps = 1.0e-20

 allocate(modary_tmp(1:par%nt))
 allocate(obsary_tmp(1:par%nt))
 allocate(resary_tmp(1:par%nt))

 modary_tmp = 0.0
 obsary_tmp = 0.0
 resary_tmp = 0.0

! if (par%amp_norm.eq.1) then
! endif

 if (par%obtype.eq.1) then

   if (par%amp_norm_trace.eq.1) then
       do ig=1,coord%ng(is)
          seis_true(:,ig) = seis_true(:,ig)/(maxval(seis_true(:,ig))+data_eps)
          seis_w(:,ig) = seis_w(:,ig)/(maxval(seis_w(:,ig))+data_eps)

          if (par%time_lapse) then
             seis_w_syn_t0(:,ig) = seis_w_syn_t0(:,ig)/(maxval(seis_w_syn_t0(:,ig))+data_eps)
             seis_true_t1(:,ig) = seis_true_t1(:,ig)/(maxval(seis_true_t1(:,ig))+data_eps)
          endif
       enddo 
   endif

   if (par%amp_norm_csg.eq.1) then
       do ig=1,coord%ng(is)
          seis_true(:,ig) = seis_true(:,ig)/(maxval(seis_true)+data_eps)
          seis_w(:,ig) = seis_w(:,ig)/(maxval(seis_w)+data_eps)

          if (par%time_lapse) then
             seis_w_syn_t0(:,ig) = seis_w_syn_t0(:,ig)/(maxval(seis_w_syn_t0)+data_eps)
             seis_true_t1(:,ig) = seis_true_t1(:,ig)/(maxval(seis_true_t1)+data_eps)
          endif
       enddo
   endif

   if (par%time_lapse) then
      seis_true_t0 = seis_true
      seis_true = seis_w_syn_t0 + (seis_true_t1 - seis_true_t0)
   endif

   shot_virtual = seis_true - seis_w 


!! this is trace-normalized misfit function
 elseif (par%obtype.eq.2) then 

   do ig=1,coord%ng(is)

       norm_mod = sqrt( sum((seis_w(:,ig) * seis_w(:,ig))) )
       norm_obs = sqrt( sum((seis_true(:,ig) * seis_true(:,ig))) )

       modary_tmp(:) = seis_w(:,ig)/(norm_mod + data_eps)
       obsary_tmp(:) = seis_true(:,ig)/(norm_obs + data_eps)

       if (par%time_lapse) then
          seis_true_t0_tmp(:) = obsary_tmp(:)

          norm_w_syn_t0 = sqrt( sum((seis_w_syn_t0(:,ig) * seis_w_syn_t0(:,ig))) )
          norm_true_t1 = sqrt( sum((seis_true_t1(:,ig) * seis_true_t1(:,ig))) )
          seis_w_syn_t0_tmp(:) = seis_w_syn_t0(:,ig)/(norm_w_syn_t0+data_eps)
          seis_true_t1_tmp(:) = seis_true_t1(:,ig)/(norm_true_t1+data_eps)

          obsary_tmp(:) = seis_w_syn_t0_tmp(:) + (seis_true_t1_tmp(:) - seis_true_t0_tmp(:))
       endif

       resary_tmp(:) = modary_tmp(:) - obsary_tmp(:)

       ! shot_virtual(:,ig) = resary_tmp(:)/(norm_mod+data_eps) + resary_tmp(:)*modary_tmp(:)/(norm_mod+data_eps)**2*seis_w(:,ig)
       shot_virtual(:,ig) = resary_tmp(:)/(norm_mod+data_eps) - modary_tmp(:)*sum(resary_tmp(:)*seis_w(:,ig))/((norm_mod+data_eps)**2)
   
       shot_virtual(:,ig) = - shot_virtual(:,ig)

       res = res + sum(resary_tmp(:)*resary_tmp(:))
   enddo 

 elseif (par%obtype.eq.3) then

   do ig=1,coord%ng(is)

       norm_mod = sqrt( sum((seis_w(:,ig) * seis_w(:,ig))) )
       norm_obs = sqrt( sum((seis_true(:,ig) * seis_true(:,ig))) )

       modary_tmp(:) = seis_w(:,ig)/(norm_mod + data_eps)
       obsary_tmp(:) = seis_true(:,ig)/(norm_obs + data_eps)

       if (par%time_lapse) then
          seis_true_t0_tmp(:) = obsary_tmp(:)

          norm_w_syn_t0 = sqrt( sum((seis_w_syn_t0(:,ig) * seis_w_syn_t0(:,ig))) )
          norm_true_t1 = sqrt( sum((seis_true_t1(:,ig) * seis_true_t1(:,ig))) )
          seis_w_syn_t0_tmp(:) = seis_w_syn_t0(:,ig)/(norm_w_syn_t0+data_eps)
          seis_true_t1_tmp(:) = seis_true_t1(:,ig)/(norm_true_t1+data_eps)

          obsary_tmp(:) = seis_w_syn_t0_tmp(:) + (seis_true_t1_tmp(:) - seis_true_t0_tmp(:))
       endif

       shot_virtual(:,ig) = modary_tmp(:) - obsary_tmp(:)

       shot_virtual(:,ig) = - shot_virtual(:,ig)

       res = res + sum((shot_virtual(:,ig)*seis_w(:,ig)))
   enddo

 elseif (par%obtype.eq.4) then  

    deltat=par%dt

    allocate(obs_s_vel(1:par%nt))
    allocate(mod_s_vel(1:par%nt))
    allocate(obs_s_vel_shift(1:par%nt))

    i_start = 1
    i_end = par%nt

    do ig=1,coord%ng(is)
          
        obsary_tmp(:) = seis_true(:,ig)
        modary_tmp(:) = seis_w(:,ig)

        call xcorr_calc(obsary_tmp,modary_tmp,par%nt,i_start,i_end,par%nlen_shift,ishift,cc_max)

        do it = 2, par%nt-1
           obs_s_vel(it) = (obsary_tmp(it+1)-obsary_tmp(it-1))/(2*deltat)
           mod_s_vel(it) = (modary_tmp(it+1)-modary_tmp(it-1))/(2*deltat)
        enddo

        obs_s_vel(1) = (obsary_tmp(2)-obsary_tmp(1))/deltat
        obs_s_vel(par%nt) = (obsary_tmp(par%nt)-obsary_tmp(par%nt-1))/deltat

        mod_s_vel(1) = (modary_tmp(2)-modary_tmp(1))/deltat
        mod_s_vel(par%nt) = (modary_tmp(par%nt)-modary_tmp(par%nt-1))/deltat

        obs_s_vel_shift = cshift(obs_s_vel, shift = ishift)

        !   norm = sum(obs_s_vel_shift * mod_s_vel) * deltat + data_eps

        shot_virtual(:,ig) = - 2.0 * obs_s_vel_shift(:) * (ishift * deltat)

        res = res + (- ishift * deltat)**2

    enddo 

    deallocate(obs_s_vel,mod_s_vel,obs_s_vel_shift)

 elseif (par%obtype.eq.5) then

    deltat=par%dt

    allocate(obs_s_vel(1:par%nt))
    allocate(mod_s_vel(1:par%nt))
    allocate(obs_s_vel_shift(1:par%nt))

    i_start = 1
    i_end = par%nt

    if (par%time_lapse==0) then 

     do ig=1,coord%ng(is)

       obsary_tmp(:) = seis_true(:,ig)
       modary_tmp(:) = seis_w(:,ig)

       call xcorr_calc_mod(modary_tmp,obsary_tmp,par%nt,i_start,i_end,par%nlen_shift,ishift,cc_max)

       do it = 2, par%nt-1
          obs_s_vel(it) = (obsary_tmp(it+1)-obsary_tmp(it-1))/(2*deltat)
          mod_s_vel(it) = (modary_tmp(it+1)-modary_tmp(it-1))/(2*deltat)
       enddo

       obs_s_vel(1) = (obsary_tmp(2)-obsary_tmp(1))/deltat
       obs_s_vel(par%nt) = (obsary_tmp(par%nt)-obsary_tmp(par%nt-1))/deltat

       mod_s_vel(1) = (modary_tmp(2)-modary_tmp(1))/deltat
       mod_s_vel(par%nt) = (modary_tmp(par%nt)-modary_tmp(par%nt-1))/deltat

       obs_s_vel_shift = cshift(obs_s_vel, shift = ishift)

  !  norm = sum(mod_s_vel * mod_s_vel) * deltat + 1.0d-32

       shot_virtual(:,ig) = mod_s_vel(:) * (ishift * deltat)

       res = res + (ishift * deltat)**2
     enddo 

    elseif (par%time_lapse==1) then

     allocate(obsary_tmp_t1(1:par%nt))
     allocate(obsary_tmp_t0(1:par%nt))
     allocate(modary_tmp_t1(1:par%nt))
     allocate(modary_tmp_t0(1:par%nt))

     allocate(modary_tmp_t1_cc(1:par%nt),modary_tmp_t0_cc(1:par%nt))
     allocate(modary_tmp_t1_cc_vel(1:par%nt),modary_tmp_t0_cc_vel(1:par%nt))

     obsary_tmp_t1 = 0.
     obsary_tmp_t0 = 0.
     modary_tmp_t1 = 0.
     modary_tmp_t0 = 0.
     modary_tmp_t1_cc = 0.
     modary_tmp_t0_cc = 0.
     modary_tmp_t1_cc_vel = 0.
     modary_tmp_t0_cc_vel = 0.

     do ig=1,coord%ng(is)

       !!! for time-lapse part
       obsary_tmp_t1(:) = seis_true_t1(:,ig)
       obsary_tmp_t0(:) = seis_true(:,ig)

       modary_tmp_t1(:) = seis_w(:,ig)
       modary_tmp_t0(:) = seis_w_syn_t0(:,ig)

       call xcorr_calc_mod(obsary_tmp_t0,obsary_tmp_t1,par%nt,i_start,i_end,par%nlen_shift,ishift_obs,cc_max_obs) ! T(d1-d2)
       call xcorr_calc_mod(modary_tmp_t0,modary_tmp_t1,par%nt,i_start,i_end,par%nlen_shift,ishift_mod,cc_max_mod) ! T(d1-d2)

       tshift_mod = ishift_mod * deltat
       tshift_obs = ishift_obs * deltat

   !    print *, 'inside gradient calc', ig, is, ishift_mod, ishift_obs
   !    print *, 'inside gradient calc', ig, is, tshift_mod, tshift_obs

       !! double-difference cc-measurement 
       ddtshift_cc = tshift_mod - tshift_obs

       res = res + (ddtshift_cc)**2

       modary_tmp_t1_cc = cshift(modary_tmp_t1, shift = ishift_mod)
       modary_tmp_t0_cc = cshift(modary_tmp_t0, shift = -ishift_mod)
     
       do it = 2, par%nt-1
          modary_tmp_t1_cc_vel(it) = (modary_tmp_t1_cc(it+1)-modary_tmp_t1_cc(it-1))/(2*deltat)
          modary_tmp_t0_cc_vel(it) = (modary_tmp_t0_cc(it+1)-modary_tmp_t0_cc(it-1))/(2*deltat)
       enddo

       modary_tmp_t1_cc_vel(1) = (modary_tmp_t1_cc(2)-modary_tmp_t1_cc(1))/deltat
       modary_tmp_t1_cc_vel(par%nt) = (modary_tmp_t1_cc(par%nt)-modary_tmp_t1_cc(par%nt-1))/deltat

       modary_tmp_t0_cc_vel(1) = (modary_tmp_t0_cc(2)-modary_tmp_t0_cc(1))/deltat
       modary_tmp_t0_cc_vel(par%nt) = (modary_tmp_t0_cc(par%nt)-modary_tmp_t0_cc(par%nt-1))/deltat

       shot_virtual(:,ig)= + ddtshift_cc * modary_tmp_t1_cc_vel(:)

       shot_virtual(:,ig)= - shot_virtual(:,ig)

     enddo 

     deallocate(obsary_tmp_t1,obsary_tmp_t0)
     deallocate(modary_tmp_t1,modary_tmp_t0)
     deallocate(modary_tmp_t1_cc, modary_tmp_t0_cc)
     deallocate(modary_tmp_t1_cc_vel,modary_tmp_t0_cc_vel)

    endif !!(par%time_lapse==1)

    deallocate(obs_s_vel,mod_s_vel,obs_s_vel_shift)


 endif

   !! (1.6) compute residual
   !! ----------------------
 if (par%obtype.eq.1) then

  do ig=1,coord%ng(is)
     do it=1,par%nt
      res=res+shot_virtual(it,ig)*shot_virtual(it,ig)
     enddo 
  enddo

 endif


 if(.True.) then
    call filename(output,par%csg_out,is,'.post.wH@')
    call write_binfile(output,seis_w,par%nt,coord%ng(is))
 endif

 !------------------------------------------------------------------------------------------

! res=sqrt(res)

 print *, 'res in compute_gradient_elastic ', res, maxval(seis_w), maxval(shot_virtual), coord%ng(is)

! if(is==5) then
 if(.False.) then
    write(is_str,'(i4.4)')  is 
    open(88888,file='shot_virtual.bin'//is_str,access='direct', form='unformatted', recl=4*par%nt)
    do ig=1,coord%ng(is)
       write(88888, rec=ig) (shot_virtual(it,ig), it=1,par%nt)
    enddo 
    close(88888)
 endif

! if(par%mute_direct.eq.1) then
!    call mute_direct(par,coord,is,c,taper,seis)
!    call mute_direct(par,coord,is,c,taper,seis_true)
! endif

 is_pre = is

 print *, 'is_pre in compute_gradient_elastic ', is_pre

 print *, 'c1_elastic_2th', c1_elastic_2th, 'par%lv', par%lv

 ! if(is.ne.is_pre) then
 !    call elastic_modeling(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,.True.,wave_pux,wave_puz,wave_pwx,wave_pwz,seis_u,seis_w)
 !    is_pre=is
 !endif

 !! Compute gradient
 !! ----------------

 ! if (save_wavefield)  then
!      staggered_elastic_fwi_kernel(is,par,coord,s,vp,vs,den,fs,fd_order,nx_pml,nz_pml,npml,damp,&
!      wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_v,energy,mig)

! if (par%adj_lagran==1) then 

  call staggered_elastic_fwi_kernel(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
       wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,shot_virtual,energy,gk)

! else 
!
!  call staggered_elastic_rtm_conven(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
!       wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,shot_virtual,energy,gk)       
!
! endif        

 !elseif (save_wavefield_boundary) then 

 !   call staggered_elastic_rtm_bvr(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,wave_boundary,shot_virtual,energy,gk) 

 !endif

 print *, 'compute_gradient_elastic ', 'test3'

 deallocate(seis_true,seis_w,seis_u)

 if (par%save_wavefield)  then

    deallocate(wave_pux,wave_puz,wave_pwx,wave_pwz)

 elseif (par%save_wavefield_boundary) then

    deallocate(wave_pux,wave_puz,wave_pwx,wave_pwz)

    deallocate(wave_boundary%u_bl, wave_boundary%u_br, wave_boundary%u_bt, wave_boundary%u_bb)
    deallocate(wave_boundary%w_bl, wave_boundary%w_br, wave_boundary%w_bt, wave_boundary%w_bb)

    deallocate(wave_boundary%last_u, wave_boundary%last_w)

!    deallocate(wave_boundary%xx_bl, wave_boundary%xx_br, wave_boundary%xx_bt, wave_boundary%xx_bb)
!    deallocate(wave_boundary%zz_bl, wave_boundary%zz_br, wave_boundary%zz_bt, wave_boundary%zz_bb)
!    deallocate(wave_boundary%xz_bl, wave_boundary%xz_br, wave_boundary%xz_bt, wave_boundary%xz_bb)

    deallocate(wave_boundary%last_xx, wave_boundary%last_zz, wave_boundary%last_xz)
 endif

 deallocate(shot_virtual)

 deallocate(modary_tmp,obsary_tmp,resary_tmp)

 if(par%time_lapse) then 
     deallocate(seis_true_t1, seis_true_t0, &
                      seis_w_syn_t0)
     deallocate(seis_w_syn_t0_tmp, seis_true_t1_tmp, &
                seis_true_t0_tmp)
 endif
              
end subroutine compute_gradient_elastic


subroutine analytical_step_length_elastic(coord,par,is,fd_order,max_offset,s,vp,vs,den,fs,nx_pml,nz_pml,npml,damp,read_data,dgg)
 implicit none

 type(acquisition), intent(in)   :: coord
 type(param),       intent(inout):: par
 integer,           intent(in)   :: is,nx_pml,nz_pml,npml,fs,fd_order,max_offset
 real,              intent(in)   :: s(:),vp(:,:),vs(:,:),den(:,:),damp(:,:)
 logical,           intent(in)   :: read_data
 double precision,  intent(out)  :: dgg

 real,allocatable :: wave_puz(:,:,:),wave_pux(:,:,:),wave_pwz(:,:,:),wave_pwx(:,:,:)
 real,allocatable :: seis_w(:,:),seis_u(:,:),seis_true(:,:)

 integer          :: kk,np,nf

 real,allocatable :: shot_res_u(:,:), shot_res_w(:,:)

 type(wave_boundary_param)::wave_boundary

 real :: max_obs,max_syn

 real :: data_eps
 real :: norm_mod, norm_obs
 real :: norm_w_syn_t0, norm_true_t1 

 real,allocatable :: modary_tmp(:), obsary_tmp(:), resary_tmp(:)

 real,allocatable :: seis_w_syn_t0_tmp(:), seis_true_t1_tmp(:), seis_true_t0_tmp(:)

 real, allocatable :: seis_true_t1(:,:), seis_true_t0(:,:), &
                      seis_w_syn_t0(:,:)

! logical :: save_wavefield, save_wavefield_boundary

 real::deltat

 integer::i_start,i_end,ishift
 real::cc_max

 real,allocatable::obs_s_vel(:),mod_s_vel(:),obs_s_vel_shift(:)

 real,allocatable::obsary_tmp_t1(:),obsary_tmp_t0(:),&
                   modary_tmp_t1(:),modary_tmp_t0(:),&
                   modary_tmp_t1_cc(:),modary_tmp_t0_cc(:),&
                   modary_tmp_t1_cc_vel(:),modary_tmp_t0_cc_vel(:)

 integer::ishift_obs,ishift_mod
 real::cc_max_obs,cc_max_mod
 real::tshift_mod,tshift_obs,ddtshift_cc

 logical::store_wave_flag

 ! save_wavefield = .false.
 !save_wavefield _boundary= .false.

 !! Initialize the output
 !! ---------------------
 dgg=0.D0

 !! Allocate and initialize the variables
 !! -------------------------------------
 allocate(seis_w(par%nt,coord%ng(is)),seis_u(par%nt,coord%ng(is)),seis_true(par%nt,coord%ng(is))); 
 seis_w=0.0; seis_u=0.0

! if (.not.par%save_wavefield) then 

 allocate(wave_pux(par%nz,par%nx,1),wave_puz(par%nz,par%nx,1),&
          wave_pwx(par%nz,par%nx,1),wave_pwz(par%nz,par%nx,1))

! endif


 !! Read observed data and generate the synthetic data
 !! -------------------------------------------------
 store_wave_flag = .false.
 call modeling_elastic_kernel(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
      wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_u,seis_w,store_wave_flag)

 call filename(output,par%csg_in,is,'.wH@')
 call read_binfile(output,seis_true,par%nt,coord%ng(is))
 !if(par%mute_direct.eq.1) call mute_direct(par,coord,is,v,taper,seis)

! if(par%time_shift_steps .ne. 0) then
!     ! positive number: shift to left (earlier time) V:(1, 2, 3, 4, 5, 6). after shift:(3, 4, 5, 6, 1, 2)
!     ! negative number: shift to right (later time)  V:(1, 2, 3, 4, 5, 6). after shift:(5, 6, 1, 2, 3, 4)
!    seis_true = cshift(seis_true, shift = par%time_shift_steps, dim=2)
!
! endif

 if (par%time_lapse) then

   print *, 'now begin time lapse'

   allocate(seis_true_t1(par%nt,coord%ng(is)), &
            seis_true_t0(par%nt,coord%ng(is)), &
            seis_w_syn_t0(par%nt,coord%ng(is)))

   seis_true_t1 = 0.0
   seis_true_t0 = 0.0
   seis_w_syn_t0 = 0.0

   allocate(seis_w_syn_t0_tmp(1:par%nt), seis_true_t1_tmp(1:par%nt))
   allocate(seis_true_t0_tmp(1:par%nt))

   seis_w_syn_t0_tmp = 0.0
   seis_true_t1_tmp = 0.0
   seis_true_t0_tmp = 0.0

   call filename(output,par%csg_t1_in,is,'.wH@')
   call read_binfile(output,seis_true_t1,par%nt,coord%ng(is))

!   if(par%time_shift_steps .ne. 0) then
!     ! positive number: shift to left (earlier time) V:(1, 2, 3, 4, 5, 6). after shift:(3, 4, 5, 6, 1, 2)
!     ! negative number: shift to right (later time)  V:(1, 2, 3, 4, 5, 6). after shift:(5, 6, 1, 2, 3, 4)
!     seis_true_t1 = cshift(seis_true_t1, shift = par%time_shift_steps, dim=2)
!
!   endif

   call filename(output,par%csg_syn_t0_in,is,'.wH@')
   call read_binfile(output,seis_w_syn_t0,par%nt,coord%ng(is))

 endif


 allocate(shot_res_u(par%nt,coord%ng(is)))
 shot_res_u = 0.0

 allocate(shot_res_w(par%nt,coord%ng(is)))
 shot_res_w = 0.0

 !------------------------------------------------------------------------------------------
  max_obs=maxval(seis_true)
  max_syn=maxval(seis_w)

  do ig=1,coord%ng(is)
  do it=1,par%nt
!       seis_true(it,ig)=seis_true(it,ig)/max_obs
!       seis_w(it,ig)=seis_w(it,ig)/max_syn
  enddo
  enddo

  if (par%time_lapse) then
!   seis_true_t1 = seis_true_t1/ maxval(seis_true_t1)
!   seis_w_syn_t0 = seis_w_syn_t0/ maxval(seis_w_syn_t0)
  endif

!  ## begin section 1: 
!  if (par%time_lapse) then
!     seis_true_t0 = seis_true
!     seis_true = 0.0
!     seis_true = seis_w_syn_t0 + (seis_true_t1 - seis_true_t0)
!  endif
! ## end section 1

  data_eps = 1.0e-20

  allocate(modary_tmp(1:par%nt))
  allocate(obsary_tmp(1:par%nt))
  allocate(resary_tmp(1:par%nt))

  modary_tmp = 0.0
  obsary_tmp = 0.0
  resary_tmp = 0.0

  dgg=0.D0
  if (par%obtype.eq.1) then

     if (par%amp_norm_trace.eq.1) then
         do ig=1,coord%ng(is)
            seis_true(:,ig) = seis_true(:,ig)/(maxval(seis_true(:,ig))+data_eps)
            seis_w(:,ig) = seis_w(:,ig)/(maxval(seis_w(:,ig))+data_eps)

            if (par%time_lapse) then
               seis_w_syn_t0(:,ig) = seis_w_syn_t0(:,ig)/(maxval(seis_w_syn_t0(:,ig))+data_eps)
               seis_true_t1(:,ig) = seis_true_t1(:,ig)/(maxval(seis_true_t1(:,ig))+data_eps)
            endif
         enddo
     endif

     if (par%amp_norm_csg.eq.1) then
         do ig=1,coord%ng(is)
            seis_true(:,ig) = seis_true(:,ig)/(maxval(seis_true)+data_eps)
            seis_w(:,ig) = seis_w(:,ig)/(maxval(seis_w)+data_eps)

            if (par%time_lapse) then
               seis_w_syn_t0(:,ig) = seis_w_syn_t0(:,ig)/(maxval(seis_w_syn_t0)+data_eps)
               seis_true_t1(:,ig) = seis_true_t1(:,ig)/(maxval(seis_true_t1)+data_eps)
            endif
         enddo
     endif

     if (par%time_lapse) then
         seis_true_t0 = seis_true
         seis_true = seis_w_syn_t0 + (seis_true_t1 - seis_true_t0)
     endif
     shot_res_w = seis_true - seis_w

  elseif (par%obtype.eq.2) then

     do ig=1,coord%ng(is)

        norm_mod = sqrt( sum((seis_w(:,ig) * seis_w(:,ig))) )
        norm_obs = sqrt( sum((seis_true(:,ig) * seis_true(:,ig))) )

        modary_tmp(:) = seis_w(:,ig)/(norm_mod + data_eps)
        obsary_tmp(:) = seis_true(:,ig)/(norm_obs + data_eps)

        if (par%time_lapse) then
            seis_true_t0_tmp(:) = obsary_tmp(:)

            norm_w_syn_t0 = sqrt( sum((seis_w_syn_t0(:,ig) * seis_w_syn_t0(:,ig))) )
            norm_true_t1 = sqrt( sum((seis_true_t1(:,ig) * seis_true_t1(:,ig))) )
            seis_w_syn_t0_tmp(:) = seis_w_syn_t0(:,ig)/(norm_w_syn_t0+data_eps)
            seis_true_t1_tmp(:) = seis_true_t1(:,ig)/(norm_true_t1+data_eps)

            obsary_tmp(:) = seis_w_syn_t0_tmp(:) + (seis_true_t1_tmp(:) - seis_true_t0_tmp(:))
        endif

        resary_tmp(:) = modary_tmp(:) - obsary_tmp(:)

        shot_res_w(:,ig) = resary_tmp(:)/(norm_mod+data_eps) + resary_tmp(:)*modary_tmp(:)/(norm_mod+data_eps)**2*seis_w(:,ig)

        shot_res_w(:,ig) = - shot_res_w(:,ig)

        dgg = dgg + sum(resary_tmp(:)*resary_tmp(:))

     enddo          

  elseif (par%obtype.eq.3) then

     do ig=1,coord%ng(is)

        norm_mod = sqrt( sum((seis_w(:,ig) * seis_w(:,ig))) )
        norm_obs = sqrt( sum((seis_true(:,ig) * seis_true(:,ig))) )

        modary_tmp(:) = seis_w(:,ig)/(norm_mod + data_eps)
        obsary_tmp(:) = seis_true(:,ig)/(norm_obs + data_eps)
        
        if (par%time_lapse) then
            seis_true_t0_tmp(:) = obsary_tmp(:)

            norm_w_syn_t0 = sqrt( sum((seis_w_syn_t0(:,ig) * seis_w_syn_t0(:,ig))) )
            norm_true_t1 = sqrt( sum((seis_true_t1(:,ig) * seis_true_t1(:,ig))) )
            seis_w_syn_t0_tmp(:) = seis_w_syn_t0(:,ig)/(norm_w_syn_t0+data_eps)
            seis_true_t1_tmp(:) = seis_true_t1(:,ig)/(norm_true_t1+data_eps)

            obsary_tmp(:) = seis_w_syn_t0_tmp(:) + (seis_true_t1_tmp(:) - seis_true_t0_tmp(:))
        endif

        shot_res_w(:,ig) = modary_tmp(:) - obsary_tmp(:)

        shot_res_w(:,ig) = - shot_res_w(:,ig)

        dgg = dgg + sum((shot_res_w(:,ig)*seis_w(:,ig)))

     enddo

  elseif (par%obtype.eq.4) then

     deltat=par%dt

     allocate(obs_s_vel(1:par%nt))
     allocate(mod_s_vel(1:par%nt))
     allocate(obs_s_vel_shift(1:par%nt))

     i_start = 1
     i_end = par%nt

     do ig=1,coord%ng(is)

        obsary_tmp(:) = seis_true(:,ig)
        modary_tmp(:) = seis_w(:,ig)

        call xcorr_calc(obsary_tmp,modary_tmp,par%nt,i_start,i_end,par%nlen_shift,ishift,cc_max)

        do it = 2, par%nt-1
           obs_s_vel(it) = (obsary_tmp(it+1)-obsary_tmp(it-1))/(2*deltat)
           mod_s_vel(it) = (modary_tmp(it+1)-modary_tmp(it-1))/(2*deltat)
        enddo

        obs_s_vel(1) = (obsary_tmp(2)-obsary_tmp(1))/deltat
        obs_s_vel(par%nt) = (obsary_tmp(par%nt)-obsary_tmp(par%nt-1))/deltat

        mod_s_vel(1) = (modary_tmp(2)-modary_tmp(1))/deltat
        mod_s_vel(par%nt) = (modary_tmp(par%nt)-modary_tmp(par%nt-1))/deltat

        obs_s_vel_shift = cshift(obs_s_vel, shift = ishift)

        !   norm = sum(obs_s_vel_shift * mod_s_vel) * deltat + data_eps

        dgg = dgg + (- ishift * deltat)**2

     enddo

     deallocate(obs_s_vel,mod_s_vel,obs_s_vel_shift)

  elseif (par%obtype.eq.5) then

     deltat=par%dt

     allocate(obs_s_vel(1:par%nt))
     allocate(mod_s_vel(1:par%nt))
     allocate(obs_s_vel_shift(1:par%nt))

     i_start = 1
     i_end = par%nt

     if (par%time_lapse==0) then

      do ig=1,coord%ng(is)

        obsary_tmp(:) = seis_true(:,ig)
        modary_tmp(:) = seis_w(:,ig)

        call xcorr_calc_mod(modary_tmp,obsary_tmp,par%nt,i_start,i_end,par%nlen_shift,ishift,cc_max)

        do it = 2, par%nt-1
           obs_s_vel(it) = (obsary_tmp(it+1)-obsary_tmp(it-1))/(2*deltat)
           mod_s_vel(it) = (modary_tmp(it+1)-modary_tmp(it-1))/(2*deltat)
        enddo

        obs_s_vel(1) = (obsary_tmp(2)-obsary_tmp(1))/deltat
        obs_s_vel(par%nt) = (obsary_tmp(par%nt)-obsary_tmp(par%nt-1))/deltat

        mod_s_vel(1) = (modary_tmp(2)-modary_tmp(1))/deltat
        mod_s_vel(par%nt) = (modary_tmp(par%nt)-modary_tmp(par%nt-1))/deltat

        obs_s_vel_shift = cshift(obs_s_vel, shift = ishift)

        dgg = dgg + (ishift * deltat)**2

      enddo

     elseif (par%time_lapse==1) then

      allocate(obsary_tmp_t1(1:par%nt))
      allocate(obsary_tmp_t0(1:par%nt))
      allocate(modary_tmp_t1(1:par%nt))
      allocate(modary_tmp_t0(1:par%nt))

      allocate(modary_tmp_t1_cc(1:par%nt),modary_tmp_t0_cc(1:par%nt))
      allocate(modary_tmp_t1_cc_vel(1:par%nt),modary_tmp_t0_cc_vel(1:par%nt))

      obsary_tmp_t1 = 0.
      obsary_tmp_t0 = 0.
      modary_tmp_t1 = 0.
      modary_tmp_t0 = 0.
      modary_tmp_t1_cc = 0.
      modary_tmp_t0_cc = 0.
      modary_tmp_t1_cc_vel = 0.
      modary_tmp_t0_cc_vel = 0.

      do ig=1,coord%ng(is)

       !!! for time-lapse part
       obsary_tmp_t1(:) = seis_true_t1(:,ig)
       obsary_tmp_t0(:) = seis_true(:,ig)

       modary_tmp_t1(:) = seis_w(:,ig)
       modary_tmp_t0(:) = seis_w_syn_t0(:,ig)

       call xcorr_calc_mod(obsary_tmp_t0,obsary_tmp_t1,par%nt,i_start,i_end,par%nlen_shift,ishift_obs,cc_max_obs) ! T(d1-d2)
       call xcorr_calc_mod(modary_tmp_t0,modary_tmp_t1,par%nt,i_start,i_end,par%nlen_shift,ishift_mod,cc_max_mod) ! T(d1-d2)

       tshift_mod = ishift_mod * deltat
       tshift_obs = ishift_obs * deltat

       !! double-difference cc-measurement
       ddtshift_cc = tshift_mod - tshift_obs

       dgg = dgg + (ddtshift_cc)**2

!       modary_tmp_t1_cc = cshift(modary_tmp_t1, shift = ishift_mod)
!       modary_tmp_t0_cc = cshift(modary_tmp_t0, shift = -ishift_mod)

!       do it = 2, par%nt-1
!          modary_tmp_t1_cc_vel(it) = (modary_tmp_t1_cc(it+1)-modary_tmp_t1_cc(it-1))/(2*deltat)
!          modary_tmp_t0_cc_vel(it) = (modary_tmp_t0_cc(it+1)-modary_tmp_t0_cc(it-1))/(2*deltat)
!       enddo

!       modary_tmp_t1_cc_vel(1) = (modary_tmp_t1_cc(2)-modary_tmp_t1_cc(1))/deltat
!       modary_tmp_t1_cc_vel(par%nt) = (modary_tmp_t1_cc(par%nt)-modary_tmp_t1_cc(par%nt-1))/deltat

!       modary_tmp_t0_cc_vel(1) = (modary_tmp_t0_cc(2)-modary_tmp_t0_cc(1))/deltat
!       modary_tmp_t0_cc_vel(par%nt) = (modary_tmp_t0_cc(par%nt)-modary_tmp_t0_cc(par%nt-1))/deltat

!       shot_virtual(:,ig)= + ddtshift_cc * modary_tmp_t1_cc_vel(:)

!       shot_virtual(:,ig)= - shot_virtual(:,ig)

     enddo

     deallocate(obsary_tmp_t1,obsary_tmp_t0)
     deallocate(modary_tmp_t1,modary_tmp_t0)
     deallocate(modary_tmp_t1_cc, modary_tmp_t0_cc)
     deallocate(modary_tmp_t1_cc_vel,modary_tmp_t0_cc_vel)

    endif !!(par%time_lapse==1)

    deallocate(obs_s_vel,mod_s_vel,obs_s_vel_shift)

  endif

  if (par%obtype.eq.1) then

  do ig=1,coord%ng(is)
  do it=1,par%nt
     dgg=dgg+shot_res_w(it,ig)*shot_res_w(it,ig)
  enddo
  enddo

  endif


!  dgg=sqrt(dgg)

  deallocate(seis_u,seis_w, seis_true)
  deallocate(shot_res_u, shot_res_w)
  deallocate(modary_tmp,obsary_tmp,resary_tmp)

  if(par%time_lapse) then
     deallocate(seis_true_t1, seis_true_t0, &
                      seis_w_syn_t0)
     deallocate(seis_w_syn_t0_tmp, seis_true_t1_tmp, &
                seis_true_t0_tmp)
  endif

!  if (.not.par%save_wavefield) then

  deallocate(wave_pux, wave_puz,&
             wave_pwx,wave_pwz)

!  endif

end subroutine analytical_step_length_elastic


end module fwi_elastic_subroutine
