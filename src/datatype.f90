! module for defining data structures
!
module datatype

 type param

  !! Parameters required for most simulations
  !! -----------------------------------------
  integer            :: nx,ny,nz,nt,free_surface,fs_type,fs_adj,fs_zero,skipshot,skiptrace,npml,first_shot,last_shot
  character(len=100) :: velfile,coordfile,fileformat,sourcefile,sourcetype,csg_in,csg_vir,tt_in,  &
                        illumfile,qualfile,snapshot_file,strategy
  real               :: dx,dt,f,frequency_max,cmin,cmax,denmin,denmax,qmin,qmax,reflmin,reflmax

  real               :: pmin, pmax, fmin, fmax, df_new

  !! Parameters required for FWI 
  !! -----------------------------
  integer             :: highpass,ic
  double precision    :: alpha,beta,gama
  character(len=100)  :: migfile,migfile_shot,reflfile,reflden

  !! ------------------------------------------
  integer             :: itermax,lsm_itermax,mute_direct,pre,n_taper,ismarine,izwb,window_size,iter,lsm_iter,n_cut
  character(len=100)  :: gradfile,logfile,csg_out,waterbottomfile,residualfile

  !! ---------------------------
  integer             :: nsearch_max,line_search,normalize,smoothgrad,smoothvel,smooth
  real                :: vmin,vmax,step,alpha_fwi,res_new,vconstraint,beta_fwi
  character(len=100)  :: velfile_out,csg_out_res

  integer             :: nlcg

  integer             :: obtype

  integer             :: fd_order

  integer             :: nt_shift

  integer             :: taper_source, taper_geophone

  integer             :: grad_illu

  integer             :: cond_wb 

  integer             :: smooth_type

  integer             :: horsmt, versmt

  integer             :: amp_norm_trace, amp_norm_csg

  !! Parameters for Time Lapse
  integer             :: time_lapse

  integer             :: nlen_shift

  integer             :: lv

  integer             :: adj_lagran
  integer             :: source_mec

  integer             :: save_wavefield, save_wavefield_boundary

  integer             :: grad_inv_phase

  integer             :: grad_twist

  character(len=100)  :: csg_t1_in, csg_syn_t0_in

  character(len=100)  :: wb_depth_file

  !! Parameters required for WQ
  !! --------------------------
  integer             :: smooth_tau
  real                :: tau_min,tau_max,alpha_wq,tau_constraint,beta_wq
  character(len=100)  :: csg_mod

  integer             :: min_num_trace

  !! Parameters for multisource modeling and migration
  !! -------------------------------------------------
  integer            :: max_delay,sg,ns_sg 
  character(len=100) :: encoding_codes,encoding_polarity

  !! Optional parameters
  !! -------------------
  integer            :: mute_data,method,exportwavefield,variable_density,nt_in,shift_wavefield,nt_out, &
                        lbfgs_start,lbfgs_end,simu_method,if_deblur,if_lsm,nscx,nscz,nfx,nfz,nrefx,nrefz,peturb_type
  real               :: dt_out,dt_in,window,xmin,xmax,offset_min,offset_max 
  character(len=100) :: densityfile,topofile,velfile_bg,vpfile,vsfile,data


 end type param

 !!!!!!!!!!!!!!!!!!!!!!
 type wave_boundary_param

   real,allocatable :: u_bl(:,:,:), u_br(:,:,:), u_bt(:,:,:), u_bb(:,:,:)
   real,allocatable :: w_bl(:,:,:), w_br(:,:,:), w_bt(:,:,:), w_bb(:,:,:)

   real,allocatable :: xx_bl(:,:,:), xx_br(:,:,:), xx_bt(:,:,:), xx_bb(:,:,:)
   real,allocatable :: zz_bl(:,:,:), zz_br(:,:,:), zz_bt(:,:,:), zz_bb(:,:,:)
   real,allocatable :: xz_bl(:,:,:), xz_br(:,:,:), xz_bt(:,:,:), xz_bb(:,:,:)

   real,allocatable :: last_u(:,:), last_w(:,:)

   real,allocatable :: last_xx(:,:), last_zz(:,:), last_xz(:,:)

 end type wave_boundary_param
 !!!!!!!!!!!!!!!!!!!!!!

 type acquisition
    integer              :: ns,ngmax
    integer, allocatable :: ng(:)
    real,    allocatable :: xs(:),zs(:),xg(:,:),zg(:,:),t(:,:)
 end type acquisition

 type acquisition3d
    integer              :: ns,ngmax
    integer, allocatable :: ng(:)
    real,    allocatable :: xs(:),ys(:),zs(:),xg(:,:),yg(:,:),zg(:,:),t(:,:)
 end type acquisition3d

 type wavelet_param
    character(len=100)   :: fileformat,inputprefix,outputprefix
    integer              :: ns,ng,nt,iref
    real                 :: dx,dt,ds,dg
 end type wavelet_param

 type boundary3d
    real, allocatable    :: top(:,:,:,:),bottom(:,:,:,:),left(:,:,:,:),right(:,:,:,:), &
                            front(:,:,:,:),back(:,:,:,:)
 end type boundary3d

 type boundary2d
    real, allocatable    :: top(:,:,:),bottom(:,:,:),left(:,:,:),right(:,:,:)
 end type boundary2d

end module datatype

