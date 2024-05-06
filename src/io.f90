! module containing I/O utilities
! many of the codes are from legacy codes

module io

use string
implicit none

interface read_sufile
  module procedure read_sufile1
  module procedure read_sufile_2d
end interface

interface read_binfile
  module procedure read_binfile_1d
  module procedure read_binfile_2d
  module procedure read_binfile_3d
end interface

interface read_binfile_intel
  module procedure read_binfile_intel_1d
  module procedure read_binfile_intel_2d
  module procedure read_binfile_intel_3d
end interface

interface write_sufile
  module procedure write_sufile_1d
  module procedure write_sufile_2d
  module procedure write_sufile_2d_coord
end interface

interface write_binfile
  module procedure write_binfile_1d
  module procedure write_binfile_2d
  module procedure write_binfile_3d
end interface

interface write_binfile_intel
  module procedure write_binfile_intel_1d
  module procedure write_binfile_intel_2d
  module procedure write_binfile_intel_3d
end interface

contains

!---------------------------------------------------------------------------------------------
! Rules for verifying the input parameters
!
! STRING:
!         - required parameter will have a default as 'n/a'
!         - optional parameter will have some default which is not 'n/a'
!
! INTEGER:
!         - required parameter will have a default as -1
!         - optional parameter will have some default which is not -1
!
! FLOATING-POINT:
!         - required parameter will have a default as -1.0
!         - optional parameter will have some default which is not -1.0
!
subroutine readparamfile(parfile, par)

use datatype
use parser
use mmi_mpi

character(len=*), intent(in) :: parfile
type(param), intent(out)     :: par
character(len=100)           :: tmp, velfile, coordfile
integer                      :: nx, nz, nt
real                         :: dx, dt, f, dt_out

if (rank.eq.0) then

  !! Parameters required for most simulations
  !! ----------------------------------------
  call readParFile(parfile, 'VEL_IN',           par%velfile,            'n/a')
  call readParFile(parfile, 'COORD_FILE',       par%coordfile,          'n/a')
  call readParFile(parfile, 'CSG_IN',           par%csg_in,             'n/a')
  call readParFile(parfile, 'CSG_VIR',          par%csg_vir,            'n/a')
  call readParFile(parfile, 'TT_IN',            par%tt_in,              'n/a')
  call readParFile(parfile, 'ILLUM_FILE',       par%illumfile,          'n/a')
  call readParFile(parfile, 'Q_IN',             par%qualfile,           'n/a')
  call readParFile(parfile, 'SNAPSHOT_FILE',    par%snapshot_file,      'n/a')
  call readParFile(parfile, 'NX',               par%nx,                    -1)
  call readParFile(parfile, 'NZ',               par%nz,                    -1)
  call readParFile(parfile, 'NY',               par%ny,                    -1)
  call readParFile(parfile, 'NT_WORK',          par%nt,                    -1)
  call readParFile(parfile, 'DX',               par%dx,                  -1.0)
  call readParFile(parfile, 'DT_WORK',          par%dt,                  -1.0)
  call readParFile(parfile, 'FILEFORMAT',       par%fileformat,           'H')
  call readParFile(parfile, 'NPML',             par%npml,                  40)
  call readParFile(parfile, 'FREQUENCY',        par%f,                    5.0)
  call readParFile(parfile, 'FREQUENCY_MAX',    par%frequency_max,                25.0)
  call readParFile(parfile, 'STRATEGY',         par%strategy,      'snapshot')
  call readParFile(parfile, 'SKIPSHOT',         par%skipshot,               0)
  call readParFile(parfile, 'SKIPTRACE',        par%skiptrace,              0)
  call readParFile(parfile, 'SOURCEFILE',       par%sourcefile,         'n/a')
  call readParFile(parfile, 'SOURCETYPE',       par%sourcetype,      'normal')

  !! ----------------
  call readParFile(parfile, 'HIGHPASS',         par%highpass,               0)
  call readParFile(parfile, 'MIG_FILE',         par%migfile,      'mig_final')
  call readParFile(parfile, 'MIG_FILE_SHOT',    par%migfile_shot,  'mig_shot')
  call readParFile(parfile, 'REFL_FILE',        par%reflfile,           'n/a')
  call readParFile(parfile, 'REFL_DEN',         par%reflden,            'n/a')
  call readParFile(parfile, 'IMAGE_CONDITION',  par%ic,                     2)

  !! Parameters for FWI
  !! ----------------------------------------
  call readParFile(parfile, 'ITERMAX',          par%itermax,              100)
  call readParFile(parfile, 'MUTE_DIRECT',      par%mute_direct,            0)
  call readParFile(parfile, 'PRECONDITION',     par%pre,                    1)
  call readParFile(parfile, 'GRAD_FILE',        par%gradfile,           'n/a')
  call readParFile(parfile, 'LOGFILE',          par%logfile,       'logs/log')
  call readParFile(parfile, 'N_TAPER',          par%n_taper,               20)
  call readParFile(parfile, 'CSG_OUT',          par%csg_out,            'csg')
  call readParFile(parfile, 'ISMARINE',         par%ismarine,               0)
  call readParFile(parfile, 'WATERBOTTOMFILE',  par%waterbottomfile,    'n/a')
  call readParFile(parfile, 'RESIDUALFILE',     par%residualfile,   'csg_res')
  call readParFile(parfile, 'IZWB',             par%izwb,                  -1)
  call readParFile(parfile, 'WINDOW_SIZE',      par%window_size,          200)
  call readParFile(parfile, 'N_CUT',            par%n_cut,                 65)
  !! FWI parameters
  !! --------------
  call readParFile(parfile, 'NSEARCH_MAX',      par%nsearch_max,           10)
  call readParFile(parfile, 'LINE_SEARCH',      par%line_search,            1)
  call readParFile(parfile, 'NORMALIZE',        par%normalize,              0)
  call readParFile(parfile, 'SMOOTHGRAD',       par%smoothgrad,             0)
  call readParFile(parfile, 'SMOOTHVEL',        par%smoothvel,              0)
  call readParFile(parfile, 'VCONSTRAINT',      par%vconstraint,         -1.0)
  call readParFile(parfile, 'VEL_OUT',          par%velfile_out,        'vel')
  call readParFile(parfile, 'VMIN',             par%vmin,              1500.0)
  call readParFile(parfile, 'VMAX',             par%vmax,              4500.0)
  call readParFile(parfile, 'STEP',             par%step,                 5.0)
  call readParFile(parfile, 'SMOOTHING_WINDOW', par%smooth,                 5)
  call readParFile(parfile, 'CSG_OUT_RES',      par%csg_out_res,        'csg')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call readParFile(parfile, 'OBTYPE',           par%obtype,        1)
  call readParFile(parfile, 'NT_SHIFT', par%nt_shift,        0)

  call readParFile(parfile, 'NLCG',          par%nlcg,        0)
  call readParFile(parfile, 'COND_WB',          par%cond_wb,        0)

  call readParFile(parfile, 'WB_DEPTH_FILE',       par%wb_depth_file,        'wb')

  call readParFile(parfile, 'SMOOTH_TYPE',          par%smooth_type,        0)

  call readParFile(parfile, 'HORSMT',          par%horsmt,        10)
  call readParFile(parfile, 'VERSMT',          par%versmt,        4)

  call readParFile(parfile, 'AMP_NORM_TRACE',          par%amp_norm_trace,        0)
  call readParFile(parfile, 'AMP_NORM_CSG',          par%amp_norm_csg,        0)

  call readParFile(parfile, 'TAPER_SOURCE',          par%taper_source,        0)
  call readParFile(parfile, 'TAPER_GEOPHONE',          par%taper_geophone,        0)

                           !grad_illu
  call readParFile(parfile, 'GRAD_ILLU',          par%grad_illu,        1)

  call readParFile(parfile, 'PMIN',          par%pmin,        0.0)
  call readParFile(parfile, 'PMAX',          par%pmax,        0.01)

  call readParFile(parfile, 'FMIN',          par%fmin,        0.1)
  call readParFile(parfile, 'FMAX',          par%fmax,        2.0)
  call readParFile(parfile, 'DF_NEW',        par%df_new,      0.1)

  !! Time-lapse parameters
  call readParFile(parfile, 'TIME_LAPSE',          par%time_lapse,        0)
  call readParFile(parfile, 'CSG_T1_IN',           par%csg_t1_in,        'csg')
  call readParFile(parfile, 'CSG_SYN_T0_IN',       par%csg_syn_t0_in,        'csg')

  call readParFile(parfile, 'NLEN_SHIFT',  par%nlen_shift,50)

  call readParFile(parfile, 'LV',  par%lv,1)

  call readParFile(parfile, 'SAVE_WAVEFIELD',  par%save_wavefield,0)
  call readParFile(parfile, 'SAVE_WAVEFIELD_BOUNDARY',  par%save_wavefield_boundary,1)

  call readParFile(parfile, 'FD_ORDER', par%fd_order,        2) 
  call readParFile(parfile, 'FREESURFACE',      par%free_surface,           0)
  call readParFile(parfile, 'ADJ_LAGRAN',  par%adj_lagran, 1)
  call readParFile(parfile, 'SOURCE_MEC',  par%source_mec, 1)
  call readParFile(parfile, 'FS_TYPE',  par%fs_type, 1)
  call readParFile(parfile, 'FS_ZERO',  par%fs_zero, 0)
  call readParFile(parfile, 'FS_ADJ',  par%fs_adj, 0)

  call readParFile(parfile, 'GRAD_INV_PHASE',  par%grad_inv_phase, 0)
  call readParFile(parfile, 'GRAD_TWIST',  par%grad_twist, 1) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  end

  !! more parameters (not used so far)
  !! -------------
  call readParFile(parfile, 'TAU_MIN',          par%tau_min,           0.0002)
  call readParFile(parfile, 'TAU_MAX',          par%tau_max,             0.07)
  call readParFile(parfile, 'SMOOTH_TAU',       par%smooth_tau,             0)
  call readParFile(parfile, 'CSG_MOD',          par%csg_mod,        'csg_mod')
  call readParFile(parfile, 'METHOD',           par%method,                 1)
  call readParFile(parfile, 'METHOD',           par%simu_method,            1)
  !! Multisource parameters
  !! ----------------------
  call readParFile(parfile, 'ENCODING_CODES',   par%encoding_codes,     'n/a')
  call readParFile(parfile, 'ENCODING_POLARITY',par%encoding_polarity,  'n/a')

  !! Optional parameters (not used so far)
  !! -------------------------------------
  call readParFile(parfile, 'VEL_BG',           par%velfile_bg,   par%velfile)
  call readParFile(parfile, 'TOPO',             par%topofile,           'n/a')
  call readParFile(parfile, 'DATA',             par%data,          'pressure')
  call readParFile(parfile, 'DENSITYFILE',      par%densityfile,        'n/a')
  call readParFile(parfile, 'VPFILE',           par%vpfile,             'n/a')
  call readParFile(parfile, 'VSFILE',           par%vsfile,             'n/a')

  ! Required Integer parameters
  call readParFile(parfile, 'NT_IN',            par%nt_in,             par%nt)
  call readParFile(parfile, 'NT_OUT',           par%nt_out,            par%nt)
  call readParFile(parfile, 'SHIFT',            par%shift_wavefield,        0)
  call readParFile(parfile, 'VARIABLE_DENSITY', par%variable_density,       0)
  call readParFile(parfile, 'MUTE_DATA',        par%mute_data,              0)
  call readParFile(parfile, 'METHOD',           par%method,                 1)
  call readParFile(parfile, 'EXPORTWAVEFIELD',  par%exportwavefield,       -1)
  call readParFile(parfile, 'METHOD',           par%simu_method,            1)
  ! Floating-point parameters
  call readParFile(parfile, 'DT_IN',            par%dt_in,             par%dt)
  call readParFile(parfile, 'DT_OUT',           par%dt_out,            par%dt)
  call readParFile(parfile, 'XMIN',             par%xmin,                 0.0)
  call readParFile(parfile, 'OFFSET_MIN',       par%offset_min,           0.0)
  call readParFile(parfile, 'WINDOW',           par%window,              -1.0)
  call readParFile(parfile, 'OFFSET_MAX',       par%offset_max,real(par%nx-1)*par%dx)
  call readParFile(parfile, 'XMAX',             par%xmax,      real(par%nx-1)*par%dx)
 
  ! deblur
  call readParFile(parfile, 'IF_DEBLUR',        par%if_deblur,              0)
  call readParFile(parfile, 'NSCX',             par%nscx,                  -1)
  call readParFile(parfile, 'NSCZ',             par%nscz,                  -1)
  call readParFile(parfile, 'NFX',              par%nfx,                   -1)
  call readParFile(parfile, 'NFZ',              par%nfz,                   -1)
  call readParFile(parfile, 'NREFX',            par%nrefx,                 -1)
  call readParFile(parfile, 'NREFZ',            par%nrefz,                 -1)
  call readParFile(parfile, 'PETURB_TYPE',      par%peturb_type,            0) 

  ! LSM
  call readParFile(parfile, 'LSM_ITERMAX',      par%lsm_itermax,           10) 
  call readParFile(parfile, 'IF_LSM',           par%if_lsm,                 0)
 
endif

!! Broadcast parameters required for all simulations
!! -------------------------------------------------
call MPI_BCAST(par%velfile,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%coordfile,100,        MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%csg_in,100,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%csg_vir,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%tt_in,100,            MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%illumfile,100,        MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%qualfile,100,         MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%snapshot_file,100,    MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nx,1,                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nz,1,                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%ny,1,                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nt,1,                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%dx,1,                      MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%dt,1,                      MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fileformat,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%npml,1,                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%f,1,                       MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%frequency_max,1,                    MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%strategy,100,         MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%skipshot,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%skiptrace,1,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%free_surface,1,         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%sourcefile,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%sourcetype,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!! ---------------------------------------
call MPI_BCAST(par%highpass,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%migfile,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%migfile_shot,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%reflfile,100,         MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%reflden,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%ic,1,                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!! Broadcast parameters required for FWI
!! ----------------------------------------------------
call MPI_BCAST(par%itermax,1,              MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%mute_direct,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%pre,1,                  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%gradfile,100,         MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%logfile,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%n_taper,1,              MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%csg_out,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%ismarine,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%waterbottomfile,100,  MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%residualfile,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%izwb,1,                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%window_size,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%csg_out_res,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%obtype,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nt_shift,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%n_cut,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%nlcg,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%cond_wb,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%fd_order,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%wb_depth_file,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%smooth_type,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%horsmt,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%versmt,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%amp_norm_trace,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%amp_norm_csg,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%taper_source,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%taper_geophone,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%grad_illu,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!!
call MPI_BCAST(par%time_lapse,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
call MPI_BCAST(par%csg_t1_in,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%csg_syn_t0_in,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%nlen_shift,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%lv,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%save_wavefield,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%save_wavefield_boundary,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%adj_lagran,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%source_mec,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fs_type,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fs_zero,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fs_adj,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%pmin,1,                MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%pmax,1,                MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fmin,1,                MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fmax,1,                MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%df_new,1,              MPI_REAL,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%grad_inv_phase,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%grad_twist,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if (par%save_wavefield == 1) then 
        par%save_wavefield_boundary = 0
endif


!! Broadcast parameters required for FWI 
!! -------------------------------------
call MPI_BCAST(par%nsearch_max,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%line_search,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%normalize,1,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%smoothgrad,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%smoothvel,1,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vconstraint,1,             MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%velfile_out,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vmin,1,                    MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vmax,1,                    MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%step,1,                    MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%smooth,1,               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!! ------------------------------------
call MPI_BCAST(par%csg_mod,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%tau_min,1,                 MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%tau_max,1,                 MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%smooth_tau,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%method,1,               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%simu_method,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!! Broadcast parameters
!! ---------------------------------------------------
call MPI_BCAST(par%encoding_codes,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%encoding_polarity,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!! Broadcast optional parameters
!! -----------------------------
call MPI_BCAST(par%velfile_bg,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%topofile,100,         MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%data,100,             MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%densityfile,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vpfile,100,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vsfile,100,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nt_in,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nt_out,1,               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%shift_wavefield,1,      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%variable_density,1,     MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%mute_data,1,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%method,1,               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%exportwavefield,1,      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%dt_in,1,                   MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%dt_out,1,                  MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%xmin,1,                    MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%xmax,1,                    MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%offset_min,1,              MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%offset_max,1,              MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%window,1,                  MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%simu_method,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%if_deblur,1,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nscx,1,                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nscz,1,                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nfx,1,                  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nfz,1,                  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nrefx,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nrefz,1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%peturb_type,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%lsm_itermax,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%if_lsm,1,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end subroutine readparamfile

!-------------------------------------------------------------------------------
subroutine readcoordfile(coordfile, coord)

use datatype
use mmi_mpi

character(len=*), intent(in)   :: coordfile
type(acquisition), intent(out) :: coord
integer                        :: i, j, is0, is, ig, ns, ng, ngmax
real                           :: xs, zs, xg, zg, t

if (rank == 0) then
  open(10,file=coordfile,form='formatted')

  ! Determine the number of sources
  ns = 0
  read(10,*,end=100) is0, ig, xs, zs, xg, zg, t
  ns = 1
  do i=2,1000000000
    read(10,*,end=100) is, ig, xs, zs, xg, zg, t
    if (is.ne.is0) then
      is0 = is
      ns = ns + 1
    endif
  enddo
  100 continue
  rewind(10)
  coord%ns = ns
  allocate(coord%ng(ns))
  allocate(coord%xs(ns))
  allocate(coord%zs(ns))

  ! Determine the number of geophones for each shot
  read(10,*,end=200) is0, ig, xs, zs, xg, zg, t
  ns = 1
  ng = 1
  do i=2,1000000000
    read(10,*,end=200) is, ig, xs, zs, xg, zg, t
    if (is.ne.is0) then
      is0 = is
      coord%ng(ns) = ng
!    write(*,*) 'ng(',ns,') = ', coord%ng(ns)
      ns = ns + 1
      ng = 1
    else
      ng = ng + 1
    endif
  enddo
  200 continue
  coord%ng(ns) = ng
  !write(*,*) 'ng(',ns,') = ', coord%ng(ns)
  rewind(10)

  ! Determine the maximum number of geophones per shot
  coord%ngmax = 0
  do is=1,coord%ns
    if (coord%ngmax < coord%ng(is)) coord%ngmax = coord%ng(is)
  enddo
  write(*,*) 'ng max = ', coord%ngmax

  allocate(coord%xg(is,coord%ngmax))
  allocate(coord%zg(is,coord%ngmax))
  allocate(coord%t(is,coord%ngmax))

  ! Read source and receiver positions
  do i=1,coord%ns
    do j=1,coord%ng(i)
      read(10,*,end=300) is, ig, xs, zs, xg, zg, t
      coord%xs(i) = xs
      coord%zs(i) = zs
      coord%xg(i,j) = xg
      coord%zg(i,j) = zg
      coord%t(i,j) = t
    enddo
  enddo
  300 continue
  close(10)

!  open(10,file='coord.txt',form='formatted')
!  do is=1,coord%ns
!    do ig=1,coord%ng(is)
!      write(10,*) is,ig,coord%xs(is),coord%zs(is),coord%xg(is,ig),coord%zg(is,ig),coord%t(is,ig)
!    enddo
!  enddo
!  close(10)
  write(*,*) 'ns = ', ns
endif

call MPI_BCAST(coord%ns,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (rank > 0) then
  allocate(coord%ng(coord%ns))
  allocate(coord%xs(coord%ns))
  allocate(coord%zs(coord%ns))
endif
call MPI_BCAST(coord%ng,coord%ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(coord%xs,coord%ns,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(coord%zs,coord%ns,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(coord%ngmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (rank > 0) then
  allocate(coord%xg(coord%ns,coord%ngmax))
  allocate(coord%zg(coord%ns,coord%ngmax))
  allocate(coord%t(coord%ns,coord%ngmax))
endif
do is=1,coord%ns
  call MPI_BCAST(coord%xg(is,:),coord%ng(is),MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(coord%zg(is,:),coord%ng(is),MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(coord%t(is,:),coord%ng(is),MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

999 continue

end subroutine readcoordfile

!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------------
subroutine readvelfile(par,c,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: c(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml
integer :: ix, iz

if (rank == 0) then
  ! Read velocity model
  open(10,file=par%velfile,access='direct',recl=i4*par%nz)
  do ix=1,par%nx
    read(10,rec=ix) c(npml+1:npml+par%nz,ix+npml)
  enddo
  close(10)

  ! Extrapolate velocity in PML regions
  do ix=1,npml
    c(npml+1:npml+par%nz,ix) = c(npml+1:npml+par%nz,npml+1)
    c(npml+1:npml+par%nz,nx_pml-npml+ix) = c(npml+1:npml+par%nz,nx_pml-npml)
  enddo
  
  do iz=1,npml
    c(iz,:) = c(npml+1,:)
    c(nz_pml-npml+iz,:) = c(nz_pml-npml,:)
  enddo

endif

do ix=1,nx_pml
  call MPI_BCAST(c(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

end subroutine readvelfile

!-------------------------------------------------------------------

subroutine readvpfile(par,c,npml,nx_pml,nz_pml,fd_order,free_surface)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: c(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml, free_surface, fd_order
integer :: ix, iz, pad_top

!pad_top=(fd_order-20)/2+1

if (rank == 0) then
  ! Read velocity model
!  if(free_surface==0)then
     
     open(10,file=par%vpfile,access='direct',recl=i4*par%nz)
     do ix=1,par%nx
        read(10,rec=ix) c(npml+1:npml+par%nz,ix+npml)
     enddo
     close(10)
     ! Extrapolate velocity in PML regions
     do ix=1,npml
        c(npml+1:npml+par%nz,ix) = c(npml+1:npml+par%nz,npml+1)
        c(npml+1:npml+par%nz,nx_pml-npml+ix) = c(npml+1:npml+par%nz,nx_pml-npml)
     enddo
     do iz=1,npml
        c(iz,:) = c(npml+1,:)
        c(nz_pml-npml+iz,:) = c(nz_pml-npml,:)
     enddo

endif
do ix=1,nx_pml
  call MPI_BCAST(c(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

end subroutine readvpfile

!---------------------------------------------------------------

subroutine readvsfile(par,c,npml,nx_pml,nz_pml,fd_order,free_surface)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: c(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml, free_surface, fd_order
integer :: ix, iz, pad_top

!pad_top=(fd_order-20)/2+1

if (rank == 0) then

  ! Read velocity model
!  if(free_surface==0)then

     open(10,file=par%vsfile,access='direct',recl=i4*par%nz)
     do ix=1,par%nx
        read(10,rec=ix) c(npml+1:npml+par%nz,ix+npml)
     enddo
     close(10)
     ! Extrapolate velocity in PML regions
     do ix=1,npml
        c(npml+1:npml+par%nz,ix) = c(npml+1:npml+par%nz,npml+1)
        c(npml+1:npml+par%nz,nx_pml-npml+ix) = c(npml+1:npml+par%nz,nx_pml-npml)
     enddo
     do iz=1,npml
        c(iz,:) = c(npml+1,:)
        c(nz_pml-npml+iz,:) = c(nz_pml-npml,:)
     enddo

endif

do ix=1,nx_pml
  call MPI_BCAST(c(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

end subroutine readvsfile
!------------------------------------------------------------------------------------

subroutine readreflfile(par,mig,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: mig(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml
integer :: ix, iz

if (rank == 0) then
  ! Read velocity model
  open(10,file=par%reflfile,access='direct',recl=i4*par%nz,convert='little_endian')
  do ix=1,par%nx
    read(10,rec=ix) mig(npml+1:npml+par%nz,ix+npml)
  enddo
  close(10)

  ! Extrapolate velocity in PML regions
  do ix=1,npml
    mig(npml+1:npml+par%nz,ix) = mig(npml+1:npml+par%nz,npml+1)
    mig(npml+1:npml+par%nz,nx_pml-npml+ix) = mig(npml+1:npml+par%nz,nx_pml-npml)
  enddo
  do iz=1,npml
    mig(iz,:) = mig(npml+1,:)
    mig(nz_pml-npml+iz,:) = mig(nz_pml-npml,:)
  enddo
endif
do ix=1,nx_pml
  call MPI_BCAST(mig(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

end subroutine readreflfile


subroutine readrefldenfile(par,mig,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: mig(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml
integer :: ix, iz

if (rank == 0) then
  ! Read velocity model
  open(10,file=par%reflden,access='direct',recl=i4*par%nz,convert='little_endian')
  do ix=1,par%nx
    read(10,rec=ix) mig(npml+1:npml+par%nz,ix+npml)
  enddo
  close(10)

  ! Extrapolate velocity in PML regions
  do ix=1,npml
    mig(npml+1:npml+par%nz,ix) = mig(npml+1:npml+par%nz,npml+1)
    mig(npml+1:npml+par%nz,nx_pml-npml+ix) = mig(npml+1:npml+par%nz,nx_pml-npml)
  enddo
  do iz=1,npml
    mig(iz,:) = mig(npml+1,:)
    mig(nz_pml-npml+iz,:) = mig(nz_pml-npml,:)
  enddo
endif
do ix=1,nx_pml
  call MPI_BCAST(mig(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo


end subroutine readrefldenfile



!-------------------------------------------------
subroutine readqualfile(par,qf,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: qf(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml
integer :: ix, iz

if (rank == 0) then
  ! Read velocity model
  open(10,file=par%qualfile,access='direct',recl=i4*par%nz,convert='little_endian')
  do ix=1,par%nx
    read(10,rec=ix) qf(npml+1:npml+par%nz,ix+npml)
  enddo
  close(10)

  ! Extrapolate Q in PML regions
  do ix=1,npml
    qf(npml+1:npml+par%nz,ix) = qf(npml+1:npml+par%nz,npml+1)
    qf(npml+1:npml+par%nz,nx_pml-npml+ix) = qf(npml+1:npml+par%nz,nx_pml-npml)
  enddo
  do iz=1,npml
    qf(iz,:) = qf(npml+1,:)
    qf(nz_pml-npml+iz,:) = qf(nz_pml-npml,:)
  enddo
endif
do ix=1,nx_pml
  call MPI_BCAST(qf(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

end subroutine readqualfile

!------------------------------------------------------------------------------------
subroutine readvelfile_new(velfile,v,nx,nz,npml)

use datatype
use global
use mmi_mpi

character(*), intent(in)  :: velfile
real,         intent(out) :: v(:,:)
integer,      intent(in)  :: nx, nz, npml
integer                   :: ix, iz, nx_pml, nz_pml

if (rank == 0) then
  nx_pml = nx+2*npml
  nz_pml = nz+2*npml

  ! Read velocity model
  open(10,file=velfile,access='direct',recl=i4*nz)
  do ix=1,nx
    read(10,rec=ix) v(npml+1:npml+nz,ix+npml)
  enddo
  close(10)

  ! Extrapolate velocity in PML regions
  do ix=1,npml
    v(npml+1:npml+nz,ix) = v(npml+1:npml+nz,npml+1)
    v(npml+1:npml+nz,nx_pml-npml+ix) = v(npml+1:npml+nz,nx_pml-npml)
  enddo
  do iz=1,npml
    v(iz,:) = v(npml+1,:)
    v(nz_pml-npml+iz,:) = v(nz_pml-npml,:)
  enddo
endif
!do ix=1,nx_pml
!  call MPI_BCAST(vp(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!enddo

end subroutine readvelfile_new



!------------------------------------------------------------------------------------
subroutine read_densityfile(par,c,den,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi
use modeling

type(param),intent(in) :: par
real, intent(in)       :: c(:,:)
real, intent(out)      :: den(:,:)
integer, intent(in)    :: npml, nx_pml, nz_pml
integer                :: ix, iz

!if (par%variable_density == 1) then
  if (par%densityfile(1:3) == 'n/a') then
    call compute_density(par%variable_density, c, den, nx_pml, nz_pml)
  else
    if (rank == 0) then
      ! Read density model
      call read_binfile_2d(par%densityfile,den(npml+1:npml+par%nz,npml+1:npml+par%nx),par%nz,par%nx)

      ! Extrapolate density in PML regions
      do ix=1,npml
        den(npml+1:npml+par%nz,ix) = den(npml+1:npml+par%nz,npml+1)
        den(npml+1:npml+par%nz,nx_pml-npml+ix) = den(npml+1:npml+par%nz,nx_pml-npml)
      enddo
      do iz=1,npml
        den(iz,:) = den(npml+1,:)
        den(nz_pml-npml+iz,:) = den(nz_pml-npml,:)
      enddo
    endif
    do ix=1,nx_pml
      call MPI_BCAST(den(:,ix),nz_pml,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    enddo
  endif
!else
!  den = 1.0
!endif

end subroutine read_densityfile



!------------------------------------------------------------------------------
subroutine read_sufile1(sufile, trace, nt, ntrace, sampling_interval)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: nt, ntrace
real, intent(out) :: sampling_interval
real, intent(out) :: trace(:,:)
integer :: nlen, i

nlen = nhead+nt
open(10,file=sufile,access='direct',recl=i4*nlen)
do i=1,ntrace
  read(10,rec=i,err=111) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass,trace(i,1:nt)
enddo
111 continue
close(10)

sampling_interval = real(dt)/1.0e6

end subroutine read_sufile1

!------------------------------------------------------------------------------
subroutine read_sufile_2d(sufile, trace, n1, n2, dd1, dd2)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1, n2
real, intent(out) :: dd1, dd2
real, intent(out) :: trace(:,:)
integer :: nlen, i

nlen = nhead+n1
open(10,file=sufile,access='direct',recl=i4*nlen)
do i=1,n2
  read(10,rec=i,err=111) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass,trace(1:n1,i)
enddo
111 continue
close(10)

dd1 = real(dt)/1.0e6
dd2 = d2
if (abs(dd2) < 1.0e-6) dd2 = 1

end subroutine read_sufile_2d

!-----------------------------------------------------------------------------------
subroutine write_sufile_1d(sufile, trace, n1, dd1)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1
real, intent(in) :: dd1, trace(:)
integer :: nlen

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = 1
dt = int(dd1*1.0e6)
d1 = dd1
f1 = 0
tracl = 0
tracr = 0

open(10,file=sufile,access='direct',recl=nlen*i4,status='replace')
tracl = tracl + 1
tracr = tracr + 1
cdp = 1
cdpt = 1
write(10,rec=1) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1)
close(10)

end subroutine write_sufile_1d

!-----------------------------------------------------------------------------------
subroutine write_sufile_2d(sufile, trace, n1, n2, dd1, dd2)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1, n2
real, intent(in) :: dd1, dd2, trace(:,:)
integer :: nlen, i

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = 1
dt = int(dd1*1.0e6)
d1 = dd1
d2 = dd2
f1 = 0
f2 = 1.0
tracl = 0
tracr = 0

open(10,file=sufile,access='direct',recl=nlen*i4,status='replace')
do i = 1,n2
  tracl = tracl + 1
  tracr = tracr + 1
  cdp = i
  cdpt = i
  write(10,rec=i) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1,i)
enddo
close(10)

end subroutine write_sufile_2d

!-----------------------------------------------------------------------------------
subroutine write_sufile_2d_offset(sufile, starting_record, is, trace, n1, n2, dd1, dd2, coord)

use global
use datatype
use su

type(acquisition), intent(in) :: coord
character(len=*), intent(in) :: sufile
integer, intent(in) :: starting_record, is, n1, n2
real, intent(in) :: dd1, dd2, trace(:,:)
integer :: nlen, i

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = is
dt = int(dd1*1.0e6)
d1 = dd1
d2 = dd2
f1 = 0
f2 = 1.0
tracl = 0
tracr = 0
scalel = -10
scalco = -10

open(10,file=sufile,access='direct',recl=nlen*i4)
do i = 1,n2
  tracl = tracl + 1
  tracr = tracr + 1
  cdp = i
  cdpt = i
  sx    = nint(coord%xs(is)*10.0)
  sy    = 0
  selev = nint(coord%zs(is)*10.0)
  gx    = nint(coord%xg(is,i)*10.0)
  gy    = 0
  gelev = nint(coord%zg(is,i)*10.0)
  
  write(10,rec=starting_record+i) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1,i)
enddo
close(10)

end subroutine write_sufile_2d_offset

!-----------------------------------------------------------------------------------
subroutine write_sufile_2d_coord(sufile, trace, n1, n2, dd1, dd2, xs, xg)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1, n2
real, intent(in) :: dd1, dd2, trace(:,:), xs, xg(n2)
integer :: nlen, i

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = 1
dt = int(dd1*1.0e6)
d1 = dd1
d2 = dd2
f1 = 0
f2 = 1.0
tracl = 0
tracr = 0
sx = nint(xs)

open(10,file=sufile,access='direct',recl=nlen*i4,status='replace')
do i = 1,n2
  tracl = tracl + 1
  tracr = tracr + 1
!  cdp = i
!  cdpt = i
  gx = nint(xg(i))
  cdp = nint(0.5*(xs+xg(i)))
  cdpt = cdp
  offset = gx - sx
  write(10,rec=i) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1,i)
enddo
close(10)

end subroutine write_sufile_2d_coord
!-----------------------------------------------------------------------------------

subroutine write_binfile_1d(filename, data, n)

use global

character(len=*), intent(in) :: filename
integer, intent(in) :: n
real, intent(in) :: data(:)
integer :: i,it

open(10,file=filename,access='direct',recl=n*i4,status='replace')
write(10,rec=1) data(1:n)
close(10)

end subroutine write_binfile_1d
!-----------------------------------------------------------------------------------

subroutine write_binfile_2d(filename, data, n1, n2)

use global

character(len=*), intent(in) :: filename
integer, intent(in) :: n1, n2
real, intent(in) :: data(:,:)
integer :: i,it

open(10,file=filename,access='direct',recl=n1*i4,status='replace')

do i = 1, n2
   write(10,rec=i)(data(it,i), it = 1, n1)
enddo
close(10)

end subroutine write_binfile_2d

!-----------------------------------------------------------------------------------
subroutine write_binfile_3d(filename, data, n1, n2, n3)

character(len=*), intent(in) :: filename
integer,          intent(in) :: n1, n2, n3
real,             intent(in) :: data(:,:,:)
integer,          parameter  :: i4=4
integer                      :: i1, i2, i3

open(10,file=filename,access='direct',recl=n1*i4,status='replace')

do i3=1,n3
  do i2=1,n2
    write(10,rec=i2+(i3-1)*n2) data(:,i2,i3)
  enddo
enddo
close(10)

end subroutine write_binfile_3d

!-----------------------------------------------------------------------------------
subroutine write_binfile_intel_1d(filename, data, n)

use global

character(len=*), intent(in) :: filename
integer, intent(in) :: n
real, intent(in) :: data(:)
integer :: i,it

open(10,file=filename,access='direct',recl=n*i4,status='replace')
write(10,rec=1) data(1:n)
close(10)

end subroutine write_binfile_intel_1d

!-----------------------------------------------------------------------------------

subroutine write_binfile_intel_2d(filename, data, n1, n2)

use global

character(len=*), intent(in) :: filename
integer, intent(in) :: n1, n2
real, intent(in) :: data(:,:)
integer :: i,it

open(10,file=filename,access='direct',recl=n1*n2*i4,status='replace')

write(10,rec=1) data
close(10)

end subroutine write_binfile_intel_2d

!-----------------------------------------------------------------------------------
subroutine write_binfile_intel_3d(filename, data, n1, n2, n3)

character(len=*), intent(in) :: filename
integer,          intent(in) :: n1, n2, n3
real,             intent(in) :: data(:,:,:)
integer,          parameter  :: i4=4
integer                      :: i1, i2, i3

open(10,file=filename,access='direct',recl=n1*n2*n3*i4,status='replace')

write(10,rec=1) data
close(10)

end subroutine write_binfile_intel_3d

!-----------------------------------------------------------------------------------
subroutine read_binfile_1d(filename, data, n)

use global

implicit none

character(len=*), intent(in) :: filename
integer, intent(in)          :: n
real                         :: data(:)

open(10,file=filename,access='direct',recl=n*i4)
read(10, rec=1) data(1:n)
close(10)

end subroutine read_binfile_1d

!-----------------------------------------------------------------------------------
subroutine read_binfile_2d(filename, data, n1, n2)

use global

implicit none

character(len=*), intent(in) :: filename
integer, intent(in) :: n1, n2
real                :: data(:,:)
integer :: i1,i

open(10,file=filename,access='direct',recl=n1*i4)

do i = 1, n2
   read(10, rec=i)(data(i1,i),i1=1,n1)
enddo
close(10)

end subroutine read_binfile_2d

!-----------------------------------------------------------------------------------
subroutine read_binfile_3d(filename, data, n1, n2, n3)

character(len=*), intent(in)  :: filename
integer,          intent(in)  :: n1, n2, n3
real,             intent(out) :: data(:,:,:)
integer,          parameter   :: i4=4
integer                       :: i1, i2, i3

open(10,file=filename,access='direct',recl=n1*i4,status='old')

do i3=1,n3
  do i2=1,n2
    read(10,rec=i2+(i3-1)*n2) data(:,i2,i3)
  enddo
enddo
close(10)

end subroutine read_binfile_3d

!-----------------------------------------------------------------------------------
subroutine read_binfile_intel_1d(filename, data, n)

use global

implicit none

character(len=*), intent(in) :: filename
integer, intent(in)          :: n
real                         :: data(:)

open(10,file=filename,access='direct',recl=n*i4)
read(10, rec=1) data(1:n)
close(10)

end subroutine read_binfile_intel_1d

!-----------------------------------------------------------------------------------
subroutine read_binfile_intel_2d(filename, data, n1, n2)

use global

implicit none

character(len=*), intent(in) :: filename
integer, intent(in) :: n1, n2
real                :: data(:,:)
integer :: i1,i

open(10,file=filename,access='direct',recl=n1*n2*i4)
read(10, rec=1) data
close(10)

end subroutine read_binfile_intel_2d

!-----------------------------------------------------------------------------------
subroutine read_binfile_intel_3d(filename, data, n1, n2, n3)

character(len=*), intent(in)  :: filename
integer,          intent(in)  :: n1, n2, n3
real,             intent(out) :: data(:,:,:)
integer,          parameter   :: i4=4
integer                       :: i1, i2, i3

open(10,file=filename,access='direct',recl=n1*n2*n3*i4,status='old')
read(10,rec=1) data
close(10)

end subroutine read_binfile_intel_3d



!-----------------------------------------------------------------------------------
subroutine clean_suheader

use su
implicit none

integer i

tracl = 0
tracr = 0
fldr = 0
tracf = 0
ep = 0
cdp = 0
cdpt = 0
trid = 0
nvs = 0
nhs = 0
duse = 0
offset = 0
gelev = 0
selev = 0
sdepth = 0
gdel = 0
sdel = 0
swdep = 0
gwdep = 0
scalel = 0
scalco = 0
sx = 0
sy = 0
gx = 0
gy = 0
counit = 0
wevel = 0
swevel = 0
sut = 0
gut = 0
sstat = 0
gstat = 0
tstat = 0
laga = 0
lagb = 0
delrt = 0
muts = 0
mute = 0
ns = 0
dt = 0
gain = 0
igc = 0
igi = 0
corr = 0
sfs = 0
sfe = 0
slen = 0
styp = 0
stas = 0
stae = 0
tatyp = 0
afilf = 0
afils = 0
nofilf = 0
nofils = 0
lcf = 0
hcf = 0
lcs = 0
hcs = 0
year = 0
day = 0
hour = 0
minute = 0
sec = 0
timbas = 0
trwf = 0
grnors = 0
grnofr = 0
grnlof = 0
gaps = 0
otrav = 0
d1 = 0.0
f1 = 0.0
d2 = 0.0
f2 = 0.0
ungpow = 0.0
unscale = 0.0
mark = 0
unass(:) = 0

end subroutine clean_suheader


!-----------------------------------------------------------------------------------
subroutine writehead(headername,binaryname,esize,data_format,&
                     n1, n2, n3, n4, n5,                     &
                     d1, d2, d3, d4, d5,                     &
                     o1, o2, o3, o4, o5,                     &
                     l1, l2, l3, l4, l5)

    character(len=*),       intent(inout) :: headername
    character(len=*),optional,intent(in)  :: binaryname
    integer,         optional,intent(in)  :: esize
    integer,         optional,intent(in)  :: n1, n2, n3, n4, n5
    real,            optional,intent(in)  :: d1, d2, d3, d4, d5
    real,            optional,intent(in)  :: o1, o2, o3, o4, o5
    character(len=*),optional,intent(in)  :: l1, l2, l3, l4, l5
    character(len=*),optional,intent(in)  :: data_format

    character(len=512) :: tmp_str
    integer            :: iunit,maxdim
    character(len=32)  :: local_data_format
    integer            :: ierr

    ierr=0

    iunit=96

    maxdim=5

    open(unit=iunit,file=trim(headername),form='formatted',position='append',iostat=ierr)

    write(iunit,'(a)') '#SEP'

    if (present(n1)) then
       write(tmp_str,*) n1
       write(iunit,'(a)') 'n1='//trim(adjustl(tmp_str))
    endif

    if (present(n2).and.(maxdim.ge.2)) then
       write(tmp_str,*) n2
       write(iunit,'(a)') 'n2='//trim(adjustl(tmp_str))
    endif

    if (present(n3).and.(maxdim.ge.3)) then
       write(tmp_str,*) n3
       write(iunit,'(a)') 'n3='//trim(adjustl(tmp_str))
    endif

    if (present(n4).and.(maxdim.ge.4)) then
       write(tmp_str,*) n4
       write(iunit,'(a)') 'n4='//trim(adjustl(tmp_str))
    endif

    if (present(n5).and.(maxdim.ge.5)) then
       write(tmp_str,*) n5
       write(iunit,'(a)') 'n5='//trim(adjustl(tmp_str))
    endif

    if (present(d1)) then
       write(tmp_str,*) d1
       write(iunit,'(a)') 'd1='//trim(adjustl(tmp_str))
    endif

    if (present(d2).and.(maxdim.ge.2)) then
       write(tmp_str,*) d2
       write(iunit,'(a)') 'd2='//trim(adjustl(tmp_str))
    endif

    if (present(d3).and.(maxdim.ge.3)) then
       write(tmp_str,*) d3
       write(iunit,'(a)') 'd3='//trim(adjustl(tmp_str))
    endif

    if (present(d4).and.(maxdim.ge.4)) then
       write(tmp_str,*) d4
       write(iunit,'(a)') 'd4='//trim(adjustl(tmp_str))
    endif

    if (present(d5).and.(maxdim.ge.5)) then
       write(tmp_str,*) d5
       write(iunit,'(a)') 'd5='//trim(adjustl(tmp_str))
    endif

    if (present(o1)) then
       write(tmp_str,*) o1
       write(iunit,'(a)') 'o1='//trim(adjustl(tmp_str))
    endif

    if (present(o2).and.(maxdim.ge.2)) then
       write(tmp_str,*) o2
       write(iunit,'(a)') 'o2='//trim(adjustl(tmp_str))
    endif

    if (present(o3).and.(maxdim.ge.3)) then
       write(tmp_str,*) o3
       write(iunit,'(a)') 'o3='//trim(adjustl(tmp_str))
    endif

    if (present(o4).and.(maxdim.ge.4)) then
       write(tmp_str,*) o4
       write(iunit,'(a)') 'o4='//trim(adjustl(tmp_str))
    endif

    if (present(o5).and.(maxdim.ge.5)) then
       write(tmp_str,*) o5
       write(iunit,'(a)') 'o5='//trim(adjustl(tmp_str))
    endif

    if (present(l1)) then
       write(iunit,'(a)') 'label1="'//trim(l1)//'"'
    endif

    if (present(l2).and.(maxdim.ge.2)) then
       write(iunit,'(a)') 'label2="'//trim(l2)//'"'
    endif

    if (present(l3).and.(maxdim.ge.3)) then
       write(iunit,'(a)') 'label3="'//trim(l3)//'"'
    endif

    if (present(l4).and.(maxdim.ge.4)) then
       write(iunit,'(a)') 'label4="'//trim(l4)//'"'
    endif

    if (present(l5).and.(maxdim.ge.5)) then
       write(iunit,'(a)') 'label5="'//trim(l5)//'"'
    endif

    write(iunit,'(a)') 'esize='//trim(adjustl("4"))

    write(iunit,'(a)') 'in='//trim(adjustl(headername))//"@"

    local_data_format="native_float"
    write(iunit,'(a)') 'data_format='//trim(adjustl(local_data_format))

    flush(iunit)

    close(iunit,iostat=ierr)

    return

  end subroutine writehead


end module io

