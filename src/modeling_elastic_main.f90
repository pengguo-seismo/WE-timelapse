
module modeling_elastic_main

 use datatype
 use io
 use math
 use mmi_mpi
 use parser
 use modeling
 use pml
 use source
 use global
 use smooth_slow
 use staggered_elastic_kernel

 implicit none

 type(param),       private              :: par
 type(acquisition), private              :: coord

 integer,           private              :: ix,iz,it,isx,isz,igx,igz,is,ig,i1,i2,&
                                            is1,is2
 real,              private              :: d1,d2,st_time,ed_time
 integer,           private, allocatable :: fs(:), izwb(:)

 double precision                        :: res,gg

 !! Medium properties
 !! -----------------
 real,              private, allocatable :: kappa(:,:),vp(:,:),vs(:,:),den(:,:),qf(:,:),tau(:,:), &
                                            tausigma(:,:),tauepsilon(:,:)

 !! --------------------------
 real,              private, allocatable :: s(:),seis_true_u(:,:),seis_true_w(:,:),seis(:,:), &
                                            energy_temp(:,:),energy(:,:),mig_temp(:,:),mig(:,:), &
                                            mig_sum(:,:),energy_sum(:,:),taper(:),wave_pux(:,:,:),wave_puz(:,:,:),wave_pwx(:,:,:),wave_pwz(:,:,:)

 contains

!-------------------------------------
subroutine modeling_elastic(parfile)

 
 character(len=*), intent(in) :: parfile
 real,allocatable :: vels(:,:)
 integer          :: fd_order

 logical :: store_wave_flag

 type(wave_boundary_param)::wave_boundary

 call start_mpi

! fd_order=22

 !! Read input parameters
 !! ---------------------
 call readparamfile(parfile, par)

 fd_order = par%fd_order
 par%lv = int(par%fd_order/2)
  
 !! PML setting
 !! -----------
    call init_pml(par%nx,par%nz,par%npml)

 !! Read acquisition geometry data
 !! ------------------------------
 call readcoordfile(par%coordfile,coord)

 !! Memory allocations
 !! ------------------
 allocate(s(par%nt),fs(nx_pml))
 allocate(vp(nz_pml,nx_pml),vs(nz_pml,nx_pml),den(nz_pml,nx_pml))
 allocate(taper(par%n_taper))
 !allocate(wave(1:par%nz,1:par%nx,1:par%nt))

 
 !! Read velocity model
 !! -------------------
 call readvpfile(par,vp,npml,nx_pml,nz_pml,fd_order,par%free_surface)
 call readvsfile(par,vs,npml,nx_pml,nz_pml,fd_order,par%free_surface)

 !! Assume density model is constant
 do ix=1,nx_pml
 do iz=1,nz_pml
    den(iz,ix)=1.0
 enddo
 enddo

!! !! Read density model
!! call read_densityfile(par,c,den,npml,nx_pml,nz_pml)

 !! Determine maximum and minimum velocities
 !! ----------------------------------------
 par%cmin=min_value(vp(iz1:iz2,ix1:ix2),par%nz,par%nx)
 par%cmax=max_value(vp(iz1:iz2,ix1:ix2),par%nz,par%nx)
 par%denmin=min_value(den(iz1:iz2,ix1:ix2),par%nz,par%nx)
 par%denmax=max_value(den(iz1:iz2,ix1:ix2),par%nz,par%nx)

 !! Setup PML damping coefficient
 !! -----------------------------
    call setup_pml(par%dx,par%cmin)

 !! Setup Ricker source wavelet 
 !! ---------------------------
 call getsource(par,s)
 !if(rank.eq.0) then
 !   call write_binfile('source.H@',s,par%nt)
 !   output=trim('source.H')
 !   call writehead(output,n1=par%nt,n2=1,n3=1,d1=par%dt,d2=1.0,d3=1.0,&
 !                  o1=0.0,o2=0.0,o3=0.0)
 !endif

 if(rank.eq.0) then
        write(*,*) 'free surface:', par%free_surface 
        write(*,*) 'npml: ',par%npml, 'nz_pml: ',nz_pml, 'nx_pml: ',nx_pml
        write(*,*) 'ns: ', coord%ns, 'ng: ',coord%ngmax
        write(*,*) "Minimum Velocity in the input model:", par%cmin
        write(*,*) "Maximum Velocity in the input model:", par%cmax
        write(*,*) "Minimum Density in the input model:", par%denmin
        write(*,*) "Maximum Density in the input model:", par%denmax
        if (par%cmax*par%dt/par%dx > 0.6) then
                write(*,*) "Stability Condition Violated !!!"
                write(*,*) "Terminating the Program"
                call flush(6)
                goto 999 
        else 
                write(*,*) "Stability Condition Satisfied"
                call flush(6)
        endif
 endif

 call readParFile(parfile,'FIRST_SHOT',par%first_shot,1)
 call readParFile(parfile,'LAST_SHOT',par%last_shot,coord%ns)

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call get_assigned(par%first_shot,par%last_shot,is1,is2)


 ! we are not going to use these arrays
 allocate(wave_pux(par%nz,par%nx,1),wave_puz(par%nz,par%nx,1),wave_pwx(par%nz,par%nx,1),wave_pwz(par%nz,par%nx,1))
 wave_pux = 0.
 wave_puz = 0.
 wave_pwx = 0.
 wave_pwz = 0.

 if(rank.eq.0) write(*,*) "Doing elastic FD modeling"
 !! Process all the shots
 !! ---------------------
 do is=is1,is2,par%skipshot+1

    !write(*,*) "Process ",rank," shot ", is
    call flush(6)
    allocate(seis_true_u(par%nt,coord%ng(is)),seis_true_w(par%nt,coord%ng(is)))

    seis_true_u=0.0; seis_true_w=0.0

    store_wave_flag = .false.

    call modeling_elastic_kernel(is,par,coord,s,vp,vs,den,par%free_surface,fd_order,nx_pml,nz_pml,npml,damp_global,&
         wave_pux,wave_puz,wave_pwx,wave_pwz,wave_boundary,seis_true_u,seis_true_w,store_wave_flag)

    !! Save the outputs
    !! ----------------
    call filename(output,par%csg_out,is,'.uH@')
    call write_binfile(output,seis_true_u(1:par%nt,1:coord%ng(is)),par%nt,coord%ng(is))
    call filename(output,par%csg_out,is,'.wH@')
    call write_binfile(output,seis_true_w(1:par%nt,1:coord%ng(is)),par%nt,coord%ng(is))

    deallocate(seis_true_u,seis_true_w)
 enddo

 print *, 'MODELING final', 'par%lv', par%lv


 999 continue
 deallocate(vp,vs,s,fs,taper,den)
 deallocate(wave_pux, wave_puz, wave_pwx, wave_pwz) 
 !deallocate(wave)
 
 
 call stop_mpi

end subroutine modeling_elastic
!---------------------------------

end module modeling_elastic_main

