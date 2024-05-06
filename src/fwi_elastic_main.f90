module fwi_elastic_main

 use modeling
 use pml
 use source
 use fwi_elastic_subroutine
 use datatype
 use io
 use math
 use mmi_mpi
 use parser
 use smooth_slow

 implicit none

 type(param),       private              :: par
 type(acquisition), private              :: coord
 logical,           private              :: store_boundary, message
 integer,           private              :: ix,iz,it,isx,isz,igx,igz,is,ig,i1,i2,&
                                            is1,is2,iter,nsearch,nsearch_max,is_pre
 real,              private              :: d1,d2,factor1,factor2
 integer,           private, allocatable :: fs(:), izwb(:)

 double precision                        :: res0,res1,res2,gg,dgg,res_tmp, &
                                            res_shot,gg_process,gg_shot, &
                                            dgg_process,dgg_shot,res_process,res_process2,res,gama,gama1,gama2

 real,              private, allocatable :: s(:),c(:,:),kappa(:,:),bc_top(:,:,:),bc_bottom(:,:,:),bc_left(:,:,:),  &
                                            bc_right(:,:,:),seis_true(:,:),seis(:,:),p_nt(:,:),p_nt_1(:,:),&
                                            energy_temp(:,:),energy(:,:),mig_temp(:,:),  &
                                            taper(:),mig_sum(:,:),energy_sum(:,:),gk(:,:),dk(:,:), &
                                            gk_process(:,:),gk_shot(:,:),energy_process(:,:),energy_shot(:,:), &
                                            seis_temp(:,:),den(:,:)

 contains

!!--------------------------------------
subroutine fwi_elastic(parfile)

 character(len=*), intent(in) :: parfile
 character(len=256)           :: file_name

 real             :: min_bound,max_bound,inner_loop,res_pre
 integer          :: fd_order,pad_top,count_num,max_offset
 real,allocatable :: vp(:,:),vs(:,:),vs_tmp(:,:)
 real,allocatable :: grad_tmp(:,:),res_all(:)
 real,allocatable :: gk_old(:,:),dk_dir(:,:)
 integer, allocatable :: wb_depth(:,:)

 double precision :: res_1st,res_2nd,res_3rd,res_4th

 call start_mpi

! fd_order=22
! pad_top=(fd_order-20)/2+1
 !! Read input parameters
 !! ---------------------
 call readparamfile(parfile, par)
 
 fd_order = par%fd_order

 pad_top = par%npml

 par%lv = int(par%fd_order/2)

 !! PML setting
 !! -----------
    call init_pml(par%nx,par%nz,par%npml)

 !! Read acquisition geometry data
 !! ------------------------------
 call readcoordfile(par%coordfile, coord)

 !! Memory allocations
 !! ------------------
 allocate(s(par%nt),fs(nx_pml))
 allocate(vp(nz_pml,nx_pml),vs(nz_pml,nx_pml),den(nz_pml,nx_pml))
 allocate(vs_tmp(nz_pml,nx_pml))

 allocate(energy(nz_pml,nx_pml))

 allocate(gk_process(nz_pml,nx_pml),gk_shot(nz_pml,nx_pml))
 allocate(gk(nz_pml,nx_pml),dk(nz_pml,nx_pml))
 allocate(energy_process(nz_pml,nx_pml),energy_shot(nz_pml,nx_pml))


 !! Read velocity model
 !! -------------------
 call readvpfile(par,vp,npml,nx_pml,nz_pml,fd_order,par%free_surface)
 call readvsfile(par,vs,npml,nx_pml,nz_pml,fd_order,par%free_surface)

 print *, 'test0 ', 'vp velocities safely read'
 print *, 'test0 ', 'vs velocities safely read'

 !! Assume constant velocity
 !! ------------------------
 do ix=1,nx_pml
 do iz=1,nz_pml
    den(iz,ix)=1.0
 enddo
 enddo
 ! !! Read density model
 ! call read_densityfile(par,c,den,npml,nx_pml,nz_pml)

 !! Determine minimum velocity
 !! --------------------------
 !! change as you wish 
 !! we need to change it as read-in parameters
 min_bound=0.
 max_bound=4000.0
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

 print *, 'test1 ', 'source wavelet safely read'

 if(rank.eq.0) then
    write(*,*) 'Free_surface:', par%free_surface
    write(*,*) 'npml: ',par%npml, 'nz_pml: ',nz_pml, 'nx_pml: ',nx_pml
    write(*,*) 'ns: ', coord%ns, 'ng: ',coord%ngmax
    write(*,*) 'dt: ', par%dt, 'dx: ', par%dx
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

 print *, 'par%first_shot', par%first_shot
 print *, 'par%last_shot', par%last_shot

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call get_assigned(par%first_shot,par%last_shot,is1,is2)

 print *, 'is1', is1, 'is2', is2

 if(rank.eq.0) then
    write(*,*) "Doing full waveform inversion for S velocity"
    call flush(6)
    file_name=trim(trim(par%logfile)//'.log')
    open(99,file=file_name)
    print *, 'file_name', file_name
 endif

 ! ****************************************!
 allocate(res_all(par%itermax))
 res_all=0.0

 ! step length
 factor1=20.0

 max_offset = int(3.0*par%nz*par%dx)
 !max_offset = 90

 if (par%nlcg) then 
    allocate(gk_old(nz_pml,nx_pml))
    gk_old = 0.
    allocate(dk_dir(nz_pml,nx_pml))
    dk_dir = 0.
 endif

 if (par%cond_wb) then
    allocate(wb_depth(2,par%nx))
    wb_depth = 0.
 endif

 dk = 0.

 inner_loop=1
 do iter=1,par%itermax

    par%iter=iter
   
    print *, 'nx_pml', nx_pml, 'nz_pml', nz_pml 

    gk_process = 0.
    energy_process = 0.
    gk = 0.
    energy = 0.

!    !$omp parallel do private(ix,iz)
!    do ix=1,nx_pml
!    do iz=1,nz_pml
!       gk_process(iz,ix)=0.0; energy_process(iz,ix)=0.0
!       gk(iz,ix)=0.0; 
!       energy(iz,ix)=0.0
!    enddo
!    enddo
!    !$omp end parallel do

    res_process=0.D0; res0=0.D0
    !! Compute the gradient for all the shots
    !! -------------------------------------- 
    do is=is1,is2,par%skipshot+1

       print *, 'is1', is1, 'is2', is2, 'inside is loop fwi_elastic_inversion'
       !    compute_gradient_elastic(coord,par,is,fd_order,max_offset,s,vp,vs,den,nx_pml,nz_pml,npml,damp,energy,gk,res)
       call compute_gradient_elastic(coord,par,is,fd_order,max_offset,s,vp,vs,den,nx_pml,nz_pml,npml,damp_global,energy_shot,gk_shot,res_shot)
       print *, 'compute_gradient_elastic finished for shot ', is, '_rank_', rank
       !$omp parallel do private(ix,iz)
       do ix=1,nx_pml
       do iz=1,nz_pml
          gk_process(iz,ix)=gk_process(iz,ix)+gk_shot(iz,ix)
          energy_process(iz,ix)=energy_process(iz,ix)+energy_shot(iz,ix)
       enddo
       enddo
       !$omp end parallel do
       res_process=res_process+res_shot
    enddo       


    !!! debug only
!    call MPI_Barrier(MPI_COMM_WORLD,ierr)
!    stop

    print *, 'fwi_elastic_inversion test1'

    !! Reduce for all the shots
    !! ------------------------
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    print *, 'loops for compute_gradient_elastic finished'

    call MPI_Allreduce(res_process,res0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(gk_process,gk,nz_pml*nx_pml,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(energy_process,energy,nz_pml*nx_pml,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

    print *, 'fwi_elastic_inversion test2'

    !! record the residual
    !! -------------------
    if(inner_loop==1)then
       res_1st=res0
    elseif(inner_loop==2)then
       res_2nd=res0
    elseif(inner_loop.GT.2)then
       res_3rd=res0
    endif

    if(rank.eq.0) then
      call filename(output,trim(par%gradfile),iter,'.raw.H@')
      call write_binfile(output,gk(pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx),par%nz,par%nx)

      call filename(output,trim(par%gradfile),iter,'.raw.all.H@')
      call write_binfile(output,gk(1:nz_pml,1:nx_pml),nz_pml,nx_pml)

      call filename(output,trim(par%gradfile),iter,'.illu.H@')
      call write_binfile(output,energy(pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx),par%nz,par%nx)

    endif

!    dk = gk

    !!!
    if (par%grad_illu==1) then 

    !! illumination of velocity and density gradient
    do ix=1,par%nx
    do iz=1,par%nz
    !   dk(iz+pad_top,ix+npml)=gk(iz+pad_top,ix+npml)/(energy(iz+pad_top,ix+npml)+0.001*maxval(abs(energy)))
       gk(iz+pad_top,ix+npml)=gk(iz+pad_top,ix+npml)/(energy(iz+pad_top,ix+npml)+0.001*maxval(abs(energy)))
    enddo
    enddo
    !dk(1:npml+11,:)=0  

    endif !!
    !!! 

    print *, 'fwi_elastic_inversion test3, ', par%grad_illu

    if(rank.eq.0) then
      call filename(output,trim(par%gradfile),iter,'.raw.illu.all.H@')
      call write_binfile(output,gk(1:nz_pml,1:nx_pml),nz_pml,nx_pml)
    endif

    if (par%taper_source) then
      
       !! Set up source position
       !! ----------------------
       if(par%free_surface==0)then
          isx=npml+int(coord%xs(is)/par%dx)
          isz=npml+int(coord%zs(is)/par%dx)
       else
          isx=npml+int(coord%xs(is)/par%dx)
          isz=pad_top+int(coord%zs(is)/par%dx)
       endif  
    
       gk(isz, isx) = 0.0   

    endif

    if (par%taper_geophone) then

       do ig=1,coord%ng(is)
         if(par%free_surface==0)then
            igx=npml+int(coord%xg(is,ig)/par%dx)
            igz=npml+int(coord%zg(is,ig)/par%dx)
         else
            igx=npml+int(coord%xg(is,ig)/par%dx)
            igz=pad_top+int(coord%zg(is,ig)/par%dx)
         endif
         gk(igz,igx) = 0.0
         gk(igz,igx) = 0.0
       enddo

    endif

    if(rank.eq.0) then
      call filename(output,trim(par%gradfile),iter,'.sg.taper.all.H@')
      call write_binfile(output,gk(1:nz_pml,1:nx_pml),nz_pml,nx_pml)
    endif

    if(par%smooth_type == 0) then 
       allocate(grad_tmp(par%nz,par%nx))
       grad_tmp=0.0
  !  call box_smooth(dk(pad_top+1:nz_pml-npml,npml+1:npml+par%nx),1,3,2,par%nz,par%nx,grad_tmp)
  !  dk(pad_top+1:nz_pml-npml,npml+1:nx_pml-npml)=grad_tmp
       call box_smooth(gk(pad_top+1:nz_pml-npml,npml+1:npml+par%nx),1,3,2,par%nz,par%nx,grad_tmp)
       gk(pad_top+1:nz_pml-npml,npml+1:nx_pml-npml)=grad_tmp
       deallocate(grad_tmp)

    ! Gaussian smoothing
    elseif(par%smooth_type == 1) then

       call grdsmth2(gk, par%horsmt, par%versmt)

    endif

    if(rank.eq.0) then
      call filename(output,trim(par%gradfile),iter,'.sg.taper.sm.all.H@')
      call write_binfile(output,gk(1:nz_pml,1:nx_pml),nz_pml,nx_pml)
    endif

    if (par%nlcg) then 

      !! use this one
      call conjugate_gradient_new(par,gk,gk_old,dk,dk_dir)

      gk_old = gk

    else

      dk = gk

    endif

    if (par%cond_wb) then 

       open(222, file=par%wb_depth_file)
       do ix=1,par%nx
         read(222, *) wb_depth(1,ix),wb_depth(2,ix)
       enddo
       close(222)

       !pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx

       do ix=1,par%nx
           dk(1:wb_depth(2,ix)+par%npml, par%npml+ix) = 0.0
       enddo

    endif

    !!!! for debug
    if (par%grad_inv_phase == 1) then 
      dk = - dk
    endif

    !! residual of previous iteration
    res_pre=res0

    print *, 'fwi_elastic_inversion test4'

    if (rank.eq.0) then
       print*,'********************************************************'
       write(*,*) 'Iter= ',par%iter,' Residual: ',res0,'loop=',int(inner_loop)
       call flush(6)
       write(99,*) 'Iter= ',par%iter,' Residual: ',res0
       call flush(99)
       print*,'max_offset = ',max_offset
    !   print*,'********************************************************'
    endif

    !! Normalize the gradient
    !! ----------------------
    dk=dk/maxval(abs(dk))
  
    if(rank.eq.0) then
      call filename(output,trim(par%gradfile),iter,'.H@')
      call write_binfile(output,dk(pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx),par%nz,par%nx)
      call filename(output,trim(par%gradfile),iter,'.all.H@')
      call write_binfile(output,dk(1:nz_pml,1:nx_pml),nz_pml,nx_pml)
      print *, 'nx_pml', nx_pml, 'nz_pml', nz_pml
    endif

    !! Trail updates 1
    !! ---------------------------------
    !$omp parallel do private(ix,iz)
    do ix=1,nx_pml
    do iz=1,nz_pml
       vs_tmp(iz,ix)=0.0
    enddo
    enddo
    !$omp end parallel do
 
    !$omp parallel do private(ix,iz)
    do ix=1,nx_pml
    do iz=1,nz_pml
       vs_tmp(iz,ix)=vs(iz,ix)-factor1*dk(iz,ix)
       if(vs_tmp(iz,ix).LT.min_bound) vs_tmp(iz,ix)=min_bound
       if(vs_tmp(iz,ix).GT.max_bound) vs_tmp(iz,ix)=max_bound
    enddo
    enddo
    !$omp end parallel do


    ! Extrapolate velocity in PML regions
    do ix=1,npml
       vs_tmp(pad_top+1:pad_top+par%nz,ix) = vs_tmp(pad_top+1:pad_top+par%nz,npml+1)
       vs_tmp(pad_top+1:pad_top+par%nz,nx_pml-npml+ix) = vs_tmp(pad_top+1:pad_top+par%nz,nx_pml-npml)
    enddo
    do iz=1,npml
       vs_tmp(nz_pml-npml+iz,:) = vs_tmp(nz_pml-npml,:)
    enddo
    do iz=1,pad_top
       vs_tmp(iz,:) = vs_tmp(pad_top+1,:)
    enddo


    dgg_process=0.D0; res1=0.D0
    do is=is1,is2,par%skipshot+1
       call analytical_step_length_elastic(coord,par,is,fd_order,max_offset,s,vp,vs_tmp,den,par%free_surface,nx_pml,nz_pml,npml,damp_global,.true.,dgg_shot)
       dgg_process=dgg_process+dgg_shot
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(dgg_process,res1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    !! Trail updates 2
    !! ----------------
    if(res1>res0)then  

      count_num=0
      do while(res1>res0)

        factor1=factor1*0.5
        !$omp parallel do private(ix,iz)
         do ix=1,nx_pml
         do iz=1,nz_pml
            vs_tmp(iz,ix)=0.0
         enddo
         enddo
        !$omp end parallel do
       
        !$omp parallel do private(ix,iz)
         do ix=1,nx_pml
         do iz=1,nz_pml
            vs_tmp(iz,ix)=vs(iz,ix)-factor1*dk(iz,ix)
            if(vs_tmp(iz,ix).LT.min_bound) vs_tmp(iz,ix)=min_bound
            if(vs_tmp(iz,ix).GT.max_bound) vs_tmp(iz,ix)=max_bound
         enddo
         enddo
        !$omp end parallel do

        ! Extrapolate velocity in PML regions
        do ix=1,npml
           vs_tmp(pad_top+1:pad_top+par%nz,ix) = vs_tmp(pad_top+1:pad_top+par%nz,npml+1)
           vs_tmp(pad_top+1:pad_top+par%nz,nx_pml-npml+ix) = vs_tmp(pad_top+1:pad_top+par%nz,nx_pml-npml)
        enddo
        do iz=1,npml
           vs_tmp(nz_pml-npml+iz,:) = vs_tmp(nz_pml-npml,:)
        enddo
        do iz=1,pad_top
           vs_tmp(iz,:) = vs_tmp(pad_top+1,:)
        enddo

        dgg_process=0.D0; res1=0.D0
        do is=is1,is2,par%skipshot+1
           call analytical_step_length_elastic(coord,par,is,fd_order,max_offset,s,vp,vs_tmp,den,par%free_surface,nx_pml,nz_pml,npml,damp_global,.true.,dgg_shot)
           dgg_process=dgg_process+dgg_shot
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(dgg_process,res1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        if(rank==0) print*,'res1= res0= factor1=',res1,res0,factor1, iter
        count_num=count_num+1 

        if(count_num.GE.5)then
                print *, 'not good count_num', count_num
                exit
        endif

     end do ! end of while do

     if(count_num.LE.5)then

       !! Update and save the image
       !! -------------------------
       !$omp parallel do private(ix,iz)
       do ix=1,nx_pml
       do iz=1,nz_pml
          vs(iz,ix)=vs_tmp(iz,ix)
       enddo
       enddo
      !$omp end parallel do
       if(rank.eq.0) then
          call filename(output,trim(par%migfile),iter,'.H@')
          call write_binfile(output,vs(pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx),par%nz,par%nx)
       endif !! endif if(rank.eq.0)
     
     endif !! (count_num.LE.5) 

   else
      
      count_num=0
      if(rank==0) print*,'direct update',res1,res0,factor1, iter
      !! Update and save the image
      !! -------------------------
      !$omp parallel do private(ix,iz)
      do ix=1,nx_pml
      do iz=1,nz_pml
         vs(iz,ix)=vs_tmp(iz,ix)
      enddo
      enddo
      !$omp end parallel do

      if(rank.eq.0) then
         call filename(output,trim(par%migfile),iter,'.H@')
         call write_binfile(output,vs(pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx),par%nz,par%nx)
      endif

   endif   

   inner_loop=inner_loop+1

!   if(res0.eq.0.0)then
!     print*,'No more offset available, problem finished!'
!     exit
!   endif

 enddo      !! End of iterations loop

 print *, 'FWI final', 'par%lv', par%lv

 call MPI_Barrier(MPI_COMM_WORLD,ierr)

 if(rank.eq.0) close(99)    !! Close the log file

 999 continue

 deallocate(s,fs)
 deallocate(vp,vs,den)
 deallocate(vs_tmp)
 deallocate(energy)
 deallocate(gk_process,gk_shot)
 deallocate(gk,dk)
 deallocate(energy_process,energy_shot)

 if (par%nlcg) deallocate(gk_old,dk_dir)
 if (par%cond_wb) deallocate(wb_depth)

 call stop_mpi

end subroutine fwi_elastic


end module fwi_elastic_main
