
     use smooth_slow

     implicit none

     real,allocatable::gk0(:,:), gk1(:,:)
     integer,allocatable::wb_depth(:,:)

     character *256::wb_depth_file, file_vp

     integer::nx, nz, npml, ii_bi
     integer::nx_pml, nz_pml
     integer::horsmt, versmt
     integer::i,j,ix

     nz = 50
     nx = 471  
     nz_pml = 150
     nx_pml = 571 
     npml = 50

     horsmt = 12
     versmt = 4

     ii_bi = 1

     allocate(gk0(1:nz_pml,1:nx_pml), gk1(1:nz_pml,1:nx_pml))

     allocate(wb_depth(2,nx))

     gk0 = 0.0
     gk1 = 0.0

     wb_depth = 0

     file_vp='grad_gorgon_1.sg.taper.all.H@'

     open(12,file=file_vp,access='direct',form='unformatted',recl=ii_bi*nz_pml)
     do i=1,nx_pml
        read(12,rec=i) (gk0(j,i), j=1,nz_pml)
     enddo 
     close(12)

     call grdsmth2(gk0, horsmt, versmt)

     open(12,file='gk0_sm',access='direct',form='unformatted',recl=ii_bi*nz_pml)
     do i=1,nx_pml
        write(12,rec=i) (gk0(j,i), j=1,nz_pml)
     enddo
     close(12)

     wb_depth_file='/scratch2/guo103/wd_lab/Gorgon_coord_file_me/wb.dat'

     open(222, file=wb_depth_file)
     do ix=1,nx
        read(222, *) wb_depth(1,ix),wb_depth(2,ix)
     enddo
     close(222)

       !pad_top+1:pad_top+par%nz,par%npml+1:par%npml+par%nx

     do ix=1,nx
        gk0(1:wb_depth(2,ix)+npml, npml+ix) = 0.0
     enddo

     open(12,file='gk0_sm_wb',access='direct',form='unformatted',recl=ii_bi*nz_pml)
     do i=1,nx_pml
        write(12,rec=i) (gk0(j,i), j=1,nz_pml)
     enddo
     close(12)

     end
