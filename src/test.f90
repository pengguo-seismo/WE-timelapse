
       use datatype
       use mmi_mpi
       use string

       implicit none
       real, allocatable::shot_virtual(:,:), seis_w(:,:)
       integer:: nt, ng, num_shot
       type(param):: par
       type(acquisition):: coord
       character is_str*4
       double precision::res
       integer::is, ig, it
       character*256::file_syn, csg_in


       num_shot = 32

       nt = 24000
!       ng = 

       par%coordfile='/scratch2/guo103/wd_lab/Gorgon_coord_file_me/coord_gorgon3924_3200.dat'

       csg_in = '/scratch2/guo103/wd_lab/sensi_syn_a/syn_0/fwi_iu_ob3_sm2/'

       call readcoordfile(par%coordfile, coord)  

       print *, 'ok1'

       res=0.

       do is = 1, num_shot

          allocate(shot_virtual(nt,coord%ng(is)))
          shot_virtual=0.0  

          allocate(seis_w(nt,coord%ng(is)))
          seis_w=0.0

          write(is_str,'(i4.4)')  is

          open(12,file='/scratch2/guo103/wd_lab/sensi_syn_a/syn_0/'//&
               'shot_virtual.bin'//is_str,access='direct', &
                form='unformatted', recl=nt)
           do ig=1, coord%ng(is)
              read(12,rec=ig) (shot_virtual(it,ig), it=1,nt)
           enddo 
          close(12)

          call filename(file_syn,csg_in,is,'.post.wH@')

         ! print *, is, trim(file_syn)

          open(12,file=file_syn,access='direct', &
                form='unformatted', recl=nt)
           do ig=1, coord%ng(is)
              read(12,rec=ig) (seis_w(it,ig), it=1,nt)
           enddo
          close(12)

          do ig=1,coord%ng(is)
             res = res + sum((shot_virtual(:,ig)*seis_w(:,ig)))
          enddo 

          deallocate(shot_virtual,seis_w)

       enddo 


       print *, sqrt(res)


       end

