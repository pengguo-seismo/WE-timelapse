
!! module for recording wave boundary values
!! this is used for source wavefield reconstruction during FWI

module wave_bound

contains      

      subroutine record_wave_boundary_uw(it, nt, u, w,&
                                         u_bl, u_br, u_bt, u_bb,&
                                         w_bl, w_br, w_bt, w_bb,&
                                         nx_pml, nz_pml, nx, nz, lv, npml)

      
       implicit none

       integer::nx_pml, nz_pml, nx, nz, lv,npml,nt
       integer::it, i, j
       
       real, intent(inout) ::u_bl(:,:,:), u_br(:,:,:), u_bt(:,:,:), u_bb(:,:,:), & 
                             w_bl(:,:,:), w_br(:,:,:), w_bt(:,:,:), w_bb(:,:,:)

       real, intent(in) :: u(:,:), w(:,:) 

       !! left & right
       do i=1,nz
          do j=1,lv
             u_bl(i,j,it)=u(i+npml,npml+j)
             w_bl(i,j,it)=w(i+npml,npml+j)
             u_br(i,j,it)=u(i+npml,nx_pml-npml-lv+j)
             w_br(i,j,it)=w(i+npml,nx_pml-npml-lv+j)

          !   u_bl(i,j,it)=u(i+npml,npml+j+1)
          !   w_bl(i,j,it)=w(i+npml,npml+j+1)
          !   u_br(i,j,it)=u(i+npml,nx_pml-npml-lv+j-1)
          !   w_br(i,j,it)=w(i+npml,nx_pml-npml-lv+j-1)
          enddo
       enddo       

       !! top & bottom
       do i=1,nx
          do j=1,lv
             u_bt(i,j,it)=u(npml+j,i+npml)
             w_bt(i,j,it)=w(npml+j,i+npml)
             u_bb(i,j,it)=u(nz_pml-npml-lv+j,i+npml)
             w_bb(i,j,it)=w(nz_pml-npml-lv+j,i+npml)

          !   u_bt(i,j,it)=u(npml+j+1,i+npml)
          !   w_bt(i,j,it)=w(npml+j+1,i+npml)
          !   u_bb(i,j,it)=u(nz_pml-npml-lv+j-1,i+npml)
          !   w_bb(i,j,it)=w(nz_pml-npml-lv+j-1,i+npml)
         enddo
       enddo

      end 

      subroutine record_wave_boundary_stress(it, nt, xx, zz, xz, &
                                             xx_bl, xx_br, xx_bt, xx_bb,&
                                             zz_bl, zz_br, zz_bt, zz_bb,&
                                             xz_bl, xz_br, xz_bt, xz_bb,&
                                             nx_pml, nz_pml, nx, nz, lv, npml)


       implicit none

       integer::nx_pml, nz_pml, nx, nz, lv,npml,nt
       integer::it, i, j

       real, intent(inout) ::xx_bl(:,:,:), xx_br(:,:,:), xx_bt(:,:,:), xx_bb(:,:,:), &
                             zz_bl(:,:,:), zz_br(:,:,:), zz_bt(:,:,:), zz_bb(:,:,:), &
                             xz_bl(:,:,:), xz_br(:,:,:), xz_bt(:,:,:), xz_bb(:,:,:)

       real, intent(in) :: xx(:,:), zz(:,:), xz(:,:)

       !! left & right
       do i=1,nz
          do j=1,lv
             xx_bl(i,j,it)=xx(i+npml,npml+j)
             zz_bl(i,j,it)=zz(i+npml,npml+j)
             xz_bl(i,j,it)=xz(i+npml,npml+j)
             xx_br(i,j,it)=xx(i+npml,nx_pml-npml-lv+j)
             zz_br(i,j,it)=zz(i+npml,nx_pml-npml-lv+j)
             xz_br(i,j,it)=xz(i+npml,nx_pml-npml-lv+j)

          !   xx_bl(i,j,it)=xx(i+npml,npml+j+1)
          !   zz_bl(i,j,it)=zz(i+npml,npml+j+1)
          !   xz_bl(i,j,it)=xz(i+npml,npml+j+1)
          !   xx_br(i,j,it)=xx(i+npml,nx_pml-npml-lv+j-1)
          !   zz_br(i,j,it)=zz(i+npml,nx_pml-npml-lv+j-1)
          !   xz_br(i,j,it)=xz(i+npml,nx_pml-npml-lv+j-1)
          enddo
       enddo

       !! top & bottom
       do i=1,nx
          do j=1,lv
             xx_bt(i,j,it)=xx(npml+j,i+npml)
             zz_bt(i,j,it)=zz(npml+j,i+npml)
             xz_bt(i,j,it)=xz(npml+j,i+npml)
             xx_bb(i,j,it)=xx(nz_pml-npml-lv+j,i+npml)
             zz_bb(i,j,it)=zz(nz_pml-npml-lv+j,i+npml)
             xz_bb(i,j,it)=xz(nz_pml-npml-lv+j,i+npml)

          !   xx_bt(i,j,it)=xx(npml+j+1,i+npml)
          !   zz_bt(i,j,it)=zz(npml+j+1,i+npml)
          !   xz_bt(i,j,it)=xz(npml+j+1,i+npml)
          !   xx_bb(i,j,it)=xx(nz_pml-npml-lv+j-1,i+npml)
          !   zz_bb(i,j,it)=zz(nz_pml-npml-lv+j-1,i+npml)
          !   xz_bb(i,j,it)=xz(nz_pml-npml-lv+j-1,i+npml)
         enddo
       enddo

      end


end module wave_bound

