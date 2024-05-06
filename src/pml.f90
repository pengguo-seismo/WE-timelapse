! module for defining PML arrays
!
module pml

implicit none
integer :: npml,nx_pml,ny_pml,nz_pml,iz1,iz2,ix1,ix2,iy1,iy2
real, allocatable :: damp(:),damp_global(:,:),damp_global_3d(:,:,:)

contains

subroutine init_pml(nx,nz,n_pml)

integer, intent(in)          :: nx, nz, n_pml

npml = n_pml
nx_pml = nx+2*npml
nz_pml = nz+2*npml
ix1 = npml+1
ix2 = npml+nx
iz1 = npml+1
iz2 = npml+nz
allocate(damp(npml))
allocate(damp_global(nz_pml,nx_pml))

end subroutine init_pml


!------------------------------------

!-------------------------------------------------------------------
subroutine setup_pml(dx, cmin)

real, intent(in) :: dx, cmin
real             :: a, xa, kappa
integer          :: ix, iz

damp_global = 0.0
a = (npml-1)*dx
kappa = 3.0*cmin*log(10000000.0)/(2.0*a)
do ix=1,npml
  xa = real(ix-1)*dx/a
  damp(ix) = kappa*xa*xa
enddo
do ix=1,npml
  do iz=1,nz_pml
    damp_global(iz,npml-ix+1) = damp(ix)
    damp_global(iz,nx_pml+ix-npml) = damp(ix)
  enddo
enddo
do iz=1,npml
  do ix=1+(npml-(iz-1)),nx_pml-(npml-(iz-1))
    damp_global(npml-iz+1,ix) = damp(iz)
    damp_global(nz_pml+iz-npml,ix) = damp(iz)
  enddo
enddo

end subroutine setup_pml
!------------------------------------



end module pml

