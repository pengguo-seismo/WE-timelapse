! module for computing density from velocity
!
module modeling

implicit none

contains

subroutine compute_density(variable_density, c, den, nx, nz)

integer, intent(in)  :: nx, nz, variable_density
real,    intent(in)  :: c(:,:)
real,    intent(out) :: den(:,:)
integer              :: ix, iz
real                 :: alpha, beta

if (variable_density == 0) then
  den = 1.0
else
  do ix=1,nx
    do iz=1,nz

      den(iz,ix) = 0.23*c(iz,ix)**0.25

    enddo
  enddo
endif

end subroutine compute_density
!------------------------------------------------------------------------
subroutine compute_density_wb(variable_density, c, den, nx, nz, npml, izwb)

integer, intent(in)  :: nx, nz, variable_density, npml, izwb(:)
real,    intent(in)  :: c(:,:)
real,    intent(out) :: den(:,:)
integer              :: ix, iz, npml2
real                 :: alpha, beta, factor

factor = 0.35
!factor = 0.23
den = 1.0
npml2 = 2*npml
if (variable_density == 1) then
  do ix=1,npml
    do iz=izwb(1)+npml+1,nz
      den(iz,ix) = factor*c(iz,ix)**0.25
    enddo
  enddo
  do ix=npml+1,nx-npml
    do iz=izwb(ix-npml)+npml+1,nz
      den(iz,ix) = factor*c(iz,ix)**0.25
    enddo
  enddo
  do ix=nx-npml+1,nx
    do iz=izwb(nx-npml2)+npml+1,nz
      den(iz,ix) = factor*c(iz,ix)**0.25
    enddo
  enddo
endif

end subroutine compute_density_wb

end module modeling

