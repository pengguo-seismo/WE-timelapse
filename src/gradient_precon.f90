module gradient_precon

 use global
 use string
 use math
 use string 
 use io
 use datatype
 use pml
 use mmi_mpi

 implicit none

! include 'fftw3.f'

 contains

!!-----------------------------------------------
subroutine precondition_gradient(par,gk,illum,dk)

 type(param),       intent(in)    :: par
 real,              intent(in)    :: gk(nz_pml,nx_pml),illum(nz_pml,nx_pml)
 real,              intent(inout) :: dk(nz_pml,nx_pml)

 integer                          :: i1, i2
 real                             :: eps

 if(par%pre.eq.1) then

    eps=max_value(illum(par%npml+1:par%npml+par%nz,&
                        par%npml+1:par%npml+par%nx),par%nz,par%nx)
    eps=0.01*eps

    !$OMP PARALLEL DO private(i1,i2)
    do i2=par%npml+1,par%npml+par%nx
    do i1=par%npml+1,par%npml+par%nz
        dk(i1,i2)=gk(i1,i2)/(illum(i1,i2)+eps)
    enddo
    enddo
    !$OMP END PARALLEL DO

 else
    dk=gk

 endif

 if(par%highpass.eq.1) call highpassfilter(par%nz,par%nx,dk(par%npml+1:par%npml+par%nz,  &
                                           par%npml+1:par%npml+par%nx))


end subroutine precondition_gradient

!!----------------------------------------------
subroutine conjugate_gradient(par,gk,gk1,dk,dk1)

 type(param),       intent(inout) :: par
 real,              intent(in)    :: gk(nz_pml,nx_pml),gk1(nz_pml,nx_pml)
 real,              intent(inout) :: dk(nz_pml,nx_pml),dk1(nz_pml,nx_pml)

 integer                          :: i1,i2
 double precision                 :: gg,dgg

 if(par%iter.gt.1) then

    gg=0.D0; dgg=0.D0
    !$omp parallel do private(i1,i2) reduction(+:gg,dgg)
    do i2=par%npml+1,par%npml+par%nx
    do i1=par%npml+1,par%npml+par%nz
        gg=gg+gk(i1,i2)*dk(i1,i2)
        dgg=dgg+gk1(i1,i2)*dk1(i1,i2)
    enddo
    enddo
    !$omp end parallel do

    if(dgg.ne.0.D0) par%beta=gg/dgg

    if(mod(par%iter,5)==0) par%beta=0.D0

    !$omp parallel do private(i1,i2)
    do i2=1,nx_pml
    do i1=1,nz_pml
        dk(i1,i2)=dk(i1,i2)+par%beta*dk1(i1,i2)
    enddo
    enddo
    !$omp end parallel do
 endif


end subroutine conjugate_gradient
!!-------------------------------

!!----------------------------------------------
subroutine conjugate_gradient_new(par,gk,gk_old,dk,dk_dir)

 type(param),       intent(inout) :: par
 real,              intent(in)    :: gk(nz_pml,nx_pml),gk_old(nz_pml,nx_pml)
 real,              intent(inout) :: dk(nz_pml,nx_pml)
 real,              intent(inout) :: dk_dir(nz_pml,nx_pml)

 integer                          :: i1,i2
 double precision                 :: gg,dgg

! if(par%iter.gt.1) then

    gg=0.D0; dgg=0.D0
    !$omp parallel do private(i1,i2) reduction(+:gg,dgg)
    do i2=par%npml+1,par%npml+par%nx
    do i1=par%npml+1,par%npml+par%nz
        gg=gg+gk(i1,i2)*gk(i1,i2)
        dgg=dgg+gk_old(i1,i2)*gk_old(i1,i2)
   !     gg=gg+dk(i1,i2)*dk(i1,i2)
   !     dgg=dgg+dk_old(i1,i2)*dk_old(i1,i2)
    enddo
    enddo
    !$omp end parallel do

    if(dgg.ne.0.D0) par%beta=gg/dgg

    if(mod(par%iter,5)==0) par%beta=0.D0

    if(par%beta<0.0) par%beta=0.0

    if(par%iter.eq.1) par%beta=0.0

    !$omp parallel do private(i1,i2)
    do i2=1,nx_pml
    do i1=1,nz_pml
    !   dk(i1,i2)=gk(i1,i2)+par%beta*dk(i1,i2)
        dk_dir(i1,i2)=-gk(i1,i2)+par%beta*dk_dir(i1,i2)
    enddo
    enddo
    !$omp end parallel do
! endif

    dk = -dk_dir

end subroutine conjugate_gradient_new
!!-------------------------------

!!this one does not work
subroutine conjugate_gradient_new0(par,dk,gk,gk_old)

 type(param),       intent(inout) :: par
 real,              intent(in)    :: gk(nz_pml,nx_pml),gk_old(nz_pml,nx_pml)
 real,              intent(inout) :: dk(nz_pml,nx_pml)

 integer                          :: i1,i2
 double precision                 :: gg,dgg

! if(par%iter.gt.1) then

    gg=0.D0; dgg=0.D0
    !$omp parallel do private(i1,i2) reduction(+:gg,dgg)
    do i2=par%npml+1,par%npml+par%nx
    do i1=par%npml+1,par%npml+par%nz
        gg=gg+gk(i1,i2)*gk(i1,i2)
        dgg=dgg+gk_old(i1,i2)*gk_old(i1,i2)
    enddo
    enddo
    !$omp end parallel do

    if(dgg.ne.0.D0) par%beta=gg/dgg

    if(mod(par%iter,5)==0) par%beta=0.D0

    if(par%beta<0.0) par%beta=0.0

    if(par%iter.eq.1) par%beta=0.0

    !$omp parallel do private(i1,i2)
    do i2=1,nx_pml
    do i1=1,nz_pml
    !   dk(i1,i2)=gk(i1,i2)+par%beta*dk(i1,i2)
        dk(i1,i2)=-gk(i1,i2)+par%beta*dk(i1,i2)
    enddo
    enddo
    !$omp end parallel do
! endif

    dk = -dk

end subroutine conjugate_gradient_new0


subroutine highpassfilter(nz,nx,mig)

 implicit none
 integer, intent(in)    :: nz,nx
 real,    intent(inout) :: mig(:,:)
 integer                :: i,j,k,n1,n2
 real, allocatable      :: mig1(:,:),mig2(:,:)

 n1=20
 n2=2
 allocate(mig1(nz,nx))
 allocate(mig2(nz,nx))
 mig2=mig
 do k=1,n1
    do j=2,nz-1
    do i=1,nx
        mig1(j,i)=0.25*mig(j-1,i)+0.5*mig(j,i)+0.25*mig(j+1,i)
    enddo
    enddo

    do i=1,nx
        mig1(1,i)=0.75*mig(1,i)+0.25*mig(2,i)
        mig1(nz,i)=0.75*mig(nz,i)+0.25*mig(nz-1,i)
    enddo

    mig=mig1
 enddo

 do k=1,n2
    do j=1,nz
    do i=2,nx-1
        mig1(j,i)=0.25*mig(j,i-1)+0.5*mig(j,i)+0.25*mig(j,i+1)
    enddo
    enddo

   do j=1,nz
        mig1(j,1)=0.75*mig(j,1)+0.25*mig(j,2)
        mig1(j,nz)=0.75*mig(j,nz)+0.25*mig(j,nz-1)
   enddo
   mig=mig1
 enddo

 mig=mig2-mig1

 deallocate(mig1,mig2)

end subroutine highpassfilter

subroutine hanning_taper_left(n1,n2,taper)

implicit none

integer, intent(in)  :: n1,n2
real,    intent(out) :: taper(n1)

integer :: kk,i
real,allocatable :: taper1(:)

allocate(taper1(n1))
taper=0.0

kk=1
do i=0,n1-1
   taper1(kk)=cos((i*1.0/(n1-1))*(3.1415926/2))**2
   kk=kk+1
enddo

do i=1,n1
   taper(i)=taper1(n1-i+1)
enddo

deallocate(taper1)

end subroutine hanning_taper_left


end module gradient_precon 
