module smooth_slow

contains


subroutine padding(A,nz,nx,win_z,win_x,B)

implicit none

real,    intent(in)  :: A(:,:)
integer, intent(in)  :: nz,nx,win_z,win_x
real,    intent(out) :: B(:,:)

integer :: ix,iz

  do ix=1,nx+2*win_x
  do iz=1,nz+2*win_z
     B(iz,ix)=0.0
  enddo
  enddo

  B(win_z+1:win_z+nz,win_x+1:win_x+nx)=A
  
  ! Extrapolate regions
  do ix=1,win_x
     B(win_z+1:win_z+nz,ix) = A(:,1)
     B(win_z+1:win_z+nz,win_x+nx+ix) = A(:,nx)
  enddo
  do iz=1,win_z
     B(iz,:) = B(win_z+1,:)
     B(win_z+nz+iz,:) = B(win_z+nz,:)
  enddo


end subroutine padding



subroutine box_smooth(s,h1,h2,h3,nz,nx,s_out)

implicit none

integer, intent(in)  :: h1,h2,h3,nz,nx
real,    intent(in)  :: s(nz,nx)
real,    intent(out) :: s_out(:,:)

integer :: nit,hh1,hh2,iter,i1,i2,i11,i12,i21,i22,iz,ix,izz,ixx
real, allocatable :: s0(:,:),s1(:,:),ss(:,:),s_temp(:,:)

allocate(s0(nz,nx),s_temp(nz,nx))

 do ix=1,nx
 do iz=1,nz
    s0(iz,ix)=s(iz,ix)
    s_out(iz,ix)=0.0
 enddo
 enddo

 nit=h3

 hh1=(h1-1)/2; hh2=(h2-1)/2

 allocate(s1(nz+2*hh1,nx+2*hh2))
 allocate(ss(h1,h2))

! print*,'3333333333'
 do iter=1,nit
   
    s_temp(:,:)=0.0 
    call padding(s0,nz,nx,hh1,hh2,s1)
!    print*,'444444'  
   
    do i2=hh2+1,hh2+nx 
    do i1=hh1+1,hh1+nz    
       i11=i1-hh1; i12=i1+hh1
       i21=i2-hh2; i22=i2+hh2
       ss=s1(i11:i12,i21:i22)   
       do ixx=1,h2
       do izz=1,h1
          s_temp(i1-hh1,i2-hh2)=s_temp(i1-hh1,i2-hh2)+ss(izz,ixx)
       enddo
       enddo
          s_temp(i1-hh1,i2-hh2)=s_temp(i1-hh1,i2-hh2)/(h1*h2)
    enddo
    enddo

    do ix=1,nx
    do iz=1,nz
       s0(iz,ix)=s_temp(iz,ix)
    enddo 
    enddo

 enddo

 do ix=1,nx
 do iz=1,nz
    s_out(iz,ix)=s0(iz,ix)
 enddo
 enddo

! print*,'4444444'
 deallocate(s0,s1,ss,s_temp)

end subroutine box_smooth


subroutine vel_smooth(vel,nz,nx,h1,h2,h3,vel2)

implicit none

real,    intent(in)  :: vel(:,:)
integer, intent(in)  :: nz,nx,h1,h2,h3
real,    intent(out) :: vel2(:,:)

integer :: ix,iz
real, allocatable :: slowness(:,:),ss2(:,:)

allocate(slowness(nz,nx),ss2(nz,nx))

 do ix=1,nx
 do iz=1,nz
    vel2(iz,ix)=0.0
    slowness(iz,ix)=1.0/vel(iz,ix)
 enddo
 enddo

! print*,'11111'
 call box_smooth(slowness,h1,h2,h3,nz,nx,ss2) 
! print*,'22222' 

 do ix=1,nx
 do iz=1,nz
    vel2(iz,ix)=1.0/ss2(iz,ix)
 enddo
 enddo

deallocate(slowness,ss2)

end subroutine vel_smooth


subroutine grdsmth2(ptr_grd, horsmt, versmt)
  implicit none
  real,intent(in out) :: ptr_grd(:,:)
 !*****************************************************
 ! GRDSMTH - Two dimensional elliptical Gaussian smoothing of the gradient
 !           using globals horsmt & versmt
 !
 !           It is the implementation of the 2D elliptical Gaussian function
 !           presented by wikipedia:
 !           http://en.wikipedia.org/wiki/Gaussian_function
 !
 ! Use size() instead of fullx and zlen in the different loops
 !
 ! Adrien ARNULF - Jan 2011
 !*****************************************************


  integer i,j
  real(kind=kind(1.e0)), allocatable :: temp(:,:),xx(:,:),yy(:,:),z(:,:)
  integer kk
  real(kind=kind(1.e0)) a,b,c,theta,sigma1,sigma2,scalesm
  integer::horsmt,versmt

  integer::nx,nz

  if (horsmt <= 0.or.versmt <= 0) then
     print *, " horsmt or versmt or both <= 0"
     return
  endif


      ! Create the two dimensional elliptical Gaussian smoothing function
      ! using globals horsmt & versmt to constrain its size.
   sigma1 = horsmt
   sigma2 = versmt
   theta = 0 ! Can become a variable if needed

   a = cos(theta)**2/2/sigma1**2 + sin(theta)**2/2/sigma2**2
   b = -sin(2*theta)/4/sigma1**2 + sin(2*theta)/4/sigma2**2
   c = sin(theta)**2/2/sigma1**2 + cos(theta)**2/2/sigma2**2

   nz = size(ptr_grd, dim=1)
   nx = size(ptr_grd, dim=2)
                                        ! Allocate local arrays
!  allocate(xx(2*horsmt+1,2*versmt+1))
!  allocate(yy(2*horsmt+1,2*versmt+1))
!  allocate(z(2*horsmt+1,2*versmt+1))
  allocate(xx(2*versmt+1,2*horsmt+1))
  allocate(yy(2*versmt+1,2*horsmt+1))
  allocate(z(2*versmt+1,2*horsmt+1))

   kk=0
   do i=1,(2*horsmt+1)
    !  xx(i,:)=-horsmt+kk
      xx(:,i)=-horsmt+kk
      kk=kk+1
   enddo

   kk=0
   do j=1,(2*versmt+1)
     ! yy(:,j)=-versmt+kk
      yy(j,:)=-versmt+kk
      kk=kk+1
   enddo

                 ! two dimensional elliptical Gaussian (z) defined within
                 ! the square(-horsmt:+horsmt,-versmt:+versmt)
                 ! size(z) = (2*horsmt+1 , 2*versmt+1)
   z = exp(-(a * xx**2 + 2*b* xx*yy + c* yy**2))

   scalesm = 1/sum(z)     ! scale factor so the filter response = 1
   z = scalesm*z
!   print *, " Size(z,1) =", size(z,1)
!   print *, " Size(z,2) =", size(z,2)


                    ! Allocate local array
   allocate(temp(nz,nx))

   print *, 'inside grdsmth2', 'nx', nx, 'nz', nz
   print *, 'inside grdsmth2', 'versmt', versmt, 'horsmt', horsmt

                ! Apply the Gaussian smoothing
   temp = ptr_grd

   print *, 'inside grdsmth2', size(temp, dim=2), size(temp, dim=1)

!   do i=versmt+1,nz-(versmt+1)
!      do j=horsmt+1,nx-(horsmt+1)
!        temp(j,i)=sum(ptr_grd(j-horsmt:j+horsmt,i-versmt:i+versmt)*z)
!      enddo
!   enddo

   do j=horsmt+1,nx-(horsmt+1)
      do i=versmt+1,nz-(versmt+1)
         temp(i,j) = sum(ptr_grd(i-versmt:i+versmt,j-horsmt:j+horsmt)*z)
      enddo 
   enddo 

   ptr_grd = temp

  deallocate(xx,yy,z,temp)

end subroutine grdsmth2



end module smooth_slow
