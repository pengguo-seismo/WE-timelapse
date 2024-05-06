module win_grad_mod
  use geometry_mod
  use dynamic_win_mod
  implicit none

  real(kind=kind(1.d0)),  allocatable :: grdla(:,:)   !! lambda gradient; see eq 3.13 Shipp thesis
  real(kind=kind(1.d0)),  allocatable :: grdmu(:,:)   !! mu gradient - (eq3.10 shipp thesis? more terms.......)

  real(kind=kind(1.d0)),  allocatable :: grdmu_a(:,:) !! sub-terms of mu-gradient calculated in s/r *crosco*
  real(kind=kind(1.d0)),  allocatable :: grdmu_b(:,:)
  real(kind=kind(1.d0)),  allocatable :: grdmu_c(:,:)

  contains

    subroutine win_grad_alloc(xlen,zlen)
      integer, intent(in) :: xlen, zlen

      allocate(grdla(xlen,zlen),grdmu(xlen,zlen))
      allocate(grdmu_a(xlen,zlen),grdmu_b(xlen,zlen),grdmu_c(xlen,zlen))

      if ( .not. (allocated(grdla) .and. allocated(grdmu)) ) then
        print '(A)','WIN_GRAD_ALLOC grdla or grdmu not allocated'
        stop 1
      endif

      if ( .not. (allocated(grdmu_a) .and. allocated(grdmu_b) .and.  &
                                                   allocated(grdmu_c)) ) then
        print '(A)','WIN_GRAD_ALLOC grdmu_[abc] not allocated'
        stop 1
      endif

    end subroutine win_grad_alloc

    subroutine win_grad_dealloc
      deallocate(grdla,grdmu)
      deallocate(grdmu_a,grdmu_b,grdmu_c)
    end subroutine win_grad_dealloc

    subroutine win_grad_zero
      grdmu = 0.0
      grdla = 0.0

      grdmu_a = 0.0
      grdmu_b = 0.0
      grdmu_c = 0.0
    end subroutine win_grad_zero

    subroutine grconv_crosco(sn, grdtop)
      integer, intent(in) :: sn
      integer, intent(in) :: grdtop
       ! Finish gradient conditions after bacmod calculation

      integer :: i,j
      real(kind=kind(1.d0)),  allocatable :: lambuf(:,:),mubuf(:,:)
      real(kind=kind(1.d0))  :: terma,termb,termc,denom
      integer :: xlen, xlen_eff, zlen
      integer :: winsta, start_eff, padval

      zlen = size(grdla,2)
      call current_wind_coords(sn,winsta,xlen, xlen_eff)

      print *, 'In grconv_crosco sn',sn,'winsta',winsta,'xlen',xlen,''

       do j=1,zlen
          do i=1,xlen
             grdla(i,j) = - grdla(i,j) / (2.0*(lambda(i,j)+mu(i,j)))
             grdla(i,j) = grdla(i,j) / (2.0*(lambda(i,j)+mu(i,j)))

             if (mu(i,j) > 0) then
                denom=(lambda(i,j)+mu(i,j))**2.0
                denom=8.0*(mu(i,j)**2.0)*denom

                terma=2.0*(lambda(i,j)*lam2mu(i,j))/denom
                terma=terma*grdmu_a(i,j)

                termb=((lambda(i,j)**2.0)+(lam2mu(i,j)**2))/denom
                termb=termb*grdmu_b(i,j)

                termc=1.0/(mu(i,j)**2.0)
                termc=termc*grdmu_c(i,j)
                grdmu(i,j) = terma - termb - termc
             else
                grdmu(i,j) = 0.0
             endif
          enddo
       enddo


                    ! Revert to effective xlen window values

      start_eff = start_pos_win(sn)
      padval    = window_padding()
      allocate(lambuf(xlen_eff,zlen),mubuf(xlen_eff,zlen))

       do i=1,xlen_eff
         lambuf(i,:) = grdla(start_eff-1+i,:)
          mubuf(i,:) = grdmu(start_eff-1+i,:)
       enddo

       grdla = 0.0
       grdmu = 0.0

       grdla(1:xlen_eff,:) = lambuf
       grdmu(1:xlen_eff,:) =  mubuf

       deallocate(lambuf,mubuf)

              ! Extend left & right by replication to write over absorbing BC columns
              ! Also extend over the bottom of the grid
              ! Top is special as the gradient will be zeroed here anyway using
              ! input parameter grdtop
        do j=grdtop,botlim-2
          do i=1,leflim+2
             grdla(i,j)=grdla(leflim+2,j)
             grdmu(i,j)=grdmu(leflim+2,j)
          enddo

          do i=riglim-2-padval,xlen_eff
             grdla(i,j)=grdla(riglim-2-padval,j)
             grdmu(i,j)=grdmu(riglim-2-padval,j)
          enddo
       enddo

       do i=1,xlen_eff
           do j=botlim-2,zlen
              grdla(i,j)=grdla(i,botlim-2)
              grdmu(i,j)=grdmu(i,botlim-2)
           enddo
      enddo

    end subroutine grconv_crosco


    subroutine grconv(sn,dt,sz,grdvp,grdvs)
      use allglobals, only:outpath,outlen
      integer, intent(in) :: sn
      real(kind=kind(1.d0)),   intent(in)  :: dt
      integer,  intent(in) :: sz
      real(kind=kind(1.d0)),  intent(out) :: grdvp(:,:), grdvs(:,:)
      ! GRCONV - gradient conversion for a single shot. Converts gradient in
      ! lambda/mu to gradients in velocity
      ! The lambda/mu gradients, grdla/grdmu, have already been left adjusted in
      !  window by subroutine bacmod but background velocities/densities winvp/s,
      ! rho have not.

      integer i,j
      real(kind=kind(1.d0)),  allocatable :: taper(:),loches(:,:)
      real(kind=kind(1.d0)) terma,termb
      integer :: start_eff
      integer :: xlen_eff, zlen
      integer :: winsta, xlen
      character*256::a1,a2
      character*256::file_bin

      zlen = size(grdvp,2)
      call current_wind_coords(sn,winsta,xlen, xlen_eff)

      print *,'grconv ', 'sn',sn, 'winsta',winsta,'xlen',xlen,'xlen_eff',xlen_eff
      print *,'zlen', zlen

      start_eff = start_pos_win(sn)

      allocate(taper(xlen_eff),loches(xlen_eff,zlen))
                          ! Initialize side taper
      taper = 1.0

      do i=1,5
         taper(i)=real(i-1)/real(5-1)
      enddo

      do i=xlen_eff-4, xlen_eff
         taper(i) = real(xlen_eff-i)/real(4)
      enddo

      grdvp = 0.0
      grdvs = 0.0

      do j=1,zlen
         loches(:,j)=sqrt(real(abs(j-sz)))
      enddo

!! this is only for debug
!      write(a1,"(I5.5)") sn
!      file_bin='gra_vp'//trim(a1)
!      open(12,file=file_bin,access='direct',form='unformatted',recl=4*nz)
!       do i=1,nx
!          write(12,rec=i) (vp(j,i),j=1,nz)
!       enddo
!      close(12)

                          ! Back in effective xlen window

      do i=1,xlen_eff
         do j=1,zlen
            rho(i-1+start_eff,j) = dt/rho(i-1+start_eff,j)
            grdvp(i,j) = 2.0*rho(i-1+start_eff,j)*winvp(i-1+start_eff,j)* &
                       grdla(i,j)
            grdvp(i,j) = grdvp(i,j)*loches(i,j)

            ! Change polarity of gradient to allow all notation to agree with Heiner's

            grdvp(i,j)=-grdvp(i,j)
            grdvp(i,j)=grdvp(i,j)*taper(i)

            terma=-4.0*rho(i-1+start_eff,j)*winvs(i-1+start_eff,j)*grdla(i,j)
            termb=2.0*rho(i-1+start_eff,j)*winvs(i-1+start_eff,j)*grdmu(i,j)
            grdvs(i,j) = terma+termb
            grdvs(i,j) = grdvs(i,j)*loches(i,j)
            grdvs(i,j) = -grdvs(i,j)
            grdvs(i,j) = grdvs(i,j)*taper(i)

            rho(i-1+start_eff,j) = dt/rho(i-1+start_eff,j)

         enddo
      enddo

      if(.false.) then

       write(a1,"(I5.5)") sn
       file_bin=outpath(1:outlen)//'/gra_vpin'//trim(a1)
       open(12,file=file_bin,access='direct',form='unformatted',recl=4*(zlen))
        do i=1,xlen_eff
           write(12,rec=i) (real(grdvp(i,j)),j=1,zlen)
        enddo
       close(12)
 
       endif

!       file_bin=outpath(1:outlen)//'/rhoin'//trim(a1)
!       open(12,file=file_bin,access='direct',form='unformatted',recl=4*(zlen-4+1))
!        do i=1,xlen_eff
!           write(12,rec=i) (real(rho(i,j)),j=4,zlen)
!        enddo
!       close(12)

!       file_bin=outpath(1:outlen)//'/loches'//trim(a1)
!       open(12,file=file_bin,access='direct',form='unformatted',recl=4*(zlen-4+1))
!        do i=1,xlen_eff
!           write(12,rec=i) (real(loches(i,j)),j=4,zlen)
!        enddo
!       close(12)

      deallocate(taper,loches)

    end subroutine grconv

end module win_grad_mod
