

subroutine xcorr_calc(dat1,dat2,npts,i1,i2,nlen,ishift,cc_max)

    implicit none

    !! dat1: obs
    !! dat2: syn
    real(kind=kind(1.e0)), intent(in) :: dat1(npts), dat2(npts)
    integer, intent(in) :: npts, i1, i2

    ! outputs:
    ! ishift = index lag (d-s) for max cross correlation
    ! cc_max = maximum of cross correlation (normalised by sqrt(synthetic*data))
    integer, intent(out) :: ishift
    real(kind=kind(1.e0)), intent(out) :: cc_max
    real(kind=kind(1.e0))::dlnA

    ! local variables
    integer :: nlen, i_time
    integer :: i_left, i_right, i, j, id_left, id_right
    real(kind=kind(1.e0)) :: cc, norm, norm_s
    real(kind=kind(1.e0)) :: SMALL_VAL

    SMALL_VAL = 1.0e-32

    ! initialise shift and cross correlation to zero
    ishift = 0
    dlnA = 0.0
    cc_max = 0.0

    ! print*,'ishift,dlnA,cc_max :',ishift,dlnA,cc_max

    if (i2 == 0 .or. i1 == 0) then
       ishift = 0
       cc_max = 0.0
       return
    endif

    if (i1.lt.0 .or. i1.gt.i2 .or. i2.gt.npts) then
        write(*,*) 'Error with window limits: i1, i2, npts ', i1, i2, npts
        return
    endif

    ! length of window (number of points, including ends)
!    nlen = i2 - i1 + 1

    ! power of synthetic signal in window
    norm_s = sqrt(sum(dat2(i1:i2)*dat2(i1:i2)))

    ! left and right limits of index (time) shift search
    ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and CC.
    !       How far to look outside, in theory, should be another parameter.
    !       However, it does not matter as much if the data and synthetics are
    !          zeroed outside the windows, as currently done in calc_criteria.
    i_left = -1*int(nlen/2.0)
    i_right = int(nlen/2.0)
    !!  i_left = -nlen
    !!  i_right = nlen

    ! i is the index to shift to be applied to DATA (d)
    do i = i_left, i_right

       ! normalization factor varies as you take different windows of d
       id_left = max(1,i1+i)      ! left index for data window
       id_right = min(npts,i2+i)  ! right index for data window
       norm_s = norm_s * sqrt(sum(dat1(id_left:id_right)*(dat1(id_left:id_right))))

       ! cc as a function of i
       cc = 0.0
       do j = i1, i2   ! loop over full window length
         if((j+i).ge.1 .and. (j+i).le.npts) cc = cc + dat2(j)*dat1(j+i)  ! d is shifted by i
       enddo
       !print*,'norm,norm_s :',norm,norm_s

       norm = maxval(dat1(i1:i2))

       ! normalized by norm of data
       if(norm > SMALL_VAL ) cc = cc/norm

       ! keeping cc-max only
       if (cc .gt. cc_max) then
          cc_max = cc
          ishift = i
          !print*,'i,cc,cc_max :',i,cc,cc_max
       endif
    enddo

end subroutine xcorr_calc



subroutine xcorr_calc_mod(dat1,dat2,npts,i1,i2,nlen,ishift,cc_max)

    implicit none

    !! dat1: syn
    !! dat2: obs
    real(kind=kind(1.e0)), intent(in) :: dat1(npts), dat2(npts)
    integer, intent(in) :: npts, i1, i2

    ! outputs:
    ! ishift = index lag (d-s) for max cross correlation
    ! cc_max = maximum of cross correlation (normalised by sqrt(synthetic*data))
    integer, intent(out) :: ishift
    real(kind=kind(1.e0)), intent(out) :: cc_max
    real(kind=kind(1.e0))::dlnA

    ! local variables
    integer :: nlen, i_time
    integer :: i_left, i_right, i, j, id_left, id_right
    real(kind=kind(1.e0)) :: cc, norm, norm_s
    real(kind=kind(1.e0)) :: SMALL_VAL

    SMALL_VAL = 1.0e-32

    ! initialise shift and cross correlation to zero
    ishift = 0
    dlnA = 0.0
    cc_max = 0.0

    ! print*,'ishift,dlnA,cc_max :',ishift,dlnA,cc_max

    if (i2 == 0 .or. i1 == 0) then
       ishift = 0
       cc_max = 0.0
       return
    endif

    if (i1.lt.0 .or. i1.gt.i2 .or. i2.gt.npts) then
        write(*,*) 'Error with window limits: i1, i2, npts ', i1, i2, npts
        return
    endif

    ! length of window (number of points, including ends)
!    nlen = i2 - i1 + 1

    ! power of synthetic signal in window
    norm_s = sqrt(sum(dat2(i1:i2)*dat2(i1:i2)))

    ! left and right limits of index (time) shift search
    ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and CC.
    !       How far to look outside, in theory, should be another parameter.
    !       However, it does not matter as much if the data and synthetics are
    !          zeroed outside the windows, as currently done in calc_criteria.
    i_left = -1*int(nlen/2.0)
    i_right = int(nlen/2.0)
    !!  i_left = -nlen
    !!  i_right = nlen

    ! i is the index to shift to be applied to DATA (d)
    do i = i_left, i_right

       ! normalization factor varies as you take different windows of d
       id_left = max(1,i1+i)      ! left index for data window
       id_right = min(npts,i2+i)  ! right index for data window
       norm = norm_s * sqrt(sum(dat1(id_left:id_right)*(dat1(id_left:id_right))))

       ! cc as a function of i
       cc = 0.0
       do j = i1, i2   ! loop over full window length
         if((j+i).ge.1 .and. (j+i).le.npts) cc = cc + dat2(j)*dat1(j+i)  ! d is shifted by i
       enddo
       !print*,'norm,norm_s :',norm,norm_s

       norm = maxval(dat2(i1:i2))

       ! normalized by norm of data
       if(norm > SMALL_VAL ) cc = cc/norm

       ! keeping cc-max only
       if (cc .gt. cc_max) then
          cc_max = cc
          ishift = i
          !print*,'i,cc,cc_max :',i,cc,cc_max
       endif
    enddo
 
end subroutine xcorr_calc_mod

