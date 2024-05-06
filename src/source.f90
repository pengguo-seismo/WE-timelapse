! module for generating source functions
!
module source

use global, only : PI
implicit none

!real, parameter :: timeshift = 1.0
real, parameter :: timeshift = 1.41421356237309504880

contains

!-------------------------------------------------------------------------
subroutine sourcegen(parfile)

use parser
use io

character(len=*), intent(in) :: parfile

character(len=100) :: sourcefile
real, allocatable  :: s(:)
integer            :: nt
real               :: dt, f

call readParFile(parfile, 'NT', nt)
call readParFile(parfile, 'DT', dt)
call readParFile(parfile, 'F', f)
call readParFile(parfile, 'SOURCEFILE', sourcefile)

allocate(s(nt))
call ricker_long(s, nt, f, dt)
call write_sufile(sourcefile, s, nt, dt)
deallocate(s)

end subroutine sourcegen

!-------------------------------------------------------------------------
subroutine getsource_is(par, is, s)

use datatype
use global
use string

type(param), intent(in) :: par
integer, intent(in)     :: is
real, intent(out)       :: s(:)
integer                 :: i
real                    :: src(par%nt)
character(len=200)      :: str

call filename(str,'../source/src',is,'.dat')
open(10,file=str,access='direct',recl=i4,err=100)
do i=1,par%nt
  read(10,rec=i,err=50) src(i)
enddo
50 continue
close(10)
s(1:par%nt-300) = src(301:par%nt)
100 continue

end subroutine getsource_is

!-------------------------------------------------------------------------
subroutine getsource(par, s)

use datatype
use global
use mmi_mpi
use io

type(param), intent(in) :: par
real, intent(out)       :: s(:)
integer                 :: i, source_is_input, shift

integer :: ii

shift = 500

if (rank == 0) then
  if (par%sourcefile(1:3) /= 'n/a') then
    s = 0.0
!    open(10,file=par%sourcefile,access='direct',recl=i4,err=100)
    open(10,file=par%sourcefile,err=100)
    do i=1,par%nt
      read(10,*,err=50) ii, s(i)
    enddo
    open(122, file='wavelet.save', access='direct', form='unformatted', recl=i4*par%nt)
    write(122, rec=1) (s(i), i=1,par%nt)
    close(122)
    50 continue
    write(*,*) 'input source has ', i,' time samples'
    if (i == 1) goto 100
    close(10)
    source_is_input = 1
    goto 200
  endif
  100 continue
  source_is_input = 0
  200 continue
  if (source_is_input.eq.0) then
    write(*,*) 'No source input'
!    call ricker_long(s,par%nt,par%f,par%dt)
     call rickerfunc_new_seri(par%dt,par%f,par%nt,s) 

    if (par%sourcetype(1:4) == 'diff') then
      write(*,*) 'diff source'
      ! Apply a time-derivative to the source to be equivalent to staggered-grid code
      do i=1,par%nt-1
        s(i) = (s(i+1)-s(i))/par%dt
      enddo
      s(par%nt) = s(par%nt-1)
    endif

    !>>>>>>>>>>>>>>>>>>> For Brazion data
!    s(shift+1:par%nt) = s(1:par%nt-shift)
!    s(1:shift) = 0.0
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  endif
endif
call MPI_BCAST(s,par%nt,MPI_REAL,0,MPI_COMM_WORLD,ierr)

end subroutine getsource

!-------------------------------------------------------------------------
subroutine source_normal(par, s)

use datatype
use global

type(param), intent(in) :: par
real, intent(out)       :: s(:)
integer                 :: i, source_is_input

open(10,file=par%sourcefile,access='direct',recl=i4,err=100)
do i=1,1000000
  read(10,rec=i,err=50) s(i)
enddo
50 continue
close(10)
source_is_input = 1
goto 200
100 continue
source_is_input = 0
200 continue
if (source_is_input.eq.0) then
  call ricker_long(s,par%nt,par%f,par%dt)
endif

end subroutine source_normal

!-------------------------------------------------------------------------
subroutine source_diff(par, s)

use datatype
use global

type(param), intent(in) :: par
real, intent(out)       :: s(:)
integer                 :: i, source_is_input

open(10,file=par%sourcefile,access='direct',recl=i4,err=100)
do i=1,1000000
  read(10,rec=i,err=50) s(i)
enddo
50 continue
close(10)
source_is_input = 1
goto 200
100 continue
source_is_input = 0
200 continue
if (source_is_input.eq.0) then
  call ricker_long(s,par%nt,par%f,par%dt)

  ! Apply a time-derivative to the source to be equivalent to staggered-grid code
  do i=1,par%nt-1
    s(i) = (s(i+1)-s(i))/par%dt
  enddo
  s(par%nt) = s(par%nt-1)
endif

end subroutine source_diff

!-------------------------------------------------------------------------
function wavelet_length(f, dt)

integer :: wavelet_length
real, intent(in) :: f, dt

wavelet_length = floor(timeshift/(f*dt)+0.5)*2+1

end function wavelet_length

!-------------------------------------------------------------------------
subroutine ricker(s, np, f, dt)
! Creaing a Ricker source wavelet 

integer, intent(in)  :: np
real,    intent(in)  :: f, dt
real,    intent(out) :: s(:)

integer :: it
real    :: pi2, const, b, tshift, tim1, tim2, u, amp, smax

tshift = timeshift/f
pi2 = sqrt(PI)/2.0
b = sqrt(6.0)/(PI*f)
const = 2.0*sqrt(6.0)/b

smax = 0.0
do it=1,np
  tim1 = real(it-1)*dt
  tim2 = tim1 - tshift
  u = const*tim2
  amp = ((u*u)/4.0-0.5)*pi2*(exp(-u*u/4.0))
  s(it) = -amp
  if (smax .lt. abs(amp)) then
    smax = abs(amp)
  endif
enddo

smax = smax*2
do it=1,np
  s(it) = s(it)/smax
enddo

end subroutine ricker


!-------------------------------------------------------------------------
subroutine ricker_long(s, nt, f, dt)
! Creating a Ricker source wavelet with the same length as the traces

integer, intent(in)  :: nt
real,    intent(in)  :: f, dt
real,    intent(out) :: s(:)
real,allocatable     :: s1(:)

integer :: it, np
real    :: pi2, const, b, tshift, tim1, tim2, u, amp, smax

np = wavelet_length(f, dt)
allocate(s1(np))
tshift = timeshift/f
!tshift = timeshift/f + 29.0*dt ! For BP data, fpeak = 3.5 Hz
pi2 = sqrt(PI)/2.0
b = sqrt(6.0)/(PI*f)
const = 2.0*sqrt(6.0)/b

smax = 0.0
do it=1,np
  tim1 = real(it-1)*dt
  tim2 = tim1 - tshift
  u = const*tim2
  amp = ((u*u)/4.0-0.5)*pi2*(exp(-u*u/4.0))
  s1(it) = -amp
  if (smax .lt. abs(amp)) then
    smax = abs(amp)
  endif
enddo

smax = smax*2
do it=1,np
  s1(it) = s1(it)/smax
enddo

do it=1,nt
   if(it.le.np) then
      s(it)=s1(it)
   else
      s(it)=0.0
   endif
enddo

end subroutine ricker_long


!! I will use this wavelet
      subroutine rickerfunc_new_seri(dt,fp,nt,ricker)

       integer length,nt,nnt,nt_fft
       integer,parameter::dp=kind(0.e0)
       real(dp)::dt,fp,tshift
       real::PI
       real(dp)::ricker(1:nt)

       integer::i,it
       real::a,t,t0,amax

       PI=4.0*atan(1.0)

       tshift=0.0

       ricker=0.

       nnt=nt

       amax=0.

       a=PI*PI*fp*fp
       t0 = 1.50 /fp

       do it=1,nt
          t = real(it-1)*dt
          ricker(it) = (1.e0 - 2.e0*a*(t-t0)**2)*exp(-a*(t-t0)**2)
          if(abs(ricker(it)).gt.amax) amax=abs(ricker(it))
       enddo

       do it=1,nt
          ricker(it)=ricker(it)/amax
       enddo

        open(12,file='ricker_wavelet',status='replace',access='direct',&
                form='unformatted',recl=4*nt)
         write(12,rec=1) (ricker(i),i=1,nt)
        close(12)

        open(12,file='ricker_wavelet.ascii')
         do i=1,nt
            write(12,*) i,ricker(i)
         enddo
        close(12)

      end


end module source

