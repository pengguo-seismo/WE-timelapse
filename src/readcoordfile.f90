
        

subroutine readcoordfile(coordfile, coord)

use datatype
use mmi_mpi

character(len=*), intent(in)   :: coordfile
type(acquisition), intent(out) :: coord
integer                        :: i, j, is0, is, ig, ns, ng, ngmax
real                           :: xs, zs, xg, zg, t

if (rank == 0) then
  open(10,file=coordfile,form='formatted')

  ! Determine the number of sources
  ns = 0
  read(10,*,end=100) is0, ig, xs, zs, xg, zg, t
  ns = 1
  do i=2,1000000000
    read(10,*,end=100) is, ig, xs, zs, xg, zg, t
    if (is.ne.is0) then
      is0 = is
      ns = ns + 1
    endif
  enddo
  100 continue
  rewind(10)
  coord%ns = ns
  allocate(coord%ng(ns))
  allocate(coord%xs(ns))
  allocate(coord%zs(ns))

  ! Determine the number of geophones for each shot
  read(10,*,end=200) is0, ig, xs, zs, xg, zg, t
  ns = 1
  ng = 1
  do i=2,1000000000
    read(10,*,end=200) is, ig, xs, zs, xg, zg, t
    if (is.ne.is0) then
      is0 = is
      coord%ng(ns) = ng
!    write(*,*) 'ng(',ns,') = ', coord%ng(ns)
      ns = ns + 1
      ng = 1
    else
      ng = ng + 1
    endif
  enddo
  200 continue
  coord%ng(ns) = ng
  !write(*,*) 'ng(',ns,') = ', coord%ng(ns)
  rewind(10)

    ! Determine the maximum number of geophones per shot
  coord%ngmax = 0
  do is=1,coord%ns
    if (coord%ngmax < coord%ng(is)) coord%ngmax = coord%ng(is)
  enddo
  write(*,*) 'ng max = ', coord%ngmax

  allocate(coord%xg(is,coord%ngmax))
  allocate(coord%zg(is,coord%ngmax))
  allocate(coord%t(is,coord%ngmax))

  ! Read source and receiver positions
  do i=1,coord%ns
    do j=1,coord%ng(i)
      read(10,*,end=300) is, ig, xs, zs, xg, zg, t
      coord%xs(i) = xs
      coord%zs(i) = zs
      coord%xg(i,j) = xg
      coord%zg(i,j) = zg
      coord%t(i,j) = t
    enddo
  enddo
  300 continue
  close(10)

!  open(10,file='coord.txt',form='formatted')
!  do is=1,coord%ns
!    do ig=1,coord%ng(is)
!      write(10,*) is,ig,coord%xs(is),coord%zs(is),coord%xg(is,ig),coord%zg(is,ig),coord%t(is,ig)
!    enddo
!  enddo
!  close(10)
  write(*,*) 'ns = ', ns
endif

end subroutine readcoordfile



