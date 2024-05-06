! Main program 
! Peng Guo
! CSIRO

program main

!use apps2d_staggered_fwi
!use apps2d_wd_staggered
!use apps2d_fwi_staggered
use fwi_elastic_main
!use apps2d_modeling_tmp
use modeling_elastic_main
use parser

implicit none

character*256 :: parfile

call readCmd('par',parfile,'parfile')
call readParFile(parfile,'jobtype',jobtype)


!if (jobtype(1:18)=='q_forward_modeling') then
!    call q_forward_modeling(parfile)


if(jobtype(1:16)=='elastic_modeling')then
  call modeling_elastic(parfile)

!elseif(jobtype(1:11)=='wave_disper')then
!  call wd_inversion(parfile)

elseif(jobtype(1:11)=='elastic_fwi')then
  call fwi_elastic(parfile)

endif


end program main

