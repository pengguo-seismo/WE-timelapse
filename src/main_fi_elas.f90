! Peng 
! peng.guo@csiro.au

program main

use fwi_elastic_main
use modeling_elastic_main
use parser

implicit none

character*256 :: parfile

call readCmd('par',parfile,'parfile')
call readParFile(parfile,'jobtype',jobtype)



if(jobtype(1:16)=='elastic_modeling')then
  call modeling_elastic(parfile)


elseif(jobtype(1:11)=='elastic_fwi')then
  call fwi_elastic(parfile)

endif


end program main

