! module defining global variables
!!!!!!!!!!!!!!

module global

implicit none

! Constants
real(4), parameter :: pi        = 3.14159265358979
real(4), parameter :: twopi     = 2.0*pi

! Type of job
character(len=200) :: jobtype, input, output

! Record length unit
integer, parameter :: i4 = 4

! 2-8 FD coefficients
real, parameter :: c1_2nd_order = -205.0/72.0
real, parameter :: c2_2nd_order = 8.0/5.0
real, parameter :: c3_2nd_order = -1.0/5.0
real, parameter :: c4_2nd_order = 8.0/315.0
real, parameter :: c5_2nd_order = -1.0/560.0
real, parameter :: c1_staggered = 9.0/8.0
real, parameter :: c2_staggered = -1.0/24.0
real, parameter :: c1_8th = -205.0/72.0
real, parameter :: c2_8th = 8.0/5.0
real, parameter :: c3_8th = -1.0/5.0
real, parameter :: c4_8th = 8.0/315.0
real, parameter :: c5_8th = -1.0/560.0

real, parameter :: c1_stag_6th = 1.171875000000000
real, parameter :: c2_stag_6th = -6.5104166666666687E-002
real, parameter :: c3_stag_6th = 4.6875000000000009E-003

real, parameter :: c1_staggered_8th = 1225.0/1024.0
real, parameter :: c2_staggered_8th = -245.0/3072.0
real, parameter :: c3_staggered_8th = 49.0/5120.0
real, parameter :: c4_staggered_8th = -5.0/7168.0

real, parameter :: c1_elastic_2th = 1.0

! 0.5%
real, parameter :: c1_elastic_4th = 1.1534
real, parameter :: c2_elastic_4th = -0.052806

real, parameter :: c1_elastic_8th = 1.19629
real, parameter :: c2_elastic_8th = -0.0797526
real, parameter :: c3_elastic_8th = 0.00957031
real, parameter :: c4_elastic_8th = -0.000697545

! Boundary record length

integer, parameter:: bc_len=5,bc_len_1=4


end module global
