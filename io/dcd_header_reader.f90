module dcd_header_reader

use iso_fortran_env, only: real32
implicit none

type :: dcd_header_info_type
    integer :: natom_found_in_dcd       = 0
    integer :: nframes_dcd_found_in_dcd = 0
end type

contains

!=========================================================
subroutine read_dcd_header_info(filename, info)

    character(*), intent(in) :: filename
    type(dcd_header_info_type), intent(out) :: info

    character(4) :: dummyc
    integer :: dummyi, i, io
    real :: dummyr

    open(10,file=trim(filename),status='old',form='unformatted',iostat=io)
    if (io /= 0) stop "ERROR opening DCD"

    read(10) dummyc, info%nframes_dcd_found_in_dcd, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(10) dummyi, dummyr
    read(10) info%natom_found_in_dcd

    close(10)

end subroutine

end module dcd_header_reader