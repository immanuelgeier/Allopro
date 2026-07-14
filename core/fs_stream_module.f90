module fs_stream_module

use project_module
use dcd_reader
implicit none

private

public :: fs_stream_type
public :: fs_init_stream
public :: fs_load_next_file
public :: fs_stream_finished

!=========================================================
type :: fs_traj_block_type
    integer :: natoms = 0
    integer :: nframes = 0
    real, allocatable :: pos(:,:,:)
end type

!=========================================================
type :: fs_stream_type
    integer :: ifile = 1
    integer :: iframe = 1
    logical :: finished = .false.

    type(fs_traj_block_type) :: block
end type

contains

subroutine fs_init_stream(stream)
    type(fs_stream_type), intent(inout) :: stream
    stream%ifile = 1
    stream%iframe = 1
    stream%finished = .false.
end subroutine

subroutine fs_load_next_file(stream, project)
    type(fs_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    character(400) :: fname

    fname = project%state(1)%rep(1)%dcd_files(stream%ifile)

    call read_dcd_into_block(fname, stream%block%pos, &
                             stream%block%natoms, stream%block%nframes)

    stream%ifile = stream%ifile + 1
end subroutine

function fs_stream_finished(stream) result(f)
    type(fs_stream_type), intent(in) :: stream
    logical :: f
    f = stream%finished
end function

end module fs_stream_module