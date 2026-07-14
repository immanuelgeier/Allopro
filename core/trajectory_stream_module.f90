module trajectory_stream_module

use project_module
use dcd_reader
implicit none

type :: traj_block_type
    integer :: nframes = 0
    integer :: natoms  = 0

    real, allocatable :: pos(:,:,:)
    character(400) :: file = ''
end type


type :: traj_stream_type
    type(traj_block_type) :: A, B

    integer :: istate = 0
    integer :: irep   = 0
    integer :: ifile  = 1

    logical :: finished = .false.
end type

contains

subroutine init_stream(stream, istate, irep)
    type(traj_stream_type), intent(inout) :: stream
    integer, intent(in) :: istate, irep

    stream%istate  = istate
    stream%irep    = irep
    stream%ifile   = 1
    stream%finished = .false.

    stream%A%nframes = 0
    stream%B%nframes = 0
end subroutine

subroutine advance(stream, project)
    type(traj_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    character(400) :: fname
    integer :: ndcd

    call swap(stream)

    ndcd = project%state(stream%istate)%rep(stream%irep)%ndcd

    if (stream%ifile > ndcd) then
        stream%finished = .true.
        stream%B%nframes = 0
        return
    end if

    fname = project%state(stream%istate)%rep(stream%irep)%dcd_files(stream%ifile)
    stream%ifile = stream%ifile + 1

    call load_into(stream%B, fname)
    stream%finished = .false.

end subroutine


subroutine load_into(block, filename)
    type(traj_block_type), intent(inout) :: block
    character(*), intent(in) :: filename

    call read_dcd_into_block(filename, block%pos, block%natoms, block%nframes)
    block%file = filename
end subroutine


!=========================================================
subroutine swap(stream)
!=========================================================
    type(traj_stream_type), intent(inout) :: stream
    type(traj_block_type) :: tmp

    tmp = stream%A
    stream%A = stream%B
    stream%B = tmp

end subroutine



!=========================================================
subroutine preload(stream, project)
!=========================================================
    type(traj_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    stream%ifile = 1
    stream%finished = .false.

    call advance(stream, project)
    call advance(stream, project)

    call swap(stream)

end subroutine


end module trajectory_stream_module