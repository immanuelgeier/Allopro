module bs_stream_module

use project_module
use dcd_reader
implicit none

private

public :: bs_stream_type
public :: bs_init_stream
public :: bs_preload
public :: bs_advance
public :: bs_get_block
public :: bs_stream_finished

!=========================================================
type :: bs_traj_block_type
    integer :: natoms = 0
    integer :: nframes = 0
    real, allocatable :: pos(:,:,:)
end type

!=========================================================
type :: bs_stream_type
    type(bs_traj_block_type), allocatable :: block(:)

    integer :: cur = 1
    integer :: nxt = 2

    integer :: istate = 0
    integer :: irep = 0

    integer :: ifile = 1
    integer :: iframe = 1

    integer :: blocksize = 32768
    logical :: finished = .false.
end type

contains

subroutine bs_init_stream(stream, istate, irep, blocksize)
    type(bs_stream_type), intent(inout) :: stream
    integer, intent(in) :: istate, irep, blocksize

    stream%istate = istate
    stream%irep = irep
    stream%blocksize = blocksize

    stream%ifile = 1
    stream%iframe = 1
    stream%cur = 1
    stream%nxt = 2
    stream%finished = .false.

    if (allocated(stream%block)) deallocate(stream%block)
    allocate(stream%block(2))
end subroutine

subroutine bs_swap(stream)
    type(bs_stream_type), intent(inout) :: stream
    integer :: tmp
    tmp = stream%cur
    stream%cur = stream%nxt
    stream%nxt = tmp
end subroutine

function bs_get_block(stream) result(b)
    type(bs_stream_type), intent(in) :: stream
    type(bs_traj_block_type) :: b
    b = stream%block(stream%cur)
end function

subroutine bs_preload(stream, project)
    type(bs_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    call bs_fill(stream, project, stream%block(stream%cur))
    if (.not. stream%finished) then
        call bs_fill(stream, project, stream%block(stream%nxt))
    end if
end subroutine

subroutine bs_advance(stream, project)
    type(bs_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    if (stream%finished) return

    call bs_swap(stream)
    call bs_fill(stream, project, stream%block(stream%nxt))
end subroutine

subroutine bs_fill(stream, project, block)
    type(bs_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project
    type(bs_traj_block_type), intent(inout) :: block

    character(400) :: fname
    integer :: natom, nframes

    block%nframes = 0

    if (stream%ifile > project%state(stream%istate)%rep(stream%irep)%ndcd) then
        stream%finished = .true.
        return
    end if

    fname = project%state(stream%istate)%rep(stream%irep)%dcd_files(stream%ifile)

    call read_dcd_into_block(fname, block%pos, natom, nframes)

    stream%ifile = stream%ifile + 1
end subroutine

function bs_stream_finished(stream) result(f)
    type(bs_stream_type), intent(in) :: stream
    logical :: f
    f = stream%finished
end function

end module bs_stream_module