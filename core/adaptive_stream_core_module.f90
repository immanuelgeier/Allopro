module adaptive_stream_core_module

use project_module
use dcd_reader
implicit none

private
public :: traj_block_type
public :: adaptive_stream_type
public :: init_stream, preload, advance, get_block, stream_finished

!=========================================================
type :: traj_block_type
    integer :: natoms  = 0
    integer :: nframes = 0
    real, allocatable :: pos(:,:,:)   ! (3,natoms,blocksize)
end type

!=========================================================
type :: adaptive_stream_type
    type(traj_block_type), allocatable :: block(:)

    integer :: cur = 1
    integer :: nxt = 2

    integer :: istate = 0
    integer :: irep   = 0

    integer :: ifile  = 1
    integer :: iframe = 1

    integer :: blocksize = 32768

    integer :: natoms = 0
    logical :: initialized = .false.

    logical :: finished = .false.
end type
contains

subroutine init_stream(stream, istate, irep, blocksize)

    type(adaptive_stream_type), intent(inout) :: stream
    integer, intent(in) :: istate, irep, blocksize

    stream%istate = istate
    stream%irep   = irep
    stream%blocksize = blocksize

    stream%ifile  = 1
    stream%iframe = 1

    stream%cur = 1
    stream%nxt = 2

    stream%finished = .false.

    if (allocated(stream%block)) deallocate(stream%block)
    allocate(stream%block(2))

end subroutine

subroutine swap(stream)

    type(adaptive_stream_type), intent(inout) :: stream
    integer :: tmp

    tmp = stream%cur
    stream%cur = stream%nxt
    stream%nxt = tmp

end subroutine


function get_block(stream) result(b)

    type(adaptive_stream_type), intent(in) :: stream
    type(traj_block_type) :: b

    b = stream%block(stream%cur)

end function


subroutine preload(stream, project)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    call fill_block(stream, project, stream%block(stream%cur))

    if (.not. stream%finished) then
        call fill_block(stream, project, stream%block(stream%nxt))
    end if

end subroutine


subroutine advance(stream, project)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    if (stream%finished) return

    call swap(stream)
    call fill_block(stream, project, stream%block(stream%nxt))

end subroutine


subroutine fill_block(stream, project, block)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project
    type(traj_block_type), intent(inout) :: block

    character(400) :: fname

    if (stream%ifile > &
        project%state(stream%istate)%rep(stream%irep)%ndcd) then

        stream%finished = .true.
        block%nframes = 0
        return

    end if

    fname = project%state(stream%istate)%rep(stream%irep)% &
            dcd_files(stream%ifile)

    if (allocated(block%pos)) then
        deallocate(block%pos)
    end if

    call read_dcd_into_block( &
        fname, &
        block%pos, &
        block%natoms, &
        block%nframes )

    stream%ifile = stream%ifile + 1

end subroutine

!=========================================================
function stream_finished(stream) result(f)
!=========================================================
    type(adaptive_stream_type), intent(in) :: stream
    logical :: f

    f = stream%finished

end function


end module adaptive_stream_core_module