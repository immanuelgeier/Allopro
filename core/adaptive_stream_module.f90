module adaptive_stream_module

use project_module
use adaptive_dcd_reader
implicit none

private
public :: traj_block_type
public :: adaptive_stream_type
public :: adaptive_init_stream
public :: adaptive_preload
public :: adaptive_advance
public :: adaptive_get_block
public :: stream_finished

!=========================================================
type :: traj_block_type
    integer :: natoms  = 0
    integer :: nframes = 0
    real, allocatable :: pos(:,:,:)
end type

!=========================================================
type :: adaptive_stream_type

    type(traj_block_type) :: block(2)

    integer :: cur = 1
    integer :: nxt = 2

    integer :: istate = 0
    integer :: irep   = 0

    integer :: ifile  = 1
    integer :: iframe = 1

    integer :: blocksize = 32768

    integer :: ignore_first = 0
    integer :: ignore_last  = 0

    logical :: finished = .false.

end type

contains
subroutine adaptive_init_stream(stream, istate, irep, blocksize)

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

    stream%block(1)%nframes = 0
    stream%block(2)%nframes = 0

end subroutine

subroutine adaptive_stream_swap(stream)

    type(adaptive_stream_type), intent(inout) :: stream
    integer :: tmp

    tmp = stream%cur
    stream%cur = stream%nxt
    stream%nxt = tmp

end subroutine

!!=========================================================
!function adaptive_get_block(stream) result(b)
!!=========================================================
!    type(adaptive_stream_type), intent(in) :: stream
!    type(traj_block_type) :: b
!
!    b = stream%block(stream%cur)
!end function

subroutine adaptive_preload(stream, project)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    call fill_block(stream, project, stream%block(stream%cur))

    if (.not. stream%finished) then
        call fill_block(stream, project, stream%block(stream%nxt))
    end if

end subroutine


!=========================================================
subroutine adaptive_advance(stream, project)
!=========================================================
    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    if (stream%finished) return

    call adaptive_stream_swap(stream)
    call fill_block(stream, project, stream%block(stream%nxt))

end subroutine

function adaptive_get_block(stream) result(b)

    type(adaptive_stream_type), intent(in) :: stream
    type(traj_block_type) :: b

    b = stream%block(stream%cur)

end function
!=========================================================
subroutine fill_block(stream, project, block)
!=========================================================
    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project
    type(traj_block_type), intent(inout) :: block

    type(dcd_info_type) :: info
    type(adaptive_dcd_type) :: reader

    character(400) :: fname
    integer :: nread

    block%nframes = 0

    do while (block%nframes < stream%blocksize)

        if (stream%ifile > project%state(stream%istate)%rep(stream%irep)%ndcd) then
            stream%finished = .true.
            exit
        end if

        fname = project%state(stream%istate)%rep(stream%irep)%dcd_files(stream%ifile)

        call adaptive_open(reader, fname)

        info%natom = reader%natom

        if (.not. allocated(block%pos)) then
            allocate(block%pos(3, info%natom, stream%blocksize))
            block%natoms = info%natom
        end if

        !call adaptive_read_frames(reader, block%pos, nread)

        block%nframes = block%nframes + nread

        call adaptive_close(reader)

        stream%ifile = stream%ifile + 1

        if (block%nframes >= stream%blocksize) exit

    end do

end subroutine

!=========================================================
function stream_finished(stream) result(f)
!=========================================================
    type(adaptive_stream_type), intent(in) :: stream
    logical :: f

    f = stream%finished
end function

end module adaptive_stream_module