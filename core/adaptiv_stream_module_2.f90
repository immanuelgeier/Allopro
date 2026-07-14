module adaptive_stream_module

use project_module
use adaptive_dcd_reader
implicit none

private

public :: traj_block_type
public :: adaptive_stream_type
public :: init_stream
public :: preload
public :: advance
public :: get_block
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

    integer :: blocksize = 0

    integer :: ignore_first = 0
    integer :: ignore_last  = 0

    logical :: finished = .false.

end type

contains

!=========================================================
subroutine init_stream(stream, istate, irep, blocksize, ignore_first, ignore_last)

    type(adaptive_stream_type), intent(inout) :: stream
    integer, intent(in) :: istate, irep, blocksize
    integer, intent(in), optional :: ignore_first, ignore_last

    stream%istate = istate
    stream%irep   = irep
    stream%blocksize = blocksize

    stream%ignore_first = 0
    stream%ignore_last  = 0

    if (present(ignore_first)) stream%ignore_first = ignore_first
    if (present(ignore_last))  stream%ignore_last  = ignore_last

    stream%ifile  = 1
    stream%iframe = 1

    stream%cur = 1
    stream%nxt = 2

    stream%block(1)%nframes = 0
    stream%block(2)%nframes = 0

    stream%finished = .false.

end subroutine

!=========================================================
subroutine swap(stream)
    type(adaptive_stream_type), intent(inout) :: stream
    integer :: tmp

    tmp = stream%cur
    stream%cur = stream%nxt
    stream%nxt = tmp
end subroutine

!=========================================================
function get_block(stream) result(b)
    type(adaptive_stream_type), intent(in) :: stream
    type(traj_block_type) :: b

    b = stream%block(stream%cur)
end function

!=========================================================
subroutine preload(stream, project)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    call fill_block(stream, project, stream%block(stream%cur))

    if (.not. stream%finished) then
        call fill_block(stream, project, stream%block(stream%nxt))
    end if

end subroutine

!=========================================================
subroutine advance(stream, project)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    if (stream%finished) return

    call swap(stream)
    call fill_block(stream, project, stream%block(stream%nxt))

end subroutine

!=========================================================
subroutine fill_block(stream, project, block)

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project
    type(traj_block_type), intent(inout) :: block

    type(dcd_info_type) :: info
    type(adaptive_dcd_type) :: reader

    character(400) :: fname
    integer :: nread, need
    integer :: first_valid, last_valid
    integer :: natom, nframes_file

    block%nframes = 0

    do while (block%nframes < stream%blocksize)

        if (stream%ifile > project%state(stream%istate)%rep(stream%irep)%ndcd) then
            stream%finished = .true.
            exit
        end if

        fname = project%state(stream%istate)%rep(stream%irep)%dcd_files(stream%ifile)

        call adaptive_open(reader, fname)

        natom        = reader%natom
        nframes_file = reader%nframes

        first_valid = stream%ignore_first + 1
        last_valid  = nframes_file - stream%ignore_last

        if (.not. allocated(block%pos)) then
            allocate(block%pos(3, natom, stream%blocksize))
            block%natoms = natom
        end if

        ! skip if already past file
        if (stream%iframe < first_valid) stream%iframe = first_valid

        if (stream%iframe > last_valid) then
            stream%ifile = stream%ifile + 1
            stream%iframe = 1
            call adaptive_close(reader)
            cycle
        end if

        need = stream%blocksize - block%nframes

        nread = min(need, last_valid - stream%iframe + 1)

        call adaptive_skip_to_frame(reader, stream%iframe)

        call adaptive_read_frames(reader, &
                                   block%pos(:,:,block%nframes+1), &
                                   nread)

        block%nframes = block%nframes + nread
        stream%iframe = stream%iframe + nread

        call adaptive_close(reader)

        if (block%nframes >= stream%blocksize) exit

        if (stream%iframe > last_valid) then
            stream%ifile = stream%ifile + 1
            stream%iframe = 1
        end if

    end do

end subroutine

!=========================================================
function stream_finished(stream) result(f)
    type(adaptive_stream_type), intent(in) :: stream
    logical :: f

    f = stream%finished
end function

end module adaptive_stream_module