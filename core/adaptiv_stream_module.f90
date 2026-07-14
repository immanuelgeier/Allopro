module adaptive_stream_module

use project_module
use dcd_reader

implicit none

!=========================================================
! DATA BLOCK
!=========================================================
type :: traj_block_type

    integer :: natoms  = 0
    integer :: nframes = 0

    real, allocatable :: pos(:,:,:)

    character(400) :: first_file = ''
    character(400) :: last_file  = ''

end type

!=========================================================
! STREAM OBJECT
!=========================================================
type :: adaptive_stream_type

    ! physical storage
    type(traj_block_type), target :: block1
    type(traj_block_type), target :: block2

    ! active pointers
    type(traj_block_type), pointer :: A => null()
    type(traj_block_type), pointer :: B => null()

    integer :: istate = 0
    integer :: irep   = 0

    integer :: ifile  = 1
    integer :: iframe = 1

    integer :: blocksize = 32768

    logical :: finished = .false.

end type

contains

!=========================================================
subroutine init_stream(stream, istate, irep, blocksize)
!=========================================================

    type(adaptive_stream_type), intent(inout) :: stream

    integer, intent(in) :: istate
    integer, intent(in) :: irep
    integer, intent(in) :: blocksize

    stream%istate = istate
    stream%irep   = irep

    stream%ifile  = 1
    stream%iframe = 1

    stream%blocksize = blocksize

    stream%A => stream%block1
    stream%B => stream%block2

    stream%finished = .false.

    stream%A%nframes = 0
    stream%B%nframes = 0

end subroutine

!=========================================================
subroutine preload(stream, project)
!=========================================================

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    call fill_block(stream, project, stream%A)

    if (.not. stream%finished) then
        call fill_block(stream, project, stream%B)
    end if

end subroutine

!=========================================================
subroutine swap(stream)
!=========================================================

    type(adaptive_stream_type), intent(inout) :: stream

    type(traj_block_type), pointer :: tmp

    tmp      => stream%A
    stream%A => stream%B
    stream%B => tmp

end subroutine

!=========================================================
subroutine advance(stream, project)
!=========================================================

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project

    if (stream%finished) return

    call swap(stream)

    call fill_block(stream, project, stream%B)

end subroutine

!=========================================================
subroutine fill_block(stream, project, block)
!=========================================================

    type(adaptive_stream_type), intent(inout) :: stream
    type(project_type), intent(in) :: project
    type(traj_block_type), intent(inout) :: block

    integer :: natom
    integer :: nframes_file
    integer :: need
    integer :: remaining

    character(400) :: fname

    block%nframes = 0

    do while (block%nframes < stream%blocksize)

        if (stream%ifile > &
            project%state(stream%istate)%rep(stream%irep)%ndcd) then

            stream%finished = .true.
            return

        end if

        fname = project%state(stream%istate)%rep(stream%irep)% &
                dcd_files(stream%ifile)

        call read_dcd_header(fname, natom, nframes_file)

        if (.not. allocated(block%pos)) then

            allocate( &
                block%pos(3, natom, stream%blocksize))

            block%natoms = natom

        end if

        remaining = nframes_file - stream%iframe + 1

        need = stream%blocksize - block%nframes

        if (remaining <= need) then

            !-----------------------------------------
            ! consume rest of DCD
            !-----------------------------------------

            call read_dcd_into_block_partial( &
                 fname, &
                 stream%iframe, &
                 nframes_file, &
                 block%pos, &
                 block%nframes + 1)

            block%nframes = block%nframes + remaining

            stream%ifile  = stream%ifile + 1
            stream%iframe = 1

        else

            !-----------------------------------------
            ! consume only part of DCD
            !-----------------------------------------

            call read_dcd_into_block_partial( &
                 fname, &
                 stream%iframe, &
                 stream%iframe + need - 1, &
                 block%pos, &
                 block%nframes + 1)

            block%nframes = stream%blocksize

            stream%iframe = stream%iframe + need

            return

        end if

    end do

end subroutine

end module adaptive_stream_module