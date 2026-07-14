module file_stream_module
use project_module
use dcd_reader
implicit none

type :: file_stream_type
    integer :: istate = 0
    integer :: irep   = 0
    integer :: ifile  = 1
    integer :: iframe = 1
    logical  :: finished = .false.
end type

contains

subroutine file_stream_init(s, istate, irep)
    type(file_stream_type), intent(inout) :: s
    integer, intent(in) :: istate, irep

    s%istate = istate
    s%irep   = irep
    s%ifile  = 1
    s%iframe = 1
    s%finished = .false.
end subroutine

subroutine file_stream_load_next(s, project, block)
    type(file_stream_type), intent(inout) :: s
    type(project_type), intent(in) :: project
    type(traj_block_type), intent(inout) :: block

    character(400) :: fname

    if (s%ifile > project%state(s%istate)%rep(s%irep)%ndcd) then
        s%finished = .true.
        return
    end if

    fname = project%state(s%istate)%rep(s%irep)%dcd_files(s%ifile)
    call read_dcd_into_block(fname, block%pos, block%natoms, block%nframes)

    s%ifile = s%ifile + 1
end subroutine

end module