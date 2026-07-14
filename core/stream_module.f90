module stream_module

use iso_fortran_env, only: real32
implicit none

!=========================================================
! STREAM STATE
!=========================================================
type :: stream_state_type

    ! ---- system size ----
    integer :: natom_full = 0
    integer :: natom_use  = 0

    ! ---- buffers (double buffering) ----
    real(real32), allocatable :: posA(:,:)
    real(real32), allocatable :: posB(:,:)

    logical :: useA = .true.

    ! ---- masks ----
    logical, allocatable :: atom_mask(:)
    integer, allocatable :: frame_idx(:)

    ! ---- iteration state ----
    integer :: current_frame = 0
    integer :: iframe_global = 0

    ! ---- trajectory info ----
    integer :: ndcd = 0
    integer :: idcd = 1
    character(400), allocatable :: dcd_files(:)

end type


contains

!=========================================================
! INIT STREAM
!=========================================================
subroutine stream_init(stream, natom, dcd_files, frame_idx, atom_mask)

    type(stream_state_type), intent(inout) :: stream
    integer, intent(in) :: natom
    character(*), intent(in) :: dcd_files(:)
    integer, intent(in) :: frame_idx(:)
    logical, intent(in) :: atom_mask(:)

    integer :: nsel

    stream%natom_full = natom
    stream%ndcd = size(dcd_files)

    allocate(stream%dcd_files(stream%ndcd))
    stream%dcd_files = dcd_files

    allocate(stream%frame_idx(size(frame_idx)))
    stream%frame_idx = frame_idx

    allocate(stream%atom_mask(natom))
    stream%atom_mask = atom_mask

    nsel = count(atom_mask)
    stream%natom_use = nsel

    allocate(stream%posA(3, nsel))
    allocate(stream%posB(3, nsel))

    stream%useA = .true.
    stream%current_frame = 0
    stream%iframe_global = 0
    stream%idcd = 1

end subroutine


!=========================================================
! NEXT FRAME
!=========================================================
subroutine stream_next(stream, pos_out, ok)

    type(stream_state_type), intent(inout) :: stream
    real(real32), intent(out) :: pos_out(:,:)
    logical, intent(out) :: ok

    real(real32), pointer :: buf(:,:)
    integer :: target_frame

    ok = .false.

    if (stream%current_frame >= size(stream%frame_idx)) return

    stream%current_frame = stream%current_frame + 1
    target_frame = stream%frame_idx(stream%current_frame)

    ! ---- select buffer ----
    if (stream%useA) then
        buf => stream%posA
    else
        buf => stream%posB
    end if

    ! ---- read frame (external DCD backend) ----
    call read_dcd_frame( &
        stream%dcd_files(stream%idcd), &
        target_frame, &
        buf )

    ! ---- apply atom mask ----
    call apply_atom_mask(buf, pos_out, stream%atom_mask)

    stream%useA = .not. stream%useA
    stream%iframe_global = stream%iframe_global + 1

    ok = .true.

end subroutine


!=========================================================
! ATOM MASK
!=========================================================
subroutine apply_atom_mask(pos_in, pos_out, mask)

    real(real32), intent(in)  :: pos_in(:,:)
    real(real32), intent(out) :: pos_out(:,:)
    logical, intent(in)       :: mask(:)

    integer :: i, j

    j = 0
    do i = 1, size(mask)
        if (mask(i)) then
            j = j + 1
            pos_out(:, j) = pos_in(:, i)
        end if
    end do

end subroutine


!=========================================================
! CLOSE (optional cleanup)
!=========================================================
subroutine stream_close(stream)

    type(stream_state_type), intent(inout) :: stream

    if (allocated(stream%posA)) deallocate(stream%posA)
    if (allocated(stream%posB)) deallocate(stream%posB)
    if (allocated(stream%atom_mask)) deallocate(stream%atom_mask)
    if (allocated(stream%frame_idx)) deallocate(stream%frame_idx)
    if (allocated(stream%dcd_files)) deallocate(stream%dcd_files)

    stream%current_frame = 0
    stream%iframe_global = 0

end subroutine


end module stream_module