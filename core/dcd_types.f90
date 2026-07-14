module dcd_types

implicit none

type :: traj_block_type
    integer :: nframes = 0
    integer :: natoms  = 0

    real(8), allocatable :: pos(:,:,:)   ! (3, natoms, frames)
    character(400) :: current_file = ''
end type traj_block_type

end module dcd_types