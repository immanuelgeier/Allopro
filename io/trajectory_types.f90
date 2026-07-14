module trajectory_types
implicit none

type :: dcd_file_info_type
    character(400) :: filename = ''
    integer :: natoms_dcd = 0
    integer :: nframes_raw = 0
end type


type :: dcd_file_view_type
    integer :: first_frame = 1
    integer :: last_frame  = 0

    integer :: nframes_valid = 0
end type


type :: traj_block_type
    integer :: natoms_loaded = 0
    integer :: nframes_loaded = 0
    real(8), allocatable :: pos(:,:,:)
end type

end module