module lifted_trajectory_io_module
use position_lift_module


    implicit none

    private
    public :: lifted_trajectory
    public :: write_trajectory
    public :: read_trajectory

    !type :: lifted_trajectory
    !    integer :: nobject
    !    integer :: ntimes
    !
    !    ! Dimensionen: (3, 2*nobject, ntimes)
    !    real(8), allocatable :: vec(:,:,:)
    !end type lifted_trajectory


contains


subroutine write_trajectory(filename, traj)
    character(len=*), intent(in) :: filename
    type(lifted_trajectory), intent(in) :: traj

    integer :: unit

    open(newunit=unit, &
         file=filename, &
         access="stream", &
         form="unformatted", &
         status="replace")

    ! Header schreiben
    write(unit) traj%nobject
    write(unit) traj%ntimes

    ! Daten schreiben
    write(unit) traj%vec

    close(unit)

end subroutine write_trajectory

subroutine write_trajectory_centroid(filename, traj)
    character(len=*), intent(in) :: filename
    type(lifted_trajectory), intent(in) :: traj

    integer :: unit

    open(newunit=unit, &
         file=filename, &
         access="stream", &
         form="unformatted", &
         status="replace")

    ! Header schreiben
    write(unit) traj%nobject
    write(unit) traj%ntimes

    ! Daten schreiben
    write(unit) traj%centroid

    close(unit)

end subroutine write_trajectory_centroid



subroutine read_trajectory(filename, traj)
    character(len=*), intent(in) :: filename
    type(lifted_trajectory), intent(out) :: traj

    integer :: unit
    integer :: n1, n2

    open(newunit=unit, &
         file=filename, &
         access="stream", &
         form="unformatted", &
         status="old")


    ! Header lesen
    read(unit) n1
    read(unit) n2
    
    write(0,*) trim(filename),2*n1,n2
    !pause
    
    traj%nobject = n1
    traj%ntimes  = n2

    ! Speicher reservieren
    allocate(traj%vec(3,2*n1,n2))

    ! Daten lesen
    read(unit) traj%vec

    close(unit)

end subroutine read_trajectory

subroutine read_trajectory_centroid(filename, traj)
    character(len=*), intent(in) :: filename
    type(lifted_trajectory), intent(out) :: traj

    integer :: unit
    integer :: n1, n2

    open(newunit=unit, &
         file=filename, &
         access="stream", &
         form="unformatted", &
         status="old")

    ! Header lesen
    read(unit) n1
    read(unit) n2

    traj%nobject = n1
    traj%ntimes  = n2

    ! Speicher reservieren
    allocate(traj%centroid(3,n1,n2))

    ! Daten lesen
    read(unit) traj%centroid

    close(unit)

end subroutine read_trajectory_centroid


end module 