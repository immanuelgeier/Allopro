module pointgroup_sym_module

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    type :: point_group_type
        character(len=20) :: name
        integer :: nsym
        real(dp), allocatable :: R(:,:,:)
        ! (3,3,nsym)
    end type
    
    type :: symmetry_operation
    character(len=8) :: label
    real(dp) :: R(3,3)
    end type

contains

subroutine make_Cn(n, axis, pg)
    use symmetry_operations
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: axis(3)
    type(point_group_type), intent(out) :: pg
    integer :: i
    real(dp) :: theta
   
    write(pg%name,'("C",I0)') n
    
    pg%nsym = n
   
    allocate(pg%R(3,3,n))
  
    do i = 0, n-1
        theta = 2.0_dp*pi*real(i,dp)/real(n,dp)
        pg%R(:,:,i+1) = rotation_matrix(axis,theta)
    end do
end subroutine make_Cn

    function get_point_group(name) result(pg)

        character(len=*), intent(in) :: name
        type(point_group_type) :: pg

        select case(trim(name))

        case("C2")
            pg%name = "C2"
            pg%nsym = 2
            allocate(pg%R(3,3,2))

            ! E
            pg%R(:,:,1) = reshape([ &
                1d0,0d0,0d0, &
                0d0,1d0,0d0, &
                0d0,0d0,1d0 ], [3,3])

            ! C2 um z
            pg%R(:,:,2) = reshape([ &
               -1d0,0d0,0d0, &
                0d0,-1d0,0d0, &
                0d0,0d0,1d0 ], [3,3])

        case("C3")
            pg%name = "C3"
            pg%nsym = 3
            allocate(pg%R(3,3,3))

            ! hier später E, C3 und C3²

        case default
            stop "Unknown point group."

        end select

    end function

end module