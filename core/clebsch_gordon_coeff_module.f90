module clebsch_gordan_module

    use iso_fortran_env, only: real64, int32

    implicit none

    private

    public :: cg_initialize
    public :: cg_finalize
    public :: clebsch_gordan
    public :: cg_size

    ! ============================================================
    ! Data type
    ! ============================================================

    type :: cg_coefficient
        integer(int32) :: l1
        integer(int32) :: m1
        integer(int32) :: l2
        integer(int32) :: m2
        integer(int32) :: l3
        integer(int32) :: m3
        real(real64)   :: value
    end type cg_coefficient

    type(cg_coefficient), allocatable :: cg_data(:)

    integer(int32) :: n_cg = 0
    integer(int32) :: lmax_cg = -1

contains

    ! ============================================================
    ! Initialize
    ! ============================================================

    subroutine cg_initialize(filename)

        character(len=*), intent(in) :: filename

        integer :: unit
        integer :: ios
        integer :: i

        character(len=512) :: line

        integer(int32) :: l1, m1
        integer(int32) :: l2, m2
        integer(int32) :: l3, m3
        real(real64)   :: value

        integer(int32) :: n

        ! --------------------------------------------------------
        ! First pass: count entries
        ! --------------------------------------------------------

        open(newunit=unit, file=filename, status="old", &
             action="read", iostat=ios)

        if (ios /= 0) then
            error stop "CG: Could not open coefficient file."
        end if

        n = 0

        do
            read(unit, '(A)', iostat=ios) line

            if (ios /= 0) exit

            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            n = n + 1
        end do

        close(unit)

        if (n == 0) then
            error stop "CG: No coefficients found."
        end if

        ! --------------------------------------------------------
        ! Allocate
        ! --------------------------------------------------------

        if (allocated(cg_data)) deallocate(cg_data)

        allocate(cg_data(n))

        n_cg = n

        ! --------------------------------------------------------
        ! Second pass: read data
        ! --------------------------------------------------------

        open(newunit=unit, file=filename, status="old", &
             action="read", iostat=ios)

        if (ios /= 0) then
            error stop "CG: Could not reopen coefficient file."
        end if

        i = 0

        do
            read(unit, '(A)', iostat=ios) line

            if (ios /= 0) exit

            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            read(line, *) l1, m1, l2, m2, l3, m3, value

            i = i + 1

            cg_data(i)%l1    = l1
            cg_data(i)%m1    = m1
            cg_data(i)%l2    = l2
            cg_data(i)%m2    = m2
            cg_data(i)%l3    = l3
            cg_data(i)%m3    = m3
            cg_data(i)%value = value

        end do

        close(unit)

        ! --------------------------------------------------------
        ! Determine lmax
        ! --------------------------------------------------------

        lmax_cg = 0

        do i = 1, n_cg
            lmax_cg = max(lmax_cg, &
                          cg_data(i)%l1, &
                          cg_data(i)%l2, &
                          cg_data(i)%l3)
        end do

        write(*,'(A,I0)') &
            "CG: loaded coefficients = ", n_cg

        write(*,'(A,I0)') &
            "CG: lmax = ", lmax_cg

    end subroutine cg_initialize


    ! ============================================================
    ! Get coefficient
    ! ============================================================

    pure function clebsch_gordan(l1, m1, l2, m2, l3, m3) result(value)

        integer, intent(in) :: l1, m1
        integer, intent(in) :: l2, m2
        integer, intent(in) :: l3, m3

        real(real64) :: value

        integer :: lo, hi, mid
        integer :: cmp

        value = 0.0_real64

        ! --------------------------------------------------------
        ! Basic selection rules
        ! --------------------------------------------------------

        if (l1 < 0 .or. l2 < 0 .or. l3 < 0) return

        if (abs(m1) > l1) return
        if (abs(m2) > l2) return
        if (abs(m3) > l3) return

        if (m1 + m2 /= m3) return

        if (l3 < abs(l1-l2)) return
        if (l3 > l1+l2) return

        ! --------------------------------------------------------
        ! Binary search
        ! --------------------------------------------------------

        lo = 1
        hi = n_cg

        do while (lo <= hi)

            mid = (lo + hi) / 2

            cmp = compare_key( &
                l1, m1, l2, m2, l3, m3, &
                cg_data(mid))

            if (cmp == 0) then

                value = cg_data(mid)%value
                return

            else if (cmp < 0) then

                hi = mid - 1

            else

                lo = mid + 1

            end if

        end do

        ! Not found -> zero
        value = 0.0_real64

    end function clebsch_gordan


    ! ============================================================
    ! Compare key
    !
    ! Returns:
    !   -1 : requested key < stored key
    !    0 : equal
    !   +1 : requested key > stored key
    ! ============================================================

    pure function compare_key(l1,m1,l2,m2,l3,m3,cg) result(cmp)

        integer, intent(in) :: l1,m1,l2,m2,l3,m3
        type(cg_coefficient), intent(in) :: cg

        integer :: cmp

        cmp = 0

        if (l1 < cg%l1) then
            cmp = -1
            return
        else if (l1 > cg%l1) then
            cmp = 1
            return
        end if

        if (m1 < cg%m1) then
            cmp = -1
            return
        else if (m1 > cg%m1) then
            cmp = 1
            return
        end if

        if (l2 < cg%l2) then
            cmp = -1
            return
        else if (l2 > cg%l2) then
            cmp = 1
            return
        end if

        if (m2 < cg%m2) then
            cmp = -1
            return
        else if (m2 > cg%m2) then
            cmp = 1
            return
        end if

        if (l3 < cg%l3) then
            cmp = -1
            return
        else if (l3 > cg%l3) then
            cmp = 1
            return
        end if

        if (m3 < cg%m3) then
            cmp = -1
            return
        else if (m3 > cg%m3) then
            cmp = 1
            return
        end if

    end function compare_key


    ! ============================================================
    ! Number of coefficients
    ! ============================================================

    pure function cg_size() result(n)

        integer :: n

        n = n_cg

    end function cg_size


    ! ============================================================
    ! Finalize
    ! ============================================================

    subroutine cg_finalize()

        if (allocated(cg_data)) deallocate(cg_data)

        n_cg = 0
        lmax_cg = -1

    end subroutine cg_finalize

end module clebsch_gordan_module