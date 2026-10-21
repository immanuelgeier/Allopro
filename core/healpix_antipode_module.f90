module healpix_antipode_module

    use healpix_types
    use pix_tools

    implicit none

    private

    public :: init_antipode_map
    public :: get_antipode
    public :: destroy_antipode_map

    integer(I4B), allocatable :: antipode_map(:)
    integer(I4B) :: map_nside = -1
    logical :: map_initialized = .false.

contains

    !------------------------------------------------------------------
    ! Build lookup table for antipodal NESTED pixels.
    !
    ! After this routine has been called, the antipodal pixel of ipix
    ! can be obtained very cheaply with
    !
    !     ipix_anti = get_antipode(ipix)
    !
    ! Pixel indices are 0-based.
    !------------------------------------------------------------------
    subroutine init_antipode_map(nside)

        integer(I4B), intent(in) :: nside

        integer(I4B) :: npix
        integer(I4B) :: ipix
        integer(I4B) :: ipix_anti
        real(DP)     :: vec(3)

        ! Number of HEALPix pixels
        npix = 12_I4B * nside * nside

        ! Reallocate if necessary
        if (allocated(antipode_map)) then
            deallocate(antipode_map)
        end if

        allocate(antipode_map(0:npix-1))

        !--------------------------------------------------------------
        ! Calculate antipodal pixel for every NESTED pixel.
        ! This is done only once.
        !--------------------------------------------------------------
        do ipix = 0, npix - 1

            call pix2vec_nest(nside, ipix, vec)

            ! Antipodal direction
            vec = -vec

            call vec2pix_nest(nside, vec, ipix_anti)

            antipode_map(ipix) = ipix_anti

        end do

        map_nside = nside
        map_initialized = .true.

    end subroutine init_antipode_map


    !------------------------------------------------------------------
    ! Return antipodal pixel from lookup table.
    !
    ! This routine should be extremely cheap compared to performing
    ! pix2vec_nest/vec2pix_nest every time.
    !------------------------------------------------------------------
    function get_antipode(ipix) result(ipix_anti)

        integer(I4B), intent(in) :: ipix
        integer(I4B)             :: ipix_anti,i

        if (.not. map_initialized) then
            error stop "healpix_antipode_module: map not initialized"
        end if
        i=ipix-1
        if (ipix < 0 .or. i > ubound(antipode_map, 1)) then
            error stop "healpix_antipode_module: invalid pixel index"
        end if

        ipix_anti = antipode_map(i)+1

    end function get_antipode


    !------------------------------------------------------------------
    ! Free lookup table.
    !------------------------------------------------------------------
    subroutine destroy_antipode_map()

        if (allocated(antipode_map)) then
            deallocate(antipode_map)
        end if

        map_nside = -1
        map_initialized = .false.

    end subroutine destroy_antipode_map

end module healpix_antipode_module

