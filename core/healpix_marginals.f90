module healpix_marginals_module
    use healpix_types
    implicit none
    private

    public :: healpix_marginals, calculate_marginals, allocate_marginals, deallocate_marginals

    type :: healpix_marginals
        integer :: npix
        integer :: nobject

        ! Zweidimensionale Marginalverteilungen (Paarbeziehungen behalten 4D)
        integer, allocatable :: xy(:,:,:,:) 
        integer, allocatable :: xz(:,:,:,:) 
        integer, allocatable :: yz(:,:,:,:) ! Wichtig für TE

        ! Eindimensionale Marginalverteilungen (Jetzt sauber entkoppelt als 2D)
        integer, allocatable :: x(:,:)    ! (npix_fine, nobject) -> Zukunft von Komponente i
        integer, allocatable :: y(:,:)    ! (npix_fine, nobject) -> Vergangenheit von Komponente i
        integer, allocatable :: z(:,:)    ! (npix_fine, nobject) -> Vergangenheit von Quellkomponente j
    end type healpix_marginals

contains

    subroutine allocate_marginals(marg, npix, nobject)
        type(healpix_marginals), intent(out) :: marg
        integer, intent(in)                  :: npix
        integer, intent(in)                  :: nobject

        marg%npix = npix
        marg%nobject   = nobject

        allocate(marg%xy(0:npix-1, 0:npix-1, nobject, nobject), source=0)
        allocate(marg%xz(0:npix-1, 0:npix-1, nobject, nobject), source=0)
        allocate(marg%yz(0:npix-1, 0:npix-1, nobject, nobject), source=0)

        allocate(marg%x(0:npix-1, nobject), source=0)
        allocate(marg%y(0:npix-1, nobject), source=0)
        allocate(marg%z(0:npix-1, nobject), source=0)
    end subroutine allocate_marginals

    subroutine deallocate_marginals(marg)
        type(healpix_marginals), intent(inout) :: marg
        if (allocated(marg%xy)) deallocate(marg%xy)
        if (allocated(marg%xz)) deallocate(marg%xz)
        if (allocated(marg%yz)) deallocate(marg%yz)
        if (allocated(marg%x))  deallocate(marg%x)
        if (allocated(marg%y))  deallocate(marg%y)
        if (allocated(marg%z))  deallocate(marg%z)
    end subroutine deallocate_marginals

    subroutine calculate_marginals(raw_counts, fine_to_sparse, npix, marg)
        integer(4), intent(in)                 :: raw_counts(:,:,:)
        integer, intent(in)                    :: fine_to_sparse(0:, 0:, 0:) 
        integer, intent(in)                    :: npix
        type(healpix_marginals), intent(inout) :: marg

        integer :: ix, iy, iz, idx_sparse
        integer :: i, j, count_val

        ! Vor der Berechnung einmal sauber nullen, um Reste zu vermeiden
        if (allocated(marg%xy)) marg%xy = 0
        if (allocated(marg%xz)) marg%xz = 0
        if (allocated(marg%yz)) marg%yz = 0
        if (allocated(marg%x))  marg%x  = 0
        if (allocated(marg%y))  marg%y  = 0
        if (allocated(marg%z))  marg%z  = 0

        ! Schleifen laufen über die HEALPix-Koordinaten (0 bis npix_fine - 1)
        do iz = 0, npix - 1
            do iy = 0, npix - 1
                do ix = 0, npix - 1
                    
                    idx_sparse = fine_to_sparse(ix, iy, iz)
                    if (idx_sparse <= 0) cycle 

                    do j = 1, marg%nobject
                        do i = 1, marg%nobject
                            count_val = raw_counts(idx_sparse, i, j)
                            if (count_val == 0) cycle

                            ! 2D-Verteilungen (Indizes verschoben um +1 wegen Fortran 1-based indexing)
                            marg%xy(ix, iy, i, j) = marg%xy(ix, iy, i, j) + count_val
                            marg%xz(ix, iz, i, j) = marg%xz(ix, iz, i, j) + count_val
                            marg%yz(iy, iz, i, j) = marg%yz(iy, iz, i, j) + count_val

                            ! 1D-Verteilungen: Getrennt nach physikalischer Komponente!
                            marg%x(ix, i) = marg%x(ix, i) + count_val  ! Zukunft von Ziel i
                            marg%y(iy, i) = marg%y(iy, i) + count_val  ! Vergangenheit von Ziel i
                            marg%z(iz, j) = marg%z(iz, j) + count_val  ! Vergangenheit von Quelle j
                        end do
                    end do

                end do
            end do
        end do
    end subroutine calculate_marginals

end module healpix_marginals_module
