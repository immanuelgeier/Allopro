module healpix_type_module
    implicit none
    private
    
    public :: refinement_lut, compressed_tau_histogram, &
              healpix_block_trajectory, healpix_block_ensemble

    ! Repräsentiert die fertig gelifteten, dichten Integer-Paare eines EINZELNEN Blocks
    type :: healpix_block_trajectory
        integer :: nobject
        integer :: npairs                          ! P_b: Anzahl valider Paare (count(mask_Y)/2)
        
        ! Real-Modus (Synchron)
        integer, allocatable :: Y_past(:,:)        ! (nobject, npairs) -> i_y
        integer, allocatable :: Y_future(:,:)      ! (nobject, npairs) -> i_x
        integer, allocatable :: Z_past_real(:,:)   ! (nobject, npairs) -> i_z (Echtlauf)
        
        ! Shuffle-Modus (Surrogat, aus fremdem Block b_shuffle)
        integer, allocatable :: Z_past_shuffle(:,:)  ! (nobject, npairs) -> i_z (Surrogat)
        integer, allocatable :: Z_future_shuffle(:,:) ! Für optionale erweiterte Analysen
    end type healpix_block_trajectory

    ! Ein simpler Container, der die gesamte Trajektorie als Liste kompakter Integer-Blöcke hält
    type :: healpix_block_ensemble
        integer :: nblocks
        type(healpix_block_trajectory), allocatable :: blocks(:) ! (nblocks)
    end type healpix_block_ensemble
    
    type :: refinement_lut
        integer :: nside_coarse
        integer :: nside_fine
        integer :: npix_coarse
        integer :: npix_fine
        integer :: sparse_size           
        logical, allocatable :: coarse_occupied(:,:,:)
        integer, allocatable :: fine_to_sparse(:,:,:)
    end type refinement_lut
        
    type :: compressed_tau_histogram
        integer(4), allocatable :: raw_counts(:,:,:)
    end type compressed_tau_histogram

end module healpix_type_module