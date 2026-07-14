module healpix_analysis_execution_module

use time_window_loader
use direction_vectors_analyse_plan_module
use position_lift_module
use healpix_types
use pix_tools

implicit none

! ======================================================================
! 1. STRUKTUREN FÜR ZEIT- UND SHUFFLE-PLANUNG (SINGLE-TAU)
! ======================================================================

type :: block_shuffle_plan
    integer :: nblocks               ! Anzahl der Blöcke in der Trajektorie
    integer :: block_length          ! Länge L eines Blocks (Frames)
    integer, allocatable :: perm(:)  ! Die blockweise Permutation für den Shuffle-Lauf
end type block_shuffle_plan

! ======================================================================
! 2. STRUKTUREN FÜR HIERARCHISCHE INDEX-KOMPRESSION (LUT)
! ======================================================================

type :: refinement_lut
    integer :: nside_coarse
    integer :: nside_fine
    integer :: npix_coarse
    integer :: npix_fine
    integer :: sparse_size           ! Gesamtanzahl allokierter feiner Bins
    
    ! Grobes Besetzungsnetz: wahr, wenn das Grob-Pixel-Tupel bei diesem tau besetzt ist
    logical, allocatable :: coarse_occupied(:,:,:)
    
    ! Look-Up-Table für die feinen Pixel (Wert: -1 oder fortlaufender Index)
    integer, allocatable :: fine_to_sparse(:,:,:)
end type refinement_lut

! Das speichereffiziente Histogramm-Array (Jetzt ohne nlags!)
type :: compressed_tau_histogram
    ! Dimensionen: (sparse_size, nobject, nobject)
    integer(4), allocatable :: raw_counts(:,:,:)
end type compressed_tau_histogram

contains

! ======================================================================
! ROTATIONS-BASIS EXTRAKTION
! ======================================================================

subroutine extract_backbone_bases(lifted, bidx)
    type(lifted_trajectory), intent(inout) :: lifted
    integer, intent(in)                    :: bidx
    integer :: t
    real(8) :: b1(3), b2(3), b3(3)

    if (.not. allocated(lifted%rot_matrix)) then
        allocate(lifted%rot_matrix(3, 3, lifted%ntimes))
    end if

    do t = 1, lifted%ntimes
        b1 = lifted%vec(:, bidx, t, 1)
        b2 = lifted%vec(:, bidx, t, 2)

        b1 = b1 / norm2(b1)
        b2 = b2 - dot_product(b1, b2) * b1
        b2 = b2 / norm2(b2)
        b3 = cross_product(b1, b2) 

        ! Speicherung als Zeilen für die spätere Multiplikation: v_rot = matmul(rot, v)
        lifted%rot_matrix(1, :) = b1
        lifted%rot_matrix(2, :) = b2
        lifted%rot_matrix(3, :) = b3
    end do
end subroutine extract_backbone_bases

! ======================================================================
! KERN-SUBROUTINEN FÜR DEN SINGLE-TAU REFINEMENT RUN
! ======================================================================

subroutine get_initial_occupancy(lifted, tau, lut, coarse_hist)
    ! STUFE 1: Bestimme Besetzung und Markov-Kerne bei Nside=1 für ein FESTES tau
    type(lifted_trajectory), intent(in) :: lifted
    integer, intent(in)                 :: tau
    type(refinement_lut), intent(inout) :: lut
    integer, intent(inout)              :: coarse_hist(0:,0:,0:) ! (12, 12, 12) Markov-Zähler
    
    integer :: t, i, j, ipix_x, ipix_y, ipix_z
    real(8) :: rot(3,3), v_rot(3)
    
    lut%nside_coarse = 1
    lut%npix_coarse  = 12
    
    if (.not. allocated(lut%coarse_occupied)) then
        allocate(lut%coarse_occupied(0:11, 0:11, 0:11))
    end if
    lut%coarse_occupied = .false.
    coarse_hist = 0
    
    ! Da wir mit Vektoren arbeiten, projizieren wir direkt in der Schleife
    do t = 1, lifted%ntimes - tau
        rot = lifted%rot_matrix(:,:,t)
        
        do i = 1, lifted%nobject
            ! Vergangenheit Ziel y(t) in Basis t
            v_rot = matmul(rot, lifted%vec(:, i, t, 1))
            call vec2pix_ring(1, v_rot, ipix_y)
            
            ! Zukunft Ziel y(t+tau) in Basis t
            v_rot = matmul(rot, lifted%vec(:, i, t + tau, 1))
            call vec2pix_ring(1, v_rot, ipix_x)
            
            do j = 1, lifted%nobject
                ! Vergangenheit Quelle z(t) in Basis t
                v_rot = matmul(rot, lifted%vec(:, j, t, 1))
                call vec2pix_ring(1, v_rot, ipix_z)
                
                ! Belegung und Grob-Histogramm (Markov-Kern) füllen
                lut%coarse_occupied(ipix_x, ipix_y, ipix_z) = .true.
                coarse_hist(ipix_x, ipix_y, ipix_z) = coarse_hist(ipix_x, ipix_y, ipix_z) + 1
            end do
        end do
    end do
end subroutine get_initial_occupancy


subroutine prepare_fine_sampling(lut, nside_fine_target)
    ! STUFE 2: Erstellung der feinen LUT exakt für die besetzten Grob-Zellen bei diesem tau
    type(refinement_lut), intent(inout) :: lut
    integer, intent(in)                 :: nside_fine_target
    
    integer :: ix, iy, iz, cx, cy, cz, refinement_factor
    
    lut%nside_fine = nside_fine_target
    lut%npix_fine  = 12 * (lut%nside_fine)**2
    refinement_factor = lut%nside_fine / lut%nside_coarse 
    
    if (allocated(lut%fine_to_sparse)) deallocate(lut%fine_to_sparse)
    allocate(lut%fine_to_sparse(0:lut%npix_fine-1, 0:lut%npix_fine-1, 0:lut%npix_fine-1))
    
    lut%fine_to_sparse = -1
    lut%sparse_size = 0
    
    do ix = 0, lut%npix_fine-1
        do iy = 0, lut%npix_fine-1
            do iz = 0, lut%npix_fine-1
                
                ! Geometrisches Mapping der feinen RING-Pixel auf das grobe Elternpixel (Nside=1)
                cx = ix / (refinement_factor**2) 
                cy = iy / (refinement_factor**2)
                cz = iz / (refinement_factor**2)
                
                if (lut%coarse_occupied(cx, cy, cz)) then
                    lut%sparse_size = lut%sparse_size + 1
                    lut%fine_to_sparse(ix, iy, iz) = lut%sparse_size
                end if
                
            end do
        end do
    end do
end subroutine prepare_fine_sampling


subroutine execute_fine_sampling(lifted, tau, shuffle_plan, lut, hist, is_shuffle)
    ! STUFE 3: Feiner Produktionslauf für ein FESTES tau unter Nutzung der LUT
    type(lifted_trajectory), intent(in) :: lifted
    integer, intent(in)                 :: tau
    type(block_shuffle_plan), intent(in):: shuffle_plan
    type(refinement_lut), intent(in)    :: lut
    type(compressed_tau_histogram), intent(inout) :: hist
    logical, intent(in)                 :: is_shuffle
    
    integer :: t, ts, ipix_x, ipix_y, ipix_z, idx_sparse, i, j
    real(8) :: rot(3,3), v_rot(3)
    
    if (.not. allocated(hist%raw_counts)) then
        allocate(hist%raw_counts(lut%sparse_size, lifted%nobject, lifted%nobject))
    end if
    hist%raw_counts = 0
    
    do t = 1, lifted%ntimes - tau
        
        ! Block-Shuffling-Index für die Quelle berechnen falls aktiv
        if (is_shuffle) then
            ! ts = ... (Berechnung via shuffle_plan%perm)
        else
            ts = t
        end if
        
        rot = lifted%rot_matrix(:,:,t) ! Basis ist IMMER vom aktuellen Zeitpunkt t
        
        do i = 1, lifted%nobject
            ! Vergangenheit Ziel y(t) in Basis t
            v_rot = matmul(rot, lifted%vec(:, i, t, 1))
            call vec2pix_ring(lut%nside_fine, v_rot, ipix_y)
            
            ! Zukunft Ziel y(t+tau) in Basis t
            v_rot = matmul(rot, lifted%vec(:, i, t + tau, 1))
            call vec2pix_ring(lut%nside_fine, v_rot, ipix_x)
            
            do j = 1, lifted%nobject
                ! Vergangenheit Quelle z(ts) ebenfalls in Basis t projiziert!
                v_rot = matmul(rot, lifted%vec(:, j, ts, 1))
                call vec2pix_ring(lut%nside_fine, v_rot, ipix_z)
                
                idx_sparse = lut%fine_to_sparse(ipix_x, ipix_y, ipix_z)
                if (idx_sparse > 0) then
                    hist%raw_counts(idx_sparse, i, j) = hist%raw_counts(idx_sparse, i, j) + 1
                end if
            end do
        end do
    end do
end subroutine execute_fine_sampling

end module healpix_analysis_execution_module