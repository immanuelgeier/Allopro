module healpix_analysis_execution_module_single_lag
    use time_window_file_mapping
    use time_window_loader 
    use healpix_sampling_plan_module
    use position_lift_module             ! <<< NEU eingebunden!
    use direction_vectors_analyse_plan_module ! Enthält project_type & cv_plan
    use healpix_types
    use pix_tools

    implicit none
    private
    public :: refinement_lut, compressed_tau_histogram, stream_and_sample_tau

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

contains
subroutine execute_healpix_analysis(project, cv_p, istate, L, mode, selector, tau)
    type(project_type), intent(in)            :: project
    type(cv_plan), intent(in)                 :: cv_p
    integer, intent(in)                       :: istate
    integer, intent(in)                       :: L
    integer, intent(in)                       :: mode
    type(replica_selector_type), intent(in)   :: selector
    integer, intent(in)                       :: tau

    ! Interne Strukturen, die die Haupt-Arbeitsschleife braucht
    type(sampling_plan)             :: plan
    type(refinement_lut)            :: lut
    type(compressed_tau_histogram)  :: hist

    ! 1. Setze die Pläne und die LUT intern auf
    call init_sampling_plan(plan, project, istate, L, mode, selector)
    call setup_refinement_lut(lut, project, istate) ! Falls du so eine Routine hast
    call allocate_compressed_histogram(hist, lut)

    ! 2. Der eigentliche Stream-Aufruf ist jetzt intern und "clean"
    ! Wir können 'is_shuffle=.false.' (Echtlauf) und '.true.' (Surrogat) nacheinander jagen
    !call stream_and_sample_tau(project, cv_p, plan, lut, hist, tau, is_shuffle=.false.)
    
    ! 3. Hier könntest du das fertige Histogramm direkt wegschreiben
    !call save_histogram(hist, istate, tau)

end subroutine execute_healpix_analysis

subroutine stream_and_sample_tau(project, cv_p, plan, lut, hist, tau, is_shuffle)
    type(project_type), intent(in)                :: project
    type(cv_plan), intent(in)                     :: cv_p
    type(sampling_plan), intent(in)               :: plan
    type(refinement_lut), intent(in)              :: lut
    type(compressed_tau_histogram), intent(inout) :: hist
    integer, intent(in)                           :: tau
    logical, intent(in)                           :: is_shuffle
    
    integer :: b_idx, target_block_id, source_block_id
    type(block_meta) :: target_meta, source_meta
    
    type(time_window_type) :: win_Y, win_Z
    logical, allocatable   :: mask_Y(:), mask_Z(:)
    type(raw_data_block)   :: raw_Y, raw_Z
    type(lifted_trajectory) :: data_Y, data_Z
    
    real(8) :: rot(3,3), v_rot(3)
    integer :: t, ipix_x, ipix_y, ipix_z, idx_sparse, i, j

    do b_idx = 1, plan%nblocks_total
        
        ! Block-Zuweisung via Sampling-Plan (wie gehabt)
        target_block_id = plan%real_run_perm(b_idx)
        target_meta     = plan%global_blocks(target_block_id)
        
        if (is_shuffle) then
            source_block_id = plan%shuffle_run_perm(b_idx)
        else
            source_block_id = target_block_id
        end if
        source_meta = plan%global_blocks(source_block_id)
        
        ! Window-Grenzen setzen
        win_Y%istate = plan%istate
        win_Y%irep   = target_meta%replica_id
        win_Y%start_time_step = target_meta%start_frame
        win_Y%end_time_step   = target_meta%start_frame + target_meta%length - 1
        
        win_Z%istate = plan%istate
        win_Z%irep   = source_meta%replica_id
        win_Z%start_time_step = source_meta%start_frame
        win_Z%end_time_step   = source_meta%start_frame + source_meta%length - 1
        
        ! Masken für genau DIESEN Block allokieren und belegen
        allocate(mask_Y(target_meta%length))
        allocate(mask_Z(source_meta%length))
        mask_Y = .false.; mask_Z = .false.
        
        ! Deine alte logische Verknüpfung: Setze True bei i und i+tau
        ! (WICHTIG: Das muss für Ziel und Quelle identisch getaktet sein, damit das Shuffling synchron bleibt!)
        do t = 1, target_meta%length
            if (t + tau > target_meta%length) exit
            mask_Y(t)       = .true.
            mask_Y(t + tau) = .true.
            
            mask_Z(t)       = .true.
            mask_Z(t + tau) = .true.
        end do
        
        ! Genauso wie in deiner alten Routine: Übergabe an den Loader inkl. Maske
        call load_raw_data_block(project, win_Y, mask_Y, raw_Y)
        call load_raw_data_block(project, win_Z, mask_Z, raw_Z)
        
        ! Konvertierung in CV-Vektoren (b1, b2)
        call lift_cv_plan(cv_p, raw_Y, data_Y)
        call lift_cv_plan(cv_p, raw_Z, data_Z)
        
        call free_raw_data_block(raw_Y)
        call free_raw_data_block(raw_Z)
        
        ! Basen-Berechnung für Y
        call extract_backbone_bases(data_Y, bidx=1) 
        
        !-----------------------------------------------------------------
        ! REUZIERTE SCHLEIFE: Da der Loader nur geladen hat, was TRUE war,
        ! laufen wir jetzt exakt die gültigen Intervalle ab.
        !-----------------------------------------------------------------
        do t = 1, target_meta%length - tau
            ! Nur rechnen, wenn an diesem Zeitschritt der t und t+tau Block aktiv geladen wurde
            if (.not. mask_Y(t)) cycle 
            
            rot = data_Y%rot_matrix(:,:,t)
            
            do i = 1, data_Y%nobject
                v_rot = matmul(rot, data_Y%vec(:, i, t, 1))
                call vec2pix_ring(lut%nside_fine, v_rot, ipix_y)
                
                v_rot = matmul(rot, data_Y%vec(:, i, t + tau, 1))
                call vec2pix_ring(lut%nside_fine, v_rot, ipix_x)
                
                do j = 1, data_Z%nobject
                    ! z(t) kommt synchron aus dem (shuffled) Z-Block!
                    v_rot = matmul(rot, data_Z%vec(:, j, t, 1))
                    call vec2pix_ring(lut%nside_fine, v_rot, ipix_z)
                    
                    idx_sparse = lut%fine_to_sparse(ipix_x, ipix_y, ipix_z)
                    if (idx_sparse > 0) then
                        hist%raw_counts(idx_sparse, i, j) = hist%raw_counts(idx_sparse, i, j) + 1
                    end if
                end do 
            end do 
        end do 
        
        deallocate(mask_Y, mask_Z)
        call free_lifted_trajectory(data_Y)
        call free_lifted_trajectory(data_Z)
        
    end do
end subroutine stream_and_sample_tau


end module healpix_analysis_execution_module_single_lag