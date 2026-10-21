module healpix_analysis_execution_module_single_lag
    use time_window_file_mapping
    use time_window_loader 
    use healpix_sampling_plan_module
    use position_lift_module             ! <<< NEU eingebunden!
    use direction_vectors_analyse_plan_module, only: cv_plan, project_type ! Enthält project_type & cv_plan
    use healpix_types
    use pix_tools
    use healpix_marginals_module
    use healpix_transfer_entropy_module
    use healpix_type_module
 
implicit none
    private
    public ::execute_healpix_analysis


contains




subroutine sample_fine_one_particle_hist(block_data, nside_max, hist)
    implicit none
    type(healpix_block_trajectory), intent(in) :: block_data
    integer, intent(in)                        :: nside_max
    type(fine_one_particle_hist), intent(inout):: hist

    integer :: i_obj, i_pair
    integer :: p_idx

    ! Allokierung beim ersten Block
    if (.not. allocated(hist%count_Y_past)) then
        hist%npix = 12 * nside_max**2
        hist%nobject = block_data%nobject
        
        allocate(hist%count_Y_past(0:hist%npix-1, hist%nobject), source=0)
        allocate(hist%count_Y_future(0:hist%npix-1, hist%nobject), source=0)
        allocate(hist%count_Z_real(0:hist%npix-1, hist%nobject), source=0)
        allocate(hist%count_Z_shuffle(0:hist%npix-1, hist%nobject), source=0)
        allocate(hist%count_Z_future_shuffle(0:hist%npix-1, hist%nobject), source=0)
    end if

    if (block_data%npairs <= 0) return

    ! Inkrementelles Aufsammeln über alle 5 Verteilungen
    do i_pair = 1, block_data%npairs
        do i_obj = 1, block_data%nobject
            
            ! 1. Y_past
            p_idx = block_data%Y_past(i_obj, i_pair)
            hist%count_Y_past(p_idx, i_obj) = hist%count_Y_past(p_idx, i_obj) + 1
            
            ! 2. Y_future
            p_idx = block_data%Y_future(i_obj, i_pair)
            hist%count_Y_future(p_idx, i_obj) = hist%count_Y_future(p_idx, i_obj) + 1
            
            ! 3. Z_past_real
            p_idx = block_data%Z_past_real(i_obj, i_pair)
            hist%count_Z_real(p_idx, i_obj) = hist%count_Z_real(p_idx, i_obj) + 1
            
            ! 4. Z_past_shuffle (falls allokiert)
            p_idx = block_data%Z_past_shuffle(i_obj, i_pair)
            hist%count_Z_shuffle(p_idx, i_obj) = hist%count_Z_shuffle(p_idx, i_obj) + 1


            ! 5. Z_future_shuffle (falls allokiert)
             p_idx = block_data%Z_future_shuffle(i_obj, i_pair)
             hist%count_Z_future_shuffle(p_idx, i_obj) = hist%count_Z_future_shuffle(p_idx, i_obj) + 1

        end do
    end do

end subroutine sample_fine_one_particle_hist


subroutine transform_trajectory_inplace(block_data, nside_fine, nside_coarse)
    implicit none
    type(healpix_block_trajectory), intent(inout) :: block_data
    integer, intent(in)                           :: nside_fine
    integer, intent(in)                           :: nside_coarse

    integer :: map_factor
    integer :: i_obj, i_pair

    ! 1. Einzige notwendige Prüfung: Wenn Auflösungen identisch, sofort raus
    if (nside_fine == nside_coarse) return
    if (block_data%npairs <= 0 .or. block_data%nobject <= 0) return

    ! Mapping-Faktor im NESTED-Schema berechnen (z.B. (4/1)**2 = 16)
    map_factor = (nside_fine / nside_coarse)**2

    ! 2. Reine 0-basierte Ganzzahldivision ohne Offsets

            
block_data%Y_past           = block_data%Y_past/ map_factor
block_data%Y_future         = block_data%Y_future/ map_factor
block_data%Z_past_real      = block_data%Z_past_real/ map_factor
block_data%Z_past_shuffle   = block_data%Z_past_shuffle/ map_factor
block_data%Z_future_shuffle = block_data%Z_future_shuffle / map_factor


end subroutine transform_trajectory_inplace


subroutine execute_healpix_analysis(project, cv_p, istate, L, mode, selector, tau)
    type(project_type), intent(in)            :: project
    type(cv_plan), intent(in)                 :: cv_p
    integer, intent(in)                       :: istate
    integer, intent(in)                       :: L
    integer, intent(in)                       :: mode
    type(replica_selector_type), intent(in)   :: selector
    integer, intent(in)                       :: tau
    type(fine_one_particle_hist)              :: hist_single_fine
    ! Lokale Variablen (vollständig deklariert)
    type(sampling_plan)             :: plan
    type(refinement_lut)            :: lut
    type(compressed_tau_histogram)  :: hist,hist_shuffle
    type(healpix_block_trajectory)  :: block_data  
    type(healpix_marginals)         :: marg_real, marg_shuffle
    integer                         :: b_idx,i,j,k
    character(len=2)                :: tau_str
    real(8), allocatable :: te_matrix(:,:),  te_matrix_shuffle(:,:)
    character*200::filename
    type(joint_entropy_storage) :: entropies,entropies_shuffel
    integer :: bidx
    integer :: nside_max,total_sparse_entries
    integer :: nside_hist_sample
    integer :: iost
    nside_max=2
    nside_hist_sample=1
    bidx=6

    
    
allocate(te_matrix(2*size(cv_p%object),2* size(cv_p%object)))
allocate(te_matrix_shuffle(2*size(cv_p%object),2* size(cv_p%object)))

! Nach der Real-Schleife und nach calculate_marginals:
if(1==1) then
    ! 1. Pläne, LUT und dichte Histogramm-Struktur aufsetzen
    call init_sampling_plan(plan, project, istate, L, mode, selector)
    call setup_refinement_lut(lut, project, istate) 
    allocate(hist_shuffle%raw_counts(lut%sparse_size, 2*size(cv_p%object), 2*size(cv_p%object)), source=0)  
    allocate        (hist%raw_counts(lut%sparse_size, 2*size(cv_p%object), 2*size(cv_p%object)), source=0)  

    ! 2. Echtlauf (Real) - Inkrementelles Stream-Sampling über alle Blöcke
    write(*,*) "=== STARTING HEALPIX  ==="
    hist%raw_counts = 0
    hist_shuffle%raw_counts = 0
    call allocate_marginals(marg_real, lut%npix_fine, 2*size(cv_p%object))
    call allocate_marginals(marg_shuffle, lut%npix_fine, 2*size(cv_p%object))
         
    do b_idx = 1,2!plan%nblocks_total
        write(*, '(A,I4,A,I4)') " -> Processing  Block ", b_idx, " of ", plan%nblocks_total
        write(0,*) 
        
        ! Einzelnen Block laden, liften und rotieren
        call precalculate_block_pair_healpix_trajectory(project, cv_p, plan,  tau, b_idx, block_data,bidx,nside_max)
        
        
        call sample_fine_one_particle_hist(block_data, nside_max, hist_single_fine)
        call transform_trajectory_inplace(block_data, nside_max, nside_hist_sample)
        
        
        
        call aggregate_histogramms(block_data, lut, hist, is_shuffle=.false.)
        call aggregate_histogramms(block_data, lut, hist_shuffle, is_shuffle=.true.)
        
       
       
    end do
    
     filename=''
     filename='entropy_short.txt'
    

     call calculate_marginals(hist%raw_counts, lut%fine_to_sparse, lut%npix_fine, marg_real)
     call calculate_all_entropies(hist%raw_counts, lut%fine_to_sparse, lut, marg_real, entropies)
    
     call save_entropies_gnuplot(trim(filename), entropies)
     !call print_1d_marginals(marg_real, "REAL RUN (TAU=" // trim(tau_str) // ")")
          
     call calculate_te_matrix(hist%raw_counts, lut%fine_to_sparse, lut, marg_real, te_matrix)
     
    do i=1,size(te_matrix,1)
        write(0,*) sum(te_matrix(i,:)),'TE SUM'
   end do
        
        
    open(1,file='test_tematrix1.txt')
    
    do i=1,size(te_matrix,1)
     write(1,*)  te_matrix(i,:)
    end do    
    close(1)
    
     call calculate_marginals(hist_shuffle%raw_counts, lut%fine_to_sparse, lut%npix_fine, marg_shuffle)
    
     !call print_1d_marginals(marg_shuffle, "SHUFFLE RUN (TAU=" // trim(tau_str) // ")")
     filename=''
     filename='entropy_shuffel_short.txt'
     	 
     	 
     	 
     call calculate_te_matrix(hist_shuffle%raw_counts, lut%fine_to_sparse, lut, marg_shuffle, te_matrix_shuffle)
     
     call calculate_all_entropies(hist%raw_counts, lut%fine_to_sparse, lut, marg_shuffle, entropies_shuffel)
     call save_entropies_gnuplot(trim(filename), entropies_shuffel)
       
     do i=1,size(te_matrix_shuffle,1)
        write(0,*) sum(te_matrix_shuffle(i,:)),'TE SUM_shuffle'
        end do
        
        
    open(1,file='test_tematrix_shuffle1.txt')
    
    do i=1,size(te_matrix,1)
     write(1,*)  te_matrix_shuffle(i,:)
    end do
    close(1)
    
    open(1,file='test_tematrix_diff1.txt')
    
    do i=1,size(te_matrix,1)
     write(1,*) te_matrix(i,:)- te_matrix_shuffle(i,:)
    end do
    close(1)
    
! =========================================================================
    ! BEIDE 1D-MARGINALIEN EXPORTIEREN (REAL & SHUFFLE) - OBJEKTGENAU
    ! =========================================================================
    open(10, file='coarse_marginals_real.txt', status='unknown', action='write')
    write(10, '(A)') "# Pixel | Object | X (Future) | Y (Past Tgt) | Z (Past Src)"
    do i = 0, 11
        do j = 1, marg_real%nobject
            write(10, '(2I8, 3I15)') i, j, marg_real%x(i, j), marg_real%y(i, j), marg_real%z(i, j)
        end do
    end do
    close(10)

    open(11, file='coarse_marginals_shuffle.txt', status='unknown', action='write')
    write(11, '(A)') "# Pixel | Object | X (Future) | Y (Past Tgt) | Z (Past Src)"
    do i = 0, 11
        do j = 1, marg_shuffle%nobject
            write(11, '(2I8, 3I15)') i, j, marg_shuffle%x(i, j), marg_shuffle%y(i, j), marg_shuffle%z(i, j)
        end do
    end do
    close(11)
    write(*,*) ">>> Exported both REAL and SHUFFLE object-specific marginals. <<<"
 
  
  filename=''
    filename='hist2'
    call save_fine_one_particle_hist(filename, hist_single_fine, iost)
 end if  
    filename=''
    filename='hist2'
    
   ! call determine_lut_size_from_hists(filename, 192, marg_shuffle%nobject, total_sparse_entries, iost)
    call check_histogram_normalization(filename, 48,2*size(cv_p%object), iost)
   
    write(*,*)  total_sparse_entries,'total_sparse_entries'

    pause
    ! 4. Cleanup
    call deallocate_marginals(marg_real)
    call deallocate_marginals(marg_shuffle)
    
    if (allocated(hist%raw_counts)) deallocate(hist%raw_counts)
    if (allocated(hist_shuffle%raw_counts)) deallocate(hist_shuffle%raw_counts)

end subroutine execute_healpix_analysis


subroutine print_1d_marginals(marg, label)
    type(healpix_marginals), intent(in) :: marg
    character(len=*), intent(in)        :: label
    
    integer :: ipix, i, j
    integer :: sum_x, sum_y, sum_z
    integer :: sum_xy, sum_xz, sum_yz
    logical :: has_data
    logical :: check_failed

    write(*,*) "========================================================="
    write(*,*) " DEBUG OUTPUT & CONSISTENCY CHECK: ", trim(label)
    write(*,*) "========================================================="
    
    has_data = .false.
    check_failed = .false.

    ! =========================================================
    ! TEIL 1: Ausgabe der globalen 1D-Marginalien pro Komponente
    ! =========================================================
    do i = 1, marg%nobject
        sum_x = sum(marg%x(:, i))
        sum_y = sum(marg%y(:, i))
        sum_z = sum(marg%z(:, i))

        ! Wenn für diese Komponente überhaupt Daten vorliegen
        if (sum_x > 0 .or. sum_y > 0 .or. sum_z > 0) then
            has_data = .true.
            write(*, '(A,I4,A,I8,A,I8,A,I8)') &
                "--- Component ID: ", i, " | Total Counts -> X:", sum_x, " | Y:", sum_y, " | Z:", sum_z

            ! Pixelweise Details zeigen
            do ipix = 0, marg%npix-1
                if (marg%x(ipix, i) > 0 .or. marg%y(ipix, i) > 0 .or. marg%z(ipix, i) > 0) then
                    write(*, '(A,I4,A,I8,A,I8,A,I8)') &
                        "    Pixel:", ipix, &
                        " | X (Future):", marg%x(ipix, i), &
                        " | Y (Past Tgt):", marg%y(ipix, i), &
                        " | Z (Past Src):", marg%z(ipix, i)
                end if
            end do
        end if
    end do

    if (.not. has_data) then
        write(*,*) "  [INFO] Alle 1D-Marginalien sind exakt Null. Keine Daten gezählt!"
        write(*,*) "========================================================="
        return
    end if

    ! =========================================================
    ! TEIL 2: Mathematische Konsistenzprüfungen (Kreuzsummen)
    ! =========================================================
    write(*,*) ""
    write(*,*) "--- RUNNING MATHEMATICAL CONSISTENCY CHECKS ---"
    
    do j = 1, marg%nobject
        do i = 1, marg%nobject
            if (i == j) cycle

            ! Hole die Gesamtsummen der 2D-Gelenkwahrscheinlichkeiten für dieses Paar
            sum_xy = sum(marg%xy(:, :, i, j))
            sum_xz = sum(marg%xz(:, :, i, j))
            sum_yz = sum(marg%yz(:, :, i, j))

            ! Check A: Symmetrie der Paarsummen
            ! Jedes gültige Triplett (X,Y,Z) muss in jedem 2D-Histogramm landen.
            ! Die Gesamtsumme von XY, XZ und YZ MUSS für das Paar (i,j) absolut identisch sein.
            if (sum_xy /= sum_xz .or. sum_xy /= sum_yz) then
                write(*, '(A,I4,A,I4,A,I8,A,I8,A,I8)') &
                    "  [CRITICAL ALERT] 2D Sum Mismatch for Pair (", i, ",", j, ")! ", &
                    "Sum(XY):", sum_xy, " | Sum(XZ):", sum_xz, " | Sum(YZ):", sum_yz
                check_failed = .true.
            end if

            ! Check B: Konsistenz zur 1D-Summe der Target-Zukunft X(i)
            ! Die Summe über alle Z-Pixel des XZ-Arrays muss kleiner/gleich der globalen 1D-Summe von X(i) sein
            if (sum_xz > sum(marg%x(:, i))) then
                write(*, '(A,I4,A,I4,A,I8,A,I8)') &
                    "  [CRITICAL ALERT] Pair counts exceed global 1D Target X limits for (", i, ",", j, ")! ", &
                    "Pair Sum(XZ):", sum_xz, " | Global 1D Sum X(i):", sum(marg%x(:, i))
                check_failed = .true.
            end if
            
            ! Check C: Konsistenz zur 1D-Summe der Source-Vergangenheit Z(j)
            ! Die Summe über alle Y-Pixel des YZ-Arrays muss kleiner/gleich der globalen 1D-Summe von Z(j) sein
            if (sum_yz > sum(marg%z(:, j))) then
                write(*, '(A,I4,A,I4,A,I8,A,I8)') &
                    "  [CRITICAL ALERT] Pair counts exceed global 1D Source Z limits for (", i, ",", j, ")! ", &
                    "Pair Sum(YZ):", sum_yz, " | Global 1D Sum Z(j):", sum(marg%z(:, j))
                check_failed = .true.
            end if
        end do
    end do

    if (.not. check_failed) then
        write(*,*) "  [PASSED] All 1D/2D conservation laws hold perfectly."
    else
        write(*,*) "  [FAILED] Conservation law violations detected. Check lookup or bounds!"
    end if
    
    write(*,*) "========================================================="
end subroutine print_1d_marginals

subroutine setup_refinement_lut(lut, project, istate)
    type(refinement_lut), intent(out) :: lut
    type(project_type), intent(in)     :: project
    integer, intent(in)                :: istate

    lut%nside_coarse = 1
    lut%nside_fine   = 1
    lut%npix_coarse  = 12
    lut%npix_fine    = 12

    lut%sparse_size = lut%npix_fine * lut%npix_fine * lut%npix_fine

    ! ALLE Dimensionen müssen von 0 bis Npix-1 laufen!
    allocate(lut%coarse_occupied(0:lut%npix_coarse-1, 0:lut%npix_coarse-1, 0:lut%npix_coarse-1))
    allocate(lut%fine_to_sparse(0:lut%npix_fine-1, 0:lut%npix_fine-1, 0:lut%npix_fine-1))

    lut%coarse_occupied = .true. 

    call rebuild_sparse_mapping(lut)
end subroutine setup_refinement_lut


subroutine allocate_compressed_histogram(hist, lut)
    type(compressed_tau_histogram), intent(out) :: hist
    type(refinement_lut), intent(in)            :: lut

    ! Das Histogramm wird exakt auf die Anzahl der erlaubten/gescannten 
    ! Pixel-Kombinationen (sparse_size) dimensioniert. 
    ! Dimension: (sparse_size, nobject, nobject)
    write(0,*) lut%sparse_size
    !pause
    allocate(hist%raw_counts(lut%sparse_size, lut%npix_coarse, lut%npix_coarse), source=0)
end subroutine allocate_compressed_histogram


! Private Hilfsroutine innerhalb des Moduls, um die IDs zuzuweisen
subroutine rebuild_sparse_mapping(lut)
    type(refinement_lut), intent(inout) :: lut
    integer :: ix, iy, iz, current_id

    lut%fine_to_sparse = 0
    current_id = 0

    ! Schleifen laufen von 0 bis Npix-1
    do iz = 0, lut%npix_fine - 1
        do iy = 0, lut%npix_fine - 1
            do ix = 0, lut%npix_fine - 1
                current_id = current_id + 1
                lut%fine_to_sparse(ix, iy, iz) = current_id
            end do
        end do
    end do
    
    lut%sparse_size = current_id
end subroutine rebuild_sparse_mapping


subroutine precalculate_block_pair_healpix_trajectory(project, cv_p, plan, tau, b_idx, block_data,bidx,nside_max)
    type(project_type), intent(in)               :: project
    type(cv_plan), intent(in)                    :: cv_p
    type(sampling_plan), intent(in)              :: plan
    !type(refinement_lut), intent(in)             :: lut
    integer, intent(in)                          :: nside_max	
    integer, intent(in)                          :: tau
    integer, intent(in)                          :: b_idx  ! <-- Welcher Block berechnet werden soll
    type(healpix_block_trajectory), intent(out)  :: block_data ! <-- Schreibt nur ein einzelnes Element
    
    integer :: target_block_id, shuffle_block_id
    type(block_meta) :: target_meta, shuffle_meta
    type(time_window) :: win_Y, win_Z_shuffle
    logical, allocatable :: mask_Y(:), mask_Z_shuffle(:)
    integer :: idx_past, idx_future
    
    type(raw_data_block) :: raw_Y, raw_Z_real, raw_Z_shuffle
    type(lifted_trajectory) :: data_Y, data_Z_real, data_Z_shuffle
    
    real(8) :: rot_Y(3,3), rot_Z_real(3,3), rot_Z_shuffle(3,3),mat_rel(3,3),R(3,3)
    real(8) :: v_cart(3), v_rot(3)
    integer :: t, p_idx, obj, ipix, n_pairs, n_obj,idx_frame
    integer, intent(in) :: bidx
    logical :: debug
    real(8) :: rot_Z_shuffle_pure(3,3)
    real(8) :: test_spalten(3,3), test_zeilen(3,3)
    type(lifted_trajectory) :: data_Z_shuffle_pure
    real(8) :: rot_Y_future(3,3), rot_Z_shuffle_future(3,3)
    real(8) :: mat_rel_past(3,3), mat_rel_future(3,3)
    
    !debug=.true.
    debug=.false.
     call generate_random_rotation_matrix(R)
    
    !für debug:
    
    n_obj = 2*size(cv_p%object)
    
    ! Block-IDs laut Sampling-Plan holen
    target_block_id  = plan%real_run_perm(b_idx)
    shuffle_block_id = plan%shuffle_run_perm(b_idx)
    
    target_meta  = plan%global_blocks(target_block_id)
    shuffle_meta = plan%global_blocks(shuffle_block_id)
    
    allocate(mask_Y(target_meta%length))
    allocate(mask_Z_shuffle(shuffle_meta%length))
    mask_Y = .false.
    mask_Z_shuffle = .false.
    
    ! Maskierung vorbereiten (wie gehabt)
    do t = 1, target_meta%length, tau
        if (t + tau > target_meta%length) exit
        mask_Y(t)       = .true.
        mask_Y(t + tau) = .true.
    end do
    
    do t = 1, shuffle_meta%length, tau
        if (t + tau > shuffle_meta%length) exit
        mask_Z_shuffle(t)       = .true.
        mask_Z_shuffle(t + tau) = .true.
    end do
    
  
    ! Zeitfenster zuweisen
    win_Y%istate          = plan%istate
    win_Y%irep            = target_meta%replica_id
    win_Y%start_time_step = target_meta%start_frame
    win_Y%end_time_step   = target_meta%start_frame + target_meta%length - 1
    
    win_Z_shuffle%istate          = plan%istate
    win_Z_shuffle%irep            = shuffle_meta%replica_id
    win_Z_shuffle%start_time_step = shuffle_meta%start_frame
    win_Z_shuffle%end_time_step   = shuffle_meta%start_frame + shuffle_meta%length - 1
    
    ! Rohdaten laden und transformieren
! Rohdaten laden und transformieren
    call load_raw_data_block(project, win_Y, mask_Y, raw_Y)
    call lift_cv_plan(cv_p, raw_Y, data_Y)
    
     do idx_frame = 1, size(data_Y%vec, 3)
        do obj = 1, size(data_Y%vec, 2)
            
            ! Drehe jeden einzelnen kartesischen Richtungsvektor global im Raum
            data_Y%vec(:, obj, idx_frame) =  data_Y%vec(:, obj, idx_frame)&
            /norm2( data_Y%vec(:, obj, idx_frame))
            
        end do
    end do
    
   
    ! Direkt kopieren statt die Datei doppelt zu lesen!
    data_Z_real = data_Y 
    
    
    ! Shuffle muss natürlich separat geladen werden, da anderes Fenster
    call load_raw_data_block(project, win_Z_shuffle, mask_Z_shuffle, raw_Z_shuffle)
    call lift_cv_plan(cv_p, raw_Z_shuffle, data_Z_shuffle)
    
      do idx_frame = 1, size(data_Z_shuffle%vec, 3)
        do obj = 1, size(data_Z_shuffle%vec, 2)
            
            ! nomrmieren
            data_Z_shuffle%vec(:, obj, idx_frame) =  data_Z_shuffle%vec(:, obj, idx_frame)&
            /norm2( data_Z_shuffle%vec(:, obj, idx_frame))
            
        end do
    end do
    
   
    data_Z_shuffle_pure = data_Z_shuffle
    ! --- DEBUG TEST: DIE GANZE TRAJEKTORIE IM RAUM VERDREHEN ---
   if(debug) then
    ! Schleife über alle Objekte (2) und alle Frames (3) der Trajektorie
    do idx_frame = 1, size(data_Z_shuffle%vec, 3)
        do obj = 1, size(data_Z_shuffle%vec, 2)
            
            ! Drehe jeden einzelnen kartesischen Richtungsvektor global im Raum
            data_Z_shuffle%vec(:, obj, idx_frame) = matmul(R, data_Z_shuffle%vec(:, obj, idx_frame))
            
        end do
    end do
    R=0._8
    end if
    ! -----------------------------------------------------------
   
    ! Krisenfeste Bestimmung der Paare anhand der echt geladenen Vektoren
    n_pairs = size(data_Y%vec, 3)-1     
    block_data%nobject = n_obj
    block_data%npairs  = n_pairs
    
    write(0,*) n_obj,'n_obj'
    !
    
    if (n_pairs < 1) then
        write(0,*) "[WARNUNG] Block ", b_idx, " hat keine validen Paare geladen!"
        deallocate(mask_Y, mask_Z_shuffle)
        return
    end if
    
    ! Lokale Arrays im Einzelblock allokieren
    allocate(block_data%Y_past          (n_obj, n_pairs))
    allocate(block_data%Y_future        (n_obj, n_pairs))
    allocate(block_data%Z_past_real     (n_obj, n_pairs))
    allocate(block_data%Z_past_shuffle  (n_obj, n_pairs))
    allocate(block_data%Z_future_shuffle(n_obj, n_pairs))
     block_data%Y_past =0
     block_data%Y_future=0
     block_data%Z_past_real=0
     block_data%Z_past_shuffle=0
     block_data%Z_future_shuffle=0

    ! HEALPix-Projektionsschleife
    write(0,*) n_pairs
	!pause
 do p_idx = 1, n_pairs
    
        idx_past   = p_idx
        idx_future = p_idx + 1
        
        ! Basen für die Vergangenheit (Zeitpunkt t und t') extrahieren
        call extract_backbone_bases(bidx, data_Y, idx_past, rot_Y)
        call extract_backbone_bases(bidx, data_Z_shuffle, idx_past, rot_Z_shuffle)
        
        ! Die einzig konsistente relative Transformation berechnen
        mat_rel =rot_Z_shuffle! matmul(transpose(rot_Y), rot_Z_shuffle)

        !mat_rel = matmul(rot_Y, transpose(rot_Z_shuffle))
        do obj = 1, n_obj
            ! --- Y PAST ---
            v_cart = data_Y%vec(:, obj, idx_past)
            v_rot  = matmul(rot_Y, v_cart)
            call vec2pix_nest(nside_max, v_rot, ipix)
            block_data%Y_past(obj, p_idx) = ipix
            
            ! --- Y FUTURE ---
            v_cart = data_Y%vec(:, obj, idx_future)
            v_rot  = matmul(rot_Y, v_cart) ! Nutzt starr die Basis der Vergangenheit
            call vec2pix_nest(nside_max, v_rot, ipix)
            block_data%Y_future(obj, p_idx) = ipix
            
            ! --- Z PAST REAL ---
            v_cart = data_Z_real%vec(:, obj, idx_past)
            v_rot  = matmul(rot_Y, v_cart) ! Identisch zu Y_past, da data_Z_real == data_Y
            call vec2pix_nest(nside_max, v_rot, ipix)
            block_data%Z_past_real(obj, p_idx) = ipix
            
            ! --- Z PAST SHUFFLE ---
            v_cart = data_Z_shuffle%vec(:, obj, idx_past)
            v_rot  = matmul(mat_rel, v_cart) 
            call vec2pix_nest(nside_max, v_rot, ipix)
            block_data%Z_past_shuffle(obj, p_idx) = ipix
            
            ! --- Z FUTURE SHUFFLE ---
            v_cart = data_Z_shuffle%vec(:, obj, idx_future)
            v_rot  = matmul(mat_rel, v_cart) ! Nutzt dieselbe relative Transformation
            call vec2pix_nest(nside_max, v_rot, ipix)
            block_data%Z_future_shuffle(obj, p_idx) = ipix
        end do
        
    end do
    
    deallocate(mask_Y, mask_Z_shuffle)
end subroutine precalculate_block_pair_healpix_trajectory



subroutine aggregate_histogramms(block_data, lut, hist, is_shuffle)
    type(healpix_block_trajectory), intent(in)   :: block_data 
    type(refinement_lut), intent(in)             :: lut
    type(compressed_tau_histogram), intent(inout) :: hist
    logical, intent(in)                          :: is_shuffle
    
    integer :: p_idx, i, j
    integer :: ipix_x, ipix_y, ipix_z, idx_sparse
    
   if (is_shuffle) then
     	 
    do p_idx = 1, block_data%npairs
        do i = 1, block_data%nobject
            ipix_y = block_data%Y_past(i, p_idx)
            ipix_x = block_data%Y_future(i, p_idx)
            
           
                do j = 1, block_data%nobject
                    ipix_z = block_data%Z_past_shuffle(j, p_idx)
                    idx_sparse = fine_to_sparse_lookup(lut, ipix_x, ipix_y, ipix_z)
                    !if (idx_sparse > 0) then
                        hist%raw_counts(idx_sparse, i, j) = hist%raw_counts(idx_sparse, i, j) + 1
                    !end if
                end do
  
        end do
        
     
    end do
    
    else
        do p_idx = 1, block_data%npairs
        do i = 1, block_data%nobject
            ipix_y = block_data%Y_past(i, p_idx)
            ipix_x = block_data%Y_future(i, p_idx)
            
       
                do j = 1, block_data%nobject
                    ipix_z = block_data%Z_past_real(j, p_idx)
                    idx_sparse = fine_to_sparse_lookup(lut, ipix_x, ipix_y, ipix_z)
                    !if (idx_sparse > 0) then
                        hist%raw_counts(idx_sparse, i, j) = hist%raw_counts(idx_sparse, i, j) + 1
                    !end if
                end do
            
        end do
        
     
    end do
    
    end if
end subroutine aggregate_histogramms




! Hilfsfunktion für den LUT Lookup 
pure integer function fine_to_sparse_lookup(lut, ix, iy, iz)
    type(refinement_lut), intent(in) :: lut
    integer, intent(in) :: ix, iy, iz
    fine_to_sparse_lookup = lut%fine_to_sparse(ix, iy, iz)
end function fine_to_sparse_lookup
!
subroutine extract_backbone_bases(bidx,data_obj, t, rot)
    type(lifted_trajectory), intent(in) :: data_obj
    integer, intent(in)                 :: t
    real(8), intent(out)                :: rot(3,3)
    real(8) :: one(3,3)
    real(8) :: v1(3),v2(3),v3(3)
    integer :: bidx
    one=0._8
    one(1,1)=1._8;one(2,2)=1._8;one(3,3)=1._8
    !write(0,*) bidx
    !pause
    ! Wir extrahieren die 3x3 Matrix für den Zeitschritt t
    ! Falls deine rot_matrix anders indiziert ist (z.B. flach), 
    ! passen wir das hier an. Hier gehen wir von (3,3,ntimes) aus:
    v1=data_obj%vec(:, 2*bidx-1, t)
    v2=data_obj%vec(:, 2*bidx  , t)
    v2=v2-dot_product(v1,v2)*v1
    v3=cross_product(v1,v2)
    
    rot(1,:) = v1
    rot(2,:) = v2
    rot(3,:) = v3
 
    !rot=transpose(rot)
end subroutine extract_backbone_bases

function cross_product(x1,x2) result(c)

implicit none

real(8), intent(in) :: x1(3), x2(3)
real(8) :: c(3)

c(1) = x1(2)*x2(3) - x1(3)*x2(2)
c(2) = x1(3)*x2(1) - x1(1)*x2(3)
c(3) = x1(1)*x2(2) - x1(2)*x2(1)

end function cross_product

subroutine save_fine_one_particle_hist(base_filename, hist, stat)
    implicit none
    character(len=*), intent(in)               :: base_filename
    type(fine_one_particle_hist), intent(in)   :: hist
    integer, intent(out), optional             :: stat

    integer :: u, ios, i_bin, idx
    character(len=256) :: full_name
    character(len=*), parameter :: suff(5) = [ &
        character(len=24) :: "_Y_past.txt", "_Y_future.txt", &
                             "_Z_real.txt", "_Z_shuffle.txt", &
                             "_Z_future_shuffle.txt" ]

    if (present(stat)) stat = 0

    do i_bin = 1, 5
        select case(i_bin)
        case(1); full_name = trim(base_filename) // suff(1)
        case(2); full_name = trim(base_filename) // suff(2)
        case(3); full_name = trim(base_filename) // suff(3)
        case(4); full_name = trim(base_filename) // suff(4)
        case(5); full_name = trim(base_filename) // suff(5)
        end select

        ! Nur das Öffnen wird noch überwacht, falls der Pfad nicht existiert
        open(newunit=u, file=trim(full_name), form='formatted', status='replace', iostat=ios)
        if (ios /= 0) then
            if (present(stat)) stat = ios
            return
        end if

        ! Reines, ungestörtes Schreiben der Matrizen (Zeile = Bins, Spalte = Objekte)
        select case(i_bin)
        case(1)
            do idx = 0, hist%npix - 1
                write(u, *) hist%count_Y_past(idx, :)
            end do
        case(2)
            do idx = 0, hist%npix - 1
                write(u, *) hist%count_Y_future(idx, :)
            end do
        case(3)
            do idx = 0, hist%npix - 1
                write(u, *) hist%count_Z_real(idx, :)
            end do
        case(4)
            do idx = 0, hist%npix - 1
                write(u, *) hist%count_Z_shuffle(idx, :)
            end do
        case(5)
            do idx = 0, hist%npix - 1
                write(u, *) hist%count_Z_future_shuffle(idx, :)
            end do
        end select

        close(u)
    end do
end subroutine save_fine_one_particle_hist

subroutine load_fine_one_particle_hist(base_filename, npix, nobject, hist, stat)
    implicit none
    character(len=*), intent(in)               :: base_filename
    integer, intent(in)                        :: npix
    integer, intent(in)                        :: nobject
    type(fine_one_particle_hist), intent(out)  :: hist
    integer, intent(out), optional             :: stat

    integer :: u, ios, i_bin, idx
    character(len=256) :: full_name
    character(len=*), parameter :: suff(5) = [ &
        character(len=24) :: "_Y_past.txt", "_Y_future.txt", &
                             "_Z_real.txt", "_Z_shuffle.txt", &
                             "_Z_future_shuffle.txt" ]

    if (present(stat)) stat = 0
    
    hist%npix = npix
    hist%nobject = nobject

    ! Allokierung prüfen wir mit stat ab
    allocate(hist%count_Y_past(0:hist%npix-1, hist%nobject), stat=ios)
    if (ios == 0) allocate(hist%count_Y_future(0:hist%npix-1, hist%nobject), stat=ios)
    if (ios == 0) allocate(hist%count_Z_real(0:hist%npix-1, hist%nobject), stat=ios)
    if (ios == 0) allocate(hist%count_Z_shuffle(0:hist%npix-1, hist%nobject), stat=ios)
    if (ios == 0) allocate(hist%count_Z_future_shuffle(0:hist%npix-1, hist%nobject), stat=ios)

    if (ios /= 0) then
        if (present(stat)) stat = ios
        return
    end if

    do i_bin = 1, 5
        select case(i_bin)
        case(1); full_name = trim(base_filename) // suff(1)
        case(2); full_name = trim(base_filename) // suff(2)
        case(3); full_name = trim(base_filename) // suff(3)
        case(4); full_name = trim(base_filename) // suff(4)
        case(5); full_name = trim(base_filename) // suff(5)
        end select

        open(newunit=u, file=trim(full_name), form='formatted', status='old', iostat=ios)
        if (ios /= 0) then
            if (present(stat)) stat = ios
            return
        end if

        select case(i_bin)
        case(1)
            do idx = 0, hist%npix - 1
                read(u, *) hist%count_Y_past(idx, :)
            end do
        case(2)
            do idx = 0, hist%npix - 1
                read(u, *) hist%count_Y_future(idx, :)
            end do
        case(3)
            do idx = 0, hist%npix - 1
                read(u, *) hist%count_Z_real(idx, :)
            end do
        case(4)
            do idx = 0, hist%npix - 1
                read(u, *) hist%count_Z_shuffle(idx, :)
            end do
        case(5)
            do idx = 0, hist%npix - 1
                read(u, *) hist%count_Z_future_shuffle(idx, :)
            end do
        end select

        close(u)
    end do
end subroutine load_fine_one_particle_hist


subroutine check_histogram_normalization(base_filename, npix, nobject, stat)
    implicit none
    character(len=*), intent(in)               :: base_filename
    integer, intent(in)                        :: npix
    integer, intent(in)                        :: nobject
    integer, intent(out), optional             :: stat

    type(fine_one_particle_hist) :: hist
    integer :: ios, i_obj, j_obj, i
    integer(8) :: sum_Y_past, sum_Y_future, sum_Z_real
    integer(8) :: sum_Z_shuffle, sum_Z_future_shuffle
    
    ! Counter für die aktiven Bins pro Objekt
    integer, allocatable :: active_Y_past_count(:)
    integer, allocatable :: active_Y_future_count(:)
    integer, allocatable :: active_Z_any_count(:)
    integer(8) :: max_theoretical_3d_bins,prod

    sum_Y_past = 0
    sum_Y_future = 0
    sum_Z_real = 0
    sum_Z_shuffle = 0
    sum_Z_future_shuffle = 0

    ! 1. Histogramme über die bereinigte Lade-Routine einlesen
    if (present(stat)) then
        call load_fine_one_particle_hist(base_filename, npix, nobject, hist, stat=ios)
        if (ios /= 0) then
            stat = ios
            return
        end if
    else
        call load_fine_one_particle_hist(base_filename, npix, nobject, hist)
    end if

    ! Hilfsarrays allokieren, um die aktiven Bins pro Objekt zu zählen
    allocate(active_Y_past_count(nobject), source=0)
    allocate(active_Y_future_count(nobject), source=0)
    allocate(active_Z_any_count(nobject), source=0)

    ! 2. Über alle Objekte loopen und aktive Bins (Counts > 0) zählen
    write(0,*) "Objekte insgesamt: ", nobject     
    
    do i_obj = 1, nobject
        do i = 0, npix-1
            if (hist%count_Y_past(i, i_obj) > 0) then
                active_Y_past_count(i_obj) = active_Y_past_count(i_obj) + 1
            end if
            if (hist%count_Y_future(i, i_obj) > 0) then
                active_Y_future_count(i_obj) = active_Y_future_count(i_obj) + 1
            end if
            ! Für Z nehmen wir an, dass es entweder im Real- ODER Shuffle-Modus auftaucht
            if (hist%count_Z_real(i, i_obj) > 0 .or. hist%count_Z_shuffle(i, i_obj) > 0) then
                active_Z_any_count(i_obj) = active_Z_any_count(i_obj) + 1
            end if
        end do
        
        write(0,*) "Objekt ", i_obj, " Aktive Bins (Past/Future/Z):", &
                   active_Y_past_count(i_obj), active_Y_future_count(i_obj), active_Z_any_count(i_obj)
        
        ! Klassische Summenbildung der Counts
        sum_Y_past           = sum_Y_past           + sum(int(hist%count_Y_past(:, i_obj), 8))
        sum_Y_future         = sum_Y_future         + sum(int(hist%count_Y_future(:, i_obj), 8))
        sum_Z_real           = sum_Z_real           + sum(int(hist%count_Z_real(:, i_obj), 8))
        sum_Z_shuffle        = sum_Z_shuffle        + sum(int(hist%count_Z_shuffle(:, i_obj), 8))
        sum_Z_future_shuffle = sum_Z_future_shuffle + sum(int(hist%count_Z_future_shuffle(:, i_obj), 8))
    end do

    ! 3. Berechnen, wie viele 3D-Kombinationen die Paarschleifen (i, j) maximal erzeugen können
    max_theoretical_3d_bins = 0
    do i_obj = 1, nobject
       
            ! Kombinationen für das Paar (i_obj, j_obj) = Bins_Y_past(i) * Bins_Y_future(i) * Bins_Z(j)
           
            prod= &
                (int(active_Y_past_count(i_obj), 8) * &
                 int(active_Y_future_count(i_obj), 8) * &
                 int(active_Z_any_count(i_obj), 8))
             
               
            max_theoretical_3d_bins = max_theoretical_3d_bins+prod
             write(0,*) prod,max_theoretical_3d_bins

    end do

    ! 4. Ergebnisse übersichtlich im Terminal ausgeben
    print *, "=== Normierungs-Check für: ", trim(base_filename), " ==="
    print *, "Gesamtsumme Y_past:           ", sum_Y_past
    print *, "Gesamtsumme Y_future:         ", sum_Y_future
    print *, "Gesamtsumme Z_real:           ", sum_Z_real
    print *, "Gesamtsumme Z_shuffle:        ", sum_Z_shuffle
    print *, "Gesamtsumme Z_future_shuffle: ", sum_Z_future_shuffle
    print *, "--------------------------------------------------------"
    print *, "Theoretisch maximale 3D-Bins aus 1D-Dichten (aufsummiert über alle Paare):", max_theoretical_3d_bins
    print *, "========================================================="
     
    ! Speicher aufräumen
    deallocate(active_Y_past_count, active_Y_future_count, active_Z_any_count)
    deallocate(hist%count_Y_past)
    deallocate(hist%count_Y_future)
    deallocate(hist%count_Z_real)
    deallocate(hist%count_Z_shuffle)
    deallocate(hist%count_Z_future_shuffle)

    if (present(stat)) stat = 0
end subroutine check_histogram_normalization

! --- HILFSROUTINE FÜR DEN DEBUG-TEST ---
subroutine generate_random_rotation_matrix(R)
    real(8), intent(out) :: R(3,3)
    real(8) :: r_ann(3)
    real(8) :: theta, phi, psi
    
    ! Erzeuge drei Pseudozufallswinkel
    call random_number(r_ann)
    theta = r_ann(1) * 3.141592653589793d0  ! 0 bis pi
    phi   = r_ann(2) * 2.d0 * 3.141592653589793d0 ! 0 bis 2pi
    psi   = r_ann(3) * 2.d0 * 3.141592653589793d0 ! 0 bis 2pi
    
    ! Eine einfache, gültige Rotationsmatrix (z.B. Euler-Winkel-Kombination)
    ! Hauptsache sie ist orthogonal (R^T * R = I)
    R(1,1) = cos(phi)*cos(psi) - sin(phi)*cos(theta)*sin(psi)
    R(1,2) = -cos(phi)*sin(psi) - sin(phi)*cos(theta)*cos(psi)
    R(1,3) = sin(phi)*sin(theta)
    
    R(2,1) = sin(phi)*cos(psi) + cos(phi)*cos(theta)*sin(psi)
    R(2,2) = -sin(phi)*sin(psi) + cos(phi)*cos(theta)*cos(psi)
    R(2,3) = -cos(phi)*sin(theta)
    
    R(3,1) = sin(theta)*sin(psi)
    R(3,2) = sin(theta)*cos(psi)
    R(3,3) = cos(theta)
end subroutine generate_random_rotation_matrix


end module healpix_analysis_execution_module_single_lag