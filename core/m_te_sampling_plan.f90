module healpix_histogram_sampling_plan_module
    use iso_fortran_env, only: real32, int32, int64
    use lifting_plan_module
    use position_lift_module
    use sparse_lut_module
    use pix_tools
    use lifted_trajectory_io_module
    use lifting_plan_execution_module
    use healpix_state_space_exploration, only: build_smoothed_backbone_basis,extract_backbone_basis
    use mod_te_io_reconstruction
    !$ use omp_lib
    implicit none
    private

    public :: sampling_plan_type, lifted_tra_int,execute_shuffle_sampling_plan
    public :: lift_trajectories,make_shuffle_te_chunk_filename
    public :: create_sampling_plan
    public :: execute_sampling_plan
    public :: free_sampling_plan

    ! ===================================================================
    ! DATENTYPEN
    ! ===================================================================

    ! Kompakte Geliftete Trajektorie im RAM (nur Integers!)
    type :: lifted_tra_int
        integer(int32), allocatable :: pix_present(:,:) ! (nobj, nframe)
        integer(int32), allocatable :: pix_lagged(:,:)  ! (nobj, nframe)
        integer :: nframe
        integer :: nobj
    end type lifted_tra_int
    
    

    ! Eine Gruppe von Chunks, die gemeinsam ins RAM-Budget passen
    type :: chunk_group_type
        integer :: n_chunks_in_group
        integer, allocatable :: chunk_ids(:)       ! globale Chunk-IDs in dieser Gruppe
        integer(int64) :: total_elements            ! Gesamtzahl Elemente (32-Bit Reals) in Gruppe
        integer(int64), allocatable :: chunk_offsets(:) ! Offset des Chunks im Gruppen-RAM-Puffer
        integer :: n_active_pairs
        integer, allocatable :: active_pair_indices(:)
    end type chunk_group_type

    ! Der Gesamte Ausführungsplan
    type :: sampling_plan_type
        real(8) :: max_ram_gb                       ! Erlaubtes RAM-Budget
        integer :: n_groups                         ! Anzahl der benötigten Durchläufe (Passes)
        type(chunk_group_type), allocatable :: groups(:)
        
        ! Fast-Lookup: Map von globaler Chunk-ID auf (group_id, local_offset)
        integer, allocatable :: chunk2group(:)
        integer(int64), allocatable :: chunk2offset(:)
    end type sampling_plan_type

contains

    ! ===================================================================
    ! PHASE 1: HEALPix Lifting (Einmaliges Pre-Processing auf Platte)
    ! ===================================================================
    subroutine lift_trajectories(lifting_plan, raw_prefix, lifted_prefix, npix, lag)
        ! Standard-Imports für deine Geometrie-Routinen hier vorausgesetzt
        type(lifting_plan_type), intent(in) :: lifting_plan
        character(len=*),        intent(in) :: raw_prefix
        character(len=*),        intent(in) :: lifted_prefix
        integer,                 intent(in) :: npix
        integer,                 intent(in) :: lag

        character(len=800) :: raw_file, lifted_file
        type(lifted_trajectory) :: raw_tra ! Deine originale 3D-Trajektorie
        integer :: istate, irep, iblock, i, t, tlag, nframe, nobj, nside, basis_index
        integer :: i1, i2, u_out
        integer(int32), allocatable :: pix_pres(:,:), pix_lag(:,:)
        real(8) :: rot_i(3,3), inv_n1, inv_n2
        real(8) :: vec_pres(3), vec_lag(3), loc_pres(3), loc_lag(3)
        real(8), allocatable :: rots(:,:,:)
        real(8):: dt,tau_mem
        nside = npix_to_nside(npix)
        basis_index = 6

        write(0,*) '-----------------------------------------------------'
        write(0,*) 'PHASE 1: Lifting Trajectories into HEALPix Space'
        write(0,*) '-----------------------------------------------------'
    ! ------------------------------------------------------------
    ! Parameters for the smoothed basis
    ! ------------------------------------------------------------

    tau_mem = 100.0_8       ! Beispielwert
    dt      = 1.0_8         ! Zeit zwischen zwei Frames


        do istate = 1, size(lifting_plan%istates)
            do irep = 1, lifting_plan%nreps(istate)
                do iblock = 1, lifting_plan%nblock(istate, irep)

                    call make_filename_lifted_tra(istate, irep, iblock, raw_prefix, raw_file)
                    call make_filename_lifted_tra(istate, irep, iblock, lifted_prefix, lifted_file)

                    write(0,*) 'Lifting: ', trim(raw_file), ' -> ', trim(lifted_file)
                    call read_trajectory(raw_file, raw_tra)
                    
                    nobj = size(raw_tra%vec, 2)
                    nframe = size(raw_tra%vec, 3)

                    allocate(pix_pres(nobj, nframe))
                    allocate(pix_lag(nobj, nframe))
                    pix_pres = 0
                    pix_lag  = 0
                    
                    allocate(rots(3,3,nframe))

                    call build_smoothed_backbone_basis( &
                    raw_tra, basis_index, tau_mem, dt, rots)

                    ! $OMP PARALLEL DO SCHEDULE(static) &
                    ! $OMP& PRIVATE(i, t, tlag, rot_i, inv_n1, inv_n2, &
                    ! $OMP&         vec_pres, vec_lag, loc_pres, loc_lag, i1, i2)
                    do i = 1, nobj
                        
                        do t = 1, nframe
                        rot_i = rots(:,:,t)
                            ! Present
                            inv_n2 = 1.0_8 / norm2(raw_tra%vec(:, i, t))
                            vec_pres = raw_tra%vec(:, i, t) * inv_n2

                            !call extract_backbone_basis(basis_index, raw_tra, t, rot_i)

                            loc_pres(1) = rot_i(1,1)*vec_pres(1) + rot_i(1,2)*vec_pres(2) + rot_i(1,3)*vec_pres(3)
                            loc_pres(2) = rot_i(2,1)*vec_pres(1) + rot_i(2,2)*vec_pres(2) + rot_i(2,3)*vec_pres(3)
                            loc_pres(3) = rot_i(3,1)*vec_pres(1) + rot_i(3,2)*vec_pres(2) + rot_i(3,3)*vec_pres(3)

                            call vec2pix_nest(nside, loc_pres, i2)
                            pix_pres(i, t) = i2

                            ! Lagged
                            tlag = t + lag
                            if (tlag <= nframe) then
                                inv_n1 = 1.0_8 / norm2(raw_tra%vec(:, i, tlag))
                                vec_lag = raw_tra%vec(:, i, tlag) * inv_n1

                                loc_lag(1) = rot_i(1,1)*vec_lag(1) + rot_i(1,2)*vec_lag(2) + rot_i(1,3)*vec_lag(3)
                                loc_lag(2) = rot_i(2,1)*vec_lag(1) + rot_i(2,2)*vec_lag(2) + rot_i(2,3)*vec_lag(3)
                                loc_lag(3) = rot_i(3,1)*vec_lag(1) + rot_i(3,2)*vec_lag(2) + rot_i(3,3)*vec_lag(3)

                                call vec2pix_nest(nside, loc_lag, i1)
                                pix_lag(i, t) = i1
                            end if
                        end do
                    end do
                    ! $OMP END PARALLEL DO

                    ! Geliftete Integer-Daten binaer auf Platte schreiben
                    open(newunit=u_out, file=trim(lifted_file), access='stream', &
                         form='unformatted', status='replace')
                    write(u_out) nobj, nframe
                    write(u_out) pix_pres
                    write(u_out) pix_lag
                    close(u_out)

                    deallocate(pix_pres, pix_lag,rots)
                end do
            end do
        end do
    end subroutine lift_trajectories

    pure function npix_to_nside(npix) result(nside)
    implicit none
    integer, intent(in) :: npix
    integer             :: nside
    
    nside = int(sqrt(real(npix, 8) / 12.0_8))
end function npix_to_nside


    ! ===================================================================
    ! PHASE 2: Planning Routine (Erstellt Gruppen basierend auf RAM-Limit)
    ! ===================================================================
    subroutine create_sampling_plan(plan, pair_map, n_chunks, max_ram_gb)
        type(sampling_plan_type), intent(out) :: plan
        type(pair_map_type),      intent(in)  :: pair_map
        integer,                  intent(in)  :: n_chunks
        real(8),                  intent(in)  :: max_ram_gb

        integer(int64) :: max_bytes_per_group, chunk_size_bytes, current_group_bytes
        integer(int64), allocatable :: chunk_sizes(:)
        integer :: ichunk, active_idx, i, j, igroup, curr_count
        integer, allocatable :: chunk_group_assignment(:)

        plan%max_ram_gb = max_ram_gb
        max_bytes_per_group = int(max_ram_gb * 1024.0_8**3, int64)

        allocate(chunk_sizes(n_chunks))
        allocate(chunk_group_assignment(n_chunks))
        allocate(plan%chunk2group(n_chunks))
        allocate(plan%chunk2offset(n_chunks))
        chunk_sizes = 0_int64

        ! 1. Exakte Elementanzahl pro Chunk aus PairMap ermitteln
        do active_idx = 1, pair_map%n_active
            i = pair_map%active_i(active_idx)
            j = pair_map%active_j(active_idx)
            ichunk = pair_map%chunk(i, j)
            if (ichunk >= 1 .and. ichunk <= n_chunks) then
                chunk_sizes(ichunk) = max(chunk_sizes(ichunk), &
                    int(pair_map%local_start(i, j), int64) + int(pair_map%size(i, j), int64) - 1_int64)
            end if
        end do

        ! 2. Chunks in Gruppen gießen
        igroup = 1
        current_group_bytes = 0_int64

        do ichunk = 1, n_chunks
            chunk_size_bytes = chunk_sizes(ichunk) * 4_int64 ! 4 Bytes pro real32

            if (chunk_size_bytes == 0_int64) cycle

            if (current_group_bytes + chunk_size_bytes > max_bytes_per_group .and. current_group_bytes > 0) then
                igroup = igroup + 1
                current_group_bytes = 0_int64
            end if

            chunk_group_assignment(ichunk) = igroup
            current_group_bytes = current_group_bytes + chunk_size_bytes
        end do

        plan%n_groups = igroup
        allocate(plan%groups(plan%n_groups))

        write(0,*) '-----------------------------------------------------'
        write(0,*) 'PHASE 2: Sampling Plan Generated'
        write(0,*) 'Total Chunks:', n_chunks
        write(0,*) 'Total Execution Passes (Groups):', plan%n_groups
        write(0,*) 'RAM Budget per Pass:', max_ram_gb, 'GB'
        write(0,*) '-----------------------------------------------------'

        ! 3. Gruppendatenstrukturen befüllen
        do igroup = 1, plan%n_groups
            curr_count = count(chunk_group_assignment == igroup)
            plan%groups(igroup)%n_chunks_in_group = curr_count
            allocate(plan%groups(igroup)%chunk_ids(curr_count))
            allocate(plan%groups(igroup)%chunk_offsets(curr_count))

            curr_count = 0
            current_group_bytes = 0_int64

            do ichunk = 1, n_chunks
                if (chunk_group_assignment(ichunk) == igroup) then
                    curr_count = curr_count + 1
                    plan%groups(igroup)%chunk_ids(curr_count) = ichunk
                    plan%groups(igroup)%chunk_offsets(curr_count) = current_group_bytes / 4_int64 ! Element-Offset

                    plan%chunk2group(ichunk) = igroup
                    plan%chunk2offset(ichunk) = current_group_bytes / 4_int64

                    current_group_bytes = current_group_bytes + chunk_sizes(ichunk) * 4_int64
                end if
            end do
            plan%groups(igroup)%total_elements = current_group_bytes / 4_int64
                ! Aktive Paare fuer diese Gruppe einsammeln (Vorgeschaltet!)
        curr_count = 0
        do active_idx = 1, pair_map%n_active
            i = pair_map%active_i(active_idx)
            j = pair_map%active_j(active_idx)
            ichunk = pair_map%chunk(i, j)
            if (chunk_group_assignment(ichunk) == igroup) then
                curr_count = curr_count + 1
            end if
        end do
        
        plan%groups(igroup)%n_active_pairs = curr_count
        allocate(plan%groups(igroup)%active_pair_indices(curr_count))
        
        curr_count = 0
        do active_idx = 1, pair_map%n_active
            i = pair_map%active_i(active_idx)
            j = pair_map%active_j(active_idx)
            ichunk = pair_map%chunk(i, j)
            if (chunk_group_assignment(ichunk) == igroup) then
                curr_count = curr_count + 1
                plan%groups(igroup)%active_pair_indices(curr_count) = active_idx
            end if
        end do

            
            
            
        end do
        

        deallocate(chunk_sizes, chunk_group_assignment)
    end subroutine create_sampling_plan


    ! ===================================================================
    ! PHASE 3: Sampling Plan Execution (Gruppenweise Ausführung)
    ! ===================================================================
! ===================================================================
subroutine execute_sampling_plan(plan, lifting_plan, hits, pair_map, &
                                 npix, lag, lifted_prefix, output_prefix)
    use iso_fortran_env, only: real32, int32, int64
    !$ use omp_lib
    implicit none

    type(sampling_plan_type),intent(in) :: plan
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(hitlist_type),      intent(in) :: hits
    type(pair_map_type),     intent(in) :: pair_map
    integer,                 intent(in) :: npix
    integer,                 intent(in) :: lag
    character(len=*),        intent(in) :: lifted_prefix
    character(len=*),        intent(in) :: output_prefix

    type(lifted_tra_int) :: tra
    character(len=800)   :: lifted_file, chunk_out_file
    real(int32), allocatable :: group_hist(:)

    integer :: igroup, istate, irep, iblock, pair_entry, active_idx, ichunk
    integer :: i, j, pair_size, u_in, u_out, idx_in_grp
    integer :: t, i1, i2, i3, nframe, nobj, nframes_to_sample
    integer :: k1, k2, k3
    
    ! Vorberechnete Multiplikatoren für das Paar (i, j)
    integer(int64) :: stride_k2, stride_k3
    integer(int64) :: global_idx, group_hist_idx, group_offset, local_start_idx
    integer(int64) :: base_pair_offset
    integer(int64) :: n_evaluations_rep1, n_hits_rep1, total_reps_count

    ! Zeitmessung
    real(8) :: t_pass_start, t_start_section
    real(8) :: t_rep1_read, t_rep1_sample, t_rep1_write
    logical :: is_rep1

    write(0,*) '-----------------------------------------------------'
    write(0,*) 'PHASE 3: Executing Sampling Plan (Parallelized)'
    write(0,*) '-----------------------------------------------------'

    ! Outer Loop: Chunk-Gruppen (Geringer RAM-Verbrauch)
    do igroup = 1, plan%n_groups
    !if(igroup<7) cycle
        t_pass_start = 0.0_8
        t_rep1_read   = 0.0_8
        t_rep1_sample = 0.0_8
        t_rep1_write  = 0.0_8
        n_evaluations_rep1 = 0_int64
        n_hits_rep1        = 0_int64
        
        !$ t_pass_start = omp_get_wtime()

        write(0,'(A,I0,A,I0,A)') '=== Processing Group Pass ', igroup, ' of ', plan%n_groups, ' ==='
        write(0,'(A,F8.2,A)')    'Allocating Group RAM: ', real(plan%groups(igroup)%total_elements * 4_int64) / 1.0e9, ' GB'
        write(0,'(A,I0)')        'Active pairs in this group: ', plan%groups(igroup)%n_active_pairs

        allocate(group_hist(plan%groups(igroup)%total_elements))
        group_hist = 0.0_real32

        ! Inner Loop 1: Trajektorienblöcke
        do istate = 1, size(lifting_plan%istates)
            do irep = 1, lifting_plan%nreps(istate)

                ! Schalter für Zeitmessung (NUR Replikat 1 messen)
                is_rep1 = (irep == 1)

                do iblock = 1, lifting_plan%nblock(istate, irep)
                
                     call make_filename_lifted_tra(istate, irep, iblock, lifted_prefix, lifted_file)
                     write(0,*) 'irep, iblock',irep, iblock
                    ! ---------------------------------------------------
                    ! 1. READ I/O
                    ! ---------------------------------------------------
           !         if (is_rep1) then
           !             !$ t_start_section = omp_get_wtime()
           !         end if
                    
                    open(newunit=u_in, file=trim(lifted_file), access='stream', &
                         form='unformatted', status='old')
                    read(u_in) nobj, nframe
                    allocate(tra%pix_present(nobj, nframe))
                    allocate(tra%pix_lagged(nobj, nframe))
                    read(u_in) tra%pix_present
                    read(u_in) tra%pix_lagged
                    close(u_in)

           !         if (is_rep1) then
           !             !$ t_rep1_read = t_rep1_read + (omp_get_wtime() - t_start_section)
           !         end if

                    ! ---------------------------------------------------
                    ! 2. SAMPLING (VOLL INLINED & STRIDE-PRECALCULATED)
           !         ! ---------------------------------------------------
           !         if (is_rep1) then
           !             !$ t_start_section = omp_get_wtime()
           !         end if

                    nframes_to_sample = nframe - lag

                    !$OMP PARALLEL DO SCHEDULE(guided) &
                    !$OMP& PRIVATE(pair_entry, active_idx, i, j, ichunk, pair_size, &
                    !$OMP&         local_start_idx, group_offset, base_pair_offset, &
                    !$OMP&         stride_k2, stride_k3, t, i1, i2, i3, k1, k2, k3, &
                    !$OMP&         global_idx, group_hist_idx) &
                    !$OMP& REDUCTION(+:n_evaluations_rep1, n_hits_rep1)
                    do pair_entry = 1, plan%groups(igroup)%n_active_pairs

                        active_idx = plan%groups(igroup)%active_pair_indices(pair_entry)

                        i = pair_map%active_i(active_idx)
                        j = pair_map%active_j(active_idx)

                        pair_size = pair_map%size(i, j)
                        if (pair_size <= 0) cycle

                        ichunk          = pair_map%chunk(i, j)
                        local_start_idx = pair_map%local_start(i, j)
                        group_offset    = plan%chunk2offset(ichunk)

                        ! --------------------------------------------------
                        ! VORBERECHNUNGEN PRO PAAR (i, j) - VOR DER T-SCHLEIFE!
                        ! --------------------------------------------------
                        ! Strides für die 3D -> 1D Indexrechnung
                        stride_k2 = int(hits%n_lagged(i), int64)
                        stride_k3 = stride_k2 * int(hits%n_present(i), int64)
                        
                        ! Basis-Offset im group_hist Array (enthält start, local_start, offset & -1_i8)
                        base_pair_offset = group_offset + local_start_idx - 1_int64

                        if (is_rep1) then
                            n_evaluations_rep1 = n_evaluations_rep1 + int(nframes_to_sample, int64)
                        end if

                        ! --------------------------------------------------
                 
                        ! --------------------------------------------------
                        do t = 1, nframes_to_sample
                            i1 = tra%pix_lagged(i, t)
                            i2 = tra%pix_present(i, t)
                            i3 = tra%pix_present(j, t)

                            k1 = hits%rank_lagged(i1, i)
                            k2 = hits%rank_present(i2, i)
                            k3 = hits%rank_present(i3, j)

                            ! Support-Check
                            if (k1 <= 0 .or. k2 <= 0 .or. k3 <= 0) cycle

                            ! Ultra-schnelle Indexrechnung
                            global_idx = int(k1, int64) + &
                                         stride_k2 * int(k2 - 1, int64) + &
                                         stride_k3 * int(k3 - 1, int64)

                            group_hist_idx = base_pair_offset + global_idx

                        
                            !$OMP ATOMIC
                            group_hist(group_hist_idx) = group_hist(group_hist_idx) + 1.
                        end do ! t

                    end do ! pair_entry
                    !$OMP END PARALLEL DO

               
                    deallocate(tra%pix_present, tra%pix_lagged)

                end do ! iblock

            end do ! irep
        end do ! istate

        ! ---------------------------------------------------
        ! 3. FLUSH / WRITE I/O
        ! ---------------------------------------------------
        !$ t_start_section = omp_get_wtime()

        write(0,*) 'Flushing finished Chunks of Group', igroup, 'to disk...'
        do idx_in_grp = 1, plan%groups(igroup)%n_chunks_in_group
            ichunk = plan%groups(igroup)%chunk_ids(idx_in_grp)
            group_offset = plan%groups(igroup)%chunk_offsets(idx_in_grp)

            call make_te_chunk_filename(output_prefix, ichunk, chunk_out_file)
            open(newunit=u_out, file=trim(chunk_out_file), access='stream', &
                 form='unformatted', status='replace')

            if (idx_in_grp < plan%groups(igroup)%n_chunks_in_group) then
                write(u_out) group_hist(group_offset + 1_int64 : plan%groups(igroup)%chunk_offsets(idx_in_grp+1))
            else
                write(u_out) group_hist(group_offset + 1_int64 : plan%groups(igroup)%total_elements)
            end if

            close(u_out)
        end do

        !$ t_rep1_write = omp_get_wtime() - t_start_section

        deallocate(group_hist)

        !$ write(0,'(A,F8.2,A)') ' Pass Write Time: ', t_rep1_write, ' s'

    end do ! igroup

end subroutine execute_sampling_plan



    subroutine free_sampling_plan(plan)
        type(sampling_plan_type), intent(inout) :: plan
        integer :: i
        if (allocated(plan%groups)) then
            do i = 1, size(plan%groups)
                if (allocated(plan%groups(i)%chunk_ids)) deallocate(plan%groups(i)%chunk_ids)
                if (allocated(plan%groups(i)%chunk_offsets)) deallocate(plan%groups(i)%chunk_offsets)
            end do
            deallocate(plan%groups)
        end if
        if (allocated(plan%chunk2group)) deallocate(plan%chunk2group)
        if (allocated(plan%chunk2offset)) deallocate(plan%chunk2offset)
    end subroutine free_sampling_plan

    
    
subroutine execute_shuffle_sampling_plan(plan, lifting_plan, hits, pair_map, &
                                          npix, lag, shuffle_block_size, &
                                          nshuffle, lifted_prefix, output_prefix)

    use iso_fortran_env, only: real32, int32, int64
    !$ use omp_lib
    implicit none

    type(sampling_plan_type), intent(in) :: plan
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(hitlist_type),      intent(in) :: hits
    type(pair_map_type),     intent(in) :: pair_map

    integer, intent(in) :: npix
    integer, intent(in) :: lag
    integer, intent(in) :: shuffle_block_size
    integer, intent(in) :: nshuffle

    character(len=*), intent(in) :: lifted_prefix
    character(len=*), intent(in) :: output_prefix

    ! ============================================================
    ! Trajektorien-Daten
    !
    ! target_tra:
    !   lagged + present für den aktuell ausgewerteten target-iblock
    !
    ! source_present:
    !   present-Daten des aktuell benötigten source-iblocks
    !
    ! Damit müssen nie mehr als zwei iblocks gleichzeitig im RAM sein.
    ! ============================================================

    type(lifted_tra_int) :: target_tra

    integer(int32), allocatable :: source_present(:,:)

    character(len=800) :: lifted_file
    character(len=800) :: chunk_out_file

    ! ============================================================
    ! Histogramm eines einzelnen Shuffle-Laufs
    ! ============================================================

    real(real32), allocatable :: shuffle_hist(:)

    ! ============================================================
    ! Globale Zeitstruktur der aktuellen Replikation
    !
    ! shuffle_map(t):
    !     Originaler globaler Frame, der an Shuffle-Zeit t steht.
    !
    ! shuffle_source_iblock(t):
    !     iblock, in dem dieser Originalframe liegt.
    !
    ! shuffle_source_local(t):
    !     lokaler Frame innerhalb dieses iblocks.
    ! ============================================================

    integer, allocatable :: shuffle_map(:)
    integer, allocatable :: shuffle_source_iblock(:)
    integer, allocatable :: shuffle_source_local(:)

    ! ============================================================
    ! Metadaten der iblocks der aktuellen Replikation
    ! ============================================================

    integer, allocatable :: block_nframe(:)
    integer, allocatable :: block_start(:)

    ! ============================================================
    ! Schleifenvariablen
    ! ============================================================

    integer :: igroup
    integer :: istate
    integer :: irep
    integer :: iblock
    integer :: source_iblock
    integer :: ishuffle
    integer :: ichunk
    integer :: idx_in_grp

    integer :: pair_entry
    integer :: active_idx
    integer :: i
    integer :: j
    integer :: pair_size

    integer :: u_in
    integer :: u_out
    integer :: ierr

    integer :: nobj
    integer :: nframe
    integer :: nobj_source
    integer :: nframe_source

    integer :: nblock
    integer :: nframes_rep
    integer :: nframes_to_sample

    integer :: t
    integer :: global_t

    integer :: i1
    integer :: i2
    integer :: i3

    integer :: k1
    integer :: k2
    integer :: k3

    integer :: source_local_t

    ! ============================================================
    ! 64-bit Histogrammindex
    ! ============================================================

    integer(int64) :: stride_k2
    integer(int64) :: stride_k3

    integer(int64) :: global_idx
    integer(int64) :: group_hist_idx
    integer(int64) :: group_offset
    integer(int64) :: local_start_idx
    integer(int64) :: base_pair_offset

    ! ============================================================
    ! Sonstiges
    ! ============================================================

    logical :: source_data_loaded
    logical :: is_rep1

    real(8) :: t_pass_start
    real(8) :: t_start_section
    real(8) :: t_rep1_write

    ! ============================================================
    ! Checks
    ! ============================================================

    if (shuffle_block_size <= 0) then
        write(0,*) 'ERROR: shuffle_block_size must be > 0'
        error stop
    end if

    if (nshuffle <= 0) then
        write(0,*) 'ERROR: nshuffle must be > 0'
        error stop
    end if

    if (lag < 0) then
        write(0,*) 'ERROR: lag must be >= 0'
        error stop
    end if

    write(0,*) '-----------------------------------------------------'
    write(0,*) 'PHASE 3S: Executing Block-Shuffle Sampling'
    write(0,*) '-----------------------------------------------------'
    write(0,'(A,I0)') 'Shuffle block size : ', shuffle_block_size
    write(0,'(A,I0)') 'Number of shuffles : ', nshuffle
    write(0,*) '-----------------------------------------------------'


    ! ============================================================
    ! Outer Loop: Histogrammgruppen
    !
    ! Genau wie in execute_sampling_plan:
    ! nur eine Gruppe gleichzeitig im RAM.
    ! ============================================================

    do igroup = 1, plan%n_groups

        t_pass_start = 0.0_8
        !$ t_pass_start = omp_get_wtime()

        write(0,'(A,I0,A,I0,A)') &
            '=== Processing Shuffle Group Pass ', &
            igroup, ' of ', plan%n_groups, ' ==='

        write(0,'(A,F8.2,A)') &
            'Allocating Shuffle Histogram RAM: ', &
            real(plan%groups(igroup)%total_elements * 4_int64) / 1.0e9, &
            ' GB'

        write(0,'(A,I0)') &
            'Active pairs in this group: ', &
            plan%groups(igroup)%n_active_pairs


        ! ========================================================
        ! Jeder Shuffle-Lauf bekommt ein eigenes Histogramm.
        !
        ! Wichtig:
        ! Dieses Histogramm wird über ALLE Replikationen akkumuliert.
        ! ========================================================

        do ishuffle = 1, nshuffle

            write(0,'(A,I0,A,I0)') &
                '-----------------------------------------------------'
            write(0,'(A,I0,A,I0)') &
                'Starting shuffle ', ishuffle, ' / ', nshuffle
            write(0,'(A,I0,A)') &
                'Histogram size: ', &
                plan%groups(igroup)%total_elements, ' entries'


            allocate(shuffle_hist( &
                plan%groups(igroup)%total_elements))

            shuffle_hist = 0.0_real32


            ! ====================================================
            ! Zustände und Replikationen
            ! ====================================================

            do istate = 1, size(lifting_plan%istates)

                do irep = 1, lifting_plan%nreps(istate)

                    is_rep1 = (irep == 1)

                    write(0,'(A,I0,A,I0,A,I0)') &
                        'Shuffle ', ishuffle, &
                        ': preparing istate=', istate, &
                        ', irep=', irep


                    ! =================================================
                    ! 1. Metadaten der gesamten Replikation einlesen
                    !
                    ! Noch keine Pixeldaten im RAM.
                    ! =================================================

                    nblock = lifting_plan%nblock(istate, irep)

                    allocate(block_nframe(nblock))
                    allocate(block_start(nblock))

                    nframes_rep = 0


                    do iblock = 1, nblock

                        call make_filename_lifted_tra( &
                            istate, irep, iblock, &
                            lifted_prefix, lifted_file)


                        open(newunit=u_in, &
                             file=trim(lifted_file), &
                             access='stream', &
                             form='unformatted', &
                             status='old', &
                             iostat=ierr)

                        if (ierr /= 0) then
                            write(0,*) 'ERROR opening: ', trim(lifted_file)
                            error stop
                        end if


                        read(u_in) nobj, nframe

                        close(u_in)


                        if (iblock == 1) then
                            block_start(iblock) = 1
                        else
                            block_start(iblock) = &
                                block_start(iblock-1) + &
                                block_nframe(iblock-1)
                        end if

                        block_nframe(iblock) = nframe

                        nframes_rep = nframes_rep + nframe

                    end do


                    if (nframes_rep <= lag) then

                        write(0,'(A,I0,A,I0,A,I0)') &
                            'WARNING: istate=', istate, &
                            ', irep=', irep, &
                            ' has only ', nframes_rep, &
                            ' frames. Skipping.'

                        deallocate(block_nframe)
                        deallocate(block_start)

                        cycle

                    end if


                    ! =================================================
                    ! 2. Neue Blockpermutation für DIESE Replikation
                    !    und DIESES Shuffle
                    !
                    ! Diese Permutation wird für alle j gemeinsam
                    ! verwendet.
                    ! =================================================

                    allocate(shuffle_map(nframes_rep))
                    allocate(shuffle_source_iblock(nframes_rep))
                    allocate(shuffle_source_local(nframes_rep))


                    call build_block_shuffle( &
                        nframes_rep, &
                        shuffle_block_size, &
                        shuffle_map)


                    ! -------------------------------------------------
                    ! Für jeden globalen Shuffle-Zeitpunkt bestimmen,
                    ! aus welchem iblock der Source-Wert kommt.
                    !
                    ! Das kostet nur 2 Integer-Arrays.
                    ! -------------------------------------------------

                    do global_t = 1, nframes_rep

                        call locate_replication_frame( &
                            shuffle_map(global_t), &
                            block_start, &
                            block_nframe, &
                            nblock, &
                            shuffle_source_iblock(global_t), &
                            shuffle_source_local(global_t))

                    end do


                    ! =================================================
                    ! 3. Target-iblocks
                    !
                    ! Der aktuelle target-iblock wird vollständig
                    ! eingelesen.
                    ! =================================================

                    do iblock = 1, nblock

                        call make_filename_lifted_tra( &
                            istate, irep, iblock, &
                            lifted_prefix, lifted_file)


                        write(0,'(A,I0,A,I0,A,I0)') &
                            'Shuffle ', ishuffle, &
                            ': target iblock ', iblock, &
                            ' / ', nblock


                        open(newunit=u_in, &
                             file=trim(lifted_file), &
                             access='stream', &
                             form='unformatted', &
                             status='old', &
                             iostat=ierr)

                        if (ierr /= 0) then
                            write(0,*) 'ERROR opening: ', trim(lifted_file)
                            error stop
                        end if


                        read(u_in) nobj, nframe


                        if (iblock == 1) then
                            if (nobj /= size(lifting_plan%istates) .and. &
                                .false.) then
                                ! Kein tatsächlicher Check.
                                ! nobj wird weiter unten geprüft.
                            end if
                        end if


                        allocate(target_tra%pix_present(nobj,nframe))
                        allocate(target_tra%pix_lagged(nobj,nframe))

                        read(u_in) target_tra%pix_present
                        read(u_in) target_tra%pix_lagged

                        close(u_in)


                        nframes_to_sample = nframe - lag


                        ! =================================================
                        ! 4. Source-iblocks
                        !
                        ! Für jeden Source-iblock wird genau einmal
                        ! pix_present geladen.
                        !
                        ! Dadurch können Shuffle-Blöcke beliebig über
                        ! iblock-Grenzen hinweg reichen.
                        ! =================================================

                        do source_iblock = 1, nblock

                            source_data_loaded = .false.


                            ! -------------------------------------------------
                            ! Prüfen, ob der aktuelle Target-iblock überhaupt
                            ! geshuffelte Sourceframes aus diesem Source-iblock
                            ! benötigt.
                            ! -------------------------------------------------

                            do t = 1, nframes_to_sample

                                global_t = &
                                    block_start(iblock) + t - 1

                                if (shuffle_source_iblock(global_t) == &
                                    source_iblock) then

                                    source_data_loaded = .true.
                                    exit

                                end if

                            end do


                            if (.not. source_data_loaded) cycle


                            ! -------------------------------------------------
                            ! Wenn Source und Target im selben iblock liegen,
                            ! können wir die bereits gelesenen Targetdaten
                            ! direkt verwenden.
                            !
                            ! Ansonsten Source-pix_present separat laden.
                            ! -------------------------------------------------

                            if (source_iblock /= iblock) then

                                call make_filename_lifted_tra( &
                                    istate, irep, source_iblock, &
                                    lifted_prefix, lifted_file)


                                open(newunit=u_in, &
                                     file=trim(lifted_file), &
                                     access='stream', &
                                     form='unformatted', &
                                     status='old', &
                                     iostat=ierr)

                                if (ierr /= 0) then
                                    write(0,*) &
                                        'ERROR opening source file: ', &
                                        trim(lifted_file)
                                    error stop
                                end if


                                read(u_in) nobj_source, nframe_source


                                if (nobj_source /= nobj) then
                                    write(0,*) &
                                        'ERROR: nobj differs between iblocks'
                                    write(0,*) &
                                        'target nobj = ', nobj
                                    write(0,*) &
                                        'source nobj = ', nobj_source
                                    error stop
                                end if


                                if (nframe_source /= &
                                    block_nframe(source_iblock)) then

                                    write(0,*) &
                                        'ERROR: inconsistent source nframe'
                                    error stop

                                end if


                                allocate(source_present( &
                                    nobj_source, nframe_source))


                                read(u_in) source_present

                                close(u_in)

                            end if


                            ! =================================================
                            ! 5. Histogramm-Sampling
                            ! =================================================

                            !$OMP PARALLEL DO SCHEDULE(guided) &
                            !$OMP& PRIVATE(pair_entry, active_idx, i, j, &
                            !$OMP&         pair_size, ichunk, local_start_idx, &
                            !$OMP&         group_offset, base_pair_offset, &
                            !$OMP&         stride_k2, stride_k3, t, global_t, &
                            !$OMP&         source_local_t, i1, i2, i3, &
                            !$OMP&         k1, k2, k3, global_idx, &
                            !$OMP&         group_hist_idx)

                            do pair_entry = 1, &
                                plan%groups(igroup)%n_active_pairs


                                active_idx = &
                                    plan%groups(igroup)% &
                                    active_pair_indices(pair_entry)


                                i = pair_map%active_i(active_idx)
                                j = pair_map%active_j(active_idx)


                                pair_size = pair_map%size(i,j)

                                if (pair_size <= 0) cycle


                                ichunk = pair_map%chunk(i,j)

                                local_start_idx = &
                                    pair_map%local_start(i,j)

                                group_offset = &
                                    plan%chunk2offset(ichunk)


                                stride_k2 = &
                                    int(hits%n_lagged(i),int64)

                                stride_k3 = &
                                    stride_k2 * &
                                    int(hits%n_present(i),int64)


                                base_pair_offset = &
                                    group_offset + &
                                    local_start_idx - 1_int64


                                ! -----------------------------------------
                                ! Target-Zeit
                                ! -----------------------------------------

                                do t = 1, nframes_to_sample

                                    global_t = &
                                        block_start(iblock) + t - 1


                                    ! -------------------------------------
                                    ! Nur die Sourceframes auswerten, die
                                    ! aus dem aktuell geladenen source_iblock
                                    ! kommen.
                                    ! -------------------------------------

                                    if (shuffle_source_iblock(global_t) /= &
                                        source_iblock) cycle


                                    source_local_t = &
                                        shuffle_source_local(global_t)


                                    ! -------------------------------------
                                    ! Target bleibt vollständig original.
                                    ! -------------------------------------

                                    i1 = &
                                        target_tra%pix_lagged(i,t)

                                    i2 = &
                                        target_tra%pix_present(i,t)


                                    ! -------------------------------------
                                    ! Source wird aus der geshuffelten
                                    ! zeitlichen Position genommen.
                                    ! -------------------------------------

                                    if (source_iblock == iblock) then

                                        i3 = &
                                            target_tra%pix_present( &
                                                j, source_local_t)

                                    else

                                        i3 = &
                                            source_present( &
                                                j, source_local_t)

                                    end if


                                    ! -------------------------------------
                                    ! Hit-Ränge
                                    ! -------------------------------------

                                    k1 = hits%rank_lagged(i1,i)
                                    k2 = hits%rank_present(i2,i)
                                    k3 = hits%rank_present(i3,j)


                                    if (k1 <= 0 .or. &
                                        k2 <= 0 .or. &
                                        k3 <= 0) cycle


                                    ! -------------------------------------
                                    ! 3D -> 1D Histogrammindex
                                    ! -------------------------------------

                                    global_idx = &
                                        int(k1,int64) + &
                                        stride_k2 * int(k2-1,int64) + &
                                        stride_k3 * int(k3-1,int64)


                                    group_hist_idx = &
                                        base_pair_offset + global_idx


                                    !$OMP ATOMIC
                                    shuffle_hist(group_hist_idx) = &
                                        shuffle_hist(group_hist_idx) + &
                                        1.0_real32

                                end do

                            end do

                            !$OMP END PARALLEL DO


                            ! -------------------------------------------------
                            ! Source iblock wieder freigeben.
                            ! -------------------------------------------------

                            if (source_iblock /= iblock) then
                                deallocate(source_present)
                            end if

                        end do ! source_iblock


                        deallocate(target_tra%pix_present)
                        deallocate(target_tra%pix_lagged)

                    end do ! iblock


                    ! =================================================
                    ! Metadaten dieser Replikation freigeben
                    ! =================================================

                    deallocate(shuffle_map)
                    deallocate(shuffle_source_iblock)
                    deallocate(shuffle_source_local)

                    deallocate(block_nframe)
                    deallocate(block_start)

                end do ! irep

            end do ! istate


            ! ========================================================
            ! 6. Histogramm dieses Shuffle-Laufs schreiben
            !
            ! Erst JETZT, nachdem alle Replikationen akkumuliert wurden.
            ! ========================================================

            write(0,'(A,I0,A)') &
                'Writing completed shuffle histogram ', &
                ishuffle, ' to disk...'


            do idx_in_grp = 1, &
                plan%groups(igroup)%n_chunks_in_group

                ichunk = &
                    plan%groups(igroup)%chunk_ids(idx_in_grp)

                group_offset = &
                    plan%groups(igroup)%chunk_offsets(idx_in_grp)


                call make_shuffle_te_chunk_filename( &
                    output_prefix, &
                    ishuffle, &
                    ichunk, &
                    chunk_out_file)


                open(newunit=u_out, &
                     file=trim(chunk_out_file), &
                     access='stream', &
                     form='unformatted', &
                     status='replace', &
                     iostat=ierr)

                if (ierr /= 0) then
                    write(0,*) &
                        'ERROR opening output: ', &
                        trim(chunk_out_file)
                    error stop
                end if


                if (idx_in_grp < &
                    plan%groups(igroup)%n_chunks_in_group) then

                    write(u_out) &
                        shuffle_hist( &
                            group_offset + 1_int64 : &
                            plan%groups(igroup)%chunk_offsets( &
                                idx_in_grp+1))

                else

                    write(u_out) &
                        shuffle_hist( &
                            group_offset + 1_int64 : &
                            plan%groups(igroup)%total_elements)

                end if


                close(u_out)

            end do


            deallocate(shuffle_hist)

        end do ! ishuffle


        !$ t_rep1_write = omp_get_wtime() - t_pass_start

        write(0,'(A,I0,A,F8.2,A)') &
            'Finished shuffle group ', &
            igroup, &
            ' in ', &
            t_rep1_write, ' s'

    end do ! igroup


end subroutine execute_shuffle_sampling_plan
   
subroutine build_block_shuffle(nframes, block_size, shuffle_map)

    use iso_fortran_env, only: int32

    implicit none

    integer, intent(in)  :: nframes
    integer, intent(in)  :: block_size
    integer, intent(out) :: shuffle_map(:)

    integer :: nblocks
    integer :: b
    integer :: k
    integer :: src_block
    integer :: dst_block

    integer :: src_start
    integer :: src_end
    integer :: dst_start
    integer :: dst_end

    integer :: offset
    integer :: src_frame
    integer :: dst_frame

    integer :: tmp
    real(8) :: r

    integer, allocatable :: permutation(:)


    if (nframes <= 0) then
        write(0,*) 'ERROR: nframes <= 0 in build_block_shuffle'
        error stop
    end if

    if (block_size <= 0) then
        write(0,*) 'ERROR: block_size <= 0 in build_block_shuffle'
        error stop
    end if

    if (size(shuffle_map) /= nframes) then
        write(0,*) 'ERROR: wrong shuffle_map size'
        error stop
    end if


    ! ------------------------------------------------------------
    ! Anzahl Blöcke
    !
    ! Der letzte Block darf kürzer sein.
    ! ------------------------------------------------------------

    nblocks = (nframes + block_size - 1) / block_size

    allocate(permutation(nblocks))


    do b = 1, nblocks
        permutation(b) = b
    end do


    ! ------------------------------------------------------------
    ! Fisher-Yates auf Blockebene
    ! ------------------------------------------------------------

    do k = nblocks, 2, -1

        call random_number(r)

        b = 1 + int(r * real(k,8))

        if (b < 1) b = 1
        if (b > k) b = k

        tmp = permutation(k)
        permutation(k) = permutation(b)
        permutation(b) = tmp

    end do


    ! ------------------------------------------------------------
    ! shuffle_map aufbauen
    !
    ! shuffle_map(dst_frame) = src_frame
    !
    ! Der Offset innerhalb des Blocks bleibt erhalten.
    ! ------------------------------------------------------------

    shuffle_map = 0


    do dst_block = 1, nblocks

        src_block = permutation(dst_block)


        src_start = (src_block - 1) * block_size + 1
        src_end   = min(src_block * block_size, nframes)

        dst_start = (dst_block - 1) * block_size + 1
        dst_end   = min(dst_block * block_size, nframes)


        do src_frame = src_start, src_end

            offset = src_frame - src_start

            dst_frame = dst_start + offset

            if (dst_frame <= dst_end) then
                shuffle_map(dst_frame) = src_frame
            end if

        end do

    end do


    deallocate(permutation)

end subroutine build_block_shuffle

subroutine locate_replication_frame(global_frame, block_start, &
                                    block_nframe, nblock, &
                                    iblock, local_frame)

    implicit none

    integer, intent(in)  :: global_frame
    integer, intent(in)  :: block_start(:)
    integer, intent(in)  :: block_nframe(:)
    integer, intent(in)  :: nblock

    integer, intent(out) :: iblock
    integer, intent(out) :: local_frame

    integer :: lo
    integer :: hi
    integer :: mid


    if (global_frame < 1) then
        write(0,*) 'ERROR: global_frame < 1'
        error stop
    end if


    lo = 1
    hi = nblock


    ! ------------------------------------------------------------
    ! Binärsuche
    ! ------------------------------------------------------------

    do while (lo <= hi)

        mid = (lo + hi) / 2


        if (global_frame < block_start(mid)) then

            hi = mid - 1


        else if (global_frame >= &
                 block_start(mid) + block_nframe(mid)) then

            lo = mid + 1


        else

            iblock = mid

            local_frame = &
                global_frame - block_start(mid) + 1

            return

        end if

    end do


    write(0,*) 'ERROR: global frame not found'
    write(0,*) 'global_frame = ', global_frame
    write(0,*) 'nblock       = ', nblock

    error stop

end subroutine locate_replication_frame

subroutine make_shuffle_te_chunk_filename(output_prefix, &
                                          ishuffle, ichunk, filename)

    implicit none

    character(len=*), intent(in)  :: output_prefix
    integer,          intent(in)  :: ishuffle
    integer,          intent(in)  :: ichunk
    character(len=*), intent(out) :: filename

    write(filename,'(A,"_shuffle",I0,"_chunk",I0)') &
        trim(output_prefix), ishuffle, ichunk

end subroutine make_shuffle_te_chunk_filename
    
end module healpix_histogram_sampling_plan_module