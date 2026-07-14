module healpix_sampling_plan_module
    use direction_vectors_analyse_plan_module ! Enthält project_type
    implicit none
    private
    
    public :: sampling_plan, block_meta, init_sampling_plan
    public :: replica_selector_type,init_replica_selector
    
    ! Die beiden Shuffling-Modi als benannte Konstanten
    integer, parameter, public :: SHUFFLE_INTRA_REPLICA = 1 
    integer, parameter, public :: SHUFFLE_INTER_REPLICA = 2 

    ! ======================================================================
    ! NEU: SELEKTOR-TYP FÜR CROSS-VALIDATION UND SUBSET-ANALYSEN
    ! ======================================================================
    type :: replica_selector_type
        integer :: nrep_total
        integer :: nrep_active
        logical, allocatable :: is_active(:) ! Dimension: (nrep_total) -> .true. falls zu nutzen
    end type replica_selector_type

    type :: block_meta
        integer :: replica_id      ! Aus welchem Replikat (1..nrep) stammt der Block
        integer :: start_frame     ! Globaler Startframe innerhalb des Replikats
        integer :: length          ! Länge des Blocks (deine Block-Länge L)
    end type block_meta

    type :: sampling_plan
        integer :: istate              
        integer :: shuffle_mode        
        integer :: block_length        
        integer :: nblocks_total       
        
        type(block_meta), allocatable :: global_blocks(:) 
        
        integer, allocatable :: real_run_perm(:)    
        integer, allocatable :: shuffle_run_perm(:) 
    end type sampling_plan

contains
! ======================================================================
    ! INITIALISIERUNG DES REPLIKAT-WÄHLERS (Standard: Alle aktiv)
    ! ======================================================================
    subroutine init_replica_selector(selector, nrep_total)
        type(replica_selector_type), intent(out) :: selector
        integer, intent(in)                      :: nrep_total
        
        selector%nrep_total  = nrep_total
        selector%nrep_active = nrep_total
        
        ! Falls die Maske bereits allokiert war (z.B. durch wiederholte Aufrufe), freigeben
        if (allocated(selector%is_active)) deallocate(selector%is_active)
        
        allocate(selector%is_active(nrep_total))
        
        ! Standardmäßig sind alle Replikate für die Analyse aktiv geschaltet
        selector%is_active = .true. 
        
    end subroutine init_replica_selector

    ! ======================================================================
    ! HAUPTROUTINE: Erweitert um den replica_selector
    ! ======================================================================
    subroutine init_sampling_plan(plan, project, istate, L, mode, selector)
        type(sampling_plan), intent(out)          :: plan
        type(project_type), intent(in)            :: project
        integer, intent(in)                       :: istate      
        integer, intent(in)                       :: L           
        integer, intent(in)                       :: mode        
        type(replica_selector_type), intent(in)   :: selector   ! <<< NEU: Selektor übergeben
        
        integer :: irep, nrep, b, idx, nblocks_this_rep, total_frames_use
        integer :: i, j, tmp, rep_start_idx, rep_end_idx
        real(8) :: rand_val

        plan%istate = istate
        plan%block_length = L
        plan%shuffle_mode = mode
        
        nrep = project%state(istate)%nrep
        
        ! Sicherheits-Check: Passt der Selektor zur Struktur?
        if (selector%nrep_total /= nrep) then
            write(0,*) "ERROR: Selektor-Dimension passt nicht zur Anzahl der Replikate im State!"
            stop
        end if
        
        ! 1. SCHRITT: Zähle Blöcke NUR von AKTIVEN Replikaten
        plan%nblocks_total = 0
        do irep = 1, nrep
            if (.not. selector%is_active(irep)) cycle ! ÜBERSPRINGEN falls inaktiv
            
            total_frames_use = project%state(istate)%rep(irep)%nframes_per_rep_use
            plan%nblocks_total = plan%nblocks_total + (total_frames_use / L)
        end do
        
        if (plan%nblocks_total == 0) then
            write(0,*) "WARNING: Keine Blöcke generiert. Alle Replikate inaktiv?"
            return
        end if
        
        ! Speicher allokieren
        allocate(plan%global_blocks(plan%nblocks_total))
        allocate(plan%real_run_perm(plan%nblocks_total))
        allocate(plan%shuffle_run_perm(plan%nblocks_total))
        
        ! 2. SCHRITT: Block-Metadaten exakt belegen (unter Berücksichtigung des Maskierungs-Filters)
        idx = 0
        do irep = 1, nrep
            if (.not. selector%is_active(irep)) cycle ! ÜBERSPRINGEN falls inaktiv
            
            total_frames_use = project%state(istate)%rep(irep)%nframes_per_rep_use
            nblocks_this_rep = total_frames_use / L
            
            do b = 1, nblocks_this_rep
                idx = idx + 1
                plan%global_blocks(idx)%replica_id = irep
                plan%global_blocks(idx)%start_frame = (b - 1) * L + 1
                plan%global_blocks(idx)%length = L
                
                plan%real_run_perm(idx)    = idx
                plan%shuffle_run_perm(idx) = idx
            end do
        end do
        
        ! 3. SCHRITT: Das Shuffling basierend auf dem gewählten Modus würfeln
        if (mode == SHUFFLE_INTER_REPLICA) then
            ! --- INTER-REPLIKAT SHUFFLE ---
            ! Alle Blöcke der AKTIVEN Replikate werden wild gemischt
            do i = plan%nblocks_total, 2, -1
                call random_number(rand_val)
                j = 1 + int(rand_val * i)
                if (j > i) j = i
                
                tmp = plan%shuffle_run_perm(i)
                plan%shuffle_run_perm(i) = plan%shuffle_run_perm(j)
                plan%shuffle_run_perm(j) = tmp
            end do
            
        else if (mode == SHUFFLE_INTRA_REPLICA) then
            ! --- INTRA-REPLIKAT SHUFFLE ---
            ! Da inaktive Replikate gar keine Blöcke in global_blocks erzeugt haben,
            ! können wir einfach sequentiell durch das generierte Array gehen.
            rep_start_idx = 1
            idx = 0
            
            do irep = 1, nrep
                if (.not. selector%is_active(irep)) cycle
                
                total_frames_use = project%state(istate)%rep(irep)%nframes_per_rep_use
                nblocks_this_rep = total_frames_use / L
                rep_end_idx = rep_start_idx + nblocks_this_rep - 1
                
                ! Fisher-Yates lokal für dieses aktive Replikat
                do i = rep_end_idx, rep_start_idx + 1, -1
                    call random_number(rand_val)
                    j = rep_start_idx + int(rand_val * (i - rep_start_idx + 1))
                    if (j > i) j = i
                    
                    tmp = plan%shuffle_run_perm(i)
                    plan%shuffle_run_perm(i) = plan%shuffle_run_perm(j)
                    plan%shuffle_run_perm(j) = tmp
                end do
                
                rep_start_idx = rep_end_idx + 1
            end do
        end if
        
   ! ======================================================================
        ! DIAGNOSTIK & TESTOUTPUT
        ! ======================================================================
        write(0,*) "========================================================"
        write(0,*) "--> SAMPLING PLAN DIAGNOSTICS FOR STATE:", istate
        write(0,*) "    Active Replicas :", count(selector%is_active), " of ", selector%nrep_total
        write(0,*) "    Total Blocks    :", plan%nblocks_total
        write(0,*) "    Block Length L  :", plan%block_length
        write(0,*) "========================================================"
        
        ! Tabellarische Ausgabe für den direkten Vergleich: Echtlauf vs. Shuffle-Lauf
        write(0,*) "  [Pos]  |  Real-Run (Block ID -> Orig-Rep)  |  Shuffle-Run (Block ID -> Orig-Rep)"
        write(0,*) "  -------+-----------------------------------+-----------------------------------"
        
        do i = 1, plan%nblocks_total
            ! Bestimme, welcher Block im Echtlauf an dieser Position ist
            j = plan%real_run_perm(i)
            
            ! Bestimme, welcher Block im Shuffle-Lauf an dieser Position gelandet ist
            tmp = plan%shuffle_run_perm(i)
            
            write(0,"(A,I5,A,I5,A,I3,A,I5,A,I3,A)") &
                "  ", i, "  |  Block ", j, " from Rep ", plan%global_blocks(j)%replica_id, &
                "  |  Block ", tmp, " from Rep ", plan%global_blocks(tmp)%replica_id, " "
        end do
        write(0,*) "========================================================"
        
    end subroutine init_sampling_plan
        


end module healpix_sampling_plan_module