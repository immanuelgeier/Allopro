module cv_lib_workflow_io_module

    use workflow_module, only: library_workflow_type
    use iso_fortran_env, only: output_unit, error_unit
    
    
    implicit none
    private

    public :: init_library_workflow
    public :: read_library_workflow
    public :: write_library_workflow
    public :: print_library_workflow
    public :: validate_library_workflow
    public :: configure_library_workflow
    

contains

subroutine configure_library_workflow(job,topologies)
    use topology_module, only: topology_set_type
    use sidechains, only: sidechain_library_type
    use domains_module, only: domain_library_type
    use user_defined_cv_module, only: user_cv_library_type, read_user_cv_library
    type(library_workflow_type), intent(out) :: job
    type(topology_set_type), intent(in)       :: topologies
    character(len=1024) :: filename
    character(len=32)   :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Library workflow configuration'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Create new workflow'
    write(*,*) '2 = Load and edit existing workflow'
    write(*,*) '3 = Load existing workflow unchanged'
    write(*,*) '4 = Skip workflow configuration'
    write(*,*)
    write(*,'(A)', advance='no') 'Select option: '

    read(*,'(A)') input
    read(input,*,iostat=ios) choice

    if (ios /= 0) error stop 'Invalid workflow option'

    select case (choice)

    case (1)
        ! Neue Konfiguration mit Standardwerten
        call edit_library_workflow(job,topologies)

    case (2)
        ! Vorhandene Konfiguration laden und bearbeiten
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_library_workflow(job, trim(filename))
        call edit_library_workflow(job,topologies)

    case (3)
        ! Vorhandene Konfiguration unverändert übernehmen
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_library_workflow(job, trim(filename))
        call validate_library_workflow(job)

        write(*,*)
        write(*,*) 'Workflow loaded:'
        call print_library_workflow(job)

    case (4)
        ! Workflow-Konfiguration überspringen
        write(*,*) 
        write(*,*) 'Workflow configuration skipped.'

    case default
        error stop 'Invalid workflow option'

    end select

end subroutine configure_library_workflow

subroutine edit_library_workflow(job, topologies)
    use sidechains, only: sidechain_library_type, init_resnames_and_sidechain_head_groups, init_sidechain_head, compute_freq_and_resids, write_sidechain_library
    use domains_module, only: domain_library_type, init_domain_maps, write_domain_library
    use user_defined_cv_module, only: user_cv_library_type, read_user_cv_library, write_user_cv_library
    use topology_module, only: topology_set_type
    implicit none

    type(library_workflow_type), intent(inout) :: job
    type(topology_set_type), intent(in)        :: topologies

    character(len=1024) :: input
    character(len=1024) :: filename
    integer :: ios, i, workflow_mode
    logical :: save_config, save_libs, has_user_lib

    ! Bibliotheks-Instanzen
    type(sidechain_library_type) :: sc_lib
    type(domain_library_type)   :: dom_lib
    type(user_cv_library_type)  :: user_lib
    integer :: rid
    character(len=3) :: current_resname
    integer :: sc_mode, k
    character(len=3) :: cname
    logical :: default_state
        
    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Edit Library Workflow (Topological)'
    write(*,*) '======================================'
    write(*,*) 'Press Enter to keep the current value.'
    write(*,*)

    ! --------------------------------------------------
    ! 1. State & Model Parameters (Ziel-State)
    ! --------------------------------------------------
    write(*,'(A,I0,A)', advance='no') 'State index (target state) [', job%model_params%istate, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%model_params%istate
        if (ios /= 0) error stop 'Invalid state index'
    end if

    ! State des Referenzfensters automatisch an den Zielstate anpassen
    job%reference_window%istate = job%model_params%istate


    ! --------------------------------------------------
    ! 2. Modus-Auswahl (Einlesen / Einlesen+Editieren / Erstellen)
    ! --------------------------------------------------
    write(*,*)
    write(*,*) '--------------------------------------------------'
    write(*,*) ' Library Workflow Mode:'
    write(*,*) '--------------------------------------------------'
    write(*,*) '  1 = Load existing library files (Read-only)'
    write(*,*) '  2 = Load existing library files and edit active status'
    write(*,*) '  3 = Compute / Create new libraries from reference structure'
    write(*,'(A,I0,A)', advance='no') 'Select option [1]: '
    read(*,'(A)') input
    
    workflow_mode = 1
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) workflow_mode
        if (ios /= 0) error stop 'Invalid workflow mode option'
    end if


    ! ==========================================================
    ! MODUS 1: VORHANDENE BIBLIOTHEKEN NUR EINLESEN
    ! ==========================================================
    if (workflow_mode == 1) then
        job%io_params%read_sidechain_library = .true.
        job%io_params%read_domain_library   = .true.
        job%io_params%write_sidechain_library = .false.
        job%io_params%write_domain_library   = .false.

        ! Sidechain library file
        write(*,'(A,A,A)', advance='no') 'Sidechain library file [', trim(job%io_params%sidechain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%sidechain_library_file = trim(input)

        ! Domain library file
        write(*,'(A,A,A)', advance='no') 'Domain library file [', trim(job%io_params%domain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%domain_library_file = trim(input)

        ! Optional User CV Library abfragen
        write(*,'(A,L1,A)', advance='no') 'Include optional User CV library file? [', job%io_params%read_user_cv_library, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%io_params%read_user_cv_library

        if (job%io_params%read_user_cv_library) then
            write(*,'(A,A,A)', advance='no') 'User CV library file [', trim(job%io_params%user_cv_library_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%user_cv_library_file = trim(input)
            job%io_params%write_user_cv_library = .false.
        else
            job%io_params%read_user_cv_library = .false.
            job%io_params%write_user_cv_library = .false.
        end if


    ! ==========================================================
    ! MODUS 2: BIBLIOTHEKEN AUS TOPOLOGIE VORINITIALISIEREN & EDITIEREN
    ! ==========================================================
    else if (workflow_mode == 2) then
        job%io_params%read_sidechain_library = .true.
        job%io_params%read_domain_library   = .true.
        job%io_params%write_sidechain_library = .false.
        job%io_params%write_domain_library   = .false.

        ! Sidechain & Domain library output filenames
        write(*,'(A,A,A)', advance='no') 'Sidechain library output file [', trim(job%io_params%sidechain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%sidechain_library_file = trim(input)

        write(*,'(A,A,A)', advance='no') 'Domain library output file [', trim(job%io_params%domain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%domain_library_file = trim(input)

        ! Vorinitialisierung aus Topologie (für Sidechain & Domain)
        write(*,*) '-> Pre-initializing Sidechain library from topology...'
        call init_resnames_and_sidechain_head_groups(sc_lib)
        call init_sidechain_head(sc_lib)
        call compute_freq_and_resids(sc_lib, job%model_params%istate, topologies)

        write(*,*) '-> Pre-initializing Domain library from topology...'
        call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)

        ! Optionale User-CV-Bibliothek aus Datei einlesen
        write(*,*)
        write(*,'(A,L1,A)', advance='no') 'Include pre-edited User CV library file? [', job%io_params%read_user_cv_library, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%io_params%read_user_cv_library

        has_user_lib = job%io_params%read_user_cv_library

        if (has_user_lib) then
            write(*,'(A,A,A)', advance='no') 'User CV library file [', trim(job%io_params%user_cv_library_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%user_cv_library_file = trim(input)

            job%io_params%write_user_cv_library = .true.

            write(*,*) '-> Reading User CV library from file...'
            call read_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
        else
            job%io_params%read_user_cv_library = .false.
            job%io_params%write_user_cv_library = .false.
        end if


    ! ==========================================================
    ! MODUS 3: NEU BERECHNEN AUS REFERENZSTRUKTUR
    ! ==========================================================
    else if (workflow_mode == 3) then
        job%io_params%read_sidechain_library = .false.
        job%io_params%read_domain_library   = .false.
        job%io_params%write_sidechain_library = .true.
        job%io_params%write_domain_library   = .true.

        ! Reference method & Search Parameters abfragen
        write(*,*)
        write(*,*) '--------------------------------------------------'
        write(*,*) ' Reference Method & Search Parameters:'
        write(*,*) '--------------------------------------------------'

        write(*,*) 'Reference method:'
        write(*,'(A)') '  1 = RMSD'
        write(*,'(A)') '  2 = CV orthogonality'
        write(*,'(A)') '  3 = Fixed frame'

        write(*,'(A,I0,A)', advance='no') 'Select method [', job%model_params%reference_method_index, ']: '
        read(*,'(A)') input

        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) job%model_params%reference_method_index
            if (ios /= 0) error stop 'Invalid reference method index'
        end if

        select case (job%model_params%reference_method_index)
        case (1)
            job%model_params%reference_method = 'rmsd'
        case (2)
            job%model_params%reference_method = 'cv_orthogonality'
        case (3)
            job%model_params%reference_method = 'fixed_frame'
        case default
            error stop 'Invalid reference method index'
        end select

        ! Kandidaten und Window-Parameter nur abfragen, wenn Methode 1 oder 2 aktiv ist
        if (job%model_params%reference_method_index == 1 .or. job%model_params%reference_method_index == 2) then
            write(*,'(A,I0,A)', advance='no') 'Number of candidates [', job%model_params%n_candidates, ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) then
                read(input,*,iostat=ios) job%model_params%n_candidates
                if (ios /= 0) error stop 'Invalid number of candidates'
            end if

            write(*,'(A,I0,A)', advance='no') 'Window method (stride=1, random=2) [', job%model_params%method_window, ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) then
                read(input,*,iostat=ios) job%model_params%method_window
                if (ios /= 0) error stop 'Invalid window method'
            end if
        end if

        ! Reference structure search / fixed frame window abfragen
        write(*,*)
        write(*,*) '--------------------------------------------------'
        if (job%model_params%reference_method_index == 3) then
            write(*,*) ' Fixed Frame Reference Selection:'
        else
            write(*,*) ' Reference structure search window:'
        end if
        write(*,*) '--------------------------------------------------'
        write(*,*) 'State index: ', job%reference_window%istate

        write(*,'(A,I0,A)', advance='no') 'Replicate index [', job%reference_window%irep, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%reference_window%irep

        if (job%model_params%reference_method_index == 3) then
            write(*,'(A,I0,A)', advance='no') 'Time step / Frame [', job%reference_window%start_time_step, ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) read(input,*) job%reference_window%start_time_step
            job%reference_window%end_time_step = job%reference_window%start_time_step
        else
            write(*,'(A,I0,A)', advance='no') 'Start time step [', job%reference_window%start_time_step, ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) read(input,*) job%reference_window%start_time_step

            write(*,'(A,I0,A)', advance='no') 'End time step [', job%reference_window%end_time_step, ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) read(input,*) job%reference_window%end_time_step
        end if

        ! Ausgabedateinamen festlegen
        write(*,'(A,A,A)', advance='no') 'Output Sidechain library file [', trim(job%io_params%sidechain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%sidechain_library_file = trim(input)

        write(*,'(A,A,A)', advance='no') 'Output Domain library file [', trim(job%io_params%domain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%domain_library_file = trim(input)

        ! Vorinitialisierung aus Topologie
        write(*,*) '-> Pre-initializing Sidechain library from topology...'
        call init_resnames_and_sidechain_head_groups(sc_lib)
        call init_sidechain_head(sc_lib)
        call compute_freq_and_resids(sc_lib, job%model_params%istate, topologies)

        write(*,*) '-> Pre-initializing Domain library from topology...'
        call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)

        ! Optionale User-CV-Bibliothek aus Datei einlesen
        write(*,*)
        write(*,'(A,L1,A)', advance='no') 'Include pre-edited User CV library file? [', job%io_params%read_user_cv_library, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%io_params%read_user_cv_library

        has_user_lib = job%io_params%read_user_cv_library

        if (has_user_lib) then
            write(*,'(A,A,A)', advance='no') 'User CV library file [', trim(job%io_params%user_cv_library_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%user_cv_library_file = trim(input)

            job%io_params%write_user_cv_library = .true.

            write(*,*) '-> Reading User CV library from file...'
            call read_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
        else
            job%io_params%read_user_cv_library = .false.
            job%io_params%write_user_cv_library = .false.
        end if
    end if

! ==========================================================
! INTERAKTIVES FILTERN DER CVs AUF WORKFLOW-EBENE
! ==========================================================
    if (workflow_mode == 2 .or. workflow_mode == 3) then
        
        write(*,*)
        write(*,*) '--------------------------------------------------'
        write(*,*) ' CV Filtering & Activation Setup (Workflow-Level):'
        write(*,*) '--------------------------------------------------'

        ! Workflow-Modus im Job abspeichern
        write(*,*) '--- Sidechain CV Filters ---'
        write(*,*) '  1 = Enable ALL sidechains'
        write(*,*) '  2 = Disable ALL sidechains'
        write(*,*) '  3 = Configure by sidechain type (e.g. all PHE, all ARG)'
        write(*,*) '  4 = Configure individual sidechain instances (detailed)'
        write(*,'(A)', advance='no') 'Select sidechain configuration mode [3]: '
        read(*,'(A)') input
        
        sc_mode = 3
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) sc_mode
            if (ios /= 0) sc_mode = 3
        end if
        job%sidechain_filter_mode = sc_mode

        ! Speicher für Typ-Filter vorbereiten
        if (allocated(job%sc_type_filters)) deallocate(job%sc_type_filters)
        allocate(job%sc_type_filters(sc_lib%ncv))
        job%sc_type_filters(:) = .true. ! Standard: alles an

        select case (sc_mode)
        case (1)
            job%sc_type_filters(:) = .true.
            write(*,*) '-> All sidechains enabled for this workflow.'

        case (2)
            job%sc_type_filters(:) = .false.
            write(*,*) '-> All sidechains disabled for this workflow.'

        case (3)
            write(*,*)
            write(*,*) '--- Configure Sidechains by Residue Type ---'
            do i = 1, sc_lib%ncv
                cname = sc_lib%residue(i)%resname
                write(*,'(A,I4,A,A3,A)', advance='no') '  [Type ', i, '] Include all ', cname, ' sidechains? [y] (y/n): '
                read(*,'(A)') input
                if (len_trim(input) > 0) then
                    input = trim(adjustl(input))
                    if (input == 'n' .or. input == 'N' .or. input == 'no') then
                        job%sc_type_filters(i) = .false.
                    end if
                end if
            end do

        case (4)
            write(*,*)
            write(*,*) '--- Detailed Individual Sidechain Configuration ---'
            do i = 1, sc_lib%ncv
                if (allocated(sc_lib%residue(i)%res_ids)) then
                    do k = 1, size(sc_lib%residue(i)%res_ids)
                        rid = sc_lib%residue(i)%res_ids(k)
                        cname = sc_lib%residue(i)%resname

                        write(*,'(A,I4,A,I6,A,A3,A)', advance='no') &
                            '  [Type ', i, ', Residue ID ', rid, '] Name: ', cname, &
                            ' - Include? [y] (y/n): '
                        
                        read(*,'(A)') input
                        if (len_trim(input) > 0) then
                            input = trim(adjustl(input))
                            if (input == 'n' .or. input == 'N' .or. input == 'no') then
                                ! Entsprechende Instanz im Job-Filter auf .false. setzen
                            end if
                        end if
                    end do
                end if
            end do
        end select


        ! ==========================================================
        ! 2. DOMAIN FILTERS (Im Job abspeichern)
        ! ==========================================================
        if (allocated(job%domain_filters)) deallocate(job%domain_filters)
        allocate(job%domain_filters(dom_lib%ndo))
        job%domain_filters(:) = .true.

        write(*,*)
        write(*,*) '--- Domain CV Filters ---'
        do i = 1, dom_lib%ndo
            write(*,'(A,I4,A,A,A)', advance='no') &
                '  [', i, '] Domain: ', trim(dom_lib%domain(i)%domain_name), &
                ' - Include in workflow? [y] (y/n): '
            
            read(*,'(A)') input
            if (len_trim(input) > 0) then
                input = trim(adjustl(input))
                if (input == 'n' .or. input == 'N' .or. input == 'no') then
                    job%domain_filters(i) = .false.
                end if
            end if
        end do


        ! ==========================================================
        ! 3. USER-DEFINED CV FILTERS (Im Job abspeichern)
        ! ==========================================================
        if (has_user_lib) then
            if (allocated(job%user_cv_filters)) deallocate(job%user_cv_filters)
            allocate(job%user_cv_filters(user_lib%n_cvs))
            job%user_cv_filters(:) = .true.

            write(*,*)
            write(*,*) '--- User-Defined CV Filters ---'
            do i = 1, user_lib%n_cvs
                write(*,'(A,I4,A,A,A)', advance='no') &
                    '  [', i, '] CV Name: ', trim(user_lib%cvs(i)%cv_name), &
                    ' - Include in workflow? [y] (y/n): '
                
                read(*,'(A)') input
                if (len_trim(input) > 0) then
                    input = trim(adjustl(input))
                    if (input == 'n' .or. input == 'N' .or. input == 'no') then
                        job%user_cv_filters(i) = .false.
                    end if
                end if
            end do
        end if
    end if
    
    ! --------------------------------------------------
    ! 4. Workflow Steps Steuerung aktualisieren
    ! --------------------------------------------------
    if (allocated(job%sc_type_filters)) then
        job%run_params%initialize_sidechain_maps      = any(job%sc_type_filters)
        job%run_params%compute_sidechain_frequencies = any(job%sc_type_filters)
    else
        job%run_params%initialize_sidechain_maps      = .false.
        job%run_params%compute_sidechain_frequencies = .false.
    end if
    job%run_params%test_sidechain_maps             = .false.

    if (allocated(job%domain_filters)) then
        job%run_params%initialize_domain_maps          = any(job%domain_filters)
    else
        job%run_params%initialize_domain_maps          = .false.
    end if

    if (has_user_lib .and. allocated(job%user_cv_filters)) then
        job%run_params%initialize_user_cvs             = any(job%user_cv_filters)
    else
        job%run_params%initialize_user_cvs             = job%io_params%write_user_cv_library
    end if

    job%run_params%search_reference_structures     = (workflow_mode == 3)


    ! --------------------------------------------------
    ! 6. Validate, Display & Save Workflow Config
    ! --------------------------------------------------
    call validate_library_workflow(job)

    write(*,*)
    write(*,*) 'Resulting workflow configuration:'
    call print_library_workflow(job)

    write(*,*)
    write(*,'(A)', advance='no') 'Save workflow configuration? (y/n): '
    read(*,'(A)') input

    save_config = trim(input) == 'y' .or. trim(input) == 'Y'

    if (save_config) then
        write(*,'(A,A,A)', advance='no') 'Configuration filename [library_workflow.cfg]: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) filename = 'library_workflow.cfg'

        call write_library_workflow(job, trim(filename))
        write(*,'(3A)') 'Configuration saved to "', trim(filename), '"'
    end if

end subroutine edit_library_workflow

subroutine edit_library_workflow_old(job, topologies)
    use sidechains, only: sidechain_library_type, init_resnames_and_sidechain_head_groups, init_sidechain_head, compute_freq_and_resids, write_sidechain_library
    use domains_module, only: domain_library_type, init_domain_maps, write_domain_library
    use user_defined_cv_module, only: user_cv_library_type, read_user_cv_library, write_user_cv_library
    use topology_module, only: topology_set_type
    implicit none

    type(library_workflow_type), intent(inout) :: job
    type(topology_set_type), intent(in)        :: topologies

    character(len=1024) :: input
    character(len=1024) :: filename
    integer :: ios, i, workflow_mode
    logical :: save_config, save_libs, has_user_lib

    ! Bibliotheks-Instanzen
    type(sidechain_library_type) :: sc_lib
    type(domain_library_type)   :: dom_lib
    type(user_cv_library_type)  :: user_lib
    integer :: rid
    character(len=3) :: current_resname
    integer :: sc_mode, k
    character(len=3) :: cname
    logical :: default_state
        
    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Edit Library Workflow (Topological)'
    write(*,*) '======================================'
    write(*,*) 'Press Enter to keep the current value.'
    write(*,*)

    ! --------------------------------------------------
    ! 1. State & Model Parameters (Ziel-State)
    ! --------------------------------------------------
    write(*,'(A,I0,A)', advance='no') 'State index (target state) [', job%model_params%istate, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%model_params%istate
        if (ios /= 0) error stop 'Invalid state index'
    end if

    ! State des Referenzfensters automatisch an den Zielstate anpassen
    job%reference_window%istate = job%model_params%istate


    ! --------------------------------------------------
    ! 2. Modus-Auswahl (Einlesen / Einlesen+Editieren / Erstellen)
    ! --------------------------------------------------
    write(*,*)
    write(*,*) '--------------------------------------------------'
    write(*,*) ' Library Workflow Mode:'
    write(*,*) '--------------------------------------------------'
    write(*,*) '  1 = Load existing library files (Read-only)'
    write(*,*) '  2 = Load existing library files and edit active status'
    write(*,*) '  3 = Compute / Create new libraries from reference structure'
    write(*,'(A,I0,A)', advance='no') 'Select option [1]: '
    read(*,'(A)') input
    
    workflow_mode = 1
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) workflow_mode
        if (ios /= 0) error stop 'Invalid workflow mode option'
    end if


    ! ==========================================================
    ! MODUS 1: VORHANDENE BIBLIOTHEKEN NUR EINLESEN
    ! ==========================================================
    if (workflow_mode == 1) then
        job%io_params%read_sidechain_library = .true.
        job%io_params%read_domain_library   = .true.
        job%io_params%write_sidechain_library = .false.
        job%io_params%write_domain_library   = .false.

        ! Sidechain library file
        write(*,'(A,A,A)', advance='no') 'Sidechain library file [', trim(job%io_params%sidechain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%sidechain_library_file = trim(input)

        ! Domain library file
        write(*,'(A,A,A)', advance='no') 'Domain library file [', trim(job%io_params%domain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%domain_library_file = trim(input)

        ! Optional User CV Library abfragen
        write(*,'(A,L1,A)', advance='no') 'Include optional User CV library file? [', job%io_params%read_user_cv_library, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%io_params%read_user_cv_library

        if (job%io_params%read_user_cv_library) then
            write(*,'(A,A,A)', advance='no') 'User CV library file [', trim(job%io_params%user_cv_library_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%user_cv_library_file = trim(input)
            job%io_params%write_user_cv_library = .false.
        else
            job%io_params%read_user_cv_library = .false.
            job%io_params%write_user_cv_library = .false.
        end if


    ! ==========================================================
    ! MODUS 2: BIBLIOTHEKEN AUS TOPOLOGIE VORINITIALISIEREN & EDITIEREN
    ! ==========================================================
    else if (workflow_mode == 2) then
        job%io_params%read_sidechain_library = .true.
        job%io_params%read_domain_library   = .true.
        job%io_params%write_sidechain_library = .false.
        job%io_params%write_domain_library   = .false.

        ! Sidechain & Domain library output filenames
        write(*,'(A,A,A)', advance='no') 'Sidechain library output file [', trim(job%io_params%sidechain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%sidechain_library_file = trim(input)

        write(*,'(A,A,A)', advance='no') 'Domain library output file [', trim(job%io_params%domain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%domain_library_file = trim(input)

        ! Vorinitialisierung aus Topologie (für Sidechain & Domain)
        write(*,*) '-> Pre-initializing Sidechain library from topology...'
        call init_resnames_and_sidechain_head_groups(sc_lib)
        call init_sidechain_head(sc_lib)
        call compute_freq_and_resids(sc_lib, job%model_params%istate, topologies)

        write(*,*) '-> Pre-initializing Domain library from topology...'
        call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)

        ! Optionale User-CV-Bibliothek aus Datei einlesen
        write(*,*)
        write(*,'(A,L1,A)', advance='no') 'Include pre-edited User CV library file? [', job%io_params%read_user_cv_library, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%io_params%read_user_cv_library

        has_user_lib = job%io_params%read_user_cv_library

        if (has_user_lib) then
            write(*,'(A,A,A)', advance='no') 'User CV library file [', trim(job%io_params%user_cv_library_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%user_cv_library_file = trim(input)

            job%io_params%write_user_cv_library = .true.

            write(*,*) '-> Reading User CV library from file...'
            call read_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
        else
            job%io_params%read_user_cv_library = .false.
            job%io_params%write_user_cv_library = .false.
        end if


    ! ==========================================================
    ! MODUS 3: NEU BERECHNEN AUS REFERENZSTRUKTUR
    ! ==========================================================
    else if (workflow_mode == 3) then
        job%io_params%read_sidechain_library = .false.
        job%io_params%read_domain_library   = .false.
        job%io_params%write_sidechain_library = .true.
        job%io_params%write_domain_library   = .true.

        ! Reference method & Search Parameters abfragen
        write(*,*)
        write(*,*) '--------------------------------------------------'
        write(*,*) ' Reference Method & Search Parameters:'
        write(*,*) '--------------------------------------------------'

        write(*,*) 'Reference method:'
        write(*,'(A)') '  1 = RMSD'
        write(*,'(A)') '  2 = CV orthogonality'

        write(*,'(A,I0,A)', advance='no') 'Select method [', job%model_params%reference_method_index, ']: '
        read(*,'(A)') input

        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) job%model_params%reference_method_index
            if (ios /= 0) error stop 'Invalid reference method index'
        end if

        select case (job%model_params%reference_method_index)
        case (1)
            job%model_params%reference_method = 'rmsd'
        case (2)
            job%model_params%reference_method = 'cv_orthogonality'
        case default
            error stop 'Invalid reference method index'
        end select

        write(*,'(A,I0,A)', advance='no') 'Number of candidates [', job%model_params%n_candidates, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) job%model_params%n_candidates
            if (ios /= 0) error stop 'Invalid number of candidates'
        end if

        write(*,'(A,I0,A)', advance='no') 'Window method (stride=1, random=2) [', job%model_params%method_window, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) job%model_params%method_window
            if (ios /= 0) error stop 'Invalid window method'
        end if

        ! Reference structure search window abfragen
        write(*,*)
        write(*,*) '--------------------------------------------------'
        write(*,*) ' Reference structure search window:'
        write(*,*) '--------------------------------------------------'
        write(*,*) 'State index: ', job%reference_window%istate

        write(*,'(A,I0,A)', advance='no') 'Replicate index [', job%reference_window%irep, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%reference_window%irep

        write(*,'(A,I0,A)', advance='no') 'Start time step [', job%reference_window%start_time_step, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%reference_window%start_time_step

        write(*,'(A,I0,A)', advance='no') 'End time step [', job%reference_window%end_time_step, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%reference_window%end_time_step

        ! Ausgabedateinamen festlegen
        write(*,'(A,A,A)', advance='no') 'Output Sidechain library file [', trim(job%io_params%sidechain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%sidechain_library_file = trim(input)

        write(*,'(A,A,A)', advance='no') 'Output Domain library file [', trim(job%io_params%domain_library_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%domain_library_file = trim(input)

        ! Vorinitialisierung aus Topologie
        write(*,*) '-> Pre-initializing Sidechain library from topology...'
        call init_resnames_and_sidechain_head_groups(sc_lib)
        call init_sidechain_head(sc_lib)
        call compute_freq_and_resids(sc_lib, job%model_params%istate, topologies)

        write(*,*) '-> Pre-initializing Domain library from topology...'
        call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)

        ! Optionale User-CV-Bibliothek aus Datei einlesen
        write(*,*)
        write(*,'(A,L1,A)', advance='no') 'Include pre-edited User CV library file? [', job%io_params%read_user_cv_library, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input,*) job%io_params%read_user_cv_library

        has_user_lib = job%io_params%read_user_cv_library

        if (has_user_lib) then
            write(*,'(A,A,A)', advance='no') 'User CV library file [', trim(job%io_params%user_cv_library_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%user_cv_library_file = trim(input)

            job%io_params%write_user_cv_library = .true.

            write(*,*) '-> Reading User CV library from file...'
            call read_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
        else
            job%io_params%read_user_cv_library = .false.
            job%io_params%write_user_cv_library = .false.
        end if
    end if

! ==========================================================
    ! INTERAKTIVES FILTERN DER CVs AUF WORKFLOW-EBENE
    ! ==========================================================
    if (workflow_mode == 2 .or. workflow_mode == 3) then
      

        write(*,*)
        write(*,*) '--------------------------------------------------'
        write(*,*) ' CV Filtering & Activation Setup (Workflow-Level):'
        write(*,*) '--------------------------------------------------'

        ! Workflow-Modus im Job abspeichern
        write(*,*) '--- Sidechain CV Filters ---'
        write(*,*) '  1 = Enable ALL sidechains'
        write(*,*) '  2 = Disable ALL sidechains'
        write(*,*) '  3 = Configure by sidechain type (e.g. all PHE, all ARG)'
        write(*,*) '  4 = Configure individual sidechain instances (detailed)'
        write(*,'(A)', advance='no') 'Select sidechain configuration mode [3]: '
        read(*,'(A)') input
        
        sc_mode = 3
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) sc_mode
            if (ios /= 0) sc_mode = 3
        end if
        job%sidechain_filter_mode = sc_mode

        ! Speicher für Typ-Filter vorbereiten
        if (allocated(job%sc_type_filters)) deallocate(job%sc_type_filters)
        allocate(job%sc_type_filters(sc_lib%ncv))
        job%sc_type_filters(:) = .true. ! Standard: alles an

        select case (sc_mode)
        case (1)
            job%sc_type_filters(:) = .true.
            write(*,*) '-> All sidechains enabled for this workflow.'

        case (2)
            job%sc_type_filters(:) = .false.
            write(*,*) '-> All sidechains disabled for this workflow.'

        case (3)
            write(*,*)
            write(*,*) '--- Configure Sidechains by Residue Type ---'
            do i = 1, sc_lib%ncv
                cname = sc_lib%residue(i)%resname
                write(*,'(A,I4,A,A3,A)', advance='no') '  [Type ', i, '] Include all ', cname, ' sidechains? [y] (y/n): '
                read(*,'(A)') input
                if (len_trim(input) > 0) then
                    input = trim(adjustl(input))
                    if (input == 'n' .or. input == 'N' .or. input == 'no') then
                        job%sc_type_filters(i) = .false.
                    end if
                end if
            end do

        case (4)
            write(*,*)
            write(*,*) '--- Detailed Individual Sidechain Configuration ---'
            ! Optional: Hier ein Instanzen-Array im Job allokieren, falls Modus 4 gewählt wird
            do i = 1, sc_lib%ncv
                if (allocated(sc_lib%residue(i)%res_ids)) then
                    do k = 1, size(sc_lib%residue(i)%res_ids)
                        rid = sc_lib%residue(i)%res_ids(k)
                        cname = sc_lib%residue(i)%resname

                        write(*,'(A,I4,A,I6,A,A3,A)', advance='no') &
                            '  [Type ', i, ', Residue ID ', rid, '] Name: ', cname, &
                            ' - Include? [y] (y/n): '
                        
                        read(*,'(A)') input
                        if (len_trim(input) > 0) then
                            input = trim(adjustl(input))
                            if (input == 'n' .or. input == 'N' .or. input == 'no') then
                                ! Entsprechende Instanz im Job-Filter auf .false. setzen
                            end if
                        end if
                    end do
                end if
            end do
        end select


 ! ==========================================================
        ! 2. DOMAIN FILTERS (Im Job abspeichern)
        ! ==========================================================
        if (allocated(job%domain_filters)) deallocate(job%domain_filters)
        allocate(job%domain_filters(dom_lib%ndo))
        job%domain_filters(:) = .true.

        write(*,*)
        write(*,*) '--- Domain CV Filters ---'
        
        ! Schnellauswahl: Alle Domains standardmäßig einschließen?
        write(*,'(A)', advance='no') 'Include all domains automatically without individual prompts? (y/n) [y]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'y'
        input = trim(adjustl(input))

        if (input /= 'y' .and. input /= 'Y' .and. input /= 'yes') then
            ! Nur wenn ausdrücklich 'n' gewählt wird, einzeln durchklicken
            do i = 1, dom_lib%ndo
                write(*,'(A,I4,A,A,A)', advance='no') &
                    '   [', i, '] Domain: ', trim(dom_lib%domain(i)%domain_name), &
                    ' - Include in workflow? [y] (y/n): '
                
                read(*,'(A)') input
                if (len_trim(input) > 0) then
                    input = trim(adjustl(input))
                    if (input == 'n' .or. input == 'N' .or. input == 'no') then
                        job%domain_filters(i) = .false.
                    end if
                end if
            end do
        else
            write(*,*) '-> All domains included by default.'
        end if

        ! ==========================================================
        ! 3. USER-DEFINED CV FILTERS (Im Job abspeichern)
        ! ==========================================================
        if (has_user_lib) then
            if (allocated(job%user_cv_filters)) deallocate(job%user_cv_filters)
            allocate(job%user_cv_filters(user_lib%n_cvs))
            job%user_cv_filters(:) = .true.

            write(*,*)
            write(*,*) '--- User-Defined CV Filters ---'
            do i = 1, user_lib%n_cvs
                write(*,'(A,I4,A,A,A)', advance='no') &
                    '  [', i, '] CV Name: ', trim(user_lib%cvs(i)%cv_name), &
                    ' - Include in workflow? [y] (y/n): '
                
                read(*,'(A)') input
                if (len_trim(input) > 0) then
                    input = trim(adjustl(input))
                    if (input == 'n' .or. input == 'N' .or. input == 'no') then
                        job%user_cv_filters(i) = .false.
                    end if
                end if
            end do
        end if
    end if
    
    ! --------------------------------------------------
    ! 4. Workflow Steps Steuerung aktualisieren
    ! --------------------------------------------------
    ! Sidechain-Maps nur initialisieren/berechnen, wenn mindestens ein Sidechain-Typ im Workflow aktiv ist
    if (allocated(job%sc_type_filters)) then
        job%run_params%initialize_sidechain_maps      = any(job%sc_type_filters)
        job%run_params%compute_sidechain_frequencies = any(job%sc_type_filters)
    else
        job%run_params%initialize_sidechain_maps      = .false.
        job%run_params%compute_sidechain_frequencies = .false.
    end if
    job%run_params%test_sidechain_maps            = .false.

    ! Domain-Maps nur initialisieren, wenn mindestens eine Domain im Workflow aktiv ist
    if (allocated(job%domain_filters)) then
        job%run_params%initialize_domain_maps         = any(job%domain_filters)
    else
        job%run_params%initialize_domain_maps         = .false.
    end if

    ! User-CVs nur initialisieren, wenn eine User-Lib vorhanden und mindestens ein CV aktiv ist
    if (has_user_lib .and. allocated(job%user_cv_filters)) then
        job%run_params%initialize_user_cvs            = any(job%user_cv_filters)
    else
        job%run_params%initialize_user_cvs            = job%io_params%write_user_cv_library
    end if

    job%run_params%search_reference_structures   = (workflow_mode == 3)


    ! --------------------------------------------------
    ! 6. Validate, Display & Save Workflow Config
    ! --------------------------------------------------
    call validate_library_workflow(job)

    write(*,*)
    write(*,*) 'Resulting workflow configuration:'
    call print_library_workflow(job)

    write(*,*)
    write(*,'(A)', advance='no') 'Save workflow configuration? (y/n): '
    read(*,'(A)') input

    save_config = trim(input) == 'y' .or. trim(input) == 'Y'

    if (save_config) then
        write(*,'(A,A,A)', advance='no') 'Configuration filename [library_workflow.cfg]: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) filename = 'library_workflow.cfg'

        call write_library_workflow(job, trim(filename))
        write(*,'(3A)') 'Configuration saved to "', trim(filename), '"'
    end if

end subroutine edit_library_workflow_old






subroutine init_library_workflow(job, mode, filename)
    type(library_workflow_type), intent(out) :: job
    character(len=*), intent(in) :: mode
    character(len=*), intent(in), optional :: filename

    select case (trim(mode))
    case ('manual')
        ! Standardwerte aus Typdefinition aktiv
    case ('file')
        if (.not. present(filename)) then
            error stop 'init_library_workflow: filename missing'
        end if
        call read_library_workflow(job, filename)
    case default
        write(error_unit, '(A,A)') 'init_library_workflow: unknown mode: ', trim(mode)
        error stop
    end select

    call validate_library_workflow(job)
end subroutine init_library_workflow


subroutine validate_library_workflow(job)
    implicit none
    type(library_workflow_type), intent(in) :: job

    if (job%model_params%istate < 1) then
        error stop 'Invalid library workflow: istate < 1'
    end if

    if (job%model_params%n_candidates < 1) then
        error stop 'Invalid library workflow: n_candidates < 1'
    end if

    if (job%model_params%method_window /= 1 .and. &
        job%model_params%method_window /= 2) then
        error stop 'Invalid library workflow: method_window'
    end if

    select case (job%model_params%reference_method_index)
    case (1)
        if (trim(job%model_params%reference_method) /= 'rmsd') then
            error stop 'Invalid library workflow: reference method mismatch'
        end if
    case (2)
        if (trim(job%model_params%reference_method) /= 'cv_orthogonality') then
            error stop 'Invalid library workflow: reference method mismatch'
        end if
    case (3)
        if (trim(job%model_params%reference_method) /= 'fixed_frame') then
            error stop 'Invalid library workflow: reference method mismatch'
        end if
    case default
        error stop 'Invalid library workflow: reference_method_index'
    end select
end subroutine validate_library_workflow


subroutine write_library_workflow(job, filename)
    type(library_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename

    integer :: unit
    integer :: ios

    call validate_library_workflow(job)

    open(newunit=unit, file=trim(filename), &
         status='replace', action='write', iostat=ios)

    if (ios /= 0) then
        write(error_unit, '(A,A)') 'Cannot open workflow file: ', trim(filename)
        error stop
    end if

    call write_library_workflow_unit(job, unit)
    close(unit, iostat=ios)
end subroutine write_library_workflow


subroutine print_library_workflow(job)
    type(library_workflow_type), intent(in) :: job
    call write_library_workflow_unit(job, output_unit)
end subroutine print_library_workflow


subroutine write_library_workflow_unit(job, unit)
    type(library_workflow_type), intent(in) :: job
    integer, intent(in) :: unit

    write(unit, '(A)') '# Library workflow configuration'
    write(unit, '(A)') 'workflow_type = library'
    write(unit, '(A)') 'format_version = 1'
    write(unit, '(A)') ''

    write(unit, '(A)') '# Model parameters'
    write(unit, '(A,I0)') 'istate = ', job%model_params%istate
    write(unit, '(A,A)') 'reference_method = ', trim(job%model_params%reference_method)
    write(unit, '(A,I0)') 'reference_method_index = ', job%model_params%reference_method_index
    write(unit, '(A,I0)') 'n_candidates = ', job%model_params%n_candidates
    write(unit, '(A,I0)') 'method_window = ', job%model_params%method_window

    write(unit, '(A)') ''
    write(unit, '(A)') '# Run parameters'
    call write_bool(unit, 'initialize_sidechain_maps', job%run_params%initialize_sidechain_maps)
    call write_bool(unit, 'compute_sidechain_frequencies', job%run_params%compute_sidechain_frequencies)
    call write_bool(unit, 'test_sidechain_maps', job%run_params%test_sidechain_maps)
    call write_bool(unit, 'initialize_domain_maps', job%run_params%initialize_domain_maps)
    call write_bool(unit, 'initialize_user_cvs', job%run_params%initialize_user_cvs)
    call write_bool(unit, 'search_reference_structures', job%run_params%search_reference_structures)

    write(unit, '(A)') ''
    write(unit, '(A)') '# I/O parameters'
    call write_bool(unit, 'read_sidechain_library', job%io_params%read_sidechain_library)
    call write_bool(unit, 'read_domain_library', job%io_params%read_domain_library)
    call write_bool(unit, 'read_user_cv_library', job%io_params%read_user_cv_library)
    call write_bool(unit, 'write_sidechain_library', job%io_params%write_sidechain_library)
    call write_bool(unit, 'write_domain_library', job%io_params%write_domain_library)
    call write_bool(unit, 'write_user_cv_library', job%io_params%write_user_cv_library)

    call write_string(unit, 'sidechain_library_file', job%io_params%sidechain_library_file)
    call write_string(unit, 'domain_library_file', job%io_params%domain_library_file)
    call write_string(unit, 'user_cv_library_file', job%io_params%user_cv_library_file)

    write(unit, '(A,I0)') 'reference_window_istate = ', job%reference_window%istate
    write(unit, '(A,I0)') 'reference_window_irep = ', job%reference_window%irep
    write(unit, '(A,I0)') 'reference_window_start_time_step = ', job%reference_window%start_time_step
    write(unit, '(A,I0)') 'reference_window_end_time_step = ', job%reference_window%end_time_step
end subroutine write_library_workflow_unit


subroutine write_bool(unit, key, value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: key
    logical, intent(in) :: value

    if (value) then
        write(unit, '(A," = .true.")') trim(key)
    else
        write(unit, '(A," = .false.")') trim(key)
    end if
end subroutine write_bool


subroutine write_string(unit, key, value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value

    write(unit, '(A," = ",A)') trim(key), trim(value)
end subroutine write_string


subroutine read_library_workflow(job, filename)
    type(library_workflow_type), intent(inout) :: job
    character(len=*), intent(in) :: filename

    character(len=2048) :: line
    character(len=128)  :: key
    character(len=1900) :: value

    integer :: unit
    integer :: ios
    integer :: line_number
    integer :: eq_pos
    integer :: format_version

    logical :: seen_type
    logical :: seen_version

    format_version = -1
    seen_type = .false.
    seen_version = .false.

    open(newunit=unit, file=trim(filename), &
         status='old', action='read', iostat=ios)

    if (ios /= 0) then
        write(error_unit, '(A,A)') 'Cannot open workflow file: ', trim(filename)
        error stop
    end if

    line_number = 0

    do
        read(unit, '(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) then
            close(unit)
            error stop 'Error reading workflow file'
        end if

        line_number = line_number + 1
        line = adjustl(line)

        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle

        eq_pos = index(line, '=')
        if (eq_pos == 0) then
            call config_error(line_number, 'Missing "=" in configuration line')
        end if

        key = adjustl(line(:eq_pos-1))
        key = trim(key)

        value = adjustl(line(eq_pos+1:))
        value = trim(value)

        if (len_trim(key) == 0) then
            call config_error(line_number, 'Empty parameter name')
        end if

        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'library') then
                call config_error(line_number, 'Wrong workflow_type')
            end if
            seen_type = .true.

        case ('format_version')
            read(value, *, iostat=ios) format_version
            call check_value_read(ios, line_number, key)
            if (format_version /= 1) then
                call config_error(line_number, 'Unsupported format_version')
            end if
            seen_version = .true.

        case ('istate')
            read(value, *, iostat=ios) job%model_params%istate
            call check_value_read(ios, line_number, key)

        case ('reference_method')
            job%model_params%reference_method = trim(value)

        case ('reference_method_index')
            read(value, *, iostat=ios) job%model_params%reference_method_index
            call check_value_read(ios, line_number, key)

        case ('n_candidates')
            read(value, *, iostat=ios) job%model_params%n_candidates
            call check_value_read(ios, line_number, key)

        case ('method_window')
            read(value, *, iostat=ios) job%model_params%method_window
            call check_value_read(ios, line_number, key)

        case ('initialize_sidechain_maps')
            read(value, *, iostat=ios) job%run_params%initialize_sidechain_maps
            call check_value_read(ios, line_number, key)

        case ('compute_sidechain_frequencies')
            read(value, *, iostat=ios) job%run_params%compute_sidechain_frequencies
            call check_value_read(ios, line_number, key)

        case ('test_sidechain_maps')
            read(value, *, iostat=ios) job%run_params%test_sidechain_maps
            call check_value_read(ios, line_number, key)

        case ('initialize_domain_maps')
            read(value, *, iostat=ios) job%run_params%initialize_domain_maps
            call check_value_read(ios, line_number, key)

        case ('initialize_user_cvs')
            read(value, *, iostat=ios) job%run_params%initialize_user_cvs
            call check_value_read(ios, line_number, key)

        case ('search_reference_structures')
            read(value, *, iostat=ios) job%run_params%search_reference_structures
            call check_value_read(ios, line_number, key)

        case ('read_sidechain_library')
            read(value, *, iostat=ios) job%io_params%read_sidechain_library
            call check_value_read(ios, line_number, key)

        case ('read_domain_library')
            read(value, *, iostat=ios) job%io_params%read_domain_library
            call check_value_read(ios, line_number, key)

        case ('read_user_cv_library')
            read(value, *, iostat=ios) job%io_params%read_user_cv_library
            call check_value_read(ios, line_number, key)

        case ('write_sidechain_library')
            read(value, *, iostat=ios) job%io_params%write_sidechain_library
            call check_value_read(ios, line_number, key)

        case ('write_domain_library')
            read(value, *, iostat=ios) job%io_params%write_domain_library
            call check_value_read(ios, line_number, key)

        case ('write_user_cv_library')
            read(value, *, iostat=ios) job%io_params%write_user_cv_library
            call check_value_read(ios, line_number, key)

        case ('sidechain_library_file')
            job%io_params%sidechain_library_file = trim(value)

        case ('domain_library_file')
            job%io_params%domain_library_file = trim(value)

        case ('user_cv_library_file')
            job%io_params%user_cv_library_file = trim(value)

        case ('reference_window_istate')
            read(value,*,iostat=ios) job%reference_window%istate
            call check_value_read(ios, line_number, key)

        case ('reference_window_irep')
            read(value,*,iostat=ios) job%reference_window%irep
            call check_value_read(ios, line_number, key)

        case ('reference_window_start_time_step')
            read(value,*,iostat=ios) job%reference_window%start_time_step
            call check_value_read(ios, line_number, key)

        case ('reference_window_end_time_step')
            read(value,*,iostat=ios) job%reference_window%end_time_step
            call check_value_read(ios, line_number, key)

        case default
            call config_error(line_number, 'Unknown parameter: '//trim(key))
        end select
    end do

    close(unit)

    if (.not. seen_type) error stop 'Workflow file has no workflow_type'
    if (.not. seen_version) error stop 'Workflow file has no format_version'
end subroutine read_library_workflow


subroutine check_value_read(ios, line_number, key)
    integer, intent(in) :: ios
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: key

    if (ios /= 0) then
        write(error_unit, '(A,I0,A,A)') &
            'Invalid value at line ', line_number, &
            ' for parameter ', trim(key)
        error stop
    end if
end subroutine check_value_read


subroutine config_error(line_number, message)
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: message

    write(error_unit, '(A,I0,2A)') &
        'Workflow configuration error at line ', &
        line_number, ': ', trim(message)
    error stop
end subroutine config_error

end module cv_lib_workflow_io_module