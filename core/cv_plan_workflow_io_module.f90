module cv_plan_workflow_io_module

    use workflow_module, only: cv_plan_workflow_type, sidechain_lib_config_type, domain_lib_config_type, user_lib_config_type
    use iso_fortran_env, only: output_unit, error_unit

    implicit none
    private

    public :: init_cv_plan_workflow
    public :: read_cv_plan_workflow
    public :: write_cv_plan_workflow
    public :: print_cv_plan_workflow
    public :: validate_cv_plan_workflow
    public :: configure_cv_plan_workflow

contains

subroutine configure_cv_plan_workflow(job)
    implicit none

    type(cv_plan_workflow_type), intent(out) :: job

    character(len=1024) :: filename
    character(len=32)   :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    CV plan workflow configuration'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Create new CV plan workflow'
    write(*,*) '2 = Load and edit existing CV plan workflow'
    write(*,*) '3 = Load existing CV plan workflow unchanged'
    write(*,*) '4 = Skip CV plan workflow configuration'
    write(*,*)
    write(*,'(A)', advance='no') 'Select option: '

    read(*,'(A)') input
    read(input,*,iostat=ios) choice

    if (ios /= 0) error stop 'Invalid CV plan workflow option'

    select case (choice)

    case (1)
        job = cv_plan_workflow_type()
        ! Standard-Arraygrößen initialisieren (1 Datei als Default)
        allocate(job%sidechain_libs(1))
        job%sidechain_libs(1)%filename = 'sidechain_lib.txt'
        job%sidechain_libs(1)%selection_mode = 1

        allocate(job%domain_libs(1))
        job%domain_libs(1)%filename = 'domain_lib.txt'
        job%domain_libs(1)%selection_mode = 1

        allocate(job%user_libs(1))
        job%user_libs(1)%filename = 'user_lib.txt'
        job%user_libs(1)%selection_mode = 1

        call edit_cv_plan_workflow(job)

    case (2)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) error stop 'No configuration filename specified'

        call read_cv_plan_workflow(job, trim(filename))
        call edit_cv_plan_workflow(job)
        job%initialized = .true.

    case (3)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) error stop 'No configuration filename specified'

        call read_cv_plan_workflow(job, trim(filename))
        call validate_cv_plan_workflow(job)

        write(*,*)
        write(*,*) 'CV plan workflow loaded:'
        call print_cv_plan_workflow(job)
        job%initialized = .true.

    case (4)
        job = cv_plan_workflow_type()
        write(*,*) 
        write(*,*) 'CV plan workflow configuration skipped.'

    case default
        error stop 'Invalid CV plan workflow option'

    end select

end subroutine configure_cv_plan_workflow


subroutine edit_cv_plan_workflow(job)
    use sidechains 
    use domains_module
    use user_defined_cv_module
    implicit none

    type(cv_plan_workflow_type), intent(inout) :: job

    character(len=1024) :: input, filename
    character(len=32)   :: cname
    integer :: ios, i, k, n, n_libs, sc_mode, rid
    logical :: save_config, has_user_lib

    ! Temporäre Bibliotheks-Objekte zum Auslesen der Namen/Strukturen
    type(sidechain_library_type) :: sc_lib
    type(domain_library_type)    :: dom_lib
    type(user_cv_library_type)   :: user_lib

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    Edit CV plan workflow'
    write(*,*) '======================================'
    write(*,*) 'Press Enter to keep the current value.'
    write(*,*)

    ! --------------------------------------------------
    ! 1. Sidechain Library Files Setup
    ! --------------------------------------------------
    n_libs = 1
    if (allocated(job%sidechain_libs)) n_libs = size(job%sidechain_libs)

    write(*,'(A,I0,A)', advance='no') 'Number of sidechain libraries [', n_libs, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) n
        if (ios /= 0 .or. n < 0) error stop 'Invalid number of sidechain libraries'
        n_libs = n
    end if

    if (.not. allocated(job%sidechain_libs) .or. size(job%sidechain_libs) /= n_libs) then
        if (allocated(job%sidechain_libs)) deallocate(job%sidechain_libs)
        if (n_libs > 0) then
            allocate(job%sidechain_libs(n_libs))
            job%sidechain_libs(1)%filename = 'sidechain_lib.txt'
            job%sidechain_libs(1)%selection_mode = 1
        end if
    end if

    do i = 1, n_libs
        write(*,'(A,I0,A,A,A)', advance='no') '   Sidechain library file #', i, ' [', trim(job%sidechain_libs(i)%filename), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            job%sidechain_libs(i)%filename = trim(input)
        end if
    end do

    ! --------------------------------------------------
    ! 2. Domain Library Files Setup
    ! --------------------------------------------------
    n_libs = 1
    if (allocated(job%domain_libs)) n_libs = size(job%domain_libs)

    write(*,*)
    write(*,'(A,I0,A)', advance='no') 'Number of domain libraries [', n_libs, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) n
        if (ios /= 0 .or. n < 0) error stop 'Invalid number of domain libraries'
        n_libs = n
    end if

    if (.not. allocated(job%domain_libs) .or. size(job%domain_libs) /= n_libs) then
        if (allocated(job%domain_libs)) deallocate(job%domain_libs)
        if (n_libs > 0) then
            allocate(job%domain_libs(n_libs))
            job%domain_libs(1)%filename = 'domain_lib.txt'
            job%domain_libs(1)%selection_mode = 1
        end if
    end if

    do i = 1, n_libs
        write(*,'(A,I0,A,A,A)', advance='no') '   Domain library file #', i, ' [', trim(job%domain_libs(i)%filename), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            job%domain_libs(i)%filename = trim(input)
        end if
    end do

    ! --------------------------------------------------
    ! 3. User-defined Library Files Setup
    ! --------------------------------------------------
    n_libs = 1
    if (allocated(job%user_libs)) n_libs = size(job%user_libs)

    write(*,*)
    write(*,'(A,I0,A)', advance='no') 'Number of user-defined libraries [', n_libs, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) n
        if (ios /= 0 .or. n < 0) error stop 'Invalid number of user libraries'
        n_libs = n
    end if

    if (.not. allocated(job%user_libs) .or. size(job%user_libs) /= n_libs) then
        if (allocated(job%user_libs)) deallocate(job%user_libs)
        if (n_libs > 0) then
            allocate(job%user_libs(n_libs))
            job%user_libs(1)%filename = 'user_lib.txt'
            job%user_libs(1)%selection_mode = 1
        end if
    end if

    do i = 1, n_libs
        write(*,'(A,I0,A,A,A)', advance='no') '   User library file #', i, ' [', trim(job%user_libs(i)%filename), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            job%user_libs(i)%filename = trim(input)
        end if
    end do

    ! ==========================================================
    ! INTERAKTIVES FILTERN DER CVs (Pro Bibliothek)
    ! ==========================================================
    write(*,*)
    write(*,*) '--------------------------------------------------'
    write(*,*) ' CV Filtering & Activation Setup (Workflow-Level):'
    write(*,*) '--------------------------------------------------'

    ! --- 1. SIDECHAIN CV FILTERS ---
    write(*,*)
    write(*,'(A,L1,A)', advance='no') 'Include sidechain CVs [', job%include_sidechain_cvs, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%include_sidechain_cvs
        if (ios /= 0) error stop 'Invalid logical value'
    end if

    if (job%include_sidechain_cvs .and. allocated(job%sidechain_libs)) then
        do i = 1, size(job%sidechain_libs)
            filename = job%sidechain_libs(i)%filename
            
            
            ios=0	
            call read_sidechain_library(trim(filename),sc_lib)
            if (ios /= 0) then
                write(*,*) 'Warning: Could not read sidechain library file "', trim(filename), '". Skipping filter.'
                cycle
            end if

            write(*,*)
            write(*,'(A,I0,A,A,A)') '--- Sidechain CV Filters for Library #', i, ' (', trim(filename), ') ---'
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
            job%sidechain_libs(i)%selection_mode = sc_mode

            if (allocated(job%sidechain_libs(i)%filters)) deallocate(job%sidechain_libs(i)%filters)
            allocate(job%sidechain_libs(i)%filters(sc_lib%ncv))
            job%sidechain_libs(i)%filters(:) = .true.

            select case (sc_mode)
            case (1)
                job%sidechain_libs(i)%filters(:) = .true.
                write(*,*) '-> All sidechains enabled for this library.'

            case (2)
                job%sidechain_libs(i)%filters(:) = .false.
                write(*,*) '-> All sidechains disabled for this library.'

            case (3)
                write(*,*)
                write(*,*) '--- Configure Sidechains by Residue Type ---'
                do k = 1, sc_lib%ncv
                    cname = sc_lib%residue(k)%resname
                    write(*,'(A,I4,A,A3,A)', advance='no') '  [Type ', k, '] Include all ', cname, ' sidechains? [y] (y/n): '
                    read(*,'(A)') input
                    if (len_trim(input) > 0) then
                        input = trim(adjustl(input))
                        if (input == 'n' .or. input == 'N' .or. input == 'no') then
                            job%sidechain_libs(i)%filters(k) = .false.
                        end if
                    end if
                end do

            case (4)
                write(*,*)
                write(*,*) '--- Detailed Individual Sidechain Configuration ---'
                do k = 1, sc_lib%ncv
                    if (allocated(sc_lib%residue(k)%res_ids)) then
                        do n = 1, size(sc_lib%residue(k)%res_ids)
                            rid = sc_lib%residue(k)%res_ids(n)
                            cname = sc_lib%residue(k)%resname

                            write(*,'(A,I4,A,I6,A,A3,A)', advance='no') &
                                '  [Type ', k, ', Residue ID ', rid, '] Name: ', cname, &
                                ' - Include? [y] (y/n): '
                            
                            read(*,'(A)') input
                            if (len_trim(input) > 0) then
                                input = trim(adjustl(input))
                                if (input == 'n' .or. input == 'N' .or. input == 'no') then
                                    ! Falls detaillierte Instanzfilter im Job gewünscht sind:
                                    ! job%sidechain_libs(i)%instance_filters(k, n) = .false.
                                end if
                            end if
                        end do
                    end if
                end do
            end select

            call clear_sidechain_library(sc_lib)
        end do
    end if

    ! --- 2. DOMAIN CV FILTERS ---
    write(*,*)
    write(*,'(A,L1,A)', advance='no') 'Include domain CVs [', job%include_domain_cvs, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%include_domain_cvs
        if (ios /= 0) error stop 'Invalid logical value'
    end if

    if (job%include_domain_cvs .and. allocated(job%domain_libs)) then
        do i = 1, size(job%domain_libs)
            filename = job%domain_libs(i)%filename

            call read_domain_library(trim(filename),dom_lib)
            if (ios /= 0) then
                write(*,*) 'Warning: Could not read domain library file "', trim(filename), '". Skipping filter.'
                cycle
            end if

            if (allocated(job%domain_libs(i)%filters)) deallocate(job%domain_libs(i)%filters)
            allocate(job%domain_libs(i)%filters(dom_lib%ndo))
            job%domain_libs(i)%filters(:) = .true.

            write(*,*)
            write(*,'(A,I0,A,A,A)') '--- Domain CV Filters for Library #', i, ' (', trim(filename), ') ---'
            write(*,'(A)', advance='no') 'Include all domains automatically without individual prompts? (y/n) [y]: '
            read(*,'(A)') input
            if (len_trim(input) == 0) input = 'y'
            input = trim(adjustl(input))

            if (input /= 'y' .and. input /= 'Y' .and. input /= 'yes') then
                job%domain_libs(i)%selection_mode = 2
                do k = 1, dom_lib%ndo
                    write(*,'(A,I4,A,A,A)', advance='no') &
                        '   [', k, '] Domain: ', trim(dom_lib%domain(k)%domain_name), &
                        ' - Include in workflow? [y] (y/n): '
                    
                    read(*,'(A)') input
                    if (len_trim(input) > 0) then
                        input = trim(adjustl(input))
                        if (input == 'n' .or. input == 'N' .or. input == 'no') then
                            job%domain_libs(i)%filters(k) = .false.
                        end if
                    end if
                end do
            else
                job%domain_libs(i)%selection_mode = 1
                write(*,*) '-> All domains included by default.'
            end if

            call clear_domain_library(dom_lib)
        end do
    end if

    ! --- 3. USER-DEFINED CV FILTERS ---
    write(*,*)
    write(*,'(A,L1,A)', advance='no') 'Include user-defined CVs [', job%include_user_cvs, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%include_user_cvs
        if (ios /= 0) error stop 'Invalid logical value'
    end if

    if (job%include_user_cvs .and. allocated(job%user_libs)) then
        do i = 1, size(job%user_libs)
            filename = job%user_libs(i)%filename

            call read_user_cv_library( trim(filename),user_lib)
            if (ios /= 0) then
                write(*,*) 'Warning: Could not read user library file "', trim(filename), '". Skipping filter.'
                cycle
            end if

            if (user_lib%n_cvs > 0) then
                if (allocated(job%user_libs(i)%filters)) deallocate(job%user_libs(i)%filters)
                allocate(job%user_libs(i)%filters(user_lib%n_cvs))
                job%user_libs(i)%filters(:) = .true.

                write(*,*)
                write(*,'(A,I0,A,A,A)') '--- User-Defined CV Filters for Library #', i, ' (', trim(filename), ') ---'
                do k = 1, user_lib%n_cvs
                    write(*,'(A,I4,A,A,A)', advance='no') &
                        '   [', k, '] CV Name: ', trim(user_lib%cvs(k)%cv_name), &
                        ' - Include in workflow? [y] (y/n): '
                    
                    read(*,'(A)') input
                    if (len_trim(input) > 0) then
                        input = trim(adjustl(input))
                        if (input == 'n' .or. input == 'N' .or. input == 'no') then
                            job%user_libs(i)%filters(k) = .false.
                        end if
                    end if
                end do
            end if

            call clear_user_cv_library(user_lib)
        end do
    end if

    ! --------------------------------------------------
    ! Validate and Save Configuration
    ! --------------------------------------------------
    call validate_cv_plan_workflow(job)

    write(*,*)
    write(*,*) 'Resulting CV plan workflow configuration:'
    call print_cv_plan_workflow(job)

    write(*,*)
    write(*,'(A)', advance='no') 'Save configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'

    save_config = (trim(input) == 'y' .or. trim(input) == 'Y' .or. trim(input) == 'yes')

    if (save_config) then
        write(*,'(A,A,A)', advance='no') 'Configuration filename [cv_plan_workflow.cfg]: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) filename = 'cv_plan_workflow.cfg'

        call write_cv_plan_workflow(job, trim(filename))
        write(*,'(3A)') 'Configuration saved to "', trim(filename), '"'
    end if

end subroutine edit_cv_plan_workflow


subroutine init_cv_plan_workflow(job, mode, filename)
    type(cv_plan_workflow_type), intent(out) :: job
    character(len=*), intent(in) :: mode
    character(len=*), intent(in), optional :: filename

    job = cv_plan_workflow_type()
    
    allocate(job%sidechain_libs(1))
    job%sidechain_libs(1)%filename = 'sidechain_lib.txt'
    job%sidechain_libs(1)%selection_mode = 1

    allocate(job%domain_libs(1))
    job%domain_libs(1)%filename = 'domain_lib.txt'
    job%domain_libs(1)%selection_mode = 1

    allocate(job%user_libs(1))
    job%user_libs(1)%filename = 'user_lib.txt'
    job%user_libs(1)%selection_mode = 1

    select case (trim(mode))
    case ('manual')
        ! Defaults behalten
    case ('file')
        if (.not. present(filename)) error stop 'init_cv_plan_workflow: filename missing'
        call read_cv_plan_workflow(job, filename)
    case default
        write(error_unit, '(A,A)') 'init_cv_plan_workflow: unknown mode: ', trim(mode)
        error stop
    end select

    call validate_cv_plan_workflow(job)
end subroutine init_cv_plan_workflow


subroutine validate_cv_plan_workflow(job)
    implicit none
    type(cv_plan_workflow_type), intent(in) :: job
    ! Zusätzliche Konsistenzprüfungen falls nötig
end subroutine validate_cv_plan_workflow


subroutine write_cv_plan_workflow(job, filename)
    type(cv_plan_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_cv_plan_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'Cannot open CV plan workflow file: ', trim(filename)
        error stop
    end if

    call write_cv_plan_workflow_unit(job, unit)
    close(unit, iostat=ios)
end subroutine write_cv_plan_workflow


subroutine print_cv_plan_workflow(job)
    type(cv_plan_workflow_type), intent(in) :: job
    call write_cv_plan_workflow_unit(job, output_unit)
end subroutine print_cv_plan_workflow


subroutine write_cv_plan_workflow_unit(job, unit)
    type(cv_plan_workflow_type), intent(in) :: job
    integer, intent(in) :: unit
    integer :: i

    write(unit, '(A)') '# CV plan workflow configuration'
    write(unit, '(A)') 'workflow_type = cv_plan'
    write(unit, '(A)') 'format_version = 1'
    write(unit, '(A)') ''

    write(unit, '(A)') '# Run parameters'
    call write_bool(unit, 'include_sidechain_cvs', job%include_sidechain_cvs)
    call write_bool(unit, 'include_domain_cvs', job%include_domain_cvs)
    call write_bool(unit, 'include_user_cvs', job%include_user_cvs)

    write(unit, '(A)') ''
    write(unit, '(A)') '# Library configurations'
    
    ! Sidechain Bibliotheken und Modi schreiben
    if (allocated(job%sidechain_libs)) then
        write(unit, '(A,I0)') 'n_sidechain_libs = ', size(job%sidechain_libs)
        do i = 1, size(job%sidechain_libs)
            write(unit, '(A,I0,A,A)') 'sidechain_library_file_', i, ' = ', trim(job%sidechain_libs(i)%filename)
            write(unit, '(A,I0,A,I0)') 'sidechain_selection_mode_', i, ' = ', job%sidechain_libs(i)%selection_mode
        end do
    end if

    write(unit, '(A)') ''
    ! Domain Bibliotheken und Modi schreiben
    if (allocated(job%domain_libs)) then
        write(unit, '(A,I0)') 'n_domain_libs = ', size(job%domain_libs)
        do i = 1, size(job%domain_libs)
            write(unit, '(A,I0,A,A)') 'domain_library_file_', i, ' = ', trim(job%domain_libs(i)%filename)
            write(unit, '(A,I0,A,I0)') 'domain_selection_mode_', i, ' = ', job%domain_libs(i)%selection_mode
        end do
    end if

    write(unit, '(A)') ''
    ! User Bibliotheken und Modi schreiben
    if (allocated(job%user_libs)) then
        write(unit, '(A,I0)') 'n_user_cv_libs = ', size(job%user_libs)
        do i = 1, size(job%user_libs)
            write(unit, '(A,I0,A,A)') 'user_library_file_', i, ' = ', trim(job%user_libs(i)%filename)
            write(unit, '(A,I0,A,I0)') 'user_selection_mode_', i, ' = ', job%user_libs(i)%selection_mode
        end do
    end if

end subroutine write_cv_plan_workflow_unit


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


subroutine read_cv_plan_workflow(job, filename)
    type(cv_plan_workflow_type), intent(inout) :: job
    character(len=*), intent(in) :: filename

    character(len=2048) :: line
    character(len=128)  :: key
    character(len=1900) :: value

    integer :: unit, ios, line_number, eq_pos, format_version, idx, n_libs
    logical :: seen_type, seen_version

    job = cv_plan_workflow_type()
    format_version = -1
    seen_type = .false.
    seen_version = .false.

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'Cannot open CV plan workflow file: ', trim(filename)
        error stop
    end if

    line_number = 0

    do
        read(unit, '(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) then
            close(unit)
            error stop 'Error reading CV plan workflow file'
        end if

        line_number = line_number + 1
        line = adjustl(line)

        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle

        eq_pos = index(line, '=')
        if (eq_pos == 0) call config_error(line_number, 'Missing "=" in configuration line')

        key = trim(adjustl(line(:eq_pos-1)))
        value = trim(adjustl(line(eq_pos+1:)))

        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'cv_plan') call config_error(line_number, 'Wrong workflow_type')
            seen_type = .true.

        case ('format_version')
            read(value, *, iostat=ios) format_version
            call check_value_read(ios, line_number, key)
            if (format_version /= 1) call config_error(line_number, 'Unsupported format_version')
            seen_version = .true.

        case ('include_sidechain_cvs')
            read(value, *, iostat=ios) job%include_sidechain_cvs
            call check_value_read(ios, line_number, key)

        case ('include_domain_cvs')
            read(value, *, iostat=ios) job%include_domain_cvs
            call check_value_read(ios, line_number, key)

        case ('include_user_cvs')
            read(value, *, iostat=ios) job%include_user_cvs
            call check_value_read(ios, line_number, key)

        case ('n_sidechain_libs')
            read(value, *, iostat=ios) n_libs
            call check_value_read(ios, line_number, key)
            if (allocated(job%sidechain_libs)) deallocate(job%sidechain_libs)
            if (n_libs >  0) allocate(job%sidechain_libs(n_libs))

        case ('n_domain_libs')
            read(value, *, iostat=ios) n_libs
            call check_value_read(ios, line_number, key)
            if (allocated(job%domain_libs)) deallocate(job%domain_libs)
            if (n_libs >  0) allocate(job%domain_libs(n_libs))

        case ('n_user_cv_libs')
            read(value, *, iostat=ios) n_libs
            call check_value_read(ios, line_number, key)
            if (allocated(job%user_libs)) deallocate(job%user_libs)
            if (n_libs >  0) allocate(job%user_libs(n_libs))

        case default
            if (index(trim(key), 'sidechain_library_file_') == 1) then
                read(key(24:), *, iostat=ios) idx
                if (ios == 0 .and. idx > 0 .and. allocated(job%sidechain_libs)) then
                    if (idx <= size(job%sidechain_libs)) job%sidechain_libs(idx)%filename = trim(value)
                end if
            else if (index(trim(key), 'sidechain_selection_mode_') == 1) then
                read(key(26:), *, iostat=ios) idx
                if (ios == 0 .and. idx > 0 .and. allocated(job%sidechain_libs)) then
                    if (idx <= size(job%sidechain_libs)) then
                        read(value, *, iostat=ios) job%sidechain_libs(idx)%selection_mode
                        call check_value_read(ios, line_number, key)
                    end if
                end if
            else if (index(trim(key), 'domain_library_file_') == 1) then
                read(key(21:), *, iostat=ios) idx
                if (ios == 0 .and. idx > 0 .and. allocated(job%domain_libs)) then
                    if (idx <= size(job%domain_libs)) job%domain_libs(idx)%filename = trim(value)
                end if
            else if (index(trim(key), 'domain_selection_mode_') == 1) then
                read(key(23:), *, iostat=ios) idx
                if (ios == 0 .and. idx > 0 .and. allocated(job%domain_libs)) then
                    if (idx <= size(job%domain_libs)) then
                        read(value, *, iostat=ios) job%domain_libs(idx)%selection_mode
                        call check_value_read(ios, line_number, key)
                    end if
                end if
            else if (index(trim(key), 'user_library_file_') == 1) then
                read(key(19:), *, iostat=ios) idx
                if (ios == 0 .and. idx > 0 .and. allocated(job%user_libs)) then
                    if (idx <= size(job%user_libs)) job%user_libs(idx)%filename = trim(value)
                end if
            else if (index(trim(key), 'user_selection_mode_') == 1) then
                read(key(21:), *, iostat=ios) idx
                if (ios == 0 .and. idx > 0 .and. allocated(job%user_libs)) then
                    if (idx <= size(job%user_libs)) then
                        read(value, *, iostat=ios) job%user_libs(idx)%selection_mode
                        call check_value_read(ios, line_number, key)
                    end if
                end if
            else
                call config_error(line_number, 'Unknown parameter: '//trim(key))
            end if
        end select
    end do

    close(unit)
    if (.not. seen_type) error stop 'CV plan workflow file has no workflow_type'
    if (.not. seen_version) error stop 'CV plan workflow file has no format_version'
end subroutine read_cv_plan_workflow


subroutine check_value_read(ios, line_number, key)
    integer, intent(in) :: ios, line_number
    character(len=*), intent(in) :: key
    if (ios /= 0) then
        write(error_unit, '(A,I0,A,A)') 'Invalid value at line ', line_number, ' for parameter ', trim(key)
        error stop
    end if
end subroutine check_value_read


subroutine config_error(line_number, message)
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: message
    write(error_unit, '(A,I0,2A)') 'CV plan workflow error at line ', line_number, ': ', trim(message)
    error stop
end subroutine config_error

end module cv_plan_workflow_io_module