module euclidean_lifting_workflow_io_module

    use workflow_module, only: euclidean_lifting_workflow_type ! Passe den Modulnamen ggf. an
    use iso_fortran_env, only: output_unit, error_unit

    implicit none
    private

    public :: init_euclidean_lifting_workflow
    public :: read_euclidean_lifting_workflow
    public :: write_euclidean_lifting_workflow
    public :: print_euclidean_lifting_workflow
    public :: validate_euclidean_lifting_workflow
    public :: configure_euclidean_lifting_workflow

contains

! ============================================================
! KONFIGURATION (Interaktiv)
! ============================================================

subroutine configure_euclidean_lifting_workflow(job)
    implicit none

    type(euclidean_lifting_workflow_type), intent(out) :: job

    character(len=1024) :: filename
    character(len=32)   :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Euclidean lifting workflow configuration'
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
        job = euclidean_lifting_workflow_type()

        call edit_euclidean_lifting_workflow(job)

    case (2)

        ! Vorhandene Konfiguration laden und bearbeiten
        write(*,'(A)', advance='no') &
            'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_euclidean_lifting_workflow(job, trim(filename))
        call edit_euclidean_lifting_workflow(job)

    case (3)

        ! Vorhandene Konfiguration unverändert übernehmen
        write(*,'(A)', advance='no') &
            'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_euclidean_lifting_workflow(job, trim(filename))
        call validate_euclidean_lifting_workflow(job)

        write(*,*)
        write(*,*) 'Workflow loaded:'
        call print_euclidean_lifting_workflow(job)

    case (4)

        ! Workflow-Konfiguration überspringen
        job = euclidean_lifting_workflow_type()

        write(*,*)
        write(*,*) 'Workflow configuration skipped.'

    case default

        error stop 'Invalid workflow option'

    end select

end subroutine configure_lifting_workflow ! Hinweis: intern konsistent mit Modul-Namen unten angepasst


! ============================================================
! EDITIEREN DER PARAMETER
! ============================================================

subroutine edit_euclidean_lifting_workflow(job)
    implicit none

    type(euclidean_lifting_workflow_type), intent(inout) :: job

    character(len=1024) :: input
    character(len=1024) :: filename
    integer :: ios
    logical :: save_config

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Edit Euclidean lifting workflow'
    write(*,*) '======================================'
    write(*,*) 'Press Enter to keep the current value.'
    write(*,*)

    ! --------------------------------------------------
    ! Model parameters
    ! --------------------------------------------------

    write(*,'(A,I0,A)', advance='no') &
        'Basis index [', job%model_params%basis_index, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%model_params%basis_index
        if (ios /= 0) error stop 'Invalid basis index'
    end if

    ! --------------------------------------------------
    ! Run parameters (Frame Selection)
    ! --------------------------------------------------

    write(*,*)
    write(*,*) 'Reference frame selection:'

    write(*,'(A,I0,A)', advance='no') &
        'Reference state index [', job%run_params%ref_istate, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%ref_istate
        if (ios /= 0) error stop 'Invalid reference state index'
    end if

    write(*,'(A,I0,A)', advance='no') &
        'Reference replicate index [', job%run_params%ref_irep, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%ref_irep
        if (ios /= 0) error stop 'Invalid reference replicate index'
    end if

    write(*,'(A,I0,A)', advance='no') &
        'Reference time step index [', job%run_params%ref_itime, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%ref_itime
        if (ios /= 0) error stop 'Invalid reference time step index'
    end if

    ! --------------------------------------------------
    ! Run parameters (Lifting Plan)
    ! --------------------------------------------------

    write(*,*)
    write(*,*) 'Lifting plan parameters:'

    write(*,'(A,I0,A)', advance='no') &
        'Offset [', job%run_params%offset, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%offset
        if (ios /= 0) error stop 'Invalid offset'
    end if

    write(*,'(A,I0,A)', advance='no') &
        'Stride [', job%run_params%stride, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%stride
        if (ios /= 0) error stop 'Invalid stride'
    end if

    write(*,'(A,I0,A)', advance='no') &
        'Block size [', job%run_params%size_block, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%size_block
        if (ios /= 0) error stop 'Invalid block size'
    end if

    ! --------------------------------------------------
    ! Run parameters (Logical Flags)
    ! --------------------------------------------------

    write(*,*)
    write(*,*) 'Execution flags (.true. / .false.):'

    write(*,'(A,L1,A)', advance='no') &
        'Test periodic boundary conditions [', job%run_params%test_pbc, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%test_pbc
        if (ios /= 0) error stop 'Invalid logical value for test_pbc'
    end if

    write(*,'(A,L1,A)', advance='no') &
        'Execute lifting [', job%run_params%execute_lifting, ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%run_params%execute_lifting
        if (ios /= 0) error stop 'Invalid logical value for execute_lifting'
    end if

    ! --------------------------------------------------
    ! I/O parameters (Output Policy)
    ! --------------------------------------------------

    write(*,*)
    write(*,*) 'I/O parameters:'

    write(*,'(A,A,A)', advance='no') &
        'Output policy [', trim(job%io_params%output_policy), ']: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        job%io_params%output_policy = trim(input)
    end if

    ! Hinweis: file_ref_type spezifische Abfragen für Centroid/Lifted Trajectory
    ! können hier ergänzt werden, je nachdem wie dein Dateireferenz-Typ definiert ist.

    ! --------------------------------------------------
    ! Validate and display
    ! --------------------------------------------------

    call validate_euclidean_lifting_workflow(job)

    write(*,*)
    write(*,*) 'Resulting workflow configuration:'
    call print_euclidean_lifting_workflow(job)

    ! --------------------------------------------------
    ! Save configuration
    ! --------------------------------------------------

    write(*,*)
    write(*,'(A)', advance='no') 'Save configuration? (y/n): '
    read(*,'(A)') input

    save_config = trim(input) == 'y' .or. trim(input) == 'Y'

    if (save_config) then

        write(*,'(A)', advance='no') &
            'Configuration filename [euclidean_lifting_workflow.cfg]: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            filename = 'euclidean_lifting_workflow.cfg'
        end if

        call write_euclidean_lifting_workflow(job, trim(filename))

        write(*,'(3A)') 'Configuration saved to "', &
            trim(filename), '"'

    end if

end subroutine edit_euclidean_lifting_workflow


! ============================================================
! INITIALISIERUNG (Modus-gesteuert)
! ============================================================

subroutine init_euclidean_lifting_workflow(job, mode, filename)

    type(euclidean_lifting_workflow_type), intent(out) :: job
    character(len=*), intent(in) :: mode
    character(len=*), intent(in), optional :: filename

    job = euclidean_lifting_workflow_type()

    select case (trim(mode))

    case ('manual')
        ! Standardwerte aktiv

    case ('file')

        if (.not. present(filename)) then
            error stop 'init_euclidean_lifting_workflow: filename missing'
        end if

        call read_euclidean_lifting_workflow(job, filename)

    case default

        write(error_unit, '(A,A)') &
            'init_euclidean_lifting_workflow: unknown mode: ', trim(mode)
        error stop

    end select

    call validate_euclidean_lifting_workflow(job)

end subroutine init_euclidean_lifting_workflow


! ============================================================
    ! VALIDIERUNG
! ============================================================

subroutine validate_euclidean_lifting_workflow(job)
    implicit none

    type(euclidean_lifting_workflow_type), intent(in) :: job

    if (job%model_params%basis_index < 1) then
        error stop 'Invalid Euclidean lifting workflow: basis_index < 1'
    end if

    if (job%run_params%ref_istate < 1) then
        error stop 'Invalid Euclidean lifting workflow: ref_istate < 1'
    end if

    if (job%run_params%ref_irep < 1) then
        error stop 'Invalid Euclidean lifting workflow: ref_irep < 1'
    end if

    if (job%run_params%ref_itime < 1) then
        error stop 'Invalid Euclidean lifting workflow: ref_itime < 1'
    end if

    if (job%run_params%stride < 1) then
        error stop 'Invalid Euclidean lifting workflow: stride < 1'
    end if

    if (job%run_params%size_block < 1) then
        error stop 'Invalid Euclidean lifting workflow: size_block < 1'
    end if

end subroutine validate_euclidean_lifting_workflow


! ============================================================
! DATEI SCHREIBEN & AUSGEBEN
! ============================================================

subroutine write_euclidean_lifting_workflow(job, filename)

    type(euclidean_lifting_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename

    integer :: unit
    integer :: ios

    call validate_euclidean_lifting_workflow(job)

    open(newunit=unit, file=trim(filename), &
         status='replace', action='write', iostat=ios)

    if (ios /= 0) then
        write(error_unit, '(A,A)') &
            'Cannot open workflow file: ', trim(filename)
        error stop
    end if

    call write_euclidean_lifting_workflow_unit(job, unit)

    close(unit, iostat=ios)

    if (ios /= 0) then
        write(error_unit, '(A,A)') &
            'Error closing workflow file: ', trim(filename)
        error stop
    end if

end subroutine write_euclidean_lifting_workflow


subroutine print_library_workflow(job) ! (Name angepasst an Konsistenz, oder print_euclidean_lifting_workflow)
    type(euclidean_lifting_workflow_type), intent(in) :: job
    call write_euclidean_lifting_workflow_unit(job, output_unit)
end subroutine print_library_workflow

subroutine print_euclidean_lifting_workflow(job)
    type(euclidean_lifting_workflow_type), intent(in) :: job
    call write_euclidean_lifting_workflow_unit(job, output_unit)
end subroutine print_euclidean_lifting_workflow


subroutine write_euclidean_lifting_workflow_unit(job, unit)

    type(euclidean_lifting_workflow_type), intent(in) :: job
    integer, intent(in) :: unit

    write(unit, '(A)') '# Euclidean lifting workflow configuration'
    write(unit, '(A)') 'workflow_type = euclidean_lifting'
    write(unit, '(A)') 'format_version = 1'
    write(unit, '(A)') ''

    write(unit, '(A)') '# Model parameters'
    write(unit, '(A,I0)') 'basis_index = ', job%model_params%basis_index

    write(unit, '(A)') ''
    write(unit, '(A)') '# Run parameters'
    write(unit, '(A,I0)') 'ref_istate = ', job%run_params%ref_istate
    write(unit, '(A,I0)') 'ref_irep = ', job%run_params%ref_irep
    write(unit, '(A,I0)') 'ref_itime = ', job%run_params%ref_itime
    write(unit, '(A,I0)') 'offset = ', job%run_params%offset
    write(unit, '(A,I0)') 'stride = ', job%run_params%stride
    write(unit, '(A,I0)') 'size_block = ', job%run_params%size_block

    call write_bool(unit, 'test_pbc', job%run_params%test_pbc)
    call write_bool(unit, 'execute_lifting', job%run_params%execute_lifting)

    write(unit, '(A)') ''
    write(unit, '(A)') '# I/O parameters'
    call write_string(unit, 'output_policy', job%io_params%output_policy)

end subroutine write_euclidean_lifting_workflow_unit


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


! ============================================================
! DATEI LESEN
! ============================================================

subroutine read_euclidean_lifting_workflow(job, filename)

    type(euclidean_lifting_workflow_type), intent(inout) :: job
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

    job = euclidean_lifting_workflow_type()

    format_version = -1
    seen_type = .false.
    seen_version = .false.

    open(newunit=unit, file=trim(filename), &
         status='old', action='read', iostat=ios)

    if (ios /= 0) then
        write(error_unit, '(A,A)') &
            'Cannot open workflow file: ', trim(filename)
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

        key = trim(adjustl(line(:eq_pos-1)))
        value = trim(adjustl(line(eq_pos+1:)))

        if (len_trim(key) == 0) then
            call config_error(line_number, 'Empty parameter name')
        end if

        select case (trim(key))

        case ('workflow_type')
            if (trim(value) /= 'euclidean_lifting') then
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

        ! Model parameters
        case ('basis_index')
            read(value, *, iostat=ios) job%model_params%basis_index
            call check_value_read(ios, line_number, key)

        ! Run parameters
        case ('ref_istate')
            read(value, *, iostat=ios) job%run_params%ref_istate
            call check_value_read(ios, line_number, key)

        case ('ref_irep')
            read(value, *, iostat=ios) job%run_params%ref_irep
            call check_value_read(ios, line_number, key)

        case ('ref_itime')
            read(value, *, iostat=ios) job%run_params%ref_itime
            call check_value_read(ios, line_number, key)

        case ('offset')
            read(value, *, iostat=ios) job%run_params%offset
            call check_value_read(ios, line_number, key)

        case ('stride')
            read(value, *, iostat=ios) job%run_params%stride
            call check_value_read(ios, line_number, key)

        case ('size_block')
            read(value, *, iostat=ios) job%run_params%size_block
            call check_value_read(ios, line_number, key)

        case ('test_pbc')
            read(value, *, iostat=ios) job%run_params%test_pbc
            call check_value_read(ios, line_number, key)

        case ('execute_lifting')
            read(value, *, iostat=ios) job%run_params%execute_lifting
            call check_value_read(ios, line_number, key)

        ! I/O parameters
        case ('output_policy')
            job%io_params%output_policy = trim(value)

        case default
            call config_error(line_number, 'Unknown parameter: '//trim(key))

        end select

    end do

    close(unit)

    if (.not. seen_type) then
        error stop 'Workflow file has no workflow_type'
    end if

    if (.not. seen_version) then
        error stop 'Workflow file has no format_version'
    end if

end subroutine read_euclidean_lifting_workflow


! ============================================================
! FEHLERBEHANDLUNG
! ============================================================

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

end module euclidean_lifting_workflow_io_module