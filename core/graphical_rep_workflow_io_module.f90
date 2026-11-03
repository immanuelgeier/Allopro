module graphical_representation_workflow_io_module

    use workflow_module, only: graphical_representation_workflow_type, workflow_config_filename, &
        graphical_representation_workflow_config_stem
    use iso_fortran_env, only: output_unit, error_unit

    implicit none
    private

    public :: init_graphical_representation_workflow
    public :: read_graphical_representation_workflow
    public :: write_graphical_representation_workflow
    public :: print_graphical_representation_workflow
    public :: validate_graphical_representation_workflow
    public :: configure_graphical_representation_workflow

contains

! ============================================================
! KONFIGURATION (Interaktiv)
! ============================================================

subroutine configure_graphical_representation_workflow(job)
    implicit none

    type(graphical_representation_workflow_type), intent(out) :: job

    character(len=1024) :: filename
    character(len=32)   :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '  Graphical representation workflow config'
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
        job = graphical_representation_workflow_type()
        call edit_graphical_representation_workflow(job)

    case (2)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_graphical_representation_workflow(job, trim(filename))
        call edit_graphical_representation_workflow(job)

    case (3)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_graphical_representation_workflow(job, trim(filename))
        call validate_graphical_representation_workflow(job)

        write(*,*)
        write(*,*) 'Workflow loaded:'
        call print_graphical_representation_workflow(job)

    case (4)
        job = graphical_representation_workflow_type()
        write(*,*)
        write(*,*) 'Workflow configuration skipped.'

    case default
        error stop 'Invalid workflow option'

    end select

end subroutine configure_graphical_representation_workflow


! ============================================================
! EDITIEREN DER PARAMETER
! ============================================================

subroutine edit_graphical_representation_workflow(job)
    implicit none

    type(graphical_representation_workflow_type), intent(inout) :: job

    character(len=2048) :: input
    character(len=2048) :: filename
    integer :: ios
    logical :: save_config

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '  Edit Graphical Representation Workflow'
    write(*,*) '======================================'
    write(*,*) 'Press Enter to keep the current value.'
    write(*,*)

    ! --------------------------------------------------
    ! 1. Files & Scripts
    ! --------------------------------------------------
    write(*,*) 'Files & Scripts:'
    write(*,'(A,A,A)', advance='no') '  PDB filename [', trim(job%params%pdb_file), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        job%params%pdb_file = trim(input)
    end if

    write(*,'(A,A,A)', advance='no') '  TCL script filename [', trim(job%params%tcl_script), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        job%params%tcl_script = trim(input)
    end if

    ! --------------------------------------------------
    ! 2. Indices & Scaling
    ! --------------------------------------------------
    write(*,'(A,I0,A)', advance='no') '  State index [', job%params%state_index, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%params%state_index
        if (ios /= 0) error stop 'Invalid state index'
    end if

    	
    write(*,*)
    write(*,*) 'Indices & Scale Parameters:'
    write(*,'(A,I0,A)', advance='no') '  Replica index [', job%params%replica_index, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%params%replica_index
        if (ios /= 0) error stop 'Invalid replica index'
    end if

    write(*,'(A,I0,A)', advance='no') '  Time index [', job%params%time_index, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%params%time_index
        if (ios /= 0) error stop 'Invalid time index'
    end if

    write(*,'(A,F0.3,A)', advance='no') '  Scale variable [', job%params%scale, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%params%scale
        if (ios /= 0) error stop 'Invalid scale value'
    end if

    ! --------------------------------------------------
    ! 3. CV Plan & State Index Integration
    ! --------------------------------------------------
    write(*,*)
    write(*,*) 'CV Plan & State Integration:'
    write(*,'(A,A,A)', advance='no') '  CV plan filename [', trim(job%params%cv_plan_file), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        job%params%cv_plan_file = trim(input)
    end if



    ! Validate and display
    call validate_graphical_representation_workflow(job)

    write(*,*)
    write(*,*) 'Resulting workflow configuration:'
    call print_graphical_representation_workflow(job)

    ! Save configuration
    write(*,*)
    write(*,'(A)', advance='no') 'Save configuration? (y/n): '
    read(*,'(A)') input

    save_config = trim(input) == 'y' .or. trim(input) == 'Y'

    if (save_config) then
        write(*,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') filename

        filename = workflow_config_filename(trim(filename), graphical_representation_workflow_config_stem)

        call write_graphical_representation_workflow(job, trim(filename))
        write(*,'(3A)') 'Configuration saved to "', trim(filename), '"'
    end if

end subroutine edit_graphical_representation_workflow


! ============================================================
! INITIALISIERUNG
! ============================================================

subroutine init_graphical_representation_workflow(job, mode, filename)
    type(graphical_representation_workflow_type), intent(out) :: job
    character(len=*), intent(in) :: mode
    character(len=*), intent(in), optional :: filename

    job = graphical_representation_workflow_type()

    select case (trim(mode))
    case ('manual')
    case ('file')
        if (.not. present(filename)) then
            error stop 'init_graphical_representation_workflow: filename missing'
        end if
        call read_graphical_representation_workflow(job, filename)
    case default
        write(error_unit, '(A,A)') 'init_graphical_representation_workflow: unknown mode: ', trim(mode)
        error stop
    end select

    call validate_graphical_representation_workflow(job)
end subroutine init_graphical_representation_workflow


! ============================================================
! VALIDIERUNG
! ============================================================

subroutine validate_graphical_representation_workflow(job)
    implicit none
    type(graphical_representation_workflow_type), intent(in) :: job

    if (len_trim(job%params%pdb_file) == 0) error stop 'Invalid workflow: pdb_file cannot be empty'
    if (len_trim(job%params%tcl_script) == 0) error stop 'Invalid workflow: tcl_script cannot be empty'
    if (job%params%replica_index < 1) error stop 'Invalid workflow: replica_index < 1'
    if (job%params%time_index < 1) error stop 'Invalid workflow: time_index < 1'
    if (len_trim(job%params%cv_plan_file) == 0) error stop 'Invalid workflow: cv_plan_file cannot be empty'
    if (job%params%state_index < 1) error stop 'Invalid workflow: state_index < 1'
end subroutine validate_graphical_representation_workflow


! ============================================================
! DATEI SCHREIBEN & AUSGEBEN
! ============================================================

subroutine write_graphical_representation_workflow(job, filename)
    type(graphical_representation_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_graphical_representation_workflow(job)

    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'Cannot open workflow file: ', trim(filename)
        error stop
    end if

    call write_graphical_representation_workflow_unit(job, unit)
    close(unit, iostat=ios)
end subroutine write_graphical_representation_workflow


subroutine print_graphical_representation_workflow(job)
    type(graphical_representation_workflow_type), intent(in) :: job
    call write_graphical_representation_workflow_unit(job, output_unit)
end subroutine print_graphical_representation_workflow


subroutine write_graphical_representation_workflow_unit(job, unit)
    type(graphical_representation_workflow_type), intent(in) :: job
    integer, intent(in) :: unit

    write(unit, '(A)') '# Graphical representation workflow configuration'
    write(unit, '(A)') 'workflow_type = graphical_representation'
    write(unit, '(A)') 'format_version = 1'
    write(unit, '(A)') ''

    write(unit, '(A)') '# Parameters'
    call write_string(unit, 'pdb_file', job%params%pdb_file)
    call write_string(unit, 'tcl_script', job%params%tcl_script)
    write(unit, '(A,I0)') 'replica_index = ', job%params%replica_index
    write(unit, '(A,I0)') 'time_index = ', job%params%time_index
    write(unit, '(A,F0.3)') 'scale = ', job%params%scale
    call write_string(unit, 'cv_plan_file', job%params%cv_plan_file)
    write(unit, '(A,I0)') 'state_index = ', job%params%state_index
end subroutine write_graphical_representation_workflow_unit


subroutine write_string(unit, key, value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value
    write(unit, '(A," = ",A)') trim(key), trim(value)
end subroutine write_string


! ============================================================
! DATEI LESEN
! ============================================================

subroutine read_graphical_representation_workflow(job, filename)
    type(graphical_representation_workflow_type), intent(inout) :: job
    character(len=*), intent(in) :: filename

    character(len=2048) :: line
    character(len=128)  :: key
    character(len=1900) :: value
    integer :: unit, ios, line_number, eq_pos, format_version
    logical :: seen_type, seen_version

    job = graphical_representation_workflow_type()
    format_version = -1
    seen_type = .false.
    seen_version = .false.

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
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

        key = trim(adjustl(line(:eq_pos-1)))
        value = trim(adjustl(line(eq_pos+1:)))

        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'graphical_representation') then
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

        case ('pdb_file')
            job%params%pdb_file = trim(value)

        case ('tcl_script')
            job%params%tcl_script = trim(value)

        case ('replica_index')
            read(value, *, iostat=ios) job%params%replica_index
            call check_value_read(ios, line_number, key)

        case ('time_index')
            read(value, *, iostat=ios) job%params%time_index
            call check_value_read(ios, line_number, key)

        case ('scale')
            read(value, *, iostat=ios) job%params%scale
            call check_value_read(ios, line_number, key)

        case ('cv_plan_file')
            job%params%cv_plan_file = trim(value)

        case ('state_index')
            read(value, *, iostat=ios) job%params%state_index
            call check_value_read(ios, line_number, key)

        case default
            call config_error(line_number, 'Unknown parameter: '//trim(key))
        end select
    end do

    close(unit)

    if (.not. seen_type) error stop 'Workflow file has no workflow_type'
    if (.not. seen_version) error stop 'Workflow file has no format_version'
end subroutine read_graphical_representation_workflow


! ============================================================
! FEHLERBEHANDLUNG
! ============================================================

subroutine check_value_read(ios, line_number, key)
    integer, intent(in) :: ios, line_number
    character(len=*), intent(in) :: key

    if (ios /=0) then
        write(error_unit, '(A,I0,A,A)') 'Invalid value at line ', line_number, ' for parameter ', trim(key)
        error stop
    end if
end subroutine check_value_read


subroutine config_error(line_number, message)
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: message

    write(error_unit, '(A,I0,2A)') 'Workflow configuration error at line ', line_number, ': ', trim(message)
    error stop
end subroutine config_error

end module graphical_representation_workflow_io_module