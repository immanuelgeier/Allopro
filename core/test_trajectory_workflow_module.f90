module test_trajectory_workflow_module
    use workflow_module, only: workflow_config_filename, test_trajectory_workflow_config_stem
    use iso_fortran_env, only: output_unit, error_unit
    implicit none

    private
    public :: test_trajectory_workflow_type, &
              edit_test_trajectory_workflow, &
              read_test_trajectory_workflow, &
              write_test_trajectory_workflow, &
              print_test_trajectory_workflow, &
              validate_test_trajectory_workflow, &
              configure_test_trajectory_workflow
              

    ! ==========================================================
    ! DATENTYP FÜR DEN TEST-TRAJEKTORIE-WORKFLOW
    ! ==========================================================
    type :: test_trajectory_workflow_type
        logical :: initialized = .false.
        
        ! Zuordnung und Zeitschritte
        integer :: istate          = 1
        integer :: irep            = 1
        integer :: start_time_step = 1
        integer :: end_time_step   = 1000
        integer :: stride          = 1

        ! Dateinamen
        character(len=1024) :: test_trajectory_filename = 'test_traj.dat'
        character(len=1024) :: cv_plan_filename         = 'cv_plan.dat'
    end type test_trajectory_workflow_type

contains

subroutine configure_test_trajectory_workflow(job)
    implicit none

    type(test_trajectory_workflow_type), intent(out) :: job

    character(len=1024) :: filename
    character(len=32)   :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    Test trajectory workflow configuration'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Create new test trajectory workflow'
    write(*,*) '2 = Load and edit existing test trajectory workflow'
    write(*,*) '3 = Load existing test trajectory workflow unchanged'
    write(*,*) '4 = Skip test trajectory workflow configuration'
    write(*,*)
    write(*,'(A)', advance='no') 'Select option: '

    read(*,'(A)') input
    read(input,*,iostat=ios) choice

    if (ios /= 0) error stop 'Invalid test trajectory workflow option'

    select case (choice)

    case (1)
        job = test_trajectory_workflow_type()
        call edit_test_trajectory_workflow(job)
        job%initialized = .true.

    case (2)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) error stop 'No configuration filename specified'

        call read_test_trajectory_workflow(job, trim(filename))
        call edit_test_trajectory_workflow(job)
        job%initialized = .true.

    case (3)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) error stop 'No configuration filename specified'

        call read_test_trajectory_workflow(job, trim(filename))
        call validate_test_trajectory_workflow(job)

        write(*,*)
        write(*,*) 'Test trajectory workflow loaded:'
        call print_test_trajectory_workflow(job)
        job%initialized = .true.

    case (4)
        job = test_trajectory_workflow_type()
        write(*,*) 
        write(*,*) 'Test trajectory workflow configuration skipped.'

    case default
        error stop 'Invalid test trajectory workflow option'

    end select

end subroutine configure_test_trajectory_workflow

    ! ==========================================================
    ! INTERAKTIVES EDITIEREN (CONFIG-MODUS)
    ! ==========================================================
    subroutine edit_test_trajectory_workflow(job)
        type(test_trajectory_workflow_type), intent(inout) :: job

        character(len=1024) :: input, filename
        integer :: ios
        logical :: save_config

        write(*,*)
        write(*,*) '======================================'
        write(*,*) '  Edit Test Trajectory Workflow       '
        write(*,*) '======================================'
        write(*,*) 'Press Enter to keep the current value.'
        write(*,*)

        ! 1. State (istate) und Replikat (irep) abfragen
        write(*,'(A,I0,A)', advance='no') 'State index (istate) [', job%istate, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) job%istate
            if (ios /= 0 .or. job%istate < 1) error stop 'Invalid istate value (must be >= 1)'
        end if

        write(*,'(A,I0,A)', advance='no') 'Replica index (irep) [', job%irep, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) job%irep
            if (ios /= 0 .or. job%irep < 1) error stop 'Invalid irep value (must be >= 1)'
        end if

        write(*,*)

        ! 2. Dateinamen abfragen
        write(*,'(A,A,A)', advance='no') 'Test trajectory output filename [', trim(job%test_trajectory_filename), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%test_trajectory_filename = trim(input)

        write(*,'(A,A,A)', advance='no') 'Input CV plan filename [', trim(job%cv_plan_filename), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%cv_plan_filename = trim(input)

        ! 3. Zeitschritte und Stride abfragen
        write(*,*)
        write(*,'(A,I0,A)', advance='no') 'Start time step [', job%start_time_step, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) job%start_time_step
            if (ios /= 0 .or. job%start_time_step < 1) error stop 'Invalid start_time_step value (must be >= 1)'
        end if

        write(*,'(A,I0,A)', advance='no') 'End time step [', job%end_time_step, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) job%end_time_step
            if (ios /= 0 .or. job%end_time_step < 1) error stop 'Invalid end_time_step value (must be >= 1)'
        end if

        write(*,'(A,I0,A)', advance='no') 'Stride (step interval) [', job%stride, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) job%stride
            if (ios /= 0 .or. job%stride < 1) error stop 'Invalid stride value (must be >= 1)'
        end if

        ! Validieren und Speichern
        call validate_test_trajectory_workflow(job)

        write(*,*)
        write(*,*) 'Resulting workflow configuration:'
        call print_test_trajectory_workflow(job)

        write(*,*)
        write(*,'(A)', advance='no') 'Save configuration? (y/n) [y]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'y'

        save_config = (trim(input) == 'y' .or. trim(input) == 'Y' .or. trim(input) == 'yes')

        if (save_config) then
            write(*,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
            read(*,'(A)') filename

            filename = workflow_config_filename(trim(filename), test_trajectory_workflow_config_stem)

            call write_test_trajectory_workflow(job, trim(filename))
            write(*,'(3A)') 'Configuration saved to "', trim(filename), '"'
        end if

    end subroutine edit_test_trajectory_workflow


    ! ==========================================================
    ! VALIDIERUNG
    ! ==========================================================
    subroutine validate_test_trajectory_workflow(job)
        type(test_trajectory_workflow_type), intent(in) :: job

        if (job%istate < 1) then
            write(error_unit, *) 'ERROR: istate must be at least 1.'
            error stop 'Invalid test trajectory workflow parameters'
        end if

        if (job%irep < 1) then
            write(error_unit, *) 'ERROR: irep must be at least 1.'
            error stop 'Invalid test trajectory workflow parameters'
        end if

        if (job%start_time_step < 1) then
            write(error_unit, *) 'ERROR: start_time_step must be at least 1.'
            error stop 'Invalid test trajectory workflow parameters'
        end if

        if (job%end_time_step < job%start_time_step) then
            write(error_unit, *) 'ERROR: end_time_step cannot be less than start_time_step.'
            error stop 'Invalid test trajectory workflow parameters'
        end if

        if (job%stride < 1) then
            write(error_unit, *) 'ERROR: stride must be at least 1.'
            error stop 'Invalid test trajectory workflow parameters'
        end if
    end subroutine validate_test_trajectory_workflow


    ! ==========================================================
    ! SCHREIBEN IN DATEI
    ! ==========================================================
    subroutine write_test_trajectory_workflow(job, filename)
        type(test_trajectory_workflow_type), intent(in) :: job
        character(len=*), intent(in) :: filename
        integer :: unit, ios

        call validate_test_trajectory_workflow(job)

        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Cannot open test trajectory workflow file for writing: ', trim(filename)
            error stop
        end if

        call write_test_trajectory_workflow_unit(job, unit)
        close(unit, iostat=ios)
    end subroutine write_test_trajectory_workflow


    ! ==========================================================
    ! AUSGABE AUF KONSOLE
    ! ==========================================================
    subroutine print_test_trajectory_workflow(job)
        type(test_trajectory_workflow_type), intent(in) :: job
        call write_test_trajectory_workflow_unit(job, output_unit)
    end subroutine print_test_trajectory_workflow


    ! Interne Einheit zum Schreiben
    subroutine write_test_trajectory_workflow_unit(job, unit)
        type(test_trajectory_workflow_type), intent(in) :: job
        integer, intent(in) :: unit

        write(unit, '(A)') '# Test trajectory workflow configuration'
        write(unit, '(A)') 'workflow_type = test_trajectory'
        write(unit, '(A)') 'format_version = 1'
        write(unit, '(A)') ''

        write(unit, '(A,I0)')   'istate                   = ', job%istate
        write(unit, '(A,I0)')   'irep                     = ', job%irep
        write(unit, '(A,A)')    'test_trajectory_filename = ', trim(job%test_trajectory_filename)
        write(unit, '(A,A)')    'cv_plan_filename         = ', trim(job%cv_plan_filename)
        write(unit, '(A,I0)')   'start_time_step          = ', job%start_time_step
        write(unit, '(A,I0)')   'end_time_step            = ', job%end_time_step
        write(unit, '(A,I0)')   'stride                   = ', job%stride
    end subroutine write_test_trajectory_workflow_unit


    ! ==========================================================
    ! EINLESEN AUS DATEI
    ! ==========================================================
    subroutine read_test_trajectory_workflow(job, filename)
        type(test_trajectory_workflow_type), intent(out) :: job
        character(len=*), intent(in) :: filename

        character(len=2048) :: line
        character(len=128)  :: key
        character(len=1900) :: value

        integer :: unit, ios, line_number, eq_pos, format_version
        logical :: seen_type, seen_version

        job = test_trajectory_workflow_type()
        format_version = -1
        seen_type = .false.
        seen_version = .false.

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Cannot open test trajectory workflow file for reading: ', trim(filename)
            error stop
        end if

        line_number = 0

        do
            read(unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) then
                close(unit)
                error stop 'Error reading test trajectory workflow file'
            end if

            line_number = line_number + 1
            line = adjustl(line)

            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle

            eq_pos = index(line, '=')
            if (eq_pos == 0) then
                write(error_unit, '(A,I0,A)') 'Error at line ', line_number, ': Missing "="'
                error stop
            end if

            key = trim(adjustl(line(:eq_pos-1)))
            value = trim(adjustl(line(eq_pos+1:)))

            select case (trim(key))
            case ('workflow_type')
                if (trim(value) /= 'test_trajectory') then
                    write(error_unit, '(A,A)') 'Wrong workflow_type in file: ', trim(value)
                    error stop
                end if
                seen_type = .true.

            case ('format_version')
                read(value, *, iostat=ios) format_version
                if (ios /= 0 .or. format_version /= 1) then
                    write(error_unit, *) 'Unsupported or invalid format_version'
                    error stop
                end if
                seen_version = .true.

            case ('istate')
                read(value, *, iostat=ios) job%istate
                if (ios /= 0) error stop 'Error reading istate'

            case ('irep')
                read(value, *, iostat=ios) job%irep
                if (ios /= 0) error stop 'Error reading irep'

            case ('test_trajectory_filename')
                job%test_trajectory_filename = trim(value)

            case ('cv_plan_filename')
                job%cv_plan_filename = trim(value)

            case ('start_time_step')
                read(value, *, iostat=ios) job%start_time_step
                if (ios /= 0) error stop 'Error reading start_time_step'

            case ('end_time_step')
                read(value, *, iostat=ios) job%end_time_step
                if (ios /= 0) error stop 'Error reading end_time_step'

            case ('stride')
                read(value, *, iostat=ios) job%stride
                if (ios /= 0) error stop 'Error reading stride'

            case default
                write(error_unit, '(A,I0,A,A)') 'Error at line ', line_number, ': Unknown parameter -> ', trim(key)
                error stop
            end select
        end do

        close(unit)

        if (.not. seen_type) error stop 'Test trajectory workflow file has no workflow_type'
        if (.not. seen_version) error stop 'Test trajectory workflow file has no format_version'

        job%initialized = .true.
        call validate_test_trajectory_workflow(job)
    end subroutine read_test_trajectory_workflow

end module test_trajectory_workflow_module