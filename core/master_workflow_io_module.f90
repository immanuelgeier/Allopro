module master_workflow_io_module
    use iso_fortran_env, only: output_unit, error_unit
    use workflow_module, only: seed_workflow_config_stem, trajectory_settings_workflow_config_stem, &
        replicate_selector_workflow_config_stem, &
        library_workflow_config_stem, cv_plan_workflow_config_stem, &
        test_trajectory_workflow_config_stem, graphical_representation_workflow_config_stem, &
        master_workflow_config_stem, euclidean_lifting_workflow_config_stem
    implicit none
    private

    public :: master_workflow_type
    public :: configure_master_workflow
    public :: read_master_workflow
    public :: write_master_workflow
    public :: print_master_workflow
    public :: validate_master_workflow

    integer, parameter, public :: master_workflow_step_type_id = 7
    integer, parameter, public :: replicate_selector_workflow_step_type_id = 8
    integer, parameter, public :: euclidean_lifting_workflow_step_type_id = 9

    ! Typ für einen einzelnen Workflow-Schritt im Master-Log/Plan
    type :: workflow_step_type
        integer           :: step_type   = 1    ! Numerische ID des Workflow-Typs
        character(len=64) :: step_name   = ''   ! z.B. 'seed', 'traj_settings', etc.
        character(len=1024):: config_file = ''   ! z.B. 'trajectory_settings.dat'
        character(len=16) :: status      = 'PENDING' ! PENDING, RUNNING, DONE, FAILED
    end type workflow_step_type

    type :: master_workflow_type
        character(len=256) :: project_name = 'default_project'
        integer :: n_steps = 0
        type(workflow_step_type), allocatable :: steps(:)
    end type master_workflow_type

contains

subroutine configure_master_workflow(master)
    type(master_workflow_type), intent(out) :: master
    character(len=1024) :: input
    character(len=1024) :: folder_path
    character(len=1024), allocatable :: available_files(:), matching_files(:)
    integer :: i, n_files, n_matching, selected_type, file_choice
    
    ! Only workflow types supported by the execution dispatch are offered.
    integer, parameter :: n_known_types = 9
    character(len=128), dimension(n_known_types), parameter :: known_step_types = [ &
        'seed                    ', &
        'traj_settings           ', &
        'cv_lib                  ', &
        'cv_plan                 ', &
        'write test trajectory   ', &
        'graphical rep vmd raw cv', &
        'master_workflow         ', &
        'replicate selector      ', &
        'euclidean lifting       ' &
    ]
    character(len=40), dimension(n_known_types), parameter :: config_stems = [ character(len=40) :: &
        seed_workflow_config_stem, &
        trajectory_settings_workflow_config_stem, &
        library_workflow_config_stem, &
        cv_plan_workflow_config_stem, &
        test_trajectory_workflow_config_stem, &
        graphical_representation_workflow_config_stem, &
        master_workflow_config_stem, &
        replicate_selector_workflow_config_stem, &
        euclidean_lifting_workflow_config_stem &
    ]

    integer :: type_choice, step_count, ios, unit_tmp
    logical :: continue_adding
    
    ! Temporäres dynamisches Array für die Schritte
    type(workflow_step_type), allocatable :: temp_steps(:)
    type(workflow_step_type), allocatable :: old_steps(:)

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Interactive Master Workflow Builder'
    write(*,*) '======================================'
    
    write(*,'(A)', advance='no') 'Project name [simulation_run]: '
    read(*,'(A)') input
    if (len_trim(input) > 0) master%project_name = trim(input)

    ! Ordner abfragen, in dem die Configs liegen
    write(*,'(A)', advance='no') 'Config directory path [.]: '
    read(*,'(A)') folder_path
    if (len_trim(folder_path) == 0) then
        folder_path = '.'
    else
        folder_path = trim(adjustl(folder_path))
    end if

    call execute_command_line('ls -1 -- ' // trim(shell_quote(trim(folder_path))) // &
        '/*.cfg > .cfg_files_tmp.list 2>/dev/null', exitstat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A,A)') 'ERROR: Cannot list .cfg files in "', trim(folder_path), '".'
        error stop 'Cannot list workflow configuration files.'
    end if
    open(newunit=unit_tmp, file='.cfg_files_tmp.list', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot read workflow configuration file list.'
    n_files = 0
    do
        read(unit_tmp, '(A)', iostat=ios) input
        if (ios /= 0) exit
        n_files = n_files + 1
    end do
    rewind(unit_tmp)
    allocate(available_files(n_files))
    do i = 1, n_files
        read(unit_tmp, '(A)', iostat=ios) available_files(i)
        if (ios /= 0) error stop 'Error reading workflow configuration file list.'
    end do
    close(unit_tmp)
    call execute_command_line('rm -f .cfg_files_tmp.list', exitstat=ios)
    if (ios /= 0) error stop 'Cannot remove temporary workflow configuration file list.'
    if (n_files == 0) then
        write(error_unit, '(A,A,A)') 'ERROR: No .cfg files found in "', trim(folder_path), '".'
        error stop 'No workflow configuration files found.'
    end if

    write(*,*)
    write(*,*) '--- Add Workflow Steps to Master Plan ---'

    step_count = 0
    continue_adding = .true.

    do while (continue_adding)
        step_count = step_count + 1
        
        ! Array dynamisch vergrößern (Reallocation-Pattern für Fortran)
        if (step_count == 1) then
            allocate(temp_steps(1))
        else
            allocate(old_steps(step_count - 1))
            old_steps = temp_steps
            deallocate(temp_steps)
            allocate(temp_steps(step_count))
            temp_steps(1:step_count-1) = old_steps
            deallocate(old_steps)
        end if

        write(*,*)
        write(*, '(A,I0,A)') '--- Configuring Step #', step_count, ' ---'
        
        do
            write(*,*) 'Select Workflow Type:'
            do i = 1, n_known_types
                write(*, '(I4, A, A)') i, ' = ', trim(known_step_types(i))
            end do

            write(*,'(A,I0,A)', advance='no') 'Choice [1]: '
            read(*,'(A)') input
            type_choice = 1
            if (len_trim(input) > 0) then
                read(input, *, iostat=ios) type_choice
                if (ios /= 0 .or. type_choice < 1 .or. type_choice > n_known_types) type_choice = 1
            end if
            selected_type = type_choice

            allocate(matching_files(n_files))
            n_matching = 0
            do i = 1, n_files
                if (is_workflow_config_file(trim(available_files(i)), trim(config_stems(selected_type)))) then
                    n_matching = n_matching + 1
                    matching_files(n_matching) = available_files(i)
                end if
            end do
            if (n_matching > 0) exit

            write(error_unit, '(A,A,A)') 'No configuration files found for workflow type "', &
                trim(known_step_types(selected_type)), '". Please choose another type.'
            deallocate(matching_files)
        end do

        write(*,*) 'Select Configuration File (.cfg):'
        do i = 1, n_matching
            write(*, '(I4, A, A)') i, ' = ', trim(matching_files(i))
        end do
        write(*,'(A,I0,A)', advance='no') 'Choice [1]: '
        read(*,'(A)') input
        file_choice = 1
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) file_choice
            if (ios /= 0 .or. file_choice < 1 .or. file_choice > n_matching) file_choice = 1
        end if
        
        ! Speichere sowohl den numerischen Typ als auch den Namen
        temp_steps(step_count)%step_type = selected_type
        temp_steps(step_count)%step_name = trim(known_step_types(selected_type))
        temp_steps(step_count)%config_file = trim(matching_files(file_choice))
        deallocate(matching_files)
        
        ! Standardstatus
        temp_steps(step_count)%status = 'PENDING'

        ! 3. Nachfragen, ob ein weiterer Schritt hinzugefügt werden soll
        write(*,*)
        write(*,'(A)', advance='no') 'Add another workflow step? (y/n) [y]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'y'
        input = trim(adjustl(input))
        continue_adding = (input == 'y' .or. input == 'Y' .or. input == 'yes')
    end do

    ! Schritte in das finale Master-Objekt übertragen
    master%n_steps = step_count
    if (allocated(master%steps)) deallocate(master%steps)
    allocate(master%steps(step_count))
    master%steps = temp_steps

    if (allocated(temp_steps)) deallocate(temp_steps)
    if (allocated(available_files)) deallocate(available_files)

    write(*,*)
    write(*,'(A,A,A)', advance='no') 'Master protocol filename [master_workflow.cfg]: '
    read(*,'(A)') input
    
    if (len_trim(input) == 0) then
        input = 'master_workflow.cfg'
    else
        input = trim(adjustl(input))
    end if

    call validate_master_workflow(master)
    call write_master_workflow(master, trim(input))
    write(*, '(3A)') '-> Master protocol successfully written to "', trim(input), '"'
    
end subroutine configure_master_workflow

function ends_with(value, suffix) result(matches)
    character(len=*), intent(in) :: value, suffix
    logical :: matches
    integer :: value_length, suffix_length

    value_length = len_trim(value)
    suffix_length = len_trim(suffix)
    matches = value_length >= suffix_length
    if (matches) matches = value(value_length-suffix_length+1:value_length) == trim(suffix)
end function ends_with

function is_workflow_config_file(filename, standard_name) result(matches)
    character(len=*), intent(in) :: filename, standard_name
    logical :: matches
    character(len=64) :: standard_filename

    standard_filename = trim(standard_name) // '.cfg'
    matches = trim(filename) == trim(standard_filename) .or. &
        ends_with(trim(filename), '/' // trim(standard_filename)) .or. &
        ends_with(trim(filename), '_' // trim(standard_filename))
end function is_workflow_config_file

function shell_quote(value) result(quoted)
    character(len=*), intent(in) :: value
    character(len=2*len(value)+2) :: quoted
    integer :: i, position

    quoted = "'"
    position = 2
    do i = 1, len_trim(value)
        if (value(i:i) == "'") then
            quoted(position:position+3) = "'\''"
            position = position + 4
        else
            quoted(position:position) = value(i:i)
            position = position + 1
        end if
    end do
    quoted(position:position) = "'"
end function shell_quote

subroutine read_master_workflow(master, filename)
    type(master_workflow_type), intent(out) :: master
    character(len=*), intent(in) :: filename
    integer :: unit, ios, i
    character(len=2048) :: line

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot open master protocol file'

    ! 1. Projekt-Name einlesen (Kommentare überspringen)
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) error stop 'Error reading project name from master protocol'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        master%project_name = trim(line)
        exit
    end do

    ! 2. Anzahl der Schritte einlesen
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) error stop 'Error reading step count from master protocol'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        read(line, *, iostat=ios) master%n_steps
        if (ios /= 0) error stop 'Invalid step count format in master protocol'
        exit
    end do

    if (master%n_steps < 1) then
        write(*,*) 'WARNING: Master workflow file contains 0 steps!'
        close(unit)
        return
    end if

    if (allocated(master%steps)) deallocate(master%steps)
    allocate(master%steps(master%n_steps))

    ! 3. Jeder Schritt besteht aus 4 aufeinanderfolgenden Werten
    do i = 1, master%n_steps
        ! step_type
        read(unit, *, iostat=ios) master%steps(i)%step_type
        if (ios /= 0) error stop 'Error reading step_type'

        ! step_name
        read(unit, '(A)', iostat=ios) master%steps(i)%step_name
        if (ios /= 0) error stop 'Error reading step_name'

        ! config_file
        read(unit, '(A)', iostat=ios) master%steps(i)%config_file
        if (ios /= 0) error stop 'Error reading config_file'

        ! status
        read(unit, '(A)', iostat=ios) master%steps(i)%status
        if (ios /= 0) error stop 'Error reading status'

        master%steps(i)%step_name   = adjustl(trim(master%steps(i)%step_name))
        master%steps(i)%config_file = adjustl(trim(master%steps(i)%config_file))
        master%steps(i)%status      = adjustl(trim(master%steps(i)%status))
    end do

    close(unit)
end subroutine read_master_workflow


subroutine write_master_workflow(master, filename)
    type(master_workflow_type), intent(in) :: master
    character(len=*), intent(in) :: filename
    integer :: unit, ios, i

    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Cannot write master protocol file'

    write(unit, '(A)') '# Master Workflow Protocol & Log'
    write(unit, '(A)') trim(master%project_name)
    write(unit, '(I0)') master%n_steps

    do i = 1, master%n_steps
        write(unit, '(I0)') master%steps(i)%step_type
        write(unit, '(A)') trim(master%steps(i)%step_name)
        write(unit, '(A)') trim(master%steps(i)%config_file)
        write(unit, '(A)') trim(master%steps(i)%status)
    end do
    close(unit)
end subroutine write_master_workflow


subroutine print_master_workflow(master, filename)
    type(master_workflow_type), intent(in) :: master
    character(len=*), intent(in) :: filename
    integer :: unit, i

    unit = output_unit
    write(unit, '(A)') '# Master Workflow Protocol & Log'
    write(unit, '(A)') trim(master%project_name)
    write(unit, '(I0)') master%n_steps

    do i = 1, master%n_steps
        write(unit, '(I0)') master%steps(i)%step_type
        write(unit, '(A)') trim(master%steps(i)%step_name)
        write(unit, '(A)') trim(master%steps(i)%config_file)
        write(unit, '(A)') trim(master%steps(i)%status)
    end do
end subroutine print_master_workflow


subroutine validate_master_workflow(master)
    type(master_workflow_type), intent(in) :: master
    integer :: i
    if (master%n_steps < 1) error stop 'Master error: No steps defined.'
    do i = 1, master%n_steps
        if (master%steps(i)%step_type < 1) error stop 'Master error: Invalid step type index.'
        if (len_trim(master%steps(i)%step_name) == 0) error stop 'Master error: Empty step name.'
        if (len_trim(master%steps(i)%config_file) == 0) error stop 'Master error: Empty config file.'
    end do
end subroutine validate_master_workflow

end module master_workflow_io_module