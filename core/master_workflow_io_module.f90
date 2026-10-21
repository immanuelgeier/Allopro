module master_workflow_io_module
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    private

    public :: master_workflow_type
    public :: configure_master_workflow
    public :: read_master_workflow
    public :: write_master_workflow
    public :: print_master_workflow
    public :: validate_master_workflow

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
    character(len=1024), allocatable :: available_files(:)
    integer :: i
    
    ! Vordefinierte, nummerierte Liste der Workflow-Typen
    integer, parameter :: n_known_types = 5
    character(len=32), dimension(n_known_types), parameter :: known_step_types = [ &
        'seed         ', &
        'traj_settings', &
        'cv_lib       ', &
        'cv_plan      ', &
        'lifting      '  &
    ]

    integer :: n_files, file_choice, type_choice
    integer :: step_count, ios, unit_tmp
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

    ! .cfg-Dateien aus dem Ordner automatisch einlesen
    call execute_command_line('ls -1 ' // trim(folder_path) // '/*.cfg > .cfg_files_tmp.list 2>/dev/null', exitstat=ios)
    
    n_files = 0
    if (ios == 0) then
        open(newunit=unit_tmp, file='.cfg_files_tmp.list', status='old', action='read', iostat=ios)
        do
            read(unit_tmp, '(A)', iostat=ios) input
            if (ios /= 0) exit
            n_files = n_files + 1
        end do
        rewind(unit_tmp)

        allocate(available_files(n_files))
        do i = 1, n_files
            read(unit_tmp, '(A)') available_files(i)
        end do
        close(unit_tmp)
        call execute_command_line('rm -f .cfg_files_tmp.list')
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
        
        ! 1. Workflow-Typ aus nummerierter Liste wählen
        write(*,*) 'Select Workflow Type:'
        do type_choice = 1, n_known_types
            write(*, '(I4, A, A)') type_choice, ' = ', trim(known_step_types(type_choice))
        end do
        
        write(*,'(A,I0,A)', advance='no') 'Choice [1]: '
        read(*,'(A)') input
        type_choice = 1
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) type_choice
            if (ios /= 0 .or. type_choice < 1 .or. type_choice > n_known_types) type_choice = 1
        end if
        
        ! Speichere sowohl den numerischen Typ als auch den Namen
        temp_steps(step_count)%step_type = type_choice
        temp_steps(step_count)%step_name = trim(known_step_types(type_choice))

        ! 2. Config-Datei aus nummerierter Ordner-Liste wählen
        if (n_files > 0) then
            write(*,*) 'Select Configuration File (.cfg):'
            do file_choice = 1, n_files
                write(*, '(I4, A, A)') file_choice, ' = ', trim(available_files(file_choice))
            end do
            
            write(*,'(A,I0,A)', advance='no') 'Choice [1]: '
            read(*,'(A)') input
            file_choice = 1
            if (len_trim(input) > 0) then
                read(input, *, iostat=ios) file_choice
                if (ios /= 0 .or. file_choice < 1 .or. file_choice > n_files) file_choice = 1
            end if
            temp_steps(step_count)%config_file = trim(available_files(file_choice))
        else
            write(*,*) 'No .cfg files found in folder. Please enter filename manually: '
            read(*,'(A)') input
            temp_steps(step_count)%config_file = trim(adjustl(input))
        end if
        
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