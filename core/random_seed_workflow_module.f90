module random_seed_workflow_module

    use random_seed_io_module, only: &
        random_seed_config_type, &
        configure_seed, &
        initialize_random_seed

    implicit none
    private

    public :: random_seed_workflow_type
    public :: configure_seed_workflow
    public :: read_seed_workflow
    public :: write_seed_workflow
    public :: run_seed_workflow

    type :: random_seed_workflow_type
        logical :: initialized = .false.
        type(random_seed_config_type) :: config
    end type random_seed_workflow_type

contains

    ! ==================================================
    ! Workflow konfigurieren
    ! ==================================================
subroutine configure_seed_workflow(job)

    implicit none

    type(random_seed_workflow_type), intent(out) :: job

    character(len=1024) :: filename
    character(len=32)   :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Random seed workflow configuration'
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
        job = random_seed_workflow_type()

        call edit_seed_workflow(job)

        job%initialized = .true.

    case (2)

        ! Vorhandene Konfiguration laden und bearbeiten
        write(*,'(A)', advance='no') &
            'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_seed_workflow(trim(filename), job)
        call edit_seed_workflow(job)

        job%initialized = .true.

    case (3)

        ! Vorhandene Konfiguration unverändert übernehmen
        write(*,'(A)', advance='no') &
            'Configuration filename: '
        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No configuration filename specified'
        end if

        call read_seed_workflow(trim(filename), job)
        call validate_seed_workflow(job)

        write(*,*)
        write(*,*) 'Workflow loaded:'
        call print_seed_workflow(job)

        job%initialized = .true.

    case (4)

        ! Workflow-Konfiguration überspringen
        job = random_seed_workflow_type()

        write(*,*)
        write(*,*) 'Workflow configuration skipped.'

    case default

        error stop 'Invalid workflow option'

    end select

    ! Konfiguration nach erfolgreicher Auswahl speichern.
    if (job%initialized) then

        write(*,'(A)', advance='no') &
            'Save workflow configuration? (y/n): '

        read(*,'(A)') input

        select case (trim(adjustl(input)))

        case ('y', 'Y', 'yes', 'YES')

            write(*,'(A)', advance='no') &
                'Configuration filename [seed_workflow.nml]: '
            read(*,'(A)') filename

            if (len_trim(filename) == 0) then
                filename = 'seed_workflow.nml'
            end if

            call write_seed_workflow(trim(filename), job)

        case ('n', 'N', 'no', 'NO')

            continue

        case default

            error stop 'Invalid save-workflow option'

        end select

    end if

end subroutine configure_seed_workflow

    ! ==================================================
    ! Workflow speichern
    ! ==================================================

    subroutine write_seed_workflow(filename, workflow)

        implicit none

        character(len=*), intent(in) :: filename
        type(random_seed_workflow_type), intent(in) :: workflow

        integer :: u, ios
        character(len=512) :: iomsg

        logical :: initialized
        logical :: use_saved_seed
        logical :: save_seed
        character(len=1024) :: seed_file

        namelist /seed_workflow/ &
            initialized, use_saved_seed, save_seed, seed_file

        initialized    = workflow%initialized
        use_saved_seed = workflow%config%use_saved_seed
        save_seed      = workflow%config%save_seed
        seed_file      = workflow%config%seed_file

        open( &
            newunit=u, &
            file=trim(filename), &
            status='replace', &
            action='write', &
            iostat=ios, &
            iomsg=iomsg)

        if (ios /= 0) then
            write(*,*) 'Could not open seed workflow file:'
            write(*,*) trim(filename)
            write(*,*) trim(iomsg)
            error stop 'write_seed_workflow: open failed'
        end if

        write(u, nml=seed_workflow, iostat=ios, iomsg=iomsg)

        if (ios /= 0) then
            close(u)
            write(*,*) trim(iomsg)
            error stop 'write_seed_workflow: write failed'
        end if

        close(u, iostat=ios, iomsg=iomsg)

        if (ios /= 0) then
            write(*,*) trim(iomsg)
            error stop 'write_seed_workflow: close failed'
        end if

    end subroutine write_seed_workflow


    ! ==================================================
    ! Workflow laden
    ! ==================================================

    subroutine read_seed_workflow(filename, workflow)

        implicit none

        character(len=*), intent(in) :: filename
        type(random_seed_workflow_type), intent(out) :: workflow

        integer :: u, ios
        character(len=512) :: iomsg

        logical :: initialized
        logical :: use_saved_seed
        logical :: save_seed
        character(len=1024) :: seed_file

        namelist /seed_workflow/ &
            initialized, use_saved_seed, save_seed, seed_file

        workflow = random_seed_workflow_type()

        initialized    = .false.
        use_saved_seed = .false.
        save_seed      = .true.
        seed_file      = 'random_seed.dat'

        open( &
            newunit=u, &
            file=trim(filename), &
            status='old', &
            action='read', &
            iostat=ios, &
            iomsg=iomsg)

        if (ios /= 0) then
            write(*,*) 'Could not open seed workflow file:'
            write(*,*) trim(filename)
            write(*,*) trim(iomsg)
            error stop 'read_seed_workflow: open failed'
        end if

        read(u, nml=seed_workflow, iostat=ios, iomsg=iomsg)

        close(u)

        if (ios /= 0) then
            write(*,*) 'Could not read seed workflow:'
            write(*,*) trim(filename)
            write(*,*) trim(iomsg)
            error stop 'read_seed_workflow: read failed'
        end if

        workflow%initialized = initialized

        workflow%config%use_saved_seed = use_saved_seed
        workflow%config%save_seed      = save_seed
        workflow%config%seed_file      = trim(seed_file)

    end subroutine read_seed_workflow


    ! ==================================================
    ! Workflow ausführen
    ! ==================================================

    subroutine run_seed_workflow(workflow)

        implicit none

        type(random_seed_workflow_type), intent(inout) :: workflow

        if (.not. workflow%initialized) then
            write(*,*) 'Seed workflow not initialized - skipping.'
            return
        end if

        call initialize_random_seed(workflow%config)

    end subroutine run_seed_workflow
subroutine edit_seed_workflow(job)

    implicit none

    type(random_seed_workflow_type), intent(inout) :: job

    character(len=1024) :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Edit random seed workflow'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Generate a new seed'
    write(*,*) '2 = Load a saved seed'
    write(*,*)
    write(*,'(A)', advance='no') 'Select option: '

    read(*,'(A)') input
    read(input,*,iostat=ios) choice

    if (ios /= 0) error stop 'Invalid seed option'

    select case (choice)

    case (1)
        job%config%use_saved_seed = .false.

    case (2)
        job%config%use_saved_seed = .true.

    case default
        error stop 'Invalid seed option'

    end select

    write(*,*)
    write(*,'(A)', advance='no') 'Save seed to file? (y/n): '

    read(*,'(A)') input

    select case (trim(adjustl(input)))

    case ('y', 'Y', 'yes', 'YES')
        job%config%save_seed = .true.

    case ('n', 'N', 'no', 'NO')
        job%config%save_seed = .false.

    case default
        error stop 'Invalid save-seed option'

    end select

    if (job%config%use_saved_seed .or. job%config%save_seed) then

        write(*,'(A)', advance='no') &
            'Seed filename [random_seed.dat]: '

        read(*,'(A)') input

        if (len_trim(input) > 0) then
            job%config%seed_file = trim(input)
        else if (len_trim(job%config%seed_file) == 0) then
            job%config%seed_file = 'random_seed.dat'
        end if

    end if

    call validate_seed_workflow(job)

end subroutine edit_seed_workflow
subroutine print_seed_workflow(job)

    implicit none

    type(random_seed_workflow_type), intent(in) :: job

    write(*,*)
    write(*,*) '--------------------------------------'
    write(*,*) 'Random seed workflow'
    write(*,*) '--------------------------------------'

    write(*,'(A,L1)') &
        'Initialized:       ', job%initialized

    if (.not. job%initialized) then
        write(*,*) 'Workflow is disabled.'
        return
    end if

    if (job%config%use_saved_seed) then
        write(*,*) 'Seed source:       saved seed'
    else
        write(*,*) 'Seed source:       newly generated seed'
    end if

    write(*,'(A,L1)') &
        'Save seed:         ', job%config%save_seed

    if (job%config%use_saved_seed .or. job%config%save_seed) then
        write(*,'(A,A)') &
            'Seed filename:     ', trim(job%config%seed_file)
    end if

    write(*,*) '--------------------------------------'

end subroutine print_seed_workflow

subroutine validate_seed_workflow(job)

    implicit none

    type(random_seed_workflow_type), intent(in) :: job

    if (.not. job%initialized) return

    if (job%config%use_saved_seed .or. job%config%save_seed) then

        if (len_trim(job%config%seed_file) == 0) then
            error stop 'Seed filename must not be empty'
        end if

    end if

end subroutine validate_seed_workflow

end module random_seed_workflow_module