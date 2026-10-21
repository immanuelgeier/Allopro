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

    subroutine configure_seed_workflow(workflow)

        implicit none

        type(random_seed_workflow_type), intent(out) :: workflow

        character(len=1024) :: filename
        character(len=32) :: input
        integer :: choice, ios

        workflow = random_seed_workflow_type()

        write(*,*)
        write(*,*) '======================================'
        write(*,*) '   Random seed workflow'
        write(*,*) '======================================'
        write(*,*)
        write(*,*) '1 = Create new workflow'
        write(*,*) '2 = Load existing workflow'
        write(*,*) '3 = Skip seed workflow'
        write(*,*)
        write(*,'(A)', advance='no') 'Select option: '

        read(*,'(A)') input
        read(input,*,iostat=ios) choice

        if (ios /= 0) error stop 'Invalid seed workflow option'

        select case (choice)

        case (1)

            ! Neue Seed-Konfiguration erstellen
            call configure_seed(workflow%config)

            workflow%initialized = .true.

        case (2)

            ! Vorhandenen Workflow laden
            write(*,'(A)', advance='no') &
                'Workflow filename [seed_workflow.nml]: '

            read(*,'(A)') filename

            if (len_trim(filename) == 0) then
                filename = 'seed_workflow.nml'
            end if

            call read_seed_workflow(workflow, trim(filename))

        case (3)

            ! Workflow überspringen
            workflow%initialized = .false.

            write(*,*)
            write(*,*) 'Seed workflow skipped.'

        case default

            error stop 'Invalid seed workflow option'

        end select

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

end module random_seed_workflow_module