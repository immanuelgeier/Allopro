module random_seed_io_module

    use iso_fortran_env, only: int64
    implicit none

    private

    public :: random_seed_config_type,configure_seed
    public :: initialize_random_seed
    public :: initialize_random_seed_from_state
    public :: save_random_seed_state

    type :: random_seed_config_type
        logical :: use_saved_seed = .false.
        logical :: save_seed      = .true.
        logical :: resume_state   = .false.

        character(len=1024) :: seed_file = 'random_seed.dat'
        character(len=1024) :: state_file = 'random_state.dat'

        integer, allocatable :: seed(:)
    end type random_seed_config_type

contains

subroutine configure_seed(seed_config)

    implicit none

    type(random_seed_config_type), intent(out) :: seed_config

    character(len=1024) :: input
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Random seed configuration'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Generate a new seed'
    write(*,*) '2 = Load a saved seed'
    write(*,*)
    write(*,'(A)', advance='no') 'Select option: '

    read(*,'(A)') input
    read(input,*,iostat=ios) choice

    if (ios /= 0) error stop 'Invalid seed option'

    ! Standardwerte
    seed_config = random_seed_config_type()
    seed_config%seed_file = 'random_seed.dat'

    select case (choice)

    case (1)
        seed_config%use_saved_seed = .false.

    case (2)
        seed_config%use_saved_seed = .true.

    case default
        error stop 'Invalid seed option'

    end select

    ! Unabhängig davon fragen, ob der Seed gespeichert werden soll
    write(*,*)
    write(*,'(A)', advance='no') 'Save seed to file? (y/n): '

    read(*,'(A)') input

    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES')
        seed_config%save_seed = .true.
    case ('n', 'N', 'no', 'NO')
        seed_config%save_seed = .false.
    case default
        error stop 'Invalid save-seed option'
    end select

    ! Dateiname abfragen, wenn geladen oder gespeichert wird
    if (seed_config%use_saved_seed .or. seed_config%save_seed) then

        write(*,'(A)', advance='no') &
            'Seed filename [random_seed.dat]: '
        read(*,'(A)') input

        if (len_trim(input) > 0) then
            seed_config%seed_file = trim(input)
        end if

    end if

end subroutine configure_seed


    subroutine initialize_random_seed(config)

        type(random_seed_config_type), intent(inout) :: config

        integer :: n
        integer :: i
        integer :: unit
        integer :: ios
        integer :: n_file

        integer(int64) :: clock64

        ! Benötigte Laenge des Seed-Vektors bestimmen.
        call random_seed(size=n)

        ! Vorhandenen Seed einlesen oder neuen erzeugen.
        if (config%use_saved_seed) then

            open(newunit=unit, file=trim(config%seed_file), &
                 status='old', action='read', iostat=ios)

            if (ios /= 0) then
                error stop 'Could not open random seed file'
            end if

            read(unit, *, iostat=ios) n_file

            if (ios /= 0) then
                close(unit)
                error stop 'Could not read seed vector size'
            end if

            if (n_file /= n) then
                close(unit)
                error stop 'Seed vector size does not match compiler'
            end if

            if (allocated(config%seed)) deallocate(config%seed)
            allocate(config%seed(n))

            read(unit, *, iostat=ios) config%seed
            close(unit)

            if (ios /= 0) then
                error stop 'Could not read seed vector'
            end if

        else

            if (allocated(config%seed)) deallocate(config%seed)
            allocate(config%seed(n))

            call system_clock(count=clock64)

            ! Neuen Seed aus dem Systemtakt ableiten.
            do i = 1, n
                config%seed(i) = int( &
                    modulo(clock64 + int(i-1, int64), &
                           int(huge(0), int64)) )
            end do

        end if

        ! Zufallsgenerator mit dem verwendeten Seed initialisieren.
        call random_seed(put=config%seed)

        ! Seed bei Bedarf fuer eine spaetere Wiederholung speichern.
        if (config%save_seed) then

            open(newunit=unit, file=trim(config%seed_file), &
                 status='replace', action='write', iostat=ios)

            if (ios /= 0) then
                error stop 'Could not open seed file for writing'
            end if

            write(unit, *) n
            write(unit, *) config%seed

            close(unit)

        end if

    end subroutine initialize_random_seed

    subroutine initialize_random_seed_from_state(config, state_filename)
        type(random_seed_config_type), intent(inout) :: config
        character(len=*), intent(in) :: state_filename

        integer, allocatable :: state(:)
        integer :: n, n_file, unit, ios, close_ios
        logical :: file_exists, seed_file_exists

        inquire(file=trim(state_filename), exist=file_exists)
        if (.not. file_exists) then
            write(*,'(3A)') 'No saved generator state at "', trim(state_filename), &
                '"; initializing from the configured seed.'
            if (config%use_saved_seed) then
                inquire(file=trim(config%seed_file), exist=seed_file_exists)
                if (.not. seed_file_exists .and. config%save_seed) then
                    write(*,'(3A)') 'No saved seed at "', trim(config%seed_file), &
                        '"; generating and saving the initial seed.'
                    config%use_saved_seed = .false.
                end if
            end if
            call initialize_random_seed(config)
            return
        end if

        call random_seed(size=n)
        open(newunit=unit, file=trim(state_filename), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Could not open random generator state file'

        read(unit, *, iostat=ios) n_file
        if (ios /= 0) then
            close(unit)
            error stop 'Could not read random generator state size'
        end if
        if (n_file /= n) then
            close(unit)
            error stop 'Random generator state size does not match compiler'
        end if

        allocate(state(n))
        read(unit, *, iostat=ios) state
        if (ios /= 0) then
            close(unit)
            error stop 'Could not read random generator state'
        end if
        close(unit, iostat=close_ios)
        if (close_ios /= 0) error stop 'Could not close random generator state file'

        call random_seed(put=state)
        write(*,'(3A)') 'Resumed random generator state from "', trim(state_filename), '".'
    end subroutine initialize_random_seed_from_state

    subroutine save_random_seed_state(state_filename)
        character(len=*), intent(in) :: state_filename

        integer, allocatable :: state(:)
        integer :: n, unit, ios, close_ios

        call random_seed(size=n)
        allocate(state(n))
        call random_seed(get=state)

        open(newunit=unit, file=trim(state_filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) error stop 'Could not open random generator state file for writing'

        write(unit, *, iostat=ios) n
        if (ios == 0) write(unit, *, iostat=ios) state
        if (ios /= 0) then
            close(unit)
            error stop 'Could not write random generator state'
        end if

        close(unit, iostat=close_ios)
        if (close_ios /= 0) error stop 'Could not close random generator state file'
    end subroutine save_random_seed_state

end module random_seed_io_module