module random_seed_io_module

    use iso_fortran_env, only: int64
    implicit none

    private

    public :: random_seed_config_type,configure_seed
    public :: initialize_random_seed

    type :: random_seed_config_type
        logical :: use_saved_seed = .false.
        logical :: save_seed      = .true.

        character(len=1024) :: seed_file = 'random_seed.dat'

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
        integer :: clock
        integer :: unit
        integer :: ios
        integer :: n_file

        integer(int64) :: clock64

        logical :: file_exists

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

end module random_seed_io_module