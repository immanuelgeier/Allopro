module trajectory_settings_module

    use project_module, only: project_type
    use topology_module, only: topologies_type

    implicit none
    private

    public :: initialize_trajectory_settings
    public :: read_trajectory_settings
    public :: write_trajectory_settings
    public :: print_trajectoy_settings

contains


! ============================================================
! INITIALIZE TRAJECTORY SETTINGS
! ============================================================

subroutine initialize_trajectory_settings(project, topologies)

    implicit none

    type(project_type), intent(inout) :: project
    type(topologies_type), intent(in) :: topologies

    character(len=1024) :: filename
    character(len=32)   :: input

    integer :: choice
    integer :: ios

    ! --------------------------------------------------------
    ! Allocate arrays and initialize default values.
    ! --------------------------------------------------------

    call set_default_trajectory_settings(project, topologies)

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Trajectory settings initialization'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Create new settings'
    write(*,*) '2 = Load and edit existing settings'
    write(*,*) '3 = Load existing settings unchanged'
    write(*,*)

    write(*,'(A)', advance='no') 'Select option: '

    read(*,'(A)') input
    read(input,*,iostat=ios) choice

    if (ios /= 0) error stop 'Invalid trajectory settings option'

    select case (choice)

    case (1)

        ! Defaults have already been initialized.
        call edit_trajectory_settings(project)

    case (2)

        ! Load existing settings and edit them.
        write(*,'(A)', advance='no') &
            'Settings filename: '

        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No settings filename specified'
        end if

        call read_trajectory_settings(trim(filename), project)

        call edit_trajectory_settings(project)

    case (3)

        ! Load existing settings without modification.
        write(*,'(A)', advance='no') &
            'Settings filename: '

        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            error stop 'No settings filename specified'
        end if

        call read_trajectory_settings(trim(filename), project)

        call validate_trajectory_settings(project)

        call print_trajectory_settings(project)

        return

    case default

        error stop 'Invalid trajectory settings option'

    end select

    call validate_trajectory_settings(project)

    call print_trajectory_settings(project)

    ! --------------------------------------------------------
    ! Optionally save the settings.
    ! --------------------------------------------------------

    write(*,*)
    write(*,'(A)', advance='no') &
        'Save trajectory settings? (y/n): '

    read(*,'(A)') input

    select case (trim(adjustl(input)))

    case ('y', 'Y', 'yes', 'YES')

        write(*,'(A)', advance='no') &
            'Settings filename [trajectory_settings.dat]: '

        read(*,'(A)') filename

        if (len_trim(filename) == 0) then
            filename = 'trajectory_settings.dat'
        end if

        call write_trajectory_settings(trim(filename), project)

    case ('n', 'N', 'no', 'NO')

        continue

    case default

        error stop 'Invalid save-settings option'

    end select

end subroutine initialize_trajectory_settings


! ============================================================
! SET DEFAULT SETTINGS
! ============================================================

subroutine set_default_trajectory_settings(project, topologies)

    implicit none

    type(project_type), intent(inout) :: project
    type(topologies_type), intent(in) :: topologies

    integer :: i, j
    integer :: ndcd

    do i = 1, project%nstate

        do j = 1, project%state(i)%nrep

            ndcd = project%state(i)%rep(j)%ndcd

            ! Remove arrays if they have already been allocated.
            if (allocated(project%state(i)%rep(j)%ignore_first_per_dcd)) &
                deallocate(project%state(i)%rep(j)%ignore_first_per_dcd)

            if (allocated(project%state(i)%rep(j)%ignore_last_per_dcd)) &
                deallocate(project%state(i)%rep(j)%ignore_last_per_dcd)

            if (allocated(project%state(i)%rep(j)%nframes_per_dcd_raw)) &
                deallocate(project%state(i)%rep(j)%nframes_per_dcd_raw)

            if (allocated(project%state(i)%rep(j)%nframes_per_dcd_use)) &
                deallocate(project%state(i)%rep(j)%nframes_per_dcd_use)

            allocate(project%state(i)%rep(j)%ignore_first_per_dcd(ndcd))
            allocate(project%state(i)%rep(j)%ignore_last_per_dcd(ndcd))

            allocate(project%state(i)%rep(j)%nframes_per_dcd_raw(ndcd))
            allocate(project%state(i)%rep(j)%nframes_per_dcd_use(ndcd))

            ! Default trajectory selection.
            project%state(i)%rep(j)%ignore_first_per_dcd = 1
            project%state(i)%rep(j)%ignore_last_per_dcd  = 1

            project%state(i)%rep(j)%offset_start = 0
            project%state(i)%rep(j)%offset_end   = 0

            ! Default atom selection.
            project%state(i)%rep(j)%nlast_natom_ignore = 0

            project%state(i)%rep(j)%nfirst_natom_to_read = &
                topologies%state(i)%topo%natom

            ! Frame counts will be determined by analyze_project_frames.
            project%state(i)%rep(j)%nframes_per_dcd_raw = 0
            project%state(i)%rep(j)%nframes_per_dcd_use = 0

        end do

    end do

end subroutine set_default_trajectory_settings


! ============================================================
! EDIT SETTINGS
! ============================================================

subroutine edit_trajectory_settings(project)

    implicit none

    type(project_type), intent(inout) :: project

    integer :: i, j, k
    character(len=128) :: label

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Edit trajectory settings'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) 'Press ENTER to keep the displayed value.'
    write(*,*)

    do i = 1, project%nstate

        write(*,*)
        write(*,*) 'State: ', i

        do j = 1, project%state(i)%nrep

            write(*,*)
            write(*,*) 'Repetition: ', j
            write(*,*) 'Number of DCD files: ', &
                project%state(i)%rep(j)%ndcd

            write(*,*)
            write(*,*) 'Trajectory offsets'

            call edit_integer('Start offset', &
                project%state(i)%rep(j)%offset_start)

            call edit_integer('End offset', &
                project%state(i)%rep(j)%offset_end)

            write(*,*)
            write(*,*) 'Atom selection'

            call edit_integer('Number of atoms to ignore at end', &
                project%state(i)%rep(j)%nlast_natom_ignore)

            call edit_integer('Number of atoms to read', &
                project%state(i)%rep(j)%nfirst_natom_to_read)

            write(*,*)
            write(*,*) 'Frame selection per DCD file'

            do k = 1, project%state(i)%rep(j)%ndcd

                write(*,*)
                write(*,*) 'DCD file: ', k

                write(label,'(A,I0,A)') &
                    'Ignore first frames, DCD ', k, ''

                call edit_integer(trim(label), &
                    project%state(i)%rep(j)%ignore_first_per_dcd(k))

                write(label,'(A,I0,A)') &
                    'Ignore last frames, DCD ', k, ''

                call edit_integer(trim(label), &
                    project%state(i)%rep(j)%ignore_last_per_dcd(k))

            end do

        end do

    end do

end subroutine edit_trajectory_settings


! ============================================================
! EDIT ONE INTEGER VALUE
! ============================================================

subroutine edit_integer(label, value)

    implicit none

    character(len=*), intent(in) :: label
    integer, intent(inout) :: value

    character(len=1024) :: input
    integer :: new_value
    integer :: ios

    write(*,'(A,I0,A)', advance='no') &
        trim(label)//' [', value, ']: '

    read(*,'(A)') input

    ! Empty input keeps the current value.
    if (len_trim(input) == 0) return

    read(input,*,iostat=ios) new_value

    if (ios /= 0) then
        error stop 'Invalid integer value'
    end if

    value = new_value

end subroutine edit_integer


! ============================================================
! VALIDATE SETTINGS
! ============================================================

subroutine validate_trajectory_settings(project)

    implicit none

    type(project_type), intent(in) :: project

    integer :: i, j, k

    do i = 1, project%nstate

        do j = 1, project%state(i)%nrep

            if (.not. allocated( &
                project%state(i)%rep(j)%ignore_first_per_dcd)) then
                error stop 'ignore_first_per_dcd is not allocated'
            end if

            if (.not. allocated( &
                project%state(i)%rep(j)%ignore_last_per_dcd)) then
                error stop 'ignore_last_per_dcd is not allocated'
            end if

            if (size(project%state(i)%rep(j)%ignore_first_per_dcd) /= &
                project%state(i)%rep(j)%ndcd) then
                error stop 'Invalid ignore_first_per_dcd size'
            end if

            if (size(project%state(i)%rep(j)%ignore_last_per_dcd) /= &
                project%state(i)%rep(j)%ndcd) then
                error stop 'Invalid ignore_last_per_dcd size'
            end if

            if (project%state(i)%rep(j)%offset_start < 0) then
                error stop 'offset_start must not be negative'
            end if

            if (project%state(i)%rep(j)%offset_end < 0) then
                error stop 'offset_end must not be negative'
            end if

            if (project%state(i)%rep(j)%nlast_natom_ignore < 0) then
                error stop 'nlast_natom_ignore must not be negative'
            end if

            if (project%state(i)%rep(j)%nfirst_natom_to_read <= 0) then
                error stop 'nfirst_natom_to_read must be positive'
            end if

            do k = 1, project%state(i)%rep(j)%ndcd

                if (project%state(i)%rep(j)%ignore_first_per_dcd(k) < 0) then
                    error stop 'ignore_first_per_dcd must not be negative'
                end if

                if (project%state(i)%rep(j)%ignore_last_per_dcd(k) < 0) then
                    error stop 'ignore_last_per_dcd must not be negative'
                end if

            end do

        end do

    end do

end subroutine validate_trajectory_settings


! ============================================================
! WRITE SETTINGS FILE
! ============================================================

subroutine write_trajectory_settings(filename, project)

    implicit none

    character(len=*), intent(in) :: filename
    type(project_type), intent(in) :: project

    integer :: unit
    integer :: ios
    integer :: i, j

    open(newunit=unit, file=trim(filename), &
         status='replace', action='write', iostat=ios)

    if (ios /= 0) then
        error stop 'Could not open trajectory settings file for writing'
    end if

    ! File format version.
    write(unit,*) 1

    ! Number of states.
    write(unit,*) project%nstate

    do i = 1, project%nstate

        ! Number of repetitions in this state.
        write(unit,*) project%state(i)%nrep

        do j = 1, project%state(i)%nrep

            ! Number of DCD files.
            write(unit,*) project%state(i)%rep(j)%ndcd

            ! Scalar settings.
            write(unit,*) &
                project%state(i)%rep(j)%offset_start, &
                project%state(i)%rep(j)%offset_end, &
                project%state(i)%rep(j)%nlast_natom_ignore, &
                project%state(i)%rep(j)%nfirst_natom_to_read

            ! Per-DCD settings.
            write(unit,*) &
                project%state(i)%rep(j)%ignore_first_per_dcd

            write(unit,*) &
                project%state(i)%rep(j)%ignore_last_per_dcd

        end do

    end do

    close(unit)

    write(*,*)
    write(*,*) 'Trajectory settings saved to: ', trim(filename)

end subroutine write_trajectory_settings


! ============================================================
! READ SETTINGS FILE
! ============================================================

subroutine read_trajectory_settings(filename, project)

    implicit none

    character(len=*), intent(in) :: filename
    type(project_type), intent(inout) :: project

    integer :: unit
    integer :: ios
    integer :: version
    integer :: nstate_file
    integer :: nrep_file
    integer :: ndcd_file
    integer :: i, j

    open(newunit=unit, file=trim(filename), &
         status='old', action='read', iostat=ios)

    if (ios /= 0) then
        error stop 'Could not open trajectory settings file'
    end if

    read(unit,*,iostat=ios) version

    if (ios /= 0) then
        close(unit)
        error stop 'Could not read trajectory settings version'
    end if

    if (version /= 1) then
        close(unit)
        error stop 'Unsupported trajectory settings file version'
    end if

    read(unit,*,iostat=ios) nstate_file

    if (ios /= 0) then
        close(unit)
        error stop 'Could not read number of states'
    end if

    if (nstate_file /= project%nstate) then
        close(unit)
        error stop 'Settings do not match project state count'
    end if

    do i = 1, project%nstate

        read(unit,*,iostat=ios) nrep_file

        if (ios /= 0) then
            close(unit)
            error stop 'Could not read number of repetitions'
        end if

        if (nrep_file /= project%state(i)%nrep) then
            close(unit)
            error stop 'Settings do not match project repetition count'
        end if

        do j = 1, project%state(i)%nrep

            read(unit,*,iostat=ios) ndcd_file

            if (ios /= 0) then
                close(unit)
                error stop 'Could not read number of DCD files'
            end if

            if (ndcd_file /= project%state(i)%rep(j)%ndcd) then
                close(unit)
                error stop 'Settings do not match DCD file count'
            end if

            read(unit,*,iostat=ios) &
                project%state(i)%rep(j)%offset_start, &
                project%state(i)%rep(j)%offset_end, &
                project%state(i)%rep(j)%nlast_natom_ignore, &
                project%state(i)%rep(j)%nfirst_natom_to_read

            if (ios /= 0) then
                close(unit)
                error stop 'Could not read trajectory scalar settings'
            end if

            read(unit,*,iostat=ios) &
                project%state(i)%rep(j)%ignore_first_per_dcd

            if (ios /= 0) then
                close(unit)
                error stop 'Could not read ignore_first_per_dcd'
            end if

            read(unit,*,iostat=ios) &
                project%state(i)%rep(j)%ignore_last_per_dcd

            if (ios /= 0) then
                close(unit)
                error stop 'Could not read ignore_last_per_dcd'
            end if

        end do

    end do

    close(unit)

    call validate_trajectory_settings(project)

    write(*,*)
    write(*,*) 'Trajectory settings loaded from: ', trim(filename)

end subroutine read_trajectory_settings


! ============================================================
! PRINT SETTINGS
! ============================================================

subroutine print_trajectory_settings(project)

    implicit none

    type(project_type), intent(in) :: project

    integer :: i, j, k

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Current trajectory settings'
    write(*,*) '======================================'

    do i = 1, project%nstate

        write(*,*)
        write(*,*) 'State: ', i

        do j = 1, project%state(i)%nrep

            write(*,*)
            write(*,*) 'Repetition: ', j

            write(*,*) 'offset_start: ', &
                project%state(i)%rep(j)%offset_start

            write(*,*) 'offset_end: ', &
                project%state(i)%rep(j)%offset_end

            write(*,*) 'nlast_natom_ignore: ', &
                project%state(i)%rep(j)%nlast_natom_ignore

            write(*,*) 'nfirst_natom_to_read: ', &
                project%state(i)%rep(j)%nfirst_natom_to_read

            do k = 1, project%state(i)%rep(j)%ndcd

                write(*,*) 'DCD ', k, &
                    ' ignore_first = ', &
                    project%state(i)%rep(j)%ignore_first_per_dcd(k), &
                    ' ignore_last = ', &
                    project%state(i)%rep(j)%ignore_last_per_dcd(k)

            end do

        end do

    end do

    write(*,*) '======================================'

end subroutine print_trajectory_settings


end module trajectory_settings_module