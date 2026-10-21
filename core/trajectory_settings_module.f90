module trajectory_settings_module

    use project_module, only: project_type
    use topology_module, only: topology_set_type

    implicit none
    private

    public :: initialize_trajectory_settings
    public :: read_trajectory_settings
    public :: write_trajectory_settings
    public :: print_trajectory_settings
    public :: validate_trajectory_settings,set_default_trajectory_settings
contains


! ============================================================
! INITIALIZE TRAJECTORY SETTINGS
! ============================================================

subroutine initialize_trajectory_settings(project, topologies)

    implicit none

    type(project_type), intent(inout) :: project
    type(topology_set_type), intent(in) :: topologies

    character(len=1024) :: filename
    character(len=32)    :: input

    integer :: choice
    integer :: ios

    ! --------------------------------------------------------
    ! Allocate arrays and initialize default values.
    ! --------------------------------------------------------

    call set_default_trajectory_settings(project, topologies)

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    Trajectory settings initialization'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) '1 = Create new settings'
    write(*,*) '2 = Load and edit existing settings'
    write(*,*) '3 = Load existing settings unchanged'
    write(*,*) '4 = Skip'
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

    case (4)

        ! Skip initialization.
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
    type(topology_set_type), intent(in) :: topologies

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
    integer :: frame_mode
    integer :: ios
    integer :: ndcd

    integer :: ref_ignore_first
    integer :: ref_ignore_last
    
    ! Variablen für globale Offsets/Selektionen in Modus 3
    integer :: ref_offset_start
    integer :: ref_offset_end
    integer :: ref_nlast_ignore
    integer :: ref_nfirst_read

    character(len=128)  :: label
    character(len=1024) :: input

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    Edit trajectory settings'
    write(*,*) '======================================'
    write(*,*)
    write(*,*) 'Press ENTER to keep the displayed value.'
    write(*,*)

    ! --------------------------------------------------------
    ! Select frame settings mode.
    ! --------------------------------------------------------

    write(*,*)
    write(*,*) 'Frame selection mode:'
    write(*,*) '1 = Set each DCD individually'
    write(*,*) '2 = Copy first DCD within each repetition'
    write(*,*) '3 = Copy first DCD to ALL repetitions and states'
    write(*,*)

    frame_mode = 1

    write(*,'(A)', advance='no') 'Select mode [1]: '
    read(*,'(A)') input

    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) frame_mode
        if (ios /= 0) then
            error stop 'Invalid frame selection mode'
        end if
    end if

    if (frame_mode < 1 .or. frame_mode > 3) then
        error stop 'Invalid frame selection mode'
    end if

    ! --------------------------------------------------------
    ! Mode 3: Edit global reference (State 1, Rep 1) once.
    ! --------------------------------------------------------

    if (frame_mode == 3) then

        if (project%nstate < 1) then
            error stop 'No states available'
        end if

        if (project%state(1)%nrep < 1) then
            error stop 'No repetitions available in state 1'
        end if

        if (project%state(1)%rep(1)%ndcd < 1) then
            error stop 'No DCD files available in reference repetition'
        end if

        write(*,*)
        write(*,*) '======================================'
        write(*,*) ' Global reference: State 1, Rep 1'
        write(*,*) '======================================'

        ! Offsets & Atom selection global edit
        call edit_integer('Start offset', project%state(1)%rep(1)%offset_start)
        call edit_integer('End offset', project%state(1)%rep(1)%offset_end)

        call edit_integer('Number of atoms to ignore at end', project%state(1)%rep(1)%nlast_natom_ignore)
        call edit_integer('Number of atoms to read', project%state(1)%rep(1)%nfirst_natom_to_read)

        ! Frame settings global edit
        call edit_integer('Ignore first frames', project%state(1)%rep(1)%ignore_first_per_dcd(1))
        call edit_integer('Ignore last frames', project%state(1)%rep(1)%ignore_last_per_dcd(1))

        ref_ignore_first = project%state(1)%rep(1)%ignore_first_per_dcd(1)
        ref_ignore_last  = project%state(1)%rep(1)%ignore_last_per_dcd(1)
        
        ref_offset_start = project%state(1)%rep(1)%offset_start
        ref_offset_end   = project%state(1)%rep(1)%offset_end
        ref_nlast_ignore = project%state(1)%rep(1)%nlast_natom_ignore
        ref_nfirst_read  = project%state(1)%rep(1)%nfirst_natom_to_read

        ! Apply immediately to all states and repetitions
        do i = 1, project%nstate
            do j = 1, project%state(i)%nrep
                ndcd = project%state(i)%rep(j)%ndcd
                if (ndcd <= 0) cycle

                project%state(i)%rep(j)%offset_start        = ref_offset_start
                project%state(i)%rep(j)%offset_end          = ref_offset_end
                project%state(i)%rep(j)%nlast_natom_ignore  = ref_nlast_ignore
                project%state(i)%rep(j)%nfirst_natom_to_read = ref_nfirst_read

                project%state(i)%rep(j)%ignore_first_per_dcd = ref_ignore_first
                project%state(i)%rep(j)%ignore_last_per_dcd  = ref_ignore_last
            end do
        end do

        write(*,*)
        write(*,*) 'Global settings applied successfully to all states and repetitions.'
        return

    end if

    ! --------------------------------------------------------
    ! Edit offsets and atom selection for Modes 1 and 2.
    ! (Only runs if frame_mode /= 3)
    ! --------------------------------------------------------

    do i = 1, project%nstate

        write(*,*)
        write(*,*) 'State: ', i

        do j = 1, project%state(i)%nrep

            ndcd = project%state(i)%rep(j)%ndcd

            write(*,*)
            write(*,*) 'Repetition: ', j
            write(*,*) 'Number of DCD files: ', ndcd

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

            ! -----------------------------------------------
            ! Frame settings for modes 1 and 2.
            ! -----------------------------------------------

            select case (frame_mode)

            case (1)

                do k = 1, ndcd

                    write(*,*)
                    write(*,*) 'DCD file: ', k

                    write(label,'(A,I0)') &
                        'Ignore first frames, DCD ', k

                    call edit_integer(trim(label), &
                        project%state(i)%rep(j)%ignore_first_per_dcd(k))

                    write(label,'(A,I0)') &
                        'Ignore last frames, DCD ', k

                    call edit_integer(trim(label), &
                        project%state(i)%rep(j)%ignore_last_per_dcd(k))

                end do

            case (2)

                if (ndcd > 0) then

                    write(*,*)
                    write(*,*) 'Reference settings: DCD 1'

                    call edit_integer('Ignore first frames', &
                        project%state(i)%rep(j)%ignore_first_per_dcd(1))

                    call edit_integer('Ignore last frames', &
                        project%state(i)%rep(j)%ignore_last_per_dcd(1))

                    if (ndcd > 1) then

                        project%state(i)%rep(j)%ignore_first_per_dcd(2:ndcd) = &
                            project%state(i)%rep(j)%ignore_first_per_dcd(1)

                        project%state(i)%rep(j)%ignore_last_per_dcd(2:ndcd) = &
                            project%state(i)%rep(j)%ignore_last_per_dcd(1)

                    end if

                end if

            end select

        end do

    end do

end subroutine edit_trajectory_settings
!====================
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
    integer :: i, j, k

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

            ! Per-DCD settings (element by element or as a clean row)
            do k = 1, project%state(i)%rep(j)%ndcd
                write(unit,*) project%state(i)%rep(j)%ignore_first_per_dcd(k), &
                              project%state(i)%rep(j)%ignore_last_per_dcd(k)
            end do

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
    integer :: i, j, k

    open(newunit=unit, file=trim(filename), &
         status='old', action='read', iostat=ios)

    if (ios /= 0) then
        error stop 'Could not open trajectory settings file'
    end if

    read(unit,*,iostat=ios) version

    if (ios /= 0 .or. version /= 1) then
        close(unit)
        error stop 'Unsupported or invalid trajectory settings file version'
    end if

    read(unit,*,iostat=ios) nstate_file
    if (ios /= 0 .or. nstate_file /= project%nstate) then
        close(unit)
        error stop 'Settings do not match project state count'
    end if

    do i = 1, project%nstate

        read(unit,*,iostat=ios) nrep_file
        if (ios /= 0 .or. nrep_file /= project%state(i)%nrep) then
            close(unit)
            error stop 'Settings do not match project repetition count'
        end if

        do j = 1, project%state(i)%nrep

            read(unit,*,iostat=ios) ndcd_file
            if (ios /= 0 .or. ndcd_file /= project%state(i)%rep(j)%ndcd) then
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

            ! Per-DCD settings einlesen
            do k = 1, project%state(i)%rep(j)%ndcd
                read(unit,*,iostat=ios) &
                    project%state(i)%rep(j)%ignore_first_per_dcd(k), &
                    project%state(i)%rep(j)%ignore_last_per_dcd(k)

                if (ios /= 0) then
                    close(unit)
                    error stop 'Could not read per-DCD ignore settings'
                end if
            end do

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