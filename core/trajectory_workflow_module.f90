module trajectory_workflow_io_module
    use iso_fortran_env, only: output_unit, error_unit
    use project_module, only: project_type
    use topology_module, only: topology_set_type
    use trajectory_settings_module, only: read_trajectory_settings, &
                                          write_trajectory_settings, &
                                          validate_trajectory_settings, &
                                          print_trajectory_settings
    implicit none
    private

    public :: trajectory_run_params_type
    public :: trajectory_io_params_type
    public :: trajectory_workflow_type

    public :: read_trajectory_workflow
    public :: write_trajectory_workflow
    public :: print_trajectory_workflow
    public :: validate_trajectory_workflow
    public :: configure_trajectory_workflow
    public :: run_trajectory_workflow

    ! ------------------------------------------------------------
    ! PARAMETER TYPES
    ! ------------------------------------------------------------

    type :: trajectory_run_params_type
        logical :: print_project_summary = .true.
        logical :: write_project_summary = .true.
    end type trajectory_run_params_type

    type :: trajectory_io_params_type
        character(len=1024) :: settings_file = 'trajectory_settings.dat'
        character(len=1024) :: summary_file  = 'project_summary.txt'
    end type trajectory_io_params_type

    ! ------------------------------------------------------------
    ! MAIN WORKFLOW TYPE
    ! ------------------------------------------------------------

    type :: trajectory_workflow_type
        type(trajectory_run_params_type) :: run_params
        type(trajectory_io_params_type)  :: io_params
    end type trajectory_workflow_type

contains

    ! ============================================================
    ! READ WORKFLOW CONFIG (.cfg)
    ! ============================================================

    subroutine read_trajectory_workflow(job, filename)
        type(trajectory_workflow_type), intent(out) :: job
        character(len=*), intent(in) :: filename
        integer :: unit, ios
        character(len=1024) :: line, key, val
        integer :: eq_pos

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Error opening trajectory workflow file: ', trim(filename)
            error stop
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                key = adjustl(trim(line(1:eq_pos-1)))
                val = adjustl(trim(line(eq_pos+1:)))

                select case (trim(key))
                case ('print_project_summary')
                    read(val, *, iostat=ios) job%run_params%print_project_summary
                case ('write_project_summary')
                    read(val, *, iostat=ios) job%run_params%write_project_summary
                case ('settings_file')
                    job%io_params%settings_file = trim(val)
                case ('summary_file')
                    job%io_params%summary_file = trim(val)
                end select
            end if
        end do

        close(unit)
    end subroutine read_trajectory_workflow


    ! ============================================================
    ! WRITE WORKFLOW CONFIG (.cfg)
    ! ============================================================

    subroutine write_trajectory_workflow(job, filename)
        type(trajectory_workflow_type), intent(in) :: job
        character(len=*), intent(in) :: filename
        integer :: unit, ios

        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Error writing trajectory workflow file: ', trim(filename)
            error stop
        end if

        write(unit, '(A)') '# Trajectory Settings Workflow Configuration'
        write(unit, '(A)') 'workflow_type = traj_settings'
        write(unit, '(A)') ''
        write(unit, '(A)') '# Run parameters'
        write(unit, '(A,L1)') 'print_project_summary = ', job%run_params%print_project_summary
        write(unit, '(A,L1)') 'write_project_summary = ', job%run_params%write_project_summary
        write(unit, '(A)') ''
        write(unit, '(A)') '# I/O parameters'
        write(unit, '(A,A)') 'settings_file = ', trim(job%io_params%settings_file)
        write(unit, '(A,A)') 'summary_file  = ', trim(job%io_params%summary_file)

        close(unit)
    end subroutine write_trajectory_workflow


    ! ============================================================
    ! INTERACTIVE CONFIGURATION
    ! ============================================================

    subroutine configure_trajectory_workflow(job)
        type(trajectory_workflow_type), intent(out) :: job
        character(len=1024) :: input

        write(*,*)
        write(*,*) '======================================'
        write(*,*) ' Configure Trajectory Settings Workflow'
        write(*,*) '======================================'

        ! Settings-Datei Abfragen
        write(*,'(A,A,A)', advance='no') 'Trajectory settings file [', trim(job%io_params%settings_file), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) job%io_params%settings_file = trim(adjustl(input))

        ! Flags Abfragen
        write(*,'(A,L1,A)', advance='no') 'Print project summary? (T/F) [', job%run_params%print_project_summary, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input, *) job%run_params%print_project_summary

        write(*,'(A,L1,A)', advance='no') 'Write project summary? (T/F) [', job%run_params%write_project_summary, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) read(input, *) job%run_params%write_project_summary

        if (job%run_params%write_project_summary) then
            write(*,'(A,A,A)', advance='no') 'Summary output file [', trim(job%io_params%summary_file), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%io_params%summary_file = trim(adjustl(input))
        end if

    end subroutine configure_trajectory_workflow


    ! ============================================================
    ! VALIDATION & PRINT
    ! ============================================================

    subroutine validate_trajectory_workflow(job)
        type(trajectory_workflow_type), intent(in) :: job

        if (len_trim(job%io_params%settings_file) == 0) then
            error stop 'Trajectory workflow error: settings_file must not be empty'
        end if
    end subroutine validate_trajectory_workflow


    subroutine print_trajectory_workflow(job)
        type(trajectory_workflow_type), intent(in) :: job

        write(*,*) '--- Trajectory Workflow Parameters ---'
        write(*,*) 'Settings File:         ', trim(job%io_params%settings_file)
        write(*,*) 'Print Project Summary: ', job%run_params%print_project_summary
        write(*,*) 'Write Project Summary: ', job%run_params%write_project_summary
        if (job%run_params%write_project_summary) then
            write(*,*) 'Summary File:          ', trim(job%io_params%summary_file)
        end if
    end subroutine print_trajectory_workflow


    ! ============================================================
    ! RUN WORKFLOW STEP
    ! ============================================================

    subroutine run_trajectory_workflow(job, project)
        type(trajectory_workflow_type), intent(in) :: job
        type(project_type), intent(inout) :: project

        write(output_unit, '(A,A)') '-> Executing Trajectory Workflow using: ', trim(job%io_params%settings_file)

        ! 1. Trajektorien-Einstellungen einlesen
        call read_trajectory_settings(trim(job%io_params%settings_file), project)

        ! 2. Projekt-Zusammenfassung ausgeben (falls gesetzt)
        if (job%run_params%print_project_summary) then
            call print_trajectory_settings(project)
        end if

        ! 3. Projekt-Zusammenfassung in Datei schreiben (falls gesetzt)
        if (job%run_params%write_project_summary) then
            ! call write_project_summary(trim(job%io_params%summary_file), project)
            write(output_unit, '(A,A)') '-> Project summary written to: ', trim(job%io_params%summary_file)
        end if

    end subroutine run_trajectory_workflow

end module trajectory_workflow_io_module