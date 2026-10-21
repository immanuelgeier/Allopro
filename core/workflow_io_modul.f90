module workflow_io_module

    use workflow_module, only: library_workflow_type

    implicit none
    private

    public :: init_library_workflow
    public :: write_library_workflow
    public :: read_library_workflow
    public :: print_library_workflow

contains

    subroutine init_library_workflow(job, mode, filename)

        type(library_workflow_type), intent(out) :: job
        character(len=*), intent(in) :: mode
        character(len=*), intent(in), optional :: filename

        ! Standardwerte der Workflow-Typen uebernehmen.
        job = library_workflow_type()

        select case (trim(mode))

        case ('manual')

            ! Die Standardwerte koennen nach dem Aufruf
            ! vom aufrufenden Programm angepasst werden.

        case ('file')

            if (.not. present(filename)) then
                error stop 'init_library_workflow: filename missing'
            end if

            call read_library_workflow(job, filename)

        case default

            error stop 'init_library_workflow: unknown mode'

        end select

    end subroutine init_library_workflow


    subroutine write_library_workflow(job, filename)

        type(library_workflow_type), intent(in) :: job
        character(len=*), intent(in) :: filename

        ! Hier werden alle Parameter unter ihren Namen
        ! in die Konfigurationsdatei geschrieben.

        ! Implementierung folgt.

    end subroutine write_library_workflow


    subroutine read_library_workflow(job, filename)

        type(library_workflow_type), intent(inout) :: job
        character(len=*), intent(in) :: filename

        ! Hier werden die gespeicherten Parameter eingelesen.
        ! Nicht angegebene Parameter behalten ihre Standardwerte.

        ! Implementierung folgt.

    end subroutine read_library_workflow


    subroutine print_library_workflow(job)

        type(library_workflow_type), intent(in) :: job

        ! Alle Einstellungen lesbar auf stdout ausgeben.

        ! Implementierung folgt.

    end subroutine print_library_workflow

end module workflow_io_module