module centroid_distance_histogram_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use workflow_module, only: workflow_config_filename, centroid_distance_histogram_workflow_config_stem
    use lifting_plan_module, only: lifting_plan_type, read_lifting_plan
    use centroid_distance_histogram_module, only: run_centroid_distance_histogram
    implicit none
    private

    public :: centroid_distance_histogram_workflow_type
    public :: configure_centroid_distance_histogram_workflow
    public :: read_centroid_distance_histogram_workflow
    public :: write_centroid_distance_histogram_workflow
    public :: validate_centroid_distance_histogram_workflow
    public :: run_centroid_distance_histogram_workflow

    type :: centroid_distance_histogram_workflow_type
        character(len=1024) :: lifting_plan_filename = 'lifting_plan.dat'
        character(len=1024) :: centroid_path = '.'
        character(len=256) :: centroid_prefix = 'centroid'
        character(len=1024) :: histogram_filename = 'centroid_distance_histogram.dat'
    end type centroid_distance_histogram_workflow_type

contains

subroutine configure_centroid_distance_histogram_workflow(job)
    type(centroid_distance_histogram_workflow_type), intent(out) :: job
    character(len=1024) :: input, config_filename
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) ' Centroid Distance Histogram Workflow'
    write(*,*) '======================================'
    write(*,*) '1 = Create new workflow'
    write(*,*) '2 = Load and edit existing workflow'
    write(*,*) '3 = Load existing workflow unchanged'
    write(*,*) '4 = Skip workflow configuration'
    write(*,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid centroid histogram workflow option'

    select case (choice)
    case (1)
        job = centroid_distance_histogram_workflow_type()
        call edit_centroid_distance_histogram_workflow(job)
    case (2)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_centroid_distance_histogram_workflow(trim(config_filename), job)
        call edit_centroid_distance_histogram_workflow(job)
    case (3)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_centroid_distance_histogram_workflow(trim(config_filename), job)
        call validate_centroid_distance_histogram_workflow(job)
        call print_centroid_distance_histogram_workflow(job)
    case (4)
        job = centroid_distance_histogram_workflow_type()
        write(*,*) 'Centroid histogram workflow configuration skipped.'
        return
    case default
        error stop 'Invalid centroid histogram workflow option'
    end select

    call validate_centroid_distance_histogram_workflow(job)
    write(*,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(*,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        config_filename = workflow_config_filename(trim(input), &
            centroid_distance_histogram_workflow_config_stem)
        call write_centroid_distance_histogram_workflow(job, trim(config_filename))
        write(*,'(3A)') 'Configuration saved to "', trim(config_filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_centroid_distance_histogram_workflow

subroutine edit_centroid_distance_histogram_workflow(job)
    type(centroid_distance_histogram_workflow_type), intent(inout) :: job
    character(len=1024) :: input

    write(*,*)
    write(*,*) 'Press Enter to keep the current value.'
    write(*,'(A,A,A)', advance='no') 'Lifting plan filename [', &
        trim(job%lifting_plan_filename), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) job%lifting_plan_filename = trim(input)

    write(*,'(A,A,A)', advance='no') 'Centroid trajectory directory [', &
        trim(job%centroid_path), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) job%centroid_path = trim(input)

    write(*,'(A,A,A)', advance='no') 'Centroid trajectory filename prefix [', &
        trim(job%centroid_prefix), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) job%centroid_prefix = trim(input)

    write(*,'(A,A,A)', advance='no') 'Histogram output filename [', &
        trim(job%histogram_filename), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) job%histogram_filename = trim(input)

    call validate_centroid_distance_histogram_workflow(job)
    call print_centroid_distance_histogram_workflow(job)
end subroutine edit_centroid_distance_histogram_workflow

subroutine validate_centroid_distance_histogram_workflow(job)
    type(centroid_distance_histogram_workflow_type), intent(in) :: job

    if (len_trim(job%lifting_plan_filename) == 0) error stop 'Lifting plan filename cannot be empty'
    if (len_trim(job%centroid_path) == 0) error stop 'Centroid trajectory path cannot be empty'
    if (len_trim(job%centroid_prefix) == 0) error stop 'Centroid trajectory prefix cannot be empty'
    if (len_trim(job%histogram_filename) == 0) error stop 'Histogram output filename cannot be empty'
end subroutine validate_centroid_distance_histogram_workflow

subroutine run_centroid_distance_histogram_workflow(job)
    type(centroid_distance_histogram_workflow_type), intent(in) :: job
    type(lifting_plan_type) :: lifting_plan
    character(len=2048) :: centroid_prefix

    call validate_centroid_distance_histogram_workflow(job)
    call read_lifting_plan(lifting_plan, trim(job%lifting_plan_filename))
    centroid_prefix = trim(job%centroid_path)
    if (len_trim(centroid_prefix) > 0) then
        if (centroid_prefix(len_trim(centroid_prefix):len_trim(centroid_prefix)) /= '/') then
            centroid_prefix = trim(centroid_prefix)//'/'
        end if
    end if
    centroid_prefix = trim(centroid_prefix)//trim(job%centroid_prefix)
    call run_centroid_distance_histogram(lifting_plan, trim(centroid_prefix), &
        trim(job%histogram_filename))
end subroutine run_centroid_distance_histogram_workflow

subroutine write_centroid_distance_histogram_workflow(job, filename)
    type(centroid_distance_histogram_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_centroid_distance_histogram_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write centroid histogram workflow: ', trim(filename)
        error stop 'Centroid histogram workflow output failed'
    end if
    write(unit,'(A)') '# Centroid distance histogram workflow'
    write(unit,'(A)') 'workflow_type = centroid_distance_histogram'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,A)') 'lifting_plan_filename = ', trim(job%lifting_plan_filename)
    write(unit,'(A,A)') 'centroid_path = ', trim(job%centroid_path)
    write(unit,'(A,A)') 'centroid_prefix = ', trim(job%centroid_prefix)
    write(unit,'(A,A)') 'histogram_filename = ', trim(job%histogram_filename)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing centroid histogram workflow file'
end subroutine write_centroid_distance_histogram_workflow

subroutine read_centroid_distance_histogram_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(centroid_distance_histogram_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1900) :: value
    integer :: unit, ios, equals_position, format_version
    logical :: saw_type, saw_version

    job = centroid_distance_histogram_workflow_type()
    saw_type = .false.
    saw_version = .false.
    format_version = -1
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read centroid histogram workflow: ', trim(filename)
        error stop 'Centroid histogram workflow input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) call config_read_error(unit, 'Error reading centroid histogram workflow')
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) call config_read_error(unit, 'Invalid workflow configuration line')
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'centroid_distance_histogram') then
                call config_read_error(unit, 'Unexpected workflow_type')
            end if
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0) call config_read_error(unit, 'Invalid format_version')
            if (format_version /= 1) call config_read_error(unit, 'Unsupported format_version')
            saw_version = .true.
        case ('lifting_plan_filename')
            job%lifting_plan_filename = trim(value)
        case ('centroid_path')
            job%centroid_path = trim(value)
        case ('centroid_prefix')
            job%centroid_prefix = trim(value)
        case ('histogram_filename')
            job%histogram_filename = trim(value)
        case default
            call config_read_error(unit, 'Unknown workflow parameter: '//trim(key))
        end select
    end do
    close(unit)
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete centroid histogram workflow file'
    call validate_centroid_distance_histogram_workflow(job)
end subroutine read_centroid_distance_histogram_workflow

subroutine print_centroid_distance_histogram_workflow(job)
    type(centroid_distance_histogram_workflow_type), intent(in) :: job
    write(output_unit,'(A,A)') 'Lifting plan: ', trim(job%lifting_plan_filename)
    write(output_unit,'(A,A)') 'Centroid path: ', trim(job%centroid_path)
    write(output_unit,'(A,A)') 'Centroid filename prefix: ', trim(job%centroid_prefix)
    write(output_unit,'(A,A)') 'Histogram output: ', trim(job%histogram_filename)
end subroutine print_centroid_distance_histogram_workflow

subroutine config_read_error(unit, message)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    close(unit)
    write(error_unit,'(A)') trim(message)
    error stop 'Centroid histogram workflow configuration error'
end subroutine config_read_error

end module centroid_distance_histogram_workflow_module
