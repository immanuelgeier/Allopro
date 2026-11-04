module c2_te_matrix_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use workflow_module, only: c2_te_matrix_workflow_type, workflow_config_filename, &
        c2_te_matrix_workflow_config_stem
    use direction_vectors_analyse_plan_module, only: cv_plan
    use cv_plan_io_module, only: read_cv_plan
    use lifting_plan_module, only: lifting_plan_type, read_lifting_plan
    use healpix_state_space_exploration, only: odf_lib_type, init_odf_lib
    use healpix_histogram_sampling_execution_module, only: run_calculate_te_matrix_chunked_C2
    implicit none
    private

    public :: c2_te_matrix_workflow_type
    public :: configure_c2_te_matrix_workflow
    public :: read_c2_te_matrix_workflow
    public :: write_c2_te_matrix_workflow
    public :: validate_c2_te_matrix_workflow
    public :: run_c2_te_matrix_workflow
    public :: check_c2_te_matrix_workflow_inputs

contains

subroutine configure_c2_te_matrix_workflow(job)
    type(c2_te_matrix_workflow_type), intent(out) :: job
    character(len=1024) :: input, config_filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' C2 Transfer-Entropy Matrix Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid C2 TE matrix workflow option'

    select case (choice)
    case (1)
        job = c2_te_matrix_workflow_type()
        call edit_c2_te_matrix_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_c2_te_matrix_workflow(trim(config_filename), job)
        if (choice == 2) then
            call edit_c2_te_matrix_workflow(job)
        else
            call print_c2_te_matrix_workflow(job)
        end if
    case (4)
        job = c2_te_matrix_workflow_type()
        write(output_unit,*) 'C2 TE matrix workflow configuration skipped.'
        return
    case default
        error stop 'Invalid C2 TE matrix workflow option'
    end select

    call validate_c2_te_matrix_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        config_filename = workflow_config_filename(trim(input), c2_te_matrix_workflow_config_stem)
        call write_c2_te_matrix_workflow(job, trim(config_filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(config_filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_c2_te_matrix_workflow

subroutine edit_c2_te_matrix_workflow(job)
    type(c2_te_matrix_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_string('CV plan filename', job%io_params%cv_plan_filename)
    call edit_string('Lifting plan filename (all listed states are processed)', &
        job%io_params%lifting_plan_filename)
    call edit_string('Active-pair matrix filename', job%io_params%pair_active_filename)
    call edit_string('Histogram directory', job%io_params%histogram_path)
    call edit_string('Histogram filename prefix', job%io_params%histogram_prefix)
    call edit_string('Output directory', job%io_params%output_path)
    call edit_string('Output filename prefix (result: prefix + te_matrix.txt)', &
        job%io_params%output_prefix)
    call edit_string('ODF directory', job%io_params%odf_path)
    call edit_string('ODF filename prefix', job%io_params%odf_prefix)
    call edit_integer('Histogram HEALPix levels', job%model_params%nlevel_target)
    call edit_integer('Maximum histogram chunk size', job%model_params%max_chunk_size)
    call edit_integer('ODF maximum lag (number of lag files minus one)', &
        job%model_params%odf_max_lag)
    call edit_integer('ODF basis index', job%model_params%odf_basis_index)
    call edit_integer('ODF HEALPix levels', job%model_params%odf_nlevel)
    call edit_logical('Read shuffled histograms (shuffle run 1)', job%model_params%shuffled)
    call edit_string('Output policy (ask/replace/rename_old)', job%io_params%output_policy)
    call validate_c2_te_matrix_workflow(job)
    call print_c2_te_matrix_workflow(job)

contains
    subroutine edit_integer(prompt, value)
        character(len=*), intent(in) :: prompt
        integer, intent(inout) :: value
        write(output_unit,'(A,A,I0,A)', advance='no') trim(prompt), ' [', value, ']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid integer workflow parameter'
    end subroutine edit_integer

    subroutine edit_logical(prompt, value)
        character(len=*), intent(in) :: prompt
        logical, intent(inout) :: value
        write(output_unit,'(A,L1,A)', advance='no') trim(prompt), value, ': '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid logical workflow parameter'
    end subroutine edit_logical

    subroutine edit_string(prompt, value)
        character(len=*), intent(in) :: prompt
        character(len=*), intent(inout) :: value
        write(output_unit,'(A,A,A)', advance='no') trim(prompt), ' [', trim(value), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) value = trim(input)
    end subroutine edit_string
end subroutine edit_c2_te_matrix_workflow

subroutine validate_c2_te_matrix_workflow(job)
    type(c2_te_matrix_workflow_type), intent(in) :: job

    if (job%model_params%nlevel_target < 1 .or. job%model_params%nlevel_target > 14) then
        error stop 'Histogram HEALPix levels must be between 1 and 14'
    end if
    if (job%model_params%max_chunk_size < 1) error stop 'Maximum chunk size must be positive'
    if (job%model_params%odf_max_lag < 0) error stop 'ODF maximum lag cannot be negative'
    if (job%model_params%odf_basis_index < 1) error stop 'ODF basis index must be positive'
    if (job%model_params%odf_nlevel < 1 .or. job%model_params%odf_nlevel > 14) then
        error stop 'ODF HEALPix levels must be between 1 and 14'
    end if
    if (len_trim(job%io_params%cv_plan_filename) == 0) error stop 'CV plan filename cannot be empty'
    if (len_trim(job%io_params%lifting_plan_filename) == 0) error stop 'Lifting plan filename cannot be empty'
    if (len_trim(job%io_params%pair_active_filename) == 0) error stop 'Pair matrix filename cannot be empty'
    if (len_trim(job%io_params%histogram_path) == 0) error stop 'Histogram path cannot be empty'
    if (len_trim(job%io_params%histogram_prefix) == 0) error stop 'Histogram prefix cannot be empty'
    if (len_trim(job%io_params%output_path) == 0) error stop 'Output path cannot be empty'
    if (len_trim(job%io_params%output_prefix) == 0) error stop 'Output prefix cannot be empty'
    if (len_trim(job%io_params%odf_path) == 0) error stop 'ODF path cannot be empty'
    if (len_trim(job%io_params%odf_prefix) == 0) error stop 'ODF prefix cannot be empty'
    select case (trim(job%io_params%output_policy))
    case ('ask', 'replace', 'rename_old')
    case default
        error stop 'Output policy must be ask, replace, or rename_old'
    end select
end subroutine validate_c2_te_matrix_workflow

subroutine run_c2_te_matrix_workflow(job)
    type(c2_te_matrix_workflow_type), intent(in) :: job
    type(cv_plan) :: cv
    type(lifting_plan_type) :: lifting_plan
    type(odf_lib_type) :: odf_lib
    logical, allocatable :: pair_active(:,:)
    integer :: nobj
    character(len=1024) :: histogram_path, odf_path

    call validate_c2_te_matrix_workflow(job)
    call read_cv_plan(cv, trim(job%io_params%cv_plan_filename))
    if (.not. allocated(cv%object)) error stop 'CV plan has no object definitions'
    if (.not. allocated(cv%term)) error stop 'CV plan has no term definitions'
    if (size(cv%object) < 1) error stop 'CV plan contains no objects'
    nobj = 2 * size(cv%object)
    call read_pair_active_matrix(trim(job%io_params%pair_active_filename), nobj, pair_active)

    call read_lifting_plan(lifting_plan, trim(job%io_params%lifting_plan_filename))
    call init_odf_lib(cv, odf_lib, lifting_plan, job%model_params%odf_max_lag, &
        job%model_params%odf_basis_index, job%model_params%odf_nlevel, nobj)

    histogram_path = ensure_directory_separator(job%io_params%histogram_path)
    odf_path = ensure_directory_separator(job%io_params%odf_path)
    call run_calculate_te_matrix_chunked_C2(cv_p=cv, odf_lib=odf_lib, lifting_plan=lifting_plan, &
        input_path=trim(histogram_path), input_prefix=trim(job%io_params%histogram_prefix), &
        output_path=trim(job%io_params%output_path), output_prefix=trim(job%io_params%output_prefix), &
        nlevel_target=job%model_params%nlevel_target, &
        max_chunk_size=job%model_params%max_chunk_size, &
        output_policy=trim(job%io_params%output_policy), &
        prefix_odf=trim(job%io_params%odf_prefix), path_odf=trim(odf_path), &
        shuffled=job%model_params%shuffled, pair_active=pair_active)
end subroutine run_c2_te_matrix_workflow

subroutine check_c2_te_matrix_workflow_inputs(job, cv_plan_filename, lifting_plan_filename, &
    cv_plan_will_be_generated, lifting_plan_will_be_generated)
    type(c2_te_matrix_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: cv_plan_filename, lifting_plan_filename
    logical, intent(in), optional :: cv_plan_will_be_generated, lifting_plan_will_be_generated
    type(cv_plan) :: cv
    type(lifting_plan_type) :: lifting_plan
    logical, allocatable :: pair_active(:,:)
    integer :: nobj, istate, ilag
    character(len=1024) :: filename, odf_path, histogram_path, output_file
    character(len=16) :: state_text, lag_text
    logical :: file_exists, cv_generated, plan_generated

    call validate_c2_te_matrix_workflow(job)
    cv_generated = .false.
    if (present(cv_plan_will_be_generated)) cv_generated = cv_plan_will_be_generated
    plan_generated = .false.
    if (present(lifting_plan_will_be_generated)) plan_generated = lifting_plan_will_be_generated

    if (.not. cv_generated) then
        inquire(file=trim(cv_plan_filename), exist=file_exists)
        if (.not. file_exists) call missing_input('CV plan', trim(cv_plan_filename))
        call read_cv_plan(cv, trim(cv_plan_filename))
        if (.not. allocated(cv%object)) error stop 'CV plan has no object definitions'
        if (.not. allocated(cv%term)) error stop 'CV plan has no term definitions'
        if (size(cv%object) < 1) error stop 'CV plan contains no objects'
        nobj = 2 * size(cv%object)
        call read_pair_active_matrix(trim(job%io_params%pair_active_filename), nobj, pair_active)
    else
        inquire(file=trim(job%io_params%pair_active_filename), exist=file_exists)
        if (.not. file_exists) call missing_input('Active-pair matrix', &
            trim(job%io_params%pair_active_filename))
    end if

    if (.not. plan_generated) then
        inquire(file=trim(lifting_plan_filename), exist=file_exists)
        if (.not. file_exists) call missing_input('Lifting plan', trim(lifting_plan_filename))
        call read_lifting_plan(lifting_plan, trim(lifting_plan_filename))
        if (size(lifting_plan%istates) == 0) error stop 'Lifting plan contains no states'
    end if

    histogram_path = ensure_directory_separator(job%io_params%histogram_path)
    inquire(file=trim(histogram_path), exist=file_exists)
    if (.not. file_exists) call missing_input('Histogram directory', trim(histogram_path))
    odf_path = ensure_directory_separator(job%io_params%odf_path)
    if (.not. plan_generated) then
        do istate = 1, size(lifting_plan%istates)
            write(state_text,'(I4.4)') lifting_plan%istates(istate)
            do ilag = 0, job%model_params%odf_max_lag
                write(lag_text,'(I3.3)') ilag
                filename = trim(odf_path)//trim(job%io_params%odf_prefix)//'state_'// &
                    trim(state_text)//'_lag_'//trim(lag_text)//'.txt'
                inquire(file=trim(filename), exist=file_exists)
                if (.not. file_exists) call missing_input('ODF file', trim(filename))
            end do
        end do
    end if
    output_file = trim(job%io_params%output_path)
    if (output_file(len_trim(output_file):len_trim(output_file)) /= '/') output_file = trim(output_file)//'/'
    output_file = trim(output_file)//trim(job%io_params%output_prefix)//'te_matrix.txt'
    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists .and. trim(job%io_params%output_policy) == 'rename_old') then
        filename = trim(output_file)//'.old'
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            write(error_unit,'(A,A)') 'Output backup target already exists: ', trim(filename)
            error stop 'C2 TE matrix output validation failed'
        end if
    end if
end subroutine check_c2_te_matrix_workflow_inputs

subroutine read_pair_active_matrix(filename, nobj, pair_active)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nobj
    logical, allocatable, intent(out) :: pair_active(:,:)
    integer :: unit, ios, i

    allocate(pair_active(nobj,nobj))
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call missing_input('Active-pair matrix', trim(filename))
    do i = 1, nobj
        read(unit,*,iostat=ios) pair_active(i,:)
        if (ios /= 0) then
            close(unit)
            write(error_unit,'(A,I0,A,A)') 'Invalid or incomplete active-pair matrix row ', i, &
                ' in ', trim(filename)
            error stop 'Active-pair matrix input failed'
        end if
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing active-pair matrix file'
end subroutine read_pair_active_matrix

subroutine missing_input(kind, filename)
    character(len=*), intent(in) :: kind, filename
    write(error_unit,'(3A)') trim(kind), ' not found or unavailable: ', trim(filename)
    error stop 'C2 TE matrix workflow input validation failed'
end subroutine missing_input

function ensure_directory_separator(path) result(directory)
    character(len=*), intent(in) :: path
    character(len=1024) :: directory
    integer :: path_length

    directory = trim(path)
    path_length = len_trim(directory)
    if (path_length == 0) error stop 'Directory path cannot be empty'
    if (directory(path_length:path_length) /= '/') directory = trim(directory)//'/'
end function ensure_directory_separator

subroutine write_c2_te_matrix_workflow(job, filename)
    type(c2_te_matrix_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_c2_te_matrix_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call missing_input('Cannot write workflow file', trim(filename))
    write(unit,'(A)') '# C2 transfer-entropy matrix workflow'
    write(unit,'(A)') 'workflow_type = c2_te_matrix'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,A)') 'cv_plan_filename = ', trim(job%io_params%cv_plan_filename)
    write(unit,'(A,A)') 'lifting_plan_filename = ', trim(job%io_params%lifting_plan_filename)
    write(unit,'(A,A)') 'pair_active_filename = ', trim(job%io_params%pair_active_filename)
    write(unit,'(A,A)') 'histogram_path = ', trim(job%io_params%histogram_path)
    write(unit,'(A,A)') 'histogram_prefix = ', trim(job%io_params%histogram_prefix)
    write(unit,'(A,A)') 'output_path = ', trim(job%io_params%output_path)
    write(unit,'(A,A)') 'output_prefix = ', trim(job%io_params%output_prefix)
    write(unit,'(A,A)') 'odf_path = ', trim(job%io_params%odf_path)
    write(unit,'(A,A)') 'odf_prefix = ', trim(job%io_params%odf_prefix)
    write(unit,'(A,I0)') 'nlevel_target = ', job%model_params%nlevel_target
    write(unit,'(A,I0)') 'max_chunk_size = ', job%model_params%max_chunk_size
    write(unit,'(A,I0)') 'odf_max_lag = ', job%model_params%odf_max_lag
    write(unit,'(A,I0)') 'odf_basis_index = ', job%model_params%odf_basis_index
    write(unit,'(A,I0)') 'odf_nlevel = ', job%model_params%odf_nlevel
    write(unit,'(A,L1)') 'shuffled = ', job%model_params%shuffled
    write(unit,'(A,A)') 'output_policy = ', trim(job%io_params%output_policy)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing C2 TE matrix workflow file'
end subroutine write_c2_te_matrix_workflow

subroutine read_c2_te_matrix_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(c2_te_matrix_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1900) :: value
    integer :: unit, ios, equals_position, format_version
    logical :: saw_type, saw_version

    job = c2_te_matrix_workflow_type()
    format_version = -1
    saw_type = .false.
    saw_version = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call missing_input('Cannot read workflow file', trim(filename))
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading C2 TE matrix workflow file'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) error stop 'Invalid C2 TE matrix workflow configuration line'
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'c2_te_matrix') error stop 'Unexpected workflow_type'
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0 .or. format_version /= 1) error stop 'Unsupported C2 TE matrix workflow version'
            saw_version = .true.
        case ('cv_plan_filename')
            job%io_params%cv_plan_filename = trim(value)
        case ('lifting_plan_filename')
            job%io_params%lifting_plan_filename = trim(value)
        case ('pair_active_filename')
            job%io_params%pair_active_filename = trim(value)
        case ('histogram_path')
            job%io_params%histogram_path = trim(value)
        case ('histogram_prefix')
            job%io_params%histogram_prefix = trim(value)
        case ('output_path')
            job%io_params%output_path = trim(value)
        case ('output_prefix')
            job%io_params%output_prefix = trim(value)
        case ('odf_path')
            job%io_params%odf_path = trim(value)
        case ('odf_prefix')
            job%io_params%odf_prefix = trim(value)
        case ('nlevel_target')
            read(value,*,iostat=ios) job%model_params%nlevel_target
            if (ios /= 0) error stop 'Invalid nlevel_target'
        case ('max_chunk_size')
            read(value,*,iostat=ios) job%model_params%max_chunk_size
            if (ios /= 0) error stop 'Invalid max_chunk_size'
        case ('odf_max_lag')
            read(value,*,iostat=ios) job%model_params%odf_max_lag
            if (ios /= 0) error stop 'Invalid odf_max_lag'
        case ('odf_basis_index')
            read(value,*,iostat=ios) job%model_params%odf_basis_index
            if (ios /= 0) error stop 'Invalid odf_basis_index'
        case ('odf_nlevel')
            read(value,*,iostat=ios) job%model_params%odf_nlevel
            if (ios /= 0) error stop 'Invalid odf_nlevel'
        case ('shuffled')
            read(value,*,iostat=ios) job%model_params%shuffled
            if (ios /= 0) error stop 'Invalid shuffled flag'
        case ('output_policy')
            job%io_params%output_policy = trim(value)
        case default
            write(error_unit,'(A,A)') 'Unknown C2 TE matrix workflow key: ', trim(key)
            error stop 'Invalid C2 TE matrix workflow configuration'
        end select
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing C2 TE matrix workflow file'
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete C2 TE matrix workflow header'
    call validate_c2_te_matrix_workflow(job)
end subroutine read_c2_te_matrix_workflow

subroutine print_c2_te_matrix_workflow(job)
    type(c2_te_matrix_workflow_type), intent(in) :: job
    character(len=32) :: shuffle_note
    shuffle_note = ''
    if (job%model_params%shuffled) shuffle_note = ' (shuffle run 1)'
    write(output_unit,'(A,A)') 'CV plan: ', trim(job%io_params%cv_plan_filename)
    write(output_unit,'(A,A)') 'Lifting plan: ', trim(job%io_params%lifting_plan_filename)
    write(output_unit,'(A,A)') 'Pair matrix: ', trim(job%io_params%pair_active_filename)
    write(output_unit,'(A,A,A,A)') 'Input histograms: ', trim(job%io_params%histogram_path), &
        trim(job%io_params%histogram_prefix), trim(shuffle_note)
    write(output_unit,'(A,A,A)') 'Output matrix: ', trim(job%io_params%output_path), &
        trim(job%io_params%output_prefix)//'te_matrix.txt'
    write(output_unit,'(A,A,A)') 'ODF files: ', trim(job%io_params%odf_path), trim(job%io_params%odf_prefix)
    write(output_unit,'(A,I0)') 'Histogram HEALPix levels: ', job%model_params%nlevel_target
    write(output_unit,'(A,I0)') 'Maximum chunk size: ', job%model_params%max_chunk_size
    write(output_unit,'(A,I0)') 'ODF max lag: ', job%model_params%odf_max_lag
    write(output_unit,'(A,I0)') 'ODF basis index: ', job%model_params%odf_basis_index
    write(output_unit,'(A,I0)') 'ODF levels: ', job%model_params%odf_nlevel
    write(output_unit,'(A,A)') 'Output policy: ', trim(job%io_params%output_policy)
end subroutine print_c2_te_matrix_workflow

end module c2_te_matrix_workflow_module
