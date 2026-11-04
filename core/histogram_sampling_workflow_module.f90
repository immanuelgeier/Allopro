module histogram_sampling_workflow_module
    use iso_fortran_env, only: output_unit, error_unit, int32
    use workflow_module, only: histogram_sampling_workflow_type, workflow_config_filename, &
        histogram_sampling_workflow_config_stem
    use direction_vectors_analyse_plan_module, only: cv_plan
    use cv_plan_io_module, only: read_cv_plan
    use lifting_plan_module, only: lifting_plan_type, read_lifting_plan
    use healpix_state_space_exploration, only: odf_lib_type, init_odf_lib
    use healpix_histogram_sampling_execution_module, only: run_te_workflow
    implicit none
    private

    public :: histogram_sampling_workflow_type
    public :: configure_histogram_sampling_workflow
    public :: read_histogram_sampling_workflow
    public :: write_histogram_sampling_workflow
    public :: validate_histogram_sampling_workflow
    public :: run_histogram_sampling_workflow
    public :: check_histogram_sampling_inputs

contains

subroutine configure_histogram_sampling_workflow(job)
    type(histogram_sampling_workflow_type), intent(out) :: job
    character(len=1024) :: input, config_filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' Histogram Sampling Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid histogram sampling workflow option'

    select case (choice)
    case (1)
        job = histogram_sampling_workflow_type()
        call edit_histogram_sampling_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_histogram_sampling_workflow(trim(config_filename), job)
        if (choice == 2) then
            call edit_histogram_sampling_workflow(job)
        else
            call print_histogram_sampling_workflow(job)
        end if
    case (4)
        job = histogram_sampling_workflow_type()
        write(output_unit,*) 'Histogram sampling workflow configuration skipped.'
        return
    case default
        error stop 'Invalid histogram sampling workflow option'
    end select

    call validate_histogram_sampling_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        config_filename = workflow_config_filename(trim(input), histogram_sampling_workflow_config_stem)
        call write_histogram_sampling_workflow(job, trim(config_filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(config_filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_histogram_sampling_workflow

subroutine edit_histogram_sampling_workflow(job)
    type(histogram_sampling_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_integer('State ID', job%model_params%state)
    call edit_string('CV plan filename', job%io_params%cv_plan_filename)
    call edit_string('Lifting plan filename', job%io_params%lifting_plan_filename)
    call edit_string('Active-pair matrix filename', job%io_params%pair_active_filename)
    call edit_string('ODF directory', job%io_params%odf_path)
    call edit_string('ODF filename prefix', job%io_params%odf_prefix)
    call edit_integer('ODF maximum lag (number of lag files minus one)', job%model_params%odf_max_lag)
    call edit_integer('ODF basis index', job%model_params%odf_basis_index)
    call edit_integer('ODF HEALPix levels', job%model_params%odf_nlevel)
    call edit_string('Lifted HEALPix trajectory directory', job%io_params%lifted_trajectory_path)
    call edit_string('Lifted HEALPix trajectory prefix', job%io_params%lifted_trajectory_prefix)
    call edit_string('Histogram output directory', job%io_params%histogram_path)
    call edit_string('Histogram output prefix', job%io_params%histogram_prefix)
    call edit_integer('Sampling lag', job%model_params%lag)
    call edit_integer('Sampling HEALPix levels', job%model_params%nlevel_target)
    call edit_integer('Maximum chunk size', job%model_params%max_chunk_size)
    call edit_real('Memory time constant tau_mem', job%model_params%tau_mem)
    call edit_real('Time step dt', job%model_params%dt)
    call edit_real('Maximum sampling-plan RAM in GB', job%model_params%max_ram_gb)
    call edit_logical('Block shuffle enabled', job%model_params%block_shuffle)
    call edit_integer('Number of shuffled runs', job%model_params%nrun_shuffle)
    call edit_integer('Shuffle block size', job%model_params%block_size_shuffle)
    call edit_integer('Backbone basis index for lifting', job%model_params%basis_index)
    call edit_string('Histogram output policy (ask/replace/new)', &
        job%io_params%histogram_output_policy)
    call validate_histogram_sampling_workflow(job)
    call print_histogram_sampling_workflow(job)

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

    subroutine edit_real(prompt, value)
        character(len=*), intent(in) :: prompt
        real(8), intent(inout) :: value
        write(output_unit,'(A,ES14.6,A)', advance='no') trim(prompt), value, ': '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid real workflow parameter'
    end subroutine edit_real

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
end subroutine edit_histogram_sampling_workflow

subroutine validate_histogram_sampling_workflow(job)
    type(histogram_sampling_workflow_type), intent(in) :: job
    integer :: nside

    if (job%model_params%state < 1) error stop 'Histogram sampling state ID must be positive'
    if (job%model_params%lag < 0) error stop 'Histogram sampling lag cannot be negative'
    if (job%model_params%basis_index < 1) error stop 'Lifting basis index must be positive'
    if (job%model_params%nlevel_target < 1 .or. job%model_params%nlevel_target > 14) then
        error stop 'Sampling HEALPix levels must be between 1 and 14'
    end if
    if (job%model_params%odf_max_lag < 0) error stop 'ODF maximum lag cannot be negative'
    if (job%model_params%odf_basis_index < 1) error stop 'ODF basis index must be positive'
    if (job%model_params%odf_nlevel < 1 .or. job%model_params%odf_nlevel > 14) then
        error stop 'ODF HEALPix levels must be between 1 and 14'
    end if
    if (job%model_params%max_chunk_size < 1) error stop 'Maximum chunk size must be positive'
    if (job%model_params%tau_mem <= 0.0_8) error stop 'tau_mem must be positive'
    if (job%model_params%dt <= 0.0_8) error stop 'dt must be positive'
    if (job%model_params%max_ram_gb <= 0.0_8) error stop 'Maximum RAM must be positive'
    if (job%model_params%nrun_shuffle < 1) error stop 'Number of shuffled runs must be positive'
    if (job%model_params%block_size_shuffle < 1) error stop 'Shuffle block size must be positive'
    if (len_trim(job%io_params%cv_plan_filename) == 0) error stop 'CV plan filename cannot be empty'
    if (len_trim(job%io_params%lifting_plan_filename) == 0) error stop 'Lifting plan filename cannot be empty'
    if (len_trim(job%io_params%pair_active_filename) == 0) error stop 'Pair matrix filename cannot be empty'
    if (len_trim(job%io_params%odf_path) == 0) error stop 'ODF path cannot be empty'
    if (len_trim(job%io_params%odf_prefix) == 0) error stop 'ODF prefix cannot be empty'
    if (len_trim(job%io_params%lifted_trajectory_path) == 0) error stop 'Lifted trajectory path cannot be empty'
    if (len_trim(job%io_params%lifted_trajectory_prefix) == 0) then
        error stop 'Lifted trajectory prefix cannot be empty'
    end if
    if (len_trim(job%io_params%histogram_path) == 0) error stop 'Histogram path cannot be empty'
    if (len_trim(job%io_params%histogram_prefix) == 0) error stop 'Histogram prefix cannot be empty'
    select case (trim(job%io_params%histogram_output_policy))
    case ('ask', 'replace', 'new')
    case default
        error stop 'Histogram output policy must be ask, replace, or new'
    end select

    nside = 2**(job%model_params%nlevel_target - 1)
    if (12 * nside * nside <= 0) error stop 'Invalid sampling HEALPix resolution'
end subroutine validate_histogram_sampling_workflow

subroutine run_histogram_sampling_workflow(job)
    type(histogram_sampling_workflow_type), intent(in) :: job
    type(cv_plan) :: cv
    type(lifting_plan_type) :: full_plan, selected_plan
    type(odf_lib_type) :: odf_lib
    logical, allocatable :: pair_active(:,:)
    integer :: nobj, state_index
    character(len=1024) :: odf_path, lifted_path, histogram_path
    character(len=1) :: unused_path, unused_prefix

    call validate_histogram_sampling_workflow(job)
    call read_cv_plan(cv, trim(job%io_params%cv_plan_filename))
    if (.not. allocated(cv%object)) error stop 'CV plan has no object definitions'
    if (size(cv%object) < 1) error stop 'CV plan contains no objects'
    nobj = 2 * size(cv%object)
    call read_pair_active_matrix(trim(job%io_params%pair_active_filename), nobj, pair_active)

    call read_lifting_plan(full_plan, trim(job%io_params%lifting_plan_filename))
    state_index = find_lifting_state(full_plan, job%model_params%state)
    call select_lifting_plan_state(full_plan, state_index, selected_plan)

    call init_odf_lib(cv, odf_lib, selected_plan, job%model_params%odf_max_lag, &
        job%model_params%odf_basis_index, job%model_params%odf_nlevel, nobj)

    odf_path = ensure_directory_separator(job%io_params%odf_path)
    lifted_path = ensure_directory_separator(job%io_params%lifted_trajectory_path)
    histogram_path = ensure_directory_separator(job%io_params%histogram_path)
    unused_path = ' '
    unused_prefix = ' '
    call run_te_workflow(odf_lib=odf_lib, lifting_plan=selected_plan, lag=job%model_params%lag, &
        basis_index=job%model_params%basis_index, prefix_vec_tra=unused_prefix, path_vec_tra=unused_path, &
        prefix_tra=trim(job%io_params%lifted_trajectory_prefix), path_tra=trim(lifted_path), &
        file_policy_tra='ask', prefix_hist=trim(job%io_params%histogram_prefix), &
        path_hist=trim(histogram_path), file_policy_hist=trim(job%io_params%histogram_output_policy), &
        nlevel_target=job%model_params%nlevel_target, &
        max_chunk_size=int(job%model_params%max_chunk_size,int32), tau_mem=job%model_params%tau_mem, &
        dt=job%model_params%dt, max_ram_gb=job%model_params%max_ram_gb, &
        block_shuffle=job%model_params%block_shuffle, nrun_shuffle=job%model_params%nrun_shuffle, &
        block_size_shuffle=job%model_params%block_size_shuffle, only_lifting=.false., &
        prefix_odf=trim(job%io_params%odf_prefix), path_odf=trim(odf_path), pair_active=pair_active)
end subroutine run_histogram_sampling_workflow

subroutine check_histogram_sampling_inputs(job, cv_plan_filename, lifting_plan_filename, &
    cv_plan_will_be_generated, lifting_plan_will_be_generated)
    type(histogram_sampling_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: cv_plan_filename, lifting_plan_filename
    logical, intent(in), optional :: cv_plan_will_be_generated, lifting_plan_will_be_generated
    type(cv_plan) :: cv
    type(lifting_plan_type) :: lifting_plan
    logical, allocatable :: pair_active(:,:)
    integer :: nobj, state_index, state_id, ilag
    character(len=1024) :: filename, odf_path
    character(len=16) :: state_text, lag_text
    logical :: file_exists, cv_generated, lifting_plan_generated

    call validate_histogram_sampling_workflow(job)
    cv_generated = .false.
    if (present(cv_plan_will_be_generated)) cv_generated = cv_plan_will_be_generated
    lifting_plan_generated = .false.
    if (present(lifting_plan_will_be_generated)) then
        lifting_plan_generated = lifting_plan_will_be_generated
    end if
    if (.not. cv_generated) then
        inquire(file=trim(cv_plan_filename), exist=file_exists)
        if (.not. file_exists) then
            write(error_unit,'(A,A)') 'CV plan file not found: ', trim(cv_plan_filename)
            error stop 'Histogram sampling dry-run failed'
        end if
    end if
    if (.not. lifting_plan_generated) then
        inquire(file=trim(lifting_plan_filename), exist=file_exists)
        if (.not. file_exists) then
            write(error_unit,'(A,A)') 'Lifting plan file not found: ', trim(lifting_plan_filename)
            error stop 'Histogram sampling dry-run failed'
        end if
    end if
    inquire(file=trim(job%io_params%pair_active_filename), exist=file_exists)
    if (.not. file_exists) then
        write(error_unit,'(A,A)') 'Active-pair matrix file not found: ', &
            trim(job%io_params%pair_active_filename)
        error stop 'Histogram sampling dry-run failed'
    end if

    if (.not. cv_generated) then
        call read_cv_plan(cv, trim(cv_plan_filename))
        if (.not. allocated(cv%object)) error stop 'CV plan has no object definitions'
        nobj = 2 * size(cv%object)
        call read_pair_active_matrix(trim(job%io_params%pair_active_filename), nobj, pair_active)
    end if
    if (.not. lifting_plan_generated) then
        call read_lifting_plan(lifting_plan, trim(lifting_plan_filename))
        state_index = find_lifting_state(lifting_plan, job%model_params%state)
        state_id = lifting_plan%istates(state_index)
    else
        state_id = job%model_params%state
    end if
    odf_path = ensure_directory_separator(job%io_params%odf_path)
    write(state_text,'(I4.4)') state_id
    do ilag = 0, job%model_params%odf_max_lag
        write(lag_text,'(I3.3)') ilag
        filename = trim(odf_path)//trim(job%io_params%odf_prefix)//'state_'// &
            trim(state_text)//'_lag_'//trim(lag_text)//'.txt'
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            write(error_unit,'(A,A)') 'ODF input file not found: ', trim(filename)
            error stop 'Histogram sampling dry-run failed'
        end if
    end do
end subroutine check_histogram_sampling_inputs

subroutine read_pair_active_matrix(filename, nobj, pair_active)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nobj
    logical, allocatable, intent(out) :: pair_active(:,:)
    integer :: unit, ios, i

    allocate(pair_active(nobj,nobj))
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read active-pair matrix: ', trim(filename)
        error stop 'Active-pair matrix input failed'
    end if
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

integer function find_lifting_state(lifting_plan, state_id) result(state_index)
    type(lifting_plan_type), intent(in) :: lifting_plan
    integer, intent(in) :: state_id
    integer :: i

    do i = 1, size(lifting_plan%istates)
        if (lifting_plan%istates(i) == state_id) then
            state_index = i
            return
        end if
    end do
    write(error_unit,'(A,I0)') 'State not present in lifting plan: ', state_id
    error stop 'Histogram sampling state is unavailable'
end function find_lifting_state

subroutine select_lifting_plan_state(source, state_index, selected)
    type(lifting_plan_type), intent(in) :: source
    integer, intent(in) :: state_index
    type(lifting_plan_type), intent(out) :: selected
    integer :: nreplicas, nblocks, block_size

    nreplicas = size(source%ireps, 2)
    nblocks = size(source%start_block, 3)
    block_size = size(source%timestep_pattern, 1)
    allocate(selected%istates(1), selected%nreps(1))
    allocate(selected%ireps(1,nreplicas), selected%nblock(1,nreplicas))
    allocate(selected%start_block(1,nreplicas,nblocks), selected%end_block(1,nreplicas,nblocks))
    allocate(selected%timestep_pattern(block_size,1,nreplicas,nblocks))
    selected%istates(1) = source%istates(state_index)
    selected%nreps(1) = source%nreps(state_index)
    selected%ireps(1,:) = source%ireps(state_index,:)
    selected%nblock(1,:) = source%nblock(state_index,:)
    selected%start_block(1,:,:) = source%start_block(state_index,:,:)
    selected%end_block(1,:,:) = source%end_block(state_index,:,:)
    selected%timestep_pattern(:,1,:,:) = source%timestep_pattern(:,state_index,:,:)
end subroutine select_lifting_plan_state

function ensure_directory_separator(path) result(directory)
    character(len=*), intent(in) :: path
    character(len=1024) :: directory
    integer :: path_length

    path_length = len_trim(path)
    directory = trim(path)
    if (directory(path_length:path_length) /= '/') directory = trim(directory)//'/'
end function ensure_directory_separator

subroutine write_histogram_sampling_workflow(job, filename)
    type(histogram_sampling_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_histogram_sampling_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write histogram sampling workflow: ', trim(filename)
        error stop 'Histogram sampling workflow output failed'
    end if
    write(unit,'(A)') '# Histogram sampling workflow'
    write(unit,'(A)') 'workflow_type = histogram_sampling'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,I0)') 'state = ', job%model_params%state
    write(unit,'(A,A)') 'cv_plan_filename = ', trim(job%io_params%cv_plan_filename)
    write(unit,'(A,A)') 'lifting_plan_filename = ', trim(job%io_params%lifting_plan_filename)
    write(unit,'(A,A)') 'pair_active_filename = ', trim(job%io_params%pair_active_filename)
    write(unit,'(A,A)') 'odf_path = ', trim(job%io_params%odf_path)
    write(unit,'(A,A)') 'odf_prefix = ', trim(job%io_params%odf_prefix)
    write(unit,'(A,I0)') 'odf_max_lag = ', job%model_params%odf_max_lag
    write(unit,'(A,I0)') 'odf_basis_index = ', job%model_params%odf_basis_index
    write(unit,'(A,I0)') 'odf_nlevel = ', job%model_params%odf_nlevel
    write(unit,'(A,A)') 'lifted_trajectory_path = ', trim(job%io_params%lifted_trajectory_path)
    write(unit,'(A,A)') 'lifted_trajectory_prefix = ', trim(job%io_params%lifted_trajectory_prefix)
    write(unit,'(A,A)') 'histogram_path = ', trim(job%io_params%histogram_path)
    write(unit,'(A,A)') 'histogram_prefix = ', trim(job%io_params%histogram_prefix)
    write(unit,'(A,I0)') 'lag = ', job%model_params%lag
    write(unit,'(A,I0)') 'basis_index = ', job%model_params%basis_index
    write(unit,'(A,I0)') 'nlevel_target = ', job%model_params%nlevel_target
    write(unit,'(A,I0)') 'max_chunk_size = ', job%model_params%max_chunk_size
    write(unit,'(A,ES24.16E3)') 'tau_mem = ', job%model_params%tau_mem
    write(unit,'(A,ES24.16E3)') 'dt = ', job%model_params%dt
    write(unit,'(A,ES24.16E3)') 'max_ram_gb = ', job%model_params%max_ram_gb
    write(unit,'(A,L1)') 'block_shuffle = ', job%model_params%block_shuffle
    write(unit,'(A,I0)') 'nrun_shuffle = ', job%model_params%nrun_shuffle
    write(unit,'(A,I0)') 'block_size_shuffle = ', job%model_params%block_size_shuffle
    write(unit,'(A,A)') 'histogram_output_policy = ', trim(job%io_params%histogram_output_policy)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing histogram sampling workflow file'
end subroutine write_histogram_sampling_workflow

subroutine read_histogram_sampling_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(histogram_sampling_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1900) :: value
    integer :: unit, ios, equals_position, format_version
    logical :: saw_type, saw_version

    job = histogram_sampling_workflow_type()
    format_version = -1
    saw_type = .false.
    saw_version = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read histogram sampling workflow: ', trim(filename)
        error stop 'Histogram sampling workflow input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) call config_read_error(unit, 'Error reading histogram sampling workflow')
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) call config_read_error(unit, 'Invalid workflow configuration line')
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'histogram_sampling') call config_read_error(unit, 'Unexpected workflow_type')
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0) call config_read_error(unit, 'Invalid format_version')
            if (format_version /= 1) call config_read_error(unit, 'Unsupported format_version')
            saw_version = .true.
        case ('state')
            read(value,*,iostat=ios) job%model_params%state
            if (ios /= 0) call config_read_error(unit, 'Invalid state')
        case ('cv_plan_filename')
            job%io_params%cv_plan_filename = trim(value)
        case ('lifting_plan_filename')
            job%io_params%lifting_plan_filename = trim(value)
        case ('pair_active_filename')
            job%io_params%pair_active_filename = trim(value)
        case ('odf_path')
            job%io_params%odf_path = trim(value)
        case ('odf_prefix')
            job%io_params%odf_prefix = trim(value)
        case ('odf_max_lag')
            read(value,*,iostat=ios) job%model_params%odf_max_lag
            if (ios /= 0) call config_read_error(unit, 'Invalid odf_max_lag')
        case ('odf_basis_index')
            read(value,*,iostat=ios) job%model_params%odf_basis_index
            if (ios /= 0) call config_read_error(unit, 'Invalid odf_basis_index')
        case ('odf_nlevel')
            read(value,*,iostat=ios) job%model_params%odf_nlevel
            if (ios /= 0) call config_read_error(unit, 'Invalid odf_nlevel')
        case ('lifted_trajectory_path')
            job%io_params%lifted_trajectory_path = trim(value)
        case ('lifted_trajectory_prefix')
            job%io_params%lifted_trajectory_prefix = trim(value)
        case ('histogram_path')
            job%io_params%histogram_path = trim(value)
        case ('histogram_prefix')
            job%io_params%histogram_prefix = trim(value)
        case ('lag')
            read(value,*,iostat=ios) job%model_params%lag
            if (ios /= 0) call config_read_error(unit, 'Invalid lag')
        case ('basis_index')
            read(value,*,iostat=ios) job%model_params%basis_index
            if (ios /= 0) call config_read_error(unit, 'Invalid basis_index')
        case ('nlevel_target')
            read(value,*,iostat=ios) job%model_params%nlevel_target
            if (ios /= 0) call config_read_error(unit, 'Invalid nlevel_target')
        case ('max_chunk_size')
            read(value,*,iostat=ios) job%model_params%max_chunk_size
            if (ios /= 0) call config_read_error(unit, 'Invalid max_chunk_size')
        case ('tau_mem')
            read(value,*,iostat=ios) job%model_params%tau_mem
            if (ios /= 0) call config_read_error(unit, 'Invalid tau_mem')
        case ('dt')
            read(value,*,iostat=ios) job%model_params%dt
            if (ios /= 0) call config_read_error(unit, 'Invalid dt')
        case ('max_ram_gb')
            read(value,*,iostat=ios) job%model_params%max_ram_gb
            if (ios /= 0) call config_read_error(unit, 'Invalid max_ram_gb')
        case ('block_shuffle')
            read(value,*,iostat=ios) job%model_params%block_shuffle
            if (ios /= 0) call config_read_error(unit, 'Invalid block_shuffle')
        case ('nrun_shuffle')
            read(value,*,iostat=ios) job%model_params%nrun_shuffle
            if (ios /= 0) call config_read_error(unit, 'Invalid nrun_shuffle')
        case ('block_size_shuffle')
            read(value,*,iostat=ios) job%model_params%block_size_shuffle
            if (ios /= 0) call config_read_error(unit, 'Invalid block_size_shuffle')
        case ('histogram_output_policy')
            job%io_params%histogram_output_policy = trim(value)
        case default
            call config_read_error(unit, 'Unknown workflow parameter: '//trim(key))
        end select
    end do
    close(unit)
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete histogram sampling workflow file'
    call validate_histogram_sampling_workflow(job)
end subroutine read_histogram_sampling_workflow

subroutine print_histogram_sampling_workflow(job)
    type(histogram_sampling_workflow_type), intent(in) :: job
    write(output_unit,'(A,I0)') 'State ID: ', job%model_params%state
    write(output_unit,'(A,A)') 'CV plan: ', trim(job%io_params%cv_plan_filename)
    write(output_unit,'(A,A)') 'Lifting plan: ', trim(job%io_params%lifting_plan_filename)
    write(output_unit,'(A,A)') 'Active-pair matrix: ', trim(job%io_params%pair_active_filename)
    write(output_unit,'(A,A)') 'ODF path: ', trim(job%io_params%odf_path)
    write(output_unit,'(A,A)') 'ODF prefix: ', trim(job%io_params%odf_prefix)
    write(output_unit,'(A,I0)') 'ODF max lag: ', job%model_params%odf_max_lag
    write(output_unit,'(A,I0)') 'ODF basis index: ', job%model_params%odf_basis_index
    write(output_unit,'(A,I0)') 'ODF levels: ', job%model_params%odf_nlevel
    write(output_unit,'(A,A)') 'Lifted trajectory path: ', trim(job%io_params%lifted_trajectory_path)
    write(output_unit,'(A,A)') 'Lifted trajectory prefix: ', trim(job%io_params%lifted_trajectory_prefix)
    write(output_unit,'(A,A)') 'Histogram path: ', trim(job%io_params%histogram_path)
    write(output_unit,'(A,A)') 'Histogram prefix: ', trim(job%io_params%histogram_prefix)
    write(output_unit,'(A,I0)') 'Sampling lag: ', job%model_params%lag
    write(output_unit,'(A,I0)') 'Sampling levels: ', job%model_params%nlevel_target
    write(output_unit,'(A,I0)') 'Maximum chunk size: ', job%model_params%max_chunk_size
    write(output_unit,'(A,ES14.6)') 'tau_mem: ', job%model_params%tau_mem
    write(output_unit,'(A,ES14.6)') 'dt: ', job%model_params%dt
    write(output_unit,'(A,ES14.6)') 'Maximum RAM (GB): ', job%model_params%max_ram_gb
    write(output_unit,'(A,L1)') 'Block shuffle: ', job%model_params%block_shuffle
    write(output_unit,'(A,I0)') 'Shuffle runs: ', job%model_params%nrun_shuffle
    write(output_unit,'(A,I0)') 'Shuffle block size: ', job%model_params%block_size_shuffle
    write(output_unit,'(A,A)') 'Histogram output policy: ', trim(job%io_params%histogram_output_policy)
end subroutine print_histogram_sampling_workflow

subroutine config_read_error(unit, message)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    close(unit)
    write(error_unit,'(A)') trim(message)
    error stop 'Histogram sampling workflow configuration read failed'
end subroutine config_read_error

end module histogram_sampling_workflow_module
