module odf_sampling_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use workflow_module, only: odf_sampling_workflow_type, workflow_config_filename, &
        odf_sampling_workflow_config_stem
    use direction_vectors_analyse_plan_module, only: cv_plan
    use cv_plan_io_module, only: read_cv_plan
    use lifting_plan_module, only: lifting_plan_type, read_lifting_plan
    use healpix_state_space_exploration, only: odf_lib_type, init_odf_lib, sample_odfs
    implicit none
    private

    public :: odf_sampling_workflow_type
    public :: configure_odf_sampling_workflow
    public :: read_odf_sampling_workflow
    public :: write_odf_sampling_workflow
    public :: validate_odf_sampling_workflow
    public :: run_odf_sampling_workflow

contains

subroutine configure_odf_sampling_workflow(job)
    type(odf_sampling_workflow_type), intent(out) :: job
    character(len=1024) :: input, config_filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' ODF Sampling Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid ODF sampling workflow option'

    select case (choice)
    case (1)
        job = odf_sampling_workflow_type()
        call edit_odf_sampling_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_odf_sampling_workflow(trim(config_filename), job)
        if (choice == 2) then
            call edit_odf_sampling_workflow(job)
        else
            call print_odf_sampling_workflow(job)
        end if
    case (4)
        job = odf_sampling_workflow_type()
        write(output_unit,*) 'ODF sampling workflow configuration skipped.'
        return
    case default
        error stop 'Invalid ODF sampling workflow option'
    end select

    call validate_odf_sampling_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        config_filename = workflow_config_filename(trim(input), odf_sampling_workflow_config_stem)
        call write_odf_sampling_workflow(job, trim(config_filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(config_filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_odf_sampling_workflow

subroutine edit_odf_sampling_workflow(job)
    type(odf_sampling_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_integer('State ID', job%model_params%state)
    call edit_string('CV plan filename', job%io_params%cv_plan_filename)
    call edit_string('Lifting plan filename', job%io_params%lifting_plan_filename)
    call edit_string('Lifted trajectory directory', job%io_params%trajectory_path)
    call edit_string('Lifted trajectory filename prefix', job%io_params%trajectory_prefix)
    call edit_integer('Maximum lag', job%model_params%max_lag)
    call edit_integer('Backbone basis index', job%model_params%basis_index)
    call edit_integer('HEALPix level count', job%model_params%nlevel)

    write(output_unit,'(A,ES14.6,A)', advance='no') 'Memory time constant tau_mem [', &
        job%model_params%tau_mem, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%model_params%tau_mem
        if (ios /= 0) error stop 'Invalid tau_mem value'
    end if

    write(output_unit,'(A,ES14.6,A)', advance='no') 'Time step dt [', job%model_params%dt, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%model_params%dt
        if (ios /= 0) error stop 'Invalid dt value'
    end if

    call edit_string('ODF output directory', job%io_params%odf_path)
    call edit_string('ODF filename prefix', job%io_params%odf_prefix)
    call edit_string('Output policy (ask/replace/new)', job%io_params%output_policy)
    call validate_odf_sampling_workflow(job)
    call print_odf_sampling_workflow(job)

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

    subroutine edit_string(prompt, value)
        character(len=*), intent(in) :: prompt
        character(len=*), intent(inout) :: value
        write(output_unit,'(A,A,A)', advance='no') trim(prompt), ' [', trim(value), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) value = trim(input)
    end subroutine edit_string
end subroutine edit_odf_sampling_workflow

subroutine validate_odf_sampling_workflow(job)
    type(odf_sampling_workflow_type), intent(in) :: job

    if (job%model_params%state < 1) error stop 'ODF sampling state ID must be positive'
    if (job%model_params%max_lag < 0) error stop 'ODF sampling maximum lag cannot be negative'
    if (job%model_params%basis_index < 1) error stop 'ODF sampling basis index must be positive'
    if (job%model_params%nlevel < 1 .or. job%model_params%nlevel > 14) then
        error stop 'ODF sampling HEALPix level count must be between 1 and 14'
    end if
    if (job%model_params%tau_mem <= 0.0_8) error stop 'ODF sampling tau_mem must be positive'
    if (job%model_params%dt <= 0.0_8) error stop 'ODF sampling dt must be positive'
    if (len_trim(job%io_params%cv_plan_filename) == 0) error stop 'CV plan filename cannot be empty'
    if (len_trim(job%io_params%lifting_plan_filename) == 0) error stop 'Lifting plan filename cannot be empty'
    if (len_trim(job%io_params%trajectory_path) == 0) error stop 'Trajectory path cannot be empty'
    if (len_trim(job%io_params%trajectory_prefix) == 0) error stop 'Trajectory prefix cannot be empty'
    if (len_trim(job%io_params%odf_path) == 0) error stop 'ODF output path cannot be empty'
    if (len_trim(job%io_params%odf_prefix) == 0) error stop 'ODF output prefix cannot be empty'
    select case (trim(job%io_params%output_policy))
    case ('ask', 'replace', 'new')
    case default
        error stop 'ODF sampling output policy must be ask, replace, or new'
    end select
end subroutine validate_odf_sampling_workflow

subroutine run_odf_sampling_workflow(job)
    type(odf_sampling_workflow_type), intent(in) :: job
    type(cv_plan) :: cv
    type(lifting_plan_type) :: full_plan, selected_plan
    type(odf_lib_type) :: odf_lib
    integer :: state_index, nobj
    character(len=1024) :: trajectory_path, odf_path

    call validate_odf_sampling_workflow(job)
    call read_cv_plan(cv, trim(job%io_params%cv_plan_filename))
    if (.not. allocated(cv%object)) error stop 'CV plan has no object definitions'
    if (size(cv%object) < 1) error stop 'CV plan must define objects for ODF sampling'
    nobj = 2 * size(cv%object)

    call read_lifting_plan(full_plan, trim(job%io_params%lifting_plan_filename))
    state_index = find_state(full_plan, job%model_params%state)
    call select_lifting_plan_state(full_plan, state_index, selected_plan)

    call init_odf_lib(cv, odf_lib, selected_plan, job%model_params%max_lag, &
        job%model_params%basis_index, job%model_params%nlevel, nobj)

    trajectory_path = trim(job%io_params%trajectory_path)
    odf_path = trim(job%io_params%odf_path)
    if (odf_path(len_trim(odf_path):len_trim(odf_path)) /= '/') then
        odf_path = trim(odf_path)//'/'
    end if

    call sample_odfs(odf_lib, selected_plan, trim(job%io_params%trajectory_prefix), &
        trim(trajectory_path), job%model_params%basis_index, job%model_params%tau_mem, &
        job%model_params%dt, trim(job%io_params%odf_prefix), trim(odf_path), &
        trim(job%io_params%output_policy), perform_analysis=.false.)
end subroutine run_odf_sampling_workflow

integer function find_state(lifting_plan, state_id) result(state_index)
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
    error stop 'ODF sampling state is unavailable'
end function find_state

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

subroutine write_odf_sampling_workflow(job, filename)
    type(odf_sampling_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_odf_sampling_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write ODF sampling workflow: ', trim(filename)
        error stop 'ODF sampling workflow output failed'
    end if
    write(unit,'(A)') '# ODF sampling workflow'
    write(unit,'(A)') 'workflow_type = odf_sampling'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,I0)') 'state = ', job%model_params%state
    write(unit,'(A,A)') 'cv_plan_filename = ', trim(job%io_params%cv_plan_filename)
    write(unit,'(A,A)') 'lifting_plan_filename = ', trim(job%io_params%lifting_plan_filename)
    write(unit,'(A,A)') 'trajectory_path = ', trim(job%io_params%trajectory_path)
    write(unit,'(A,A)') 'trajectory_prefix = ', trim(job%io_params%trajectory_prefix)
    write(unit,'(A,I0)') 'max_lag = ', job%model_params%max_lag
    write(unit,'(A,I0)') 'basis_index = ', job%model_params%basis_index
    write(unit,'(A,I0)') 'nlevel = ', job%model_params%nlevel
    write(unit,'(A,ES24.16E3)') 'tau_mem = ', job%model_params%tau_mem
    write(unit,'(A,ES24.16E3)') 'dt = ', job%model_params%dt
    write(unit,'(A,A)') 'odf_path = ', trim(job%io_params%odf_path)
    write(unit,'(A,A)') 'odf_prefix = ', trim(job%io_params%odf_prefix)
    write(unit,'(A,A)') 'output_policy = ', trim(job%io_params%output_policy)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing ODF sampling workflow file'
end subroutine write_odf_sampling_workflow

subroutine read_odf_sampling_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(odf_sampling_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1900) :: value
    integer :: unit, ios, equals_position, format_version
    logical :: saw_type, saw_version

    job = odf_sampling_workflow_type()
    format_version = -1
    saw_type = .false.
    saw_version = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read ODF sampling workflow: ', trim(filename)
        error stop 'ODF sampling workflow input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) call odf_config_read_error(unit, 'Error reading ODF sampling workflow')
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) call odf_config_read_error(unit, 'Invalid workflow configuration line')
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'odf_sampling') call odf_config_read_error(unit, 'Unexpected workflow_type')
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid format_version')
            if (format_version /= 1) call odf_config_read_error(unit, 'Unsupported format_version')
            saw_version = .true.
        case ('state')
            read(value,*,iostat=ios) job%model_params%state
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid state')
        case ('cv_plan_filename')
            job%io_params%cv_plan_filename = trim(value)
        case ('lifting_plan_filename')
            job%io_params%lifting_plan_filename = trim(value)
        case ('trajectory_path')
            job%io_params%trajectory_path = trim(value)
        case ('trajectory_prefix')
            job%io_params%trajectory_prefix = trim(value)
        case ('max_lag')
            read(value,*,iostat=ios) job%model_params%max_lag
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid max_lag')
        case ('basis_index')
            read(value,*,iostat=ios) job%model_params%basis_index
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid basis_index')
        case ('nlevel')
            read(value,*,iostat=ios) job%model_params%nlevel
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid nlevel')
        case ('tau_mem')
            read(value,*,iostat=ios) job%model_params%tau_mem
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid tau_mem')
        case ('dt')
            read(value,*,iostat=ios) job%model_params%dt
            if (ios /= 0) call odf_config_read_error(unit, 'Invalid dt')
        case ('odf_path')
            job%io_params%odf_path = trim(value)
        case ('odf_prefix')
            job%io_params%odf_prefix = trim(value)
        case ('output_policy')
            job%io_params%output_policy = trim(value)
        case default
            call odf_config_read_error(unit, 'Unknown workflow parameter: '//trim(key))
        end select
    end do
    close(unit)
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete ODF sampling workflow file'
    call validate_odf_sampling_workflow(job)
end subroutine read_odf_sampling_workflow

subroutine print_odf_sampling_workflow(job)
    type(odf_sampling_workflow_type), intent(in) :: job
    write(output_unit,'(A,I0)') 'State ID: ', job%model_params%state
    write(output_unit,'(A,A)') 'CV plan: ', trim(job%io_params%cv_plan_filename)
    write(output_unit,'(A,A)') 'Lifting plan: ', trim(job%io_params%lifting_plan_filename)
    write(output_unit,'(A,A)') 'Lifted trajectory directory: ', trim(job%io_params%trajectory_path)
    write(output_unit,'(A,A)') 'Lifted trajectory prefix: ', trim(job%io_params%trajectory_prefix)
    write(output_unit,'(A,I0)') 'Maximum lag: ', job%model_params%max_lag
    write(output_unit,'(A,I0)') 'Basis index: ', job%model_params%basis_index
    write(output_unit,'(A,I0)') 'HEALPix levels: ', job%model_params%nlevel
    write(output_unit,'(A,ES14.6)') 'tau_mem: ', job%model_params%tau_mem
    write(output_unit,'(A,ES14.6)') 'dt: ', job%model_params%dt
    write(output_unit,'(A,A)') 'ODF output directory: ', trim(job%io_params%odf_path)
    write(output_unit,'(A,A)') 'ODF filename prefix: ', trim(job%io_params%odf_prefix)
    write(output_unit,'(A,A)') 'Output policy: ', trim(job%io_params%output_policy)
end subroutine print_odf_sampling_workflow

subroutine odf_config_read_error(unit, message)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    close(unit)
    write(error_unit,'(A)') trim(message)
    error stop 'ODF sampling workflow configuration read failed'
end subroutine odf_config_read_error

end module odf_sampling_workflow_module
