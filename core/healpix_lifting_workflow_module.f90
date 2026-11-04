module healpix_lifting_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use workflow_module, only: healpix_lifting_workflow_type, workflow_config_filename, &
        healpix_lifting_workflow_config_stem
    use lifting_plan_module, only: lifting_plan_type, read_lifting_plan
    use healpix_histogram_sampling_plan_module, only: lift_trajectories
    use iso_fortran_env, only: int64
    implicit none
    private

    public :: healpix_lifting_workflow_type
    public :: configure_healpix_lifting_workflow
    public :: read_healpix_lifting_workflow
    public :: write_healpix_lifting_workflow
    public :: validate_healpix_lifting_workflow
    public :: run_healpix_lifting_workflow

contains

subroutine configure_healpix_lifting_workflow(job)
    type(healpix_lifting_workflow_type), intent(out) :: job
    character(len=1024) :: input, config_filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' HEALPix Trajectory Lifting Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid HEALPix lifting workflow option'

    select case (choice)
    case (1)
        job = healpix_lifting_workflow_type()
        call edit_healpix_lifting_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_healpix_lifting_workflow(trim(config_filename), job)
        if (choice == 2) then
            call edit_healpix_lifting_workflow(job)
        else
            call print_healpix_lifting_workflow(job)
        end if
    case (4)
        job = healpix_lifting_workflow_type()
        write(output_unit,*) 'HEALPix lifting workflow configuration skipped.'
        return
    case default
        error stop 'Invalid HEALPix lifting workflow option'
    end select

    call validate_healpix_lifting_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        config_filename = workflow_config_filename(trim(input), healpix_lifting_workflow_config_stem)
        call write_healpix_lifting_workflow(job, trim(config_filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(config_filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_healpix_lifting_workflow

subroutine edit_healpix_lifting_workflow(job)
    type(healpix_lifting_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_integer('State ID', job%model_params%state)
    call edit_string('Lifting plan filename', job%io_params%lifting_plan_filename)
    call edit_string('Input trajectory directory', job%io_params%trajectory_path)
    call edit_string('Input trajectory filename prefix', job%io_params%trajectory_prefix)
    call edit_string('HEALPix lifted output directory', job%io_params%lifted_path)
    call edit_string('HEALPix lifted output filename prefix', job%io_params%lifted_prefix)
    call edit_integer('HEALPix pixel count (npix)', job%model_params%npix)
    call edit_integer('Trajectory lag', job%model_params%lag)
    call edit_integer('Backbone basis index', job%model_params%basis_index)

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
    call edit_string('Output policy (ask/replace/new)', job%io_params%output_policy)
    call validate_healpix_lifting_workflow(job)
    call print_healpix_lifting_workflow(job)

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
end subroutine edit_healpix_lifting_workflow

subroutine validate_healpix_lifting_workflow(job)
    type(healpix_lifting_workflow_type), intent(in) :: job
    integer :: nside

    if (job%model_params%state < 1) error stop 'HEALPix lifting state ID must be positive'
    if (job%model_params%npix < 12) error stop 'HEALPix npix must be at least 12'
    nside = int(sqrt(real(job%model_params%npix, 8) / 12.0_8))
    if (12_int64 * int(nside, int64)**2 /= int(job%model_params%npix, int64)) then
        error stop 'npix must equal 12*nside**2'
    end if
    if (iand(nside, nside - 1) /= 0) error stop 'HEALPix nside must be a power of two'
    if (job%model_params%lag < 0) error stop 'HEALPix lifting lag cannot be negative'
    if (job%model_params%basis_index < 1) error stop 'HEALPix lifting basis index must be positive'
    if (job%model_params%tau_mem <= 0.0_8) error stop 'HEALPix lifting tau_mem must be positive'
    if (job%model_params%dt <= 0.0_8) error stop 'HEALPix lifting dt must be positive'
    if (len_trim(job%io_params%lifting_plan_filename) == 0) error stop 'Lifting plan filename cannot be empty'
    if (len_trim(job%io_params%trajectory_path) == 0) error stop 'Input trajectory path cannot be empty'
    if (len_trim(job%io_params%trajectory_prefix) == 0) error stop 'Input trajectory prefix cannot be empty'
    if (len_trim(job%io_params%lifted_path) == 0) error stop 'Lifted output path cannot be empty'
    if (len_trim(job%io_params%lifted_prefix) == 0) error stop 'Lifted output prefix cannot be empty'
    select case (trim(job%io_params%output_policy))
    case ('ask', 'replace', 'new')
    case default
        error stop 'HEALPix lifting output policy must be ask, replace, or new'
    end select
end subroutine validate_healpix_lifting_workflow

subroutine run_healpix_lifting_workflow(job)
    type(healpix_lifting_workflow_type), intent(in) :: job
    type(lifting_plan_type) :: full_plan, selected_plan
    integer :: state_index
    character(len=1024) :: trajectory_path, lifted_path

    call validate_healpix_lifting_workflow(job)
    call read_lifting_plan(full_plan, trim(job%io_params%lifting_plan_filename))
    state_index = find_lifting_state(full_plan, job%model_params%state)
    call select_lifting_plan_state(full_plan, state_index, selected_plan)

    trajectory_path = ensure_directory_separator(job%io_params%trajectory_path)
    lifted_path = ensure_directory_separator(job%io_params%lifted_path)
    call lift_trajectories(selected_plan, trim(trajectory_path), trim(job%io_params%trajectory_prefix), &
        trim(lifted_path), trim(job%io_params%lifted_prefix), job%model_params%npix, &
        job%model_params%lag, job%model_params%basis_index, job%model_params%tau_mem, &
        job%model_params%dt, trim(job%io_params%output_policy))
end subroutine run_healpix_lifting_workflow

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
    error stop 'HEALPix lifting state is unavailable'
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

subroutine write_healpix_lifting_workflow(job, filename)
    type(healpix_lifting_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_healpix_lifting_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write HEALPix lifting workflow: ', trim(filename)
        error stop 'HEALPix lifting workflow output failed'
    end if
    write(unit,'(A)') '# HEALPix trajectory lifting workflow'
    write(unit,'(A)') 'workflow_type = healpix_lifting'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,I0)') 'state = ', job%model_params%state
    write(unit,'(A,A)') 'lifting_plan_filename = ', trim(job%io_params%lifting_plan_filename)
    write(unit,'(A,A)') 'trajectory_path = ', trim(job%io_params%trajectory_path)
    write(unit,'(A,A)') 'trajectory_prefix = ', trim(job%io_params%trajectory_prefix)
    write(unit,'(A,A)') 'lifted_path = ', trim(job%io_params%lifted_path)
    write(unit,'(A,A)') 'lifted_prefix = ', trim(job%io_params%lifted_prefix)
    write(unit,'(A,I0)') 'npix = ', job%model_params%npix
    write(unit,'(A,I0)') 'lag = ', job%model_params%lag
    write(unit,'(A,I0)') 'basis_index = ', job%model_params%basis_index
    write(unit,'(A,ES24.16E3)') 'tau_mem = ', job%model_params%tau_mem
    write(unit,'(A,ES24.16E3)') 'dt = ', job%model_params%dt
    write(unit,'(A,A)') 'output_policy = ', trim(job%io_params%output_policy)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing HEALPix lifting workflow file'
end subroutine write_healpix_lifting_workflow

subroutine read_healpix_lifting_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(healpix_lifting_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1900) :: value
    integer :: unit, ios, equals_position, format_version
    logical :: saw_type, saw_version

    job = healpix_lifting_workflow_type()
    format_version = -1
    saw_type = .false.
    saw_version = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read HEALPix lifting workflow: ', trim(filename)
        error stop 'HEALPix lifting workflow input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) call config_read_error(unit, 'Error reading HEALPix lifting workflow')
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) call config_read_error(unit, 'Invalid workflow configuration line')
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'healpix_lifting') call config_read_error(unit, 'Unexpected workflow_type')
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0) call config_read_error(unit, 'Invalid format_version')
            if (format_version /= 1) call config_read_error(unit, 'Unsupported format_version')
            saw_version = .true.
        case ('state')
            read(value,*,iostat=ios) job%model_params%state
            if (ios /= 0) call config_read_error(unit, 'Invalid state')
        case ('lifting_plan_filename')
            job%io_params%lifting_plan_filename = trim(value)
        case ('trajectory_path')
            job%io_params%trajectory_path = trim(value)
        case ('trajectory_prefix')
            job%io_params%trajectory_prefix = trim(value)
        case ('lifted_path')
            job%io_params%lifted_path = trim(value)
        case ('lifted_prefix')
            job%io_params%lifted_prefix = trim(value)
        case ('npix')
            read(value,*,iostat=ios) job%model_params%npix
            if (ios /= 0) call config_read_error(unit, 'Invalid npix')
        case ('lag')
            read(value,*,iostat=ios) job%model_params%lag
            if (ios /= 0) call config_read_error(unit, 'Invalid lag')
        case ('basis_index')
            read(value,*,iostat=ios) job%model_params%basis_index
            if (ios /= 0) call config_read_error(unit, 'Invalid basis_index')
        case ('tau_mem')
            read(value,*,iostat=ios) job%model_params%tau_mem
            if (ios /= 0) call config_read_error(unit, 'Invalid tau_mem')
        case ('dt')
            read(value,*,iostat=ios) job%model_params%dt
            if (ios /= 0) call config_read_error(unit, 'Invalid dt')
        case ('output_policy')
            job%io_params%output_policy = trim(value)
        case default
            call config_read_error(unit, 'Unknown workflow parameter: '//trim(key))
        end select
    end do
    close(unit)
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete HEALPix lifting workflow file'
    call validate_healpix_lifting_workflow(job)
end subroutine read_healpix_lifting_workflow

subroutine print_healpix_lifting_workflow(job)
    type(healpix_lifting_workflow_type), intent(in) :: job
    write(output_unit,'(A,I0)') 'State ID: ', job%model_params%state
    write(output_unit,'(A,A)') 'Lifting plan: ', trim(job%io_params%lifting_plan_filename)
    write(output_unit,'(A,A)') 'Input trajectory directory: ', trim(job%io_params%trajectory_path)
    write(output_unit,'(A,A)') 'Input trajectory prefix: ', trim(job%io_params%trajectory_prefix)
    write(output_unit,'(A,A)') 'Lifted output directory: ', trim(job%io_params%lifted_path)
    write(output_unit,'(A,A)') 'Lifted output prefix: ', trim(job%io_params%lifted_prefix)
    write(output_unit,'(A,I0)') 'HEALPix pixel count: ', job%model_params%npix
    write(output_unit,'(A,I0)') 'Lag: ', job%model_params%lag
    write(output_unit,'(A,I0)') 'Basis index: ', job%model_params%basis_index
    write(output_unit,'(A,ES14.6)') 'tau_mem: ', job%model_params%tau_mem
    write(output_unit,'(A,ES14.6)') 'dt: ', job%model_params%dt
    write(output_unit,'(A,A)') 'Output policy: ', trim(job%io_params%output_policy)
end subroutine print_healpix_lifting_workflow

subroutine config_read_error(unit, message)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    close(unit)
    write(error_unit,'(A)') trim(message)
    error stop 'HEALPix lifting workflow configuration read failed'
end subroutine config_read_error

end module healpix_lifting_workflow_module
