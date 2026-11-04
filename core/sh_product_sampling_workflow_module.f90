module sh_product_sampling_workflow_module
    use iso_fortran_env, only: output_unit, error_unit, real64, int64
    use workflow_module, only: sh_product_sampling_workflow_type, workflow_config_filename, &
        sh_product_sampling_workflow_config_stem
    use lifting_plan_module, only: lifting_plan_type, read_lifting_plan
    use lifted_trajectory_io_module, only: lifted_trajectory, read_trajectory
    use lifting_plan_execution_module, only: make_filename_lifted_tra
    use replicate_selector_workflow_module, only: replica_selector_type, read_replica_selector_file
    use sh_lifting_module, only: sh_lift_type, sh_lift_initialize, sh_lift_finalize, &
        sh_lift_vectors_block, sh_lm_index
    use output_policy_module, only: get_output_prefix, check_output_policy, OUTPUT_WRITE, &
        OUTPUT_SKIP, OUTPUT_COMPUTE_ONLY, OUTPUT_ABORT
    implicit none
    private

    integer, parameter :: SH_MEAN_SINGLE = 1
    integer, parameter :: SH_MEAN_PAIR = 2
    integer, parameter :: SH_MEAN_TRIPLE = 3
    integer, parameter :: SH_MEAN_SINGLE_FUTURE = 4
    integer, parameter :: SH_MEAN_PAIR_ZERO_LAG = 5

    type :: sh_product_mean_record_type
        integer :: kind = 0
        integer :: channel1 = 0, channel2 = 0, channel3 = 0
        integer :: l1 = 0, m1 = 0, l2 = 0, m2 = 0, l3 = 0, m3 = 0
        real(real64) :: mean_re = 0.0_real64, mean_im = 0.0_real64
        integer(int64) :: n_samples = 0_int64
    end type sh_product_mean_record_type

    type :: sh_product_means_type
        integer :: format_version = 4
        integer :: state_index = 0
        integer :: lmax = -1
        integer :: lag = 0
        integer :: n_channels = 0
        integer(int64) :: n_records = 0_int64
        type(sh_product_mean_record_type), allocatable :: records(:)
    end type sh_product_means_type

    public :: sh_product_sampling_workflow_type
    public :: SH_MEAN_SINGLE, SH_MEAN_PAIR, SH_MEAN_TRIPLE
    public :: SH_MEAN_SINGLE_FUTURE, SH_MEAN_PAIR_ZERO_LAG
    public :: sh_product_mean_record_type
    public :: sh_product_means_type
    public :: configure_sh_product_sampling_workflow
    public :: read_sh_product_sampling_workflow
    public :: write_sh_product_sampling_workflow
    public :: validate_sh_product_sampling_workflow
    public :: run_sh_product_sampling_workflow
    public :: check_sh_product_sampling_inputs
    public :: read_sh_product_means
    public :: write_sh_product_means

contains

subroutine configure_sh_product_sampling_workflow(job)
    type(sh_product_sampling_workflow_type), intent(out) :: job
    character(len=1024) :: input, filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' Spherical-Harmonic Product Sampling'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid SH product sampling workflow option'

    select case (choice)
    case (1)
        job = sh_product_sampling_workflow_type()
        call edit_sh_product_sampling_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename
        if (len_trim(filename) == 0) error stop 'No configuration filename specified'
        call read_sh_product_sampling_workflow(trim(filename), job)
        if (choice == 2) then
            call edit_sh_product_sampling_workflow(job)
        else
            call print_sh_product_sampling_workflow(job)
        end if
    case (4)
        job = sh_product_sampling_workflow_type()
        write(output_unit,*) 'SH product sampling workflow configuration skipped.'
        return
    case default
        error stop 'Invalid SH product sampling workflow option'
    end select

    call validate_sh_product_sampling_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        filename = workflow_config_filename(trim(input), sh_product_sampling_workflow_config_stem)
        call write_sh_product_sampling_workflow(job, trim(filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_sh_product_sampling_workflow

subroutine edit_sh_product_sampling_workflow(job)
    type(sh_product_sampling_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios

    write(output_unit,*)
    write(output_unit,*) 'The state and active replicates are read from the selector file.'
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_string('Replicate selector filename', job%selector_filename)
    call edit_string('Lifting plan filename', job%lifting_plan_filename)
    call edit_string('Active-pair matrix filename', job%pair_active_filename)
    call edit_string('Euclidean lifted trajectory directory', job%trajectory_path)
    call edit_string('Euclidean lifted trajectory prefix', job%trajectory_prefix)
    call edit_integer('Maximum spherical-harmonic degree L', job%lmax)
    call edit_integer('Sampling lag in frames (within each block)', job%lag)
    call edit_string('Result directory', job%output_path)
    call edit_string('Result filename prefix', job%output_prefix)
    call edit_string('Output policy (ask/replace/new)', job%output_policy)
    call validate_sh_product_sampling_workflow(job)
    call print_sh_product_sampling_workflow(job)

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
end subroutine edit_sh_product_sampling_workflow

subroutine validate_sh_product_sampling_workflow(job)
    type(sh_product_sampling_workflow_type), intent(in) :: job

    if (job%lmax < 0 .or. job%lmax > 64) error stop 'SH maximum degree must be between 0 and 64'
    if (job%lag < 0) error stop 'SH product sampling lag cannot be negative'
    if (len_trim(job%selector_filename) == 0) error stop 'Replicate selector filename cannot be empty'
    if (len_trim(job%lifting_plan_filename) == 0) error stop 'Lifting plan filename cannot be empty'
    if (len_trim(job%pair_active_filename) == 0) error stop 'Active-pair matrix filename cannot be empty'
    if (len_trim(job%trajectory_path) == 0) error stop 'Trajectory path cannot be empty'
    if (len_trim(job%trajectory_prefix) == 0) error stop 'Trajectory prefix cannot be empty'
    if (len_trim(job%output_path) == 0) error stop 'Output path cannot be empty'
    if (len_trim(job%output_prefix) == 0) error stop 'Output prefix cannot be empty'
    select case (trim(job%output_policy))
    case ('ask', 'replace', 'new')
    case default
        error stop 'SH result output policy must be ask, replace, or new'
    end select
end subroutine validate_sh_product_sampling_workflow

subroutine run_sh_product_sampling_workflow(job)
    type(sh_product_sampling_workflow_type), intent(in) :: job
    type(lifting_plan_type) :: lifting_plan
    type(replica_selector_type) :: selector
    type(sh_lift_type) :: sh
    type(lifted_trajectory) :: trajectory
    logical, allocatable :: pair_active(:,:)
    integer :: state_index, n_channels, irep, iblock, unit, ios, output_action
    integer(int64) :: n_samples, n_records
    character(len=2048) :: output_filename
    character(len=1024) :: output_directory
    character(len=2048) :: output_name
    character(len=2048) :: trajectory_prefix, trajectory_filename

    call validate_sh_product_sampling_workflow(job)
    call read_replica_selector_file(trim(job%selector_filename), selector, state_index)
    call read_lifting_plan(lifting_plan, trim(job%lifting_plan_filename))
    call validate_selection(selector, state_index, lifting_plan)
    call check_sh_product_sampling_inputs(job)

    irep = first_active_replica(selector)
    iblock = first_block(lifting_plan, state_index, irep)
    trajectory_prefix = join_path(job%trajectory_path, job%trajectory_prefix)
    call make_filename_lifted_tra(state_index, irep, iblock, trim(trajectory_prefix), trajectory_filename)
    call read_trajectory(trim(trajectory_filename), trajectory)
    n_channels = 2 * trajectory%nobject
    if (n_channels < 1) error stop 'Trajectory contains no direction channels'
    call read_pair_active_matrix(trim(job%pair_active_filename), n_channels, pair_active)

    n_samples = count_lagged_samples(lifting_plan, selector, state_index, &
        trim(job%trajectory_path), trim(job%trajectory_prefix), job%lag)
    if (n_samples < 1_int64) error stop 'Selected blocks contain no samples for the configured lag'

    output_directory = trim(job%output_path)
    output_name = trim(job%output_prefix)
    if (output_directory(len_trim(output_directory):len_trim(output_directory)) /= '/') then
        output_directory = trim(output_directory)//'/'
    end if
    output_name = get_output_prefix(trim(output_name), trim(job%output_policy))
    output_filename = trim(output_directory)//trim(output_name)//'_means.txt'
    call check_output_policy(trim(output_filename), trim(job%output_policy), output_action)
    if (output_action == OUTPUT_ABORT) error stop 'SH product sampling output policy aborted workflow'
    if (output_action == OUTPUT_SKIP) return

    n_records = count_output_records(pair_active, job%lmax)
    write(output_unit,'(A,I0,A,I0,A,I0,A)') 'Sampling ', n_samples, ' lagged frames (lag ', &
        job%lag, ') into ', n_records, ' selected SH records.'
    call sh_lift_initialize(sh, 2 * job%lmax)
    if (output_action == OUTPUT_WRITE) then
        open(newunit=unit, file=trim(output_filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Cannot write SH product means: ', trim(output_filename)
            error stop 'SH product means output failed'
        end if
        call write_means_header(unit, state_index, job%lmax, job%lag, n_channels, n_records)
    else
        unit = 0
    end if

    call sample_single_means(unit, output_action, lifting_plan, selector, state_index, &
        trim(job%trajectory_path), trim(job%trajectory_prefix), sh, job%lmax, &
        job%lag, pair_active, n_channels, n_samples)
    call sample_pair_means(unit, output_action, lifting_plan, selector, state_index, &
        trim(job%trajectory_path), trim(job%trajectory_prefix), sh, job%lmax, &
        job%lag, pair_active, n_channels, n_samples)
    call sample_triple_means(unit, output_action, lifting_plan, selector, state_index, &
        trim(job%trajectory_path), trim(job%trajectory_prefix), sh, job%lmax, &
        job%lag, pair_active, n_channels, n_samples)

    if (output_action == OUTPUT_WRITE) then
        close(unit, iostat=ios)
        if (ios /= 0) error stop 'Error closing SH product means file'
        write(output_unit,'(3A)') 'SH product means saved to "', trim(output_filename), '"'
    end if
    call sh_lift_finalize(sh)
end subroutine run_sh_product_sampling_workflow

subroutine check_sh_product_sampling_inputs(job, selector_will_be_generated, plan_will_be_generated)
    type(sh_product_sampling_workflow_type), intent(in) :: job
    logical, intent(in), optional :: selector_will_be_generated, plan_will_be_generated
    type(replica_selector_type) :: selector
    type(lifting_plan_type) :: lifting_plan
    type(lifted_trajectory) :: trajectory
    logical, allocatable :: pair_active(:,:)
    integer :: state_index, irep, iblock, n_channels
    integer(int64) :: n_samples
    logical :: selector_generated, plan_generated, exists

    call validate_sh_product_sampling_workflow(job)
    selector_generated = .false.
    plan_generated = .false.
    if (present(selector_will_be_generated)) selector_generated = selector_will_be_generated
    if (present(plan_will_be_generated)) plan_generated = plan_will_be_generated

    inquire(file=trim(job%selector_filename), exist=exists)
    if (.not. exists .and. .not. selector_generated) then
        call missing_file('Replicate selector', trim(job%selector_filename))
    end if
    inquire(file=trim(job%lifting_plan_filename), exist=exists)
    if (.not. exists .and. .not. plan_generated) then
        call missing_file('Lifting plan', trim(job%lifting_plan_filename))
    end if
    inquire(file=trim(job%pair_active_filename), exist=exists)
    if (.not. exists) call missing_file('Active-pair matrix', trim(job%pair_active_filename))
    if (selector_generated .or. plan_generated) return

    call read_replica_selector_file(trim(job%selector_filename), selector, state_index)
    call read_lifting_plan(lifting_plan, trim(job%lifting_plan_filename))
    call validate_selection(selector, state_index, lifting_plan)
    call check_selected_trajectory_files(lifting_plan, selector, state_index, &
        trim(job%trajectory_path), trim(job%trajectory_prefix))
    irep = first_active_replica(selector)
    iblock = first_block(lifting_plan, state_index, irep)
    block
        character(len=2048) :: prefix, filename
        prefix = join_path(job%trajectory_path, job%trajectory_prefix)
        call make_filename_lifted_tra(state_index, irep, iblock, trim(prefix), filename)
        call read_trajectory(trim(filename), trajectory)
    end block
    if (.not. allocated(trajectory%vec)) error stop 'Euclidean trajectory block has no vectors'
    n_channels = size(trajectory%vec,2)
    if (n_channels /= 2*trajectory%nobject) error stop 'Trajectory direction-channel count is inconsistent'
    call read_pair_active_matrix(trim(job%pair_active_filename), n_channels, pair_active)
    if (.not. any(pair_active)) error stop 'Active-pair matrix selects no channel combinations'
    n_samples = count_lagged_samples(lifting_plan, selector, state_index, &
        trim(job%trajectory_path), trim(job%trajectory_prefix), job%lag)
    if (n_samples < 1_int64) error stop 'Selected blocks contain no samples for the configured lag'
end subroutine check_sh_product_sampling_inputs

subroutine validate_selection(selector, state_index, lifting_plan)
    type(replica_selector_type), intent(in) :: selector
    integer, intent(in) :: state_index
    type(lifting_plan_type), intent(in) :: lifting_plan

    if (.not. allocated(selector%is_active)) error stop 'Replicate selector has no activation array'
    if (.not. allocated(lifting_plan%istates) .or. .not. allocated(lifting_plan%nreps) .or. &
        .not. allocated(lifting_plan%nblock)) error stop 'Lifting plan is incomplete'
    if (selector%nrep_active < 1) error stop 'Replicate selector has no active replicates'
    if (size(lifting_plan%nreps) /= size(lifting_plan%istates) .or. &
        size(lifting_plan%nblock,1) /= size(lifting_plan%istates)) then
        error stop 'Lifting plan state arrays are inconsistent'
    end if
    if (state_index < 1 .or. state_index > size(lifting_plan%istates)) then
        error stop 'Replicate selector state is outside the lifting plan'
    end if
    if (lifting_plan%istates(state_index) /= state_index) then
        error stop 'Replicate selector state does not match the lifting plan'
    end if
    if (lifting_plan%nreps(state_index) /= selector%nrep_total) then
        error stop 'Replicate selector and lifting plan have different replica counts'
    end if
    if (size(selector%is_active) /= selector%nrep_total) then
        error stop 'Replicate selector array size is inconsistent'
    end if
    if (size(lifting_plan%nblock,2) < selector%nrep_total) then
        error stop 'Lifting plan replica array is smaller than the selector'
    end if
    if (count(selector%is_active) /= selector%nrep_active) then
        error stop 'Replicate selector active-replica count is inconsistent'
    end if
    block
        integer :: irep
        do irep = 1, selector%nrep_total
            if (selector%is_active(irep) .and. lifting_plan%nblock(state_index,irep) < 1) then
                error stop 'An active replica has no lifting blocks in the plan'
            end if
        end do
    end block
end subroutine validate_selection

subroutine check_selected_trajectory_files(lifting_plan, selector, state_index, trajectory_path, trajectory_prefix)
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(replica_selector_type), intent(in) :: selector
    integer, intent(in) :: state_index
    character(len=*), intent(in) :: trajectory_path, trajectory_prefix
    integer :: irep, iblock
    character(len=2048) :: prefix, filename
    logical :: exists

    prefix = join_path(trajectory_path, trajectory_prefix)
    do irep = 1, selector%nrep_total
        if (.not. selector%is_active(irep)) cycle
        do iblock = 1, lifting_plan%nblock(state_index,irep)
            call make_filename_lifted_tra(state_index, irep, iblock, trim(prefix), filename)
            inquire(file=trim(filename), exist=exists)
            if (.not. exists) call missing_file('Euclidean lifted trajectory', trim(filename))
        end do
    end do
end subroutine check_selected_trajectory_files

integer(int64) function count_lagged_samples(lifting_plan, selector, state_index, trajectory_path, &
    trajectory_prefix, lag) result(n_samples)
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(replica_selector_type), intent(in) :: selector
    integer, intent(in) :: state_index, lag
    character(len=*), intent(in) :: trajectory_path, trajectory_prefix
    integer :: irep, iblock, unit, ios, nobject, ntimes
    character(len=2048) :: prefix, filename

    n_samples = 0_int64
    prefix = join_path(trajectory_path, trajectory_prefix)
    do irep = 1, selector%nrep_total
        if (.not. selector%is_active(irep)) cycle
        do iblock = 1, lifting_plan%nblock(state_index,irep)
            call make_filename_lifted_tra(state_index, irep, iblock, trim(prefix), filename)
            open(newunit=unit, file=trim(filename), access='stream', form='unformatted', &
                status='old', action='read', iostat=ios)
            if (ios /= 0) call missing_file('Cannot open trajectory', trim(filename))
            read(unit, iostat=ios) nobject
            if (ios == 0) read(unit, iostat=ios) ntimes
            close(unit)
            if (ios /= 0) call missing_file('Invalid trajectory header', trim(filename))
            if (nobject < 1 .or. ntimes < 1) call missing_file('Empty trajectory block', trim(filename))
            if (ntimes > lag) n_samples = n_samples + int(ntimes-lag, int64)
        end do
    end do
end function count_lagged_samples

subroutine read_pair_active_matrix(filename, n_channels, pair_active)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n_channels
    logical, allocatable, intent(out) :: pair_active(:,:)
    integer :: unit, ios, channel
    character(len=2048) :: line

    allocate(pair_active(n_channels,n_channels))
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call missing_file('Cannot read active-pair matrix', trim(filename))
    do channel = 1, n_channels
        read(unit,*,iostat=ios) pair_active(channel,:)
        if (ios /= 0) then
            close(unit)
            write(error_unit,'(A,I0,A,A)') 'Invalid or incomplete active-pair matrix row ', &
                channel, ' in ', trim(filename)
            error stop 'Active-pair matrix input failed'
        end if
    end do
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading active-pair matrix'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        close(unit)
        error stop 'Active-pair matrix contains more rows than trajectory channels'
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing active-pair matrix'
end subroutine read_pair_active_matrix

integer function first_active_replica(selector) result(irep)
    type(replica_selector_type), intent(in) :: selector
    do irep = 1, selector%nrep_total
        if (selector%is_active(irep)) return
    end do
    error stop 'No active replicate in selector'
end function first_active_replica

integer function first_block(lifting_plan, state_index, irep) result(iblock)
    type(lifting_plan_type), intent(in) :: lifting_plan
    integer, intent(in) :: state_index, irep
    if (lifting_plan%nblock(state_index,irep) < 1) error stop 'Active replicate has no lifting blocks'
    iblock = 1
end function first_block

subroutine load_sh_channels(filename, channels, n_channels, sh, Y_re, Y_im)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: channels(:)
    integer, intent(in) :: n_channels
    type(sh_lift_type), intent(in) :: sh
    real(real64), allocatable, intent(out) :: Y_re(:,:,:), Y_im(:,:,:)
    type(lifted_trajectory) :: trajectory
    real(real64), allocatable :: vectors(:,:,:)
    integer :: i, nframes

    call read_trajectory(filename, trajectory)
    if (.not. allocated(trajectory%vec)) error stop 'Euclidean trajectory block has no vectors'
    if (size(trajectory%vec,1) /= 3 .or. size(trajectory%vec,2) /= n_channels .or. &
        size(trajectory%vec,2) /= 2*trajectory%nobject .or. &
        size(trajectory%vec,3) /= trajectory%ntimes) then
        error stop 'Euclidean trajectory block has inconsistent dimensions'
    end if
    if (maxval(channels) > size(trajectory%vec,2)) then
        error stop 'Requested SH channel is outside the Euclidean trajectory'
    end if
    nframes = trajectory%ntimes
    allocate(vectors(3,size(channels),nframes))
    do i = 1, size(channels)
        vectors(:,i,:) = trajectory%vec(:,channels(i),:)
    end do
    allocate(Y_re(sh%nlm,size(channels),nframes), Y_im(sh%nlm,size(channels),nframes))
    call sh_lift_vectors_block(sh, vectors, Y_re, Y_im)
end subroutine load_sh_channels

subroutine sample_single_means(unit, output_action, lifting_plan, selector, state_index, &
    trajectory_path, trajectory_prefix, sh, lmax, lag, pair_active, n_channels, n_samples)
    integer, intent(in) :: unit, output_action, state_index, lmax, lag, n_channels
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(replica_selector_type), intent(in) :: selector
    character(len=*), intent(in) :: trajectory_path, trajectory_prefix
    type(sh_lift_type), intent(in) :: sh
    logical, intent(in) :: pair_active(:,:)
    integer(int64), intent(in) :: n_samples
    integer :: channel, irep, iblock, frame, lm, l, m
    integer, allocatable :: channels(:)
    logical, allocatable :: selected_channel(:)
    real(real64), allocatable :: Y_re(:,:,:), Y_im(:,:,:)
    real(real64), allocatable :: current_re(:), current_im(:), future_re(:), future_im(:)
    character(len=2048) :: prefix, filename

    allocate(channels(1), selected_channel(n_channels), &
        current_re(sh%nlm), current_im(sh%nlm), &
        future_re(sh%nlm), future_im(sh%nlm))
    selected_channel = any(pair_active,dim=1) .or. any(pair_active,dim=2)
    prefix = join_path(trajectory_path, trajectory_prefix)
    do channel = 1, n_channels
        if (.not. selected_channel(channel)) cycle
        channels(1) = channel
        current_re = 0.0_real64
        current_im = 0.0_real64
        future_re = 0.0_real64
        future_im = 0.0_real64
        do irep = 1, selector%nrep_total
            if (.not. selector%is_active(irep)) cycle
            do iblock = 1, lifting_plan%nblock(state_index,irep)
                call make_filename_lifted_tra(state_index, irep, iblock, trim(prefix), filename)
                call load_sh_channels(trim(filename), channels, n_channels, sh, Y_re, Y_im)
                do frame = 1, size(Y_re,3)-lag
                    do l = 0, 2*lmax
                        do m = -l, l
                            lm = sh_lm_index(l,m)
                            current_re(lm) = current_re(lm) + Y_re(lm,1,frame)
                            current_im(lm) = current_im(lm) + Y_im(lm,1,frame)
                            future_re(lm) = future_re(lm) + Y_re(lm,1,frame+lag)
                            future_im(lm) = future_im(lm) + Y_im(lm,1,frame+lag)
                        end do
                    end do
                end do
                deallocate(Y_re, Y_im)
            end do
        end do
        do l = 0, 2*lmax
            do m = -l, l
                lm = sh_lm_index(l,m)
                call write_mean_record(unit, output_action, SH_MEAN_SINGLE, channel, 0, 0, &
                    l, m, 0, 0, 0, 0, current_re(lm)/real(n_samples,real64), &
                    current_im(lm)/real(n_samples,real64), n_samples)
                call write_mean_record(unit, output_action, SH_MEAN_SINGLE_FUTURE, channel, 0, 0, &
                    l, m, 0, 0, 0, 0, future_re(lm)/real(n_samples,real64), &
                    future_im(lm)/real(n_samples,real64), n_samples)
            end do
        end do
    end do
end subroutine sample_single_means

subroutine sample_pair_means(unit, output_action, lifting_plan, selector, state_index, &
    trajectory_path, trajectory_prefix, sh, lmax, lag, pair_active, n_channels, n_samples)
    integer, intent(in) :: unit, output_action, state_index, lmax, lag, n_channels
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(replica_selector_type), intent(in) :: selector
    character(len=*), intent(in) :: trajectory_path, trajectory_prefix
    type(sh_lift_type), intent(in) :: sh
    logical, intent(in) :: pair_active(:,:)
    integer(int64), intent(in) :: n_samples
    integer :: channel1, channel2, irep, iblock, frame, lm1, lm2, l1, m1, l2, m2
    integer :: channels(2), n_lm, n_lm_first
    logical, allocatable :: lag_pairs(:,:)
    real(real64), allocatable :: Y_re(:,:,:), Y_im(:,:,:)
    real(real64), allocatable :: lag_sum_re(:,:), lag_sum_im(:,:)
    real(real64), allocatable :: zero_sum_re(:,:), zero_sum_im(:,:)
    real(real64) :: product_re, product_im
    character(len=2048) :: prefix, filename

    n_lm = sh%nlm
    n_lm_first = (lmax+1)**2
    allocate(lag_pairs(n_channels,n_channels))
    lag_pairs = pair_active
    do channel1 = 1, n_channels
        if (any(pair_active(channel1,:))) lag_pairs(channel1,channel1) = .true.
    end do
    allocate(lag_sum_re(n_lm_first,n_lm), lag_sum_im(n_lm_first,n_lm), &
        zero_sum_re(n_lm_first,n_lm), zero_sum_im(n_lm_first,n_lm))
    prefix = join_path(trajectory_path, trajectory_prefix)
    do channel1 = 1, n_channels
        do channel2 = 1, n_channels
            if (.not. lag_pairs(channel1,channel2) .and. .not. pair_active(channel1,channel2)) cycle
            channels = [channel1, channel2]
            lag_sum_re = 0.0_real64
            lag_sum_im = 0.0_real64
            zero_sum_re = 0.0_real64
            zero_sum_im = 0.0_real64
            do irep = 1, selector%nrep_total
                if (.not. selector%is_active(irep)) cycle
                do iblock = 1, lifting_plan%nblock(state_index,irep)
                    call make_filename_lifted_tra(state_index, irep, iblock, trim(prefix), filename)
                    call load_sh_channels(trim(filename), channels, n_channels, sh, Y_re, Y_im)
                    do frame = 1, size(Y_re,3)-lag
                        do l1 = 0, lmax
                            do m1 = -l1, l1
                                lm1 = sh_lm_index(l1,m1)
                                do l2 = 0, 2*lmax
                                    do m2 = -l2, l2
                                        lm2 = sh_lm_index(l2,m2)
                                        if (lag_pairs(channel1,channel2)) then
                                            product_re = Y_re(lm1,1,frame+lag)*Y_re(lm2,2,frame) - &
                                                Y_im(lm1,1,frame+lag)*Y_im(lm2,2,frame)
                                            product_im = Y_re(lm1,1,frame+lag)*Y_im(lm2,2,frame) + &
                                                Y_im(lm1,1,frame+lag)*Y_re(lm2,2,frame)
                                            lag_sum_re(lm1,lm2) = lag_sum_re(lm1,lm2) + product_re
                                            lag_sum_im(lm1,lm2) = lag_sum_im(lm1,lm2) + product_im
                                        end if
                                        if (pair_active(channel1,channel2)) then
                                            product_re = Y_re(lm1,1,frame)*Y_re(lm2,2,frame) - &
                                                Y_im(lm1,1,frame)*Y_im(lm2,2,frame)
                                            product_im = Y_re(lm1,1,frame)*Y_im(lm2,2,frame) + &
                                                Y_im(lm1,1,frame)*Y_re(lm2,2,frame)
                                            zero_sum_re(lm1,lm2) = zero_sum_re(lm1,lm2) + product_re
                                            zero_sum_im(lm1,lm2) = zero_sum_im(lm1,lm2) + product_im
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                    deallocate(Y_re, Y_im)
                end do
            end do
            if (lag_pairs(channel1,channel2)) then
                do l1 = 0, lmax
                    do m1 = -l1, l1
                        lm1 = sh_lm_index(l1,m1)
                        do l2 = 0, 2*lmax
                            do m2 = -l2, l2
                                lm2 = sh_lm_index(l2,m2)
                                call write_mean_record(unit, output_action, SH_MEAN_PAIR, &
                                    channel1, channel2, 0, l1, m1, l2, m2, 0, 0, &
                                    lag_sum_re(lm1,lm2)/real(n_samples,real64), &
                                    lag_sum_im(lm1,lm2)/real(n_samples,real64), n_samples)
                            end do
                        end do
                    end do
                end do
            end if
            if (pair_active(channel1,channel2)) then
                do l1 = 0, lmax
                    do m1 = -l1, l1
                        lm1 = sh_lm_index(l1,m1)
                        do l2 = 0, 2*lmax
                            do m2 = -l2, l2
                                lm2 = sh_lm_index(l2,m2)
                                call write_mean_record(unit, output_action, SH_MEAN_PAIR_ZERO_LAG, &
                                    channel1, channel2, 0, l1, m1, l2, m2, 0, 0, &
                                    zero_sum_re(lm1,lm2)/real(n_samples,real64), &
                                    zero_sum_im(lm1,lm2)/real(n_samples,real64), n_samples)
                            end do
                        end do
                    end do
                end do
            end if
        end do
    end do
end subroutine sample_pair_means

subroutine sample_triple_means(unit, output_action, lifting_plan, selector, state_index, &
    trajectory_path, trajectory_prefix, sh, lmax, lag, pair_active, n_channels, n_samples)
    integer, intent(in) :: unit, output_action, state_index, lmax, lag, n_channels
    type(lifting_plan_type), intent(in) :: lifting_plan
    type(replica_selector_type), intent(in) :: selector
    character(len=*), intent(in) :: trajectory_path, trajectory_prefix
    type(sh_lift_type), intent(in) :: sh
    logical, intent(in) :: pair_active(:,:)
    integer(int64), intent(in) :: n_samples
    integer :: channel1, channel2, irep, iblock, frame
    integer :: l1, m1, l2, m2, l3, m3, lm1, lm2, lm3, index
    integer :: channels(3)
    integer(int64) :: n_components
    real(real64), allocatable :: Y_re(:,:,:), Y_im(:,:,:), sum_re(:), sum_im(:)
    real(real64) :: product_re, product_im, partial_re, partial_im
    character(len=2048) :: prefix, filename

    n_components = count_triple_components(lmax)
    if (n_components > int(huge(index),int64)) error stop 'Too many SH triple products to allocate'
    allocate(sum_re(int(n_components)), sum_im(int(n_components)))
    prefix = join_path(trajectory_path, trajectory_prefix)
    do channel1 = 1, n_channels
        do channel2 = 1, n_channels
            if (.not. pair_active(channel1,channel2)) cycle
                channels = [channel1, channel1, channel2]
                sum_re = 0.0_real64
                sum_im = 0.0_real64
                do irep = 1, selector%nrep_total
                    if (.not. selector%is_active(irep)) cycle
                    do iblock = 1, lifting_plan%nblock(state_index,irep)
                        call make_filename_lifted_tra(state_index, irep, iblock, trim(prefix), filename)
                        call load_sh_channels(trim(filename), channels, n_channels, sh, Y_re, Y_im)
                        do frame = 1, size(Y_re,3)-lag
                            index = 0
                            do l1 = 0, lmax
                                do m1 = -l1, l1
                                    lm1 = sh_lm_index(l1,m1)
                                    do l2 = 0, lmax
                                        do m2 = -l2, l2
                                            lm2 = sh_lm_index(l2,m2)
                                            m3 = -m1 - m2
                                            do l3 = abs(l1-l2), l1+l2
                                                if (abs(m3) > l3) cycle
                                                index = index + 1
                                                lm3 = sh_lm_index(l3,m3)
                                                partial_re = Y_re(lm1,1,frame+lag)*Y_re(lm2,2,frame) - &
                                                    Y_im(lm1,1,frame+lag)*Y_im(lm2,2,frame)
                                                partial_im = Y_re(lm1,1,frame+lag)*Y_im(lm2,2,frame) + &
                                                    Y_im(lm1,1,frame+lag)*Y_re(lm2,2,frame)
                                                product_re = partial_re*Y_re(lm3,3,frame) - &
                                                    partial_im*Y_im(lm3,3,frame)
                                                product_im = partial_re*Y_im(lm3,3,frame) + &
                                                    partial_im*Y_re(lm3,3,frame)
                                                sum_re(index) = sum_re(index) + product_re
                                                sum_im(index) = sum_im(index) + product_im
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                            if (index /= n_components) error stop 'Internal SH triple component count mismatch'
                        end do
                        deallocate(Y_re, Y_im)
                    end do
                end do

                index = 0
                do l1 = 0, lmax
                    do m1 = -l1, l1
                        do l2 = 0, lmax
                            do m2 = -l2, l2
                                m3 = -m1 - m2
                                do l3 = abs(l1-l2), l1+l2
                                    if (abs(m3) > l3) cycle
                                    index = index + 1
                                    call write_mean_record(unit, output_action, SH_MEAN_TRIPLE, &
                                        channel1, channel1, channel2, l1, m1, l2, m2, l3, m3, &
                                        sum_re(index)/real(n_samples,real64), &
                                        sum_im(index)/real(n_samples,real64), n_samples)
                                end do
                            end do
                        end do
                    end do
                end do
        end do
    end do
end subroutine sample_triple_means

subroutine write_mean_record(unit, output_action, kind, channel1, channel2, channel3, &
    l1, m1, l2, m2, l3, m3, mean_re, mean_im, n_samples)
    integer, intent(in) :: unit, output_action, kind
    integer, intent(in) :: channel1, channel2, channel3, l1, m1, l2, m2, l3, m3
    real(real64), intent(in) :: mean_re, mean_im
    integer(int64), intent(in) :: n_samples

    if (output_action /= OUTPUT_WRITE) return
    write(unit,*) kind, channel1, channel2, channel3, l1, m1, l2, m2, l3, m3, &
        mean_re, mean_im, n_samples
end subroutine write_mean_record

integer(int64) function count_output_records(pair_active, lmax) result(n_records)
    logical, intent(in) :: pair_active(:,:)
    integer, intent(in) :: lmax
    integer(int64) :: n_lm, n_lm_full, n_pair_components, n_triple_components
    integer(int64) :: n_single, n_lag_pair, n_zero_pair, n_triple, record_limit
    integer(int64) :: n_active_channels, n_active_pairs
    integer :: channel
    logical, allocatable :: lag_pairs(:,:)

    if (size(pair_active,1) /= size(pair_active,2)) error stop 'Active-pair matrix must be square'
    n_active_channels = int(count(any(pair_active,dim=1) .or. any(pair_active,dim=2)),int64)
    n_active_pairs = int(count(pair_active),int64)
    allocate(lag_pairs(size(pair_active,1),size(pair_active,2)))
    lag_pairs = pair_active
    do channel = 1, size(pair_active,1)
        if (any(pair_active(channel,:))) lag_pairs(channel,channel) = .true.
    end do
    n_lag_pair = int(count(lag_pairs),int64)
    n_lm = int(lmax+1,int64)**2
    n_lm_full = int(2*lmax+1,int64)**2
    n_pair_components = n_lm*n_lm_full
    n_triple_components = count_triple_components(lmax)
    record_limit = int(huge(0),int64)
    n_single = checked_record_product(2_int64*n_active_channels, n_lm_full, record_limit)
    n_zero_pair = checked_record_product(n_active_pairs, n_pair_components, record_limit)
    n_lag_pair = checked_record_product(n_lag_pair, n_pair_components, record_limit)
    n_triple = checked_record_product(n_active_pairs, n_triple_components, record_limit)
    if (n_single > record_limit-n_zero_pair .or. n_single+n_zero_pair > record_limit-n_lag_pair .or. &
        n_single+n_zero_pair+n_lag_pair > record_limit-n_triple) then
        error stop 'SH product output exceeds the supported record limit; reduce L or channel count'
    end if
    n_records = n_single + n_zero_pair + n_lag_pair + n_triple
end function count_output_records

integer(int64) function count_triple_components(lmax) result(n_components)
    integer, intent(in) :: lmax
    integer :: l1, l2, m3, l3_min, multiplicity
    integer :: m1_min, m1_max
    n_components = 0_int64
    do l1 = 0, lmax
        do l2 = 0, lmax
            do m3 = -l1-l2, l1+l2
                m1_min = max(-l1, m3-l2)
                m1_max = min(l1, m3+l2)
                multiplicity = max(0, m1_max-m1_min+1)
                l3_min = max(abs(l1-l2), abs(m3))
                if (multiplicity > 0 .and. l3_min <= l1+l2) then
                    n_components = n_components + int(multiplicity,int64) * &
                        int(l1+l2-l3_min+1,int64)
                end if
            end do
        end do
    end do
end function count_triple_components

integer(int64) function checked_record_product(a, b, limit) result(product)
    integer(int64), intent(in) :: a, b, limit
    if (a < 0_int64 .or. b < 0_int64) error stop 'Invalid SH product record-count factor'
    if (b /= 0_int64 .and. a > limit / b) then
        error stop 'SH product output exceeds the supported record limit; reduce L or channel count'
    end if
    product = a*b
end function checked_record_product

subroutine write_means_header(unit, state_index, lmax, lag, n_channels, n_records)
    integer, intent(in) :: unit, state_index, lmax, lag, n_channels
    integer(int64), intent(in) :: n_records
    write(unit,'(A)') '# Lagged spherical-harmonic product means'
    write(unit,'(A)') 'format_version = 4'
    write(unit,'(A,I0)') 'state_index = ', state_index
    write(unit,'(A,I0)') 'lmax = ', lmax
    write(unit,'(A,I0)') 'lag = ', lag
    write(unit,'(A,I0)') 'n_channels = ', n_channels
    write(unit,'(A,I0)') 'n_records = ', n_records
    write(unit,'(A)') '# kind channels(3) l1 m1 l2 m2 l3 m3 mean_re mean_im n_samples'
    write(unit,'(A)') '# pair: Y(channel1,l1,m1,t+lag) * Y(channel2,l2,m2,t)'
    write(unit,'(A)') '# triple: Y(channel1,l1,m1,t+lag) * Y(channel2,l2,m2,t) * Y(channel3,l3,m3,t)'
    write(unit,'(A)') '# triple modes obey m1 + m2 + m3 = 0'
    write(unit,'(A)') '# selected pair (i,j) maps to triple channels (i,i,j)'
    write(unit,'(A)') '# kind 1=single(t), 4=single(t+lag), 2=pair(lag), 5=pair(t), 3=triple'
end subroutine write_means_header

subroutine write_sh_product_means(results, filename)
    type(sh_product_means_type), intent(in) :: results
    character(len=*), intent(in) :: filename
    integer :: unit, ios, i

    if (.not. allocated(results%records)) error stop 'SH means records are not allocated'
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call missing_file('Cannot write SH product means', trim(filename))
    call write_means_header(unit, results%state_index, results%lmax, results%lag, &
        results%n_channels, int(size(results%records),int64))
    do i = 1, size(results%records)
        write(unit,*) results%records(i)%kind, results%records(i)%channel1, &
            results%records(i)%channel2, results%records(i)%channel3, &
            results%records(i)%l1, results%records(i)%m1, results%records(i)%l2, &
            results%records(i)%m2, results%records(i)%l3, results%records(i)%m3, &
            results%records(i)%mean_re, results%records(i)%mean_im, results%records(i)%n_samples
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing SH product means file'
end subroutine write_sh_product_means

subroutine read_sh_product_means(filename, results)
    character(len=*), intent(in) :: filename
    type(sh_product_means_type), intent(out) :: results
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1800) :: value
    integer :: unit, ios, equals_position, i, alloc_status
    integer :: format_version
    integer(int64) :: n_records
    logical :: saw_version, saw_state, saw_lmax, saw_lag, saw_channels, saw_records

    n_records = -1_int64
    format_version = -1
    saw_version = .false.
    saw_state = .false.
    saw_lmax = .false.
    saw_lag = .false.
    saw_channels = .false.
    saw_records = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call missing_file('Cannot read SH product means', trim(filename))
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading SH product means header'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) exit
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0 .or. format_version < 1 .or. format_version > 4) &
                error stop 'Unsupported SH product means format'
            saw_version = .true.
        case ('state_index')
            read(value,*,iostat=ios) results%state_index
            if (ios /= 0) error stop 'Invalid SH means state index'
            saw_state = .true.
        case ('lmax')
            read(value,*,iostat=ios) results%lmax
            if (ios /= 0) error stop 'Invalid SH means lmax'
            saw_lmax = .true.
        case ('lag')
            read(value,*,iostat=ios) results%lag
            if (ios /= 0) error stop 'Invalid SH means lag'
            saw_lag = .true.
        case ('n_channels')
            read(value,*,iostat=ios) results%n_channels
            if (ios /= 0) error stop 'Invalid SH means channel count'
            saw_channels = .true.
        case ('n_records')
            read(value,*,iostat=ios) n_records
            if (ios /= 0) error stop 'Invalid SH means record count'
            saw_records = .true.
        case default
            error stop 'Unexpected SH product means header key'
        end select
    end do
    backspace(unit)
    if (.not. saw_version .or. .not. saw_state .or. .not. saw_lmax .or. .not. saw_channels .or. &
        .not. saw_records .or. &
        (format_version >= 2 .and. .not. saw_lag)) then
        error stop 'SH product means header is incomplete'
    end if
    results%format_version = format_version
    if (n_records < 0_int64 .or. n_records > int(huge(i),int64)) then
        error stop 'SH product means record count is invalid or too large to load'
    end if
    results%n_records = n_records
    allocate(results%records(int(n_records)), stat=alloc_status)
    if (alloc_status /= 0) error stop 'Cannot allocate SH product means records'
    do i = 1, int(n_records)
        read(unit,*,iostat=ios) results%records(i)%kind, results%records(i)%channel1, &
            results%records(i)%channel2, results%records(i)%channel3, &
            results%records(i)%l1, results%records(i)%m1, results%records(i)%l2, &
            results%records(i)%m2, results%records(i)%l3, results%records(i)%m3, &
            results%records(i)%mean_re, results%records(i)%mean_im, results%records(i)%n_samples
        if (ios /= 0) error stop 'Incomplete SH product means data'
    end do
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading SH product means data'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        error stop 'Unexpected data after SH product means records'
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing SH product means input'
end subroutine read_sh_product_means

subroutine write_sh_product_sampling_workflow(job, filename)
    type(sh_product_sampling_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_sh_product_sampling_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call missing_file('Cannot write SH product workflow', trim(filename))
    write(unit,'(A)') '# Spherical-harmonic product sampling workflow'
    write(unit,'(A)') 'workflow_type = sh_product_sampling'
    write(unit,'(A)') 'format_version = 3'
    write(unit,'(A,I0)') 'lmax = ', job%lmax
    write(unit,'(A,I0)') 'lag = ', job%lag
    write(unit,'(A,A)') 'selector_filename = ', trim(job%selector_filename)
    write(unit,'(A,A)') 'lifting_plan_filename = ', trim(job%lifting_plan_filename)
    write(unit,'(A,A)') 'pair_active_filename = ', trim(job%pair_active_filename)
    write(unit,'(A,A)') 'trajectory_path = ', trim(job%trajectory_path)
    write(unit,'(A,A)') 'trajectory_prefix = ', trim(job%trajectory_prefix)
    write(unit,'(A,A)') 'output_path = ', trim(job%output_path)
    write(unit,'(A,A)') 'output_prefix = ', trim(job%output_prefix)
    write(unit,'(A,A)') 'output_policy = ', trim(job%output_policy)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing SH product workflow file'
end subroutine write_sh_product_sampling_workflow

subroutine read_sh_product_sampling_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(sh_product_sampling_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1800) :: value
    integer :: unit, ios, equals_position, version
    logical :: saw_type, saw_version, saw_lag

    job = sh_product_sampling_workflow_type()
    version = -1
    saw_type = .false.
    saw_version = .false.
    saw_lag = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call missing_file('Cannot read SH product workflow', trim(filename))
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading SH product workflow'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) error stop 'Invalid SH product workflow configuration line'
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'sh_product_sampling') error stop 'Unexpected SH workflow type'
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) version
            if (ios /= 0 .or. version < 1 .or. version > 3) &
                error stop 'Unsupported SH product workflow version'
            saw_version = .true.
        case ('lmax')
            read(value,*,iostat=ios) job%lmax
            if (ios /= 0) error stop 'Invalid SH lmax'
        case ('lag')
            read(value,*,iostat=ios) job%lag
            if (ios /= 0) error stop 'Invalid SH sampling lag'
            saw_lag = .true.
        case ('selector_filename')
            job%selector_filename = trim(value)
        case ('lifting_plan_filename')
            job%lifting_plan_filename = trim(value)
        case ('pair_active_filename')
            job%pair_active_filename = trim(value)
        case ('trajectory_path')
            job%trajectory_path = trim(value)
        case ('trajectory_prefix')
            job%trajectory_prefix = trim(value)
        case ('output_path')
            job%output_path = trim(value)
        case ('output_prefix')
            job%output_prefix = trim(value)
        case ('output_policy')
            job%output_policy = trim(value)
        case default
            write(error_unit,'(A,A)') 'Unknown SH product workflow key: ', trim(key)
            error stop 'Invalid SH product workflow configuration'
        end select
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing SH product workflow file'
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete SH product workflow header'
    if (version >= 2 .and. .not. saw_lag) error stop 'SH product workflow is missing lag'
    call validate_sh_product_sampling_workflow(job)
end subroutine read_sh_product_sampling_workflow

subroutine print_sh_product_sampling_workflow(job)
    type(sh_product_sampling_workflow_type), intent(in) :: job
    write(output_unit,'(A,I0)') 'Maximum SH degree L: ', job%lmax
    write(output_unit,'(A,I0)') 'Sampling lag in frames: ', job%lag
    write(output_unit,'(A,A)') 'Replicate selector: ', trim(job%selector_filename)
    write(output_unit,'(A,A)') 'Lifting plan: ', trim(job%lifting_plan_filename)
    write(output_unit,'(A,A)') 'Active-pair matrix: ', trim(job%pair_active_filename)
    write(output_unit,'(A,A,A)') 'Euclidean vectors: ', trim(job%trajectory_path), &
        '/'//trim(job%trajectory_prefix)
    write(output_unit,'(A,A,A)') 'Result file: ', trim(job%output_path), &
        '/'//trim(job%output_prefix)//'_means.txt'
    write(output_unit,'(A,A)') 'Output policy: ', trim(job%output_policy)
end subroutine print_sh_product_sampling_workflow

subroutine missing_file(kind, filename)
    character(len=*), intent(in) :: kind, filename
    write(error_unit,'(3A)') trim(kind), ' missing or unreadable: ', trim(filename)
    error stop 'SH product sampling I/O failed'
end subroutine missing_file

function join_path(directory, filename) result(path)
    character(len=*), intent(in) :: directory, filename
    character(len=2048) :: path
    integer :: n

    n = len_trim(directory)
    if (n == 0) error stop 'Directory path cannot be empty'
    if (directory(n:n) == '/') then
        path = trim(directory)//trim(filename)
    else
        path = trim(directory)//'/'//trim(filename)
    end if
end function join_path

end module sh_product_sampling_workflow_module
