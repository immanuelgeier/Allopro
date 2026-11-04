module odf_analysis_workflow_module
    use iso_fortran_env, only: output_unit, error_unit, real64
    use workflow_module, only: odf_analysis_workflow_type, workflow_config_filename, &
        odf_analysis_workflow_config_stem
    use direction_vectors_analyse_plan_module, only: cv_plan
    use cv_plan_io_module, only: read_cv_plan
    use lifting_plan_module, only: lifting_plan_type
    use odf_sampling_workflow_module, only: odf_sampling_workflow_type, read_odf_sampling_workflow
    use healpix_state_space_exploration, only: odf_lib_type, odf_level_analysis_type, &
        init_odf_lib, read_odf_lib, compute_odf_entropy_levels, write_entropy, &
        write_odf_analysis_text, downsample_hist
    implicit none
    private

    public :: odf_analysis_workflow_type
    public :: configure_odf_analysis_workflow
    public :: read_odf_analysis_workflow
    public :: write_odf_analysis_workflow
    public :: validate_odf_analysis_workflow
    public :: check_odf_analysis_workflow_inputs
    public :: run_odf_analysis_workflow

    type :: odf_source_type
        type(odf_sampling_workflow_type) :: workflow
        type(cv_plan) :: cv
        type(odf_lib_type) :: library
        character(len=1024) :: cv_plan_filename
        character(len=1024) :: odf_path
    end type odf_source_type

contains

subroutine configure_odf_analysis_workflow(job)
    type(odf_analysis_workflow_type), intent(out) :: job
    character(len=1024) :: input, config_filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' ODF Analysis Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid ODF analysis workflow option'

    select case (choice)
    case (1)
        job = odf_analysis_workflow_type()
        call edit_odf_analysis_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') config_filename
        if (len_trim(config_filename) == 0) error stop 'No configuration filename specified'
        call read_odf_analysis_workflow(trim(config_filename), job)
        if (choice == 2) call edit_odf_analysis_workflow(job)
        if (choice == 3) call print_odf_analysis_workflow(job)
    case (4)
        job = odf_analysis_workflow_type()
        write(output_unit,*) 'ODF analysis workflow configuration skipped.'
        return
    case default
        error stop 'Invalid ODF analysis workflow option'
    end select

    call validate_odf_analysis_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        config_filename = workflow_config_filename(trim(input), odf_analysis_workflow_config_stem)
        call write_odf_analysis_workflow(job, trim(config_filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(config_filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_odf_analysis_workflow

subroutine edit_odf_analysis_workflow(job)
    type(odf_analysis_workflow_type), intent(inout) :: job
    type(odf_source_type) :: source_a, source_b
    character(len=1024) :: input
    integer :: ios

    write(output_unit,*)
    write(output_unit,*) 'ODF library and model parameters are imported from ODF sampling workflows.'
    call edit_string('ODF sampling workflow A configuration', job%sampling_workflow_a)
    call load_odf_source(trim(job%sampling_workflow_a), source_a)
    call print_source_summary('A', source_a)

    call edit_logical('Write entropy analysis', job%analyze_entropy)
    call edit_logical('Write maxima/basin analysis', job%analyze_details)
    call edit_logical('Compare two ODFs', job%compare_odfs)
    call edit_string('Entropy output filename', job%entropy_filename)
    call edit_string('ODF detail output filename', job%analysis_filename)
    call edit_string('ODF comparison output filename', job%comparison_filename)

    if (job%compare_odfs) then
        call edit_string('ODF sampling workflow B configuration', job%sampling_workflow_b)
        call load_odf_source(trim(job%sampling_workflow_b), source_b)
        call print_source_summary('B', source_b)
        call print_available_lags(source_a%workflow, 'A')
        call print_available_lags(source_b%workflow, 'B')
        call print_cv_objects(source_a%cv, 'A')
        call print_cv_objects(source_b%cv, 'B')

        call edit_integer('Lag in ODF library A', job%lag_a)
        call edit_integer('CV object index in library A', job%cv_index_a)
        call edit_component('CV component in library A (1 or 2)', job%cv_component_a)
        call edit_integer('Lag in ODF library B', job%lag_b)
        call edit_integer('CV object index in library B', job%cv_index_b)
        call edit_component('CV component in library B (1 or 2)', job%cv_component_b)
    end if
    call validate_odf_analysis_workflow(job)
    call print_odf_analysis_workflow(job)

contains
    subroutine edit_string(prompt, value)
        character(len=*), intent(in) :: prompt
        character(len=*), intent(inout) :: value
        write(output_unit,'(A,A,A)',advance='no') trim(prompt), ' [', trim(value), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) value = trim(input)
    end subroutine edit_string

    subroutine edit_logical(prompt, value)
        character(len=*), intent(in) :: prompt
        logical, intent(inout) :: value
        character(len=8) :: current
        current = 'false'
        if (value) current = 'true'
        write(output_unit,'(A,A,A)',advance='no') trim(prompt), ' [', trim(current), ']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        select case (trim(adjustl(input)))
        case ('y','Y','yes','YES','true','.true.')
            value = .true.
        case ('n','N','no','NO','false','.false.')
            value = .false.
        case default
            error stop 'Enter yes or no for the analysis selection'
        end select
    end subroutine edit_logical

    subroutine edit_integer(prompt, value)
        character(len=*), intent(in) :: prompt
        integer, intent(inout) :: value
        write(output_unit,'(A,I0,A)',advance='no') trim(prompt), value, ']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid integer selection'
    end subroutine edit_integer

    subroutine edit_component(prompt, value)
        character(len=*), intent(in) :: prompt
        integer, intent(inout) :: value
        write(output_unit,'(A,I0,A)',advance='no') trim(prompt), value, ']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0 .or. value < 1 .or. value > 2) error stop 'CV component must be 1 or 2'
    end subroutine edit_component
end subroutine edit_odf_analysis_workflow

subroutine validate_odf_analysis_workflow(job)
    type(odf_analysis_workflow_type), intent(in) :: job

    if (.not. job%analyze_entropy .and. .not. job%analyze_details .and. .not. job%compare_odfs) then
        error stop 'Enable at least one ODF analysis or comparison'
    end if
    if (len_trim(job%sampling_workflow_a) == 0) error stop 'ODF sampling workflow A is required'
    if (job%compare_odfs) then
        if (len_trim(job%sampling_workflow_b) == 0) error stop 'ODF sampling workflow B is required for comparison'
        if (job%lag_a < 0 .or. job%lag_b < 0) error stop 'ODF comparison lags must be nonnegative'
        if (job%cv_index_a < 1 .or. job%cv_index_b < 1) error stop 'CV object indices must be positive'
        if (job%cv_component_a < 1 .or. job%cv_component_a > 2 .or. &
            job%cv_component_b < 1 .or. job%cv_component_b > 2) error stop 'CV component indices must be 1 or 2'
        if (len_trim(job%comparison_filename) == 0) error stop 'Comparison output filename cannot be empty'
    end if
    if (job%analyze_entropy .and. len_trim(job%entropy_filename) == 0) &
        error stop 'Entropy output filename cannot be empty'
    if (job%analyze_details .and. len_trim(job%analysis_filename) == 0) &
        error stop 'ODF detail output filename cannot be empty'
    if (job%analyze_entropy .and. job%analyze_details) then
        if (trim(job%entropy_filename) == trim(job%analysis_filename)) &
            error stop 'Entropy and detail output filenames must be different'
    end if
    if (job%compare_odfs) then
        if ((job%analyze_entropy .and. trim(job%comparison_filename) == trim(job%entropy_filename)) .or. &
            (job%analyze_details .and. trim(job%comparison_filename) == trim(job%analysis_filename))) then
            error stop 'ODF output filenames must be different'
        end if
    end if
end subroutine validate_odf_analysis_workflow

subroutine check_odf_analysis_workflow_inputs(job, skip_a_cv, skip_a_odfs, skip_b_cv, skip_b_odfs)
    type(odf_analysis_workflow_type), intent(in) :: job
    logical, intent(in), optional :: skip_a_cv, skip_a_odfs, skip_b_cv, skip_b_odfs
    type(odf_sampling_workflow_type) :: source_a, source_b
    character(len=1024) :: cv_a, path_a, cv_b, path_b
    logical :: skip

    call validate_odf_analysis_workflow(job)
    call read_odf_sampling_workflow(trim(job%sampling_workflow_a), source_a)
    call resolve_source_paths(job%sampling_workflow_a, source_a, cv_a, path_a)
    if (job%compare_odfs .and. job%lag_a > source_a%model_params%max_lag) &
        error stop 'Comparison lag A exceeds workflow A maximum lag'
    skip = .false.
    if (present(skip_a_cv)) skip = skip_a_cv
    if (.not. skip) call check_cv_file(cv_a)
    skip = .false.
    if (present(skip_a_odfs)) skip = skip_a_odfs
    if (job%analyze_entropy .or. job%analyze_details) then
        if (.not. skip) call check_source_odfs(source_a, path_a, 0, source_a%model_params%max_lag)
    else if (job%compare_odfs .and. .not. skip) then
        call check_source_odfs(source_a, path_a, 0, source_a%model_params%max_lag)
    end if

    if (job%compare_odfs) then
        call read_odf_sampling_workflow(trim(job%sampling_workflow_b), source_b)
        call resolve_source_paths(job%sampling_workflow_b, source_b, cv_b, path_b)
        if (job%lag_b > source_b%model_params%max_lag) &
            error stop 'Comparison lag B exceeds workflow B maximum lag'
        skip = .false.
        if (present(skip_b_cv)) skip = skip_b_cv
        if (.not. skip) call check_cv_file(cv_b)
        skip = .false.
        if (present(skip_b_odfs)) skip = skip_b_odfs
        if (.not. skip) call check_source_odfs(source_b, path_b, 0, source_b%model_params%max_lag)
    end if
end subroutine check_odf_analysis_workflow_inputs

subroutine run_odf_analysis_workflow(job)
    type(odf_analysis_workflow_type), intent(in) :: job
    type(odf_source_type) :: source_a, source_b
    integer :: object_a, object_b, common_level

    call validate_odf_analysis_workflow(job)
    call check_odf_analysis_workflow_inputs(job)
    call load_odf_source(trim(job%sampling_workflow_a), source_a)

    if (job%analyze_entropy .or. job%analyze_details .or. job%compare_odfs) then
        call compute_odf_entropy_levels(source_a%library)
    end if
    if (job%analyze_entropy) call write_entropy(source_a%library, trim(job%entropy_filename))
    if (job%analyze_details) call write_odf_analysis_text(source_a%library, trim(job%analysis_filename))

    if (job%compare_odfs) then
        call load_odf_source(trim(job%sampling_workflow_b), source_b)
        call validate_comparison_objects(job, source_a, source_b, object_a, object_b)
        call compute_odf_entropy_levels(source_b%library)
        common_level = min(source_a%workflow%model_params%nlevel, source_b%workflow%model_params%nlevel)
        call write_comparison_report(job, source_a, source_b, object_a, object_b, common_level)
    end if
end subroutine run_odf_analysis_workflow

subroutine load_odf_source(config_filename, source)
    character(len=*), intent(in) :: config_filename
    type(odf_source_type), intent(out) :: source
    type(lifting_plan_type) :: plan
    integer :: nobj

    call read_odf_sampling_workflow(trim(config_filename), source%workflow)
    call resolve_source_paths(config_filename, source%workflow, &
        source%cv_plan_filename, source%odf_path)
    call read_cv_plan(source%cv, trim(source%cv_plan_filename))
    if (.not. allocated(source%cv%object)) error stop 'CV plan has no object definitions'
    if (size(source%cv%object) < 1) error stop 'CV plan contains no CV objects'
    nobj = 2 * size(source%cv%object)
    allocate(plan%istates(1))
    plan%istates(1) = source%workflow%model_params%state
    call init_odf_lib(source%cv, source%library, plan, source%workflow%model_params%max_lag, &
        source%workflow%model_params%basis_index, source%workflow%model_params%nlevel, nobj)
    call read_odf_lib(source%library, trim(source%workflow%io_params%odf_prefix), trim(source%odf_path))
end subroutine load_odf_source

subroutine resolve_source_paths(config_filename, source, cv_plan_filename, odf_path)
    character(len=*), intent(in) :: config_filename
    type(odf_sampling_workflow_type), intent(in) :: source
    character(len=*), intent(out) :: cv_plan_filename, odf_path

    cv_plan_filename = resolve_relative_to_config(config_filename, source%io_params%cv_plan_filename)
    odf_path = resolve_relative_to_config(config_filename, source%io_params%odf_path)
    if (len_trim(odf_path) == 0) error stop 'ODF library path cannot be empty'
    if (odf_path(len_trim(odf_path):len_trim(odf_path)) /= '/') then
        odf_path = trim(odf_path)//'/'
    end if
end subroutine resolve_source_paths

function resolve_relative_to_config(config_filename, configured_path) result(resolved_path)
    character(len=*), intent(in) :: config_filename, configured_path
    character(len=1024) :: resolved_path, candidate
    integer :: slash
    logical :: exists

    if (len_trim(configured_path) == 0) error stop 'Empty path in ODF sampling workflow'
    if (configured_path(1:1) == '/') then
        resolved_path = trim(configured_path)
        return
    end if
    slash = scan(trim(config_filename), '/', back=.true.)
    if (slash == 0) then
        resolved_path = trim(configured_path)
    else
        candidate = config_filename(:slash)//trim(configured_path)
        inquire(file=trim(candidate), exist=exists)
        if (exists) then
            resolved_path = candidate
        else
            inquire(file=trim(configured_path), exist=exists)
            if (exists) then
                resolved_path = trim(configured_path)
            else
                resolved_path = candidate
            end if
        end if
    end if
end function resolve_relative_to_config

subroutine check_cv_file(filename)
    character(len=*), intent(in) :: filename
    logical :: exists
    inquire(file=trim(filename), exist=exists)
    if (.not. exists) then
        write(error_unit,'(A,A)') 'ERROR: CV plan file not found: ', trim(filename)
        error stop 'ODF analysis input validation failed'
    end if
end subroutine check_cv_file

subroutine check_source_odfs(source, path, first_lag, last_lag)
    type(odf_sampling_workflow_type), intent(in) :: source
    character(len=*), intent(in) :: path
    integer, intent(in) :: first_lag, last_lag
    character(len=4096) :: filename
    integer :: lag, unit, ios, close_ios

    do lag = first_lag, last_lag
        write(filename,'(A,A,"state_",I4.4,"_lag_",I3.3,".txt")') trim(path), &
            trim(source%io_params%odf_prefix), source%model_params%state, lag
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'ERROR: ODF file not found or unreadable: ', trim(filename)
            error stop 'ODF analysis input validation failed'
        end if
        close(unit, iostat=close_ios)
        if (close_ios /= 0) error stop 'Error closing ODF file during validation'
    end do
end subroutine check_source_odfs

subroutine validate_comparison_objects(job, source_a, source_b, object_a, object_b)
    type(odf_analysis_workflow_type), intent(in) :: job
    type(odf_source_type), intent(in) :: source_a, source_b
    integer, intent(out) :: object_a, object_b
    integer :: nobj_a, nobj_b

    if (job%lag_a > source_a%workflow%model_params%max_lag .or. &
        job%lag_b > source_b%workflow%model_params%max_lag) error stop 'ODF comparison lag is out of range'
    if (job%cv_index_a > size(source_a%cv%object) .or. job%cv_index_b > size(source_b%cv%object)) then
        error stop 'ODF comparison CV object index is out of range'
    end if
    object_a = 2 * (job%cv_index_a - 1) + job%cv_component_a
    object_b = 2 * (job%cv_index_b - 1) + job%cv_component_b
    nobj_a = size(source_a%library%odf(1)%hist, 2)
    nobj_b = size(source_b%library%odf(1)%hist, 2)
    if (object_a > nobj_a .or. object_b > nobj_b) error stop 'ODF comparison histogram column is out of range'
end subroutine validate_comparison_objects

subroutine write_comparison_report(job, source_a, source_b, object_a, object_b, common_level)
    type(odf_analysis_workflow_type), intent(in) :: job
    type(odf_source_type), intent(in) :: source_a, source_b
    integer, intent(in) :: object_a, object_b, common_level
    integer :: unit, ios, level
    integer, allocatable :: hist_a(:), hist_b(:)
    real(real64) :: squared_difference, entropy_a, entropy_b

    open(newunit=unit, file=trim(job%comparison_filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write ODF comparison report: ', trim(job%comparison_filename)
        error stop 'ODF comparison output failed'
    end if
    write(unit,'(A)') '# Selected ODF comparison across ODF sampling workflow libraries'
    write(unit,'(A,A)') '# Workflow A: ', trim(job%sampling_workflow_a)
    write(unit,'(A,A)') '# ODF A: ', trim(source_a%library%odf(job%lag_a+1)%odf_name)
    write(unit,'(A,A)') '# CV A: ', trim(source_a%cv%object(job%cv_index_a)%name(job%cv_component_a))
    write(unit,'(A,A)') '# Workflow B: ', trim(job%sampling_workflow_b)
    write(unit,'(A,A)') '# ODF B: ', trim(source_b%library%odf(job%lag_b+1)%odf_name)
    write(unit,'(A,A)') '# CV B: ', trim(source_b%cv%object(job%cv_index_b)%name(job%cv_component_b))
    write(unit,'(A,I0)') '# Comparison HEALPix levels: ', common_level

    do level = 1, common_level
        call downsample_hist(source_a%library%odf(job%lag_a+1)%hist(:,object_a), &
            source_a%workflow%model_params%nlevel, level, hist_a)
        call downsample_hist(source_b%library%odf(job%lag_b+1)%hist(:,object_b), &
            source_b%workflow%model_params%nlevel, level, hist_b)
        squared_difference = sum((real(hist_a, real64) - real(hist_b, real64))**2)
        deallocate(hist_a, hist_b)
        entropy_a = source_a%library%odf(job%lag_a+1)%analysis(object_a,level)%entropy
        entropy_b = source_b%library%odf(job%lag_b+1)%analysis(object_b,level)%entropy
        write(unit,'(/,A,I0)') 'LEVEL ', level
        write(unit,'(A,ES24.16E3)') '  PIXEL_COUNT_SQUARED_DIFFERENCE ', squared_difference
        write(unit,'(A,3(1X,ES24.16E3))') '  ENTROPY A / B / A-B ', entropy_a, entropy_b, entropy_a-entropy_b
        call write_level_analysis(unit, source_a%library%odf(job%lag_a+1)%analysis(object_a,level), 'A')
        call write_level_analysis(unit, source_b%library%odf(job%lag_b+1)%analysis(object_b,level), 'B')
    end do
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing ODF comparison report'
end subroutine write_comparison_report

subroutine write_level_analysis(unit, analysis, label)
    integer, intent(in) :: unit
    type(odf_level_analysis_type), intent(in) :: analysis
    character(len=*), intent(in) :: label
    integer :: i

    write(unit,'(A,A,A,I0)',advance='no') '  ', label, ' MAXIMA ', analysis%nmax
    if (analysis%nmax > 0 .and. allocated(analysis%maxima_pix)) then
        write(unit,'(A)',advance='no') ' PIX:'
        do i = 1, analysis%nmax
            write(unit,'(1X,I0)',advance='no') analysis%maxima_pix(i)
        end do
        write(unit,'(A)',advance='no') ' VALUES:'
        do i = 1, analysis%nmax
            write(unit,'(1X,I0)',advance='no') analysis%maxima_value(i)
        end do
    end if
    write(unit,*)
    write(unit,'(A,A,A,I0)',advance='no') '  ', label, ' BASINS ', analysis%nbasin
    if (analysis%nbasin > 0 .and. allocated(analysis%basin_pix)) then
        write(unit,'(A)',advance='no') ' PIX:'
        do i = 1, analysis%nbasin
            write(unit,'(1X,I0)',advance='no') analysis%basin_pix(i)
        end do
        write(unit,'(A)',advance='no') ' COUNTS:'
        do i = 1, analysis%nbasin
            write(unit,'(1X,I0)',advance='no') analysis%basin_counts(i)
        end do
    end if
    write(unit,*)
end subroutine write_level_analysis

subroutine print_source_summary(label, source)
    character(len=*), intent(in) :: label
    type(odf_source_type), intent(in) :: source
    write(output_unit,'(A,A,A,I0)') 'Library ', label, ' state: ', source%workflow%model_params%state
    write(output_unit,'(A,A,A,I0)') 'Library ', label, ' maximum lag: ', source%workflow%model_params%max_lag
    write(output_unit,'(A,A,A,I0)') 'Library ', label, ' HEALPix levels: ', source%workflow%model_params%nlevel
    write(output_unit,'(A,A,A,A)') 'Library ', label, ' CV plan: ', trim(source%cv_plan_filename)
    write(output_unit,'(A,A,A,A)') 'Library ', label, ' ODF path: ', trim(source%odf_path)
end subroutine print_source_summary

subroutine print_available_lags(source, label)
    type(odf_sampling_workflow_type), intent(in) :: source
    character(len=*), intent(in) :: label
    write(output_unit,'(A,A,A,I0,A,I0)') 'Library ', label, ' lag range: 0-', &
        source%model_params%max_lag, ', state ', source%model_params%state
end subroutine print_available_lags

subroutine print_cv_objects(cv, label)
    type(cv_plan), intent(in) :: cv
    character(len=*), intent(in) :: label
    integer :: i
    write(output_unit,'(A,A)') 'CV objects in library ', label
    do i = 1, size(cv%object)
        write(output_unit,'(I5,2A)') i, ' = ', trim(cv%object(i)%name(1))//' / '//trim(cv%object(i)%name(2))
    end do
end subroutine print_cv_objects

subroutine write_odf_analysis_workflow(job, filename)
    type(odf_analysis_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    call validate_odf_analysis_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write ODF analysis workflow: ', trim(filename)
        error stop 'ODF analysis workflow output failed'
    end if
    write(unit,'(A)') '# ODF analysis workflow'
    write(unit,'(A)') 'workflow_type = odf_analysis'
    write(unit,'(A)') 'format_version = 2'
    write(unit,'(A,L1)') 'analyze_entropy = ', job%analyze_entropy
    write(unit,'(A,L1)') 'analyze_details = ', job%analyze_details
    write(unit,'(A,L1)') 'compare_odfs = ', job%compare_odfs
    write(unit,'(A,A)') 'sampling_workflow_a = ', trim(job%sampling_workflow_a)
    write(unit,'(A,A)') 'sampling_workflow_b = ', trim(job%sampling_workflow_b)
    write(unit,'(A,I0)') 'lag_a = ', job%lag_a
    write(unit,'(A,I0)') 'lag_b = ', job%lag_b
    write(unit,'(A,I0)') 'cv_index_a = ', job%cv_index_a
    write(unit,'(A,I0)') 'cv_component_a = ', job%cv_component_a
    write(unit,'(A,I0)') 'cv_index_b = ', job%cv_index_b
    write(unit,'(A,I0)') 'cv_component_b = ', job%cv_component_b
    write(unit,'(A,A)') 'entropy_filename = ', trim(job%entropy_filename)
    write(unit,'(A,A)') 'analysis_filename = ', trim(job%analysis_filename)
    write(unit,'(A,A)') 'comparison_filename = ', trim(job%comparison_filename)
    close(unit, iostat=ios)
    if (ios /= 0) error stop 'Error closing ODF analysis workflow file'
end subroutine write_odf_analysis_workflow

subroutine read_odf_analysis_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(odf_analysis_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1900) :: value
    integer :: unit, ios, equals_position, format_version
    logical :: saw_type, saw_version

    job = odf_analysis_workflow_type()
    format_version = -1
    saw_type = .false.
    saw_version = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read ODF analysis workflow: ', trim(filename)
        error stop 'ODF analysis workflow input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) call config_read_error(unit, 'Error reading ODF analysis workflow')
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line, '=')
        if (equals_position == 0) call config_read_error(unit, 'Invalid workflow configuration line')
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'odf_analysis') call config_read_error(unit, 'Unexpected workflow_type')
            saw_type = .true.
        case ('format_version')
            read(value,*,iostat=ios) format_version
            if (ios /= 0 .or. format_version /= 2) call config_read_error(unit, 'Unsupported format_version')
            saw_version = .true.
        case ('analyze_entropy')
            read(value,*,iostat=ios) job%analyze_entropy
            if (ios /= 0) call config_read_error(unit, 'Invalid analyze_entropy')
        case ('analyze_details')
            read(value,*,iostat=ios) job%analyze_details
            if (ios /= 0) call config_read_error(unit, 'Invalid analyze_details')
        case ('compare_odfs')
            read(value,*,iostat=ios) job%compare_odfs
            if (ios /= 0) call config_read_error(unit, 'Invalid compare_odfs')
        case ('sampling_workflow_a')
            job%sampling_workflow_a = trim(value)
        case ('sampling_workflow_b')
            job%sampling_workflow_b = trim(value)
        case ('lag_a')
            read(value,*,iostat=ios) job%lag_a
            if (ios /= 0) call config_read_error(unit, 'Invalid lag_a')
        case ('lag_b')
            read(value,*,iostat=ios) job%lag_b
            if (ios /= 0) call config_read_error(unit, 'Invalid lag_b')
        case ('cv_index_a')
            read(value,*,iostat=ios) job%cv_index_a
            if (ios /= 0) call config_read_error(unit, 'Invalid cv_index_a')
        case ('cv_component_a')
            read(value,*,iostat=ios) job%cv_component_a
            if (ios /= 0) call config_read_error(unit, 'Invalid cv_component_a')
        case ('cv_index_b')
            read(value,*,iostat=ios) job%cv_index_b
            if (ios /= 0) call config_read_error(unit, 'Invalid cv_index_b')
        case ('cv_component_b')
            read(value,*,iostat=ios) job%cv_component_b
            if (ios /= 0) call config_read_error(unit, 'Invalid cv_component_b')
        case ('entropy_filename')
            job%entropy_filename = trim(value)
        case ('analysis_filename')
            job%analysis_filename = trim(value)
        case ('comparison_filename')
            job%comparison_filename = trim(value)
        case default
            call config_read_error(unit, 'Unknown workflow parameter: '//trim(key))
        end select
    end do
    close(unit)
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete ODF analysis workflow file'
    call validate_odf_analysis_workflow(job)
end subroutine read_odf_analysis_workflow

subroutine print_odf_analysis_workflow(job)
    type(odf_analysis_workflow_type), intent(in) :: job
    write(output_unit,'(A,L1)') 'Entropy analysis: ', job%analyze_entropy
    write(output_unit,'(A,L1)') 'Detailed analysis: ', job%analyze_details
    write(output_unit,'(A,L1)') 'Compare ODFs: ', job%compare_odfs
    write(output_unit,'(A,A)') 'ODF sampling workflow A: ', trim(job%sampling_workflow_a)
    if (job%compare_odfs) then
        write(output_unit,'(A,A)') 'ODF sampling workflow B: ', trim(job%sampling_workflow_b)
        write(output_unit,'(A,2(1X,I0))') 'Selected lags A/B:', job%lag_a, job%lag_b
        write(output_unit,'(A,4(1X,I0))') 'CV object/component A/B:', &
            job%cv_index_a, job%cv_component_a, job%cv_index_b, job%cv_component_b
    end if
end subroutine print_odf_analysis_workflow

subroutine config_read_error(unit, message)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    close(unit)
    write(error_unit,'(A)') trim(message)
    error stop 'ODF analysis workflow configuration error'
end subroutine config_read_error

end module odf_analysis_workflow_module
