module te_plot_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use workflow_module, only: te_plot_workflow_type, workflow_config_filename, &
        te_plot_workflow_config_stem
    use te_io_module, only: te_plot_data_type, read_te_plot_data
    use te_plot_module, only: write_te_summary_table, fill_pdb_occ_and_beta_sidechains
    use output_policy_module, only: check_output_policy, OUTPUT_WRITE, OUTPUT_SKIP, &
        OUTPUT_COMPUTE_ONLY, OUTPUT_ABORT
    use project_module, only: project_type
    use topology_module, only: topology_set_type
    use direction_vectors_analyse_plan_module, only: cv_plan
    use cv_plan_io_module, only: read_cv_plan
    use position_lift_module, only: lifted_trajectory, lift_cv_plan
    use time_window_file_mapping, only: time_window
    use time_window_loader, only: raw_data_block, load_raw_data_block
    use io_pdb_module, only: write_subgroup_pdb
    implicit none
    private

    public :: te_plot_workflow_type
    public :: configure_te_plot_workflow
    public :: read_te_plot_workflow
    public :: write_te_plot_workflow
    public :: validate_te_plot_workflow
    public :: check_te_plot_workflow_inputs
    public :: run_te_plot_workflow

contains

subroutine configure_te_plot_workflow(job)
    type(te_plot_workflow_type), intent(out) :: job
    character(len=1024) :: input, filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' TE Matrix Plot Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)',advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid TE plot workflow option'

    select case(choice)
    case(1)
        job = te_plot_workflow_type()
        call edit_te_plot_workflow(job)
    case(2,3)
        write(output_unit,'(A)',advance='no') 'Configuration filename: '
        read(*,'(A)') filename
        if (len_trim(filename) == 0) error stop 'No TE plot workflow filename specified'
        call read_te_plot_workflow(trim(filename),job)
        if (choice == 2) then
            call edit_te_plot_workflow(job)
        else
            call print_te_plot_workflow(job)
        end if
    case(4)
        job = te_plot_workflow_type()
        write(output_unit,*) 'TE plot workflow configuration skipped.'
        return
    case default
        error stop 'Invalid TE plot workflow option'
    end select

    call validate_te_plot_workflow(job)
    write(output_unit,'(A)',advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case(trim(adjustl(input)))
    case('y','Y','yes','YES','Yes')
        write(output_unit,'(A)',advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        filename = workflow_config_filename(trim(input),te_plot_workflow_config_stem)
        call write_te_plot_workflow(job,trim(filename))
        write(output_unit,'(3A)') 'Configuration saved to "',trim(filename),'"'
    case('n','N','no','NO','No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_te_plot_workflow

subroutine edit_te_plot_workflow(job)
    type(te_plot_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios, i, n_files

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_string('Raw TE matrix filename',job%raw_matrix_filename)
    n_files = job%n_shuffle_files
    write(output_unit,'(A,I0,A)',advance='no') 'Number of shuffle matrix files [',n_files,']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) n_files
        if (ios /= 0 .or. n_files < 1) error stop 'Invalid number of shuffle matrices'
        if (allocated(job%shuffle_filenames)) deallocate(job%shuffle_filenames)
        allocate(job%shuffle_filenames(n_files))
        job%n_shuffle_files = n_files
        do i = 1, n_files
            write(output_unit,'(A,I0,A)',advance='no') 'Shuffle matrix ',i,' filename: '
            read(*,'(A)') job%shuffle_filenames(i)
        end do
    else if (n_files > 0 .and. .not. allocated(job%shuffle_filenames)) then
        error stop 'Shuffle filenames are missing'
    end if
    if (n_files > 0) then
        do i = 1, n_files
            write(output_unit,'(A,I0,2A)',advance='no') 'Shuffle ',i,' [', &
                trim(job%shuffle_filenames(i))//']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) job%shuffle_filenames(i) = trim(input)
        end do
    end if

    call edit_logical('Write summary table',job%write_summary)
    if (job%write_summary) call edit_string('Summary filename',job%summary_filename)
    call edit_logical('Write PDB output',job%write_pdb)
    if (job%write_pdb) then
        call edit_string('CV plan filename',job%cv_plan_filename)
        call edit_integer('State index',job%state_index)
        call edit_integer('Replicate index',job%replicate_index)
        call edit_integer('Frame index',job%frame_index)
        call edit_integer('First CV target to include (offset)',job%cv_target_offset)
        call edit_real('PDB value scale',job%pdb_scale)
        call edit_real('Absolute B-factor threshold for atom selection',job%pdb_threshold)
        call edit_string('PDB output filename',job%pdb_filename)
    end if
    call edit_string('Output policy (ask/replace/new)',job%output_policy)
    call validate_te_plot_workflow(job)
    call print_te_plot_workflow(job)

contains
    subroutine edit_string(prompt,value)
        character(len=*), intent(in) :: prompt
        character(len=*), intent(inout) :: value
        write(output_unit,'(A,A,A,A,A)',advance='no') trim(prompt),' [',trim(value),']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) value = trim(input)
    end subroutine edit_string

    subroutine edit_integer(prompt,value)
        character(len=*), intent(in) :: prompt
        integer, intent(inout) :: value
        write(output_unit,'(A,A,I0,A)',advance='no') trim(prompt),' [',value,']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid integer TE plot setting'
    end subroutine edit_integer

    subroutine edit_real(prompt,value)
        character(len=*), intent(in) :: prompt
        real(8), intent(inout) :: value
        write(output_unit,'(A,A,ES12.4,A)',advance='no') trim(prompt),' [',value,']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid real TE plot setting'
    end subroutine edit_real

    subroutine edit_logical(prompt,value)
        character(len=*), intent(in) :: prompt
        logical, intent(inout) :: value
        write(output_unit,'(A,L1,A)',advance='no') trim(prompt),' [',value,']: '
        read(*,'(A)') input
        if (len_trim(input) == 0) return
        read(input,*,iostat=ios) value
        if (ios /= 0) error stop 'Invalid logical TE plot setting'
    end subroutine edit_logical
end subroutine edit_te_plot_workflow

subroutine validate_te_plot_workflow(job)
    type(te_plot_workflow_type), intent(in) :: job
    integer :: i

    if (len_trim(job%raw_matrix_filename) == 0) error stop 'Raw TE matrix filename cannot be empty'
    if (job%n_shuffle_files < 1) error stop 'At least one shuffle matrix is required'
    if (.not. allocated(job%shuffle_filenames)) error stop 'Shuffle matrix filenames are not allocated'
    if (size(job%shuffle_filenames) /= job%n_shuffle_files) &
        error stop 'Shuffle filename count does not match configuration'
    do i = 1, job%n_shuffle_files
        if (len_trim(job%shuffle_filenames(i)) == 0) error stop 'Shuffle matrix filename cannot be empty'
    end do
    if (.not. job%write_summary .and. .not. job%write_pdb) &
        error stop 'At least one TE plot output must be enabled'
    if (job%write_summary .and. len_trim(job%summary_filename) == 0) &
        error stop 'TE summary filename cannot be empty'
    if (job%write_pdb) then
        if (len_trim(job%cv_plan_filename) == 0) error stop 'CV plan filename cannot be empty'
        if (len_trim(job%pdb_filename) == 0) error stop 'PDB filename cannot be empty'
        if (job%state_index < 1 .or. job%replicate_index < 1 .or. job%frame_index < 1) &
            error stop 'PDB state, replicate and frame indices must be positive'
        if (job%cv_target_offset < 0) error stop 'CV target offset cannot be negative'
        if (job%pdb_scale <= 0.0_8) error stop 'PDB scale must be positive'
        if (job%pdb_threshold < 0.0_8) error stop 'PDB threshold cannot be negative'
    end if
    select case(trim(job%output_policy))
    case('ask','replace','new')
    case default
        error stop 'TE plot output policy must be ask, replace, or new'
    end select
end subroutine validate_te_plot_workflow

subroutine check_te_plot_workflow_inputs(job, cv_plan_will_be_generated)
    type(te_plot_workflow_type), intent(in) :: job
    logical, intent(in), optional :: cv_plan_will_be_generated
    type(te_plot_data_type) :: data
    type(cv_plan) :: plan
    logical :: plan_generated

    call validate_te_plot_workflow(job)
    call read_te_plot_data(trim(job%raw_matrix_filename),job%shuffle_filenames,data)
    plan_generated = .false.
    if (present(cv_plan_will_be_generated)) plan_generated = cv_plan_will_be_generated
    if (job%write_pdb .and. .not. plan_generated) then
        call read_cv_plan(plan,trim(job%cv_plan_filename))
        if (.not. allocated(plan%object)) error stop 'CV plan contains no objects'
        if (data%matrix_size /= 2*size(plan%object)) &
            error stop 'TE matrix size does not match the CV plan dimensions'
    end if
end subroutine check_te_plot_workflow_inputs

subroutine run_te_plot_workflow(job, project, topologies)
    type(te_plot_workflow_type), intent(in) :: job
    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(te_plot_data_type) :: data
    type(cv_plan) :: plan
    type(raw_data_block) :: block
    type(lifted_trajectory) :: lifted
    type(time_window) :: window
    logical, allocatable :: mask(:), group(:)
    real(4), allocatable :: atom_occupancy(:), atom_beta(:)
    real(8), allocatable :: positions(:,:)
    integer :: summary_action, pdb_action, atom_count

    call validate_te_plot_workflow(job)
    call read_te_plot_data(trim(job%raw_matrix_filename),job%shuffle_filenames,data)
    summary_action = OUTPUT_SKIP
    pdb_action = OUTPUT_SKIP
    if (job%write_summary) then
        call check_output_policy(trim(job%summary_filename),trim(job%output_policy),summary_action)
        if (summary_action == OUTPUT_ABORT) error stop 'TE summary output policy aborted workflow'
    end if
    if (job%write_pdb) then
        call check_output_policy(trim(job%pdb_filename),trim(job%output_policy),pdb_action)
        if (pdb_action == OUTPUT_ABORT) error stop 'TE PDB output policy aborted workflow'
    end if
    if (pdb_action == OUTPUT_WRITE .or. pdb_action == OUTPUT_COMPUTE_ONLY) then
        if (data%matrix_size < 1) error stop 'TE matrix is empty'
        if (job%state_index > project%nstate) error stop 'PDB state index is outside the project'
        if (job%replicate_index > project%state(job%state_index)%nrep) &
            error stop 'PDB replicate index is outside the selected state'
        if (job%frame_index > project%state(job%state_index)%rep(job%replicate_index)%n_frames_use_total) &
            error stop 'PDB frame index is outside the selected replicate'

        call read_cv_plan(plan,trim(job%cv_plan_filename))
        window%istate = job%state_index
        window%irep = job%replicate_index
        window%start_time_step = job%frame_index
        window%end_time_step = job%frame_index
        allocate(mask(1))
        mask = .true.
        call load_raw_data_block(project,window,mask,block)
        call lift_cv_plan(plan,block,lifted)

        atom_count = size(block%xyz,2)
        if (size(topologies%state(job%state_index)%topo%atom_id) /= atom_count) &
            error stop 'PDB topology and selected frame have different atom counts'
        if (size(lifted%vec,2) /= data%matrix_size) &
            error stop 'TE matrix size does not match the CV plan dimensions'
        allocate(positions(3,atom_count),atom_occupancy(atom_count),atom_beta(atom_count),group(atom_count))
        positions = block%xyz(:,:,1)
        call fill_pdb_occ_and_beta_sidechains(plan,data%antisymmetric_column_sum,job%cv_target_offset, &
            job%pdb_scale,atom_occupancy,atom_beta)
        group = abs(real(atom_beta,8)) > job%pdb_threshold
    end if

    if (summary_action == OUTPUT_WRITE) call write_te_summary_table(job%summary_filename,data)
    if (pdb_action == OUTPUT_WRITE) then
        call write_subgroup_pdb(trim(job%pdb_filename),topologies%state(job%state_index)%topo, &
            real(positions,4),atom_occupancy,atom_beta,group)
    end if
    if (summary_action == OUTPUT_COMPUTE_ONLY .or. pdb_action == OUTPUT_COMPUTE_ONLY) then
        write(output_unit,*) 'TE plot outputs calculated without writing by output policy.'
    end if
end subroutine run_te_plot_workflow

subroutine write_te_plot_workflow(job,filename)
    type(te_plot_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios, i

    call validate_te_plot_workflow(job)
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) error stop 'Cannot write TE plot workflow configuration'
    write(unit,'(A)') '# TE matrix plotting workflow'
    write(unit,'(A)') 'workflow_type = te_plot'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,A)') 'raw_matrix_filename = ',trim(job%raw_matrix_filename)
    write(unit,'(A,I0)') 'n_shuffle_files = ',job%n_shuffle_files
    do i = 1, job%n_shuffle_files
        write(unit,'(A,A)') 'shuffle_filename = ',trim(job%shuffle_filenames(i))
    end do
    write(unit,'(A,L1)') 'write_summary = ',job%write_summary
    write(unit,'(A,A)') 'summary_filename = ',trim(job%summary_filename)
    write(unit,'(A,L1)') 'write_pdb = ',job%write_pdb
    write(unit,'(A,A)') 'cv_plan_filename = ',trim(job%cv_plan_filename)
    write(unit,'(A,I0)') 'state_index = ',job%state_index
    write(unit,'(A,I0)') 'replicate_index = ',job%replicate_index
    write(unit,'(A,I0)') 'frame_index = ',job%frame_index
    write(unit,'(A,I0)') 'cv_target_offset = ',job%cv_target_offset
    write(unit,'(A,ES24.16E3)') 'pdb_scale = ',job%pdb_scale
    write(unit,'(A,ES24.16E3)') 'pdb_threshold = ',job%pdb_threshold
    write(unit,'(A,A)') 'pdb_filename = ',trim(job%pdb_filename)
    write(unit,'(A,A)') 'output_policy = ',trim(job%output_policy)
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing TE plot workflow configuration'
end subroutine write_te_plot_workflow

subroutine read_te_plot_workflow(filename,job)
    character(len=*), intent(in) :: filename
    type(te_plot_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1800) :: value
    integer :: unit, ios, equals_position, version, shuffle_index
    logical :: saw_type, saw_version, saw_shuffle_count

    job = te_plot_workflow_type()
    version = -1
    shuffle_index = 0
    saw_type = .false.
    saw_version = .false.
    saw_shuffle_count = .false.
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read TE plot workflow: ',trim(filename)
        error stop 'TE plot workflow input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading TE plot workflow configuration'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line,'=')
        if (equals_position == 0) error stop 'Invalid TE plot workflow configuration line'
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case(trim(key))
        case('workflow_type')
            if (trim(value) /= 'te_plot') error stop 'Unexpected TE plot workflow type'
            saw_type = .true.
        case('format_version')
            read(value,*,iostat=ios) version
            if (ios /= 0 .or. version /= 1) error stop 'Unsupported TE plot workflow version'
            saw_version = .true.
        case('raw_matrix_filename')
            job%raw_matrix_filename = trim(value)
        case('n_shuffle_files')
            read(value,*,iostat=ios) job%n_shuffle_files
            if (ios /= 0 .or. job%n_shuffle_files < 1) error stop 'Invalid shuffle filename count'
            if (allocated(job%shuffle_filenames)) deallocate(job%shuffle_filenames)
            allocate(job%shuffle_filenames(job%n_shuffle_files))
            shuffle_index = 0
            saw_shuffle_count = .true.
        case('shuffle_filename')
            if (.not. saw_shuffle_count) error stop 'Shuffle count must precede shuffle filenames'
            shuffle_index = shuffle_index + 1
            if (shuffle_index > job%n_shuffle_files) error stop 'Too many shuffle filenames'
            job%shuffle_filenames(shuffle_index) = trim(value)
        case('write_summary')
            read(value,*,iostat=ios) job%write_summary
            if (ios /= 0) error stop 'Invalid write_summary setting'
        case('summary_filename')
            job%summary_filename = trim(value)
        case('write_pdb')
            read(value,*,iostat=ios) job%write_pdb
            if (ios /= 0) error stop 'Invalid write_pdb setting'
        case('cv_plan_filename')
            job%cv_plan_filename = trim(value)
        case('state_index')
            read(value,*,iostat=ios) job%state_index
            if (ios /= 0) error stop 'Invalid TE plot state index'
        case('replicate_index')
            read(value,*,iostat=ios) job%replicate_index
            if (ios /= 0) error stop 'Invalid TE plot replicate index'
        case('frame_index')
            read(value,*,iostat=ios) job%frame_index
            if (ios /= 0) error stop 'Invalid TE plot frame index'
        case('cv_target_offset')
            read(value,*,iostat=ios) job%cv_target_offset
            if (ios /= 0) error stop 'Invalid TE plot CV target offset'
        case('pdb_scale')
            read(value,*,iostat=ios) job%pdb_scale
            if (ios /= 0) error stop 'Invalid TE plot PDB scale'
        case('pdb_threshold')
            read(value,*,iostat=ios) job%pdb_threshold
            if (ios /= 0) error stop 'Invalid TE plot PDB threshold'
        case('pdb_filename')
            job%pdb_filename = trim(value)
        case('output_policy')
            job%output_policy = trim(value)
        case default
            error stop 'Unknown TE plot workflow configuration key'
        end select
    end do
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing TE plot workflow configuration'
    if (.not. saw_type .or. .not. saw_version .or. .not. saw_shuffle_count) &
        error stop 'Incomplete TE plot workflow configuration'
    if (shuffle_index /= job%n_shuffle_files) error stop 'Missing shuffle filenames'
    call validate_te_plot_workflow(job)
end subroutine read_te_plot_workflow

subroutine print_te_plot_workflow(job)
    type(te_plot_workflow_type), intent(in) :: job
    integer :: i

    write(output_unit,'(A,A)') 'Raw matrix: ',trim(job%raw_matrix_filename)
    write(output_unit,'(A,I0)') 'Shuffle matrices: ',job%n_shuffle_files
    do i = 1, job%n_shuffle_files
        write(output_unit,'(A,A)') '  ',trim(job%shuffle_filenames(i))
    end do
    write(output_unit,'(A,L1,A,A)') 'Write summary: ',job%write_summary, &
        ' -> ',trim(job%summary_filename)
    write(output_unit,'(A,L1)') 'Write PDB: ',job%write_pdb
    if (job%write_pdb) then
        write(output_unit,'(A,A)') 'CV plan: ',trim(job%cv_plan_filename)
        write(output_unit,'(A,3(I0,1X))') 'State/replicate/frame: ', &
            job%state_index,job%replicate_index,job%frame_index
        write(output_unit,'(A,A)') 'PDB file: ',trim(job%pdb_filename)
    end if
end subroutine print_te_plot_workflow

end module te_plot_workflow_module
