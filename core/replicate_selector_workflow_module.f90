module replicate_selector_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use project_module, only: project_type
    use workflow_module, only: workflow_config_filename, replicate_selector_workflow_config_stem
    implicit none
    private

    public :: replica_selector_type
    public :: init_replica_selector
    public :: replicate_selector_workflow_type
    public :: configure_replicate_selector_workflow
    public :: read_replicate_selector_workflow
    public :: write_replicate_selector_workflow
    public :: validate_replicate_selector_workflow
    public :: run_replicate_selector_workflow
    public :: read_replica_selector_file

    type :: replica_selector_type
        integer :: nrep_total = 0
        integer :: nrep_active = 0
        logical, allocatable :: is_active(:)
    end type replica_selector_type

    type :: replicate_selector_workflow_type
        logical :: initialized = .false.
        integer :: state_index = 1
        character(len=1024) :: selector_file = 'replicate_selector.dat'
        logical, allocatable :: active_replicates(:)
    end type replicate_selector_workflow_type

contains

subroutine init_replica_selector(selector, nrep_total)
    type(replica_selector_type), intent(out) :: selector
    integer, intent(in) :: nrep_total

    selector%nrep_total = nrep_total
    selector%nrep_active = nrep_total
    allocate(selector%is_active(nrep_total))
    selector%is_active = .true.
end subroutine init_replica_selector

subroutine configure_replicate_selector_workflow(job, project)
    type(replicate_selector_workflow_type), intent(out) :: job
    type(project_type), intent(in) :: project
    character(len=1024) :: input, filename
    integer :: choice, ios

    write(*,*)
    write(*,*) '======================================'
    write(*,*) ' Replicate Selector Workflow'
    write(*,*) '======================================'
    write(*,*) '1 = Create new workflow'
    write(*,*) '2 = Load and edit existing workflow'
    write(*,*) '3 = Load existing workflow unchanged'
    write(*,*) '4 = Skip workflow configuration'
    write(*,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid replicate selector workflow option'

    select case (choice)
    case (1)
        job = replicate_selector_workflow_type()
        call edit_replicate_selector_workflow(job, project)
        job%initialized = .true.

    case (2)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename
        if (len_trim(filename) == 0) error stop 'No configuration filename specified'
        call read_replicate_selector_workflow(trim(filename), job)
        call edit_replicate_selector_workflow(job, project)
        job%initialized = .true.

    case (3)
        write(*,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename
        if (len_trim(filename) == 0) error stop 'No configuration filename specified'
        call read_replicate_selector_workflow(trim(filename), job)
        call validate_replicate_selector_workflow(job, project)
        call print_replicate_selector_workflow(job)
        job%initialized = .true.

    case (4)
        job = replicate_selector_workflow_type()
        write(*,*) 'Replicate selector workflow configuration skipped.'

    case default
        error stop 'Invalid replicate selector workflow option'
    end select

    if (.not. job%initialized) return

    write(*,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(*,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') filename
        filename = workflow_config_filename(trim(filename), replicate_selector_workflow_config_stem)
        call write_replicate_selector_workflow(job, trim(filename))
        write(*,'(3A)') 'Configuration saved to "', trim(filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_replicate_selector_workflow

subroutine edit_replicate_selector_workflow(job, project)
    type(replicate_selector_workflow_type), intent(inout) :: job
    type(project_type), intent(in) :: project
    character(len=1024) :: input
    integer :: i, ios, nrep, selection_mode
    logical :: reset_selection

    if (project%nstate < 1 .or. .not. allocated(project%state)) then
        error stop 'Cannot configure replicate selector: project has no states'
    end if

    write(*,*)
    write(*,*) 'Available states and replicate counts:'
    do i = 1, project%nstate
        write(*,'(I4,A,A,A,I0,A)') i, ' = ', trim(project%state(i)%name), &
            ' (', project%state(i)%nrep, ' replicates)'
    end do
    write(*,'(A,I0,A)', advance='no') 'Select state [', job%state_index, ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) job%state_index
        if (ios /= 0) error stop 'Invalid state index'
    end if
    if (job%state_index < 1 .or. job%state_index > project%nstate) then
        error stop 'Selected state index is outside the project'
    end if

    nrep = project%state(job%state_index)%nrep
    if (nrep < 1 .or. .not. allocated(project%state(job%state_index)%rep)) then
        error stop 'Selected project state has no available replicates'
    end if
    reset_selection = .not. allocated(job%active_replicates)
    if (allocated(job%active_replicates)) then
        reset_selection = size(job%active_replicates) /= nrep
    end if
    if (reset_selection) then
        if (allocated(job%active_replicates)) deallocate(job%active_replicates)
        allocate(job%active_replicates(nrep))
        job%active_replicates = .true.
    end if

    write(*,'(A,A,A,I0,A)') 'State "', trim(project%state(job%state_index)%name), &
        '" provides ', nrep, ' replicates.'
    write(*,'(A,A,A)', advance='no') 'Selector output filename [', trim(job%selector_file), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) job%selector_file = trim(input)
    if (len_trim(job%selector_file) == 0) error stop 'Selector output filename cannot be empty'

    write(*,*)
    write(*,*) '1 = Activate all replicates'
    write(*,*) '2 = Select replicates individually'
    write(*,'(A)', advance='no') 'Selection mode [1]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) then
        selection_mode = 1
    else
        read(input,*,iostat=ios) selection_mode
        if (ios /= 0) error stop 'Invalid replicate selection mode'
    end if

    select case (selection_mode)
    case (1)
        job%active_replicates = .true.
    case (2)
        do i = 1, nrep
            if (job%active_replicates(i)) then
                write(*,'(A,I0,A)', advance='no') 'Activate replicate ', i, '? [y]: '
            else
                write(*,'(A,I0,A)', advance='no') 'Activate replicate ', i, '? [n]: '
            end if
            read(*,'(A)') input
            if (len_trim(input) == 0) cycle
            select case (trim(adjustl(input)))
            case ('y', 'Y', 'yes', 'YES', 'Yes')
                job%active_replicates(i) = .true.
            case ('n', 'N', 'no', 'NO', 'No')
                job%active_replicates(i) = .false.
            case default
                error stop 'Invalid replicate activation; enter y or n'
            end select
        end do
    case default
        error stop 'Invalid replicate selection mode'
    end select

    call validate_replicate_selector_workflow(job, project)
end subroutine edit_replicate_selector_workflow

subroutine run_replicate_selector_workflow(job, project, selector)
    type(replicate_selector_workflow_type), intent(in) :: job
    type(project_type), intent(in) :: project
    type(replica_selector_type), intent(out) :: selector
    integer :: output_unit_number, ios

    call validate_replicate_selector_workflow(job, project)
    call init_replica_selector(selector, project%state(job%state_index)%nrep)
    selector%is_active = job%active_replicates
    selector%nrep_active = count(selector%is_active)

    open(newunit=output_unit_number, file=trim(job%selector_file), status='replace', &
        action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write replicate selector file: ', trim(job%selector_file)
        error stop 'Replicate selector output failed'
    end if

    write(output_unit_number,'(A,I0)') 'state_index = ', job%state_index
    write(output_unit_number,'(A,I0)') 'n_replicates = ', selector%nrep_total
    write(output_unit_number,'(A,I0)') 'n_active = ', selector%nrep_active
    write(output_unit_number,'(A)',advance='no') 'active_replicates ='
    write(output_unit_number,'(*(1X,L1))') selector%is_active
    close(output_unit_number, iostat=ios)
    if (ios /= 0) error stop 'Error closing replicate selector output file'

    write(*,'(A,I0,A,I0,A)') 'Initialized selector with ', selector%nrep_active, &
        ' active of ', selector%nrep_total, ' replicates.'
    write(*,'(3A)') 'Replicate selector written to "', trim(job%selector_file), '"'
end subroutine run_replicate_selector_workflow

subroutine read_replica_selector_file(filename, selector, state_index)
    character(len=*), intent(in) :: filename
    type(replica_selector_type), intent(out) :: selector
    integer, intent(out), optional :: state_index
    character(len=2048) :: line, key, value
    integer :: input_unit, ios, equals_position
    integer :: file_state, nrep, nactive
    logical, allocatable :: active(:)
    logical :: saw_state, saw_nrep, saw_active

    file_state = 0
    nrep = 0
    nactive = -1
    saw_state = .false.
    saw_nrep = .false.
    saw_active = .false.

    open(newunit=input_unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot open replicate selector file: ', trim(filename)
        error stop 'Replicate selector input failed'
    end if

    do
        read(input_unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) error stop 'Error reading replicate selector file'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

        equals_position = index(line, '=')
        if (equals_position == 0) error stop 'Invalid line in replicate selector file'
        key = adjustl(trim(line(:equals_position-1)))
        value = adjustl(trim(line(equals_position+1:)))
        select case (trim(key))
        case ('state_index')
            read(value,*,iostat=ios) file_state
            if (ios /= 0) error stop 'Invalid state_index in replicate selector file'
            saw_state = .true.
        case ('n_replicates')
            read(value,*,iostat=ios) nrep
            if (ios /= 0) error stop 'Invalid n_replicates in selector file'
            if (nrep < 1) error stop 'Invalid n_replicates in selector file'
            if (allocated(active)) deallocate(active)
            allocate(active(nrep))
            saw_nrep = .true.
        case ('n_active')
            read(value,*,iostat=ios) nactive
            if (ios /= 0) error stop 'Invalid n_active in replicate selector file'
        case ('active_replicates')
            if (.not. saw_nrep) error stop 'n_replicates must precede active_replicates'
            read(value,*,iostat=ios) active
            if (ios /= 0) error stop 'Invalid active_replicates in selector file'
            saw_active = .true.
        end select
    end do
    close(input_unit)

    if (.not. saw_state .or. .not. saw_nrep .or. .not. saw_active) then
        error stop 'Replicate selector file is missing required values'
    end if

    call init_replica_selector(selector, nrep)
    selector%is_active = active
    selector%nrep_active = count(active)
    if (nactive >= 0 .and. nactive /= selector%nrep_active) then
        error stop 'n_active does not match active_replicates in selector file'
    end if
    if (present(state_index)) state_index = file_state
end subroutine read_replica_selector_file

subroutine read_replicate_selector_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(replicate_selector_workflow_type), intent(out) :: job
    character(len=2048) :: line, key, value
    integer :: input_unit, ios, equals_position, nrep
    logical :: saw_type, saw_state, saw_filename, saw_nrep, saw_active

    nrep = 0
    saw_type = .false.
    saw_state = .false.
    saw_filename = .false.
    saw_nrep = .false.
    saw_active = .false.

    open(newunit=input_unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot open replicate selector workflow: ', trim(filename)
        error stop 'Replicate selector workflow input failed'
    end if

    do
        read(input_unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) error stop 'Error reading replicate selector workflow'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

        equals_position = index(line, '=')
        if (equals_position == 0) error stop 'Invalid line in replicate selector workflow'
        key = adjustl(trim(line(:equals_position-1)))
        value = adjustl(trim(line(equals_position+1:)))
        select case (trim(key))
        case ('workflow_type')
            if (trim(value) /= 'replicate_selector') error stop 'Unexpected workflow_type'
            saw_type = .true.
        case ('state_index')
            read(value,*,iostat=ios) job%state_index
            if (ios /= 0) error stop 'Invalid state_index in workflow'
            saw_state = .true.
        case ('selector_file')
            job%selector_file = trim(value)
            saw_filename = .true.
        case ('n_replicates')
            read(value,*,iostat=ios) nrep
            if (ios /= 0) error stop 'Invalid n_replicates in workflow'
            if (nrep < 1) error stop 'Invalid n_replicates in workflow'
            if (allocated(job%active_replicates)) deallocate(job%active_replicates)
            allocate(job%active_replicates(nrep))
            saw_nrep = .true.
        case ('active_replicates')
            if (.not. saw_nrep) error stop 'n_replicates must precede active_replicates'
            read(value,*,iostat=ios) job%active_replicates
            if (ios /= 0) error stop 'Invalid active_replicates in workflow'
            saw_active = .true.
        end select
    end do
    close(input_unit)

    if (.not. saw_type .or. .not. saw_state .or. .not. saw_filename .or. &
        .not. saw_nrep .or. .not. saw_active) then
        error stop 'Replicate selector workflow is missing required values'
    end if
    job%initialized = .true.
end subroutine read_replicate_selector_workflow

subroutine write_replicate_selector_workflow(job, filename)
    type(replicate_selector_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: output_unit_number, ios

    if (.not. allocated(job%active_replicates)) error stop 'Replicate selector workflow has no selection'
    open(newunit=output_unit_number, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write replicate selector workflow: ', trim(filename)
        error stop 'Replicate selector workflow output failed'
    end if

    write(output_unit_number,'(A)') '# Replicate selector workflow configuration'
    write(output_unit_number,'(A)') 'workflow_type = replicate_selector'
    write(output_unit_number,'(A,I0)') 'state_index = ', job%state_index
    write(output_unit_number,'(A,A)') 'selector_file = ', trim(job%selector_file)
    write(output_unit_number,'(A,I0)') 'n_replicates = ', size(job%active_replicates)
    write(output_unit_number,'(A)',advance='no') 'active_replicates ='
    write(output_unit_number,'(*(1X,L1))') job%active_replicates
    close(output_unit_number, iostat=ios)
    if (ios /= 0) error stop 'Error closing replicate selector workflow'
end subroutine write_replicate_selector_workflow

subroutine validate_replicate_selector_workflow(job, project)
    type(replicate_selector_workflow_type), intent(in) :: job
    type(project_type), intent(in) :: project

    if (project%nstate < 1 .or. .not. allocated(project%state)) then
        error stop 'Replicate selector project has no states'
    end if
    if (size(project%state) < project%nstate) error stop 'Project state array is inconsistent'
    if (.not. job%initialized .and. .not. allocated(job%active_replicates)) then
        error stop 'Replicate selector workflow is not configured'
    end if
    if (job%state_index < 1 .or. job%state_index > project%nstate) then
        error stop 'Replicate selector state is outside the project'
    end if
    if (project%state(job%state_index)%nrep < 1 .or. &
        .not. allocated(project%state(job%state_index)%rep)) then
        error stop 'Selected state has no available replicates'
    end if
    if (size(project%state(job%state_index)%rep) /= project%state(job%state_index)%nrep) then
        error stop 'Project replicate array is inconsistent'
    end if
    if (.not. allocated(job%active_replicates)) then
        error stop 'Replicate selector workflow has no active-replicate mask'
    end if
    if (size(job%active_replicates) /= project%state(job%state_index)%nrep) then
        error stop 'Replicate selector mask size does not match selected state'
    end if
    if (count(job%active_replicates) == 0) then
        error stop 'At least one replicate must be active'
    end if
    if (len_trim(job%selector_file) == 0) error stop 'Replicate selector filename cannot be empty'
end subroutine validate_replicate_selector_workflow

subroutine print_replicate_selector_workflow(job)
    type(replicate_selector_workflow_type), intent(in) :: job
    integer :: i

    write(output_unit,'(A,I0)') 'State index: ', job%state_index
    write(output_unit,'(A,A)') 'Selector file: ', trim(job%selector_file)
    if (allocated(job%active_replicates)) then
        write(output_unit,'(A,I0,A,I0)') 'Active replicates: ', count(job%active_replicates), &
            ' of ', size(job%active_replicates)
        do i = 1, size(job%active_replicates)
            write(output_unit,'(A,I0,A,L1)') '  Replicate ', i, ' active = ', job%active_replicates(i)
        end do
    end if
end subroutine print_replicate_selector_workflow

end module replicate_selector_workflow_module
