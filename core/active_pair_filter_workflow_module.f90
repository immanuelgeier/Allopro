module active_pair_filter_workflow_module
    use iso_fortran_env, only: output_unit, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use workflow_module, only: active_pair_filter_workflow_type, workflow_config_filename, &
        active_pair_filter_workflow_config_stem
    use te_io_module, only: read_square_matrix_auto
    use direction_vectors_analyse_plan_module, only: cv_plan
    use cv_plan_io_module, only: read_cv_plan
    use output_policy_module, only: check_output_policy, OUTPUT_WRITE, OUTPUT_SKIP, &
        OUTPUT_COMPUTE_ONLY, OUTPUT_ABORT
    implicit none
    private

    public :: active_pair_filter_workflow_type
    public :: configure_active_pair_filter_workflow
    public :: read_active_pair_filter_workflow
    public :: write_active_pair_filter_workflow
    public :: validate_active_pair_filter_workflow
    public :: check_active_pair_filter_workflow_inputs
    public :: run_active_pair_filter_workflow

contains

subroutine configure_active_pair_filter_workflow(job)
    type(active_pair_filter_workflow_type), intent(out) :: job
    character(len=1024) :: input, filename
    integer :: choice, ios, n_possible

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' Active CV-Pair Filter Workflow'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)',advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid active-pair filter workflow option'

    select case(choice)
    case(1)
        job = active_pair_filter_workflow_type()
        call edit_active_pair_filter_workflow(job)
    case(2,3)
        write(output_unit,'(A)',advance='no') 'Configuration filename: '
        read(*,'(A)') filename
        if (len_trim(filename) == 0) error stop 'No active-pair filter workflow filename specified'
        call read_active_pair_filter_workflow(trim(filename),job)
        if (choice == 2) then
            call edit_active_pair_filter_workflow(job)
        else
            call print_active_pair_filter_workflow(job)
            n_possible = check_active_pair_filter_workflow_inputs(job)
            if (job%n_pairs + job%n_random_pairs > n_possible) then
                write(error_unit,'(A,I0,A,I0,A,I0)') 'Configured ranked pairs: ',job%n_pairs, &
                    ', random pairs: ',job%n_random_pairs,', available: ',n_possible
                error stop 'Active-pair filter configuration is invalid'
            end if
            write(output_unit,'(A,I0)') 'Available unordered CV-object pairs: ',n_possible
        end if
    case(4)
        job = active_pair_filter_workflow_type()
        write(output_unit,*) 'Active-pair filter workflow configuration skipped.'
        return
    case default
        error stop 'Invalid active-pair filter workflow option'
    end select

    call validate_active_pair_filter_workflow(job)
    write(output_unit,'(A)',advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case(trim(adjustl(input)))
    case('y','Y','yes','YES','Yes')
        write(output_unit,'(A)',advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        filename = workflow_config_filename(trim(input),active_pair_filter_workflow_config_stem)
        call write_active_pair_filter_workflow(job,trim(filename))
        write(output_unit,'(3A)') 'Configuration saved to "',trim(filename),'"'
    case('n','N','no','NO','No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_active_pair_filter_workflow

subroutine edit_active_pair_filter_workflow(job)
    type(active_pair_filter_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios, n_possible, choice

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_string('TE matrix filename',job%matrix_filename)
    call edit_string('CV plan filename',job%cv_plan_filename)
    call edit_string('Active-pair matrix output filename',job%active_pair_filename)
    call edit_string('Output policy (ask/replace/new)',job%output_policy)
    write(output_unit,*) 'Matrix layout: 1 = one value per CV object, 2 = two CV components per object'
    choice = 2
    if (trim(job%matrix_layout) == 'one_per_object') choice = 1
    write(output_unit,'(A,I0,A)',advance='no') 'Select matrix layout [',choice,']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) choice
        if (ios /= 0 .or. choice < 1 .or. choice > 2) error stop 'Invalid pair-score matrix layout'
        if (choice == 1) then
            job%matrix_layout = 'one_per_object'
            job%score_mode = 'sum'
        else
            job%matrix_layout = 'two_per_object'
        end if
    end if
    call validate_active_pair_filter_workflow(job)
    n_possible = check_active_pair_filter_workflow_inputs(job)
    write(output_unit,'(A,I0)') 'Available unordered CV-object pairs: ',n_possible

    if (trim(job%matrix_layout) == 'two_per_object') then
        write(output_unit,*) 'Pair score: 1 = sum, 2 = sum of squares, 3 = largest absolute entry'
        choice = 1
        select case(trim(job%score_mode))
        case('sum_squares')
            choice = 2
        case('largest_abs')
            choice = 3
        end select
        write(output_unit,'(A,I0,A)',advance='no') 'Select pair score [',choice,']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) choice
            if (ios /= 0 .or. choice < 1 .or. choice > 3) error stop 'Invalid pair-score mode'
            select case(choice)
            case(1)
                job%score_mode = 'sum'
            case(2)
                job%score_mode = 'sum_squares'
            case(3)
                job%score_mode = 'largest_abs'
            end select
        end if
    end if
    call edit_string('Select pairs with largest or smallest scores',job%selection)
    call edit_integer('Number of pairs to select',job%n_pairs)
    call edit_integer('Additional random pairs to activate (0 for none)',job%n_random_pairs)
    if (job%n_random_pairs > 0) then
        call edit_string('File for additional random pair indices',job%random_pair_indices_filename)
    end if
    call validate_active_pair_filter_workflow(job)
    if (job%n_pairs + job%n_random_pairs > n_possible) &
        error stop 'Ranked and random pair counts exceed available CV-object pairs'
    call print_active_pair_filter_workflow(job)

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
        if (ios /= 0) error stop 'Invalid active-pair filter integer'
    end subroutine edit_integer
end subroutine edit_active_pair_filter_workflow

subroutine validate_active_pair_filter_workflow(job)
    type(active_pair_filter_workflow_type), intent(in) :: job

    if (len_trim(job%matrix_filename) == 0) error stop 'Pair-score matrix filename cannot be empty'
    if (len_trim(job%cv_plan_filename) == 0) error stop 'CV plan filename cannot be empty'
    if (len_trim(job%active_pair_filename) == 0) error stop 'Active-pair matrix filename cannot be empty'
    if (trim(job%active_pair_filename) == trim(job%matrix_filename) .or. &
        trim(job%active_pair_filename) == trim(job%cv_plan_filename)) &
        error stop 'Active-pair output must not overwrite an input file'
    select case(trim(job%matrix_layout))
    case('one_per_object','two_per_object')
    case default
        error stop 'Matrix layout must be one_per_object or two_per_object'
    end select
    select case(trim(job%score_mode))
    case('sum','sum_squares','largest_abs')
    case default
        error stop 'Pair score mode must be sum, sum_squares, or largest_abs'
    end select
    if (trim(job%matrix_layout) == 'one_per_object' .and. trim(job%score_mode) /= 'sum') &
        error stop 'One-value-per-object matrices support only the directional sum score'
    select case(trim(job%selection))
    case('largest','smallest')
    case default
        error stop 'Pair selection must be largest or smallest'
    end select
    if (job%n_pairs < 1) error stop 'Number of selected pairs must be positive'
    if (job%n_random_pairs < 0) error stop 'Number of random additional pairs cannot be negative'
    if (job%n_random_pairs > 0) then
        if (len_trim(job%random_pair_indices_filename) == 0) &
            error stop 'Random pair indices filename cannot be empty'
        if (trim(job%random_pair_indices_filename) == trim(job%active_pair_filename)) &
            error stop 'Random pair indices and active-pair matrix need different filenames'
        if (trim(job%random_pair_indices_filename) == trim(job%matrix_filename) .or. &
            trim(job%random_pair_indices_filename) == trim(job%cv_plan_filename)) &
            error stop 'Random pair index output must not overwrite an input file'
    end if
    select case(trim(job%output_policy))
    case('ask','replace','new')
    case default
        error stop 'Active-pair output policy must be ask, replace, or new'
    end select
end subroutine validate_active_pair_filter_workflow

integer function check_active_pair_filter_workflow_inputs(job,cv_plan_will_be_generated) result(n_possible)
    type(active_pair_filter_workflow_type), intent(in) :: job
    logical, intent(in), optional :: cv_plan_will_be_generated
    type(cv_plan) :: plan
    real(8), allocatable :: matrix(:,:)
    integer :: n_objects, matrix_stride
    logical :: plan_generated

    call validate_active_pair_filter_workflow(job)
    plan_generated = .false.
    if (present(cv_plan_will_be_generated)) plan_generated = cv_plan_will_be_generated
    matrix_stride = 2
    if (trim(job%matrix_layout) == 'one_per_object') matrix_stride = 1
    if (plan_generated) then
        call read_square_matrix_auto(trim(job%matrix_filename),matrix)
        if (.not. all(ieee_is_finite(matrix))) error stop 'Pair-score matrix contains non-finite values'
        if (mod(size(matrix,1),matrix_stride) /= 0) &
            error stop 'Pair-score matrix dimension does not match the selected layout'
        n_objects = size(matrix,1)/matrix_stride
    else
        call read_cv_plan(plan,trim(job%cv_plan_filename))
        if (.not. allocated(plan%object)) error stop 'CV plan has no object definitions'
        n_objects = size(plan%object)
        call read_square_matrix_auto(trim(job%matrix_filename),matrix)
        if (.not. all(ieee_is_finite(matrix))) error stop 'Pair-score matrix contains non-finite values'
        if (size(matrix,1) /= matrix_stride*n_objects) then
            write(error_unit,'(A,I0,A,I0)') 'Pair-score matrix dimension is ',size(matrix,1), &
                ', expected ',matrix_stride*n_objects
            error stop 'Pair-score matrix and CV plan dimensions do not match'
        end if
    end if
    if (n_objects < 2) error stop 'CV plan must contain at least two objects'
    n_possible = n_objects*(n_objects-1)/2
end function check_active_pair_filter_workflow_inputs

subroutine run_active_pair_filter_workflow(job)
    type(active_pair_filter_workflow_type), intent(in) :: job
    type(cv_plan) :: plan
    real(8), allocatable :: matrix(:,:), scores(:)
    integer, allocatable :: pair_i(:), pair_j(:)
    logical, allocatable :: pair_active(:,:)
    integer :: n_objects, n_possible, i, j, pair_index, matrix_stride, swap_index, temp_index
    integer :: output_action, random_output_action, unit, random_unit, ios
    real(8) :: random_value, temp_score

    n_possible = check_active_pair_filter_workflow_inputs(job)
    if (job%n_pairs + job%n_random_pairs > n_possible) then
        write(error_unit,'(A,I0,A,I0,A,I0)') 'Requested ',job%n_pairs,' ranked plus ', &
            job%n_random_pairs,' random pairs, but only ',n_possible,' unordered pairs are available'
        error stop 'Requested active-pair count exceeds the available pairs'
    end if
    call read_cv_plan(plan,trim(job%cv_plan_filename))
    call read_square_matrix_auto(trim(job%matrix_filename),matrix)
    n_objects = size(plan%object)
    matrix_stride = 2
    if (trim(job%matrix_layout) == 'one_per_object') matrix_stride = 1
    allocate(scores(n_possible),pair_i(n_possible),pair_j(n_possible))
    pair_index = 0
    do i = 1, n_objects-1
        do j = i+1, n_objects
            pair_index = pair_index + 1
            pair_i(pair_index) = i
            pair_j(pair_index) = j
            if (matrix_stride == 1) then
                scores(pair_index) = matrix(i,j) + matrix(j,i)
            else
                select case(trim(job%score_mode))
                case('sum')
                    scores(pair_index) = sum(matrix(2*i-1:2*i,2*j-1:2*j)) + &
                        sum(matrix(2*j-1:2*j,2*i-1:2*i))
                case('sum_squares')
                    scores(pair_index) = sum(matrix(2*i-1:2*i,2*j-1:2*j)**2) + &
                        sum(matrix(2*j-1:2*j,2*i-1:2*i)**2)
                case('largest_abs')
                    scores(pair_index) = max(maxval(abs(matrix(2*i-1:2*i,2*j-1:2*j))), &
                        maxval(abs(matrix(2*j-1:2*j,2*i-1:2*i))))
                end select
            end if
            if (.not. ieee_is_finite(scores(pair_index))) &
                error stop 'Pair-score matrix contains a non-finite value'
        end do
    end do
    call sort_pair_scores(scores,pair_i,pair_j,trim(job%selection) == 'largest')

    write(output_unit,'(A,I0,A,I0,A)') 'Selecting ',job%n_pairs,' of ',n_possible, &
        ' unordered CV-object pairs.'
    write(output_unit,'(A,A,A,A,A,ES16.8,A,ES16.8)') 'Matrix layout: ',trim(job%matrix_layout), &
        ', score: ',trim(job%score_mode),', selected score range: ', &
        minval(scores(1:job%n_pairs)),' to ', &
        maxval(scores(1:job%n_pairs))

    if (job%n_random_pairs > 0) then
        do i = n_possible,job%n_pairs+2,-1
            call random_number(random_value)
            swap_index = job%n_pairs+1+int(random_value*real(i-job%n_pairs,8))
            temp_index = pair_i(i)
            pair_i(i) = pair_i(swap_index)
            pair_i(swap_index) = temp_index
            temp_index = pair_j(i)
            pair_j(i) = pair_j(swap_index)
            pair_j(swap_index) = temp_index
            temp_score = scores(i)
            scores(i) = scores(swap_index)
            scores(swap_index) = temp_score
        end do
        write(output_unit,'(A,I0)') 'Additional unique random pairs: ',job%n_random_pairs
    end if

    allocate(pair_active(2*n_objects,2*n_objects))
    pair_active = .false.
    do pair_index = 1,job%n_pairs+job%n_random_pairs
        i = pair_i(pair_index)
        j = pair_j(pair_index)
        pair_active(2*i-1:2*i,2*j-1:2*j) = .true.
        pair_active(2*j-1:2*j,2*i-1:2*i) = .true.
    end do

    call check_output_policy(trim(job%active_pair_filename),trim(job%output_policy),output_action)
    if (output_action == OUTPUT_ABORT) error stop 'Active-pair matrix output policy aborted workflow'
    random_output_action = OUTPUT_WRITE
    if (job%n_random_pairs > 0) then
        call check_output_policy(trim(job%random_pair_indices_filename),trim(job%output_policy), &
            random_output_action)
        if (random_output_action == OUTPUT_ABORT) &
            error stop 'Random pair indices output policy aborted workflow'
    end if
    if (output_action == OUTPUT_SKIP .or. random_output_action == OUTPUT_SKIP) then
        write(output_unit,*) 'Active-pair outputs skipped by output policy.'
        return
    end if
    if (output_action == OUTPUT_COMPUTE_ONLY .or. random_output_action == OUTPUT_COMPUTE_ONLY) then
        write(output_unit,*) 'Active-pair matrix computed but not written by output policy.'
        return
    end if

    open(newunit=unit,file=trim(job%active_pair_filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write active-pair matrix: ',trim(job%active_pair_filename)
        error stop 'Active-pair matrix output failed'
    end if
    if (job%n_random_pairs > 0) then
        open(newunit=random_unit,file=trim(job%random_pair_indices_filename),status='replace', &
            action='write',iostat=ios)
        if (ios /= 0) then
            close(unit)
            write(error_unit,'(A,A)') 'Cannot write random pair index list: ', &
                trim(job%random_pair_indices_filename)
            error stop 'Random pair index output failed'
        end if
        write(random_unit,'(A)') '# Additional random CV-object pairs (1-based object indices)'
        write(random_unit,'(A,I0)') 'n_random_pairs = ',job%n_random_pairs
        do pair_index = job%n_pairs+1,job%n_pairs+job%n_random_pairs
            write(random_unit,'(I0,1X,I0)',iostat=ios) pair_i(pair_index),pair_j(pair_index)
            if (ios /= 0) then
                close(unit)
                close(random_unit)
                error stop 'Error writing random pair index list'
            end if
        end do
        close(random_unit,iostat=ios)
        if (ios /= 0) then
            close(unit)
            error stop 'Error closing random pair index list'
        end if
    end if
    do i = 1, size(pair_active,1)
        write(unit,'(*(L1,1X))',iostat=ios) pair_active(i,:)
        if (ios /= 0) then
            close(unit)
            error stop 'Error writing active-pair matrix'
        end if
    end do
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing active-pair matrix'
    write(output_unit,'(3A)') 'Active-pair matrix saved to "',trim(job%active_pair_filename),'"'
    if (job%n_random_pairs > 0) &
        write(output_unit,'(3A)') 'Random pair indices saved to "',trim(job%random_pair_indices_filename),'"'
end subroutine run_active_pair_filter_workflow

subroutine sort_pair_scores(scores,pair_i,pair_j,largest_first)
    real(8), intent(inout) :: scores(:)
    integer, intent(inout) :: pair_i(:),pair_j(:)
    logical, intent(in) :: largest_first

    call quicksort(1,size(scores))

contains
    recursive subroutine quicksort(left,right)
        integer, intent(in) :: left,right
        integer :: i,j,temp_i,temp_j
        real(8) :: pivot_score,temp_score
        integer :: pivot_i,pivot_j

        if (left >= right) return
        i = left
        j = right
        pivot_score = scores((left+right)/2)
        pivot_i = pair_i((left+right)/2)
        pivot_j = pair_j((left+right)/2)
        do
            do while (i <= right)
                if (.not. pair_precedes(scores(i),pair_i(i),pair_j(i), &
                    pivot_score,pivot_i,pivot_j,largest_first)) exit
                i = i+1
            end do
            do while (j >= left)
                if (.not. pair_precedes(pivot_score,pivot_i,pivot_j, &
                    scores(j),pair_i(j),pair_j(j),largest_first)) exit
                j = j-1
            end do
            if (i > j) exit
            temp_score = scores(i)
            scores(i) = scores(j)
            scores(j) = temp_score
            temp_i = pair_i(i)
            pair_i(i) = pair_i(j)
            pair_i(j) = temp_i
            temp_j = pair_j(i)
            pair_j(i) = pair_j(j)
            pair_j(j) = temp_j
            i = i+1
            j = j-1
            if (i > j) exit
        end do
        if (left < j) call quicksort(left,j)
        if (i < right) call quicksort(i,right)
    end subroutine quicksort

    logical function pair_precedes(score_a,i_a,j_a,score_b,i_b,j_b,descending)
        real(8), intent(in) :: score_a,score_b
        integer, intent(in) :: i_a,j_a,i_b,j_b
        logical, intent(in) :: descending

        if (score_a < score_b) then
            pair_precedes = .not. descending
        else if (score_a > score_b) then
            pair_precedes = descending
        else
            if (i_a == i_b) then
                pair_precedes = j_a < j_b
            else
                pair_precedes = i_a < i_b
            end if
        end if
    end function pair_precedes
end subroutine sort_pair_scores

subroutine write_active_pair_filter_workflow(job,filename)
    type(active_pair_filter_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit,ios

    call validate_active_pair_filter_workflow(job)
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) error stop 'Cannot write active-pair filter workflow configuration'
    write(unit,'(A)') '# Active CV-pair filter workflow'
    write(unit,'(A)') 'workflow_type = active_pair_filter'
    write(unit,'(A)') 'format_version = 2'
    write(unit,'(A,A)') 'matrix_filename = ',trim(job%matrix_filename)
    write(unit,'(A,A)') 'cv_plan_filename = ',trim(job%cv_plan_filename)
    write(unit,'(A,A)') 'matrix_layout = ',trim(job%matrix_layout)
    write(unit,'(A,A)') 'score_mode = ',trim(job%score_mode)
    write(unit,'(A,A)') 'selection = ',trim(job%selection)
    write(unit,'(A,I0)') 'n_pairs = ',job%n_pairs
    write(unit,'(A,I0)') 'n_random_pairs = ',job%n_random_pairs
    write(unit,'(A,A)') 'active_pair_filename = ',trim(job%active_pair_filename)
    write(unit,'(A,A)') 'random_pair_indices_filename = ',trim(job%random_pair_indices_filename)
    write(unit,'(A,A)') 'output_policy = ',trim(job%output_policy)
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing active-pair filter configuration'
end subroutine write_active_pair_filter_workflow

subroutine read_active_pair_filter_workflow(filename,job)
    character(len=*), intent(in) :: filename
    type(active_pair_filter_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1800) :: value
    integer :: unit,ios,equals_position,version
    logical :: saw_type,saw_version

    job = active_pair_filter_workflow_type()
    saw_type = .false.
    saw_version = .false.
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot read active-pair filter workflow: ',trim(filename)
        error stop 'Active-pair filter configuration input failed'
    end if
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading active-pair filter configuration'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line,'=')
        if (equals_position == 0) error stop 'Invalid active-pair filter configuration line'
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case(trim(key))
        case('workflow_type')
            if (trim(value) /= 'active_pair_filter') error stop 'Unexpected active-pair workflow type'
            saw_type = .true.
        case('format_version')
            version = 0
            read(value,*,iostat=ios) version
            if (ios /= 0 .or. (version /= 1 .and. version /= 2)) &
                error stop 'Unsupported active-pair filter workflow version'
            saw_version = .true.
        case('matrix_filename')
            job%matrix_filename = trim(value)
        case('cv_plan_filename')
            job%cv_plan_filename = trim(value)
        case('matrix_layout')
            job%matrix_layout = trim(value)
        case('score_mode')
            job%score_mode = trim(value)
        case('selection')
            job%selection = trim(value)
        case('n_pairs')
            read(value,*,iostat=ios) job%n_pairs
            if (ios /= 0) error stop 'Invalid selected pair count'
        case('n_random_pairs')
            read(value,*,iostat=ios) job%n_random_pairs
            if (ios /= 0) error stop 'Invalid random pair count'
        case('active_pair_filename')
            job%active_pair_filename = trim(value)
        case('random_pair_indices_filename')
            job%random_pair_indices_filename = trim(value)
        case('output_policy')
            job%output_policy = trim(value)
        case default
            error stop 'Unknown active-pair filter workflow configuration key'
        end select
    end do
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing active-pair filter configuration'
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete active-pair filter configuration'
    call validate_active_pair_filter_workflow(job)
end subroutine read_active_pair_filter_workflow

subroutine print_active_pair_filter_workflow(job)
    type(active_pair_filter_workflow_type), intent(in) :: job

    write(output_unit,'(A,A)') 'Pair-score matrix: ',trim(job%matrix_filename)
    write(output_unit,'(A,A)') 'CV plan: ',trim(job%cv_plan_filename)
    write(output_unit,'(A,A)') 'Matrix layout: ',trim(job%matrix_layout)
    write(output_unit,'(A,A)') 'Pair score: ',trim(job%score_mode)
    write(output_unit,'(A,A,A,I0)') 'Selection: ',trim(job%selection),' ',job%n_pairs
    write(output_unit,'(A,A)') 'Active-pair output: ',trim(job%active_pair_filename)
    write(output_unit,'(A,I0)') 'Additional random pairs: ',job%n_random_pairs
    if (job%n_random_pairs > 0) &
        write(output_unit,'(A,A)') 'Random pair indices: ',trim(job%random_pair_indices_filename)
    write(output_unit,'(A,A)') 'Output policy: ',trim(job%output_policy)
end subroutine print_active_pair_filter_workflow

end module active_pair_filter_workflow_module
