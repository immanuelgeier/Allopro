module sh_invariant_contraction_workflow_module
    use iso_fortran_env, only: output_unit, error_unit, real64, int64
    use workflow_module, only: sh_invariant_contraction_workflow_type, workflow_config_filename, &
        sh_invariant_contraction_workflow_config_stem
    use sh_product_sampling_workflow_module, only: sh_product_means_type, sh_product_mean_record_type, &
        SH_MEAN_SINGLE, SH_MEAN_PAIR, SH_MEAN_TRIPLE, SH_MEAN_SINGLE_FUTURE, &
        SH_MEAN_PAIR_ZERO_LAG, read_sh_product_means
    use output_policy_module, only: get_output_prefix, check_output_policy, OUTPUT_WRITE, &
        OUTPUT_SKIP, OUTPUT_COMPUTE_ONLY, OUTPUT_ABORT
    implicit none
    private

    public :: sh_invariant_contraction_workflow_type
    public :: configure_sh_invariant_contraction_workflow
    public :: read_sh_invariant_contraction_workflow
    public :: write_sh_invariant_contraction_workflow
    public :: validate_sh_invariant_contraction_workflow
    public :: check_sh_invariant_contraction_inputs
    public :: run_sh_invariant_contraction_workflow
    public :: wigner_3j

contains

subroutine configure_sh_invariant_contraction_workflow(job)
    type(sh_invariant_contraction_workflow_type), intent(out) :: job
    character(len=1024) :: input, filename
    integer :: choice, ios

    write(output_unit,*)
    write(output_unit,*) '======================================'
    write(output_unit,*) ' SH Wigner-3j Invariant Contraction'
    write(output_unit,*) '======================================'
    write(output_unit,*) '1 = Create new workflow'
    write(output_unit,*) '2 = Load and edit existing workflow'
    write(output_unit,*) '3 = Load existing workflow unchanged'
    write(output_unit,*) '4 = Skip workflow configuration'
    write(output_unit,'(A)', advance='no') 'Select option: '
    read(*,'(A)') input
    read(input,*,iostat=ios) choice
    if (ios /= 0) error stop 'Invalid SH invariant contraction workflow option'

    select case (choice)
    case (1)
        job = sh_invariant_contraction_workflow_type()
        call edit_sh_invariant_contraction_workflow(job)
    case (2, 3)
        write(output_unit,'(A)', advance='no') 'Configuration filename: '
        read(*,'(A)') filename
        if (len_trim(filename) == 0) error stop 'No configuration filename specified'
        call read_sh_invariant_contraction_workflow(trim(filename), job)
        if (choice == 2) then
            call edit_sh_invariant_contraction_workflow(job)
        else
            call print_sh_invariant_contraction_workflow(job)
        end if
    case (4)
        job = sh_invariant_contraction_workflow_type()
        write(output_unit,*) 'SH invariant contraction configuration skipped.'
        return
    case default
        error stop 'Invalid SH invariant contraction workflow option'
    end select

    call validate_sh_invariant_contraction_workflow(job)
    write(output_unit,'(A)', advance='no') 'Save workflow configuration? (y/n) [y]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'y'
    select case (trim(adjustl(input)))
    case ('y', 'Y', 'yes', 'YES', 'Yes')
        write(output_unit,'(A)', advance='no') 'Workflow file prefix (blank for no prefix): '
        read(*,'(A)') input
        filename = workflow_config_filename(trim(input), sh_invariant_contraction_workflow_config_stem)
        call write_sh_invariant_contraction_workflow(job, trim(filename))
        write(output_unit,'(3A)') 'Configuration saved to "', trim(filename), '"'
    case ('n', 'N', 'no', 'NO', 'No')
        continue
    case default
        error stop 'Invalid answer; enter y or n'
    end select
end subroutine configure_sh_invariant_contraction_workflow

subroutine edit_sh_invariant_contraction_workflow(job)
    type(sh_invariant_contraction_workflow_type), intent(inout) :: job
    character(len=1024) :: input
    integer :: ios, n_triplets, i

    write(output_unit,*)
    write(output_unit,*) 'Press Enter to keep the current value.'
    call edit_string('SH product means filename', job%means_filename)
    call edit_string('Output directory', job%output_path)
    call edit_string('Output filename prefix', job%output_prefix)
    call edit_string('Contractions (all/selected)', job%contraction_mode)
    if (trim(job%contraction_mode) == 'selected') then
        n_triplets = job%n_selected_l_triplets
        write(output_unit,'(A,I0,A)', advance='no') 'Number of selected (l1,l2,l3) triplets [', &
            n_triplets, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) n_triplets
            if (ios /= 0 .or. n_triplets < 1) error stop 'Invalid number of SH contractions'
            if (allocated(job%selected_l_triplets)) deallocate(job%selected_l_triplets)
            allocate(job%selected_l_triplets(3,n_triplets))
            job%n_selected_l_triplets = n_triplets
            do i = 1, n_triplets
                write(output_unit,'(A,I0,A)', advance='no') 'Enter l1 l2 l3 for contraction ', i, ': '
                read(*,'(A)') input
                read(input,*,iostat=ios) job%selected_l_triplets(:,i)
                if (ios /= 0) error stop 'Invalid (l1,l2,l3) contraction'
            end do
        else if (n_triplets > 0 .and. .not. allocated(job%selected_l_triplets)) then
            error stop 'Selected contraction triplets are missing'
        end if
    end if
    call edit_string('Output policy (ask/replace/new)', job%output_policy)
    call validate_sh_invariant_contraction_workflow(job)
    call print_sh_invariant_contraction_workflow(job)

contains
    subroutine edit_string(prompt, value)
        character(len=*), intent(in) :: prompt
        character(len=*), intent(inout) :: value
        write(output_unit,'(A,A,A)', advance='no') trim(prompt), ' [', trim(value), ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) value = trim(input)
    end subroutine edit_string
end subroutine edit_sh_invariant_contraction_workflow

subroutine validate_sh_invariant_contraction_workflow(job)
    type(sh_invariant_contraction_workflow_type), intent(in) :: job

    if (len_trim(job%means_filename) == 0) error stop 'SH product means filename cannot be empty'
    if (len_trim(job%output_path) == 0) error stop 'SH invariant output path cannot be empty'
    if (len_trim(job%output_prefix) == 0) error stop 'SH invariant output prefix cannot be empty'
    select case (trim(job%contraction_mode))
    case ('all')
        if (job%n_selected_l_triplets /= 0) &
            error stop 'All-contractions mode cannot contain selected l triplets'
    case ('selected')
        if (job%n_selected_l_triplets < 1) error stop 'At least one selected l triplet is required'
        if (.not. allocated(job%selected_l_triplets)) error stop 'Selected l triplets are not allocated'
        if (size(job%selected_l_triplets,1) /= 3 .or. &
            size(job%selected_l_triplets,2) /= job%n_selected_l_triplets) then
            error stop 'Selected l-triplet array has an invalid shape'
        end if
    case default
        error stop 'Contraction mode must be all or selected'
    end select
    select case (trim(job%output_policy))
    case ('ask', 'replace', 'new')
    case default
        error stop 'SH invariant output policy must be ask, replace, or new'
    end select
end subroutine validate_sh_invariant_contraction_workflow

subroutine run_sh_invariant_contraction_workflow(job)
    type(sh_invariant_contraction_workflow_type), intent(in) :: job
    type(sh_product_means_type) :: means
    integer, allocatable :: record_hash(:)
    logical, allocatable :: active_pairs(:,:)
    integer, allocatable :: l_triplets(:,:)
    integer :: n_triplets, i, j, l1, l2, l3, m1, m2, m3, index
    integer :: record_index, unit, ios, output_action, pair_count
    complex(real64) :: raw_value, connected_value, triple_mean
    complex(real64) :: pair_ab, pair_ac, pair_bc, single_a, single_b, single_c
    real(real64) :: coefficient
    character(len=2048) :: output_directory, output_name, output_filename

    call validate_sh_invariant_contraction_workflow(job)
    call read_sh_product_means(trim(job%means_filename), means)
    call validate_means_for_connected_contraction(means)
    call collect_active_pairs(means, active_pairs)
    call build_record_hash(means%records, record_hash)
    call resolve_l_triplets(job, means%lmax, l_triplets)
    n_triplets = size(l_triplets,2)
    pair_count = count(active_pairs)

    output_directory = trim(job%output_path)
    if (output_directory(len_trim(output_directory):len_trim(output_directory)) /= '/') then
        output_directory = trim(output_directory)//'/'
    end if
    output_name = get_output_prefix(trim(job%output_prefix), trim(job%output_policy))

    do index = 1, n_triplets
        l1 = l_triplets(1,index)
        l2 = l_triplets(2,index)
        l3 = l_triplets(3,index)
        write(output_filename,'(A,A,"_l1_",I0,"_l2_",I0,"_l3_",I0,".txt")') &
            trim(output_directory), trim(output_name), l1, l2, l3
        call check_output_policy(trim(output_filename), trim(job%output_policy), output_action)
        if (output_action == OUTPUT_ABORT) error stop 'SH invariant output policy aborted workflow'
        if (output_action == OUTPUT_SKIP) cycle
        if (output_action == OUTPUT_WRITE) then
            open(newunit=unit, file=trim(output_filename), status='replace', action='write', iostat=ios)
            if (ios /= 0) call output_error('Cannot write SH invariant matrix', trim(output_filename))
            write(unit,'(A)') '# Wigner-3j contracted SH product means'
            write(unit,'(A)') 'format_version = 1'
            write(unit,'(A,I0)') 'state_index = ', means%state_index
            write(unit,'(A,I0)') 'lag = ', means%lag
            write(unit,'(A,I0)') 'l1 = ', l1
            write(unit,'(A,I0)') 'l2 = ', l2
            write(unit,'(A,I0)') 'l3 = ', l3
            write(unit,'(A)') '# channel_i channel_j raw_re raw_im connected_re connected_im'
        else
            unit = 0
        end if

        do i = 1, means%n_channels
            do j = 1, means%n_channels
                if (.not. active_pairs(i,j)) cycle
                raw_value = (0.0_real64,0.0_real64)
                connected_value = (0.0_real64,0.0_real64)
                do m1 = -l1, l1
                    do m2 = -l2, l2
                        m3 = -m1-m2
                        if (abs(m3) > l3) cycle
                        coefficient = wigner_3j(l1,l2,l3,m1,m2,m3)
                        if (coefficient == 0.0_real64) cycle
                        record_index = find_record(record_hash, means%records, SH_MEAN_TRIPLE, &
                            i, i, j, l1, m1, l2, m2, l3, m3)
                        triple_mean = record_value(means%records(record_index))
                        pair_ab = fetch_record_value(record_hash, means%records, SH_MEAN_PAIR, &
                            i, i, 0, l1, m1, l2, m2, 0, 0)
                        pair_ac = fetch_record_value(record_hash, means%records, SH_MEAN_PAIR, &
                            i, j, 0, l1, m1, l3, m3, 0, 0)
                        pair_bc = fetch_record_value(record_hash, means%records, SH_MEAN_PAIR_ZERO_LAG, &
                            i, j, 0, l2, m2, l3, m3, 0, 0)
                        single_a = fetch_record_value(record_hash, means%records, SH_MEAN_SINGLE_FUTURE, &
                            i, 0, 0, l1, m1, 0, 0, 0, 0)
                        single_b = fetch_record_value(record_hash, means%records, SH_MEAN_SINGLE, &
                            i, 0, 0, l2, m2, 0, 0, 0, 0)
                        single_c = fetch_record_value(record_hash, means%records, SH_MEAN_SINGLE, &
                            j, 0, 0, l3, m3, 0, 0, 0, 0)
                        raw_value = raw_value + coefficient*triple_mean
                        connected_value = connected_value + coefficient * (triple_mean - &
                            pair_ab*single_c - pair_ac*single_b - pair_bc*single_a + &
                            2.0_real64*single_a*single_b*single_c)
                    end do
                end do
                if (output_action == OUTPUT_WRITE) then
                    write(unit,'(2(I0,1X),4(ES24.16E3,1X))',iostat=ios) i, j, &
                        real(raw_value,real64), aimag(raw_value), &
                        real(connected_value,real64), aimag(connected_value)
                    if (ios /= 0) error stop 'Error writing SH invariant matrix row'
                end if
            end do
        end do
        if (output_action == OUTPUT_WRITE) then
            close(unit, iostat=ios)
            if (ios /= 0) error stop 'Error closing SH invariant matrix file'
            write(output_unit,'(A,A)') 'Wrote contracted matrix: ', trim(output_filename)
        end if
    end do
    write(output_unit,'(A,I0,A,I0)') 'Contracted ', n_triplets, ' l-triplets over ', &
        pair_count, ' selected channel pairs.'
end subroutine run_sh_invariant_contraction_workflow

subroutine check_sh_invariant_contraction_inputs(job)
    type(sh_invariant_contraction_workflow_type), intent(in) :: job
    type(sh_product_means_type) :: means
    logical, allocatable :: active_pairs(:,:)
    integer, allocatable :: record_hash(:), l_triplets(:,:)

    call validate_sh_invariant_contraction_workflow(job)
    call read_sh_product_means(trim(job%means_filename), means)
    call validate_means_for_connected_contraction(means)
    call collect_active_pairs(means, active_pairs)
    call build_record_hash(means%records, record_hash)
    call resolve_l_triplets(job, means%lmax, l_triplets)
end subroutine check_sh_invariant_contraction_inputs

subroutine validate_means_for_connected_contraction(means)
    type(sh_product_means_type), intent(in) :: means
    integer :: i
    logical :: has_future, has_zero_lag

    if (.not. allocated(means%records)) error stop 'SH means records are not loaded'
    if (means%format_version /= 4) then
        error stop 'SH means use an obsolete mode convention; rerun SH product sampling'
    end if
    has_future = .false.
    has_zero_lag = .false.
    do i = 1, size(means%records)
        has_future = has_future .or. means%records(i)%kind == SH_MEAN_SINGLE_FUTURE
        has_zero_lag = has_zero_lag .or. means%records(i)%kind == SH_MEAN_PAIR_ZERO_LAG
    end do
    if (.not. has_future .or. .not. has_zero_lag) then
        error stop 'Connected contraction requires SH means generated by the updated lagged sampler'
    end if
end subroutine validate_means_for_connected_contraction

subroutine collect_active_pairs(means, active_pairs)
    type(sh_product_means_type), intent(in) :: means
    logical, allocatable, intent(out) :: active_pairs(:,:)
    integer :: i, c1, c2, c3

    allocate(active_pairs(means%n_channels,means%n_channels))
    active_pairs = .false.
    do i = 1, size(means%records)
        if (means%records(i)%kind /= SH_MEAN_TRIPLE) cycle
        c1 = means%records(i)%channel1
        c2 = means%records(i)%channel2
        c3 = means%records(i)%channel3
        if (c1 < 1 .or. c1 > means%n_channels .or. c2 /= c1 .or. &
            c3 < 1 .or. c3 > means%n_channels) error stop 'Invalid SH triple channel tuple'
        active_pairs(c1,c3) = .true.
    end do
    if (.not. any(active_pairs)) error stop 'SH means contain no selected channel pairs'
end subroutine collect_active_pairs

subroutine resolve_l_triplets(job, lmax, triplets)
    type(sh_invariant_contraction_workflow_type), intent(in) :: job
    integer, intent(in) :: lmax
    integer, allocatable, intent(out) :: triplets(:,:)
    integer :: l1, l2, l3, n, index

    if (trim(job%contraction_mode) == 'selected') then
        triplets = job%selected_l_triplets
        do index = 1, size(triplets,2)
            l1 = triplets(1,index)
            l2 = triplets(2,index)
            l3 = triplets(3,index)
            call validate_l_triplet(l1,l2,l3,lmax)
            if (any(all(triplets(:,1:index-1) == spread(triplets(:,index),2,index-1),dim=1))) then
                error stop 'Duplicate selected l triplet'
            end if
        end do
        return
    end if

    n = 0
    do l1 = 0, lmax
        do l2 = 0, lmax
            n = n + 2*min(l1,l2) + 1
        end do
    end do
    allocate(triplets(3,n))
    index = 0
    do l1 = 0, lmax
        do l2 = 0, lmax
            do l3 = abs(l1-l2), l1+l2
                index = index + 1
                triplets(:,index) = [l1,l2,l3]
            end do
        end do
    end do
end subroutine resolve_l_triplets

subroutine validate_l_triplet(l1,l2,l3,lmax)
    integer, intent(in) :: l1,l2,l3,lmax
    if (l1 < 0 .or. l1 > lmax .or. l2 < 0 .or. l2 > lmax) then
        error stop 'Selected l1/l2 exceeds the sampled SH maximum degree'
    end if
    if (l3 < abs(l1-l2) .or. l3 > l1+l2) error stop 'Selected l triplet violates triangle rule'
end subroutine validate_l_triplet

subroutine build_record_hash(records, table)
    type(sh_product_mean_record_type), intent(in) :: records(:)
    integer, allocatable, intent(out) :: table(:)
    integer(int64) :: requested_size
    integer :: i, slot, table_size

    if (size(records) > (huge(0)-1)/2) error stop 'Too many SH means records to index'
    requested_size = 2_int64*int(size(records),int64)+1_int64
    table_size = int(requested_size)
    allocate(table(table_size))
    table = 0
    do i = 1, size(records)
        slot = hash_slot(records(i)%kind, records(i)%channel1, records(i)%channel2, &
            records(i)%channel3, records(i)%l1, records(i)%m1, records(i)%l2, records(i)%m2, &
            records(i)%l3, records(i)%m3, table_size)
        do
            if (table(slot) == 0) then
                table(slot) = i
                exit
            end if
            if (same_record_key(records(i), records(table(slot)))) then
                error stop 'Duplicate SH means record key'
            end if
            slot = slot + 1
            if (slot > table_size) slot = 1
        end do
    end do
end subroutine build_record_hash

integer function find_record(table, records, kind, c1, c2, c3, l1, m1, l2, m2, l3, m3) result(index)
    integer, intent(in) :: table(:), kind, c1, c2, c3, l1, m1, l2, m2, l3, m3
    type(sh_product_mean_record_type), intent(in) :: records(:)
    integer :: slot, attempts

    slot = hash_slot(kind,c1,c2,c3,l1,m1,l2,m2,l3,m3,size(table))
    do attempts = 1, size(table)
        index = table(slot)
        if (index == 0) error stop 'Required SH mean is missing for connected contraction'
        if (same_key(records(index),kind,c1,c2,c3,l1,m1,l2,m2,l3,m3)) return
        slot = slot + 1
        if (slot > size(table)) slot = 1
    end do
    error stop 'SH means index search failed'
end function find_record

function fetch_record_value(table, records, kind, c1, c2, c3, l1, m1, l2, m2, l3, m3) result(value)
    integer, intent(in) :: table(:), kind, c1, c2, c3, l1, m1, l2, m2, l3, m3
    type(sh_product_mean_record_type), intent(in) :: records(:)
    complex(real64) :: value
    integer :: index

    index = find_record(table,records,kind,c1,c2,c3,l1,m1,l2,m2,l3,m3)
    value = record_value(records(index))
end function fetch_record_value

function record_value(record) result(value)
    type(sh_product_mean_record_type), intent(in) :: record
    complex(real64) :: value
    value = cmplx(record%mean_re,record%mean_im,real64)
end function record_value

integer function hash_slot(kind,c1,c2,c3,l1,m1,l2,m2,l3,m3,table_size) result(slot)
    integer, intent(in) :: kind,c1,c2,c3,l1,m1,l2,m2,l3,m3,table_size
    integer(int64), parameter :: modulus = 2147483629_int64
    integer(int64) :: hash
    integer :: values(10), i

    values = [kind,c1,c2,c3,l1,m1+256,l2,m2+256,l3,m3+256]
    hash = 17_int64
    do i = 1, size(values)
        hash = modulo(hash*131_int64+int(values(i),int64),modulus)
    end do
    slot = int(modulo(hash,int(table_size,int64))) + 1
end function hash_slot

logical function same_record_key(a,b)
    type(sh_product_mean_record_type), intent(in) :: a,b
    same_record_key = same_key(a,b%kind,b%channel1,b%channel2,b%channel3, &
        b%l1,b%m1,b%l2,b%m2,b%l3,b%m3)
end function same_record_key

logical function same_key(record,kind,c1,c2,c3,l1,m1,l2,m2,l3,m3)
    type(sh_product_mean_record_type), intent(in) :: record
    integer, intent(in) :: kind,c1,c2,c3,l1,m1,l2,m2,l3,m3
    same_key = record%kind == kind .and. record%channel1 == c1 .and. &
        record%channel2 == c2 .and. record%channel3 == c3 .and. &
        record%l1 == l1 .and. record%m1 == m1 .and. record%l2 == l2 .and. &
        record%m2 == m2 .and. record%l3 == l3 .and. record%m3 == m3
end function same_key

real(real64) function wigner_3j(j1,j2,j3,m1,m2,m3) result(value)
    integer, intent(in) :: j1,j2,j3,m1,m2,m3
    integer :: z, z_min, z_max, phase
    real(real64) :: log_prefactor, log_denominator, term

    value = 0.0_real64
    if (min(j1,j2,j3) < 0) return
    if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3) return
    if (m1+m2+m3 /= 0) return
    if (j3 < abs(j1-j2) .or. j3 > j1+j2) return

    log_prefactor = 0.5_real64 * ( &
        log_factorial(j1+j2-j3) + log_factorial(j1-j2+j3) + &
        log_factorial(-j1+j2+j3) - log_factorial(j1+j2+j3+1) + &
        log_factorial(j1+m1) + log_factorial(j1-m1) + &
        log_factorial(j2+m2) + log_factorial(j2-m2) + &
        log_factorial(j3+m3) + log_factorial(j3-m3))
    z_min = max(0,j2-j3-m1,j1-j3+m2)
    z_max = min(j1+j2-j3,j1-m1,j2+m2)
    if (z_min > z_max) return
    do z = z_min, z_max
        log_denominator = log_factorial(z) + log_factorial(j1+j2-j3-z) + &
            log_factorial(j1-m1-z) + log_factorial(j2+m2-z) + &
            log_factorial(j3-j2+m1+z) + log_factorial(j3-j1-m2+z)
        term = exp(log_prefactor-log_denominator)
        if (modulo(z,2) == 1) term = -term
        value = value + term
    end do
    phase = j1-j2-m3
    if (modulo(phase,2) /= 0) value = -value
end function wigner_3j

real(real64) function log_factorial(n) result(value)
    integer, intent(in) :: n
    if (n < 0) error stop 'Negative factorial argument in Wigner-3j calculation'
    value = log_gamma(real(n+1,real64))
end function log_factorial

subroutine write_sh_invariant_contraction_workflow(job, filename)
    type(sh_invariant_contraction_workflow_type), intent(in) :: job
    character(len=*), intent(in) :: filename
    integer :: unit, ios, i

    call validate_sh_invariant_contraction_workflow(job)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call output_error('Cannot write SH invariant workflow',trim(filename))
    write(unit,'(A)') '# SH Wigner-3j invariant contraction workflow'
    write(unit,'(A)') 'workflow_type = sh_invariant_contraction'
    write(unit,'(A)') 'format_version = 1'
    write(unit,'(A,A)') 'means_filename = ',trim(job%means_filename)
    write(unit,'(A,A)') 'output_path = ',trim(job%output_path)
    write(unit,'(A,A)') 'output_prefix = ',trim(job%output_prefix)
    write(unit,'(A,A)') 'contraction_mode = ',trim(job%contraction_mode)
    write(unit,'(A,I0)') 'n_selected_l_triplets = ',job%n_selected_l_triplets
    if (allocated(job%selected_l_triplets)) then
        do i = 1, size(job%selected_l_triplets,2)
            write(unit,*) job%selected_l_triplets(:,i)
        end do
    end if
    write(unit,'(A,A)') 'output_policy = ',trim(job%output_policy)
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing SH invariant workflow file'
end subroutine write_sh_invariant_contraction_workflow

subroutine read_sh_invariant_contraction_workflow(filename, job)
    character(len=*), intent(in) :: filename
    type(sh_invariant_contraction_workflow_type), intent(out) :: job
    character(len=2048) :: line
    character(len=128) :: key
    character(len=1800) :: value
    integer :: unit, ios, equals_position, version, i
    logical :: saw_type, saw_version, reading_triplets

    job = sh_invariant_contraction_workflow_type()
    saw_type = .false.
    saw_version = .false.
    reading_triplets = .false.
    i = 0
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios /= 0) call output_error('Cannot read SH invariant workflow',trim(filename))
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) exit
        if (ios /= 0) error stop 'Error reading SH invariant workflow'
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        equals_position = index(line,'=')
        if (equals_position == 0) then
            if (.not. reading_triplets) error stop 'Invalid SH invariant workflow line'
            i = i+1
            if (i > job%n_selected_l_triplets) error stop 'Too many selected l triplets'
            read(line,*,iostat=ios) job%selected_l_triplets(:,i)
            if (ios /= 0) error stop 'Invalid selected l triplet in workflow'
            cycle
        end if
        key = trim(adjustl(line(:equals_position-1)))
        value = trim(adjustl(line(equals_position+1:)))
        select case(trim(key))
        case('workflow_type')
            if (trim(value) /= 'sh_invariant_contraction') error stop 'Unexpected SH invariant workflow type'
            saw_type = .true.
        case('format_version')
            read(value,*,iostat=ios) version
            if (ios /= 0 .or. version /= 1) error stop 'Unsupported SH invariant workflow version'
            saw_version = .true.
        case('means_filename')
            job%means_filename = trim(value)
        case('output_path')
            job%output_path = trim(value)
        case('output_prefix')
            job%output_prefix = trim(value)
        case('contraction_mode')
            job%contraction_mode = trim(value)
        case('n_selected_l_triplets')
            read(value,*,iostat=ios) job%n_selected_l_triplets
            if (ios /= 0 .or. job%n_selected_l_triplets < 0) error stop 'Invalid selected l-triplet count'
            if (allocated(job%selected_l_triplets)) deallocate(job%selected_l_triplets)
            if (job%n_selected_l_triplets > 0) allocate(job%selected_l_triplets(3,job%n_selected_l_triplets))
            i = 0
            reading_triplets = job%n_selected_l_triplets > 0
        case('output_policy')
            job%output_policy = trim(value)
            reading_triplets = .false.
        case default
            error stop 'Unknown SH invariant workflow configuration key'
        end select
    end do
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing SH invariant workflow input'
    if (.not. saw_type .or. .not. saw_version) error stop 'Incomplete SH invariant workflow header'
    if (i /= job%n_selected_l_triplets) then
        if (job%n_selected_l_triplets > 0) error stop 'Missing selected l triplets'
    end if
    call validate_sh_invariant_contraction_workflow(job)
end subroutine read_sh_invariant_contraction_workflow

subroutine print_sh_invariant_contraction_workflow(job)
    type(sh_invariant_contraction_workflow_type), intent(in) :: job
    write(output_unit,'(A,A)') 'SH product means: ',trim(job%means_filename)
    write(output_unit,'(A,A)') 'Contractions: ',trim(job%contraction_mode)
    write(output_unit,'(A,I0)') 'Selected l triplets: ',job%n_selected_l_triplets
    write(output_unit,'(A,A,A)') 'Output prefix: ',trim(job%output_path),'/'//trim(job%output_prefix)
end subroutine print_sh_invariant_contraction_workflow

subroutine output_error(message, filename)
    character(len=*), intent(in) :: message, filename
    write(error_unit,'(A,A)') trim(message)//': ',trim(filename)
    error stop 'SH invariant contraction I/O failed'
end subroutine output_error

end module sh_invariant_contraction_workflow_module
