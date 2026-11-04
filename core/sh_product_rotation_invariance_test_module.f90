module sh_product_rotation_invariance_test_module
    use iso_fortran_env, only: real64, output_unit
    use workflow_module, only: sh_product_sampling_workflow_type, &
        sh_invariant_contraction_workflow_type
    use sh_product_sampling_workflow_module, only: run_sh_product_sampling_workflow
    use lifted_trajectory_io_module, only: lifted_trajectory, write_trajectory
    use lifting_plan_module, only: lifting_plan_type, write_lifting_plan
    use lifting_plan_execution_module, only: make_filename_lifted_tra
    use sh_invariant_contraction_workflow_module, only: run_sh_invariant_contraction_workflow
    implicit none
    private

    public :: run_sh_product_rotation_invariance_test

    integer, parameter :: test_lmax = 1
    integer, parameter :: n_samples = 23
    real(real64), parameter :: tolerance = 1.0e-10_real64

contains

subroutine run_sh_product_rotation_invariance_test()
    type(sh_invariant_contraction_workflow_type) :: job
    real(real64) :: rotation(3,3)
    complex(real64) :: raw_base, connected_base, raw_rotated, connected_rotated
    complex(real64) :: raw_base_l2, connected_base_l2, raw_rotated_l2, connected_rotated_l2
    character(len=64) :: test_prefix
    character(len=128) :: base_means_file, rotated_means_file
    character(len=160) :: base_output_file, rotated_output_file
    character(len=160) :: base_output_l2_file, rotated_output_l2_file
    character(len=24) :: clock_text
    integer :: clock_count

    call system_clock(count=clock_count)
    write(clock_text,'(I0)') clock_count
    test_prefix = 'sh_rotation_test_'//trim(clock_text)
    base_means_file = trim(test_prefix)//'_base_means.txt'
    rotated_means_file = trim(test_prefix)//'_rotated_means.txt'
    base_output_file = trim(test_prefix)//'_base_l1_1_l2_1_l3_1.txt'
    rotated_output_file = trim(test_prefix)//'_rotated_l1_1_l2_1_l3_1.txt'
    base_output_l2_file = trim(test_prefix)//'_base_l1_1_l2_1_l3_2.txt'
    rotated_output_l2_file = trim(test_prefix)//'_rotated_l1_1_l2_1_l3_2.txt'
    call make_rotation(rotation)

    call run_sampling(identity_matrix(),trim(test_prefix)//'_base',base_means_file)
    call configure_test_job(job, trim(base_means_file), trim(test_prefix)//'_base')
    call run_sh_invariant_contraction_workflow(job)
    call read_invariant_value(trim(base_output_file), raw_base, connected_base)
    call read_invariant_value(trim(base_output_l2_file), raw_base_l2, connected_base_l2)

    call run_sampling(rotation,trim(test_prefix)//'_rotated',rotated_means_file)
    call configure_test_job(job, trim(rotated_means_file), trim(test_prefix)//'_rotated')
    call run_sh_invariant_contraction_workflow(job)
    call read_invariant_value(trim(rotated_output_file), raw_rotated, connected_rotated)
    call read_invariant_value(trim(rotated_output_l2_file), raw_rotated_l2, connected_rotated_l2)

    if (abs(raw_base) <= tolerance) error stop 'SH rotation test generated a trivial raw invariant'
    if (abs(connected_base) <= tolerance) &
        error stop 'SH rotation test generated a trivial connected invariant'
    call assert_invariant(raw_base, raw_rotated, 'raw')
    call assert_invariant(connected_base, connected_rotated, 'connected')
    if (abs(raw_base_l2) <= tolerance) error stop 'SH rotation test generated a trivial L=2 raw invariant'
    if (abs(connected_base_l2) <= tolerance) &
        error stop 'SH rotation test generated a trivial L=2 connected invariant'
    call assert_invariant(raw_base_l2, raw_rotated_l2, 'L=2 raw')
    call assert_invariant(connected_base_l2, connected_rotated_l2, 'L=2 connected')

    call delete_file(trim(base_means_file))
    call delete_file(trim(rotated_means_file))
    call delete_file(trim(base_output_file))
    call delete_file(trim(rotated_output_file))
    call delete_file(trim(base_output_l2_file))
    call delete_file(trim(rotated_output_l2_file))
    call cleanup_sampling_inputs(trim(test_prefix)//'_base')
    call cleanup_sampling_inputs(trim(test_prefix)//'_rotated')
    write(output_unit,'(A,ES12.4)') 'SH product rotation-invariance test passed; tolerance = ', tolerance
end subroutine run_sh_product_rotation_invariance_test

subroutine run_sampling(rotation, prefix, means_filename)
    real(real64), intent(in) :: rotation(3,3)
    character(len=*), intent(in) :: prefix
    character(len=*), intent(out) :: means_filename
    type(sh_product_sampling_workflow_type) :: job
    type(lifting_plan_type) :: plan
    type(lifted_trajectory) :: trajectory
    character(len=256) :: selector_filename, plan_filename, pair_filename
    character(len=256) :: trajectory_prefix, trajectory_filename
    integer :: unit, frame

    selector_filename = trim(prefix)//'_selector.txt'
    plan_filename = trim(prefix)//'_lifting_plan.txt'
    pair_filename = trim(prefix)//'_active_pairs.txt'
    trajectory_prefix = trim(prefix)//'_trajectory'
    means_filename = trim(prefix)//'_products_means.txt'

    open(newunit=unit,file=trim(selector_filename),status='replace',action='write')
    write(unit,'(A)') 'state_index = 1'
    write(unit,'(A)') 'n_replicates = 1'
    write(unit,'(A)') 'n_active = 1'
    write(unit,'(A)') 'active_replicates = T'
    close(unit)

    allocate(plan%istates(1),plan%nreps(1),plan%ireps(1,1),plan%nblock(1,1))
    allocate(plan%start_block(1,1,1),plan%end_block(1,1,1))
    allocate(plan%timestep_pattern(n_samples+1,1,1,1))
    plan%istates = 1
    plan%nreps = 1
    plan%ireps = 1
    plan%nblock = 1
    plan%start_block = 1
    plan%end_block = n_samples+1
    plan%timestep_pattern = .true.
    call write_lifting_plan(plan,trim(plan_filename))

    open(newunit=unit,file=trim(pair_filename),status='replace',action='write')
    write(unit,*) .false., .true.
    write(unit,*) .false., .false.
    close(unit)

    trajectory%nobject = 1
    trajectory%ntimes = n_samples+1
    allocate(trajectory%vec(3,2,trajectory%ntimes))
    do frame = 1, trajectory%ntimes
        trajectory%vec(:,1,frame) = matmul(rotation,synthetic_direction(frame,1))
        trajectory%vec(:,2,frame) = matmul(rotation,synthetic_direction(frame,2))
    end do
    call make_filename_lifted_tra(1,1,1,'./'//trim(trajectory_prefix),trajectory_filename)
    call write_trajectory(trim(trajectory_filename),trajectory)

    job%lmax = test_lmax
    job%lag = 1
    job%selector_filename = selector_filename
    job%lifting_plan_filename = plan_filename
    job%pair_active_filename = pair_filename
    job%trajectory_path = '.'
    job%trajectory_prefix = trajectory_prefix
    job%output_path = '.'
    job%output_prefix = trim(prefix)//'_products'
    job%output_policy = 'replace'
    call run_sh_product_sampling_workflow(job)
end subroutine run_sampling

subroutine cleanup_sampling_inputs(prefix)
    character(len=*), intent(in) :: prefix
    character(len=256) :: trajectory_prefix, trajectory_filename

    call delete_file(trim(prefix)//'_selector.txt')
    call delete_file(trim(prefix)//'_lifting_plan.txt')
    call delete_file(trim(prefix)//'_active_pairs.txt')
    trajectory_prefix = './'//trim(prefix)//'_trajectory'
    call make_filename_lifted_tra(1,1,1,trim(trajectory_prefix),trajectory_filename)
    call delete_file(trim(trajectory_filename))
end subroutine cleanup_sampling_inputs

subroutine configure_test_job(job, means_filename, output_prefix)
    type(sh_invariant_contraction_workflow_type), intent(out) :: job
    character(len=*), intent(in) :: means_filename, output_prefix

    job%means_filename = means_filename
    job%output_path = '.'
    job%output_prefix = output_prefix
    job%contraction_mode = 'selected'
    job%n_selected_l_triplets = 2
    allocate(job%selected_l_triplets(3,2))
    job%selected_l_triplets(:,1) = [1,1,1]
    job%selected_l_triplets(:,2) = [1,1,2]
    job%output_policy = 'replace'
end subroutine configure_test_job

function synthetic_direction(sample, role) result(direction)
    integer, intent(in) :: sample, role
    real(real64) :: direction(3)
    real(real64) :: z, phi, theta, pi

    pi = acos(-1.0_real64)
    z = 2.0_real64*modulo(real(17*sample+29*role,real64)*0.6180339887498949_real64, &
        1.0_real64)-1.0_real64
    phi = 2.0_real64*pi*modulo(real(13*sample+7*role,real64)* &
        0.4142135623730950_real64,1.0_real64)
    theta = acos(max(-1.0_real64,min(1.0_real64,z)))
    direction = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
end function synthetic_direction

function identity_matrix() result(matrix)
    real(real64) :: matrix(3,3)
    matrix = 0.0_real64
    matrix(1,1) = 1.0_real64
    matrix(2,2) = 1.0_real64
    matrix(3,3) = 1.0_real64
end function identity_matrix

subroutine make_rotation(matrix)
    real(real64), intent(out) :: matrix(3,3)
    real(real64) :: axis(3), cosine, sine, one_minus_cosine
    integer :: i, j

    axis = [1.0_real64,2.0_real64,3.0_real64]/sqrt(14.0_real64)
    cosine = cos(0.731_real64)
    sine = sin(0.731_real64)
    one_minus_cosine = 1.0_real64-cosine
    do i = 1, 3
        do j = 1, 3
            matrix(i,j) = one_minus_cosine*axis(i)*axis(j)
        end do
        matrix(i,i) = matrix(i,i)+cosine
    end do
    matrix(1,2) = matrix(1,2)-axis(3)*sine
    matrix(2,1) = matrix(2,1)+axis(3)*sine
    matrix(1,3) = matrix(1,3)+axis(2)*sine
    matrix(3,1) = matrix(3,1)-axis(2)*sine
    matrix(2,3) = matrix(2,3)-axis(1)*sine
    matrix(3,2) = matrix(3,2)+axis(1)*sine
end subroutine make_rotation

subroutine read_invariant_value(filename, raw_value, connected_value)
    character(len=*), intent(in) :: filename
    complex(real64), intent(out) :: raw_value, connected_value
    character(len=2048) :: line
    integer :: unit, ios, channel1, channel2
    real(real64) :: raw_re, raw_im, connected_re, connected_im

    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios /= 0) error stop 'Cannot read SH rotation-test output'
    do
        read(unit,'(A)',iostat=ios) line
        if (ios < 0) error stop 'SH rotation-test output has no matrix row'
        if (ios /= 0) error stop 'Error reading SH rotation-test output'
        read(line,*,iostat=ios) channel1, channel2, raw_re, raw_im, connected_re, connected_im
        if (ios /= 0) cycle
        if (channel1 == 1 .and. channel2 == 2) exit
    end do
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing SH rotation-test output'
    raw_value = cmplx(raw_re,raw_im,real64)
    connected_value = cmplx(connected_re,connected_im,real64)
end subroutine read_invariant_value

subroutine assert_invariant(original, rotated, name)
    complex(real64), intent(in) :: original, rotated
    character(len=*), intent(in) :: name
    real(real64) :: difference, allowed_error

    difference = abs(original-rotated)
    allowed_error = tolerance*max(1.0_real64,abs(original),abs(rotated))
    if (difference > allowed_error) then
        write(output_unit,'(A,A,A,ES16.8,A,ES16.8)') 'Rotation changed ',trim(name), &
            ' invariant by ',difference,' (allowed ',allowed_error,')'
        error stop 'SH product rotation-invariance test failed'
    end if
end subroutine assert_invariant

subroutine delete_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios
    logical :: exists

    inquire(file=trim(filename),exist=exists)
    if (.not. exists) return
    open(newunit=unit,file=trim(filename),status='old',iostat=ios)
    if (ios /= 0) error stop 'Cannot open SH rotation-test temporary file for cleanup'
    close(unit,status='delete',iostat=ios)
    if (ios /= 0) error stop 'Cannot delete SH rotation-test temporary file'
end subroutine delete_file

end module sh_product_rotation_invariance_test_module
