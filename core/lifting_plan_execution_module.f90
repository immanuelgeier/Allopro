module lifting_plan_execution_module

use lifting_plan_module
use position_lift_module
use direction_vectors_analyse_plan_module
use project_module
use lifted_trajectory_io_module
use time_window_loader
use refstruc_pbc_module
contains


subroutine execute_lifting_plan( &
    cv_p, &
    lifting_plan, &
    project, &
    ref, &
    path_lifted_tra, &
    prefix_lifted_tra, &
    output_policy)

    use output_policy_module

    implicit none

    type(lifting_plan_type) :: lifting_plan
    type(time_window)      :: window
    type(raw_data_block)   :: block
    type(lifted_trajectory) :: lifted_tra

    type(cv_plan), intent(in) :: cv_p
    type(project_type)       :: project
    type(refstruc_type), intent(in) :: ref

    character(len=*), intent(in) :: path_lifted_tra
    character(len=*), intent(in) :: prefix_lifted_tra
    character(len=*), intent(in) :: output_policy

    integer :: irep, istate, iblock
    integer :: output_action

    logical, allocatable :: mask(:)

    character(len=800) :: filename
    character(len=800) :: prefix
    character(len=800) :: effective_prefix


    ! ------------------------------------------------------------
    ! Construct prefix including path
    ! ------------------------------------------------------------

    prefix = trim(path_lifted_tra)//'/'//trim(prefix_lifted_tra)

    effective_prefix = get_output_prefix( &
        prefix, output_policy)


    ! ------------------------------------------------------------
    ! Loop over lifting plan
    ! ------------------------------------------------------------

    do istate = 1, size(lifting_plan%istates)

        do irep = 1, lifting_plan%nreps(istate)

            do iblock = 1, lifting_plan%nblock(istate,irep)


                ! ------------------------------------------------
                ! Construct output filename FIRST
                ! ------------------------------------------------

                call make_filename_lifted_tra( &
                    istate, &
                    irep, &
                    iblock, &
                    effective_prefix, &
                    filename)


                ! ------------------------------------------------
                ! Check output policy BEFORE loading raw data
                ! ------------------------------------------------

                call check_output_policy( &
                    filename, &
                    output_policy, &
                    output_action)


                select case (output_action)

                case (OUTPUT_SKIP)

                    write(0,*) 'Skipping lifting:'
                    write(0,*) '  ', trim(filename)

                    cycle


                case (OUTPUT_ABORT)

                    write(0,*) 'Aborting lifting workflow.'
                    return


                case (OUTPUT_WRITE)

                    write(0,*) 'Lifting -> ', trim(filename)


                case (OUTPUT_COMPUTE_ONLY)

                    write(0,*) 'Computing without writing:'
                    write(0,*) '  ', trim(filename)


                case default

                    write(0,*) 'ERROR: Invalid output action.'
                    return

                end select


                ! ------------------------------------------------
                ! Load raw data
                ! ------------------------------------------------

                call init_time_window_from_lifting_plan( &
                    lifting_plan, &
                    istate, &
                    irep, &
                    iblock, &
                    window, &
                    mask)

                write(0,*) 'loading block'

                call load_raw_data_block( &
                    project, &
                    window, &
                    mask, &
                    block)


                ! ------------------------------------------------
                ! Unwrap
                ! ------------------------------------------------

                call unwrap_xyz_block( &
                    ref, &
                    block%xyz, &
                    block%boxl)


                ! ------------------------------------------------
                ! Lift CV
                ! ------------------------------------------------

                call lift_cv_plan( &
                    cv_p, &
                    block, &
                    lifted_tra)


                ! ------------------------------------------------
                ! Write result
                ! ------------------------------------------------

                if (output_action == OUTPUT_WRITE) then

                    write(0,*) 'Writing: ', trim(filename)

                    call write_trajectory( &
                        filename, &
                        lifted_tra)

                else if (output_action == OUTPUT_COMPUTE_ONLY) then

                    write(0,*) &
                        'Output skipped: ', trim(filename)

                end if


                if (allocated(lifted_tra%vec)) then
                    deallocate(lifted_tra%vec)
                endif

            end do
        end do
    end do

end subroutine execute_lifting_plan

subroutine execute_lifting_plan_centroid( &
    cv_p, &
    lifting_plan, &
    project, &
    ref, &
    path_lifted_tra, &
    prefix_lifted_tra, &
    path_centroid_tra, &
    prefix_centroid_tra, &
    output_policy)

    use output_policy_module

    implicit none

    type(lifting_plan_type)  :: lifting_plan
    type(time_window)        :: window
    type(raw_data_block)     :: block

    type(lifted_trajectory)  :: lifted_tra
    type(lifted_trajectory)  :: centroid_tra

    type(cv_plan), intent(in) :: cv_p
    type(project_type)        :: project
    type(refstruc_type), intent(in) :: ref

    character(len=*), intent(in) :: path_lifted_tra
    character(len=*), intent(in) :: prefix_lifted_tra

    character(len=*), intent(in) :: path_centroid_tra
    character(len=*), intent(in) :: prefix_centroid_tra

    character(len=*), intent(in) :: output_policy

    integer :: irep, istate, iblock
    integer :: output_action_lifted
    integer :: output_action_centroid

    logical, allocatable :: mask(:)

    character(len=800) :: filename_lifted
    character(len=800) :: filename_centroid

    character(len=800) :: prefix_lifted
    character(len=800) :: prefix_centroid

    character(len=800) :: effective_prefix_lifted
    character(len=800) :: effective_prefix_centroid


    ! ------------------------------------------------------------
    ! Construct prefixes including paths
    ! ------------------------------------------------------------

    prefix_lifted = &
        trim(path_lifted_tra)//'/'//trim(prefix_lifted_tra)

    prefix_centroid = &
        trim(path_centroid_tra)//'/'//trim(prefix_centroid_tra)


    effective_prefix_lifted = &
        get_output_prefix(prefix_lifted, output_policy)

    effective_prefix_centroid = &
        get_output_prefix(prefix_centroid, output_policy)


    ! ------------------------------------------------------------
    ! Loop over lifting plan
    ! ------------------------------------------------------------

    do istate = 1, size(lifting_plan%istates)

        do irep = 1, lifting_plan%nreps(istate)

            do iblock = 1, lifting_plan%nblock(istate,irep)


                ! ------------------------------------------------
                ! Construct BOTH output filenames FIRST
                ! ------------------------------------------------

                call make_filename_lifted_tra( &
                    istate, &
                    irep, &
                    iblock, &
                    effective_prefix_lifted, &
                    filename_lifted)

                call make_filename_lifted_tra( &
                    istate, &
                    irep, &
                    iblock, &
                    effective_prefix_centroid, &
                    filename_centroid)


                ! ------------------------------------------------
                ! Check BOTH output policies BEFORE computation
                ! ------------------------------------------------

                call check_output_policy( &
                    filename_lifted, &
                    output_policy, &
                    output_action_lifted)

                call check_output_policy( &
                    filename_centroid, &
                    output_policy, &
                    output_action_centroid)


                ! ------------------------------------------------
                ! Abort if either output requests abort
                ! ------------------------------------------------

                if (output_action_lifted == OUTPUT_ABORT .or. &
                    output_action_centroid == OUTPUT_ABORT) then

                    write(0,*) 'Aborting lifting workflow.'
                    return

                endif


                ! ------------------------------------------------
                ! If BOTH outputs can be skipped, do not calculate
                ! ------------------------------------------------

                if (output_action_lifted == OUTPUT_SKIP .and. &
                    output_action_centroid == OUTPUT_SKIP) then

                    write(0,*) 'Skipping block:'
                    write(0,*) '  ', trim(filename_lifted)
                    write(0,*) '  ', trim(filename_centroid)

                    cycle

                endif


                ! ------------------------------------------------
                ! Report actions
                ! ------------------------------------------------

                if (output_action_lifted == OUTPUT_WRITE) then

                    write(0,*) 'Lifting -> ', trim(filename_lifted)

                else if (output_action_lifted == OUTPUT_COMPUTE_ONLY) then

                    write(0,*) 'Computing lifted trajectory without writing:'
                    write(0,*) '  ', trim(filename_lifted)

                else if (output_action_lifted == OUTPUT_SKIP) then

                    write(0,*) 'Skipping lifted trajectory:'
                    write(0,*) '  ', trim(filename_lifted)

                endif


                if (output_action_centroid == OUTPUT_WRITE) then

                    write(0,*) 'Centroid -> ', trim(filename_centroid)

                else if (output_action_centroid == OUTPUT_COMPUTE_ONLY) then

                    write(0,*) 'Computing centroid trajectory without writing:'
                    write(0,*) '  ', trim(filename_centroid)

                else if (output_action_centroid == OUTPUT_SKIP) then

                    write(0,*) 'Skipping centroid trajectory:'
                    write(0,*) '  ', trim(filename_centroid)

                endif


                ! ------------------------------------------------
                ! Load raw data
                ! ------------------------------------------------

                call init_time_window_from_lifting_plan( &
                    lifting_plan, &
                    istate, &
                    irep, &
                    iblock, &
                    window, &
                    mask)

                write(0,*) 'loading block'

                call load_raw_data_block( &
                    project, &
                    window, &
                    mask, &
                    block)


                ! ------------------------------------------------
                ! Unwrap
                ! ------------------------------------------------

                call unwrap_xyz_block( &
                    ref, &
                    block%xyz, &
                    block%boxl)


                ! ------------------------------------------------
                ! Calculate weighted trajectory AND centroid
                ! ------------------------------------------------

                call lift_cv_plan_all( &
                    cv_p, &
                    block, &
                    lifted_tra, &
                    centroid_tra)


                ! ------------------------------------------------
                ! Write lifted trajectory
                ! ------------------------------------------------

                if (output_action_lifted == OUTPUT_WRITE) then

                    write(0,*) 'Writing: ', trim(filename_lifted)

                    call write_trajectory( &
                        filename_lifted, &
                        lifted_tra)

                endif


                ! ------------------------------------------------
                ! Write centroid trajectory
                ! ------------------------------------------------

                if (output_action_centroid == OUTPUT_WRITE) then

                    write(0,*) 'Writing: ', trim(filename_centroid)

                    call write_trajectory( &
                        filename_centroid, &
                        centroid_tra)

                endif


                ! ------------------------------------------------
                ! Free memory for next block
                ! ------------------------------------------------

                if (allocated(lifted_tra%vec)) then
                    deallocate(lifted_tra%vec)
                endif

                if (allocated(centroid_tra%vec)) then
                    deallocate(centroid_tra%vec)
                endif

            end do
        end do
    end do

end subroutine execute_lifting_plan_centroid

subroutine make_filename_lifted_tra(istate,irep,iblock,prefix,filename)

   implicit none

   integer, intent(in) :: istate
   integer, intent(in) :: irep
   integer, intent(in) :: iblock

   character(*), intent(in)  :: prefix
   character(*), intent(out) :: filename

   write(filename,'(A,"_state",I4.4,"_rep",I4.4,"_block",I5.5,".bin")') &
        trim(prefix), istate, irep, iblock

end subroutine make_filename_lifted_tra




end module 