module output_policy_module

    implicit none

    integer, parameter :: OUTPUT_WRITE       = 1
    integer, parameter :: OUTPUT_SKIP        = 2
    integer, parameter :: OUTPUT_COMPUTE_ONLY = 3
    integer, parameter :: OUTPUT_ABORT       = 4

    ! ask_state:
    ! 0 = normal fragen
    ! 1 = alle vorhandenen Outputs überspringen
    ! 2 = abort
    integer, save :: ask_state = 0

contains


    ! =========================================================
    ! Construct effective output prefix
    !
    ! replace / ask -> unchanged
    ! new            -> append "_new"
    ! =========================================================

    function get_output_prefix(prefix, policy) result(effective_prefix)

        implicit none

        character(len=*), intent(in) :: prefix
        character(len=*), intent(in) :: policy

        character(len=2048) :: effective_prefix

        select case (trim(policy))

        case ('replace', 'ask')

            effective_prefix = trim(prefix)

        case ('new')

            effective_prefix = trim(prefix)//'_new'
            write(0,*) 'changed prefix to:',trim(prefix)//'_new'

        case default

            write(0,*) 'ERROR: Unknown output policy: ', trim(policy)
            write(0,*) 'Valid policies: replace, ask, new'
            error stop

        end select

    end function get_output_prefix


    ! =========================================================
    ! Decide what to do BEFORE calculation
    !
    ! OUTPUT_WRITE:
    !     calculate + write
    !
    ! OUTPUT_SKIP:
    !     don't calculate
    !
    ! OUTPUT_COMPUTE_ONLY:
    !     calculate, don't write
    !
    ! OUTPUT_ABORT:
    !     abort workflow
    ! =========================================================

    subroutine check_output_policy( &
        output_file, output_policy, action)

        implicit none

        character(len=*), intent(in) :: output_file
        character(len=*), intent(in) :: output_policy

        integer, intent(out) :: action

        logical :: file_exists
        character(len=10) :: answer

        action = OUTPUT_WRITE

        inquire( &
            file   = trim(output_file), &
            exist  = file_exists )

        ! -----------------------------------------------------
        ! No existing output -> just calculate and write
        ! -----------------------------------------------------

        if (.not. file_exists) return


        ! -----------------------------------------------------
        ! Apply policy
        ! -----------------------------------------------------

        select case (trim(output_policy))


        case ('replace')

            write(0,*) 'Existing output will be replaced:'
            write(0,*) trim(output_file)

            action = OUTPUT_WRITE


        case ('new')

            ! This means that *_new already exists.
            ! We do NOT silently overwrite it.

            write(0,*) 'Output file already exists:'
            write(0,*) trim(output_file)
            write(0,*) 'Policy "new" cannot overwrite it.'

            action = OUTPUT_SKIP


        case ('ask')

            ! -------------------------------------------------
            ! Persistent "skip all"
            ! -------------------------------------------------

            if (ask_state == 1) then
                action = OUTPUT_SKIP
                return
            end if

            if (ask_state == 2) then
                action = OUTPUT_ABORT
                return
            end if


            ! -------------------------------------------------
            ! Ask user
            ! -------------------------------------------------

            write(0,*)
            write(0,*) 'Output file already exists:'
            write(0,*) trim(output_file)
            write(0,*)
            write(0,*) '[o] overwrite and calculate'
            write(0,*) '[s] skip calculation'
            write(0,*) '[c] calculate but do not save'
            write(0,*) '[a] skip all remaining'
            write(0,*) '[q] abort workflow'
            write(0,*)

            do

                write(0,'(A)',advance='no') &
                    'Choice [o/s/c/a/q]: '

                read(*,'(A)') answer

                select case (answer(1:1))


                case ('o','O')

                    action = OUTPUT_WRITE
                    exit


                case ('s','S')

                    action = OUTPUT_SKIP
                    exit


                case ('c','C')

                    action = OUTPUT_COMPUTE_ONLY
                    exit


                case ('a','A')

                    ask_state = 1
                    action = OUTPUT_SKIP
                    exit


                case ('q','Q')

                    ask_state = 2
                    action = OUTPUT_ABORT
                    exit


                case default

                    write(0,*) 'Invalid choice.'
                    write(0,*) &
                        'Please enter o, s, c, a or q.'

                end select

            end do


        case default

            write(0,*) 'ERROR: Unknown output policy: ', &
                       trim(output_policy)

            write(0,*) &
                'Valid policies: replace, ask, new'

            action = OUTPUT_ABORT

        end select

    end subroutine check_output_policy


    ! =========================================================
    ! Reset persistent interactive state
    ! =========================================================

    subroutine reset_output_policy()

        implicit none

        ask_state = 0

    end subroutine reset_output_policy


end module output_policy_module