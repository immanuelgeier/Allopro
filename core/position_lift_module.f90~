module position_lift_module
    use time_window_file_mapping
    use time_window_loader
    use direction_vectors_analyse_plan_module

    implicit none
    private
    public :: lifted_trajectory, lift_cv_plan

    type :: lifted_trajectory
        integer :: nobject
        integer :: ntimes
        ! Dimensionen: (3, 2 * nobject, ntimes)
        ! Die 2. Dimension liegt abwechselnd im Speicher:
        ! 2*obj - 1  -> b1 (comp=1)
        ! 2*obj      -> b2 (comp=2)
        real(8), allocatable :: vec(:,:,:)
        !real(8), allocatable :: rot_matrix(:,:,:)
    end type lifted_trajectory

contains

    subroutine lift_cv_plan(plan, block, lifted)
        type(cv_plan), intent(in)              :: plan
        type(raw_data_block), intent(in)       :: block
        type(lifted_trajectory), intent(out)   :: lifted

        integer :: nobject
        integer :: t, k, obj, comp, atom, target_idx

        nobject = size(plan%object)
        lifted%ntimes = size(block%xyz, 3)
        lifted%nobject = nobject

        ! 2. Dimension verdoppelt für das verschachtelte (interleaved) Layout
        allocate(lifted%vec(3, 2 * nobject, lifted%ntimes))
        lifted%vec = 0.d0

        ! ----------------------------------------
        ! Frames & CV Akkumulation
        ! ----------------------------------------
        do t = 1, lifted%ntimes
            do k = 1, size(plan%term)
                atom = plan%term(k)%atom
                obj  = plan%term(k)%target
                comp = plan%term(k)%comp ! 1 für b1, 2 für b2

                ! Mathematisch exakte Verschachtelung:
                ! Wenn comp = 1 (b1) -> 2 * obj - 1
                ! Wenn comp = 2 (b2) -> 2 * obj
                target_idx = 2 * obj - (2 - comp)

                lifted%vec(:, target_idx, t) = &
                    lifted%vec(:, target_idx, t) + &
                    block%xyz(:, atom, t) * plan%term(k)%coeff
            end do
        end do
        
        !! Rotationsmatrix-Speicher vorbereiten (3x3 pro Zeitschritt)
        !allocate(lifted%rot_matrix(3, 3, lifted%ntimes))
        !lifted%rot_matrix = 0.d0

    end subroutine lift_cv_plan

end module position_lift_module
