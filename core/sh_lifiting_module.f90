module sh_lifting_module

    use iso_fortran_env, only: real64, int32

    implicit none

    private

    public :: sh_lift_type
    public :: sh_lift_initialize
    public :: sh_lift_finalize
    public :: sh_lift_block
    public :: sh_lm_index

    ! ============================================================
    ! SH lifting object
    ! ============================================================

    type :: sh_lift_type

        integer(int32) :: lmax = -1
        integer(int32) :: nlm  = 0

        ! Normalization constants
        !
        ! norm(l,m) for m >= 0
        !
        real(real64), allocatable :: norm(:,:)

    end type sh_lift_type


contains


    ! ============================================================
    ! Initialization
    ! ============================================================

    subroutine sh_lift_initialize(sh, lmax)

        type(sh_lift_type), intent(out) :: sh
        integer, intent(in)              :: lmax

        integer :: l, m
        real(real64) :: pi

        if (lmax < 0) then
            error stop "SH: lmax must be >= 0"
        end if

        sh%lmax = lmax
        sh%nlm  = (lmax + 1)**2

        allocate(sh%norm(0:lmax,0:lmax))

        sh%norm = 0.0_real64

        pi = acos(-1.0_real64)

        ! --------------------------------------------------------
        ! Condon-Shortley convention:
        !
        ! Y_l^m(theta,phi)
        !
        ! = N_lm P_l^m(cos(theta)) exp(i m phi)
        !
        ! with
        !
        ! N_lm =
        ! sqrt[(2l+1)/(4 pi) * (l-m)!/(l+m)!]
        !
        ! The Condon-Shortley phase is contained in P_l^m.
        ! --------------------------------------------------------

        do l = 0, lmax
            do m = 0, l

                sh%norm(l,m) = sqrt( &
                    real(2*l + 1, real64) / (4.0_real64*pi) * &
                    factorial_ratio(l,m) )

            end do
        end do

    end subroutine sh_lift_initialize


    ! ============================================================
    ! Factorial ratio
    !
    ! (l-m)! / (l+m)!
    !
    ! Computed recursively to avoid factorial overflow.
    ! ============================================================

    pure function factorial_ratio(l,m) result(r)

        integer, intent(in) :: l, m
        real(real64) :: r

        integer :: k

        r = 1.0_real64

        if (m == 0) return

        do k = l-m+1, l+m
            r = r / real(k, real64)
        end do

    end function factorial_ratio


    ! ============================================================
    ! SH index
    !
    ! Maps (l,m) to a flat Fortran array index.
    !
    ! For each l:
    !
    !   m = -l ... +l
    !
    ! Index:
    !
    !   idx = l^2 + l + m + 1
    !
    ! Hence
    !
    !   l=0: index 1
    !   l=1: indices 2..4
    !   l=2: indices 5..9
    !   ...
    ! ============================================================

    pure function sh_lm_index(l,m) result(idx)

        integer, intent(in) :: l, m
        integer :: idx

        idx = l*l + l + m + 1

    end function sh_lm_index


    ! ============================================================
    ! Lift one block of trajectories
    !
    ! Input:
    !
    !   r1(3,nframe)
    !   r2(3,nframe)
    !   r3(3,nframe)
    !
    ! Output:
    !
    !   Y1_re(nlm,nframe)
    !   Y1_im(nlm,nframe)
    !
    !   Y2_re(...)
    !   Y2_im(...)
    !
    !   Y3_re(...)
    !   Y3_im(...)
    !
    ! All vectors are internally normalized.
    !
    ! The output contains
    !
    !   Y_l^m
    !
    ! for
    !
    !   0 <= l <= lmax
    !   -l <= m <= l
    !
    ! ============================================================

    subroutine sh_lift_block( &
        sh, &
        r1, r2, r3, &
        Y1_re, Y1_im, &
        Y2_re, Y2_im, &
        Y3_re, Y3_im)

        type(sh_lift_type), intent(in) :: sh

        real(real64), intent(in) :: r1(:,:)
        real(real64), intent(in) :: r2(:,:)
        real(real64), intent(in) :: r3(:,:)

        real(real64), intent(out) :: Y1_re(:,:)
        real(real64), intent(out) :: Y1_im(:,:)

        real(real64), intent(out) :: Y2_re(:,:)
        real(real64), intent(out) :: Y2_im(:,:)

        real(real64), intent(out) :: Y3_re(:,:)
        real(real64), intent(out) :: Y3_im(:,:)

        integer :: nframe
        integer :: t

        ! --------------------------------------------------------
        ! Basic checks
        ! --------------------------------------------------------

        if (size(r1,1) /= 3 .or. &
            size(r2,1) /= 3 .or. &
            size(r3,1) /= 3) then

            error stop "SH: trajectory must have shape (3,nframe)"

        end if

        nframe = size(r1,2)

        if (size(r2,2) /= nframe .or. &
            size(r3,2) /= nframe) then

            error stop "SH: trajectory frame counts differ"

        end if

        if (size(Y1_re,1) /= sh%nlm .or. &
            size(Y1_im,1) /= sh%nlm .or. &
            size(Y2_re,1) /= sh%nlm .or. &
            size(Y2_im,1) /= sh%nlm .or. &
            size(Y3_re,1) /= sh%nlm .or. &
            size(Y3_im,1) /= sh%nlm) then

            error stop "SH: wrong output dimensions"

        end if

        if (size(Y1_re,2) /= nframe .or. &
            size(Y2_re,2) /= nframe .or. &
            size(Y3_re,2) /= nframe) then

            error stop "SH: wrong output frame dimension"

        end if

        ! --------------------------------------------------------
        ! Parallelize over trajectory
        ! --------------------------------------------------------

        !$omp parallel do schedule(static)
        do t = 1, nframe

            call sh_lift_single( &
                sh, r1(:,t), &
                Y1_re(:,t), Y1_im(:,t))

            call sh_lift_single( &
                sh, r2(:,t), &
                Y2_re(:,t), Y2_im(:,t))

            call sh_lift_single( &
                sh, r3(:,t), &
                Y3_re(:,t), Y3_im(:,t))

        end do
        !$omp end parallel do

    end subroutine sh_lift_block


    ! ============================================================
    ! Lift one vector
    ! ============================================================

    subroutine sh_lift_single( &
        sh, r, &
        Y_re, Y_im)

        type(sh_lift_type), intent(in) :: sh

        real(real64), intent(in)  :: r(3)

        real(real64), intent(out) :: Y_re(:)
        real(real64), intent(out) :: Y_im(:)

        integer :: l, m
        integer :: idx

        real(real64) :: x, y, z
        real(real64) :: rho
        real(real64) :: rnorm

        real(real64) :: cos_theta

        real(real64) :: pmm
        real(real64) :: pmmp1
        real(real64) :: pll

        real(real64) :: cphi
        real(real64) :: sphi

        real(real64) :: cos_mphi
        real(real64) :: sin_mphi

        real(real64) :: phase_re
        real(real64) :: phase_im

        real(real64) :: factor
        real(real64) :: p_prev
        real(real64) :: p_curr

        ! --------------------------------------------------------
        ! Cartesian coordinates
        ! --------------------------------------------------------

        x = r(1)
        y = r(2)
        z = r(3)

        rnorm = sqrt(x*x + y*y + z*z)

        if (rnorm <= tiny(1.0_real64)) then
            error stop "SH: zero-length vector encountered"
        end if

        x = x / rnorm
        y = y / rnorm
        z = z / rnorm

        ! --------------------------------------------------------
        ! cos(theta)
        ! --------------------------------------------------------

        cos_theta = max(-1.0_real64, &
                        min(1.0_real64, z))

        ! --------------------------------------------------------
        ! exp(i phi)
        !
        ! For rho = sqrt(x^2+y^2)
        !
        ! exp(i phi) = (x + i y)/rho
        !
        ! At the poles phi is irrelevant.
        ! We choose phi = 0 there.
        ! --------------------------------------------------------

        rho = sqrt(x*x + y*y)

        if (rho > 100.0_real64*tiny(1.0_real64)) then

            cphi = x / rho
            sphi = y / rho

        else

            cphi = 1.0_real64
            sphi = 0.0_real64

        end if

        ! --------------------------------------------------------
        ! Clear output
        ! --------------------------------------------------------

        Y_re = 0.0_real64
        Y_im = 0.0_real64

        ! ========================================================
        ! m = 0
        !
        ! P_l^0
        ! ========================================================

        p_prev = 1.0_real64

        idx = sh_lm_index(0,0)

        Y_re(idx) = sh%norm(0,0)
        Y_im(idx) = 0.0_real64

        do l = 1, sh%lmax

            if (l == 1) then

                p_curr = cos_theta

            else

                pll = ( &
                    real(2*l-1,real64)*cos_theta*p_curr - &
                    real(l-1,real64)*p_prev ) / &
                    real(l,real64)

                p_prev = p_curr
                p_curr = pll

            end if

            idx = sh_lm_index(l,0)

            Y_re(idx) = sh%norm(l,0) * p_curr
            Y_im(idx) = 0.0_real64

        end do


        ! ========================================================
        ! m > 0
        ! ========================================================

        ! exp(i*0*phi)
        cos_mphi = 1.0_real64
        sin_mphi = 0.0_real64

        ! --------------------------------------------------------
        ! P_m^m starting from P_0^0
        !
        ! P_m^m(x)
        ! =
        ! (-1)^m (2m-1)!! (1-x^2)^(m/2)
        !
        ! We generate it recursively.
        ! --------------------------------------------------------

        pmm = 1.0_real64

        do m = 1, sh%lmax

            pmm = -real(2*m-1,real64) * rho * pmm

            ! ----------------------------------------------------
            ! Update exp(i m phi)
            ! ----------------------------------------------------

            phase_re = cos_mphi*cphi - sin_mphi*sphi
            phase_im = sin_mphi*cphi + cos_mphi*sphi

            cos_mphi = phase_re
            sin_mphi = phase_im

            ! ----------------------------------------------------
            ! l = m
            ! ----------------------------------------------------

            idx = sh_lm_index(m,m)

            Y_re(idx) = sh%norm(m,m) * pmm * cos_mphi
            Y_im(idx) = sh%norm(m,m) * pmm * sin_mphi

            ! ----------------------------------------------------
            ! Negative m:
            !
            ! Y_l^{-m}
            ! =
            ! (-1)^m conj(Y_l^m)
            !
            ! ----------------------------------------------------

            idx = sh_lm_index(m,-m)

            if (mod(m,2) == 0) then

                Y_re(idx) =  Y_re(sh_lm_index(m,m))
                Y_im(idx) = -Y_im(sh_lm_index(m,m))

            else

                Y_re(idx) = -Y_re(sh_lm_index(m,m))
                Y_im(idx) =  Y_im(sh_lm_index(m,m))

            end if

            ! ----------------------------------------------------
            ! l > m
            !
            ! Associated Legendre recurrence:
            !
            ! P_l^m(x) =
            ! [ (2l-1)x P_{l-1}^m
            !   - (l+m-1)P_{l-2}^m ]/(l-m)
            ! ----------------------------------------------------

            if (m < sh%lmax) then

                pmmp1 = cos_theta * &
                        real(2*m+1,real64) * pmm

                ! l = m+1

                idx = sh_lm_index(m+1,m)

                Y_re(idx) = sh%norm(m+1,m) * pmmp1 * cos_mphi
                Y_im(idx) = sh%norm(m+1,m) * pmmp1 * sin_mphi

                ! negative m

                idx = sh_lm_index(m+1,-m)

                if (mod(m,2) == 0) then

                    Y_re(idx) =  Y_re(sh_lm_index(m+1,m))
                    Y_im(idx) = -Y_im(sh_lm_index(m+1,m))

                else

                    Y_re(idx) = -Y_re(sh_lm_index(m+1,m))
                    Y_im(idx) =  Y_im(sh_lm_index(m+1,m))

                end if

                p_prev = pmm
                p_curr = pmmp1

                do l = m+2, sh%lmax

                    pll = ( &
                        real(2*l-1,real64)*cos_theta*p_curr - &
                        real(l+m-1,real64)*p_prev ) / &
                        real(l-m,real64)

                    idx = sh_lm_index(l,m)

                    Y_re(idx) = sh%norm(l,m) * pll * cos_mphi
                    Y_im(idx) = sh%norm(l,m) * pll * sin_mphi

                    ! negative m

                    idx = sh_lm_index(l,-m)

                    if (mod(m,2) == 0) then

                        Y_re(idx) =  Y_re(sh_lm_index(l,m))
                        Y_im(idx) = -Y_im(sh_lm_index(l,m))

                    else

                        Y_re(idx) = -Y_re(sh_lm_index(l,m))
                        Y_im(idx) =  Y_im(sh_lm_index(l,m))

                    end if

                    p_prev = p_curr
                    p_curr = pll

                end do

            end if

        end do

    end subroutine sh_lift_single


    ! ============================================================
    ! Finalize
    ! ============================================================

    subroutine sh_lift_finalize(sh)

        type(sh_lift_type), intent(inout) :: sh

        if (allocated(sh%norm)) deallocate(sh%norm)

        sh%lmax = -1
        sh%nlm  = 0

    end subroutine sh_lift_finalize


end module sh_lifting_module