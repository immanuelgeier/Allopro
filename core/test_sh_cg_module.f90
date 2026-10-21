module test_sh_cg_module
use iso_fortran_env, only: real64
use sh_lifting_module

type(sh_lift_type) :: sh

real(real64), allocatable :: r1(:,:), r2(:,:), r3(:,:)
real(real64), allocatable :: Y1_re(:,:), Y1_im(:,:)
real(real64), allocatable :: Y2_re(:,:), Y2_im(:,:)
real(real64), allocatable :: Y3_re(:,:), Y3_im(:,:)

integer :: nframe
integer :: lmax
integer :: nlm

lmax = 20
nframe = 10000

call sh_lift_initialize(sh, lmax)

nlm = sh%nlm

allocate(r1(3,nframe))
allocate(r2(3,nframe))
allocate(r3(3,nframe))

allocate(Y1_re(nlm,nframe))
allocate(Y1_im(nlm,nframe))
allocate(Y2_re(nlm,nframe))
allocate(Y2_im(nlm,nframe))
allocate(Y3_re(nlm,nframe))
allocate(Y3_im(nlm,nframe))

! r1, r2, r3 befüllen ...

call sh_lift_block( &
    sh, r1, r2, r3, &
    Y1_re, Y1_im, &
    Y2_re, Y2_im, &
    Y3_re, Y3_im)

call sh_lift_finalize(sh)

end module test_sh_cg_module