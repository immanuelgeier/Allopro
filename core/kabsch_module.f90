module kabsch_module

contains
subroutine kabsch_rotation(n, ref, xyz, rot, rmsd)
implicit none

integer, intent(in) :: n
real(8), intent(in)  :: ref(3,n)
real(8), intent(in)  :: xyz(3,n)
real(8), intent(out) :: rot(3,3)
real(8), intent(out), optional :: rmsd

real(8) :: cref(3), cxyz(3)
real(8) :: X(3,n), Y(3,n)
real(8) :: H(3,3)
real(8) :: U(3,3), VT(3,3), V(3,3)
real(8) :: S(3)
real(8) :: work(100)
real(8) :: D(3,3)
real(8) :: tmp(3)
real(8) :: err
integer :: i, info, lwork
integer :: sign_det

! ------------------------------------------------------------
! CENTROIDS
! ------------------------------------------------------------

cref = sum(ref,2) / dble(n)
cxyz = sum(xyz,2) / dble(n)

! ------------------------------------------------------------
! CENTERED COORDINATES
!
! X = current structure
! Y = reference structure
! ------------------------------------------------------------

do i=1,n
    X(:,i) = xyz(:,i) - cxyz
    Y(:,i) = ref(:,i) - cref
enddo

! ------------------------------------------------------------
! COVARIANCE MATRIX
!
! H = X Y^T
!
! We want:
!
!       rot X ~= Y
!
! ------------------------------------------------------------

H = matmul(X, transpose(Y))

! ------------------------------------------------------------
! SVD
!
! H = U S VT
! ------------------------------------------------------------

lwork = 100

call dgesvd('A','A',3,3,H,3,S,U,3,VT,3,work,lwork,info)

if (info /= 0) then
    write(0,*) 'KABSCH: DGESVD failed, info=',info
    stop
endif

V = transpose(VT)

! ------------------------------------------------------------
! Ensure a PROPER rotation:
!
! det(R) = +1
!
! R = V D U^T
! ------------------------------------------------------------

D = 0.d0
D(1,1) = 1.d0
D(2,2) = 1.d0
D(3,3) = 1.d0

if (det3(V) * det3(U) < 0.d0) then
    D(3,3) = -1.d0
endif

rot = matmul(V,matmul(D,transpose(U)))

! ------------------------------------------------------------
! OPTIONAL RMSD
! ------------------------------------------------------------

if (present(rmsd)) then

    err = 0.d0

    do i=1,n
        tmp = matmul(rot,X(:,i)) - Y(:,i)
        err = err + dot_product(tmp,tmp)
    enddo

    rmsd = sqrt(err/dble(n))

endif

end subroutine kabsch_rotation


real(8) function det3(A)
    implicit none

    real(8), intent(in) :: A(3,3)

    det3 = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
         - A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
         + A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))

end function det3

end module