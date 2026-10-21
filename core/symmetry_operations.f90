module symmetry_operations

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    real(dp), parameter :: pi = acos(-1.0_dp)

contains


!=========================================================
! Rotation um beliebige Achse
!
! Rodrigues-Formel:
!
! R = cos(theta) I
!   + (1-cos(theta)) n*n^T
!   + sin(theta) [n]x
!
!=========================================================
function rotation_matrix(axis, theta) result(R)

    real(dp), intent(in) :: axis(3)
    real(dp), intent(in) :: theta

    real(dp) :: R(3,3)

    real(dp) :: n(3)
    real(dp) :: c, s, v


    n = axis / sqrt(sum(axis**2))

    c = cos(theta)
    s = sin(theta)
    v = 1.0_dp - c


    R(1,1) = c + n(1)*n(1)*v
    R(2,2) = c + n(2)*n(2)*v
    R(3,3) = c + n(3)*n(3)*v


    R(1,2) = n(1)*n(2)*v - n(3)*s
    R(2,1) = n(1)*n(2)*v + n(3)*s

    R(1,3) = n(1)*n(3)*v + n(2)*s
    R(3,1) = n(1)*n(3)*v - n(2)*s

    R(2,3) = n(2)*n(3)*v - n(1)*s
    R(3,2) = n(2)*n(3)*v + n(1)*s


end function rotation_matrix



!=========================================================
! Spiegelung an Ebene
!
! Ebene definiert durch Normalenvektor n
!
! sigma = I - 2 n n^T
!
!=========================================================
function reflection_matrix(normal) result(R)

    real(dp), intent(in) :: normal(3)

    real(dp) :: R(3,3)

    real(dp) :: n(3)


    n = normal / sqrt(sum(normal**2))


    R = 0.0_dp

    R(1,1)=1.0_dp
    R(2,2)=1.0_dp
    R(3,3)=1.0_dp


    R = R - 2.0_dp * outer_product(n,n)


end function reflection_matrix



!=========================================================
! Inversion
!
! x -> -x
!
!=========================================================
function inversion_matrix() result(R)

    real(dp) :: R(3,3)

    R = 0.0_dp

    R(1,1)=-1.0_dp
    R(2,2)=-1.0_dp
    R(3,3)=-1.0_dp

end function inversion_matrix



!=========================================================
! kleine Hilfsfunktion
! äußeres Produkt
!=========================================================
function outer_product(a,b) result(M)

    real(dp), intent(in) :: a(3), b(3)

    real(dp) :: M(3,3)

    integer :: i,j


    do i=1,3
        do j=1,3
            M(i,j)=a(i)*b(j)
        end do
    end do

end function outer_product


end module symmetry_operations