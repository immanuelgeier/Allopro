module vmf_module

use iso_fortran_env, only: real64

implicit none

contains


function vmf_norm(kappa) result(c)

real(real64),intent(in)::kappa
real(real64)::c

c = kappa/(4.0_real64*pi*sinh(kappa))

end function



function vmf_density(z,mu,kappa) result(rho)

real(real64),intent(in)::z(3)
real(real64),intent(in)::mu(3)
real(real64),intent(in)::kappa

real(real64)::rho
real(real64)::c

c=vmf_norm(kappa)

rho=c*exp(kappa*dot_product(z,mu))


end function



function vmf_entropy_exact(kappa) result(h)

real(real64),intent(in)::kappa

real(real64)::h
real(real64)::A,c

c=vmf_norm(kappa)

A=1.0_real64/tanh(kappa)-1.0_real64/kappa


h=-log(c)-kappa*A


end function




subroutine integrate_pixel_vmf( &
        vertices, &
        nvert, &
        mu, &
        kappa, &
        prob)

real(real64),intent(in)::vertices(3,nvert)
integer,intent(in)::nvert

real(real64),intent(in)::mu(3)
real(real64),intent(in)::kappa

real(real64),intent(out)::prob

integer::i
integer::nq

real(real64)::area
real(real64)::z(3)
real(real64)::rho_sum


!
! placeholder:
! replace by exact HEALPix polygon area
!

area=spherical_polygon_area(vertices,nvert)



!
! simple quadrature:
! use triangle fan around center
!

rho_sum=0.0_real64

do i=2,nvert-1


    z = vertices(:,1)+ &
        vertices(:,i)+ &
        vertices(:,i+1)

    z=z/sqrt(dot_product(z,z))


    rho_sum=rho_sum+vmf_density(z,mu,kappa)


end do


nq=nvert-2

prob=area*rho_sum/real(nq,real64)



end subroutine




subroutine entropy_healpix(prob,omega,h)

real(real64),intent(in)::prob(:)
real(real64),intent(in)::omega(:)

real(real64),intent(out)::h


integer::i


h=0.0_real64


do i=1,size(prob)

    if(prob(i)>0.0_real64)then

        h=h-prob(i)* &
        log(prob(i)/omega(i))

    endif

enddo


end subroutine



end module