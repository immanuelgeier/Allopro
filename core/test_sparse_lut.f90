program test_sparse_lut
implicit none
integer :: nside,npix,no
integer, allocatable :: odf1(:,:), odf1(:,:), odf2(:,:) 

no=10

nside=1
npix=12*(4**(nside-1))

allocate(odf1(npix,no),odf2(npix,no),odf3(npix,no))

contains

subroutine init_sparse_random_odf(ncount,npi,n)
implicit none
integer :: ncount,i,j,k,npi,n
real(8) :: rn
odf1=0._8
odf2=0._8
odf3=0._8

do j=1,n
  do i=1,ncount
    call random_number(rn)
    npi=floor(int(rn*dble(npi))
    odf1(npi,j)=odf1(npi,j)+1
    call random_number(rn)
    npi=floor(int(rn*dble(npi))
    odf2(npi,j)=odf2(npi,j)+1
    call random_number(rn)
    npi=floor(int(rn*dble(npi))
    odf3(npi,j)=odf3(npi,j)+1
  end do
end do


end subroutine init sparse_random_odf


end program test_sparse_lut