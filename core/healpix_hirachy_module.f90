module healpix_transfer_entropie_module

use healpix_types
use pix_tools
use tra_lift_module, only: lifted_trajectory,normalize_and_rotate_vec

contains
! 
!type lifted_trajectory
!    integer :: nobject
!    integer :: ntimes
!
!    ! b1/b2 Vektoren
!    ! dimension:
!    ! (3,nobject,ntstep,2)
!    ! letzter Index:
!    ! 1 = b1
!    ! 2 = b2
!    real(8), allocatable :: vec(:,:,:,:)
!
!end type

type :: bin_occupancy

logical,allocatable :: bocc(:,:,:,:,:)
!bocc(npix,npix,npix,nobj,nobj)
end type

type :: healpix_pair_trajectory_block

integer,allocatable :: pix(:,:,:)
!pix(npix(nside),nobj,2*npair)

end type






end module healpix_transfer_entropie_module
