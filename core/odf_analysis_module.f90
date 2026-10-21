module odf_analysis_module

use healpix_types
use pix_tools

implicit none

contains

!===========================================================
! Count local maxima in a HEALPix histogram
!===========================================================

subroutine find_local_maxima(hist,nside,nmax,max_pix,max_value)

implicit none

integer,intent(in)  :: hist(:)
integer,intent(in)  :: nside

integer,intent(out) :: nmax
integer,intent(out) :: max_pix(:)
integer,intent(out) :: max_value(:)

integer :: npix
integer :: ipix
integer :: nneigh
integer :: neigh(8)
integer :: i

logical :: is_max


npix=size(hist)

nmax=0

do ipix=0,npix-1

    if(hist(ipix+1)==0) cycle

    call neighbours_nest(nside,ipix,neigh,nneigh)

    is_max=.true.

    do i=1,nneigh

        if(neigh(i)<0) cycle

        if(hist(neigh(i)+1) > hist(ipix+1)) then
            is_max=.false.
            exit
        endif

    enddo

    if(is_max) then

        nmax=nmax+1

        if(nmax<=size(max_pix)) then
            max_pix(nmax)=ipix
            max_value(nmax)=hist(ipix+1)
        endif

    endif

enddo

end subroutine find_local_maxima


!===========================================================
! Sort maxima by height (descending)
!===========================================================

subroutine sort_maxima(nmax,max_pix,max_value)

implicit none

integer,intent(in)    :: nmax
integer,intent(inout) :: max_pix(:)
integer,intent(inout) :: max_value(:)

integer :: i,j
integer :: tmp

do i=1,nmax-1
do j=i+1,nmax

    if(max_value(j)>max_value(i)) then

        tmp=max_value(i)
        max_value(i)=max_value(j)
        max_value(j)=tmp

        tmp=max_pix(i)
        max_pix(i)=max_pix(j)
        max_pix(j)=tmp

    endif

enddo
enddo

end subroutine sort_maxima


!===========================================================
! Return only the N largest maxima
!===========================================================

subroutine largest_maxima(hist,nside,nkeep,nfound,pix,value)

implicit none

integer,intent(in)  :: hist(:)
integer,intent(in)  :: nside
integer,intent(in)  :: nkeep

integer,intent(out) :: nfound
integer,intent(out) :: pix(nkeep)
integer,intent(out) :: value(nkeep)

integer :: npix

integer,allocatable :: all_pix(:)
integer,allocatable :: all_value(:)


npix=size(hist)

allocate(all_pix(npix))
allocate(all_value(npix))

call find_local_maxima( &
    hist, &
    nside, &
    nfound, &
    all_pix, &
    all_value)

call sort_maxima( &
    nfound, &
    all_pix, &
    all_value)

nfound=min(nfound,nkeep)

pix(1:nfound)=all_pix(1:nfound)
value(1:nfound)=all_value(1:nfound)

if(nfound<nkeep) then
    pix(nfound+1:nkeep)=0
    value(nfound+1:nkeep)=0
endif

deallocate(all_pix)
deallocate(all_value)

end subroutine largest_maxima


subroutine find_hist_basins(hist,nside, &
                            basin, &
                            maxima, &
                            basin_counts, &
                            nmax)

implicit none

integer, intent(in) :: hist(:)
integer, intent(in) :: nside

integer, allocatable, intent(out) :: basin(:)
integer, allocatable, intent(out) :: maxima(:)
integer, allocatable, intent(out) :: basin_counts(:)

integer, intent(out) :: nmax

integer :: npix
integer :: i,j,k
integer :: maxpix
integer :: next
integer :: current
integer :: nneigh
integer :: neighbors(8)

integer, allocatable :: flow(:)
integer, allocatable :: root(:)


npix=size(hist)


if(npix /= 12*nside*nside) then
    write(0,*) "ERROR: inconsistent HEALPix grid"
    write(0,*) "npix =",npix
    write(0,*) "nside=",nside
    stop
endif


allocate(flow(npix))
allocate(root(npix))
allocate(basin(npix))

flow=0
root=0
basin=0


!
! --------------------------------------------------
! 1) Gradient flow
! --------------------------------------------------
!
! Each occupied pixel points to its highest neighbour.
!

do i=1,npix


    !
    ! Empty pixels are ignored
    !
    if(hist(i)==0) then
        flow(i)=i
        cycle
    endif


    maxpix=i


    neighbors=-1

    call neighbours_nest( &
        nside, &
        i-1, &
        neighbors, &
        nneigh )


    do k=1,nneigh


        j=neighbors(k)+1


        if(j < 1 .or. j > npix) then
            write(0,*) "bad neighbour"
            write(0,*) "pixel=",i
            write(0,*) "neighbor=",neighbors(k)
            stop
        endif


        !
        ! move uphill
        !
        if(hist(j)>hist(maxpix)) then

            maxpix=j

        elseif(hist(j)==hist(maxpix)) then

            !
            ! deterministic plateau handling
            !
            if(j < maxpix) maxpix=j

        endif


    enddo


    flow(i)=maxpix


enddo



!
! --------------------------------------------------
! 2) Follow gradient to maxima
! --------------------------------------------------
!

do i=1,npix


    current=i


    do

        next=flow(current)


        if(next==current) exit


        current=next


    enddo


    root(i)=current


enddo



!
! --------------------------------------------------
! 3) Collect maxima
! --------------------------------------------------
!

allocate(maxima(npix))

nmax=0


do i=1,npix


    if(hist(i)>0 .and. root(i)==i) then

        nmax=nmax+1

        maxima(nmax)=i

    endif


enddo


if(nmax==0) then

    write(0,*) "WARNING: no maxima found"

    allocate(basin_counts(0))

    deallocate(flow)
    deallocate(root)

    return

endif


call resize_integer_array(maxima,nmax)



!
! --------------------------------------------------
! 4) Assign basin IDs
! --------------------------------------------------
!

do i=1,npix


    if(hist(i)==0) cycle


    do j=1,nmax


        if(root(i)==maxima(j)) then

            basin(i)=j

            exit

        endif


    enddo


enddo



!
! check
!

do i=1,npix

    if(hist(i)>0 .and. basin(i)==0) then

        write(0,*) "ERROR: pixel without basin"
        write(0,*) "pixel=",i
        write(0,*) "hist=",hist(i)
        write(0,*) "root=",root(i)

        stop

    endif

enddo



!
! --------------------------------------------------
! 5) Sum histogram counts per basin
! --------------------------------------------------
!

allocate(basin_counts(nmax))

basin_counts=0


do i=1,npix

    if(basin(i)>0) then

        basin_counts(basin(i)) = &
            basin_counts(basin(i)) + hist(i)

    endif

enddo

maxima=maxima-1

deallocate(flow)
deallocate(root)


end subroutine find_hist_basins

subroutine resize_integer_array(a,n)

implicit none

integer,allocatable,intent(inout) :: a(:)
integer,intent(in) :: n

integer,allocatable :: tmp(:)


allocate(tmp(n))

tmp=a(1:n)

deallocate(a)

allocate(a(n))

a=tmp

deallocate(tmp)


end subroutine resize_integer_array




end module odf_analysis_module