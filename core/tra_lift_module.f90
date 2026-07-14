module tra_lift_module

use topology_module
use time_window_file_mapping
use time_window_loader
use sidechains
use domains_module
use healpix_types
use pix_tools
   
implicit none

type shuffled_time_map

    integer, allocatable :: index(:)

end type shuffled_time_map

type healpix_marginals

    ! Zweidimensionale Marginalverteilungen
    !
    ! x = p_i(t+tau)      : zukünftiger Zustand des Zielobjektes i
    ! y = p_i(t)          : vergangener Zustand des Zielobjektes i
    ! z = p_j(t)          : vergangener Zustand des Quellobjektes j
    ! Dimensionen:
    ! xy(npix,npix,nob,nob)
    !   xy(ipix_x,ipix_y,target_component,source_component)
    !   Summe über z:
    !   C(x,y) = sum_z C(x,y,z)
    integer, allocatable :: xy(:,:,:,:)
    ! xz(npix,npix,nob,nob)
    !   Summe über y:
    !   C(x,z) = sum_y C(x,y,z)
    !   Wird z.B. für alternative Informationsmaße benötigt.
    integer, allocatable :: xz(:,:,:,:)
    ! yz(npix,npix,nob,nob)
    !   Summe über x:
    !   C(y,z) = sum_x C(x,y,z)
    !   Dies ist die wichtige Randverteilung für TE:
    !   p(y_t,z_t)
    integer, allocatable :: yz(:,:,:,:)
    ! Eindimensionale Marginalverteilungen
    ! x(npix,nob,nob)
    !   C(x) = sum_y sum_z C(x,y,z)
    !   zukünftiger Zustand:
    !   p_i(t+tau)
    integer, allocatable :: x(:,:,:)
    ! y(npix,nob,nob)
    !   C(y) = sum_x sum_z C(x,y,z)
    !   vergangener Zustand des Zielobjektes:
    !   p_i(t)
    !   Wird in der TE-Formel benötigt.
    integer, allocatable :: y(:,:,:)
    ! z(npix,nob,nob)
    !   C(z) = sum_x sum_y C(x,y,z)
    !   vergangener Zustand des Quellobjektes:
    !   p_j(t)
    integer, allocatable :: z(:,:,:)

end type healpix_marginals


type lagged_healpix_density

integer, allocatable :: raw_counts(:,:,:,:,:)
logical, allocatable :: occ(:,:,:,:,:)

end type

type scalar_block_mean

real(8), allocatable :: mean_length(:,:)
real(8), allocatable :: mean_squared_length(:,:)
real(8), allocatable :: mean_length_cross(:)
end type

type lifted_trajectory
    integer :: nobject
    integer :: ntimes

    ! b1/b2 Vektoren
    ! dimension:
    ! (3,nobject,ntstep,2)
    ! letzter Index:
    ! 1 = b1
    ! 2 = b2
    real(8), allocatable :: vec(:,:,:,:)

end type

type cv_term
    integer :: atom      ! absoluter Atomindex
    integer :: target    ! welche CV?
    integer :: comp      ! 1=b1, 2=b2
    real(8) :: coeff
end type

type cv_object
    integer :: first_term
    integer :: last_term
    integer :: object_type   ! domain / residue
    integer :: object_id
end type

type cv_plan
    type(cv_term), allocatable   :: term(:)
    type(cv_object), allocatable :: object(:)
end type

contains

subroutine generate_shuffle(ntimes,shuffle)

implicit none

integer,intent(in) :: ntimes
type(shuffled_time_map),intent(out) :: shuffle

integer :: i,j,tmp
real(8) :: r


allocate(shuffle%index(ntimes))


do i=1,ntimes
    shuffle%index(i)=i
end do


! Fisher-Yates shuffle

do i=ntimes,2,-1

    call random_number(r)

    j=1+int(r*i)

    tmp=shuffle%index(i)
    shuffle%index(i)=shuffle%index(j)
    shuffle%index(j)=tmp

end do


end subroutine generate_shuffle


subroutine compute_transfer_entropy(nside,hdens,marg,TE)

implicit none

integer,intent(in) :: nside


type(lagged_healpix_density),intent(in) :: hdens
type(healpix_marginals),intent(in) :: marg

real(8),allocatable,intent(out) :: TE(:,:)

integer :: npix,nob
integer :: i,j,ix,iy,iz
integer :: a,b,c
integer :: total
real(8) :: pxyz,pxy,pyz,py
real(8) :: logterm


npix=size(hdens%raw_counts,1)
nob= size(hdens%raw_counts,4)


allocate(TE(nob,nob))

TE=0.0d0


do i=1,nob
do j=1,nob


   total=sum(hdens%raw_counts(:,:,:,i,j))


   if(total==0) cycle


   do ix=0,npix-1
   do iy=0,npix-1
   do iz=0,npix-1


      a=hdens%raw_counts(ix,iy,iz,i,j)


      if(a==0) cycle


      ! Wahrscheinlichkeiten

      pxyz = real(a,8)/real(total,8)

      pxy = real(marg%xy(ix,iy,i,j),8)/real(total,8)

      pyz = real(marg%yz(iy,iz,i,j),8)/real(total,8)

      py  = real(marg%y(iy,i,j),8)/real(total,8)


      if(pxy>0.d0 .and. pyz>0.d0 .and. py>0.d0) then

         logterm = log( pxyz*py/(pxy*pyz) )

         TE(i,j)=TE(i,j)+pxyz*logterm

      end if


   end do
   end do
   end do


end do
end do

end subroutine compute_transfer_entropy



subroutine marginalize_healpix_density(nside,lifted,hdens,marg)

implicit none

integer,intent(in) :: nside

type(lifted_trajectory),intent(in) :: lifted
type(lagged_healpix_density),intent(in) :: hdens
type(healpix_marginals),intent(out) :: marg

integer :: npix,nob
integer :: i,j,ipix1,ipix2,ipix3
integer :: a,b,c


npix=12*nside**2
nob=2*lifted%nobject


allocate(marg%xy(0:npix-1,0:npix-1,nob,nob))
allocate(marg%xz(0:npix-1,0:npix-1,nob,nob))
allocate(marg%yz(0:npix-1,0:npix-1,nob,nob))

allocate(marg%x(0:npix-1,nob,nob))
allocate(marg%y(0:npix-1,nob,nob))
allocate(marg%z(0:npix-1,nob,nob))


marg%xy=0
marg%xz=0
marg%yz=0

marg%x=0
marg%y=0
marg%z=0


do i=1,nob
do j=1,nob

   do ipix1=0,npix-1
   do ipix2=0,npix-1
   do ipix3=0,npix-1


      if(.not.hdens%occ(ipix1,ipix2,ipix3,i,j)) cycle


      a=hdens%raw_counts(ipix1,ipix2,ipix3,i,j)


      ! xy = sum über z
      marg%xy(ipix1,ipix2,i,j)= &
      marg%xy(ipix1,ipix2,i,j)+a


      ! xz = sum über y
      marg%xz(ipix1,ipix3,i,j)= &
      marg%xz(ipix1,ipix3,i,j)+a


      ! yz = sum über x
      marg%yz(ipix2,ipix3,i,j)= &
      marg%yz(ipix2,ipix3,i,j)+a


      ! Einzelmarginale

      marg%x(ipix1,i,j)= &
      marg%x(ipix1,i,j)+a

      marg%y(ipix2,i,j)= &
      marg%y(ipix2,i,j)+a

      marg%z(ipix3,i,j)= &
      marg%z(ipix3,i,j)+a


   end do
   end do
   end do

end do
end do

end subroutine marginalize_healpix_density



subroutine test_healpix_density(tstart,tend,big_interval,tau,plan,project,istate,irep,nside)
implicit none
integer,intent(in) :: big_interval,tau,tstart,tend,istate,irep,nside
type(lifted_trajectory)                   :: lifted
logical,allocatable                       :: window_mask(:)
type(project_type), intent(in)            :: project
type(raw_data_block)                      :: block
type(cv_plan), intent(in)                 :: plan
type(lagged_healpix_density),save:: hdens,hdens_shuffel
type(time_window) :: window
type(scalar_block_mean) :: scalar_stats
type(healpix_marginals) :: marg,marg_shuffel
real(8),allocatable     ::te(:,:),te_shuffel(:,:)
real(8),allocatable     ::te_old(:,:),te_shuffel_old(:,:)
integer :: ntstep,nt,i,j,k,ii,jj,ipix1,ipix2,ipix3,maxi,c,ki,kj,total
real(8) :: neff
integer :: nob
integer :: n,nocc
real(8) :: occupancy
integer :: unit_te
integer :: unit_shuffle
integer :: unit_diff
integer,allocatable :: order(:)

character(len=32)  :: replicate_name
character(len=512) :: filename_te
character(len=512) :: filename_shuffle
character(len=512) :: filename_diff

unit_te      = 20
unit_shuffle = 21
unit_diff    = 22

write(replicate_name,'("replicate_",I0)') irep


filename_te = trim(project%rootdir)//'/output/'// &
              trim(replicate_name)//'_te_matrix.txt'

filename_shuffle = trim(project%rootdir)//'/output/'// &
                   trim(replicate_name)//'_te_shuffle_matrix.txt'

filename_diff = trim(project%rootdir)//'/output/'// &
                trim(replicate_name)//'_te_difference.txt'


open(unit_te, file=trim(filename_te), &
     status="replace", action="write")

open(unit_shuffle, file=trim(filename_shuffle), &
     status="replace", action="write")

open(unit_diff, file=trim(filename_diff), &
     status="replace", action="write")
!logical,allocatable                       ::nocc(:,:)

ntstep=tend-tstart+1

window%istate =istate
window%irep   =irep
window%start_time_step =tstart
window%end_time_step   =tend
allocate(window_mask(ntstep))

window_mask=.false.

do i=1,ntstep,big_interval
if(i+tau>ntstep) exit
window_mask(i)    =.true.
window_mask(i+tau)=.true.
end do

nt=count(window_mask)
write(0,*) nt

call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)
  
call lift_cv_plan(plan, block, lifted)
call normalize_and_rotate_vec(lifted,8,scalar_stats)

if(irep==1) then
	
call healpix_density_block_spaced_pairs(nside,lifted,hdens,.false.) 
write(0,*) allocated(	hdens%raw_counts),irep
call marginalize_healpix_density(nside,lifted,hdens,marg) 
write(0,*) allocated(	hdens%raw_counts),irep
call healpix_density_block_shuffle(nside,lifted,hdens_shuffel,.false.)
write(0,*) allocated(	hdens_shuffel%raw_counts),irep
call marginalize_healpix_density(nside,lifted,hdens_shuffel,marg_shuffel) 
write(0,*) allocated(	hdens%raw_counts),irep
else
	write(0,*) allocated(	hdens%raw_counts),allocated(	hdens_shuffel%raw_counts)
!pause
	
call healpix_density_block_spaced_pairs(nside,lifted,hdens,.true.)    
call marginalize_healpix_density(nside,lifted,hdens,marg)  
call healpix_density_block_shuffle(nside,lifted,hdens_shuffel,.true.)
call marginalize_healpix_density(nside,lifted,hdens_shuffel,marg_shuffel) 


end if


!do i=1,size(hdens%raw_counts,4)
!do j=1,size(hdens%raw_counts,4)
!
!total=sum(hdens%raw_counts(:,:,:,i,j))
!
!write(*,*) i,j,total, &
!sum(marg%xy(:,:,i,j)), &
!sum(marg%yz(:,:,i,j)), &
!sum(marg%y(:,i,j))
!
!enddo
!enddo

call compute_transfer_entropy(nside,hdens,marg,te)

call compute_transfer_entropy(nside,hdens_shuffel,marg_shuffel,te_shuffel)



write(0,*) maxval(te),'maxval te',maxloc(te)

do i=1,size(te,1)
   write(0,*) i,maxval(te(i,:)),maxval(te_shuffel(i,:))
end do

!pause
write(0,*) "hdens shape =",shape(hdens%occ)
write(0,*) "block xyz =",shape(block%xyz)


!
! Statistik der einzelnen TE-Paare
!




nob=size(hdens%raw_counts,4)

write(0,*) "i j N Neff Nocc Nocc/N TE"
call make_object_order(plan, order)
allocate(        te_old(nob,nob))
allocate(te_shuffel_old(nob,nob))

te_old=te
te_shuffel_old=te_shuffel
do i = 1, nob

   do j = 1, nob

      ! alter Index
      ii = order(i)
      jj = order(j)

      ! Anzahl aller Übergänge für dieses Paar
      n = sum(hdens%raw_counts(:,:,:,ii,jj))

      if (n == 0) then
         te(i,j)         = 0.d0
         te_shuffel(i,j) = 0.d0
      else

         nocc = count(hdens%occ(:,:,:,ii,jj))

         neff = real(n,8)**2 / &
                sum(real(hdens%raw_counts(:,:,:,ii,jj),8)**2)

         occupancy = real(nocc,8) / real(n,8)

         ! Falls te bereits in alter Reihenfolge berechnet wurde:
         te(i,j)         = te_old(ii,jj)
         te_shuffel(i,j) = te_shuffel_old(ii,jj)

         write(0,'(2I5,4F12.4,1X,F12.6,F12.6)') &
                i, j, &
                real(n,8), &
                neff, &
                real(nocc,8), &
                occupancy, &
                te(i,j), te_shuffel(i,j)

      end if

   end do

   write(unit_te,*)      te(i,:)
   write(unit_shuffle,*) te_shuffel(i,:)
   write(unit_diff,*)    te(i,:) - te_shuffel(i,:)

end do


end subroutine test_healpix_density



subroutine healpix_density_block_spaced_pairs(nside,lifted,hdens,accumulate)

implicit none

integer, intent(in) :: nside
logical, intent(in) :: accumulate

type(lifted_trajectory), intent(in)    :: lifted
type(lagged_healpix_density), intent(inout) :: hdens

integer :: npix,nob
integer :: ipix1,ipix2,ipix3
integer :: t,i,j,ki,kj

real(8) :: x1(3),x2(3),x3(3)

nob  = 2*lifted%nobject
npix = 12*nside**2

!--------------------------------------------------
! Nur beim ersten Aufruf allokieren/initialisieren
!--------------------------------------------------
if (.not.accumulate) then

    if (allocated(hdens%occ))        deallocate(hdens%occ)
    if (allocated(hdens%raw_counts)) deallocate(hdens%raw_counts)

    allocate(hdens%occ(0:npix-1,0:npix-1,0:npix-1,nob,nob))
    allocate(hdens%raw_counts(0:npix-1,0:npix-1,0:npix-1,nob,nob))

    hdens%occ        = .false.
    hdens%raw_counts = 0

else

    ! Sicherheitsprüfung
    if (.not.allocated(hdens%raw_counts)) then
        stop "healpix_density_block_spaced_pairs: accumulate=.true. aber hdens nicht initialisiert."
    endif

endif

!--------------------------------------------------
! Histogramm aufbauen
!--------------------------------------------------

do ki=1,2

    do t=1,lifted%ntimes-1,2

        do i=1,lifted%nobject

            x1 = lifted%vec(:,i,t+1,ki)
            x2 = lifted%vec(:,i,t  ,ki)

            call vec2pix_ring(nside,x1,ipix1)
            call vec2pix_ring(nside,x2,ipix2)

            do j=1,lifted%nobject

                do kj=1,2

                    x3 = lifted%vec(:,j,t,kj)

                    call vec2pix_ring(nside,x3,ipix3)

                    hdens%raw_counts(ipix1,ipix2,ipix3, &
                                     2*(i-1)+ki, &
                                     2*(j-1)+kj) = &
                    hdens%raw_counts(ipix1,ipix2,ipix3, &
                                     2*(i-1)+ki, &
                                     2*(j-1)+kj) + 1

                    hdens%occ(ipix1,ipix2,ipix3, &
                              2*(i-1)+ki, &
                              2*(j-1)+kj) = .true.

                end do

            end do

        end do

    end do

end do

end subroutine healpix_density_block_spaced_pairs


subroutine test_lift_cv_plan(plan,project)
implicit none
type(lifted_trajectory)                :: lifted
logical,allocatable                    :: window_mask(:)
type(project_type), intent(in)         :: project
type(raw_data_block)                   :: block
type(cv_plan), intent(in)              :: plan
type(time_window) :: window
integer :: ntstep,i,t


window%istate =1
window%irep   =1
window%start_time_step =1
window%end_time_step   =10
ntstep=window%end_time_step-window%start_time_step+1

allocate(window_mask(ntstep)) 
window_mask=.false.
do i=1,ntstep,1
window_mask(i)=.true.
end do


write(0,*) window_mask
!pause

call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)
lifted%ntimes=ntstep

!allocate(lifted%vec(3,size(plan%object),count(window_mask),2))

call lift_cv_plan(plan, block, lifted)

do i=1,size(block%xyz,2)
do t=1,size(block%xyz,3)
write(0,*) block%xyz(:,i,t)
end do
end do

!pause
do i=1,lifted%nobject
do t=1,lifted%ntimes
write(0,*)lifted%vec(:,i,t,1)
end do
end do
write(0,*) lifted%nobject
end subroutine



subroutine lift_cv_plan(plan, block, lifted)

implicit none

type(cv_plan), intent(in)              :: plan
type(raw_data_block), intent(in)       :: block
type(lifted_trajectory), intent(out)   :: lifted

integer :: nobject
integer :: ntstep

integer :: t,k,obj,comp,atom

real(8) :: norm
real(8) :: v(3)
!----------------------------------------
! Dimensionen
!----------------------------------------
nobject=size(plan%object)
lifted%ntimes=size(block%xyz,3)
lifted%nobject=nobject
!write(0,*) shape(block%xyz)
!write(0,*) lifted%ntimes,'ntimes in lifting routine'
!pause
allocate(lifted%vec(3,nobject,lifted%ntimes,2))

lifted%vec=0.d0


!----------------------------------------
! Frames
!----------------------------------------
do t=1,lifted%ntimes

    !------------------------------------
    ! CV accumulation
    !------------------------------------
    do k=1,size(plan%term)

        atom = plan%term(k)%atom
        obj  = plan%term(k)%target
        comp = plan%term(k)%comp

        lifted%vec(:,obj,t,comp) = &
        lifted%vec(:,obj,t,comp) + &
        block%xyz(:,atom,t)*plan%term(k)%coeff
    enddo
    !!------------------------------------
    !! Normierung und Länge speichern
    !!------------------------------------
    !do obj=1,nobject
    !    do comp=1,2
    !        v = lifted%vec(:,obj,t,comp)
    !        norm=sqrt(dot_product(v,v))
    !        lifted%length(comp,obj,t)=norm
    !        if(norm > tiny(1.d0)) then
    !            lifted%vec(:,obj,t,comp)=v/norm
    !        else
    !            lifted%vec(:,obj,t,comp)=0.d0
    !        endif
    !    enddo
    !enddo
enddo

end subroutine


subroutine make_object_order(plan, order)

    implicit none

    type(cv_plan), intent(in) :: plan
    integer, allocatable, intent(out) :: order(:)

    integer :: i, j, n, tmp
    integer, allocatable :: obj_order(:)

    !----------------------------------------
    ! Anzahl der Objekte (Residuen + Domains)
    !----------------------------------------
    n = size(plan%object)

    allocate(obj_order(n))

    do i = 1, n
        obj_order(i) = i
    end do

    !----------------------------------------
    ! Objekte nach Typ und ID sortieren
    !----------------------------------------
    do i = 1, n-1
        do j = i+1, n

            if (plan%object(obj_order(j))%object_type < &
                plan%object(obj_order(i))%object_type) then

                tmp = obj_order(i)
                obj_order(i) = obj_order(j)
                obj_order(j) = tmp

            elseif (plan%object(obj_order(j))%object_type == &
                    plan%object(obj_order(i))%object_type) then

                if (plan%object(obj_order(j))%object_id < &
                    plan%object(obj_order(i))%object_id) then

                    tmp = obj_order(i)
                    obj_order(i) = obj_order(j)
                    obj_order(j) = tmp

                end if
            end if

        end do
    end do

    !----------------------------------------
    ! Objekt-Permutation -> CV-Permutation
    !
    ! Reihenfolge:
    ! Obj1_CV1
    ! Obj1_CV2
    ! Obj2_CV1
    ! Obj2_CV2
    ! ...
    !----------------------------------------
    allocate(order(2*n))

    do i = 1, n

        order(2*i-1) = 2*obj_order(i)-1
        order(2*i  ) = 2*obj_order(i)

    end do

    deallocate(obj_order)

end subroutine make_object_order



subroutine make_cv_plan(topologies,sc_lib,dom_lib,plan)
implicit none

type(topology_set_type), intent(in)      :: topologies
type(sidechain_library_type), intent(in) :: sc_lib
type(domain_library_type), intent(in)    :: dom_lib
type(cv_plan), intent(out)               :: plan

integer :: nterm,nobj
integer :: i,j,k,ncvtype,ires,idom
integer :: atom
integer :: term_idx,obj_idx
integer :: istate

!---------------------------------------------------
! 1) Erst Größe bestimmen
!---------------------------------------------------

nterm=0
nobj=0

! Domains
do idom=1,dom_lib%ndo

    nobj=nobj+1
    nterm=nterm+2*dom_lib%domain(idom)%natom_keep

enddo


! Residuen
do ncvtype=1,size(sc_lib%residue)

    if (.not.allocated(sc_lib%residue(ncvtype)%cv)) cycle

    do j=1,sc_lib%residue(ncvtype)%freq

        nobj=nobj+1

        nterm=nterm+2*sc_lib%residue(ncvtype)%natom_head

    enddo

enddo


allocate(plan%term(nterm))
allocate(plan%object(nobj))


!---------------------------------------------------
! 2) Terms erzeugen
!---------------------------------------------------

term_idx=0
obj_idx=0


!===================================================
! Domains
!===================================================

do idom=1,dom_lib%ndo

    obj_idx=obj_idx+1

    plan%object(obj_idx)%object_type=1
    plan%object(obj_idx)%object_id=idom

    plan%object(obj_idx)%first_term=term_idx+1
    ! CV1
    do i=1,dom_lib%domain(idom)%natom_keep

        term_idx=term_idx+1

        plan%term(term_idx)%atom = &
            dom_lib%domain(idom)%map(i)+ &
            dom_lib%domain(idom)%start_atom-1

        plan%term(term_idx)%target=obj_idx
        plan%term(term_idx)%comp=1
        plan%term(term_idx)%coeff= &
            dom_lib%domain(idom)%cv(1,i)
    enddo
    ! CV2
    do i=1,dom_lib%domain(idom)%natom_keep

        term_idx=term_idx+1

        plan%term(term_idx)%atom = &
            dom_lib%domain(idom)%map(i)+ &
            dom_lib%domain(idom)%start_atom-1

        plan%term(term_idx)%target=obj_idx
        plan%term(term_idx)%comp=2
        plan%term(term_idx)%coeff= &
            dom_lib%domain(idom)%cv(2,i)
    enddo
    plan%object(obj_idx)%last_term=term_idx

enddo
!===================================================
! Residuen
!===================================================
istate=1
do ncvtype=1,size(sc_lib%residue)
    if (.not.allocated(sc_lib%residue(ncvtype)%cv)) cycle
    do j=1,sc_lib%residue(ncvtype)%freq
        ires=sc_lib%residue(ncvtype)%res_ids(j)
        obj_idx=obj_idx+1
        plan%object(obj_idx)%object_type=2
        plan%object(obj_idx)%object_id=ires
        plan%object(obj_idx)%first_term=term_idx+1

        ! CV1
        do i=1,sc_lib%residue(ncvtype)%natom_head
            term_idx=term_idx+1

            atom = sc_lib%residue(ncvtype)%map(i) + &
                   topologies%state(istate)%topo% &
                   acid_start_list(ires)-1

            plan%term(term_idx)%atom=atom
            plan%term(term_idx)%target=obj_idx
            plan%term(term_idx)%comp=1
            plan%term(term_idx)%coeff= &
                 sc_lib%residue(ncvtype)%cv(1,i)
        enddo
        ! CV2
        do i=1,sc_lib%residue(ncvtype)%natom_head

            term_idx=term_idx+1


            atom = sc_lib%residue(ncvtype)%map(i) + &
                   topologies%state(istate)%topo% &
                   acid_start_list(ires)-1


            plan%term(term_idx)%atom=atom
            plan%term(term_idx)%target=obj_idx
            plan%term(term_idx)%comp=2
            plan%term(term_idx)%coeff= &
                 sc_lib%residue(ncvtype)%cv(2,i)
        enddo
        plan%object(obj_idx)%last_term=term_idx
    enddo
enddo


end subroutine

subroutine healpix_density_block_shuffle(nside,lifted,hdens,accumulate)

implicit none

integer, intent(in) :: nside
logical, intent(in) :: accumulate

type(lifted_trajectory), intent(in) :: lifted
type(lagged_healpix_density), intent(inout) :: hdens

integer :: ipix1,ipix2,ipix3
integer :: npix,nob
integer :: t,ts
integer :: i,j,ki,kj

real(8) :: x1(3),x2(3),x3(3)

type(shuffled_time_map) :: shuffle

nob  = 2*lifted%nobject
npix = 12*nside**2

!--------------------------------------------------
! Nur beim ersten Aufruf allokieren
!--------------------------------------------------
if (.not.accumulate) then

    if (allocated(hdens%occ))        deallocate(hdens%occ)
    if (allocated(hdens%raw_counts)) deallocate(hdens%raw_counts)

    allocate(hdens%occ(0:npix-1,0:npix-1,0:npix-1,nob,nob))
    allocate(hdens%raw_counts(0:npix-1,0:npix-1,0:npix-1,nob,nob))

    hdens%occ        = .false.
    hdens%raw_counts = 0

else

    if (.not.allocated(hdens%raw_counts)) then
        stop "healpix_density_block_shuffle: accumulate=.true. aber hdens nicht initialisiert."
    endif

endif

!--------------------------------------------------
! Neue Zufallspermutation für dieses Replikat
!--------------------------------------------------
call generate_shuffle(lifted%ntimes,shuffle)

!--------------------------------------------------
! Histogramm aufbauen
!--------------------------------------------------

do ki = 1,2

    do t = 1,lifted%ntimes-1,2

        ts = shuffle%index(t)

        do i = 1,lifted%nobject

            x1 = lifted%vec(:,i,t+1,ki)
            x2 = lifted%vec(:,i,t  ,ki)

            call vec2pix_ring(nside,x1,ipix1)
            call vec2pix_ring(nside,x2,ipix2)

            do j = 1,lifted%nobject

                do kj = 1,2

                    x3 = lifted%vec(:,j,ts,kj)

                    call vec2pix_ring(nside,x3,ipix3)

                    hdens%raw_counts( &
                        ipix1,ipix2,ipix3, &
                        2*(i-1)+ki, &
                        2*(j-1)+kj ) = &
                    hdens%raw_counts( &
                        ipix1,ipix2,ipix3, &
                        2*(i-1)+ki, &
                        2*(j-1)+kj ) + 1

                    hdens%occ( &
                        ipix1,ipix2,ipix3, &
                        2*(i-1)+ki, &
                        2*(j-1)+kj ) = .true.

                end do

            end do

        end do

    end do

end do

end subroutine healpix_density_block_shuffle



subroutine test_lifting(topologies,sc_lib,dom_lib,project,bidx1,bidx2)
implicit none

type(time_window) :: window
type(project_type), intent(in)           ::project
type(sidechain_library_type), intent(in) ::sc_lib
type(domain_library_type) , intent(in)   ::dom_lib
type(topology_set_type)   , intent(in)   ::topologies
type(raw_data_block) :: block
logical,allocatable  ::activ_res(:)
logical,allocatable  ::activ_dom(:)
logical,allocatable  ::window_mask(:) 
integer :: ntstep
integer, intent(in) :: bidx1,bidx2
real(8),allocatable :: bc11(:),bc12(:),bc21(:),bc22(:)
integer,allocatable :: m1(:),m2(:)
integer :: ds1,ds2, iatom1,iatom2,iatom,ires,s1,s2,i,j,k,t,idom,ncvtype 
real(8) :: b11(3),b12(3),b21(3),b22(3),b13(3),b23(3),b1(3),b2(3),b1_res(3),b2_res(3)
real(8) :: v1(3,2),v2(3,2),l11,l12,l21,l22,or1(3,3),or2(3,3),rel(3,3),rel_last(3,3),shift_rel(3,3)
real(8) :: lr1,lr2

allocate(bc11(dom_lib%domain(bidx1)%natom_keep))
allocate(bc12(dom_lib%domain(bidx1)%natom_keep))
allocate(m1(dom_lib%domain(bidx1)%natom_keep))
bc11=dom_lib%domain(bidx1)%cv(1,:)
bc12=dom_lib%domain(bidx1)%cv(2,:)
m1=dom_lib%domain(bidx1)%map


allocate(bc21(dom_lib%domain(bidx2)%natom_keep))
allocate(bc22(dom_lib%domain(bidx2)%natom_keep))
allocate(m2(dom_lib%domain(bidx2)%natom_keep))
bc21=dom_lib%domain(bidx2)%cv(1,:)
bc22=dom_lib%domain(bidx2)%cv(2,:)
m2=dom_lib%domain(bidx2)%map

window%istate =1
window%irep   =1
window%start_time_step =1
window%end_time_step   =100000
ntstep=window%end_time_step-window%start_time_step+1

allocate(window_mask(ntstep)) 
window_mask=.false.
do i=1,ntstep,1000
window_mask(i)=.true.
end do

call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)


ds1=dom_lib%domain(bidx1)%start_atom-1
s1=dom_lib%domain(bidx1)%natom_keep

ds2=dom_lib%domain(bidx2)%start_atom-1
s2=dom_lib%domain(bidx2)%natom_keep

rel_last=0._8
rel=0._8
do t=1,count(window_mask)
v1=0._8
v2=0._8
b1=0._8
b2=0._8
idom=bidx1
    do i=1,dom_lib%domain(idom)%natom_keep     
        iatom1=dom_lib%domain(idom)%map(i)+dom_lib%domain(idom)%start_atom-1
        b1=b1+block%xyz(:,iatom1,t)*dom_lib%domain(idom)%cv(1,i)
        b2=b2+block%xyz(:,iatom1,t)*dom_lib%domain(idom)%cv(2,i)        
    end do

          l11=dot_product(b1,b1)
          l12=dot_product(b2,b2)
          l11=sqrt(l11)
          l12=sqrt(l12)
          b11=b1/l11
          b12=b2/l12
        or1(:,1)=b11
        or1(:,2)=b12-dot_product(b11,b12)*b11
        or1(:,3)=cross_product(b11,b12)
        

b1=0._8
b2=0._8
idom=bidx2
  do i=1,dom_lib%domain(idom)%natom_keep     
        iatom1=dom_lib%domain(idom)%map(i)+dom_lib%domain(idom)%start_atom-1
        b1=b1+block%xyz(:,iatom1,t)*dom_lib%domain(idom)%cv(1,i)
        b2=b2+block%xyz(:,iatom1,t)*dom_lib%domain(idom)%cv(2,i)        
    end do
          l21=dot_product(b1,b1)
          l22=dot_product(b2,b2)
          l21=sqrt(l21)
          l22=sqrt(l22)
          b21=b1/l21
          b22=b2/l22
        or2(:,1)=b21
        or2(:,2)=b22
        or2(:,3)=cross_product(b21,b22)
        rel_last=rel
        rel=matmul(or2,transpose(or1))
        shift_rel=matmul(rel_last,transpose(rel))
        !write(0,*) rel(1,:)
        !write(0,*) rel(2,:)
        !write(0,*) rel(3,:)
        !write(0,*)
        	
        write(0,*) shift_rel(1,:)
        write(0,*) shift_rel(2,:)
        write(0,*) shift_rel(3,:)
        write(0,*)
        
    do ncvtype = 1, size(sc_lib%residue)

    if (.not. allocated(sc_lib%residue(ncvtype)%cv)) cycle

    do j = 1, sc_lib%residue(ncvtype)%freq

        ires = sc_lib%residue(ncvtype)%res_ids(j)

        b1_res = 0.d0
        b2_res = 0.d0

        do i = 1, sc_lib%residue(ncvtype)%natom_head

            iatom = sc_lib%residue(ncvtype)%map(i) + &
                    topologies%state(window%istate)%topo%acid_start_list(ires) - 1

            b1_res = b1_res + block%xyz(:,iatom,t) * &
                               sc_lib%residue(ncvtype)%cv(1,i)

            b2_res = b2_res + block%xyz(:,iatom,t) * &
                               sc_lib%residue(ncvtype)%cv(2,i)
                               
      
 
         end do
           lr1=dot_product(b1_res,b1_res)
           lr1=sqrt(lr1)
           b1_res= b1_res/lr1
           
           lr2=dot_product(b2_res,b2_res)
           lr2=sqrt(lr2)
           b2_res= b2_res/lr2
         
         !write(0,*) sc_lib%residue(ncvtype)%cv(1,:),'cv(1,:)'
         !write(0,*) sc_lib%residue(ncvtype)%cv(2,:),'cv(2,:)'
         
         write(0,*) dot_product(b1_res,b2_res)!,b1_res,b2_res
         b1_res=matmul(transpose(or1),b1_res)
         b2_res=matmul(transpose(or1),b2_res)
         write(0,*) dot_product(b1_res,b2_res)
         write(0,*) 
    end do

end do
end do


end subroutine




subroutine test_simple_vector_cv_statisics(topologies,sc_lib,dom_lib,project,stats_b,stats_s)
implicit none

type(time_window) :: window
type(project_type), intent(in)           ::project
type(sidechain_library_type), intent(in) ::sc_lib
type(domain_library_type) , intent(in)   ::dom_lib
type(topology_set_type)   , intent(in)   ::topologies

logical,allocatable  ::activ_res(:)
logical,allocatable  ::activ_dom(:)
logical,allocatable  ::window_mask(:) 
integer :: ntstep
real(8),allocatable :: stats_b(:),stats_s(:)

!call sidechain_cv_inventur(topologies,1,sc_lib)


window%istate =1
window%irep   =1
window%start_time_step =1
window%end_time_step   =1000
ntstep=window%end_time_step-window%start_time_step+1

allocate(activ_dom(dom_lib%ndo))
allocate(window_mask(ntstep))
activ_dom=.true.
window_mask=.false.
window_mask=.true.
!window_mask(ntstep-10:ntstep)=.true.
call simple_vector_cv_statisics(topologies,window,window_mask,activ_res,activ_dom,sc_lib,dom_lib,project,stats_b,stats_s)

end subroutine




subroutine simple_vector_cv_statisics(topologies,window,window_mask,activ_res&
,activ_dom,sc_lib,dom_lib,project,stats_b,stats_s)
implicit none
type(topology_set_type), intent(in)      ::topologies
type(project_type), intent(in)           ::project
type(domain_library_type) , intent(in)   ::dom_lib
type(sidechain_library_type), intent(in) ::sc_lib
type(time_window),intent(in) :: window
type(raw_data_block) :: block
type(window_map)     :: map

logical,allocatable, intent(in) :: activ_res(:),activ_dom(:)
logical,allocatable, intent(in)  :: window_mask(:) 
integer :: time_step,i,iatom,idom,ntstep
real(8) :: b1(3),b2(3),s1(3),s2(3),lb1,lb2,b12,basis(3,3),t1(3),t2(3),t3(3),t4(3),t5(3)
real(8) :: b1_res(3),b2_res(3),s1_res(3),s2_res(3),lb1_res,lb2_res,b12_res,basis_res(3,3)
real(8) ::q1(3),q2(3)
real(8),allocatable :: l1diffsq(:),l2diffsq(:),b12diffsq(:),l1s,l2s,b12s
real(8),allocatable :: l1diffsq_res(:),l2diffsq_res(:),b12diffsq_res(:),l1s_res,l2s_res,b12s_res
real(8),allocatable :: stats_b(:),stats_s(:)
integer :: ndomain,nres,ncv,ncvtype,j, ires


!call get_window_map(window, project, map)
call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)
    
    ndomain=dom_lib%ndo
    
allocate( l1diffsq(ndomain))
allocate( l2diffsq(ndomain))
allocate(b12diffsq(ndomain))

nres=topologies%state(window%istate)%topo%nres
allocate( l1diffsq_res(nres))
allocate( l2diffsq_res(nres))
allocate(b12diffsq_res(nres))

ntstep=count(window_mask)

 l1diffsq=0._8
 l2diffsq=0._8
b12diffsq=0._8

l1diffsq_res = 0.d0
l2diffsq_res = 0.d0
b12diffsq_res = 0.d0

do time_step=1,ntstep




    do idom=1,dom_lib%ndo
    
    b1=0._8
    b2=0._8
    s1=0._8
    s2=0._8
    
        do i=1,dom_lib%domain(idom)%natom_keep
        
        iatom=dom_lib%domain(idom)%map(i)+dom_lib%domain(idom)%start_atom-1
        b1=b1+block%xyz(:,iatom,time_step)*dom_lib%domain(idom)%cv(1,i)
        b2=b2+block%xyz(:,iatom,time_step)*dom_lib%domain(idom)%cv(2,i)
        
        end do
        !write(0,*) block%xyz(:,iatom,time_step),'#############'
        !write(0,*)shape(block%xyz)
    lb1=dot_product(b1,b1)
    lb2=dot_product(b2,b2)
    lb1=sqrt(lb1)
    lb2=sqrt(lb2)
    b1=b1/lb1
    b2=b2/lb2
    b12=dot_product(b1,b2)
    !write(0,*) lb1,lb2,b12
    
    l1s =1._8-lb1
    l2s =1._8-lb1
    b12s=b12
             
     l1diffsq(idom)= l1diffsq(idom)+l1s *l1s 
     l2diffsq(idom)= l2diffsq(idom)+l2s *l2s 
    b12diffsq(idom)=b12diffsq(idom)+b12s*b12s
    
    basis(1,:)=b1
    basis(2,:)=b2
    basis(3,:)=cross_product(b1,b2)
    
    
    end do
    


do ncvtype = 1, size(sc_lib%residue)

    if (.not. allocated(sc_lib%residue(ncvtype)%cv)) cycle

    do j = 1, sc_lib%residue(ncvtype)%freq

        ires = sc_lib%residue(ncvtype)%res_ids(j)

        b1_res = 0.d0
        b2_res = 0.d0

        do i = 1, sc_lib%residue(ncvtype)%natom_head

            iatom = sc_lib%residue(ncvtype)%map(i) + &
                    topologies%state(window%istate)%topo%acid_start_list(ires) - 1

            b1_res = b1_res + block%xyz(:,iatom,time_step) * &
                               sc_lib%residue(ncvtype)%cv(1,i)

            b2_res = b2_res + block%xyz(:,iatom,time_step) * &
                               sc_lib%residue(ncvtype)%cv(2,i)

         !write(0,*) topologies%state(window%istate)%topo%atom_id(iatom),sc_lib%residue(ncvtype)%map(i)
        
         
       
         
         end do
         
         

!write(0,*) dot_product(b1_res,b2_res)       
!write(0,*) b1_res, b2_res,sc_lib%residue(ncvtype)%resname
        lb1_res = sqrt(dot_product(b1_res,b1_res))
        lb2_res = sqrt(dot_product(b2_res,b2_res))
!write(0,*)sc_lib%residue(ncvtype)%cv(1,:)
!write(0,*)sc_lib%residue(ncvtype)%cv(2,:)
!pause
        !if (lb1_res <= tiny(1.d0) .or. lb2_res <= tiny(1.d0)) cycle

        b1_res = b1_res / lb1_res
        b2_res = b2_res / lb2_res

        b12_res = dot_product(b1_res,b2_res)
        write(0,*) lb1_res,lb2_res,b12_res,'llc'&
        ,topologies%state(window%istate)%topo%resnames(ires)&
        ,ires,sc_lib%residue(ncvtype)%resname
        
        l1s_res  = 1.d0 - lb1_res
        l2s_res  = 1.d0 - lb2_res
        b12s_res = b12_res

        l1diffsq_res(ires)  = l1diffsq_res(ires)  + l1s_res*l1s_res
        l2diffsq_res(ires)  = l2diffsq_res(ires)  + l2s_res*l2s_res
        b12diffsq_res(ires) = b12diffsq_res(ires) + b12s_res*b12s_res

        basis_res(1,:) = b1_res
        basis_res(2,:) = b2_res
        basis_res(3,:) = cross_product(b1_res,b2_res)

    end do

end do
 
!pause
end do



 l1diffsq= l1diffsq/dble(ntstep)
 l2diffsq= l2diffsq/dble(ntstep)
b12diffsq=b12diffsq/dble(ntstep)





write(0,*)
write(0,*) 'local domain basis instability scores'
write(0,*)
do idom=1,ndomain
!write(0,*) dom_lib%domain(idom)%domain_name,l1diffsq(idom),l2diffsq(idom),b12diffsq(idom)
!write(0,*) dom_lib%domain(idom)%natom_keep/3,-log(l1diffsq(idom)),-log(l2diffsq(idom)),-log(b12diffsq(idom))
!write(0,*) idom,-log(l1diffsq(idom))-log(l2diffsq(idom))-log(b12diffsq(idom))
end do


 l1diffsq_res= l1diffsq_res/dble(ntstep)
 l2diffsq_res= l2diffsq_res/dble(ntstep)
b12diffsq_res=b12diffsq_res/dble(ntstep)


i=0
write(0,*)
write(0,*) 'local sidechain basis instability scores'
write(0,*)
do idom=2,topologies%state(window%istate)%topo%nres-1
if(topologies%state(window%istate)%topo%resnames(idom)=='GLY'.or. topologies%state(window%istate)%topo%resnames(idom)=='ALA') cycle
i=i+1
write(0,*) topologies%state(window%istate)%topo%resnames(idom)
write(0,*) l1diffsq_res(idom),l2diffsq_res(idom),b12diffsq_res(idom)
!write(0,*) -log(l1diffsq_res(idom)),-log(l2diffsq_res(idom)),-log(b12diffsq_res(idom))
!write(0,*) idom,-log(l1diffsq_res(idom))-log(l2diffsq_res(idom))-log(b12diffsq_res(idom))
!write(0,*)
end do

!write(0,*) ndomain,'???',dom_lib%ndo
if(allocated(stats_b)) deallocate(stats_b)
allocate(stats_b(dom_lib%ndo))


if(allocated(stats_s)) deallocate(stats_s)
allocate(stats_s(i))
i=0
do ncvtype = 1, size(sc_lib%residue)

    if (.not. allocated(sc_lib%residue(ncvtype)%cv)) cycle

    do j = 1, sc_lib%residue(ncvtype)%freq
        i=i+1
        ires = sc_lib%residue(ncvtype)%res_ids(j)
        
        stats_s(i)=b12diffsq_res(ires)

  end do
end do
        
        
        
stats_b=b12diffsq

!pause
end subroutine


subroutine sidechain_cv_inventur(topologies,istate,sc_lib)

type(sidechain_library_type), intent(in) ::sc_lib
type(topology_set_type), intent(in) :: topologies
integer, intent(in) :: istate
integer:: nres,ires,ncvtype,icv,freq,i,j,k,iatom


do ncvtype=1,size(sc_lib%residue)
   if( allocated(sc_lib%residue(ncvtype)%cv) .eqv. .false.) cycle
     write(0,*) allocated(sc_lib%residue(ncvtype)%cv),ncvtype
     write(0,*) sc_lib%residue(ncvtype)%cv(1,:)
     write(0,*) sc_lib%residue(ncvtype)%cv(2,:)
    do j=1,sc_lib%residue(ncvtype)%freq
    ires=sc_lib%residue(ncvtype)%res_ids(j)
    
    write(0,*) topologies%state(istate)%topo%resnames(ires),ires,sc_lib%residue(ncvtype)%freq
    
        do i=1,sc_lib%residue(ncvtype)%natom_head
        iatom=sc_lib%residue(ncvtype)%map(i)&
             +topologies%state(istate)%topo%acid_start_list(ires)-1
        write(0,*) iatom,topologies%state(istate)%topo%atom_id(iatom)
       
        end do
        
    end do
end do

!pause
end subroutine




function cross_product(x1,x2) result(c)

implicit none

real(8), intent(in) :: x1(3), x2(3)
real(8) :: c(3)

c(1) = x1(2)*x2(3) - x1(3)*x2(2)
c(2) = x1(3)*x2(1) - x1(1)*x2(3)
c(3) = x1(1)*x2(2) - x1(2)*x2(1)

end function cross_product

subroutine normalize_and_rotate_vec(lifted,bidx,scalar_stats)

implicit none

type(scalar_block_mean)              :: scalar_stats
type(lifted_trajectory), intent(inout) :: lifted

integer :: bidx,t,i,j

real(8) :: b1(3),b2(3),b3(3)
real(8) :: l1,l2
real(8) :: rot(3,3)
real(8) :: v(3),l,c

allocate(scalar_stats%mean_length(2,lifted%nobject))
allocate(scalar_stats%mean_squared_length(2,lifted%nobject))
allocate(scalar_stats%mean_length_cross(lifted%nobject))

do t=1,lifted%ntimes

    b1=lifted%vec(:,bidx,t,1)
    b2=lifted%vec(:,bidx,t,2)

    l1=norm2(b1)
    l2=norm2(b2)

    b1=b1/l1
    b2=b2/l2
    b2=b2-dot_product(b1,b2)*b1
    b3=cross_product(b1,b2)

    scalar_stats%mean_squared_length(1,bidx)= &
        scalar_stats%mean_squared_length(1,bidx)+l1*l1

    scalar_stats%mean_length(1,bidx)= &
        scalar_stats%mean_length(1,bidx)+l1

    scalar_stats%mean_squared_length(2,bidx)= &
        scalar_stats%mean_squared_length(2,bidx)+l2*l2

    scalar_stats%mean_length(2,bidx)= &
        scalar_stats%mean_length(2,bidx)+l2

    rot(:,1)=b1
    rot(:,2)=b2
    rot(:,3)=b3

    scalar_stats%mean_length_cross(bidx)=norm2(b3)

    do i=1,2
        do j=1,lifted%nobject

            v=lifted%vec(:,j,t,i)
            l=norm2(v)

            scalar_stats%mean_length(2,bidx)=l
            scalar_stats%mean_length(2,bidx)=l*l

            lifted%vec(:,j,t,i)=matmul(rot,v)

        end do
    end do

end do

end subroutine

!N   HIS X  86 N   HIS X  71
!H   HIS X  86 H   HIS X  71
!CA  HIS X  86 CA  HIS X  71
!HA  HIS X  86 HA  HIS X  71
!CB  HIS X  86 CB  HIS X  71
!HB1 HIS X  86 HB1 HIS X  71
!HB2 HIS X  86 HB2 HIS X  71
!CG  HIS X  86 CG  HIS X  71
!ND1 HIS X  86 ND1 HIS X  71
!CD2 HIS X  86 HD1 HIS X  71
!HD2 HIS X  86 CD2 HIS X  71
!CE1 HIS X  86 HD2 HIS X  71
!HE1 HIS X  86 CE1 HIS X  71
!NE2 HIS X  86 HE1 HIS X  71
!HE2 HIS X  86 NE2 HIS X  71
!C   HIS X  86 C   HIS X  71
!O   HIS X  86 O   HIS X  71
end module tra_lift_module