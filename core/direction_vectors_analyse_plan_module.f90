module direction_vectors_analyse_plan_module

use topology_module
use sidechains
use domains_module
use user_defined_cv_module
use project_module
use time_window_file_mapping
use time_window_loader
!use sidechainhead_refstruc

implicit none


type cv_term
    integer :: atom      ! absoluter Atomindex
    integer :: target    ! welche CV?
    integer :: comp      ! 1=b1, 2=b2
    real(8) :: coeff
    logical :: has_c2_sym=.false.
end type

type cv_object
    integer :: first_term
    integer :: last_term
    integer :: object_type   ! domain / residue
    integer :: object_id
    integer :: resid=0
    integer :: domid=0
    character*4 ::point_group
    character(len=128) :: name(2)
    integer :: refid
end type

type cv_plan
    type(cv_term), allocatable   :: term(:)
    type(cv_object), allocatable :: object(:)
end type


contains



subroutine destroy_cv_plan(plan)
    type(cv_plan), intent(inout) :: plan

    if (allocated(plan%term))   deallocate(plan%term)
    if (allocated(plan%object)) deallocate(plan%object)

end subroutine destroy_cv_plan

subroutine get_atoms_sidechain(istate,sc_lib,topologies,icv,num_freq,atoms)
integer , intent(in)  :: icv,num_freq,istate
integer  :: ires,nh,i

integer, allocatable,intent(out):: atoms(:)
type(topology_set_type), intent(in)      :: topologies
type(sidechain_library_type), intent(inout) :: sc_lib

if(num_freq>sc_lib%residue(icv)%freq) then
	write(0,*) 'error in residue counting'
	stop
end if
nh=sc_lib%residue(icv)%natom_head
allocate(atoms(nh))
ires=sc_lib%residue(icv)%res_ids(num_freq)

do i=1,nh
          atoms(i) = &
                    sc_lib%residue(icv)%map(i) + &
                    topologies%state(istate)%topo% &
                    acid_start_list(ires) - 1
                    
end do                    
                    
end subroutine get_atoms_sidechain


            
subroutine get_atoms_domaine(istate,dom_lib,topologies,icv,atoms)
integer , intent(in)  :: icv,istate
integer  :: ires,nh,i

integer, allocatable,intent(out):: atoms(:)
type(topology_set_type), intent(in)      :: topologies
type(domain_library_type), intent(inout) :: dom_lib


nh=dom_lib%domain(icv)%natom_keep
allocate(atoms(nh))

do i=1,nh
          atoms(i)  =  dom_lib%domain(icv)%map(i) + &
            dom_lib%domain(icv)%start_atom - 1
                    
end do                    
                    
end subroutine get_atoms_domaine

subroutine get_atoms_user_cv(istate,user_lib,topologies,icv,atoms)
integer , intent(in)  :: icv,istate
integer  :: ires,nh,i

integer, allocatable,intent(out):: atoms(:)
type(topology_set_type), intent(in)      :: topologies
type(user_cv_library_type), intent(inout) :: user_lib


nh=user_lib%cvs(icv)%natom
allocate(atoms(nh))

do i=1,nh
          atoms(i)  =  user_lib%cvs(icv)%map(i) + &
            user_lib%cvs(icv)%start_atom - 1
                    
end do                    
                    
end subroutine get_atoms_user_cv
            
            
            
subroutine add_cv_pair_sidechain(plan, object_type, resid, refid, &
                       name, point_group, atoms, cv, has_c2_sym)

implicit none

type(cv_plan), intent(inout) :: plan
real(8), intent(in) :: cv(:,:)

integer, intent(in) :: object_type
integer, intent(in) :: resid
integer, intent(in) :: refid

character(len=*), intent(in) :: name(2)
character(len=*), intent(in) :: point_group

integer, intent(in) :: atoms(:)

logical, intent(in) :: has_c2_sym

integer :: old_nobj
integer :: old_nterm
integer :: n
integer :: i
integer :: object_id

type(cv_object), allocatable :: object_tmp(:)
type(cv_term),   allocatable :: term_tmp(:)

n = size(atoms)

if (size(cv,2) /= n) then
    write(0,*) 'ERROR add_cv_pair_sidechain: inconsistent sizes'
    stop
endif

! ---------------------------------------------------------
! Alte Größen
! ---------------------------------------------------------

if (allocated(plan%object)) then
    old_nobj = size(plan%object)
else
    old_nobj = 0
endif

if (allocated(plan%term)) then
    old_nterm = size(plan%term)
else
    old_nterm = 0
endif

object_id = old_nobj + 1

! ---------------------------------------------------------
! Neue Arrays erzeugen
! ---------------------------------------------------------

allocate(object_tmp(old_nobj+1))
allocate(term_tmp(old_nterm+2*n))

! ---------------------------------------------------------
! Alte Daten kopieren
! ---------------------------------------------------------

if (old_nobj > 0) then
    object_tmp(1:old_nobj) = plan%object
endif

if (old_nterm > 0) then
    term_tmp(1:old_nterm) = plan%term
endif

! ---------------------------------------------------------
! Neues Object
! ---------------------------------------------------------

object_tmp(old_nobj+1)%object_type = object_type
object_tmp(old_nobj+1)%object_id   = object_id
object_tmp(old_nobj+1)%resid       = resid
object_tmp(old_nobj+1)%refid       = refid

object_tmp(old_nobj+1)%first_term = old_nterm + 1
object_tmp(old_nobj+1)%last_term  = old_nterm + 2*n

object_tmp(old_nobj+1)%point_group = point_group
object_tmp(old_nobj+1)%name        = name

! ---------------------------------------------------------
! CV1
! ---------------------------------------------------------

do i=1,n

    term_tmp(old_nterm+i)%atom       = atoms(i)
    term_tmp(old_nterm+i)%target     = old_nobj+1
    term_tmp(old_nterm+i)%comp       = 1
    term_tmp(old_nterm+i)%coeff      = cv(1,i)
    term_tmp(old_nterm+i)%has_c2_sym = .false.

enddo

! ---------------------------------------------------------
! CV2
! ---------------------------------------------------------

do i=1,n

    term_tmp(old_nterm+n+i)%atom       = atoms(i)
    term_tmp(old_nterm+n+i)%target     = old_nobj+1
    term_tmp(old_nterm+n+i)%comp       = 2
    term_tmp(old_nterm+n+i)%coeff      = cv(2,i)
    term_tmp(old_nterm+n+i)%has_c2_sym = has_c2_sym

enddo

! ---------------------------------------------------------
! Arrays ersetzen
! ---------------------------------------------------------

call move_alloc(object_tmp,plan%object)
call move_alloc(term_tmp,plan%term)

end subroutine add_cv_pair_sidechain


subroutine add_domain_cv_pair(plan, domid, name, &
                              point_group, atoms, cv)

implicit none

type(cv_plan), intent(inout) :: plan

integer, intent(in) :: domid

character(len=*), intent(in) :: name(2)
character(len=*), intent(in) :: point_group

integer, intent(in) :: atoms(:)
real(8), intent(in) :: cv(:,:)

integer :: old_nobj
integer :: old_nterm
integer :: n
integer :: i
integer :: object_id

type(cv_object), allocatable :: object_tmp(:)
type(cv_term),   allocatable :: term_tmp(:)

n = size(atoms)

if (size(cv,2) /= n) then
    write(0,*) 'ERROR add_domain_cv_pair: inconsistent sizes'
    stop
endif

! ---------------------------------------------------------
! Alte Größen
! ---------------------------------------------------------

if (allocated(plan%object)) then
    old_nobj = size(plan%object)
else
    old_nobj = 0
endif

if (allocated(plan%term)) then
    old_nterm = size(plan%term)
else
    old_nterm = 0
endif

object_id = old_nobj + 1

! ---------------------------------------------------------
! Neue Arrays
! ---------------------------------------------------------

allocate(object_tmp(old_nobj+1))
allocate(term_tmp(old_nterm+2*n))

! ---------------------------------------------------------
! Alte Daten kopieren
! ---------------------------------------------------------

if (old_nobj > 0) then
    object_tmp(1:old_nobj) = plan%object
endif

if (old_nterm > 0) then
    term_tmp(1:old_nterm) = plan%term
endif

! ---------------------------------------------------------
! Neues Object
! ---------------------------------------------------------

object_tmp(old_nobj+1)%object_type = 1
object_tmp(old_nobj+1)%object_id   = object_id
object_tmp(old_nobj+1)%domid       = domid
object_tmp(old_nobj+1)%resid       = 0
object_tmp(old_nobj+1)%refid       = 0

object_tmp(old_nobj+1)%first_term = old_nterm + 1
object_tmp(old_nobj+1)%last_term  = old_nterm + 2*n

object_tmp(old_nobj+1)%point_group = point_group
object_tmp(old_nobj+1)%name        = name

! ---------------------------------------------------------
! CV1
! ---------------------------------------------------------

do i=1,n

    term_tmp(old_nterm+i)%atom       = atoms(i)
    term_tmp(old_nterm+i)%target     = old_nobj+1
    term_tmp(old_nterm+i)%comp       = 1
    term_tmp(old_nterm+i)%coeff      = cv(1,i)
    term_tmp(old_nterm+i)%has_c2_sym = .false.

enddo

! ---------------------------------------------------------
! CV2
! ---------------------------------------------------------

do i=1,n

    term_tmp(old_nterm+n+i)%atom       = atoms(i)
    term_tmp(old_nterm+n+i)%target     = old_nobj+1
    term_tmp(old_nterm+n+i)%comp       = 2
    term_tmp(old_nterm+n+i)%coeff      = cv(2,i)
    term_tmp(old_nterm+n+i)%has_c2_sym = .false.

enddo

! ---------------------------------------------------------
! Arrays ersetzen
! ---------------------------------------------------------

call move_alloc(object_tmp,plan%object)
call move_alloc(term_tmp,plan%term)

end subroutine add_domain_cv_pair

subroutine make_cv_names(names, basename, id)

    implicit none
            
    character(len=*), intent(out) :: names(2)
    character(len=*), intent(in)  :: basename
    integer, intent(in)           :: id

    write(names(1),'(A,"_",I0,"_1")') trim(basename), id
    write(names(2),'(A,"_",I0,"_2")') trim(basename), id

end subroutine make_cv_names
subroutine make_domain_names(names, basename)

    implicit none

    character(len=*), intent(out) :: names(2)
    character(len=*), intent(in)  :: basename

    names(1) = trim(basename)//'_1'
    names(2) = trim(basename)//'_2'

end subroutine make_domain_names

subroutine add_user_cv_pair(plan, istate, name, atoms, cv)

implicit none

type(cv_plan), intent(inout) :: plan

integer, intent(in) :: istate

character(len=*), intent(in) :: name(2)

integer, intent(in) :: atoms(:)
real(8), intent(in) :: cv(:,:)

integer :: old_nobj
integer :: old_nterm
integer :: n
integer :: i
integer :: object_id

type(cv_object), allocatable :: object_tmp(:)
type(cv_term),   allocatable :: term_tmp(:)

n = size(atoms)

if (size(cv,2) /= n) then
    write(0,*) 'ERROR add_user_cv_pair: inconsistent sizes'
    stop
endif

! ---------------------------------------------------------
! Alte Größen
! ---------------------------------------------------------

if (allocated(plan%object)) then
    old_nobj = size(plan%object)
else
    old_nobj = 0
endif

if (allocated(plan%term)) then
    old_nterm = size(plan%term)
else
    old_nterm = 0
endif

object_id = old_nobj + 1

! ---------------------------------------------------------
! Neue Arrays
! ---------------------------------------------------------

allocate(object_tmp(old_nobj+1))
allocate(term_tmp(old_nterm+2*n))

! ---------------------------------------------------------
! Alte Daten kopieren
! ---------------------------------------------------------

if (old_nobj > 0) then
    object_tmp(1:old_nobj) = plan%object
endif

if (old_nterm > 0) then
    term_tmp(1:old_nterm) = plan%term
endif

! ---------------------------------------------------------
! Neues Object (object_type = 3 für User-CVs, analog zu Domains/Sidechains)
! ---------------------------------------------------------

object_tmp(old_nobj+1)%object_type = 3
object_tmp(old_nobj+1)%object_id   = object_id
!object_tmp(old_nobj+1)%istate      = istate
object_tmp(old_nobj+1)%resid       = 0
object_tmp(old_nobj+1)%refid       = 0

object_tmp(old_nobj+1)%first_term = old_nterm + 1
object_tmp(old_nobj+1)%last_term  = old_nterm + 2*n

object_tmp(old_nobj+1)%point_group = 'none'
object_tmp(old_nobj+1)%name        = name

! ---------------------------------------------------------
! CV1
! ---------------------------------------------------------

do i=1,n

    term_tmp(old_nterm+i)%atom       = atoms(i)
    term_tmp(old_nterm+i)%target     = old_nobj+1
    term_tmp(old_nterm+i)%comp       = 1
    term_tmp(old_nterm+i)%coeff      = cv(1,i)
    term_tmp(old_nterm+i)%has_c2_sym = .false.

enddo

! ---------------------------------------------------------
! CV2
! ---------------------------------------------------------

do i=1,n

    term_tmp(old_nterm+n+i)%atom       = atoms(i)
    term_tmp(old_nterm+n+i)%target     = old_nobj+1
    term_tmp(old_nterm+n+i)%comp       = 2
    term_tmp(old_nterm+n+i)%coeff      = cv(2,i)
    term_tmp(old_nterm+n+i)%has_c2_sym = .false.

enddo

! ---------------------------------------------------------
! Arrays ersetzen
! ---------------------------------------------------------

call move_alloc(object_tmp,plan%object)
call move_alloc(term_tmp,plan%term)

end subroutine add_user_cv_pair


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
    call make_domain_names( &
    plan%object(obj_idx)%name, &
    dom_lib%domain(idom)%domain_name)

    plan%object(obj_idx)%first_term=term_idx+1
    plan%object(obj_idx)%domid=idom
    
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
    plan%object(obj_idx)%point_group='none'
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
        plan%object(obj_idx)%resid=ires
        plan%object(obj_idx)%first_term=term_idx+1
        
        call make_cv_names( &
        plan%object(obj_idx)%name, &
        sc_lib%residue(ncvtype)%resname, &
        ires)
        
        if(sc_lib%residue(ncvtype)%has_sym) then
            plan%object(obj_idx)%point_group='C2  '
        else
        	plan%object(obj_idx)%point_group='none'
         end if

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
        plan%term(obj_idx)%has_c2_sym=sc_lib%residue(ncvtype)%has_sym
    enddo
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

subroutine get_cv_name_list(plan,names)

    implicit none

    type(cv_plan), intent(in) :: plan
    character(len=128), allocatable, intent(out) :: names(:)

    integer :: i
    integer :: nobject

    nobject = size(plan%object)

    allocate(names(2*nobject))

    do i=1,nobject

        names(2*i-1) = trim(plan%object(i)%name(1))
        names(2*i)   = trim(plan%object(i)%name(2))

    enddo

end subroutine get_cv_name_list

subroutine  get_cv_residue_order(plan,order)

    implicit none

    type(cv_plan), intent(in) :: plan
    character(len=128) name
    integer,allocatable,intent(out):: order(:)
    integer :: i,k,nd,ns,ml
    integer :: nobject
    integer,allocatable :: resnums(:),perm(:)
    
    nobject = size(plan%object)

    allocate(order(2*nobject))
    nd=0
    ns=0
    do i=1,nobject

        name= trim(plan%object(i)%name(1))
       
        write(0,*)plan%object(i)%object_type&
        , plan%object(i)%name(1),plan%object(i)%object_type,i,plan%object(i)%domid,plan%object(i)%resid
        if(plan%object(i)%object_type==1) then
        	k=i
        order(2*i-1)=2*k-1
        order(2*i)=2*k
        nd=nd+1
        else
        ns=ns+1	
        
        end if
    enddo
    write(0,*) order(1:2*nd)
    pause
    
    allocate(resnums(ns),perm(ns))
      do i=nd+1,nobject

        name= trim(plan%object(i)%name(1))
       
        write(0,*)plan%object(i)%object_type&
        , plan%object(i)%name(1),plan%object(i)%object_type,i,plan%object(i)%domid,plan%object(i)%resid
         resnums(i-nd)=plan%object(i)%resid
    enddo
    
    do i=1,ns
    ml=minloc(resnums,1)
    perm(i)=ml
    resnums(ml)=huge(ml)
    end do
    
    do i=nd+1,nobject
    k=perm(i-nd)+nd
     order(2*i-1)=2*k-1
     order(2*i)=2*k
    end do
    
    do i=1,nobject
    write(0,*) order(2*i-1), order(2*i)
    end do
    
    pause

end subroutine get_cv_residue_order

end module