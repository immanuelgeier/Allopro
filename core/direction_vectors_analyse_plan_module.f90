module direction_vectors_analyse_plan_module


use topology_module
use sidechains
use domains_module
use project_module

implicit none


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


end module