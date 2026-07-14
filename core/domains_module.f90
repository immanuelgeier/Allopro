module domains_module

use project_module
use topology_module
implicit none

private


public :: init_domain_maps
public :: domain_library_type

!=========================================================
type :: domain_cv_type

    character(3) :: domain_name = ''

    integer :: natom_all  = 0
    integer :: natom_keep = 0

    character(4) :: dir1_atom_1 = 'none'
    character(4) :: dir1_atom_2 = 'none'

    character(4) :: dir2_atom_1 = 'none'
    character(4) :: dir2_atom_2 = 'none'

    real(8) :: d1_ref(3) = 0.d0
    real(8) :: d2_ref(3) = 0.d0

    real(8), allocatable :: ref_xyz(:,:)
    real(8), allocatable :: cv(:,:)

    integer, allocatable :: map(:) !4*nres

    character(4), allocatable :: atom_names(:)
    integer :: nres_dom
    integer,allocatable :: res_ids(:)
    integer :: start_atom
    integer :: end_atom
end type

type :: domain_library_type
    integer :: ndo
    type(domain_cv_type), allocatable :: domain(:)

end type


contains

!=========================================================
! BUILD RESIDUE + ATOM MAPS
!=========================================================
subroutine init_domain_maps(state,domain_lib,topologies)

    implicit none

    type(state_topology_type), intent(inout)  :: state
    type(domain_library_type), intent(inout)  :: domain_lib
    type(topology_set_type)  , intent(in)     :: topologies


    integer :: i, j, res,iatom,ca
    integer :: nres, natom,nres_do
    integer :: first_res, last_res
    integer :: astart, aend

    nres  = state%topo%nres
    natom = state%topo%natom

    !====================================================
    ! residue -> domain map
    !====================================================
    if (allocated(state%domain_key)) deallocate(state%domain_key)
    allocate(state%domain_key(nres))
    state%domain_key = 0
    if(allocated(domain_lib%domain)) deallocate(domain_lib%domain)
    allocate(domain_lib%domain(state%ndomain))
    
    domain_lib%ndo=state%ndomain
    !====================================================
    ! LOOP OVER DOMAINS
    !====================================================
    do i = 1, state%ndomain
         write(0,*) 'Domain:',state%domain(i)%name
        !------------------------------------------------
        ! DEFINE RESIDUE RANGE (THIS IS THE KEY FIX)
        !------------------------------------------------
        first_res = state%domain(i)%first_res
        last_res  = state%domain(i)%last_res

        state%domain(i)%nres = last_res - first_res + 1

        !------------------------------------------------
        ! allocate and fill resid
        !------------------------------------------------
        if (allocated(state%domain(i)%resid)) then
            deallocate(state%domain(i)%resid)
        end if

        allocate(state%domain(i)%resid(state%domain(i)%nres))

        do j = 1, state%domain(i)%nres
            state%domain(i)%resid(j) = first_res + j - 1
        end do
        
         if (allocated(domain_lib%domain(i)%res_ids)) then
            deallocate(domain_lib%domain(i)%res_ids)
        end if

        allocate(domain_lib%domain(i)%res_ids(state%domain(i)%nres))

        do j = 1, state%domain(i)%nres
            domain_lib%domain(i)%res_ids(j) = first_res + j - 1
        end do

        !------------------------------------------------
        ! fill domain_key
        !------------------------------------------------
        do j = 1, state%domain(i)%nres

            res = state%domain(i)%resid(j)

            if (res >= 1 .and. res <= nres) then
                state%domain_key(res) = i
            end if

        end do

    end do

    !====================================================
    ! ATOM BOUNDARIES FROM TOPOLOGY
    !====================================================
    do i = 1, state%ndomain
        
        domain_lib%domain(i)%domain_name=state%domain(i)%name

        first_res = state%domain(i)%resid(1)
        last_res  = state%domain(i)%resid(state%domain(i)%nres)

        state%domain(i)%first_atom = &
            state%topo%acid_start_list(first_res)
        
        domain_lib%domain(i)%start_atom=state%domain(i)%first_atom
            
        state%domain(i)%last_atom  = &
            state%topo%acid_start_list(last_res) + &
            state%topo%acid_size(last_res) - 1
            
        domain_lib%domain(i)%end_atom=state%domain(i)%last_atom
            
            domain_lib%domain(i)%natom_all&
            =state%domain(i)%last_atom&
            -state%domain(i)%first_atom+1
            domain_lib%domain(i)%natom_keep&
            =state%domain(i)%nres*3
            domain_lib%domain(i)%nres_dom=state%domain(i)%nres
            allocate(domain_lib%domain(i)%map(domain_lib%domain(i)%natom_keep))
            ca=0
            do iatom=state%domain(i)%first_atom,state%domain(i)%last_atom
                   if(trim(state%topo%atom_id(iatom))=='N') then
                   ca=ca+1
                   domain_lib%domain(i)%map(ca)=iatom-state%domain(i)%first_atom+1
                   end if
                   if(trim(state%topo%atom_id(iatom))=='CA') then
                   ca=ca+1
                   domain_lib%domain(i)%map(ca)=iatom-state%domain(i)%first_atom+1
                   end if
                   if(trim(state%topo%atom_id(iatom))=='C') then
                   ca=ca+1
                   domain_lib%domain(i)%map(ca)=iatom-state%domain(i)%first_atom+1
                   end if
            end do
            
            !do iatom=1,domain_lib%domain(i)%natom_keep
            !ca=domain_lib%domain(i)%map(iatom)+state%domain(i)%first_atom-1
            !write(0,*) domain_lib%domain(i)%map,state%topo%atom_id(ca)
            !end do
    end do
    
    write(0,*) 'initialisation of domain maps done'
    !pause
end subroutine init_domain_maps


end module domains_module