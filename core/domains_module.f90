module domains_module

use project_module
use topology_module
implicit none

private


public :: init_domain_maps
public :: domain_library_type
public :: read_domain_library
public :: write_domain_library
public :: clear_domain_library
!=========================================================
type :: domain_cv_type
    

    logical :: active
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

subroutine clear_domain_library(dom_lib)
    implicit none
    type(domain_library_type), intent(inout) :: dom_lib
    integer :: i

    if (allocated(dom_lib%domain)) then
        do i = 1, size(dom_lib%domain)
            if (allocated(dom_lib%domain(i)%ref_xyz))    deallocate(dom_lib%domain(i)%ref_xyz)
            if (allocated(dom_lib%domain(i)%cv))         deallocate(dom_lib%domain(i)%cv)
            if (allocated(dom_lib%domain(i)%map))        deallocate(dom_lib%domain(i)%map)
            if (allocated(dom_lib%domain(i)%atom_names)) deallocate(dom_lib%domain(i)%atom_names)
            if (allocated(dom_lib%domain(i)%res_ids))    deallocate(dom_lib%domain(i)%res_ids)
        end do
        deallocate(dom_lib%domain)
    end if
    
    dom_lib%ndo = 0
end subroutine clear_domain_library
!=========================================================
!! BUILD RESIDUE + ATOM MAPS
!!=========================================================
!subroutine init_domain_maps(state,domain_lib,topologies)
!
!    implicit none
!
!    type(state_topology_type), intent(inout)  :: state
!    type(domain_library_type), intent(inout)  :: domain_lib
!    type(topology_set_type)  , intent(in)     :: topologies
!
!
!    integer :: i, j, res,iatom,ca
!    integer :: nres, natom,nres_do
!    integer :: first_res, last_res
!    integer :: astart, aend
!
!    nres  = state%topo%nres
!    natom = state%topo%natom
!
!    !====================================================
!    ! residue -> domain map
!    !====================================================
!    if (allocated(state%domain_key)) deallocate(state%domain_key)
!    allocate(state%domain_key(nres))
!    state%domain_key = 0
!    if(allocated(domain_lib%domain)) deallocate(domain_lib%domain)
!    allocate(domain_lib%domain(state%ndomain))
!    
!    domain_lib%ndo=state%ndomain
!    !====================================================
!    ! LOOP OVER DOMAINS
!    !====================================================
!    i = 1, state%ndomain
!      write(0,*) 'Domain:',state%domain(i)%name
!    !------------------------------------------------
!    ! DEFINE RESIDUE RANGE (THIS IS THE KEY FIX)
!    !------------------------------------------------
!    first_res = state%domain(i)%first_res
!    last_res  = state%domain(i)%last_res
!   
!    state%domain(i)%nres = last_res - first_res + 1
!   
!    !------------------------------------------------
!    ! allocate and fill resid
!    !------------------------------------------------
!    if (allocated(state%domain(i)%resid)) then
!        deallocate(state%domain(i)%resid)
!    end if
!   
!    allocate(state%domain(i)%resid(state%domain(i)%nres))
!   
!    do j = 1, state%domain(i)%nres
!        state%domain(i)%resid(j) = first_res + j - 1
!    end do
!    
!         if (allocated(domain_lib%domain(i)%res_ids)) then
!            deallocate(domain_lib%domain(i)%res_ids)
!        end if
!
!        allocate(domain_lib%domain(i)%res_ids(state%domain(i)%nres))
!
!        do j = 1, state%domain(i)%nres
!            domain_lib%domain(i)%res_ids(j) = first_res + j - 1
!        end do
!
!        !------------------------------------------------
!        ! fill domain_key
!        !------------------------------------------------
!        do j = 1, state%domain(i)%nres
!
!            res = state%domain(i)%resid(j)
!
!            if (res >= 1 .and. res <= nres) then
!                state%domain_key(res) = i
!            end if
!
!        end do
!
!    end do
!
!    !====================================================
!    ! ATOM BOUNDARIES FROM TOPOLOGY
!    !====================================================
!    do i = 1, state%ndomain
!        
!        domain_lib%domain(i)%domain_name=state%domain(i)%name
!
!        first_res = state%domain(i)%resid(1)
!        last_res  = state%domain(i)%resid(state%domain(i)%nres)
!
!        state%domain(i)%first_atom = &
!            state%topo%acid_start_list(first_res)
!        
!        domain_lib%domain(i)%start_atom=state%domain(i)%first_atom
!            
!        state%domain(i)%last_atom  = &
!            state%topo%acid_start_list(last_res) + &
!            state%topo%acid_size(last_res) - 1
!            
!        domain_lib%domain(i)%end_atom=state%domain(i)%last_atom
!            
!            domain_lib%domain(i)%natom_all&
!            =state%domain(i)%last_atom&
!            -state%domain(i)%first_atom+1
!            domain_lib%domain(i)%natom_keep&
!            =state%domain(i)%nres*3
!            domain_lib%domain(i)%nres_dom=state%domain(i)%nres
!            allocate(domain_lib%domain(i)%map(domain_lib%domain(i)%natom_keep))
!            ca=0
!            do iatom=state%domain(i)%first_atom,state%domain(i)%last_atom
!                   if(trim(state%topo%atom_id(iatom))=='N') then
!                   ca=ca+1
!                   domain_lib%domain(i)%map(ca)=iatom-state%domain(i)%first_atom+1
!                   end if
!                   if(trim(state%topo%atom_id(iatom))=='CA') then
!                   ca=ca+1
!                   domain_lib%domain(i)%map(ca)=iatom-state%domain(i)%first_atom+1
!                   end if
!                   if(trim(state%topo%atom_id(iatom))=='C') then
!                   ca=ca+1
!                   domain_lib%domain(i)%map(ca)=iatom-state%domain(i)%first_atom+1
!                   end if
!            end do
!            
!            !do iatom=1,domain_lib%domain(i)%natom_keep
!            !ca=domain_lib%domain(i)%map(iatom)+state%domain(i)%first_atom-1
!            !write(0,*) domain_lib%domain(i)%map,state%topo%atom_id(ca)
!            !end do
!    end do
!    
!    write(0,*) 'initialisation of domain maps done'
!    !pause
!end subroutine init_domain_maps
!
! ==========================================================
! 1. Topologie-Daten für Domänen im State finalisieren
! ==========================================================
!subroutine finalize_state_topology_domains(state)
!    implicit none
!
!    type(state_topology_type), intent(inout) :: state
!
!    integer :: i, j, res, iatom
!    integer :: nres, first_res, last_res
!
!    nres = state%topo%nres
!
!    ! Residue -> Domain Map initialisieren
!    if (allocated(state%domain_key)) deallocate(state%domain_key)
!    allocate(state%domain_key(nres))
!    state%domain_key = 0
!
!    ! Loop über alle Domänen zur Bestimmung der Residuenbereiche und Atomgrenzen
!    do i = 1, state%ndomain
!        write(0,*) 'Finalizing Domain Topology:', state%domain(i)%name
!
!        first_res = state%domain(i)%first_res
!        last_res  = state%domain(i)%last_res
!       
!        state%domain(i)%nres = last_res - first_res + 1
!       
!        if (allocated(state%domain(i)%resid)) deallocate(state%domain(i)%resid)
!        allocate(state%domain(i)%resid(state%domain(i)%nres))
!       
!        do j = 1, state%domain(i)%nres
!            state%domain(i)%resid(j) = first_res + j - 1
!        end do
!
!        ! domain_key füllen
!        do j = 1, state%domain(i)%nres
!            res = state%domain(i)%resid(j)
!            if (res >= 1 .and. res <= nres) then
!                state%domain_key(res) = i
!            end if
!        end do
!
!        ! Atomgrenzen festlegen
!        first_res = state%domain(i)%resid(1)
!        last_res  = state%domain(i)%resid(state%domain(i)%nres)
!
!        state%domain(i)%first_atom = state%topo%acid_start_list(first_res)
!        state%domain(i)%last_atom  = state%topo%acid_start_list(last_res) + &
!                                     state%topo%acid_size(last_res) - 1
!    end do
!
!    write(0,*) 'State topology domain finalization done'
!end subroutine finalize_state_topology_domains


! ==========================================================
! 2. Domain-Library aus finalisiertem State initialisieren
! ==========================================================
subroutine init_domain_maps(state, domain_lib)
    implicit none

    type(state_topology_type), intent(in)    :: state
    type(domain_library_type), intent(inout) :: domain_lib

    integer :: i, j, iatom, ca
    integer :: first_res, last_res

    if (allocated(domain_lib%domain)) deallocate(domain_lib%domain)
    allocate(domain_lib%domain(state%ndomain))
    
    domain_lib%ndo = state%ndomain

    do i = 1, state%ndomain
        !write(0,*) 'Initializing Domain Library Entry:', state%domain(i)%name

        domain_lib%domain(i)%domain_name = state%domain(i)%name
        domain_lib%domain(i)%nres_dom    = state%domain(i)%nres

        ! Residuen-IDs kopieren
        if (allocated(domain_lib%domain(i)%res_ids)) deallocate(domain_lib%domain(i)%res_ids)
        allocate(domain_lib%domain(i)%res_ids(state%domain(i)%nres))

        do j = 1, state%domain(i)%nres
            domain_lib%domain(i)%res_ids(j) = state%domain(i)%resid(j)
        end do

        ! Atom-Boundaries übernehmen
        domain_lib%domain(i)%start_atom  = state%domain(i)%first_atom
        domain_lib%domain(i)%end_atom    = state%domain(i)%last_atom
        domain_lib%domain(i)%natom_all   = state%domain(i)%last_atom - state%domain(i)%first_atom + 1
        domain_lib%domain(i)%natom_keep  = state%domain(i)%nres * 3

        ! Backbone-Map (N, CA, C) erstellen
        allocate(domain_lib%domain(i)%map(domain_lib%domain(i)%natom_keep))
        ca = 0
        do iatom = state%domain(i)%first_atom, state%domain(i)%last_atom
            if (trim(state%topo%atom_id(iatom)) == 'N' .or. &
                trim(state%topo%atom_id(iatom)) == 'CA' .or. &
                trim(state%topo%atom_id(iatom)) == 'C') then
                ca = ca + 1
                domain_lib%domain(i)%map(ca) = iatom - state%domain(i)%first_atom + 1
            end if
        end do
    end do

   ! write(0,*) 'Initialisation of domain library done'
end subroutine init_domain_maps

subroutine init_domain_maps_old(state,domain_lib,topologies)

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
end subroutine init_domain_maps_old

subroutine write_domain_library(filename, lib)

    character(*), intent(in) :: filename
    type(domain_library_type), intent(in) :: lib

    integer :: i,j,u,n

    open(newunit=u,file=filename,status="replace",action="write")

    write(u,'(A,I0)') "NDOMAINS ",lib%ndo

    do i=1,lib%ndo

        write(u,*)
        write(u,'(A)') "DOMAIN"

        write(u,'(A,1X,A)') "NAME",trim(lib%domain(i)%domain_name)

        write(u,'(A,I0)') "NATOM_ALL ",lib%domain(i)%natom_all
        write(u,'(A,I0)') "NATOM_KEEP ",lib%domain(i)%natom_keep
        write(u,'(A,I0)') "NRES ",lib%domain(i)%nres_dom

        write(u,'(A,I0)') "START_ATOM ",lib%domain(i)%start_atom
        write(u,'(A,I0)') "END_ATOM ",lib%domain(i)%end_atom

        write(u,'(A,1X,A,1X,A)') "DIR1", &
             trim(lib%domain(i)%dir1_atom_1), &
             trim(lib%domain(i)%dir1_atom_2)

        write(u,'(A,1X,A,1X,A)') "DIR2", &
             trim(lib%domain(i)%dir2_atom_1), &
             trim(lib%domain(i)%dir2_atom_2)

        write(u,'(A)') "D1_REF"
        write(u,'(3F24.16)') lib%domain(i)%d1_ref

        write(u,'(A)') "D2_REF"
        write(u,'(3F24.16)') lib%domain(i)%d2_ref
        write(u,'(A)') "MAP"
        n=size(lib%domain(i)%map)
        write(u,*) n
        write(u,*) (lib%domain(i)%map(j), j=1,n)
        write(u,'(A)') "RES_IDS"
        n=lib%domain(i)%nres_dom
        write(u,*) size(lib%domain(i)%res_ids)
        write(u,*) lib%domain(i)%res_ids
        write(u,'(A)') "REF_XYZ"
        write(u,*) size(lib%domain(i)%ref_xyz,2)

        do j=1,size(lib%domain(i)%ref_xyz,2)
            write(u,'(3F24.16)') lib%domain(i)%ref_xyz(:,j)
        end do

        
         if (allocated(lib%domain(i)%cv)) then
            write(u,'(A)') "CV"
            write(u,*) size(lib%domain(i)%cv,2)

            do j=1,size(lib%domain(i)%cv,2)
                write(u,'(3F24.16)') lib%domain(i)%cv(:,j)
            end do
        endif

                write(u,'(A)') "END_DOMAIN"
    end do

    close(u)

end subroutine
   
subroutine read_domain_library(filename, lib)

    character(*), intent(in)  :: filename
    type(domain_library_type), intent(out) :: lib

    integer :: u
    integer :: i,j,n,nresdom

    character(32) :: key

    open(newunit=u,file=filename,status="old",action="read")

    !-----------------------------
    ! Anzahl Domains
    !-----------------------------
    read(u,*) key, lib%ndo
    allocate(lib%domain(lib%ndo))

    do i=1,lib%ndo

        ! Leerzeile
        read(u,*)

        ! DOMAIN
        read(u,*) key
       ! write(0,*) key

        ! NAME
        read(u,*) key, lib%domain(i)%domain_name
!write(0,*) key
        ! NATOM_ALL
        read(u,*) key, lib%domain(i)%natom_all
!write(0,*) key
        ! NATOM_KEEP
        read(u,*) key, lib%domain(i)%natom_keep
!write(0,*) key
        ! NRES
        read(u,*) key, lib%domain(i)%nres_dom
!write(0,*) key
        ! START_ATOM
        read(u,*) key, lib%domain(i)%start_atom
!write(0,*) key
        ! END_ATOM
         read(u,*) key, lib%domain(i)%end_atom
!write(0,*) key
        ! DIR1
        read(u,*) key, lib%domain(i)%dir1_atom_1, &
                  lib%domain(i)%dir1_atom_2
!write(0,*) key
        ! DIR2
        read(u,*) key, lib%domain(i)%dir2_atom_1, &
                  lib%domain(i)%dir2_atom_2
!write(0,*) key
        ! D1_REF
        read(u,*) key
        read(u,*) lib%domain(i)%d1_ref
!write(0,*) key
        ! D2_REF
        read(u,*) key
        read(u,*) lib%domain(i)%d2_ref
!write(0,*) key
        !-----------------------------
        ! MAP
        !-----------------------------
        read(u,*) key
        read(u,*) n
!write(0,*) key
        allocate(lib%domain(i)%map(n))
        read(u,*) lib%domain(i)%map
        
        !-----------------------------
        ! RES_IDS
        !-----------------------------
        read(u,*) key
        read(u,*) nresdom
!write(0,*) n,nresdom,'n',lib%domain(i)%map,key
!pause
        allocate(lib%domain(i)%res_ids(nresdom))
        read(u,*) lib%domain(i)%res_ids !Bad integer


        !-----------------------------
        ! REF_XYZ
        !-----------------------------
        read(u,*) key
        read(u,*) n

        allocate(lib%domain(i)%ref_xyz(3,n))

        do j=1,n
            read(u,*) lib%domain(i)%ref_xyz(:,j)
        end do

        !-----------------------------
        ! CV
        !-----------------------------
        read(u,*) key

        if(trim(key)=="CV") then

            read(u,*) n

            allocate(lib%domain(i)%cv(2,n))

            do j=1,n
                read(u,*) lib%domain(i)%cv(:,j)
            end do

            read(u,*) key      ! END_DOMAIN

        else
            ! Falls keine CV gespeichert wurde,
            ! ist key bereits END_DOMAIN.
        endif

    end do

    close(u)

end subroutine

end module domains_module