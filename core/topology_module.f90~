module topology_module

implicit none

!====================================================
! SINGLE TOPOLOGY
!====================================================
type :: topology_type


! atom level
integer :: natom = 0

character(4), allocatable :: atom_id(:)
character(3), allocatable :: acid_id(:)
character(1), allocatable :: strain_id(:)

integer, allocatable :: acid_number(:)
integer, allocatable :: res_key(:)

! residue level
integer :: nres = 0

integer, allocatable :: acid_size(:)
integer, allocatable :: acid_num_list(:)
integer, allocatable :: acid_start_list(:)

character(3), allocatable :: acid_id_list(:)
character(3), allocatable :: resnames(:)

! geometry
real(4), allocatable :: refx(:)
real(4), allocatable :: refy(:)
real(4), allocatable :: refz(:)

! mapping
integer, allocatable :: residx(:)


end type topology_type

!====================================================
! ONE TOPOLOGY PER STATE
!====================================================
type :: state_topology_type


character(400) :: pdb_file    = ''
character(400) :: domain_file = ''
integer :: ndomain = 0

type(topology_type) :: topo
type(domain_type), allocatable :: domain(:)

! residue -> domain map
integer, allocatable :: domain_key(:)

end type state_topology_type

!====================================================
! COLLECTION OF STATE TOPOLOGIES
!====================================================
type :: topology_set_type


integer :: nstate = 0

type(state_topology_type), allocatable :: state(:)


end type topology_set_type


type :: domain_type

character(40) :: name = ''

    integer :: nres = 0
    integer :: first_res = 0
    integer :: last_res  = 0
    integer :: first_atom = 0
    integer :: last_atom  = 0
integer, allocatable :: resid(:)

end type domain_type




contains

subroutine read_domains(filename,state)

implicit none

character(*),intent(in) :: filename
type(state_topology_type),intent(inout) :: state

integer :: i,j
integer :: first,last
character(40) :: name

open(10,file=filename,status='old',action='read')

read(10,*) state%ndomain

allocate(state%domain(state%ndomain))
allocate(state%domain_key(state%topo%nres))

state%domain_key = 0

do i=1,state%ndomain

   read(10,*) name, first, last

   state%domain(i)%name      = trim(name)
   state%domain(i)%first_res = first
   state%domain(i)%last_res  = last
   state%domain(i)%nres      = last-first+1

   allocate(state%domain(i)%resid(state%domain(i)%nres))

   do j=1,state%domain(i)%nres

      state%domain(i)%resid(j) = first+j-1

   end do

   !----------------------------------------
   ! residue -> domain map
   !----------------------------------------
   do j=1,state%topo%nres

      if(state%topo%acid_number(state%topo%acid_start_list(j)) >= first .and. &
         state%topo%acid_number(state%topo%acid_start_list(j)) <= last) then

         state%domain_key(j)=i

      end if

   end do

end do

close(10)

!----------------------------------------
! determine atom ranges
!----------------------------------------
do i=1,state%ndomain

   do j=1,state%topo%nres

      if(state%topo%acid_number(state%topo%acid_start_list(j)) == &
         state%domain(i)%first_res) then

         state%domain(i)%first_atom = state%topo%acid_start_list(j)

      end if

      if(state%topo%acid_number(state%topo%acid_start_list(j)) == &
         state%domain(i)%last_res) then

         state%domain(i)%last_atom = &
            state%topo%acid_start_list(j) + &
            state%topo%acid_size(j) - 1

      end if

   end do

end do

end subroutine


!====================================================
! ALLOCATE TOPOLOGY SET
!====================================================
subroutine allocate_topology_set(topologies, nstate)

type(topology_set_type), intent(inout) :: topologies
integer, intent(in) :: nstate

topologies%nstate = nstate

if (allocated(topologies%state)) then
    deallocate(topologies%state)
end if

allocate(topologies%state(nstate))


end subroutine allocate_topology_set

!====================================================
! ALLOCATE ATOM LEVEL ARRAYS
!====================================================
subroutine allocate_topology(topo)

type(topology_type), intent(inout) :: topo

allocate(topo%atom_id(topo%natom))
allocate(topo%acid_id(topo%natom))
allocate(topo%strain_id(topo%natom))

allocate(topo%acid_number(topo%natom))
allocate(topo%res_key(topo%natom))

allocate(topo%refx(topo%natom))
allocate(topo%refy(topo%natom))
allocate(topo%refz(topo%natom))


end subroutine allocate_topology

!====================================================
! BUILD RESIDUE SIZE + KEY MAP
!====================================================
subroutine make_acid_size_list(topo)

type(topology_type), intent(inout) :: topo

integer :: i
integer :: res
integer :: last_res
integer :: count

if (allocated(topo%acid_size))       deallocate(topo%acid_size)
if (allocated(topo%acid_num_list))   deallocate(topo%acid_num_list)
if (allocated(topo%acid_id_list))    deallocate(topo%acid_id_list)
if (allocated(topo%acid_start_list)) deallocate(topo%acid_start_list)
if (allocated(topo%resnames))        deallocate(topo%resnames)

res = 1
last_res = topo%acid_number(1)

do i = 2, topo%natom

    if (topo%acid_number(i) /= last_res) then
        res = res + 1
        last_res = topo%acid_number(i)
    end if

end do

topo%nres = res

allocate(topo%acid_size(res))
allocate(topo%acid_start_list(res))
allocate(topo%resnames(res))

topo%acid_size = 0

res = 1
count = 1

last_res = topo%acid_number(1)

topo%res_key(1) = 1

do i = 2, topo%natom

    if (topo%acid_number(i) /= last_res) then

        topo%acid_size(res) = count

        res = res + 1
        count = 0

        last_res = topo%acid_number(i)

    end if

    count = count + 1

    topo%res_key(i) = res

end do

topo%acid_size(res) = count


end subroutine make_acid_size_list

!====================================================
! BUILD RESIDUE START LIST
!====================================================
subroutine make_acid_start_list(topo)

type(topology_type), intent(inout) :: topo

integer :: i
integer :: counter

counter = 1

do i = 1, topo%nres

    topo%acid_start_list(i) = counter

    counter = counter + topo%acid_size(i)

end do

do i = 1, topo%nres

    topo%resnames(i) = &
        topo%acid_id(topo%acid_start_list(i))

end do


end subroutine make_acid_start_list

!====================================================
! BUILD RESIDUE INDEX MAP
!====================================================
subroutine make_residx(topo)


type(topology_type), intent(inout) :: topo

integer :: i
integer :: j

if (allocated(topo%residx)) then
    deallocate(topo%residx)
end if

allocate(topo%residx(topo%natom))

topo%residx = 0

do i = 1, topo%nres

    do j = topo%acid_start_list(i), &
           topo%acid_start_list(i) + topo%acid_size(i) - 1

        topo%residx(j) = i

    end do

end do


end subroutine make_residx

!====================================================
! COMPLETE TOPOLOGY FINALIZATION
!====================================================
subroutine finalize_topology(topo)

type(topology_type), intent(inout) :: topo

call make_acid_size_list(topo)
call make_acid_start_list(topo)
call make_residx(topo)

end subroutine finalize_topology

end module topology_module
