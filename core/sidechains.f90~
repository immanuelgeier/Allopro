module sidechains

use topology_module

implicit none

private

public :: residue_cv_type
public :: sidechain_library_type
public :: test_sidechain_head_maps
public :: compute_freq_and_resids
public :: init_resnames_and_sidechain_head_groups
public :: init_sidechain_head
public :: make_sc_res_maps
public :: read_sidechain_library
public :: write_sidechain_library
public :: clear_sidechain_library

!=========================================================
type :: residue_cv_type
    
    logical :: active
    character(3) :: resname = ''

    integer :: natom_all  = 0
    integer :: natom_head = 0

    character(4) :: dir1_atom_1 = '    '
    character(4) :: dir1_atom_2 = '    '

    character(4) :: dir2_atom_1 = '    '
    character(4) :: dir2_atom_2 = '    '

    real(8) :: d1_ref(3) = 0.d0
    real(8) :: d2_ref(3) = 0.d0

    real(8), allocatable :: ref_xyz(:,:)
    real(8), allocatable :: cv(:,:),cv_id(:,:)

    integer, allocatable :: map(:)

    character(4), allocatable :: atom_names(:)
    integer :: freq=0
    integer,allocatable :: res_ids(:)
    logical             :: has_sym=.false.
    integer,allocatable :: perms(:,:),order(:)

end type

!=========================================================
type :: sidechain_library_type
    integer :: ncv
    type(residue_cv_type), allocatable :: residue(:)

end type

contains
subroutine clear_sidechain_library(sc_lib)
    implicit none
    type(sidechain_library_type), intent(inout) :: sc_lib
    integer :: i

    if (allocated(sc_lib%residue)) then
        do i = 1, size(sc_lib%residue)
            if (allocated(sc_lib%residue(i)%ref_xyz))    deallocate(sc_lib%residue(i)%ref_xyz)
            if (allocated(sc_lib%residue(i)%cv))         deallocate(sc_lib%residue(i)%cv)
            if (allocated(sc_lib%residue(i)%cv_id))      deallocate(sc_lib%residue(i)%cv_id)
            if (allocated(sc_lib%residue(i)%map))        deallocate(sc_lib%residue(i)%map)
            if (allocated(sc_lib%residue(i)%atom_names)) deallocate(sc_lib%residue(i)%atom_names)
            if (allocated(sc_lib%residue(i)%res_ids))    deallocate(sc_lib%residue(i)%res_ids)
            if (allocated(sc_lib%residue(i)%perms))      deallocate(sc_lib%residue(i)%perms)
            if (allocated(sc_lib%residue(i)%order))      deallocate(sc_lib%residue(i)%order)
        end do
        deallocate(sc_lib%residue)
    end if
    
    sc_lib%ncv = 0
end subroutine clear_sidechain_library



subroutine init_resnames_and_sidechain_head_groups(lib)

    implicit none

    type(sidechain_library_type), intent(inout) :: lib
   
    lib%ncv=19
    allocate(lib%residue(lib%ncv))
!=========================================================
! PHE
!=========================================================
lib%residue(1)%resname = 'PHE'
allocate(lib%residue(1)%atom_names(6))
lib%residue(1)%atom_names = (/ &
    'CG  ','CD1 ','CE1 ', &
    'CD2 ','CE2 ','CZ  ' /)
    
lib%residue(1)%dir1_atom_1="CG  "
lib%residue(1)%dir1_atom_2="CZ  "
lib%residue(1)%dir2_atom_1="CD1 "
lib%residue(1)%dir2_atom_2="CD2 "
lib%residue(1)%has_sym=.true.
allocate(lib%residue(1)%perms(6,1))
allocate(lib%residue(1)%cv_id(2,6))
lib%residue(1)%cv_id=0._8
lib%residue(1)%cv_id(1,1)=-1._8/3._8 !CG
lib%residue(1)%cv_id(1,6)= 1._8/3._8 !CZ

lib%residue(1)%cv_id(1,2)=-1._8/6._8 !CD1
lib%residue(1)%cv_id(1,3)= 1._8/6._8 !CD2

lib%residue(1)%cv_id(1,4)=-1._8/6._8 !CE1
lib%residue(1)%cv_id(1,5)= 1._8/6._8 !CE2

lib%residue(1)%cv_id(2,2)=-sqrt(3._8)/3._8 !CD1
lib%residue(1)%cv_id(2,4)= sqrt(3._8)/3._8 !CD2
lib%residue(1)%cv_id(2,3)=-sqrt(3._8)/3._8 !CE1
lib%residue(1)%cv_id(2,5)= sqrt(3._8)/3._8 !CDE2

!=========================================================
! HIS (HIE)
!=========================================================
lib%residue(2)%resname = 'HIE'
allocate(lib%residue(2)%atom_names(5))
lib%residue(2)%atom_names = (/ &
    'CG  ','ND1 ','CE1 ', &
    'CD2 ','NE2 ' /)
    
lib%residue(2)%dir1_atom_1="CG  "
lib%residue(2)%dir1_atom_2="NE2 " 
lib%residue(2)%dir2_atom_1="ND1 "
lib%residue(2)%dir2_atom_2="CD2 "

!=========================================================
! THR
!=========================================================
lib%residue(3)%resname = 'THR'
allocate(lib%residue(3)%atom_names(3))
lib%residue(3)%atom_names = (/ &
    'CB  ','OG1 ','CG2 ' /)

lib%residue(3)%dir1_atom_1="CB  "
lib%residue(3)%dir1_atom_2="OG1 "
lib%residue(3)%dir2_atom_1="CB  "
lib%residue(3)%dir2_atom_2="CG2 "

!=========================================================
! ARG
!=========================================================
lib%residue(4)%resname = 'ARG'
allocate(lib%residue(4)%atom_names(4))
lib%residue(4)%atom_names = (/ &
    'CZ  ','NE  ','NH1 ','NH2 ' /)

lib%residue(4)%dir1_atom_1="CZ  "
lib%residue(4)%dir1_atom_2="NE  "
lib%residue(4)%dir2_atom_1="NH1 "
lib%residue(4)%dir2_atom_2="NH2 "
lib%residue(4)%has_sym=.true.
allocate(lib%residue(4)%perms(4,1))
allocate(lib%residue(4)%cv_id(2,4))
lib%residue(4)%perms(1,1)=1
lib%residue(4)%perms(2,1)=2
lib%residue(4)%perms(3,1)=4
lib%residue(4)%perms(4,1)=3  
lib%residue(4)%cv_id=0._8
lib%residue(4)%cv_id(1,1)=-1._8/2._8
lib%residue(4)%cv_id(1,2)= 1._8/2._8
lib%residue(4)%cv_id(2,3)=-1._8/2._8
lib%residue(4)%cv_id(2,4)= 1._8/2._8
!=========================================================
! TYR
!=========================================================
lib%residue(5)%resname = 'TYR'
allocate(lib%residue(5)%atom_names(6))
lib%residue(5)%atom_names = (/ &
    'CG  ','CD1 ','CE1 ', &
    'CD2 ','CE2 ','CZ  ' /)

lib%residue(5)%dir1_atom_1="CG  "
lib%residue(5)%dir1_atom_2="CZ  "
lib%residue(5)%dir2_atom_1="CD1 "
lib%residue(5)%dir2_atom_2="CD2 " 
lib%residue(5)%has_sym=.true.
allocate(lib%residue(5)%perms(6,1))
lib%residue(5)%perms(1,1)=1
lib%residue(5)%perms(2,1)=4
lib%residue(5)%perms(3,1)=5
lib%residue(5)%perms(4,1)=2
lib%residue(5)%perms(5,1)=3
lib%residue(5)%perms(6,1)=6
allocate(lib%residue(5)%cv_id(2,6))
lib%residue(5)%cv_id=0._8
lib%residue(5)%cv_id(1,1)=-1._8/3._8 !CG
lib%residue(5)%cv_id(1,6)= 1._8/3._8 !CZ

lib%residue(5)%cv_id(1,2)=-1._8/6._8 !CD1
lib%residue(5)%cv_id(1,3)= 1._8/6._8 !CD2

lib%residue(5)%cv_id(1,4)=-1._8/6._8 !CE1
lib%residue(5)%cv_id(1,5)= 1._8/6._8 !CE2

lib%residue(5)%cv_id(2,2)=-sqrt(3._8)/3._8 !CD1
lib%residue(5)%cv_id(2,4)= sqrt(3._8)/3._8 !CD2
lib%residue(5)%cv_id(2,3)=-sqrt(3._8)/3._8 !CE1
lib%residue(5)%cv_id(2,5)= sqrt(3._8)/3._8 !CDE2


!=========================================================
! ASP
!=========================================================
lib%residue(6)%resname = 'ASP'
allocate(lib%residue(6)%atom_names(4))
lib%residue(6)%atom_names = (/ &
    'CB  ','CG  ','OD1 ','OD2 ' /)

lib%residue(6)%dir1_atom_1="CG  "
lib%residue(6)%dir1_atom_2="CB  "
lib%residue(6)%dir2_atom_1="OD1 "
lib%residue(6)%dir2_atom_2="OD2 "
lib%residue(6)%has_sym=.true.
allocate(lib%residue(6)%perms(4,1))
lib%residue(6)%perms(1,1)=1
lib%residue(6)%perms(2,1)=2
lib%residue(6)%perms(3,1)=4
lib%residue(6)%perms(4,1)=3
allocate(lib%residue(6)%cv_id(2,4))
lib%residue(6)%cv_id=0._8
lib%residue(6)%cv_id(1,1)=-1._8/2._8
lib%residue(6)%cv_id(1,2)= 1._8/2._8
lib%residue(6)%cv_id(2,3)=-1._8/2._8
lib%residue(6)%cv_id(2,4)= 1._8/2._8


!=========================================================
! ASN
!=========================================================
lib%residue(7)%resname = 'ASN'
allocate(lib%residue(7)%atom_names(4))
lib%residue(7)%atom_names = (/ &
    'CB  ','CG  ','OD1 ','ND2 ' /)
    
lib%residue(7)%dir1_atom_1="CB  "
lib%residue(7)%dir1_atom_2="CG "
lib%residue(7)%dir2_atom_1="OD1  "
lib%residue(7)%dir2_atom_2="ND2 "

!=========================================================
! GLN
!=========================================================
lib%residue(8)%resname = 'GLN'
allocate(lib%residue(8)%atom_names(4))
lib%residue(8)%atom_names = (/ &
    'CG  ','CD  ','OE1 ','NE2 ' /)

lib%residue(8)%dir1_atom_1="CD  "
lib%residue(8)%dir1_atom_2="OE1 " 
lib%residue(8)%dir2_atom_1="CD  "
lib%residue(8)%dir2_atom_2="NE2 "
           
!=========================================================
! GLU
!=========================================================
lib%residue(9)%resname = 'GLU'
allocate(lib%residue(9)%atom_names(4))
lib%residue(9)%atom_names = (/ &
    'CG  ','CD  ','OE1 ','OE2 ' /)

lib%residue(9)%dir1_atom_1="CD  "
lib%residue(9)%dir1_atom_2="OE1 "
lib%residue(9)%dir2_atom_1="CD  "
lib%residue(9)%dir2_atom_2="OE2 "
lib%residue(9)%has_sym=.true.
allocate(lib%residue(9)%perms(4,1))
lib%residue(9)%perms(1,1)=1
lib%residue(9)%perms(2,1)=2
lib%residue(9)%perms(3,1)=4
lib%residue(9)%perms(4,1)=3
allocate(lib%residue(9)%cv_id(2,4))
lib%residue(9)%cv_id=0._8
lib%residue(9)%cv_id(1,1)=-1._8/2._8
lib%residue(9)%cv_id(1,2)= 1._8/2._8
lib%residue(9)%cv_id(2,3)=-1._8/2._8
lib%residue(9)%cv_id(2,4)= 1._8/2._8


!=========================================================
! ILE
!=========================================================
lib%residue(10)%resname = 'ILE'
allocate(lib%residue(10)%atom_names(4))
lib%residue(10)%atom_names = (/ &
    'CB  ','CG1 ','CG2 ','CD  ' /)

lib%residue(10)%dir1_atom_1="CB  "
lib%residue(10)%dir1_atom_2="CD  "
lib%residue(10)%dir2_atom_1="CG1 "
lib%residue(10)%dir2_atom_2="CG2 "    

!=========================================================
! VAL
!=========================================================
lib%residue(11)%resname = 'VAL'
allocate(lib%residue(11)%atom_names(3))
lib%residue(11)%atom_names = (/ &
    'CB  ','CG1 ','CG2 ' /)
    
lib%residue(11)%dir1_atom_1="CB  "
lib%residue(11)%dir1_atom_2="CG1 "
lib%residue(11)%dir2_atom_1="CG1 "
lib%residue(11)%dir2_atom_2="CG2 "     

!=========================================================
! LYS
!=========================================================
lib%residue(12)%resname = 'LYS'
allocate(lib%residue(12)%atom_names(3))
lib%residue(12)%atom_names = (/ &
    'CD  ','CE  ','NZ  ' /)

lib%residue(12)%dir1_atom_1="CD  "
lib%residue(12)%dir1_atom_2="NZ  "
lib%residue(12)%dir2_atom_1="CD  "
lib%residue(12)%dir2_atom_2="CE  "      

!=========================================================
! LEU
!=========================================================
lib%residue(13)%resname = 'LEU'
allocate(lib%residue(13)%atom_names(4))
lib%residue(13)%atom_names = (/ &
    'CB  ','CG  ','CD1 ','CD2 ' /)
    
lib%residue(13)%dir1_atom_1="CB  "
lib%residue(13)%dir1_atom_2="CG  "
lib%residue(13)%dir2_atom_1="CD1 "
lib%residue(13)%dir2_atom_2="CD2 " 

!=========================================================
! MET
!=========================================================
lib%residue(14)%resname = 'MET'
allocate(lib%residue(14)%atom_names(4))
lib%residue(14)%atom_names = (/ &
    'CB  ','CG  ','SD  ','CE  ' /)

        
lib%residue(14)%dir1_atom_1="CB  "
lib%residue(14)%dir1_atom_2="CG  "
lib%residue(14)%dir2_atom_1="SD  "
lib%residue(14)%dir2_atom_2="CE  " 

!=========================================================
! CYS
!=========================================================
lib%residue(15)%resname = 'CYS'
allocate(lib%residue(15)%atom_names(3))
lib%residue(15)%atom_names = (/ &
    'CA  ','CB  ','SG  ' /)
    
lib%residue(15)%dir1_atom_1="CA  "
lib%residue(15)%dir1_atom_2="CB  "
lib%residue(15)%dir2_atom_1="CB  "
lib%residue(15)%dir2_atom_2="SG  " 


!=========================================================
! SER
!=========================================================
lib%residue(16)%resname = 'SER'
allocate(lib%residue(16)%atom_names(3))
lib%residue(16)%atom_names = (/ &
    'CA  ','CB  ','OG  ' /)

lib%residue(16)%dir1_atom_1="CB  "
lib%residue(16)%dir1_atom_2="OG  "
lib%residue(16)%dir2_atom_1="CA  "
lib%residue(16)%dir2_atom_2="CB  " 

!=========================================================
! PRO
!=========================================================
lib%residue(17)%resname = 'PRO'
allocate(lib%residue(17)%atom_names(5))
lib%residue(17)%atom_names = (/ &
    'N   ','CA  ','CB  ','CG  ','CD  ' /)

lib%residue(17)%dir1_atom_1="N   "
lib%residue(17)%dir1_atom_2="CG  "
lib%residue(17)%dir2_atom_1="CD  "
lib%residue(17)%dir2_atom_2="CG  " 

!=========================================================
! TRP
!=========================================================
lib%residue(18)%resname = 'TRP'
allocate(lib%residue(18)%atom_names(9))
lib%residue(18)%atom_names = (/ &
    'CG  ','CD1 ','NE1 ', &
    'CE2 ','CD2 ','CE3 ', &
    'CZ2 ','CZ3 ','CH2 ' /)
    
lib%residue(18)%dir1_atom_1='CD1 '
lib%residue(18)%dir1_atom_2='CZ3 '
lib%residue(18)%dir2_atom_1='CE2 '
lib%residue(18)%dir2_atom_2='CE3 '

!=========================================================
! HID
!=========================================================
lib%residue(19)%resname = 'HID'
allocate(lib%residue(19)%atom_names(5))
lib%residue(19)%atom_names = (/ &
    'CG  ','ND1 ','CE1 ', &
    'CD2 ','NE2 ' /)
    
lib%residue(19)%dir1_atom_1="CG  "
lib%residue(19)%dir1_atom_2="NE2 " 
lib%residue(19)%dir2_atom_1="ND1 "
lib%residue(19)%dir2_atom_2="CD2 "


end subroutine     



subroutine init_sidechain_head(lib)

    implicit none

    type(sidechain_library_type), intent(inout) :: lib

    integer :: i

    do i = 1, size(lib%residue)

        lib%residue(i)%natom_head = &
            size(lib%residue(i)%atom_names)
            

    end do

end subroutine

subroutine make_sc_res_maps(lib, istate, topologies)

    implicit none

    type(sidechain_library_type), intent(inout) :: lib
    integer, intent(in) :: istate
    type(topology_set_type), intent(in) :: topologies

    integer :: i,j,k,l
    integer :: ind_res,nres
    integer :: res_start_atom,natom_res

    character(3) :: resname
    character(3) :: cv_resname

    logical :: has_hd1,has_he2
    logical, allocatable :: his_used(:)

    nres = topologies%state(istate)%topo%nres

    allocate(his_used(nres))
    his_used = .false.

    do j = 1,size(lib%residue)

        resname = lib%residue(j)%resname
        
        if(resname == 'HID' .or. resname == 'HIE') resname = 'HIS'

        ind_res = -1

        do i = 2,nres-1

            !==================================================
            ! Normale Aminosäuren
            !==================================================
            if(resname /= 'HIS') then

                if(trim(topologies%state(istate)%topo%resnames(i)) == trim(resname)) then
                    ind_res = i
                    exit
                end if

            !==================================================
            ! Histidine
            !==================================================
            else

                if(trim(topologies%state(istate)%topo%resnames(i)) /= 'HIS') cycle
                if(his_used(i)) cycle

                res_start_atom = topologies%state(istate)%topo%acid_start_list(i)
                natom_res      = topologies%state(istate)%topo%acid_size(i)

                has_hd1 = .false.
                has_he2 = .false.

                do k = res_start_atom,res_start_atom+natom_res-1

                    if(trim(topologies%state(istate)%topo%atom_id(k)) == 'HD1') &
                        has_hd1 = .true.

                    if(trim(topologies%state(istate)%topo%atom_id(k)) == 'HE2') &
                        has_he2 = .true.

                end do

                if(has_hd1 .and. .not. has_he2) then
                    cv_resname = 'HID'
                else if(has_he2 .and. .not. has_hd1) then
                    cv_resname = 'HIE'
                else
                    cv_resname = 'HIP'
                end if

                if(trim(lib%residue(j)%resname) == trim(cv_resname)) then
                    ind_res = i
                    his_used(i) = .true.
                    exit
                end if

            end if

        end do
        


        if(ind_res < 0) cycle
        
  !f (trim(lib%residue(j)%resname) == 'PRO') then
  ! write(0,*) 'Library residue: ', j, trim(lib%residue(j)%resname)
  ! write(0,*) 'Selected topology residue: ', ind_res
  ! write(0,*) 'Topology residue name: ', &
  !     trim(topologies%state(istate)%topo%resnames(ind_res))
  !
  ! res_start_atom = topologies%state(istate)%topo%acid_start_list(ind_res)
  ! natom_res      = topologies%state(istate)%topo%acid_size(ind_res)
  !
  ! do k = res_start_atom,res_start_atom+natom_res-1
  !     write(0,*) trim(topologies%state(istate)%topo%atom_id(k))
  ! end do
  !end if

        res_start_atom = topologies%state(istate)%topo%acid_start_list(ind_res)
        natom_res      = topologies%state(istate)%topo%acid_size(ind_res)

        if(allocated(lib%residue(j)%map)) deallocate(lib%residue(j)%map)
        allocate(lib%residue(j)%map(lib%residue(j)%natom_head))
        lib%residue(j)%map = -1

        do k = res_start_atom,res_start_atom+natom_res-1
            do l = 1,lib%residue(j)%natom_head

                if(trim(adjustl(topologies%state(istate)%topo%atom_id(k))) == &
                   trim(adjustl(lib%residue(j)%atom_names(l)))) then

                    lib%residue(j)%map(l) = k-res_start_atom+1

                end if

            end do
        end do
    
   !if (trim(lib%residue(j)%resname) == 'PRO') then
   !write(0,*) 'PRO map:'
   !do l = 1,lib%residue(j)%natom_head
   !    write(0,*) trim(lib%residue(j)%atom_names(l)), &
   !               lib%residue(j)%map(l)
   !               
   !end do
   !end if

    end do


    deallocate(his_used)
    
    

do i=1,lib%ncv
if(allocated(lib%residue(i)%map)) then
	allocate(lib%residue(i)%order(lib%residue(i)%natom_head))
call argsort_int( lib%residue(i)%map,lib%residue(i)%order)
!write(0,*) lib%residue(i)%resname
!write(0,*)lib%residue(i)%map(:),'map'
!write(0,*) lib%residue(i)%order,'order'
end if
end do

!pause


end subroutine

subroutine argsort_int(a, idx)
   implicit none

   integer, intent(in)  :: a(:)
   integer, intent(out) :: idx(size(a))

   integer :: i, j, tmp, n

   n = size(a)

   ! Indexvektor initialisieren
   do i = 1, n
      idx(i) = i
   end do

   ! Einfacher Insertion Sort auf den Indizes
   do i = 2, n
      tmp = idx(i)
      j = i - 1

    do while (j >= 1)
       if (a(idx(j)) <= a(tmp)) exit
       idx(j+1) = idx(j)
       j = j - 1
    end do

idx(j+1) = tmp

      idx(j+1) = tmp
   end do

end subroutine argsort_int


logical function is_hydrogen(atomname)
    character(*), intent(in) :: atomname

    is_hydrogen = (atomname(1:1) == 'H')
end function




subroutine test_sidechain_head_maps(istate,lib,topologies)
implicit none
integer :: i,j,k,istate
integer ::res_start_atom,natom_res
type(sidechain_library_type), intent(in) :: lib
type(topology_set_type), intent(in) :: topologies

do i=1,size(lib%residue)

if(lib%residue(i)%freq==0) cycle

res_start_atom = &
    topologies%state(istate)%topo%acid_start_list(lib%residue(i)%res_ids(1))

natom_res = &
    topologies%state(istate)%topo%acid_size(lib%residue(i)%res_ids(1))
 
    
    write(0,*)i, lib%residue(i)%resname
    if(allocated(lib%residue(i)%map)) then
    	write(0,*) lib%residue(i)%atom_names
    	
        do j=1,size(lib%residue(i)%map) 
            write(0,*) lib%residue(i)%map(j),topologies%state(istate)%topo%atom_id(lib%residue(i)%map(j)+res_start_atom-1)
        end do
        write(0,*)
    end if
end do

end subroutine


subroutine compute_freq_and_resids(lib, istate, topologies)

    implicit none

    type(sidechain_library_type), intent(inout) :: lib
    integer, intent(in) :: istate
    type(topology_set_type), intent(in) :: topologies

    integer :: i,j,k
    integer :: nres
    integer :: res_start_atom,natom_res

    character(3) :: resname
    character(3) :: his_type

    logical :: has_hd1,has_he2
    logical, allocatable :: visited(:)

    nres = topologies%state(istate)%topo%nres

    allocate(visited(nres))
    visited = .false.

    do i = 1,size(lib%residue)

        resname = lib%residue(i)%resname

        lib%residue(i)%freq = 0

        if(allocated(lib%residue(i)%res_ids)) deallocate(lib%residue(i)%res_ids)
        allocate(lib%residue(i)%res_ids(nres))

        do j = 2,nres-1

            if(visited(j)) cycle

            !--------------------------------------------------
            ! normale Aminosäuren
            !--------------------------------------------------
            if(resname /= 'HID' .and. resname /= 'HIE') then

                if(trim(topologies%state(istate)%topo%resnames(j)) /= trim(resname)) cycle

            !--------------------------------------------------
            ! Histidine klassifizieren
            !--------------------------------------------------
            else

                if(trim(topologies%state(istate)%topo%resnames(j)) /= 'HIS') cycle

                res_start_atom = topologies%state(istate)%topo%acid_start_list(j)
                natom_res      = topologies%state(istate)%topo%acid_size(j)

                has_hd1 = .false.
                has_he2 = .false.

                do k = res_start_atom,res_start_atom+natom_res-1

                    if(trim(topologies%state(istate)%topo%atom_id(k)) == 'HD1') has_hd1 = .true.
                    if(trim(topologies%state(istate)%topo%atom_id(k)) == 'HE2') has_he2 = .true.

                end do

                if(has_hd1 .and. .not. has_he2) then
                    his_type = 'HID'
                else if(has_he2 .and. .not. has_hd1) then
                    his_type = 'HIE'
                else
                    his_type = 'HIP'
                end if

                if(trim(his_type) /= trim(resname)) cycle

            end if

            lib%residue(i)%freq = lib%residue(i)%freq + 1
            lib%residue(i)%res_ids(lib%residue(i)%freq) = j
            visited(j) = .true.

        end do

        if(lib%residue(i)%freq > 0) then
            lib%residue(i)%res_ids = &
                lib%residue(i)%res_ids(1:lib%residue(i)%freq)
        else
            deallocate(lib%residue(i)%res_ids)
        end if

    end do

    deallocate(visited)

end subroutine



subroutine write_sidechain_library(filename, lib)

    implicit none

    character(*), intent(in) :: filename
    type(sidechain_library_type), intent(in) :: lib

    integer :: u,i,j
    logical :: res_is_in_protein

    open(newunit=u,file=filename,status="replace",action="write")

    write(u,'(A,I0)') "NCV ",lib%ncv

    do i=1,lib%ncv
        write(u,'(A,I0)') "FREQ ",lib%residue(i)%freq
        if(lib%residue(i)%freq==0) cycle

        write(u,'(A,1X,A)') "RESNAME",trim(lib%residue(i)%resname)

        write(u,'(A,I0)') "NATOM_ALL ",lib%residue(i)%natom_all
        write(u,'(A,I0)') "NATOM_HEAD ",lib%residue(i)%natom_head
        
        write(u,'(A,L1)') "HAS_SYM ",lib%residue(i)%has_sym

        write(u,'(A,1X,A,1X,A)') "DIR1", &
             trim(lib%residue(i)%dir1_atom_1), &
             trim(lib%residue(i)%dir1_atom_2)

        write(u,'(A,1X,A,1X,A)') "DIR2", &
             trim(lib%residue(i)%dir2_atom_1), &
             trim(lib%residue(i)%dir2_atom_2)

        write(u,'(A)') "D1_REF"
        write(u,'(3F24.16)') lib%residue(i)%d1_ref

        write(u,'(A)') "D2_REF"
        write(u,'(3F24.16)') lib%residue(i)%d2_ref

        write(u,'(A)') "MAP"
        write(u,*) size(lib%residue(i)%map)
        write(u,*) lib%residue(i)%map

        write(u,'(A)') "RES_IDS"
        write(u,*) size(lib%residue(i)%res_ids)
        write(u,*) lib%residue(i)%res_ids

        write(u,'(A)') "REF_XYZ"
        write(u,*) size(lib%residue(i)%ref_xyz,2)

        do j=1,size(lib%residue(i)%ref_xyz,2)
            write(u,'(3F24.16)') lib%residue(i)%ref_xyz(:,j)
        end do

        if(allocated(lib%residue(i)%cv)) then
            write(u,'(A)') "CV"
            write(u,*) size(lib%residue(i)%cv,1),size(lib%residue(i)%cv,2)

            do j=1,size(lib%residue(i)%cv,2)
                write(u,'(3F24.16)') lib%residue(i)%cv(:,j)
            end do
        endif

        if(allocated(lib%residue(i)%cv_id)) then
            write(u,'(A)') "CV_ID"
            write(u,*) size(lib%residue(i)%cv_id,1),size(lib%residue(i)%cv_id,2)

            do j=1,size(lib%residue(i)%cv_id,2)
                write(u,'(3F24.16)') lib%residue(i)%cv_id(:,j)
            end do
        endif

        if(allocated(lib%residue(i)%order)) then
            write(u,'(A)') "ORDER"
            write(u,*) size(lib%residue(i)%order)
            write(u,*) lib%residue(i)%order
        endif

        if(allocated(lib%residue(i)%perms)) then
            write(u,'(A)') "PERMS"
            write(u,*) size(lib%residue(i)%perms,1),size(lib%residue(i)%perms,2)

            do j=1,size(lib%residue(i)%perms,2)
                write(u,*) lib%residue(i)%perms(:,j)
            end do
        endif

        write(u,'(A)') "END_RESIDUE"

    end do

    close(u)

end subroutine


subroutine read_sidechain_library(filename,lib)

    implicit none

    character(*),intent(in) :: filename
    type(sidechain_library_type),intent(out) :: lib
    logical :: res_is_in_protein     
    integer :: u,i,j,n,n1,n2
    character(32) :: key

    open(newunit=u,file=filename,status="old",action="read")

    read(u,*) key,lib%ncv
    allocate(lib%residue(lib%ncv))
    !read(u,*)
    do i=1,lib%ncv
        read(u,*) key,lib%residue(i)%freq
        if(lib%residue(i)%freq==0) cycle
        !read(u,*)
        !read(u,*) key
        

        read(u,*) key,lib%residue(i)%resname
        read(u,*) key,lib%residue(i)%natom_all
        read(u,*) key,lib%residue(i)%natom_head
        
        read(u,*) key,lib%residue(i)%has_sym

        read(u,*) key, &
             lib%residue(i)%dir1_atom_1, &
             lib%residue(i)%dir1_atom_2

        read(u,*) key, &
             lib%residue(i)%dir2_atom_1, &
             lib%residue(i)%dir2_atom_2

        read(u,*) key
        read(u,*) lib%residue(i)%d1_ref

        read(u,*) key
        read(u,*) lib%residue(i)%d2_ref

        read(u,*) key
        read(u,*) n
        allocate(lib%residue(i)%map(n))
        read(u,*) lib%residue(i)%map

        read(u,*) key
        read(u,*) n
        allocate(lib%residue(i)%res_ids(n))
        read(u,*) lib%residue(i)%res_ids

        read(u,*) key
        read(u,*) n
        allocate(lib%residue(i)%ref_xyz(3,n))

        do j=1,n
            read(u,*) lib%residue(i)%ref_xyz(:,j)
        end do

        read(u,*) key

        if(trim(key)=="CV") then

            read(u,*) n1,n2
            allocate(lib%residue(i)%cv(n1,n2))

            do j=1,n2
                read(u,*) lib%residue(i)%cv(:,j)
            end do

            read(u,*) key

        endif

        if(trim(key)=="CV_ID") then

            read(u,*) n1,n2
            allocate(lib%residue(i)%cv_id(n1,n2))

            do j=1,n2
                read(u,*) lib%residue(i)%cv_id(:,j)
            end do

            read(u,*) key

        endif

        if(trim(key)=="ORDER") then

            read(u,*) n
            allocate(lib%residue(i)%order(n))
            read(u,*) lib%residue(i)%order

            read(u,*) key

        endif

        if(trim(key)=="PERMS") then

            read(u,*) n1,n2
            allocate(lib%residue(i)%perms(n1,n2))

            do j=1,n2
                read(u,*) lib%residue(i)%perms(:,j)
            end do

            read(u,*) key

        endif

        ! key sollte jetzt END_RESIDUE sein

    end do

    close(u)

end subroutine
                
                   
end module sidechains