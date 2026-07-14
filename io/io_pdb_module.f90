module io_pdb_module

use topology_module
implicit none

contains

subroutine count_atoms_pdb(pdb_name, natom)

character(*), intent(in) :: pdb_name
integer, intent(out)     :: natom

integer :: io
character(200) :: line
character(6) :: record

natom = 0

open(1, file=trim(pdb_name), status='old', action='read', iostat=io)

if (io /= 0) then
    write(0,*) 'error opening file:', trim(pdb_name)
    stop
end if

do
    read(1,'(A)', iostat=io) line
    if (io /= 0) exit

    record = line(1:6)

    if (record == 'END   ') exit
    if (record == 'REMARK') cycle

    if (record == 'ATOM  ' .or. record == 'HETATM') then
        natom = natom + 1
    end if

end do

close(1)

write(0,*) 'atoms found:', natom

end subroutine count_atoms_pdb


subroutine read_pdb(pdb_name, topo)

character(*), intent(in) :: pdb_name
type(topology_type), intent(inout) :: topo

integer :: io_read, i
real(4) :: x,y,z,occ,beta
integer(8) :: atom_number, res_number
character(6) :: remark
character(4) :: atom_specifier
character(3) :: acid
character(1) :: strain
character(100) :: headerline
integer :: natom

open(1,file=trim(pdb_name),status='old',action='read')

read(1,*) headerline

natom = topo%natom   ! <-- WICHTIG

do i=1,natom
   ! write(0,*) i,natom
    read(1,*) remark, atom_number, atom_specifier, acid, strain, &
              res_number, x, y, z, occ, beta

    topo%atom_id(i)     = atom_specifier
    topo%acid_id(i)     = acid
    topo%strain_id(i)   = strain
    topo%acid_number(i) = res_number

    topo%refx(i) = x
    topo%refy(i) = y
    topo%refz(i) = z

end do

close(1)

end subroutine read_pdb

end module io_pdb_module