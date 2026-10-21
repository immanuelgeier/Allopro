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

subroutine write_subgroup_pdb(outfile,topo,pos,occ,beta,group)
    implicit none

    ! --- Inputs ---
    type(topology_type), intent(in) :: topo
    character(len=*), intent(in) :: outfile
    logical, intent(in) :: group(:)

    ! --- Locals ---
    integer :: i, natom_total, atom_counter
    integer(8) :: atom_number, res_number
    real(4) :: x, y, z
    real(4),intent(in) ::  occ(:), beta(:)
    character(6) :: remark
    character(4) :: atom_specifier
    character(3) :: acid
    character(1) :: strain
    character(100) :: format_pdb='(A6,I5,1X,A4,1X,A3,1X,A1,I4,3X,3F8.3,2F6.2,10X,A2,A2)'
    real(4), intent(in) :: pos(:,:)

    natom_total = size(group)

    open(10, file=trim(outfile), status='replace', action='write')

    ! Optional Header
    write(10,'(A)') 'REMARK   Subgroup PDB generated'

    atom_counter = 0

    do i = 1, natom_total
        if (group(i)) then
            atom_counter = atom_counter + 1

            remark = 'ATOM  '
            atom_number   = atom_counter
            atom_specifier = topo%atom_id(i)
            acid           = topo%acid_id(i)
            strain         = topo%strain_id(i)
            res_number     = topo%acid_number(i)

    !topo%atom_id(i)     = atom_specifier
    !topo%acid_id(i)     = acid
    !topo%strain_id(i)   = strain
    !topo%acid_number(i) = res_number
            
            
            x = pos(1,i)
            y = pos(2,i)
            z = pos(3,i)

            

            write(10,fmt=format_pdb) remark, atom_number, atom_specifier, &
                 acid, strain, res_number, x, y, z, occ(i), beta(i)
                 !remark,atom_number,atom_specifier,acid,strain,res_number,x,y,z,occ,beta
        end if
    end do

    write(10,'(A)') 'END'

    close(10)

    write(0,*) 'Subgroup PDB written:', trim(outfile), ' atoms:', atom_counter

end subroutine write_subgroup_pdb

end module io_pdb_module