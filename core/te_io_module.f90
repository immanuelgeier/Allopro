module te_io_module

implicit none
private

public :: te_plot_data_type
public :: read_te_plot_data
public :: read_square_matrix_auto
public :: delta_te_matrix, anti_sym_te, te_matrix, te_shuflle_set, te_shuffle_mean
public :: te_shuffle_std, row_sum, col_sum, antisym_row_sum, antisym_col_sum
public :: shuflle_names, raw_te_matrix_file, nshuflle, msize, load_te_data

type :: te_plot_data_type
    integer :: matrix_size = 0
    integer :: n_shuffles = 0
    real(8), allocatable :: raw_matrix(:,:)
    real(8), allocatable :: shuffle_mean(:,:)
    real(8), allocatable :: delta_matrix(:,:)
    real(8), allocatable :: antisymmetric_matrix(:,:)
    real(8), allocatable :: row_sum(:)
    real(8), allocatable :: column_sum(:)
    real(8), allocatable :: antisymmetric_row_sum(:)
    real(8), allocatable :: antisymmetric_column_sum(:)
end type te_plot_data_type

real(8),allocatable ::delta_te_matrix(:,:)&
,anti_sym_te(:,:), te_matrix(:,:),te_shuflle_set(:,:,:),te_shuffle_mean(:,:),te_shuffle_std(:,:)
real(8),allocatable :: row_sum(:),col_sum(:),antisym_row_sum(:),antisym_col_sum(:)
character*1024, allocatable :: shuflle_names(:)
character*1024 :: raw_te_matrix_file
integer :: nshuflle,msize


contains

subroutine read_te_plot_data(raw_filename, shuffle_filenames, data)
character(len=*), intent(in) :: raw_filename
character(len=*), intent(in) :: shuffle_filenames(:)
type(te_plot_data_type), intent(out) :: data
real(8), allocatable :: matrix(:,:)
integer :: i, n

if (len_trim(raw_filename) == 0) error stop 'Raw TE matrix filename cannot be empty'
if (size(shuffle_filenames) < 1) error stop 'At least one shuffle matrix is required'

call read_square_matrix_auto(raw_filename, matrix)
n = size(matrix,1)
data%matrix_size = n
data%n_shuffles = size(shuffle_filenames)
allocate(data%raw_matrix(n,n), data%shuffle_mean(n,n), data%delta_matrix(n,n), &
    data%antisymmetric_matrix(n,n), data%row_sum(n), data%column_sum(n), &
    data%antisymmetric_row_sum(n), data%antisymmetric_column_sum(n))
data%raw_matrix = matrix
data%shuffle_mean = 0.0d0

do i = 1, size(shuffle_filenames)
    if (len_trim(shuffle_filenames(i)) == 0) then
        error stop 'Shuffle matrix filename cannot be empty'
    end if
    call read_square_matrix_auto(trim(shuffle_filenames(i)), matrix)
    if (size(matrix,1) /= n) then
        write(0,'(A,A,A,I0,A,I0)') 'Shuffle matrix "', trim(shuffle_filenames(i)), &
            '" has size ', size(matrix,1), ', expected ', n
        error stop 'TE matrix dimensions do not match'
    end if
    data%shuffle_mean = data%shuffle_mean + matrix
end do
data%shuffle_mean = data%shuffle_mean / real(size(shuffle_filenames),8)
data%delta_matrix = data%raw_matrix - data%shuffle_mean
data%antisymmetric_matrix = data%delta_matrix - transpose(data%delta_matrix)

do i = 1, n
    data%row_sum(i) = sum(data%delta_matrix(i,:)) - data%delta_matrix(i,i)
    data%column_sum(i) = sum(data%delta_matrix(:,i)) - data%delta_matrix(i,i)
    data%antisymmetric_row_sum(i) = sum(data%antisymmetric_matrix(i,:))
    data%antisymmetric_column_sum(i) = sum(data%antisymmetric_matrix(:,i))
end do
end subroutine read_te_plot_data

subroutine read_square_matrix_auto(filename, matrix)
character(len=*), intent(in) :: filename
real(8), allocatable, intent(out) :: matrix(:,:)
character(len=32768) :: line
integer :: unit, ios, rows, columns, fields, row

open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
if (ios /= 0) then
    write(0,'(A,A)') 'Cannot open TE matrix: ', trim(filename)
    error stop 'TE matrix input failed'
end if

rows = 0
columns = 0
do
    read(unit,'(A)',iostat=ios) line
    if (ios < 0) exit
    if (ios /= 0) then
        close(unit)
        error stop 'Error reading TE matrix dimensions'
    end if
    call strip_matrix_comment(line)
    if (len_trim(line) == 0) cycle
    fields = count_matrix_fields(line)
    if (fields < 1) cycle
    if (columns == 0) columns = fields
    if (fields /= columns) then
        close(unit)
        write(0,'(A,A)') 'Inconsistent number of columns in TE matrix: ', trim(filename)
        error stop 'TE matrix is not rectangular'
    end if
    rows = rows + 1
end do

if (rows < 1 .or. rows /= columns) then
    close(unit)
    write(0,'(A,A,A,I0,A,I0)') 'TE matrix "', trim(filename), '" is not square (', &
        rows, ' rows, ', columns, ' columns)'
    error stop 'TE matrix must be non-empty and square'
end if

allocate(matrix(rows,columns))
rewind(unit)
row = 0
do
    read(unit,'(A)',iostat=ios) line
    if (ios < 0) exit
    if (ios /= 0) then
        close(unit)
        error stop 'Error reading TE matrix values'
    end if
    call strip_matrix_comment(line)
    if (len_trim(line) == 0) cycle
    row = row + 1
    read(line,*,iostat=ios) matrix(row,:)
    if (ios /= 0) then
        close(unit)
        write(0,'(A,I0,A,A)') 'Invalid TE matrix values on row ', row, ' in ', trim(filename)
        error stop 'TE matrix contains invalid values'
    end if
end do
close(unit,iostat=ios)
if (ios /= 0) error stop 'Error closing TE matrix input'
end subroutine read_square_matrix_auto

subroutine strip_matrix_comment(line)
character(len=*), intent(inout) :: line
integer :: comment
comment = index(line,'#')
if (comment > 0) line(comment:) = ' '
end subroutine strip_matrix_comment

integer function count_matrix_fields(line) result(n_fields)
character(len=*), intent(in) :: line
integer :: i, char_code
logical :: in_field
n_fields = 0
in_field = .false.
do i = 1, len_trim(line)
    char_code = iachar(line(i:i))
    if (line(i:i) == ',' .or. line(i:i) == ' ' .or. char_code == 9) then
        in_field = .false.
    else if (.not. in_field) then
        n_fields = n_fields + 1
        in_field = .true.
    end if
end do
end function count_matrix_fields

subroutine load_te_data
integer:: i
msize=162
nshuflle=10
allocate(shuflle_names(nshuflle))
call init_file_names()!shuflle_names,raw_te_matrix_file)
do i=1,nshuflle
write(0,*)  trim(shuflle_names(i))

end do


call init_matrizes(msize,nshuflle)

allocate(delta_te_matrix(msize,msize))

te_shuffle_mean(:,:)=sum(te_shuflle_set,3)

delta_te_matrix=te_matrix-te_shuffle_mean/dble(nshuflle)


anti_sym_te=delta_te_matrix-transpose(delta_te_matrix)


end subroutine load_te_data



subroutine init_file_names!(shuflle_names,raw_te_matrix_file)
implicit none
!character*1024, allocatable ,intent(out):: shuflle_names(:)
!character*1024 ,intent(out):: raw_te_matrix_file

raw_te_matrix_file='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/C2_te_l100_te_matrix_saved.txt'
shuflle_names(1)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/C2_te_l100_te_matrix.txt'
shuflle_names(2)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s2_C2_te_l100_te_matrix.txt'
shuflle_names(3)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s3_C2_te_l100_te_matrix.txt'
shuflle_names(4)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s4_C2_te_l100_te_matrix.txt'
shuflle_names(5)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s5_C2_te_l100_te_matrix.txt'
shuflle_names(6)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s6_C2_te_l100_te_matrix.txt'
shuflle_names(7)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s7_C2_te_l100_te_matrix.txt'
shuflle_names(8)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s8_C2_te_l100_te_matrix.txt'
shuflle_names(9)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s9_C2_te_l100_te_matrix.txt'
shuflle_names(10)='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/s10_C2_te_l100_te_matrix.txt'

end subroutine init_file_names



subroutine init_matrizes(msize,nshuflle)
implicit none
integer ,intent(in)             :: msize,nshuflle
integer :: i

allocate(te_matrix(msize,msize))
allocate(anti_sym_te(msize,msize))
allocate(te_shuflle_set(msize,msize,nshuflle))
allocate(te_shuffle_mean(msize,msize) )
allocate(te_shuffle_std(msize,msize)  )
allocate(row_sum(msize)              )
allocate(col_sum(msize)              )
allocate(antisym_row_sum(msize)      )
allocate(antisym_col_sum(msize)      )


do i=1,nshuflle
call read_square_matrix(shuflle_names(i),msize,te_matrix)
te_shuflle_set(:,:,i)=te_matrix
end do
call read_square_matrix(raw_te_matrix_file,msize,te_matrix)

end subroutine init_matrizes

subroutine read_square_matrix(filename,msize,matrix)
integer ,intent(in)             :: msize
real(8),allocatable,intent(inout) :: matrix(:,:)
character(len=*),intent(in)      :: filename
integer :: i
open(1,file=trim(filename))

do i=1,msize
read(1,*)matrix(i,:)
end do
close(1)

end subroutine read_square_matrix

end module te_io_module
