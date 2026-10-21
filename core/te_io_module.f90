module te_io_module



implicit none
real(8),allocatable ::delta_te_matrix(:,:)&
,anti_sym_te(:,:), te_matrix(:,:),te_shuflle_set(:,:,:),te_shuffle_mean(:,:),te_shuffle_std(:,:)
real(8),allocatable :: row_sum(:),col_sum(:),antisym_row_sum(:),antisym_col_sum(:)
character*1024, allocatable :: shuflle_names(:)
character*1024 :: raw_te_matrix_file
integer :: nshuflle,msize


contains

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


