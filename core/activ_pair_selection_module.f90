module activ_pair_selection_module
use te_io_module
use direction_vectors_analyse_plan_module
use sparse_lut_module
implicit  none

contains







subroutine get_reduced_matrix(cv_p,size,matrix_active_final)
real(8),dimension(size/2,size/2) :: mred,sred,sqred,matrix_red
integer :: size,i,j
integer, allocatable :: order(:)

type(cv_plan), intent(in) :: cv_p
real    :: matrix_sorted(size,size)
logical,allocatable :: matrix_active(:,:)
logical,allocatable,intent(inout)::matrix_active_final(:,:)

call load_te_data()

call get_cv_residue_order(cv_p,order)


!do i=1,size
!do j=1,size
!matrix_sorted(order(i),order(j))=delta_te_matrix(i,j)
!end do
!end do

do i=1,size
do j=1,size
matrix_sorted(i,j)=delta_te_matrix(order(i),order(j))
end do
end do
allocate(matrix_active(size,size))


delta_te_matrix=matrix_sorted



!call maxval_group_reduction(size,delta_te_matrix,mred)
!call ranked_pairs_sym_mat(300,cv_p,size/2,mred,matrix_red,matrix_active)
open(1,file='delta_te_max_sorted.txt')
do i=1,size
write(1,*) delta_te_matrix(i,:)
end do
close(1)
!!do i=1,size/2
!!write(0,*) count(matrix_active(i,:)),'häufigkeit'
!!end do
!!pause
!do i=1,size/2
!write(1,*) matrix_red(i,:)
!end do
!close(1)

call sum_group_reduction(size,delta_te_matrix,sred)
call maxval_group_reduction(size,delta_te_matrix,mred)
call ranked_pairs_sym_mat(300,cv_p,size/2,sred,matrix_red,matrix_active)

write(0,*) count(matrix_active)
call add_new_random_pairs(size/2,matrix_active,100)
write(0,*) shape(matrix_active_final)

matrix_active_final=matrix_active
!do i=1,size/2
!do j=1,size/2
!matrix_active_final(2*i+1:2*i,2*j-1:2*j)=matrix_active(i,j)
!matrix_active_final(2*j+1:2*j,2*i-1:2*i)=matrix_active(i,j)
!end do
!end do


write(0,*) shape(matrix_active_final)
write(0,*) count(matrix_active_final),'nach routine'
pause
open(1,file='top_delta_te_sum.txt')
do i=1,size/2
write(1,*) matrix_red(i,:)
end do
close(1)

!stop



!return
call sort_pairs_sym_mat(cv_p,size/2,sred)

call sum_squares_group_reduction(size,delta_te_matrix,sqred)

call sum_group_reduction(size,delta_te_matrix,sred)
call maxval_group_reduction(size,delta_te_matrix,mred)
call ranked_pairs_sym_mat(300,cv_p,size/2,sqred,matrix_red,matrix_active)


open(1,file='top_delta_te_sum_sq.txt')
do i=1,size/2
write(1,*) matrix_red(i,:)
end do
close(1)




open(1,file='delta_te_matrix.txt')
do i=1,size
write(1,*) real(delta_te_matrix(i,:),4)
end do
close(1)

open(1,file='delta_te_matrix_antisym.txt')
do i=1,size
write(1,*) real(anti_sym_te(i,:),4)
end do
close(1)

open(1,file='delta_te_matrix_antisym_rsd.txt')
do i=1,size-4
write(1,*) i,real(anti_sym_te(i,i+1),4),real(anti_sym_te(i,i+2),4)&
,real(anti_sym_te(i,i+3),4),real(anti_sym_te(i,i+4),4),sum(anti_sym_te(:,i))
end do
close(1)

end subroutine get_reduced_matrix

subroutine add_new_random_pairs(hs,matrix_active,nnew)
integer, intent(in):: nnew,hs
logical,allocatable, intent(inout):: matrix_active(:,:)
integer::i,j,nfound
real(8):: rn

nfound=0
do while(nfound<nnew)
call random_number(rn)
i=floor(dble(hs)*rn)+1
call random_number(rn)
j=floor(dble(hs)*rn)+1
if(i==j) cycle

if(matrix_active(i,j) .or. matrix_active(j,i)) then
  cycle
else
  nfound=nfound+1	
  matrix_active(i,j)=.true.
  matrix_active(j,i)=.true.
end if
end do

write(0,*) count(matrix_active),'in routine'
pause
end subroutine add_new_random_pairs


subroutine ranked_pairs_sym_mat(rank,cv_p,hs,matrix,matrix_red,matrix_active)
integer, intent(in)   :: hs
real(8), intent(in) :: matrix(hs,hs)
real(8), intent(out) ::matrix_red(hs,hs)
logical,allocatable, intent(inout) ::matrix_active(:,:)
real(8),allocatable :: vec(:)
integer :: npair,ipair,i,j,k,ml,rank
integer,allocatable :: key(:,:),perm(:)
type(cv_plan), intent(in) :: cv_p
npair=(hs-1)*hs
npair=npair/2

allocate(vec(npair))
allocate(key(2,npair))
allocate(perm(npair))
ipair=0
do i=1,hs
do j=i+1,hs
ipair=ipair+1
vec(ipair)=matrix(i,j)
key(1,ipair)=i
key(2,ipair)=j
end do
end do
npair=ipair
 do ipair=1,npair
    ml=maxloc(vec,1)
    perm(ipair)=ml
    vec(ml)=-huge(1.)
  end do

do ipair=1,npair
i=key(1,perm(ipair))
j=key(2,perm(ipair))
write(0,*) ipair,i,j,matrix(i,j),trim(cv_p%object(i)%name(mod(i,2)+1)),trim(cv_p%object(j)%name(mod(j,2)+1))

end do
  
write(0,*) minloc(matrix)
write(0,*) maxloc(matrix)
matrix_active=.false.
matrix_red=0._8
do ipair=1,rank
i=key(1,perm(ipair))
j=key(2,perm(ipair))
matrix_red(i,j)=matrix(i,j)
matrix_red(j,i)=matrix(i,j)
matrix_active(i,j)=.true.
matrix_active(j,i)=.true.
end do

end subroutine ranked_pairs_sym_mat

subroutine sort_pairs_sym_mat(cv_p,hs,matrix)
integer, intent(in)   :: hs
real(8), intent(in) :: matrix(hs,hs)
real(8),allocatable :: vec(:)
integer :: npair,ipair,i,j,k,ml
integer,allocatable :: key(:,:),perm(:)
type(cv_plan), intent(in) :: cv_p
npair=(hs-1)*hs
npair=npair/2

allocate(vec(npair))
allocate(key(2,npair))
allocate(perm(npair))
ipair=0
do i=1,hs
do j=i+1,hs
ipair=ipair+1
vec(ipair)=matrix(i,j)
key(1,ipair)=i
key(2,ipair)=j
end do
end do
npair=ipair
 do ipair=1,npair
    ml=maxloc(vec,1)
    perm(ipair)=ml
    vec(ml)=-huge(1.)
  end do

open(1,file='delta_te_order.txt')  
do ipair=1,npair
i=key(1,perm(ipair))
j=key(2,perm(ipair))

write(0,*) ipair,i,j,matrix(i,j),trim(cv_p%object(i)%name(mod(i,2)+1)),trim(cv_p%object(j)%name(mod(j,2)+1))
write(1,*) ipair,i,j,matrix(i,j),trim(cv_p%object(i)%name(mod(i,2)+1)),trim(cv_p%object(j)%name(mod(j,2)+1))

end do
close(1)

write(0,*) minloc(matrix)
write(0,*) maxloc(matrix)



end subroutine sort_pairs_sym_mat



subroutine maxval_group_reduction(size,matrix,matrix_red)
implicit none

integer, intent(in) :: size
real(8), intent(in) :: matrix(size,size)
real(8), intent(out):: matrix_red(size/2,size/2)
real(8) :: block(2,2,2)
integer :: i,j

do i=1,size/2
do j=1,size/2
block(1,:,:)=matrix(2*i-1:2*i,2*j-1:2*j)
block(2,:,:)=matrix(2*j-1:2*j,2*i-1:2*i)
matrix_red(i,j)=maxval(block)
end do
end do

open(1,file='maxval_group_te_matrix.txt')
do i=1,size/2
write(1,*) real(matrix_red(i,:),4)
end do
close(1)

end subroutine maxval_group_reduction


subroutine sum_group_reduction(size,matrix,matrix_red)
implicit none

integer, intent(in) :: size
real(8), intent(in) :: matrix(size,size)
real(8), intent(out):: matrix_red(size/2,size/2)
real(8) :: block(2,2,2)
integer :: i,j

do i=1,size/2
do j=1,size/2
block(1,:,:)=matrix(2*i-1:2*i,2*j-1:2*j)
block(2,:,:)=matrix(2*j-1:2*j,2*i-1:2*i)

matrix_red(i,j)=sum(block)
end do
end do

open(1,file='sum_group_te_matrix.txt')
do i=1,size/2
write(1,*) real(matrix_red(i,:),4)
end do
close(1)

end subroutine sum_group_reduction

subroutine sum_squares_group_reduction(size,matrix,matrix_red)
implicit none

integer, intent(in) :: size
real(8), intent(in) :: matrix(size,size)
real(8), intent(out):: matrix_red(size/2,size/2)
real(8) :: block(2,2,2)
integer :: i,j

do i=1,size/2
do j=1,size/2
block(1,:,:)=matrix(2*i-1:2*i,2*j-1:2*j)
block(2,:,:)=matrix(2*j-1:2*j,2*i-1:2*i)
block=block*block
matrix_red(i,j)=sum(block)
end do
end do



open(1,file='sum_squared_group_te_matrix.txt')
do i=1,size/2
write(1,*) real(matrix_red(i,:),4)
end do
close(1)

end subroutine sum_squares_group_reduction


subroutine matrix_freshold(size,matrix,x)
integer, intent(in) :: size
real(8), intent(inout) ::x
real(8), intent(inout) :: matrix(size,size)
integer:: i,j
do i=1,size
do j=1,size
if(delta_te_matrix(i,j)<x) delta_te_matrix(i,j)=0._8
end do
end do
end subroutine matrix_freshold






end module activ_pair_selection_module