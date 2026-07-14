module domain_refstruc

use project_module
use topology_module
use sidechains
use time_window_file_mapping
use time_window_loader
use domains_module




contains

subroutine ref_from_frames_domains(time,irep,istate,lib,topologies,project)
implicit none
integer :: istate,irep,time_step,i,j,k,nt,ntstep,time,ido
type(time_window):: window
type(window_map)  :: map
type(raw_data_block) :: block
type(project_type), intent(in) :: project
type(domain_library_type), intent(inout) ::lib

character*4          :: atom_name
logical,allocatable  :: window_mask(:) 
integer, allocatable :: ref_time_steps(:)
logical, allocatable :: atom_groups(:)
logical, allocatable :: activ_dom(:)
integer :: method
integer :: icv,iatom,nhead,na,ng,ires,ndom,idom
real(8),allocatable :: coords(:,:,:),struc(:,:) 
integer :: ncv_all,ncoords,ngroup
integer,allocatable :: g_start(:),g_size(:)
real(8) :: v11(3),v12(3),v21(3),v22(3),d1(3),d2(3)
type(topology_set_type) :: topologies




window%start_time_step=time
window%end_time_step=time
window%istate=1
window%irep  =1
ntstep=1
ndom=topologies%state(window%istate)%ndomain
!write(0,*) 'ndom',ndom
!pause
allocate(activ_dom(ndom))
!if(allocated(lib%domain)) deallocate(lib%domain)
activ_dom=.false.

activ_dom=.true.


!pause
!call print_window_map(window, map)
allocate(window_mask(ntstep))
window_mask=.true.
!write(0,*) 'size test',size(window_mask),window%end_time_step-window%start_time_step+1
!pause
call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)

    nt=count(window_mask)
    write(0,*) block%istate
    
    do ido=1,lib%ndo
        !write(0,*) lib%domain(ido)%domain_name
        !write(0,*) 'ido',ido,allocated(lib%domain(ido)%ref_xyz)
       
        	
        	
        	
        	
       if(allocated(lib%domain(ido)%ref_xyz)) deallocate(lib%domain(ido)%ref_xyz)	
        allocate(lib%domain(ido)%ref_xyz(3,lib%domain(ido)%natom_keep))
       if(allocated(lib%domain(ido)%cv)) deallocate(lib%domain(ido)%cv)
        allocate(lib%domain(ido)%cv(2,lib%domain(ido)%natom_keep)) 
        
        
        
        
        do i=1,lib%domain(ido)%natom_keep
             
           
            iatom=lib%domain(ido)%map(i)+lib%domain(ido)%start_atom-1
            lib%domain(ido)%ref_xyz(:,i)=block%xyz(:,iatom,1)
            !write(0,*) topologies%state(istate)%topo%atom_id(iatom)
        end do
        call direction_vecs_ref_inv_gyration&
        (lib%domain(ido)%natom_keep,lib%domain(ido)%ref_xyz,.true.,d1,d2,lib%domain(ido)%cv)
    end do
    
    deallocate(block%xyz)

end subroutine




subroutine direction_vecs_ref_inv_gyration(n,ref,use_axes,d1_orig,d2_orig,cvec)
implicit none

integer,intent(in) :: n
real(8),intent(in) :: d1_orig(3),d2_orig(3)
real(8)            :: d1(3),d2(3)
real(8),allocatable,intent(in)  :: ref(:,:)
real(8),allocatable,intent(inout) :: cvec(:,:)

real(8) :: cog(3), ri(3)
real(8) :: c(n), c2(n)
real(8) :: v1(3), v2(3)

real(8) :: G(3,3)
real(8) :: U(3,3), VT(3,3), S(3)
real(8) :: work(100)

real(8) :: x1(3), x2(3)
real(8) :: tmp1(3), tmp2(3)

integer :: info, lwork, i, j
logical :: use_axes

!do i=1,n
!write(0,*) ref(:,i)
!end do
!pause
!orthogonalisation of d2

! ----------------------------------------
! CENTROID
! ----------------------------------------
cog = sum(ref,2) / dble(n)

! ----------------------------------------
! GYRATION TENSOR
! ----------------------------------------
G = 0.d0

do i=1,n
   ri = ref(:,i) - cog
   do j=1,3
      G(j,1) = G(j,1) + ri(j)*ri(1)
      G(j,2) = G(j,2) + ri(j)*ri(2)
      G(j,3) = G(j,3) + ri(j)*ri(3)
   end do
end do

! ----------------------------------------
! SVD
! ----------------------------------------
lwork = 100
call dgesvd('A','A',3,3,G,3,S,U,3,VT,3,work,lwork,info)

if (info /= 0) then
   write(0,*) 'DGESVD failed, info=',info
   stop
end if

! ----------------------------------------
! build x = V S^-1 U^T d
! ----------------------------------------



if (use_axes) then

    ! Eigenvektoren des symmetrischen Gyrationstensors
    d1 = U(:,1)
    d2 = U(:,2)

else

    d1 = d1_orig
    d2 = d2_orig
    d2 = d2-dot_product(d1,d2)/dot_product(d1,d1)*d1


end if
! U^T d
tmp1 = matmul(transpose(U), d1)
tmp2 = matmul(transpose(U), d2)

! S^-1
do i=1,3
   if (S(i) > 1d-12) then
      tmp1(i) = tmp1(i) / S(i)
      tmp2(i) = tmp2(i) / S(i)
   else
      tmp1(i) = 0.d0
      tmp2(i) = 0.d0
   end if
end do

! V = VT^T
x1 = matmul(transpose(VT), tmp1)
x2 = matmul(transpose(VT), tmp2)

! ----------------------------------------
! PROJECT ON ATOMS
! ----------------------------------------
do i=1,n
   ri = ref(:,i) - cog

   c(i)  = dot_product(ri, x1)
   c2(i) = dot_product(ri, x2)
end do

! ----------------------------------------
! CENTERING
! ----------------------------------------
c  = c  - sum(c)/dble(n)
c2 = c2 - sum(c2)/dble(n)

! ----------------------------------------
! RECONSTRUCTION CHECK
! ----------------------------------------
v1 = 0.d0
v2 = 0.d0

do i=1,n
   ri = ref(:,i) - cog
   v1 = v1 + c(i)*ri
   v2 = v2 + c2(i)*ri
end do
!
!write(0,*) 'v1',v1
!write(0,*) 'v2',v2
!
!write(0,*) 'd1',d1
!write(0,*) 'd2',d2
!
!write(0,*)'dot(v1,d1)',dot_product(v1,d1/dot_product(d1,d1))&
!         ,'dot(v2,d2)',dot_product(v2,d2/dot_product(d2,d2))
!
!write(0,*) 'c1',c
!write(0,*) 'c2',c2
!
!pause

! ----------------------------------------
! OUTPUT
! ----------------------------------------
do i=1,n
   cvec(1,i) = c(i)
   cvec(2,i) = c2(i)
end do

end subroutine



end module domain_refstruc