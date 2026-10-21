module lifting_plan_module

use project_module
use time_window_file_mapping
use time_window_loader

implicit  none

type :: lifting_plan_type

integer,allocatable :: istates (:)  
integer,allocatable :: ireps   (:,:)
integer,allocatable :: nreps   (:)  
integer,allocatable :: nblock  (:,:) 
integer,allocatable :: start_block(:,:,:)
integer,allocatable :: end_block(:,:,:)
logical,allocatable::  timestep_pattern(:,:,:,:)

end type


contains

subroutine init_time_window_from_lifting_plan(lifting_plan,istate,irep,iblock,window,mask)
integer,intent(in) :: istate,irep,iblock
type(lifting_plan_type),intent(in) ::lifting_plan
type(time_window), intent(out) ::window
logical, allocatable, intent(out) :: mask(:)

window%istate         =lifting_plan%istates(istate)
window%irep           =lifting_plan%ireps(istate,irep)
window%start_time_step=lifting_plan%start_block(istate,irep,iblock)
window%end_time_step  =lifting_plan%end_block(istate,irep,iblock)

allocate(mask(window%end_time_step-window%start_time_step+1))
mask=lifting_plan%timestep_pattern(:,istate,irep,iblock)
end subroutine

subroutine init_lifting_plan_all(lifting_plan,offset,stride,size_block,project)

   implicit none
   type(project_type),     intent(in)  :: project
   type(lifting_plan_type),intent(out) :: lifting_plan
   integer, intent(in) :: offset
   integer, intent(in) :: stride
   integer, intent(in) :: size_block
   integer :: i,j,k
   integer :: max_rep,max_block
   integer :: nframes
   !----------------------------------------
   ! determine maximum dimensions
   !----------------------------------------

   max_rep   = 0
   max_block = 0
 
   do i=1,project%nstate

      max_rep = max(max_rep,project%state(i)%nrep)

      do j=1,project%state(i)%nrep

         nframes = project%state(i)%rep(j)%n_frames_use_total - offset

         if (nframes > 0) then
            max_block = max(max_block,nframes/size_block)
         endif

      enddo
   enddo

   !----------------------------------------
   ! allocate arrays
   !----------------------------------------

   allocate(lifting_plan%istates(project%nstate))
   allocate(lifting_plan%nreps(project%nstate))

   allocate(lifting_plan%ireps(project%nstate,max_rep))
   allocate(lifting_plan%nblock(project%nstate,max_rep))

   allocate(lifting_plan%start_block(project%nstate,max_rep,max_block))
   allocate(lifting_plan%end_block(project%nstate,max_rep,max_block))

   allocate(lifting_plan%timestep_pattern(size_block,project%nstate,max_rep,max_block))

   !----------------------------------------
   ! initialize
   !----------------------------------------

   lifting_plan%istates = 0
   lifting_plan%nreps   = 0
   lifting_plan%ireps   = 0
   lifting_plan%nblock  = 0

   lifting_plan%start_block = 0
   lifting_plan%end_block   = 0

   lifting_plan%timestep_pattern = .false.

   do i=1,project%nstate
      lifting_plan%istates(i) = i
      lifting_plan%nreps(i)   = project%state(i)%nrep
   enddo

   !----------------------------------------
   ! fill block information
   !----------------------------------------

   do i=1,project%nstate

      do j=1,project%state(i)%nrep

         lifting_plan%ireps(i,j) = j

         nframes = project%state(i)%rep(j)%n_frames_use_total - offset

         if (nframes <= 0) cycle

         lifting_plan%nblock(i,j) = nframes/size_block

         do k=1,lifting_plan%nblock(i,j)

            lifting_plan%start_block(i,j,k) = &
                 offset + (k-1)*size_block + 1

            lifting_plan%end_block(i,j,k) = &
                 offset + k*size_block

         enddo

      enddo

   enddo

   !----------------------------------------
   ! build stride mask
   !----------------------------------------
   do k=1,max_block
      do j=1,max_rep
         do i=1,project%nstate

            lifting_plan%timestep_pattern(:,i,j,k) = .false.
          
            lifting_plan%timestep_pattern(1:size_block:stride,i,j,k) = .true.

         enddo
      enddo
   enddo
  
end subroutine init_lifting_plan_all

subroutine init_lifting_plan_state(istate,lifting_plan,offset,stride,size_block,project)

   implicit none
   type(project_type),     intent(in)  :: project
   type(lifting_plan_type),intent(out) :: lifting_plan
   integer, intent(in) :: istate
   integer, intent(in) :: offset
   integer, intent(in) :: stride
   integer, intent(in) :: size_block
   integer :: i,j,k
   integer :: max_rep,max_block
   integer :: nframes
   !----------------------------------------
   ! determine maximum dimensions
   !----------------------------------------

   max_rep   = 0
   max_block = 0
   
   if(istate>1 .or. istate > project%nstate) then
   	   write(0,*) 'Invalid state number'
   	   write(0,*) 'Valid states >= 1 <= ',project%nstate
   end if
 
   do i=istate,istate

      max_rep = max(max_rep,project%state(i)%nrep)

      do j=1,project%state(i)%nrep

         nframes = project%state(i)%rep(j)%n_frames_use_total - offset

         if (nframes > 0) then
            max_block = max(max_block,nframes/size_block)
         endif

      enddo
   enddo

   !----------------------------------------
   ! allocate arrays
   !----------------------------------------

   allocate(lifting_plan%istates(project%nstate))
   allocate(lifting_plan%nreps(project%nstate))

   allocate(lifting_plan%ireps(project%nstate,max_rep))
   allocate(lifting_plan%nblock(project%nstate,max_rep))

   allocate(lifting_plan%start_block(project%nstate,max_rep,max_block))
   allocate(lifting_plan%end_block(project%nstate,max_rep,max_block))

   allocate(lifting_plan%timestep_pattern(size_block,project%nstate,max_rep,max_block))

   !----------------------------------------
   ! initialize
   !----------------------------------------

   lifting_plan%istates = 0
   lifting_plan%nreps   = 0
   lifting_plan%ireps   = 0
   lifting_plan%nblock  = 0

   lifting_plan%start_block = 0
   lifting_plan%end_block   = 0

   lifting_plan%timestep_pattern = .false.

   do i=1,project%nstate
      lifting_plan%istates(i) = i
      lifting_plan%nreps(i)   = project%state(i)%nrep
   enddo

   !----------------------------------------
   ! fill block information
   !----------------------------------------

   do i=1,project%nstate

      do j=1,project%state(i)%nrep

         lifting_plan%ireps(i,j) = j

         nframes = project%state(i)%rep(j)%n_frames_use_total - offset

         if (nframes <= 0) cycle

         lifting_plan%nblock(i,j) = nframes/size_block

         do k=1,lifting_plan%nblock(i,j)

            lifting_plan%start_block(i,j,k) = &
                 offset + (k-1)*size_block + 1

            lifting_plan%end_block(i,j,k) = &
                 offset + k*size_block

         enddo

      enddo

   enddo

   !----------------------------------------
   ! build stride mask
   !----------------------------------------
   do k=1,max_block
      do j=1,max_rep
         do i=1,project%nstate

            lifting_plan%timestep_pattern(:,i,j,k) = .false.
          
            lifting_plan%timestep_pattern(1:size_block:stride,i,j,k) = .true.

         enddo
      enddo
   enddo
  
end subroutine init_lifting_plan_state


end module