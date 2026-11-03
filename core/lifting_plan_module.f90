module lifting_plan_module

use project_module
use time_window_file_mapping
use time_window_loader
use iso_fortran_env, only: error_unit

implicit  none
private

public :: lifting_plan_type
public :: init_time_window_from_lifting_plan
public :: init_lifting_plan_all
public :: init_lifting_plan_state
public :: write_lifting_plan
public :: read_lifting_plan

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

subroutine write_lifting_plan(lifting_plan, filename)
   type(lifting_plan_type), intent(in) :: lifting_plan
   character(len=*), intent(in) :: filename
   integer :: unit, ios, nstates, nreplicas, nblocks, block_size

   call validate_lifting_plan(lifting_plan)
   nstates = size(lifting_plan%istates)
   nreplicas = size(lifting_plan%ireps, 2)
   nblocks = size(lifting_plan%start_block, 3)
   block_size = size(lifting_plan%timestep_pattern, 1)

   open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(error_unit, '(A,A)') 'Cannot write lifting plan file: ', trim(filename)
      error stop 'Lifting plan output failed'
   end if

   write(unit, '(A)', iostat=ios) '# Lifting plan data'
   if (ios == 0) write(unit, *, iostat=ios) 1, nstates, nreplicas, nblocks, block_size
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%istates
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%nreps
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%ireps
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%nblock
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%start_block
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%end_block
   if (ios == 0) write(unit, *, iostat=ios) lifting_plan%timestep_pattern
   if (ios /= 0) then
      close(unit)
      write(error_unit, '(A,A)') 'Error writing lifting plan file: ', trim(filename)
      error stop 'Lifting plan output failed'
   end if
   close(unit, iostat=ios)
   if (ios /= 0) error stop 'Error closing lifting plan file'
end subroutine write_lifting_plan

subroutine read_lifting_plan(lifting_plan, filename)
   type(lifting_plan_type), intent(out) :: lifting_plan
   character(len=*), intent(in) :: filename
   integer :: unit, ios, format_version, nstates, nreplicas, nblocks, block_size
   character(len=256) :: comment

   open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(error_unit, '(A,A)') 'Cannot read lifting plan file: ', trim(filename)
      error stop 'Lifting plan input failed'
   end if

   read(unit, '(A)', iostat=ios) comment
   if (ios /= 0) call lifting_plan_read_error(unit, 'Missing lifting plan header')
   read(unit, *, iostat=ios) format_version, nstates, nreplicas, nblocks, block_size
   if (ios /= 0) call lifting_plan_read_error(unit, 'Invalid lifting plan dimensions')
   if (format_version /= 1) call lifting_plan_read_error(unit, 'Unsupported lifting plan format version')
   if (nstates < 1 .or. nreplicas < 1 .or. nblocks < 0 .or. block_size < 1) then
      call lifting_plan_read_error(unit, 'Invalid lifting plan dimensions')
   end if

   allocate(lifting_plan%istates(nstates), lifting_plan%nreps(nstates))
   allocate(lifting_plan%ireps(nstates, nreplicas), lifting_plan%nblock(nstates, nreplicas))
   allocate(lifting_plan%start_block(nstates, nreplicas, nblocks))
   allocate(lifting_plan%end_block(nstates, nreplicas, nblocks))
   allocate(lifting_plan%timestep_pattern(block_size, nstates, nreplicas, nblocks))

   read(unit, *, iostat=ios) lifting_plan%istates
   if (ios == 0) read(unit, *, iostat=ios) lifting_plan%nreps
   if (ios == 0) read(unit, *, iostat=ios) lifting_plan%ireps
   if (ios == 0) read(unit, *, iostat=ios) lifting_plan%nblock
   if (ios == 0) read(unit, *, iostat=ios) lifting_plan%start_block
   if (ios == 0) read(unit, *, iostat=ios) lifting_plan%end_block
   if (ios == 0) read(unit, *, iostat=ios) lifting_plan%timestep_pattern
   if (ios /= 0) call lifting_plan_read_error(unit, 'Incomplete lifting plan data')
   close(unit, iostat=ios)
   if (ios /= 0) error stop 'Error closing lifting plan file'

   call validate_lifting_plan(lifting_plan)
end subroutine read_lifting_plan

subroutine validate_lifting_plan(lifting_plan)
   type(lifting_plan_type), intent(in) :: lifting_plan
   integer :: nstates, nreplicas, nblocks, block_size

   if (.not. allocated(lifting_plan%istates) .or. .not. allocated(lifting_plan%nreps) .or. &
       .not. allocated(lifting_plan%ireps) .or. .not. allocated(lifting_plan%nblock) .or. &
       .not. allocated(lifting_plan%start_block) .or. .not. allocated(lifting_plan%end_block) .or. &
       .not. allocated(lifting_plan%timestep_pattern)) then
      error stop 'Lifting plan is missing allocated arrays'
   end if

   nstates = size(lifting_plan%istates)
   if (nstates < 1) error stop 'Lifting plan has no states'
   if (size(lifting_plan%nreps) /= nstates) error stop 'Lifting plan state arrays have inconsistent sizes'
   if (size(lifting_plan%ireps, 1) /= nstates .or. size(lifting_plan%nblock, 1) /= nstates .or. &
       size(lifting_plan%start_block, 1) /= nstates .or. size(lifting_plan%end_block, 1) /= nstates .or. &
       size(lifting_plan%timestep_pattern, 2) /= nstates) then
      error stop 'Lifting plan arrays have inconsistent state dimensions'
   end if
   nreplicas = size(lifting_plan%ireps, 2)
   if (nreplicas < 1 .or. size(lifting_plan%nblock, 2) /= nreplicas .or. &
       size(lifting_plan%start_block, 2) /= nreplicas .or. &
       size(lifting_plan%end_block, 2) /= nreplicas .or. &
       size(lifting_plan%timestep_pattern, 3) /= nreplicas) then
      error stop 'Lifting plan arrays have inconsistent replicate dimensions'
   end if
   nblocks = size(lifting_plan%start_block, 3)
   if (size(lifting_plan%end_block, 3) /= nblocks .or. &
       size(lifting_plan%timestep_pattern, 4) /= nblocks) then
      error stop 'Lifting plan arrays have inconsistent block dimensions'
   end if
   block_size = size(lifting_plan%timestep_pattern, 1)
   if (block_size < 1) error stop 'Lifting plan has an invalid block size'
   if (any(lifting_plan%nreps < 0) .or. any(lifting_plan%nreps > nreplicas)) then
      error stop 'Lifting plan has invalid replicate counts'
   end if
   if (any(lifting_plan%nblock < 0) .or. any(lifting_plan%nblock > nblocks)) then
      error stop 'Lifting plan has invalid block counts'
   end if
end subroutine validate_lifting_plan

subroutine lifting_plan_read_error(unit, message)
   integer, intent(in) :: unit
   character(len=*), intent(in) :: message
   close(unit)
   write(error_unit, '(A)') trim(message)
   error stop 'Lifting plan input failed'
end subroutine lifting_plan_read_error

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
   
   if (istate < 1 .or. istate > project%nstate) then
      error stop 'Invalid state number for lifting plan'
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
   enddo
   lifting_plan%nreps(istate) = project%state(istate)%nrep

   !----------------------------------------
   ! fill block information
   !----------------------------------------

   do i=istate,istate

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