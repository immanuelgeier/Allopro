module analysis_plan_module

use project_module
implicit none

!=========================================================
! ANALYSIS STEP
!=========================================================
type :: analysis_step_type


character(64) :: name = ''


end type analysis_step_type

!=========================================================
! ONE FILE TASK
!=========================================================
type :: task_type

integer :: istate = 0
integer :: irep   = 0
integer :: ifile  = 0

character(400) :: filename = ''

integer :: natom  = 0
integer :: nframe = 0

integer :: first_frame = 1
integer :: last_frame  = 0

logical :: active = .true.


end type task_type

!=========================================================
! ONE REPLICATE
!=========================================================
type :: replicate_plan_type

    integer :: istate
    integer :: irep

    integer :: nfile_task = 0
    type(file_task_type), allocatable :: file_task(:)

    integer :: nblock_task = 0
    type(block_task_type), allocatable :: block_task(:)

end type


!=========================================================
! GLOBAL SETTINGS
!=========================================================
type :: analysis_settings_type

    logical :: dry_run = .false.

    logical :: use_streaming = .true.

    integer :: ignore_first_frame = 1
    integer :: ignore_last_frame  = 1

    integer :: block_size = 32768

end type

!=========================================================
! COMPLETE ANALYSIS PLAN
!=========================================================
type :: analysis_plan_type


type(analysis_settings_type) :: settings

integer :: nstep = 0
type(analysis_step_type), allocatable :: step(:)

integer :: nreplicate = 0
type(replicate_plan_type), allocatable :: replicate(:)
  
integer :: mode = 1
    ! 1 = file-based
    ! 2 = streaming blocks

end type analysis_plan_type


type :: file_task_type

    integer :: istate
    integer :: irep
    integer :: ifile

    character(400) :: filename

    integer :: natom
    integer :: nframe

    integer :: first_frame
    integer :: last_frame

end type

type :: block_task_type

    integer :: iblock

    integer :: first_global_frame
    integer :: last_global_frame

    integer :: nframe

    integer :: istate
    integer :: irep

end type

contains

!=========================================================
! ADD ANALYSIS STEP
!=========================================================
subroutine add_step(plan, step_name)

implicit none

type(analysis_plan_type), intent(inout) :: plan
character(*), intent(in) :: step_name

type(analysis_step_type), allocatable :: tmp(:)
integer :: n

n = plan%nstep + 1

allocate(tmp(n))

if (plan%nstep > 0) then
    tmp(1:plan%nstep) = plan%step
    deallocate(plan%step)
end if

tmp(n)%name = trim(step_name)

call move_alloc(tmp, plan%step)

plan%nstep = n


end subroutine add_step

!=========================================================
! PRINT PLAN
!=========================================================
subroutine print_plan(plan)


implicit none

type(analysis_plan_type), intent(in) :: plan

integer :: i

write(0,*)
write(0,*) '================================'
write(0,*) 'ANALYSIS PLAN'
write(0,*) '================================'

write(0,*) 'dry run:', plan%settings%dry_run

write(0,*)
write(0,*) 'analysis steps:'

do i=1,plan%nstep
    write(0,*) '  ',trim(plan%step(i)%name)
end do

write(0,*)


end subroutine print_plan

end module analysis_plan_module
