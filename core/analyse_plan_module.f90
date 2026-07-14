module analysis_plan_module

use project_module
implicit none

!=========================================================
! ANALYSIS STEP
!=========================================================
type :: analysis_step_type

```
character(64) :: name = ''
```

end type analysis_step_type

!=========================================================
! ONE FILE TASK
!=========================================================
type :: task_type

```
integer :: istate = 0
integer :: irep   = 0
integer :: ifile  = 0

character(400) :: filename = ''

integer :: natom  = 0
integer :: nframe = 0

integer :: first_frame = 1
integer :: last_frame  = 0

logical :: active = .true.
```

end type task_type

!=========================================================
! ONE REPLICATE
!=========================================================
type :: replicate_plan_type

```
integer :: istate = 0
integer :: irep   = 0

integer :: ntask = 0

type(task_type), allocatable :: task(:)
```

end type replicate_plan_type

!=========================================================
! GLOBAL SETTINGS
!=========================================================
type :: analysis_settings_type

```
logical :: dry_run = .false.

integer :: ignore_first_frame = 1
integer :: ignore_last_frame  = 1

integer :: block_size = 2
```

end type analysis_settings_type

!=========================================================
! COMPLETE ANALYSIS PLAN
!=========================================================
type :: analysis_plan_type

```
type(analysis_settings_type) :: settings

integer :: nstep = 0
type(analysis_step_type), allocatable :: step(:)

integer :: nreplicate = 0
type(replicate_plan_type), allocatable :: replicate(:)
```

end type analysis_plan_type

contains

!=========================================================
! ADD ANALYSIS STEP
!=========================================================
subroutine add_step(plan, step_name)

```
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
```

end subroutine add_step

!=========================================================
! PRINT PLAN
!=========================================================
subroutine print_plan(plan)

```
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
```

end subroutine print_plan

end module analysis_plan_module
