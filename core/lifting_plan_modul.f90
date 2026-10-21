module lifting_plan_module
implicit  none

use project_module

type :: lifting_plan_type

integer :: istate
integer :: irep
integer :: nblock
integer :: size_block
integer :: ntimes
integer :: stride
logical, allocatable :: timestep_pattern(:)

end type


contains

subroutine init_lifiting_plan(irep,istate,stride,size_block,project)
implicit none
integer :: istate,stride,size_block

write(0,*) project%state(istate)%rep(irep)%n_total_use


end subroutine

end module lifting_plan_module