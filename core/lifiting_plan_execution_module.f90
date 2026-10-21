module lifting_plan_execution_module

use lifting_plan_moduel
use position_lift_module
use direction_vectors_analyse_plan_module
use project_module
use trajectory_io
use time_window_loader

contains


subroutine execute_lifting_plan(cv_p,lifting_plan,project)
implicit none
type(lifting_plan_typ)   :: lifting_plan
type(time_window)        :: time_window
type(raw_data_block)     :: block
type(lifted_trajectory)  :: lifted_tra
type(cv_plan), intent(in):: cv_p
integer :: irep,istate,iblock
logical, allocatable :: mask(:)
character*800 :: filename
character*800 :: prefix

do istate=1,size(lifting_plan%istates)
do irep=1,lifting_plan%nreps(istate)
do iblock=1,lifting_plan%nblock(istate,irep)

call init_time_window_from_lifting_plan(lifting_plan,istate,irep,iblock,time_window,mask)

call load_raw_data_block(project, time_window, mask, block)
call lift_cv_plan(cv_p, block, lifted_tra)

call make_filename_lifted_tra(istate,irep,iblock,prefix,filename)

prefix=''
prefix='lifted_tra'
call write_trajectory(filename, traj)

end do
end do
end do


end subroutine











module lifting_execution_module