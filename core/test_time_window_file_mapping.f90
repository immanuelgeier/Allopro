module test_time_window_file_mapping

use project_module
use time_window_file_mapping

implicit none

contains

!=========================================================
subroutine print_window_map(window, map)

    implicit none

    type(time_window), intent(in) :: window
    type(window_map),  intent(in) :: map

    integer :: i

    write(0,*)
    write(0,*) '================================'
    write(0,*) 'WINDOW MAP'
    write(0,*) '================================'

    write(0,*) 'state :', window%istate
    write(0,*) 'rep   :', window%irep

    write(0,*)
    write(0,*) 'window start:', window%start_time_step
    write(0,*) 'window end  :', window%end_time_step

    write(0,*)
    write(0,*) 'first dcd:', map%istart_file
    write(0,*) 'last  dcd:', map%iend_file

    write(0,*)
    write(0,*) 'local start time_step in first file:', &
                 map%start_timestep_first_file

    write(0,*) 'local last time_step in last file:', &
                 map%end_timestep_last_file

    write(0,*)
    write(0,*) 'nsegment:', map%nsegment
    write(0,*)

    do i = 1, map%nsegment

        write(0,*) '--------------------------------'
        write(0,*) 'segment:', i
        write(0,*) 'state   :', map%segment(i)%istate
        write(0,*) 'rep     :', map%segment(i)%irep
        write(0,*) 'file    :', map%segment(i)%ifile
        write(0,*) 't_start :', map%segment(i)%start_time_step_global
        write(0,*) 't_end   :', map%segment(i)%end_time_step_global
        write(0,*) 'start frame in file:', map%segment(i)%start_frame_in_file
        write(0,*) 'end frame in file  :', map%segment(i)%end_frame_in_file
        write(0,*) 'total raw frames :', map%segment(i)%nframes_raw
        write(0,*) 'use all  :', map%segment(i)%use_all

    end do

end subroutine


!=========================================================
subroutine print_window_dcd_usage(project, window)

    implicit none

    type(project_type), intent(in) :: project
    type(time_window),  intent(in) :: window

    integer :: t
    integer :: ifile
    integer :: istep

    integer :: ndcd
    integer, allocatable :: counts(:)

    ndcd = project%state(window%istate)%rep(window%irep)%ndcd
    allocate(counts(ndcd))
    counts = 0

    do t = window%start_time_step, window%end_time_step

        call find_dcd( &
            project, &
            window%istate, &
            window%irep, &
            t, &
            ifile, istep)

        if (ifile > 0) counts(ifile) = counts(ifile) + 1

    end do

    write(0,*)
    write(0,*) '================================'
    write(0,*) 'WINDOW DCD USAGE'
    write(0,*) '================================'

    write(0,*) 'state :', window%istate
    write(0,*) 'rep   :', window%irep
    write(0,*)

    do ifile = 1, ndcd

        if (counts(ifile) > 0) then
            write(0,'(A,I6,A,I12)') &
                'dcd ', ifile, ' frames ', counts(ifile)
        end if

    end do

    write(0,*)

    deallocate(counts)

end subroutine


!=========================================================
subroutine test_window_map(project)

    implicit none

    type(project_type), intent(in) :: project

    type(time_window) :: window
    type(window_map)  :: map

    integer :: istate,i
    integer :: irep,iseg,ntstep,cstep
    logical, allocatable :: time_mask(:),frame_to_keep(:)

    write(0,*)
    write(0,*) '================================'
    write(0,*) 'WINDOW MAP TEST'
    write(0,*) '================================'

    do istate = 1, min(1, project%nstate)

        do irep = 1, min(1, project%state(istate)%nrep)

            window%istate = istate
            window%irep   = irep

            window%start_time_step = 1
            window%end_time_step   = 100000
            ntstep=window%end_time_step-window%start_time_step+1
           
            allocate(time_mask(ntstep))
            
            time_mask=.false.
            do i=1,ntstep
                if(mod(i,10)==0) time_mask(i)=.true.
            end do
            
            call get_window_map(window, project, map)

            call print_window_map(window, map)

            call print_window_dcd_usage(project, window)
            
            cstep=0
            do iseg=1,map%nsegment
            call build_frame_to_keep_mask(iseg,project,window,map,time_mask,frame_to_keep)
            cstep=cstep+count(frame_to_keep)
            write(0,*) 'iseg',iseg
            write(0,*) 'size frame_to_keep_mask' ,size(frame_to_keep)
            write(0,*) 'count frame_to_keep_mask',count(frame_to_keep)
            write(0,*) 'consume complete file',map%segment(iseg)%use_all
            write(0,*)
            end do
            write(0,*) 'total number of timestep',cstep
            if(cstep/=count(time_mask)) then
            	write(0,*) 'miss match with numer of timesteps'
            	stop
            end if
        end do

    end do

end subroutine

!=========================================================
end module test_time_window_file_mapping