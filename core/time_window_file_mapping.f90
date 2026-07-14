module time_window_file_mapping

use project_module
implicit none

!=========================================================
type :: time_window

    integer :: istate
    integer :: irep

    integer :: start_time_step
    integer :: end_time_step

end type

!=========================================================
type :: window_segment

    integer :: istate
    integer :: irep
    integer :: ifile

    ! local frame indices inside file
    integer :: start_frame_in_file
    integer :: end_frame_in_file

    ! global time indices (IMPORTANT for VAR/FFT alignment)
    integer :: start_time_step_global
    integer :: end_time_step_global

    integer :: nframes_raw
    integer :: nframes_use
    integer :: nframes_keep
    
    logical :: use_all
end type

!=========================================================
type :: window_map

    integer :: istart_file
    integer :: iend_file

    integer :: start_timestep_first_file
    integer :: end_timestep_last_file

    integer :: nsegment
    type(window_segment), allocatable :: segment(:)

end type


contains

!=========================================================
subroutine find_dcd( &
    project, &
    istate, &
    irep, &
    time_step, &
    ifile, &
    time_step_file )

    implicit none

    type(project_type), intent(in) :: project

    integer, intent(in)  :: istate
    integer, intent(in)  :: irep
    integer, intent(in)  :: time_step

    integer, intent(out) :: ifile
    integer, intent(out) :: time_step_file

    integer :: i
    integer :: tstart, tend

    ifile          = -1
    time_step_file = -1
    
    !write(0,*)istate,irep,'test'

    !pause
    !------------------------------------------------------
    ! sanity check
    !------------------------------------------------------
    if ( time_step < 1 .or. &
         time_step > project%state(istate)%rep(irep)%nframes_per_rep_use ) then

        write(0,*) 'find_dcd error'
        write(0,*) 'time_step:', time_step
        write(0,*) 'state    :', istate
        write(0,*) 'rep      :', irep
        stop
    end if

    !------------------------------------------------------
    ! global start shifted by offset
    !------------------------------------------------------
    tstart = 1 - project%state(istate)%rep(irep)%offset_start

    do i = 1, project%state(istate)%rep(irep)%ndcd

        tend = tstart + &
               project%state(istate)%rep(irep)%nframes_per_dcd_use(i) - 1

        if (time_step >= tstart .and. time_step <= tend) then

            ifile = i

            ! local frame index in file
            time_step_file = time_step - tstart + 1

            return
        end if

        tstart = tend + 1
    end do

    write(0,*) 'find_dcd internal error'
    write(0,*) 'time_step=', time_step
    stop

end subroutine


subroutine get_window_map(window, project, map)

    implicit none

    type(time_window),  intent(in)  :: window
    type(project_type), intent(in)  :: project
    type(window_map),   intent(out) :: map

    integer :: istate, irep
    integer :: i, ifile

    integer :: file_use_start, file_use_end
    integer :: seg_use_start, seg_use_end

    integer :: rel_start, rel_end
    integer :: raw_offset

    istate = window%istate
    irep   = window%irep

    !======================================================
    ! locate boundary files (USE space)
    !======================================================
    !write(0,*) istate,irep,'test2'
    !pause
    call find_dcd(project, istate, irep, &
        window%start_time_step, &
        map%istart_file, map%start_timestep_first_file)

    call find_dcd(project, istate, irep, &
        window%end_time_step, &
        map%iend_file, map%end_timestep_last_file)

    !======================================================
    ! allocate segments
    !======================================================
    map%nsegment = map%iend_file - map%istart_file + 1
    allocate(map%segment(map%nsegment))

    !======================================================
    ! build segments
    !======================================================
    do i = 1, map%nsegment

        ifile = map%istart_file + i - 1

        map%segment(i)%istate = istate
        map%segment(i)%irep   = irep
        map%segment(i)%ifile  = ifile

        !--------------------------------------------------
        ! USE-space file boundaries
        !--------------------------------------------------
        call compute_file_global_range( &
            project, istate, irep, ifile, &
            file_use_start, file_use_end)

        !--------------------------------------------------
        ! clip in USE space
        !--------------------------------------------------
        seg_use_start = max(file_use_start, window%start_time_step)
        seg_use_end   = min(file_use_end,   window%end_time_step)

        !--------------------------------------------------
        ! relative position inside file (USE space)
        !--------------------------------------------------
        rel_start = seg_use_start - file_use_start
        rel_end   = seg_use_end   - file_use_start

        !--------------------------------------------------
        ! convert to RAW frame index
        !--------------------------------------------------
        raw_offset = project%state(istate)%rep(irep)%ignore_first_per_dcd(ifile)

        map%segment(i)%start_frame_in_file = raw_offset + 1 + rel_start
        map%segment(i)%end_frame_in_file   = raw_offset + 1 + rel_end

        map%segment(i)%nframes_raw = project%state(istate)%rep(irep)%nframes_per_dcd_raw(ifile)
        map%segment(i)%nframes_use = project%state(istate)%rep(irep)%nframes_per_dcd_use(ifile)

        !--------------------------------------------------
        ! GLOBAL TIME (correct)
        !--------------------------------------------------
        map%segment(i)%start_time_step_global = seg_use_start
        map%segment(i)%end_time_step_global   = seg_use_end
        if( map%segment(i)%end_frame_in_file&
        == map%segment(i)%nframes_raw-project%state(istate)%rep(irep)%ignore_last_per_dcd(ifile)) then
        map%segment(i)%use_all=.true.
        else
        map%segment(i)%use_all=.false.	
        end if
        
    end do

end subroutine

subroutine compute_file_global_range(project, istate, irep, ifile, gstart, gend)

    implicit none

    type(project_type), intent(in) :: project
    integer, intent(in) :: istate, irep, ifile
    integer, intent(out) :: gstart, gend

    integer :: i
    integer :: t

    t = 1 - project%state(istate)%rep(irep)%offset_start

    do i = 1, ifile - 1
        t = t + project%state(istate)%rep(irep)%nframes_per_dcd_use(i)
    end do

    gstart = t
    gend   = t + project%state(istate)%rep(irep)%nframes_per_dcd_use(ifile) - 1

end subroutine


subroutine build_frame_to_keep_mask(iseg,project,window,map,window_mask,frame_to_keep)


    implicit none
    integer                             :: iseg
    type(window_map),      intent(in)   :: map
    type(project_type),    intent(in)   :: project
    type(time_window),  intent(in)      :: window
    logical, intent(in) :: window_mask(:)
    
    logical, allocatable, intent(out) :: frame_to_keep(:)

    integer :: nraw
    integer :: nuse
    integer :: i,f
    integer :: tglobal
    
    
    if(size(window_mask)==window%end_time_step-window%start_time_step+1) then
    else
    	 write(0,*) 'missmatach in number of time_steps in window ?'
    write(0,*) size(window_mask),window%end_time_step-window%start_time_step+1 
    stop
    end if
    

    nraw = map%segment(iseg)%nframes_raw
    nuse = map%segment(iseg)%nframes_raw
    if(map%segment(iseg)%use_all) then
       allocate(frame_to_keep(nraw))
    else
       i=map%segment(iseg)%end_time_step_global+map%segment(iseg)%start_frame_in_file&
       -map%segment(iseg)%start_time_step_global
        allocate(frame_to_keep(i))
    end if
    
    frame_to_keep=.false.
    

    frame_to_keep(1:map%segment(iseg)%start_frame_in_file-1)=.false.
    f=map%segment(iseg)%start_frame_in_file-1
    write(0,*) f
    write(0,*) 1+map%segment(iseg)%end_time_step_global-map%segment(iseg)%start_time_step_global
    do i = 1, 1+map%segment(iseg)%end_time_step_global-map%segment(iseg)%start_time_step_global
        f=f+1
        tglobal = &
            map%segment(iseg)%start_time_step_global+i-1
        !write(0,*) f,tglobal
        !pause
        frame_to_keep(f) = window_mask(i)
        
    end do
    !write(0,*)size(frame_to_keep) 
     
    !do i=1,size(frame_to_keep) 
    !write(0,*) frame_to_keep(i)
    !end do
    
  !    do iseg = 1, map%nsegment
  !    
  !        write(0,*) map%segment(iseg)%start_time_step_global-window%start_time_step-1,'!!'
  !        !map%segment(iseg)%nframes_keep=count(window_mask(iseg))
  !          pause
  !  end do
  

end subroutine




!=========================================================
end module time_window_file_mapping