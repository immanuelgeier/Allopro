module time_window_loader

use project_module
use time_window_file_mapping
use adaptive_dcd_reader

implicit none

private

public :: raw_data_block
public :: load_raw_data_block

!=========================================================
type :: raw_data_block

    integer :: natom
    integer :: nframe

    integer :: istate
    integer :: irep

    integer :: start_time_step
    integer :: end_time_step

    real, allocatable :: xyz(:,:,:)

end type

!=========================================================
contains
!

!=========================================================
subroutine load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)

    type(project_type), intent(in) :: project
    type(time_window),  intent(in) :: window

    logical, intent(in) :: window_mask(:)

    type(raw_data_block), intent(out) :: block

    type(window_map) :: map

    integer :: natom
    integer :: nframe_total

    integer :: iseg
    integer :: iframe_out
    integer :: n_keep,nt_start,nt_end,n_keep_sum

    !-----------------------------------------------------
    ! mapping
    !-----------------------------------------------------
    !write(0,*) window%istate,'test3'
    call get_window_map(window, project, map)

    natom = project%state(window%istate) &
                   %rep(window%irep) &
                   %nfirst_natom_to_read

    !-----------------------------------------------------
    ! total output frames
    !-----------------------------------------------------
    nframe_total = 0
    n_keep=0

    do iseg = 1, map%nsegment
        nframe_total = nframe_total + &
                       map%segment(iseg)%nframes_use
       nt_start=map%segment(iseg)%start_time_step_global-window%start_time_step+1
       nt_end  =map%segment(iseg)%end_time_step_global-window%start_time_step+1
       n_keep=n_keep+ count(window_mask(nt_start:nt_end))          
    end do
    
    
    !write(0,*) n_keep,'mmmmmmmmmmmmmmmmmm'
    !pause

    allocate(block%xyz(3,natom,n_keep))
    block%xyz=-42._8
    

    block%natom = natom
    block%nframe = nframe_total

    block%istate = window%istate
    block%irep   = window%irep

    block%start_time_step = window%start_time_step
    block%end_time_step   = window%end_time_step

    iframe_out = 1

    !-----------------------------------------------------
    ! load segments
    !-----------------------------------------------------
    n_keep_sum=0
    n_keep=0
    do iseg = 1, map%nsegment
        !write(0,*) map%nsegment,'nesegments'
       ! pause
       nt_start=map%segment(iseg)%start_time_step_global-window%start_time_step+1
       nt_end  =map%segment(iseg)%end_time_step_global-window%start_time_step+1
       
	   n_keep_sum=n_keep_sum+n_keep
       n_keep=count(window_mask( nt_start:nt_end))
       
       
       write(0,*) n_keep_sum+1,n_keep_sum+n_keep,'range',iseg
        call load_segment( &
            project, &
            window, &
            map, &
            iseg, &
            window_mask, &
            block%xyz(:,:, &
            n_keep_sum+1:n_keep_sum+n_keep ))
            

            
          
        iframe_out = iframe_out + &
                     map%segment(iseg)%nframes_use

    end do

end subroutine


!=========================================================
subroutine load_segment( &
    project, &
    window ,&
    map, &
    iseg, &
    window_mask, &
    xyz)

    type(project_type), intent(in) :: project
    type(window_map),   intent(in) :: map
    type(time_window) :: window

    integer, intent(in) :: iseg

    logical, intent(in) :: window_mask(:)

    real, intent(inout) :: xyz(:,:,:)

    logical, allocatable :: frame_to_keep(:)

    type(adaptive_dcd_type) :: reader

    character(400) :: filename

    integer :: nread

    !-----------------------------------------------------
    ! segment mask
    !-----------------------------------------------------
    call build_frame_to_keep_mask( &
        iseg, &
        project, &
        window, &
        map, &
        window_mask, &
        frame_to_keep)

    ! write(0,*)'frame_to_keep', frame_to_keep  
     !pause
    !-----------------------------------------------------
    ! open dcd
    !-----------------------------------------------------
    filename = &
        project%state(map%segment(iseg)%istate) &
               %rep(map%segment(iseg)%irep) &
               %dcd_files(map%segment(iseg)%ifile)

    call adaptive_open(reader, filename)

    !write(0,*) nread
    !-----------------------------------------------------
    ! read selected frames
    !-----------------------------------------------------
    call read_selected_frames( &
        reader, &
        frame_to_keep, &
        xyz, &
        nread)

    !-----------------------------------------------------
    ! sanity check
    !-----------------------------------------------------
    if (count( frame_to_keep) /= nread) then

        write(0,*) 'ERROR load_segment'
        write(0,*) 'segment =', iseg
        write(0,*) 'expected=', count( frame_to_keep) 
        write(0,*) 'read    =', nread
        pause
        !stop

    end if

    call adaptive_close(reader)

    deallocate(frame_to_keep)

end subroutine


!=========================================================
subroutine read_selected_frames( &
    reader, &
    frame_to_keep, &
    xyz, &
    nread)

    type(adaptive_dcd_type), intent(inout) :: reader

    logical, intent(in) :: frame_to_keep(:)

    real, intent(inout) :: xyz(:,:,:)

    integer, intent(out) :: nread

    integer :: iframe
    integer :: iframe_out

    integer :: natom_read

    real :: box

    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: z(:)

    natom_read = size(xyz,2)

    allocate(x(reader%natom))
    allocate(y(reader%natom))
    allocate(z(reader%natom))

    iframe_out = 0

    do iframe = 1, size(frame_to_keep)

        read(reader%unit) box
        read(reader%unit) x
        read(reader%unit) y
        read(reader%unit) z
        
        if (frame_to_keep(iframe)) then
        
            iframe_out = iframe_out + 1

            xyz(1,:,iframe_out) = x(1:natom_read)
            xyz(2,:,iframe_out) = y(1:natom_read)
            xyz(3,:,iframe_out) = z(1:natom_read)
            !write(0,*) iframe,'frame_to_keep(iframe)',frame_to_keep(iframe),iframe_out
            !write(0,*) x(1:5)
        end if
        

    end do

    nread = iframe_out

    deallocate(x)
    deallocate(y)
    deallocate(z)

end subroutine

!=========================================================
end module time_window_loader