module trajectory_module


contains

subroutine write_dcd(filename, natom, nframes,delta)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: natom, nframes
    real, intent(in) :: tra(3,natom,nframes)
   ! real, intent(in) :: pos_x(nframes, natom), pos_y(nframes, natom), pos_z(nframes, natom)
    integer*8, intent(in) :: delta

    character*4 :: coor
    integer :: i, n
    integer :: two, twentyfour, zero, peroff
    character*80 :: string1, string2
    integer :: io_open

    ! Header Variablen
    integer :: ifirst, nfreq, total

    ! Setze Header-Werte
    coor = 'CORD'
    ifirst = 0
    nfreq = 1
    total = nframes
    peroff = 0
    twentyfour = 24
    zero = 0
    two = 2
    string1 = 'Welcome to the heaven!'
    string2 = 'Go to hell!'

    ! --- Datei öffnen ---
    open(10, file=trim(filename), form='unformatted', status='replace', iostat=io_open)
    if (io_open /= 0) then
        print *, 'Error opening file for writing:', filename
        stop
    end if

    ! --- Header schreiben ---
    write(10) coor, nframes, ifirst, nfreq, total, zero, zero, zero, zero, zero, delta, peroff, zero, zero, zero, zero, zero, zero, zero, twentyfour
    write(10) two, string1, string2
    write(10) natom

    ! --- Frames schreiben ---
    write(0,*) shape(pos_x),'shape pos_x'	
    do n = 1, nframes
        write(10) (tra(1,n,i), i=1,natom)
        write(10) (tra(2,n,i), i=1,natom)
        write(10) (tra(3,n,i), i=1,natom)
    end do

    close(10)
    print *, 'DCD written successfully:', filename
    
end subroutine write_dcd

subroutine read_dcd_start_to_end_no_box(filename, start_frame, end_frame)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: start_frame, end_frame

    integer :: i, j, natom_header, nframes_header, nframes_to_read
    integer :: start_local, end_local, io_open
    real :: d(6), dummyr
    character*4 :: dummyc
    integer :: dummyi
    real, allocatable :: tmp_x(:), tmp_y(:), tmp_z(:)

    ! --- DCD Datei öffnen ---
    open(10, file=trim(filename), status='old', form='unformatted', iostat=io_open)
    if (io_open /= 0) then
        print *, "Error opening DCD file: ", trim(filename)
        stop
    end if

    ! --- Header lesen ---
    read(10) dummyc, nframes_header, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(10) dummyi, dummyr
    read(10) natom_header

    ! Validierung Start/End
    start_local = max(1, start_frame)
    end_local   = min(nframes_header, end_frame)
    if (start_local > end_local) then
        print *, "Error: start_frame > end_frame"
        stop
    end if
    nframes_to_read = end_local - start_local + 1
    print *, "Reading frames ", start_local, " to ", end_local, " (total ", nframes_to_read, ")"

    ! --- Globale Arrays allokieren ---
    if (allocated(pos_x)) deallocate(pos_x)
    if (allocated(pos_y)) deallocate(pos_y)
    if (allocated(pos_z)) deallocate(pos_z)
    if (allocated(box))   deallocate(box)

    allocate(pos_x(nframes_to_read, natom_header))
    allocate(pos_y(nframes_to_read, natom_header))
    allocate(pos_z(nframes_to_read, natom_header))
    allocate(box(nframes_to_read, 3))

    ! Temporäre Arrays für Einlesen eines Frames
    allocate(tmp_x(natom_header))
    allocate(tmp_y(natom_header))
    allocate(tmp_z(natom_header))

    ! --- Frames einlesen ---
    do i = 1, nframes_header
       ! read(10) (d(j), j=1,6)            ! Box info
        read(10) tmp_x
        read(10) tmp_y
        read(10) tmp_z

        if (i >= start_local .and. i <= end_local) then
            pos_x(i-start_local+1,:) = tmp_x
            pos_y(i-start_local+1,:) = tmp_y
            pos_z(i-start_local+1,:) = tmp_z
            box(i-start_local+1,1) = d(1)
            box(i-start_local+1,2) = d(3)
            box(i-start_local+1,3) = d(6)
        end if
    end do

    ! Aufräumen temporäre Arrays
    deallocate(tmp_x, tmp_y, tmp_z)

    close(10)
    print *, "Finished reading DCD: ",trim(filename), nframes_to_read, " frames, ", natom_header, " atoms."
end subroutine read_dcd_start_to_end_no_box



end module trajectory_module