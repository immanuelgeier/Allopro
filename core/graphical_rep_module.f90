module graphical_rep_module

use project_module
use topology_module
use sidechains
use domains_module
use time_window_file_mapping
use time_window_loader
use direction_vectors_analyse_plan_module
use position_lift_module
use io_pdb_module
use te_io_module
implicit none

type :: frame2plot
    integer :: istate,irep,itime
end type


contains

subroutine make_vmd_line_graph(project,topologies,plan,frame)

    type(topology_set_type), intent(in) :: topologies
    type(project_type), intent(in) :: project
    type(frame2plot) :: frame
    type(cv_plan), intent(in) :: plan

    type(raw_data_block) :: block
    type(time_window) :: window
    type(lifted_trajectory) :: lifted

    real(8), allocatable :: vecs(:,:), centroids(:,:)

    logical, allocatable :: mask(:)
    real(4), allocatable :: occ(:), beta(:),pos(:,:)
    logical, allocatable :: group(:)
    real(8), allocatable ::si(:,:)
    integer :: n_objecte
    integer :: i,mind(2)

    character(len=400) :: outfile
    character(len=400) :: vmdfile

    real(8) :: start(3), ende(3)
    real(8) :: scale,val,max_te,x
    character(len=128), allocatable:: names(:)
    real(8) :: rgb(3)
    integer :: j,k
    
    integer :: nobs
    integer :: obs_i, obs_j
    integer :: obj_i, obj_j
    integer :: nplot
    integer :: t

    integer, allocatable :: obs_to_obj(:)
    integer, allocatable :: obs_to_comp(:)
    logical, allocatable :: obs_is_sidechain(:)
    integer :: color_id
integer, parameter :: ncolors = 64


color_id = 100 + int(x * real(ncolors - 1,8))

    !------------------------------------------
    ! Set time window
    !------------------------------------------

    window%istate          = frame%istate
    window%irep            = frame%irep
    window%start_time_step = frame%itime
    window%end_time_step   = frame%itime

    allocate(mask(1))
    mask = .true.


    !------------------------------------------
    ! Load trajectory frame
    !------------------------------------------

    call load_raw_data_block( &
         project, &
         window, &
         mask, &
         block)

    call load_te_data()
    !------------------------------------------
    ! Get direction vectors
    !------------------------------------------

    call lift_cv_plan(plan, block, lifted)

    allocate(vecs(3,size(lifted%vec,2)))
    vecs = lifted%vec(:,:,1)


    !------------------------------------------
    ! Get centroid positions
    !------------------------------------------

    call lift_cv_plan_centroid(plan, block, lifted)

    allocate(centroids(3,size(lifted%vec,2)))
    centroids = lifted%vec(:,:,1)


    !------------------------------------------
    ! Original atom positions for PDB
    !------------------------------------------

    allocate(pos(3,size(block%xyz,2)))
    pos = block%xyz(:,:,1)

    allocate(occ(size(pos,2)))
    allocate(beta(size(pos,2)))
    allocate(group(size(pos,2)))

    group = .true.
    occ   = 0.0
    beta  = 0.0

    call lift_cv_plan_occ(plan,block,lifted, occ)

    !------------------------------------------
    ! Write PDB
    !------------------------------------------

    outfile = 'test_pdb_graph.pdb'

    call write_subgroup_pdb( &
         outfile, &
         topologies%state(frame%istate)%topo, &
         pos, &
         occ, &
         beta, &
         group)


    !------------------------------------------
    ! Create VMD script
    !------------------------------------------

    vmdfile = 'pvg.tcl'

    call make_vmd_plot_script(vmdfile)


    !------------------------------------------
    ! Draw vectors
    !------------------------------------------

    scale = 1.0d0
    !allocate(si(size(delta_te_matrix,1),size(delta_te_matrix,2)))
    !si=delta_te_matrix
    
    
    !------------------------------------------------------------
    ! Find and plot 100 strongest TE signals
    !------------------------------------------------------------

   ! allocate(si(size(delta_te_matrix,1),size(delta_te_matrix,2)))
   ! si = delta_te_matrix
   nobs = size(delta_te_matrix,1)

    allocate(obs_to_obj(nobs))
    allocate(obs_to_comp(nobs))
    allocate(obs_is_sidechain(nobs))

    obs_to_obj       = 0
    obs_to_comp      = 0
    obs_is_sidechain = .false.


    !------------------------------------------------------------
    ! Build mapping:
    !
    ! observable index -> group/object + component
    !
    ! target = group number
    ! comp   = 1 or 2
    !------------------------------------------------------------

    do t = 1, size(plan%term)

        if (plan%term(t)%target < 1 .or. &
            plan%term(t)%target > size(plan%object)) cycle


        obs_i = 2 * (plan%term(t)%target - 1) + &
                plan%term(t)%comp

        if (obs_i < 1 .or. obs_i > nobs) cycle

        obs_to_obj(obs_i)  = plan%term(t)%target
        obs_to_comp(obs_i) = plan%term(t)%comp

        if (plan%object(plan%term(t)%target)%object_type == &
            2) then

            obs_is_sidechain(obs_i) = .true.

        end if

    end do


allocate(si(nobs,nobs))
si = delta_te_matrix


nplot = 0

do while (nplot < 250)

    mind = maxloc(si)

    obs_i = mind(1)
    obs_j = mind(2)
    


    val = si(obs_i,obs_j)
   

    if (abs(val) == 0.0d0) exit

    ! Eintrag entfernen
    si(obs_i,obs_j) = 0.0d0
    si(obs_j,obs_i) = 0.0d0
   if(obs_i<=18 .or. obs_j<=18) cycle
   if( obs_is_sidechain(obs_i) .eqv. .false.) cycle
   if( obs_is_sidechain(obs_j) .eqv. .false.) cycle
    
       
     if(obs_i==obs_j) cycle
     if((obs_i+1)/2==(obs_j+1)/2) cycle
     
     if(nplot==0) max_te=val
     
     write(0,*) val,obs_i,obs_j,max_te
     !pause
    ! ---------------------------------------------------------
    ! WICHTIG:
    ! CV-Index == Centroid-Index
    ! ---------------------------------------------------------

    if (obs_i < 1 .or. obs_i > size(centroids,2)) cycle
    if (obs_j < 1 .or. obs_j > size(centroids,2)) cycle

    start = centroids(:,obs_i)
    ende  = centroids(:,obs_j)

    ! ---------------------------------------------------------
    ! Debug
    ! ---------------------------------------------------------

    write(0,'("TE line ",I4, &
              ": CV ",I5," -> ",I5, &
              " group ",I5," -> ",I5, &
              " comp ",I1," -> ",I1, &
              " dist2=",ES14.6, &
              " TE=",ES14.6)') &
         nplot+1, &
         obs_i, obs_j, &
         (obs_i+1)/2, (obs_j+1)/2, &
         mod(obs_i-1,2)+1, &
         mod(obs_j-1,2)+1, &
         sqrt(dot_product(start-ende,start-ende)), &
         val

    ! ---------------------------------------------------------
    ! Farbe
    ! ---------------------------------------------------------

    !------------------------------------------------------------
! Color: white (TE=0) -> deep red (TE=max_te)
!------------------------------------------------------------

if (max_te > 0.0d0) then

    x = val / max_te

    ! Sicherheitsbegrenzung
    x = max(0.0d0, min(1.0d0, x))

    rgb(1) = 1.0d0
    rgb(2) = 1.0d0 - x
    rgb(3) = 1.0d0 - x

else

    rgb = (/ 1.0d0, 1.0d0, 1.0d0 /)

end if

nplot = nplot + 1

call add_line_to_vmd_plot_script_w( &
     vmdfile, &
     start, &
     ende, &
x)
     
end do

write(0,'("Plotted ",I5," TE lines")') nplot

deallocate(si)
    deallocate(obs_to_obj)
    deallocate(obs_to_comp)
    deallocate(obs_is_sidechain)
    

    !------------------------------------------

    deallocate(mask)
    deallocate(pos)
    deallocate(vecs)
    deallocate(centroids)
    deallocate(occ)
    deallocate(beta)
    deallocate(group)

    call  get_cv_name_list(plan,names)
    open(1,file='cv_name_list.txt')
    do i=1,size(names)
     write(0,*) i,trim(names(i))
     write(1,*) i,trim(names(i))
    end do
    close(1)
    
end subroutine make_vmd_line_graph

!!------------------------------
!subroutine generate_pdb_file(pos,filename,occ_data,beta_data)
!implicit none
!type( frame2plot) :: frame
!character(len=*) :: filename
!
!
!
!end subroutine generate_pdb_file
!------------------------------
subroutine get_pos_vecs_and_centroids(project,topologies,plan,frame)

    type(topology_set_type), intent(in) :: topologies
    type(project_type), intent(in) :: project
    type(frame2plot) :: frame
    type(cv_plan), intent(in) :: plan

    type(raw_data_block) :: block
    type(time_window) :: window
    type(lifted_trajectory) :: lifted

    real(8), allocatable :: vecs(:,:), centroids(:,:)

    logical, allocatable :: mask(:)
    real(4), allocatable :: occ(:), beta(:),pos(:,:)
    logical, allocatable :: group(:)

    integer :: n_objecte
    integer :: i

    character(len=400) :: outfile
    character(len=400) :: vmdfile

    real(8) :: start(3), ende(3)
    real(8) :: scale,x
    character(len=128), allocatable:: names(:)

    !------------------------------------------
    ! Set time window
    !------------------------------------------

    window%istate          = frame%istate
    window%irep            = frame%irep
    window%start_time_step = frame%itime
    window%end_time_step   = frame%itime

    allocate(mask(1))
    mask = .true.


    !------------------------------------------
    ! Load trajectory frame
    !------------------------------------------

    call load_raw_data_block( &
         project, &
         window, &
         mask, &
         block)


    !------------------------------------------
    ! Get direction vectors
    !------------------------------------------

    call lift_cv_plan(plan, block, lifted)

    allocate(vecs(3,size(lifted%vec,2)))
    vecs = lifted%vec(:,:,1)


    !------------------------------------------
    ! Get centroid positions
    !------------------------------------------

    call lift_cv_plan_centroid(plan, block, lifted)

    allocate(centroids(3,size(lifted%vec,2)))
    centroids = lifted%vec(:,:,1)


    !------------------------------------------
    ! Original atom positions for PDB
    !------------------------------------------

    allocate(pos(3,size(block%xyz,2)))
    pos = block%xyz(:,:,1)

    allocate(occ(size(pos,2)))
    allocate(beta(size(pos,2)))
    allocate(group(size(pos,2)))

    group = .true.
    occ   = 0.0
    beta  = 0.0

    call lift_cv_plan_occ(plan,block,lifted, occ)

    !------------------------------------------
    ! Write PDB
    !------------------------------------------

    outfile = 'test_pdb.pdb'

    call write_subgroup_pdb( &
         outfile, &
         topologies%state(frame%istate)%topo, &
         pos, &
         occ, &
         beta, &
         group)


    !------------------------------------------
    ! Create VMD script
    !------------------------------------------

    vmdfile = 'plot_vectors.tcl'

    call make_vmd_plot_script(vmdfile)


    !------------------------------------------
    ! Draw vectors
    !------------------------------------------

    scale = 5.0d0

    do i = 1, size(centroids,2)
        write(0,*) centroids(:,i)
        start = centroids(:,i)

        ende = start + scale * vecs(:,i)

        call add_line_to_vmd_plot_script( &
             vmdfile, &
             start, &
             ende, &
             1, &       ! color
             3)         ! line width

    end do


    !------------------------------------------
    ! Cleanup
    !------------------------------------------

    deallocate(mask)
    deallocate(pos)
    deallocate(vecs)
    deallocate(centroids)
    deallocate(occ)
    deallocate(beta)
    deallocate(group)

    call  get_cv_name_list(plan,names)
    open(1,file='cv_name_list.txt')
    do i=1,size(names)
     write(0,*) i,trim(names(i))
     write(1,*) i,trim(names(i))
    end do
    close(1)
    
end subroutine get_pos_vecs_and_centroids




subroutine add_line_to_vmd_plot_script(filename,start,ende,color,linewi)
    implicit none

    character(len=*), intent(in) :: filename
    real(8), intent(in) :: start(3), ende(3)
    integer, intent(in) :: color, linewi

    integer :: unit

    open(newunit=unit, file=filename, status='old', &
         action='write', position='append', form='formatted')

    ! Set color and line width
    write(unit,'("draw color ",I0)') color
   !write(unit,'("draw width ",I0)') linewi

    ! Draw line from start to end
    write(unit,'("draw line {",3(ES20.12,1X),"} {",3(ES20.12,1X),"}")') &
         start(1),start(2),start(3),ende(1),ende(2),ende(3)

    close(unit)

end subroutine add_line_to_vmd_plot_script
subroutine make_vmd_plot_script(filename)
    implicit none

    character(len=*), intent(in) :: filename
    integer :: unit
    integer :: i
    real(8) :: x
    real(8) :: r, g, b

    open(newunit=unit, file=filename, status='replace', &
         action='write', form='formatted')

    write(unit,'(A)') '# VMD visualization script generated by Fortran'
    write(unit,'(A)') 'draw delete all'
    write(unit,'(A)') ''

    ! ------------------------------------------------------------
    ! Define the color palette ONCE (IDs 100 to 599)
    ! 100 = white, 599 = red
    ! ------------------------------------------------------------
    do i = 0, 499
        x = real(i,8) / 499.0d0

        r = 1.0d0
        g = 1.0d0 - x
        b = 1.0d0 - x

        write(unit, &
            '("color change rgb ",I0," ",F10.6," ",F10.6," ",F10.6)') &
            100+i, r, g, b
    end do

    write(unit,'(A)') ''
    write(unit,'(A)') '# --- Zylinder folgen ---'
    write(unit,'(A)') ''

    close(unit)

end subroutine make_vmd_plot_script

subroutine add_line_to_vmd_plot_script_w(filename, start, ende, linewi)
    implicit none

    character(len=*), intent(in) :: filename
    real(8), intent(in) :: start(3), ende(3)
    !real(8), intent(in) :: val      ! Hier vorher "rgb(3)", jetzt ein klarer Skalar (z.B. Wert von 0 bis 1)
    real(8), intent(in) :: linewi

    integer :: unit
    integer :: color_id
    real(8) :: radius
    real(8) :: intensity

    radius = 0.5d0 * linewi

    ! intensity soll von 0 bis 1 laufen
    intensity =1.! max(0.0d0, min(1.0d0, val))

    ! Mappe intensity [0, 1] sauber auf die IDs 100 bis 599
    color_id = 100 + int(intensity * 499.0d0 + 0.5d0)  ! +0.5d0 für sauberes Runden

    if (color_id < 100) color_id = 100
    if (color_id > 599) color_id = 599

    open(newunit=unit, file=filename, status='old', &
         action='write', position='append', form='formatted')

    write(unit, '("graphics top color ",I0)') color_id

    write(unit, '("graphics top cylinder {",3(ES20.12,1X),"} {", &
                3(ES20.12,1X),"} radius ",ES12.5)') &
          start(1), start(2), start(3), &
          ende(1), ende(2), ende(3), &
          radius

    close(unit)

end subroutine add_line_to_vmd_plot_script_w

subroutine add_line_to_vmd_plot_script_rgb(filename, start, ende, val, linewi)
    implicit none

    character(len=*), intent(in) :: filename
    real(8), intent(in) :: start(3), ende(3)
    real(8), intent(in) :: val      ! Hier vorher "rgb(3)", jetzt ein klarer Skalar (z.B. Wert von 0 bis 1)
    integer, intent(in) :: linewi

    integer :: unit
    integer :: color_id
    real(8) :: radius
    real(8) :: intensity

    radius = 0.05d0 * real(linewi, 8)

    ! intensity soll von 0 bis 1 laufen
    intensity = max(0.0d0, min(1.0d0, val))

    ! Mappe intensity [0, 1] sauber auf die IDs 100 bis 599
    color_id = 100 + int(intensity * 499.0d0 + 0.5d0)  ! +0.5d0 für sauberes Runden

    if (color_id < 100) color_id = 100
    if (color_id > 599) color_id = 599

    open(newunit=unit, file=filename, status='old', &
         action='write', position='append', form='formatted')

    write(unit, '("graphics top color ",I0)') color_id

    write(unit, '("graphics top cylinder {",3(ES20.12,1X),"} {", &
                3(ES20.12,1X),"} radius ",ES12.5)') &
          start(1), start(2), start(3), &
          ende(1), ende(2), ende(3), &
          radius

    close(unit)

end subroutine add_line_to_vmd_plot_script_rgb

end module graphical_rep_module