program main

use system_module
use project_module
use topology_module
use analysis_plan_module
use time_window_file_mapping
use test_time_window_file_mapping
use time_window_loader
use sidechains
use sidechainhead_refstruc
use domains_module
use domain_refstruc
use user_defined_cv_module
use direction_vectors_analyse_plan_module	
use sparse_lut_module
use lifting_plan_module
use lifting_plan_execution_module
use healpix_state_space_exploration
use te_plot_module
use graphical_rep_module
use refstruc_pbc_module
use healpix_histogram_sampling_execution_module
use centroid_distance_histogram_module
use activ_pair_selection_module
use ref_struc_selection_module
use random_seed_io_module, only: random_seed_config_type, &
                                 initialize_random_seed,configure_seed
use workflow_module  
use cv_lib_workflow_io_module
use cv_plan_workflow_io_module
use random_seed_workflow_module
use trajectory_settings_module
use euclidean_lifting_workflow_io_module
use trajectory_workflow_io_module
use test_trajectory_workflow_module
use graphical_representation_workflow_io_module
use replicate_selector_workflow_module

implicit none

type(project_type)           :: project
type(topology_set_type)      :: topologies
type(analysis_plan_type)     :: plan
type(sidechain_library_type) :: sc_lib
type(domain_library_type)    :: dom_lib
type(user_cv_library_type)   :: user_lib
type(cv_plan)                :: cv_plan1
type(lifting_plan_type)      :: lifting_plan
type (frame2plot)            :: frame2p
type (refstruc_type)         :: ref1
type(odf_lib_type)           :: odf_lib
type(time_window)            :: window1,window2
type(random_seed_config_type) :: seed_config
type(library_workflow_type)   :: cv_lib_workflow_initial
type(library_workflow_type)   :: cv_lib_workflow_use
type(cv_plan_workflow_type)   :: cv_plan_workflow1
type(euclidean_lifting_workflow_type):: euclidian_lifting_workflow
type(trajectory_workflow_type) :: traj_workflow
type(test_trajectory_workflow_type) :: write_test_tra
type(graphical_representation_workflow_type) :: plot_vectors_vmd
type(replicate_selector_workflow_type) :: replicate_selector_workflow
!integer :: irep,istate

logical :: found
integer :: n, clock,i
logical, allocatable :: matrix_active_final(:,:)

type(random_seed_workflow_type) ::seed_job, initial_seed_workflow,prodcution_seed_workflow
logical :: run


!call configure_euclidean_lifting_workflow(euclidian_lifting_workflow)
!stop
call run_workflow_driver()
!call configure_cv_plan_workflow(cv_plan_workflow1)

stop
!call configure_sub_workflows()
	


run=.false.
run=.true.

write(0,*)
write(0,*) 'INITIAL SEED WORKFLOW CONFIGURATION'

call configure_seed_workflow(initial_seed_workflow)
if(run) then
call run_seed_workflow(initial_seed_workflow)
end if

write(0,*)
write(0,*) 'PRODUCTION SEED WORKFLOW CONFIGURATION'

call configure_seed_workflow(prodcution_seed_workflow)
if(run) then
call run_seed_workflow(prodcution_seed_workflow)
end if


call initialize_system(project,topologies)

call configure_trajectory_workflow(traj_workflow,project,topologies)

!call initialize_trajectory_settings(project, topologies)

stop
if(1==1) then
call analyze_project_frames(project)
call print_project_summary(project)
end if



call init_resnames_and_sidechain_head_groups(sc_lib)
call init_sidechain_head(sc_lib)


call  create_multiple_user_cv_libraries(topologies)

write(0,*)!--------------------------------------------
write(0,*) 'INITIAL CV-LIBRARY WORKFLOW CONFIGURATION '
write(0,*)!--------------------
call configure_library_workflow(cv_lib_workflow_initial,topologies)




if(run) then
call run_library_workflow( &
    cv_lib_workflow_initial,topologies, sc_lib, dom_lib, user_lib)

end if

write(0,*)!---------------------------------------------
write(0,*) 'PRODUCTION CV-PLAN WORKFLOW CONFIGURATION '
write(0,*)!-------------------- ------------------------ 

call configure_cv_plan_workflow(cv_plan_workflow1,topologies)



stop


write(0,*) 'EUCLIDIAN LIFTING WORKFLWO CONFUGURATION:'
call configure_euclidean_lifting_workflow(euclidian_lifting_workflow)



if(run) then
call run_euclidean_lifting_workflow(euclidian_lifting_workflow, project, topologies)
end if


stop
	





!call te_rep_workflow(cv_plan1,project,topologies)
frame2p%istate=1;frame2p%irep=1;frame2p%itime=1
call make_vmd_line_graph(project,topologies,cv_plan1,frame2p)

!call get_pos_vecs_and_centroids(project,topologies,cv_plan1,frame2p)
stop





allocate( matrix_active_final(162,162))
matrix_active_final=.false.
call get_reduced_matrix(cv_plan1,162,matrix_active_final)
write(0,*) count(matrix_active_final),'in main'
stop
!open(1,file='activ_pairs')
!do i=1,162
!write(1,*) matrix_active_final(i,:)
!end do
!close(1)

open(1,file='activ_pairs')
do i=1,162
read(1,*) matrix_active_final(i,:)
write(0,*) matrix_active_final(i,:)
end do
close(1)
!
!stop
!call get_cv_residue_order(cv_plan1)
!stop



!call te_rep_workflow(cv_plan1,project,topologies)

!stop

!frame2p%istate=1;frame2p%irep=1;frame2p%itime=4510551
call prepare_pbc_correction(1,1,1,project,ref1)
    
!frame2p%istate=1;frame2p%irep=1;frame2p%itime=4510551
!
!call test_pbc_correction(project,topologies,frame2p,ref1)
!
!
!
!call get_pos_vecs_and_centroids(project,topologies, cv_plan1,frame2p)!,pos,vecs,centroids)

!stop

call init_lifting_plan_all(lifting_plan,0,10,1000000,project)



write(0,*) 'execute euclidian lifting'

if(1==2) then
call execute_lifting_plan_centroid(&
 cv_p               =cv_plan1&
,lifting_plan           =lifting_plan&
,project                =project&
,ref                    =ref1&
,path_lifted_tra        ='/media/pbuser/T7/lifted_trajectories/'&
,prefix_centroid_tra    ='cvec_tra'&
,path_centroid_tra      ='/media/pbuser/T7/lifted_trajectories/'&
,prefix_lifted_tra      ='vec_tra'&
,output_policy          ='ask')

end if

write(0,*) 



if(1==2) then
call run_centroid_distance_histogram(lifting_plan, '/media/pbuser/T7/lifted_trajectories/cvec_tra')
end if


!call lift_trajectories(lifting_plan, &
!path_vec_tra     ='/media/pbuser/T7/lifted_trajectories',&
!prefix_vec_tra   ='vec_tra',&
!path_tra         ='/media/pbuser/T7/lifted_trajectories',&
!prefix_tra       ='htra_b9_l20_m100',&
!npix, lag,basis_index,tau_mem,dt,file_policy_tra)



  
 if(1==1) then
call run_te_workflow( &
    odf_lib            = odf_lib, &
    lifting_plan       = lifting_plan, &
    lag                = 10, &
    basis_index        = 9, &
    prefix_vec_tra     = 'vec_tra', &
    path_vec_tra       = '/media/pbuser/T7/lifted_trajectories/', &
    prefix_tra         = 'htra_b9_l100_m100_e4', &
    path_tra           = '/media/pbuser/T7/lifted_trajectories/', &
    file_policy_tra    = 'replace', &
    prefix_hist        = 'hte_hist_b9_l200_m100_e4', &
    path_hist          = '/media/pbuser/T7/histograms/', &
    file_policy_hist   = 'replace', &
    nlevel_target      = 4, &
    max_chunk_size     = 2*10**8, &
    tau_mem            = 100._8, &
    dt                 = 1._8, &
    max_ram_gb         = 8.0_8, &
    block_shuffle      = .false., &
    nrun_shuffle       = 1, &
    block_size_shuffle = 10000, &
    only_lifting=.true. ,&
    prefix_odf ='odf_',&
    path_odf   ='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/',&
    pair_active= matrix_active_final)
 

 end if
 
 

call init_odf_lib(cv_plan1,odf_lib, lifting_plan,&
max_lag=20,&
basis_index=9,&
nlevel=6,&
nobj=162)



if(1==2) then
call sample_odfs_C2( &
    odf_lib, &
    lifting_plan, &
    prefix_tra='vec_tra', &
    path_tra='/media/pbuser/T7/lifted_trajectories/', &
    basis_index=9, &
    tau_mem            = 100._8, &
    dt                 = 1._8,&
    prefix_odf ='odf_',&
    path_odf   ='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/',&
    output_policy= 'replace')
    end if


 if(1==2) then
call run_te_workflow( &
    odf_lib            = odf_lib, &
    lifting_plan       = lifting_plan, &
    lag                = 10, &
    basis_index        = 9, &
    prefix_vec_tra     = 'vec_tra', &
    path_vec_tra       = '/media/pbuser/T7/lifted_trajectories/', &
    prefix_tra         = 'htra_b9_l100_m100_e4', &
    path_tra           = '/media/pbuser/T7/lifted_trajectories', &
    file_policy_tra    = 'ask', &
    prefix_hist        = 'hte_hist_b9_l100_m100_e4', &
    path_hist          = '/media/pbuser/T7/histograms/', &
    file_policy_hist   = 'ask', &
    nlevel_target      = 4, &
    max_chunk_size     = 2*10**8, &
    tau_mem            = 100._8, &
    dt                 = 1._8, &
    max_ram_gb         = 8.0_8, &
    block_shuffle      = .false., &
    nrun_shuffle       = 1, &
    block_size_shuffle = 10000, &
    only_lifting=.false.,&
    prefix_odf ='odf_',&
    path_odf   ='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/',&
    pair_active=matrix_active_final)


 end if
 

 !Random kontroll set muss gespeichert werden
call run_calculate_te_matrix_chunked_C2(& 
    cv_p           = cv_plan1,&
    odf_lib        = odf_lib, &
    lifting_plan   = lifting_plan, &
    input_path     = '/media/pbuser/T7/histograms/', &
    input_prefix   = 'hte_hist_b9_l100_m100_e4', &
    output_path    = '/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/', &
    output_prefix  = 'C2_te_l100_e4', &
    nlevel_target  = 4, &
    max_chunk_size = 2*10**8,&
    output_policy  = 'ask',& 
    prefix_odf ='odf_',&
    path_odf   ='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/',&
    shuffled=.false.,&
     pair_active=matrix_active_final)
 


stop

call run_te_pair(odf_lib, lifting_plan, 77, 78,found,    prefix_odf ='odf_',&
    path_odf   ='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/')
 


stop	
call test_mapping_real(odf_lib,    prefix_odf ='odf_',&
     path_odf   ='/home/pbuser/Schreibtisch/workspace_mdgraph/project/output/3lnx/')

stop
write(0,*) 'start sampling'
call sample_rare_events(odf_lib,lifting_plan,10)

!
	
stop
!! 1. Selector für den State initialisieren (z.B. State 1 hat 10 Replikate)
!call init_replica_selector(rep_selector, nrep_total=project%state(1)%nrep)
!
!! 2. Selektives Ausschalten für Cross-Validation (z.B. die letzten beiden weglassen)
!rep_selector%is_active(3:10) = .false.
!
!! 3. Plan mit dem konfigurierten Selector initialisieren 
!call init_sampling_plan(train_plan, project, istate=1, L=1000000, &
!                        mode=SHUFFLE_INTER_REPLICA, selector=rep_selector)
                        
!call execute_healpix_analysis( &
!        project      = project, &
!        cv_p         = cv_plan1, &
!        istate       = 1, &
!        L            = 100000, &
!        mode         = SHUFFLE_INTER_REPLICA, &
!        selector     = rep_selector, &
!        tau          = 100)
!stop
!call test_lift_cv_plan(cv_plan1,project)
!test_healpix_density(tstart,tend,big_interval,tau,plan,project)
!do irep=1,10
!call test_healpix_density(1,1000000,1000,100,cv_plan1,project,1,irep,1)
!end do
!call test_lifting(topologies,sc_lib,dom_lib,project,6,7)
write(0,*) 'all done'
contains


subroutine run_graphical_representation_workflow(job,project,topologies)
    use workflow_module, only: graphical_representation_workflow_type
    use cv_plan_io_module, only: read_cv_plan
    implicit none
   
    type(graphical_representation_workflow_type), intent(in) :: job
    type(frame2plot) :: frame2pl
    
    
    type(project_type)      ,intent(in):: project
    type(topology_set_type) ,intent(in):: topologies
    type(cv_plan)          :: plan
    real(8)                :: scale
    character(len=2048) :: outfile
    character(len=2048) :: vmdfile

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '  Running Graphical Representation'
    write(*,*) '======================================'

    ! 1. Validierung vor dem Ausführen sicherstellen
    call validate_graphical_representation_workflow(job)
    
    call read_cv_plan(plan,job%params%cv_plan_file)
    
    ! 2. Zeitpunkt / Zustand für den Frame setzen (gemäß deiner Vorgabe)
    frame2pl%istate = job%params%state_index   ! 
    frame2pl%irep   = job%params%replica_index ! Replikat-Index
    frame2pl%itime  = job%params%time_index    ! Zeit-Index


    scale      = job%params%scale            ! Skalierungsfaktor (Default = 1.0)
    vmdfile    = job%params%tcl_script       ! TCL-Skript für VMD
    outfile    = job%params%pdb_file 


    write(*,*) '  -> State Index  : ', frame2pl%istate
    write(*,*) '  -> Replica Index: ', frame2pl%irep
    write(*,*) '  -> Time Index   : ', frame2pl%itime
    write(*,*) '  -> Scale        : ', scale
    write(*,*)

    !
     call get_pos_vecs_and_centroids(project, topologies, plan, frame2pl, scale, outfile, vmdfile)
     
    call destroy_cv_plan(plan) 
     
    write(*,*) 'Graphical representation workflow completed successfully.'

end subroutine run_graphical_representation_workflow


subroutine run_test_trajectory_workflow(job, project)
    use test_trajectory_workflow_module
    use cv_plan_io_module
    use position_lift_module

    implicit none

    type(test_trajectory_workflow_type), intent(in) :: job
   
    type(time_window) :: window
    type(raw_data_block) :: block
    type(cv_plan) :: plan_to_test
    type(project_type) :: project
    type(lifted_trajectory) :: lifted     
    logical, allocatable :: window_mask(:)
    
    integer :: ntstep, i, ios
    
    ! Platzhalter für die gelifteten Daten und das Centroid (entsprechend deiner Routinen anpassen)
    ! type(lifted_data_type) :: lifted
    ! type(centroid_type) :: centroid

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    Running Test Trajectory Workflow'
    write(*,*) '======================================'

    ! ---------------------------------------------------------
    ! 1. Zeitfenster für den Datenblock definieren & allokieren
    ! ---------------------------------------------------------
    window%istate          = job%istate
    window%irep            = job%irep
    window%start_time_step = job%start_time_step
    window%end_time_step   = job%end_time_step

    ntstep = window%end_time_step - window%start_time_step + 1
    if (ntstep <= 0) then
        error stop 'run_test_trajectory_workflow: Invalid time window range (end_time_step < start_time_step)'
    end if

    allocate(window_mask(ntstep))
    window_mask = .true.

    !  Stride auf die Maske anwenden, falls Schritte übersprungen werden sollen
    if (job%stride > 1) then
        do i = 1, ntstep
            if (mod(i - 1, job%stride) /= 0) then
                window_mask(i) = .false.
            end if
        end do
    end if

    write(*,'(A,I0,A,I0,A,I0,A,I0)') '-> Loading raw data block: istate=', window%istate, &
               ', irep=', window%irep, ', steps [', window%start_time_step, ' : ', window%end_time_step, ']'

    call load_raw_data_block(project, window, window_mask, block)

    ! ---------------------------------------------------------
    ! 2. CV Plan aus Datei einlesen
    ! ---------------------------------------------------------
    write(*,'(3A)') '-> Reading CV plan from "', trim(job%cv_plan_filename), '"...'
    call read_cv_plan(plan_to_test, trim(job%cv_plan_filename))
    write(*,*) '-> CV plan successfully loaded.'

    ! ---------------------------------------------------------
    ! 3. CV Plan "liften" (lift_cv_plan_all)
    ! ---------------------------------------------------------
    write(*,*) '-> Executing lift_cv_plan_all...'
    
    ! Aufruf deiner Lifting-Routine (Passe die Argumente 'lifted' und 'centroid' an deine Typen an)
     call  lift_euclidian_cv_plan(plan_to_test, block, lifted)

    ! ---------------------------------------------------------
    ! 4. Trajektorie / Ergebnisse speichern
    ! ---------------------------------------------------------
    write(*,'(3A)') '-> Saving resulting trajectory to "', trim(job%test_trajectory_filename), '"...'
    
    ! Beispiel für das Schreiben der Ergebnisse in die Zieldatei:
     call write_lifted_trajectory(trim(job%test_trajectory_filename), lifted)

    write(*,*) 'Test trajectory workflow execution completed successfully.'

    ! Aufräumen
    if (allocated(window_mask)) deallocate(window_mask)
    ! Falls plan%term / plan%object beim Einlesen alloziert wurden und hier freigegeben werden müssen:
    ! if (allocated(plan%term)) deallocate(plan%term)
    ! if (allocated(plan%object)) deallocate(plan%object)
!call destroy_plan(plan_to_test)

end subroutine run_test_trajectory_workflow

subroutine write_lifted_trajectory(filename, lifted)
    use iso_fortran_env, only: error_unit
    implicit none

    character(len=*), intent(in) :: filename
    type(lifted_trajectory), intent(in) :: lifted

    integer :: unit, ios, t, obj, idx
    real(8), allocatable :: row_data(:)
    character(len=32768) :: header_line
    character(len=64)    :: temp_str

    ! 1. Datei öffnen
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'Cannot open trajectory output file for writing: ', trim(filename)
        error stop
    end if

    ! 2. Header schreiben
    write(unit, '(A)') '# Lifted trajectory data'
    write(unit, '(A,I0)') '# Number of objects: ', lifted%nobject
    write(unit, '(A,I0)') '# Number of time steps: ', lifted%ntimes

    header_line = '# Col  1: TimeStep'
    do obj = 1, lifted%nobject
        write(temp_str, '(I0)') obj
        header_line = trim(header_line) // ' | Obj_' // trim(temp_str) // &
                      ' [b1_x b1_y b1_z b2_x b2_y b2_z Cent_x Cent_y Cent_z]'
    end do
    write(unit, '(A)') trim(header_line)

    ! 3. Temporäres Array für eine Zeile allokieren (9 Werte pro Objekt)
    allocate(row_data(9 * lifted%nobject))

! 4. Daten zeilenweise über das temporäre Array schreiben
    do t = 1, lifted%ntimes
        do obj = 1, lifted%nobject
            ! Startindex für dieses Objekt im Zeilen-Array berechnen
            idx = (obj - 1) * 9 + 1
            
            ! b1 (3 Komponenten: x, y, z)
            row_data(idx   : idx+2) = lifted%vec(:, 2*obj - 1, t)
            !write(0,*) lifted%vec(:, 2*obj - 1, t)
            !pause
            ! b2 (3 Komponenten: x, y, z)
            row_data(idx+3 : idx+5) = lifted%vec(:, 2*obj,     t)
            
            ! Centroid (3 Komponenten: x, y, z) -> Korrigiert auf 'obj' statt '2*obj - 1'
            row_data(idx+6 : idx+8) = lifted%centroid(:, obj, t)
        end do

        ! Zeile sauber ausgeben: TimeStep gefolgt vom gesamten Zeilen-Array
        write(unit, *) t, row_data
    end do

    ! 5. Aufräumen und schließen
    deallocate(row_data)
    close(unit, iostat=ios)
end subroutine write_lifted_trajectory




subroutine run_cv_plan_workflow(job, plan)
    use sidechains 
    use domains_module
    use user_defined_cv_module
    use cv_plan_io_module
    implicit none

    type(cv_plan_workflow_type), intent(in)    :: job
    type(cv_plan),               intent(out)   :: plan

    character(len=1024) :: filename
    character(len=32)   :: names(2)
    integer :: i, k, n, rid, ios, sc_mode, idx
    logical :: include_it, type_include
    integer, allocatable :: atoms(:)

    ! Temporäre Bibliotheks-Objekte zum Auslesen der Daten
    type(sidechain_library_type) :: sc_lib
    type(domain_library_type)    :: dom_lib
    type(user_cv_library_type)   :: user_lib

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '    Running CV Plan Workflow'
    write(*,*) '======================================'

    call print_cv_plan_workflow(job)
    !pause
    
    if (job%read_existing_cv_plan) then
        write(*,'(3A)') '-> Loading existing CV plan from "', trim(job%existing_cv_plan_filename), '"...'
        call read_cv_plan(plan, trim(job%existing_cv_plan_filename))
        write(*,*) '-> Existing CV plan successfully loaded.'
        return
    end if

    ! ---------------------------------------------------------
    ! 1. SIDECHAIN CVs verarbeiten (Unterstützt Typ- & Instanz-Filter)
    ! ---------------------------------------------------------
    if (job%include_sidechain_cvs .and. allocated(job%sidechain_libs)) then
        write(*,*) 'Processing sidechain libraries...'
        
        do i = 1, size(job%sidechain_libs)
            filename = job%sidechain_libs(i)%filename
            sc_mode = job%sidechain_libs(i)%selection_mode
            
            call read_sidechain_library(trim(filename), sc_lib)
            
            if (allocated(sc_lib%residue)) then
                idx = 0 ! Fortlaufender Index für Modus 4 (Instanz-Filter)
                
                do k = 1, sc_lib%ncv
                    ! Standard-Typ-Filter für Modus 1, 2, 3
                    type_include = .true.
                    if (sc_mode /= 4 .and. allocated(job%sidechain_libs(i)%filters)) then
                        if (k <= size(job%sidechain_libs(i)%filters)) then
                            type_include = job%sidechain_libs(i)%filters(k)
                        end if
                    end if

                    if (allocated(sc_lib%residue(k)%res_ids)) then
                        do n = 1, size(sc_lib%residue(k)%res_ids)
                            rid = sc_lib%residue(k)%res_ids(n)
                            
                            ! Wenn Modus 4 aktiv ist, prüfen wir den individuellen Instanz-Filter
                            if (sc_mode == 4) then
                                idx = idx + 1
                                include_it = .true.
                                if (allocated(job%sidechain_libs(i)%filters)) then
                                    if (idx <= size(job%sidechain_libs(i)%filters)) then
                                        include_it = job%sidechain_libs(i)%filters(idx)
                                    end if
                                end if
                            else
                                ! Für Modus 1-3 erbt die Instanz den Status des Residuen-Typs
                                include_it = type_include
                            end if

                            if (include_it) then
                                call make_cv_names(names, sc_lib%residue(k)%resname, rid)
                                call get_atoms_sidechain(job%istate, sc_lib, topologies, k, n, atoms)
                                call add_cv_pair_sidechain( &
                                    plan         = plan, &
                                    object_type  = 2, &  
                                    resid        = rid, &
                                    refid        = 0, &  
                                    name         = names, &
                                    point_group  = 'none', &  
                                    atoms        = atoms, &
                                    cv           = sc_lib%residue(k)%cv, &    
                                    has_c2_sym   = sc_lib%residue(k)%has_sym &
                                )
                            end if
                        end do
                    end if
                end do
            end if

            call clear_sidechain_library(sc_lib)
        end do
    end if

    ! ---------------------------------------------------------
    ! 2. DOMAIN CVs verarbeiten
    ! ---------------------------------------------------------
write(0,*) job%domain_libs(1)%filters

!pause
if (job%include_domain_cvs .and. allocated(job%domain_libs)) then
        write(*,*) 'Processing domain libraries...'
        write(*,*)  trim(filename)
        do i = 1, size(job%domain_libs)
            filename = job%domain_libs(i)%filename
            
            call read_domain_library(trim(filename), dom_lib)
            
            if (allocated(dom_lib%domain)) then
                do k = 1, dom_lib%ndo
                    include_it = .true.
                    
                    if (allocated(job%domain_libs(i)%filters)) then
                        if (k <= size(job%domain_libs(i)%filters)) then
                           
                            include_it = job%domain_libs(i)%filters(k)
                            write(0,*) 'Domain Index:', k, 'Name:', trim(dom_lib%domain(k)%domain_name), 'Filter-Val:', include_it
                        end if
                    end if

                    if (include_it) then
                        call make_domain_names(names, dom_lib%domain(k)%domain_name)
                        call get_atoms_domaine(job%istate, dom_lib, topologies, k, atoms)
                        call add_domain_cv_pair( &
                            plan         = plan, &
                            domid        = k, &
                            name         = names, &
                            point_group  = 'none', &
                            atoms        = atoms, &
                            cv           = dom_lib%domain(k)%cv &
                        )
                    end if
                end do
            end if

            call clear_domain_library(dom_lib)
        end do
    end if
    ! ---------------------------------------------------------
    ! 3. USER-DEFINED CVs verarbeiten
    ! ---------------------------------------------------------
    if (job%include_user_cvs .and. allocated(job%user_libs)) then
        write(*,*) 'Processing user-defined libraries...'
        
        do i = 1, size(job%user_libs)
            filename = job%user_libs(i)%filename
            
            call read_user_cv_library(trim(filename), user_lib)
            
            if (user_lib%n_cvs > 0 .and. allocated(user_lib%cvs)) then
                do k = 1, user_lib%n_cvs
                    include_it = .true.
                    if (allocated(job%user_libs(i)%filters)) then
                        if (k <= size(job%user_libs(i)%filters)) then
                            include_it = job%user_libs(i)%filters(k)
                        end if
                    end if

                    if (include_it) then
                        names(1) = trim(user_lib%cvs(k)%cv_name) // '_1'
                        names(2) = trim(user_lib%cvs(k)%cv_name) // '_2'

                        call add_user_cv_pair( &
                            plan   = plan, &
                            istate = user_lib%cvs(k)%istate, &
                            name   = names, &
                            atoms  = user_lib%cvs(k)%atom_indices, &
                            cv     = user_lib%cvs(k)%cv_values       &
                        )
                    end if
                end do
            end if

            call clear_user_cv_library(user_lib)
        end do
    end if

    write(*,*) 'CV plan workflow execution completed successfully.'
    if (.not. job%read_existing_cv_plan) then
        write(*,*) 'Written CV plan to: ', trim(job%cv_plan_filename)
        call write_cv_plan(plan, trim(job%cv_plan_filename))
    end if

end subroutine run_cv_plan_workflow


subroutine configure_sub_workflows()
    use iso_fortran_env, only: output_unit, error_unit
    
    implicit none

    character(len=32) :: input
    integer :: choice, ios
    logical :: running

    running = .true.

    do while (running)
        write(output_unit, *)
        write(output_unit, *) '======================================'
        write(output_unit, *) '      SUB-WORKFLOW CONFIGURATION MENU '
        write(output_unit, *) '======================================'
        write(output_unit, *) '1 = Configure Seed Workflow'
        write(output_unit, *) '2 = Configure Trajectory Settings'
        write(output_unit, *) '3 = Configure User-Defined CV Libraries'
        write(output_unit, *) '4 = Configure CV-Library Workflow'
        write(output_unit, *) '5 = Configure CV-Plan Workflow'
        write(output_unit, *) '6 = Configure Test-CV-Trajectory Workflow'
        write(output_unit, *) '7 = Configure Gaphical-Rep-VMD-CV Workflow'
        write(output_unit, *) '8 = Configure Replicate Selector Workflow'
        write(output_unit, *) '9 = Configure Euclidean Lifting Workflow'
        write(output_unit, *) '0 = Return to Main Menu'
        write(output_unit, *)
        write(output_unit, '(A)', advance='no') 'Select option [0]: '

        read(*,'(A)') input
        if (len_trim(input) == 0) then
            choice = 0
        else
            read(input, *, iostat=ios) choice
            if (ios /= 0) then
                write(error_unit, *) 'Invalid option. Please enter a number.'
                cycle
            end if
        end if

        select case (choice)
        case (0)
            running = .false.
            write(output_unit, *) '-> Exiting sub-workflow menu and returning to main driver.'

        case (1)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '           SEED WORKFLOW CONFIGURATION            '
            write(output_unit, *) '--------------------------------------------------'
            call configure_seed_workflow(initial_seed_workflow)

        case (2)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '         TRAJECTORY SETTINGS CONFIGURATION         '
            write(output_unit, *) '--------------------------------------------------'

            call configure_trajectory_workflow(traj_workflow,project,topologies)
           

        case (3)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '        USER CV LIBRARY CONFIGURATION             '
            write(output_unit, *) '--------------------------------------------------'
            call configure_user_cv_library(topologies)

        case (4)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '       CV-LIBRARY WORKFLOW CONFIGURATION          '
            write(output_unit, *) '--------------------------------------------------'
            call configure_library_workflow(cv_lib_workflow_use, topologies)

        case (5)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '        CV-PLAN WORKFLOW CONFIGURATION            '
            write(output_unit, *) '--------------------------------------------------'
            call configure_cv_plan_workflow(cv_plan_workflow1,topologies)
            
            
        case(6)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '        TEST CV TRAJECTORY WORKFLOW CONFIGURATION            '
            write(output_unit, *) '--------------------------------------------------'
            call configure_test_trajectory_workflow(write_test_tra)
        
         case(7)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '        GRAPHICAL REP VMD WORKFLOW CONFIGURATION            '
            write(output_unit, *) '--------------------------------------------------'
            call configure_graphical_representation_workflow(plot_vectors_vmd)

         case (8)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '       REPLICATE SELECTOR WORKFLOW CONFIGURATION'
            write(output_unit, *) '--------------------------------------------------'
            call configure_replicate_selector_workflow(replicate_selector_workflow, project)

         case (9)
            write(output_unit, *)
            write(output_unit, *) '--------------------------------------------------'
            write(output_unit, *) '       EUCLIDEAN LIFTING WORKFLOW CONFIGURATION'
            write(output_unit, *) '--------------------------------------------------'
            call configure_euclidean_lifting_workflow(euclidian_lifting_workflow)

         case default
         end select
    end do

    write(output_unit, *) '-> All interactive sub-workflow configurations completed.'
end subroutine configure_sub_workflows



subroutine run_workflow_driver()
    use iso_fortran_env, only: output_unit, error_unit
    use master_workflow_io_module, only: master_workflow_type, configure_master_workflow
    implicit none

    character(len=1024) :: input
    character(len=1024) :: master_cfg_filename
    integer :: main_choice, ios
    type(master_workflow_type) :: master
    character(len=1024) :: master_stack(64)
    character(len=1024) :: generated_selector_files(256)
    integer :: generated_selector_count
    logical :: running

    ! Initialize project and topologies once, before entering the menu loop.
    call initialize_system(project, topologies)

    running = .true.
    do while (running)
        write(output_unit, *)
        write(output_unit, *) '======================================'
        write(output_unit, *) '    PROTEIN WORKFLOW DRIVER SYSTEM    '
        write(output_unit, *) '======================================'
        write(output_unit, *) '1 = Configure sub-workflows interactively'
        write(output_unit, *) '2 = Configure master workflow only (using existing configs)'
        write(output_unit, *) '3 = Validate, Dry-Run and Execute master workflow (Batch)'
        write(output_unit, *) '4 = Quit'
        write(output_unit, *)
        write(output_unit, '(A)', advance='no') 'Select option [3]: '

        read(*,'(A)') input
        main_choice = 3
        if (len_trim(input) > 0) then
            read(input,*,iostat=ios) main_choice
            if (ios /= 0) then
                write(error_unit, '(A)') 'Invalid main option. Please enter a number from 1 to 4.'
                cycle
            end if
        end if

        select case (main_choice)

        case (1)
            write(output_unit, *) '-> Entering Sub-Workflow Configuration Mode...'
            call configure_sub_workflows()

        case (2)
            write(output_unit, *) '-> Entering Master-Workflow Configuration Mode...'
            call configure_master_workflow(master)
            write(output_unit, *) 'Master protocol configuration completed.'

        case (3)
            write(output_unit, *) '-> Entering Master Execution Mode...'

            write(output_unit, '(A)', advance='no') &
                'Master protocol filename to execute [master_workflow.cfg]: '
            read(*,'(A)') input
            if (len_trim(input) == 0) then
                master_cfg_filename = 'master_workflow.cfg'
            else
                master_cfg_filename = trim(adjustl(input))
            end if

            write(output_unit, *) '-> Running Dry-Run consistency check...'
            master_stack = ''
            generated_selector_files = ''
            generated_selector_count = 0
            call validate_master_tree(trim(master_cfg_filename), master_stack, 0, &
                generated_selector_files, generated_selector_count)
            write(output_unit, *) '-> Dry-run passed successfully! All configuration files are present.'

            write(output_unit, *) '-> Starting Execution Loop...'
            master_stack = ''
            call execute_master_tree(trim(master_cfg_filename), master_stack, 0)
            write(output_unit, *) 'All workflow steps completed successfully!'

        case (4)
            running = .false.
            write(output_unit, *) 'Exiting workflow driver.'

        case default
            write(output_unit, '(A,I0)') 'Invalid main option: ', main_choice
        end select
    end do
end subroutine run_workflow_driver

recursive subroutine validate_master_tree(master_filename, master_stack, depth, &
    generated_selector_files, generated_selector_count)
    use iso_fortran_env, only: error_unit
    use master_workflow_io_module, only: master_workflow_type, read_master_workflow, &
        validate_master_workflow, master_workflow_step_type_id, replicate_selector_workflow_step_type_id, &
        euclidean_lifting_workflow_step_type_id
    use replicate_selector_workflow_module, only: replicate_selector_workflow_type, &
        read_replicate_selector_workflow, validate_replicate_selector_workflow
    use euclidean_lifting_workflow_io_module, only: read_euclidean_lifting_workflow, &
        validate_euclidean_lifting_workflow
    implicit none

    character(len=*), intent(in) :: master_filename
    character(len=1024), intent(in) :: master_stack(:)
    integer, intent(in) :: depth
    character(len=1024), intent(inout) :: generated_selector_files(:)
    integer, intent(inout) :: generated_selector_count

    type(master_workflow_type) :: master
    character(len=1024) :: current_master, step_config, next_master_stack(size(master_stack))
    character(len=1024) :: selector_output
    integer :: i, generated_index
    logical :: file_exists, selector_will_be_generated

    current_master = normalize_workflow_path(master_filename)
    if (depth >= size(master_stack)) then
        write(error_unit, '(A)') 'ERROR: Maximum nested master workflow depth exceeded.'
        error stop 'Master workflow nesting is too deep.'
    end if
    do i = 1, depth
        if (trim(master_stack(i)) == trim(current_master)) then
            write(error_unit, '(A,A)') 'ERROR: Circular master workflow reference: ', trim(current_master)
            error stop 'Circular master workflow reference.'
        end if
    end do
    inquire(file=trim(current_master), exist=file_exists)
    if (.not. file_exists) then
        write(error_unit, '(A,A)') 'ERROR: Master workflow file not found: ', trim(current_master)
        error stop 'Missing master workflow file.'
    end if

    next_master_stack = master_stack
    next_master_stack(depth + 1) = current_master
    call read_master_workflow(master, trim(current_master))
    call validate_master_workflow(master)

    do i = 1, master%n_steps
        if (master%steps(i)%step_type < 1 .or. &
            master%steps(i)%step_type > euclidean_lifting_workflow_step_type_id) then
            write(error_unit, '(A,I0,A,A)') 'ERROR: Unknown workflow step type ', &
                master%steps(i)%step_type, ' in ', trim(current_master)
            error stop 'Unknown workflow step type.'
        end if

        step_config = resolve_master_step_path(current_master, trim(master%steps(i)%config_file))
        inquire(file=trim(step_config), exist=file_exists)
        if (.not. file_exists) then
            write(error_unit, '(A,A,A,A,A)') 'ERROR: Configuration file "', trim(step_config), &
                '" referenced by "', trim(current_master), '" was not found.'
            error stop 'Dry-run failed due to missing configuration files.'
        end if

        if (master%steps(i)%step_type == master_workflow_step_type_id) then
            call validate_master_tree(trim(step_config), next_master_stack, depth + 1, &
                generated_selector_files, generated_selector_count)
        else if (master%steps(i)%step_type == replicate_selector_workflow_step_type_id) then
            block
                type(replicate_selector_workflow_type) :: selector_job
                call read_replicate_selector_workflow(trim(step_config), selector_job)
                call validate_replicate_selector_workflow(selector_job, project)
                selector_output = resolve_master_step_path(step_config, trim(selector_job%selector_file))
                if (generated_selector_count >= size(generated_selector_files)) then
                    error stop 'Too many replicate selector outputs in master workflow.'
                end if
                generated_selector_count = generated_selector_count + 1
                generated_selector_files(generated_selector_count) = selector_output
            end block
        else if (master%steps(i)%step_type == euclidean_lifting_workflow_step_type_id) then
            block
                type(euclidean_lifting_workflow_type) :: lifting_job
                character(len=1024) :: data_file
                logical :: data_file_exists

                call read_euclidean_lifting_workflow(lifting_job, trim(step_config))
                call validate_euclidean_lifting_workflow(lifting_job)
                if (lifting_job%run_params%execute_lifting) then
                    data_file = resolve_master_step_path(step_config, &
                        trim(lifting_job%io_params%cv_plan_filename))
                    inquire(file=trim(data_file), exist=data_file_exists)
                    if (.not. data_file_exists) then
                        write(error_unit, '(A,A)') 'ERROR: CV plan file not found: ', trim(data_file)
                        error stop 'Dry-run failed due to missing CV plan file.'
                    end if
                    data_file = resolve_master_step_path(step_config, &
                        trim(lifting_job%io_params%selector_file))
                    inquire(file=trim(data_file), exist=data_file_exists)
                    selector_will_be_generated = .false.
                    do generated_index = 1, generated_selector_count
                        if (trim(generated_selector_files(generated_index)) == trim(data_file)) then
                            selector_will_be_generated = .true.
                            exit
                        end if
                    end do
                    if (.not. data_file_exists .and. .not. selector_will_be_generated) then
                        write(error_unit, '(A,A)') 'ERROR: Replicate selector file not found: ', trim(data_file)
                        error stop 'Dry-run failed due to missing replicate selector file.'
                    end if
                end if
            end block
        end if
    end do
end subroutine validate_master_tree

recursive subroutine execute_master_tree(master_filename, master_stack, depth)
    use iso_fortran_env, only: output_unit, error_unit
    use master_workflow_io_module, only: master_workflow_type, read_master_workflow, &
        validate_master_workflow, master_workflow_step_type_id, replicate_selector_workflow_step_type_id, &
        euclidean_lifting_workflow_step_type_id
    use random_seed_workflow_module, only: random_seed_workflow_type, read_seed_workflow, run_seed_workflow
    use trajectory_workflow_io_module, only: trajectory_workflow_type, read_trajectory_workflow, &
        run_trajectory_workflow
    use cv_lib_workflow_io_module, only: read_library_workflow, validate_library_workflow
    use replicate_selector_workflow_module, only: replicate_selector_workflow_type, &
        read_replicate_selector_workflow, validate_replicate_selector_workflow, &
        run_replicate_selector_workflow
    use replicate_selector_workflow_module, only: replica_selector_type
    implicit none

    character(len=*), intent(in) :: master_filename
    character(len=1024), intent(in) :: master_stack(:)
    integer, intent(in) :: depth

    type(master_workflow_type) :: master
    character(len=1024) :: current_master, step_config, next_master_stack(size(master_stack))
    integer :: i

    current_master = normalize_workflow_path(master_filename)
    if (depth >= size(master_stack)) then
        write(error_unit, '(A)') 'ERROR: Maximum nested master workflow depth exceeded.'
        error stop 'Master workflow nesting is too deep.'
    end if
    do i = 1, depth
        if (trim(master_stack(i)) == trim(current_master)) then
            write(error_unit, '(A,A)') 'ERROR: Circular master workflow reference: ', trim(current_master)
            error stop 'Circular master workflow reference.'
        end if
    end do

    next_master_stack = master_stack
    next_master_stack(depth + 1) = current_master
    call read_master_workflow(master, trim(current_master))
    call validate_master_workflow(master)

    do i = 1, master%n_steps
        step_config = resolve_master_step_path(current_master, trim(master%steps(i)%config_file))
        write(output_unit, '(A,A,A,I0,A)') '>>> Executing Step: ', trim(master%steps(i)%step_name), &
            ' [type ', master%steps(i)%step_type, ']'

        select case (master%steps(i)%step_type)
        case (1)
            block
                type(random_seed_workflow_type) :: seed_job
                call read_seed_workflow(trim(step_config), seed_job)
                call run_seed_workflow(seed_job)
            end block

        case (2)
            block
                type(trajectory_workflow_type) :: traj_job
                call read_trajectory_workflow(traj_job, trim(step_config))
                call run_trajectory_workflow(traj_job, project, topologies)
            end block

        case (3)
            block
                type(library_workflow_type) :: cv_job
                call read_library_workflow(cv_job, trim(step_config))
                call validate_library_workflow(cv_job)
                call run_library_workflow(cv_job, topologies, sc_lib, dom_lib, user_lib)
            end block

        case (4)
            block
                type(cv_plan_workflow_type) :: cv_plan_job
                type(cv_plan) :: generated_plan
                call read_cv_plan_workflow(cv_plan_job, trim(step_config))
                call run_cv_plan_workflow(cv_plan_job, generated_plan)
                write(output_unit, *) '-> CV Plan successfully generated from workflow.'
            end block

        case (5)
            block
                type(test_trajectory_workflow_type) :: test_job
                call read_test_trajectory_workflow(test_job, trim(step_config))
                call run_test_trajectory_workflow(test_job, project)
            end block

        case (6)
            block
                type(graphical_representation_workflow_type) :: graphical_job
                call read_graphical_representation_workflow(graphical_job, trim(step_config))
                call run_graphical_representation_workflow(graphical_job, project, topologies)
            end block

        case (replicate_selector_workflow_step_type_id)
            block
                type(replicate_selector_workflow_type) :: selector_job
                type(replica_selector_type) :: selector
                call read_replicate_selector_workflow(trim(step_config), selector_job)
                call validate_replicate_selector_workflow(selector_job, project)
                selector_job%selector_file = resolve_master_step_path(step_config, &
                    trim(selector_job%selector_file))
                call run_replicate_selector_workflow(selector_job, project, selector)
            end block

        case (euclidean_lifting_workflow_step_type_id)
            block
                type(euclidean_lifting_workflow_type) :: lifting_job
                call read_euclidean_lifting_workflow(lifting_job, trim(step_config))
                call validate_euclidean_lifting_workflow(lifting_job)
                lifting_job%io_params%cv_plan_filename = resolve_master_step_path( &
                    step_config, trim(lifting_job%io_params%cv_plan_filename))
                lifting_job%io_params%selector_file = resolve_master_step_path( &
                    step_config, trim(lifting_job%io_params%selector_file))
                lifting_job%io_params%lifting_plan_filename = resolve_master_step_path( &
                    step_config, trim(lifting_job%io_params%lifting_plan_filename))
                call run_euclidean_lifting_workflow(lifting_job, project, topologies)
            end block

        case (master_workflow_step_type_id)
            call execute_master_tree(trim(step_config), next_master_stack, depth + 1)

        case default
            write(error_unit, '(A,I0)') 'ERROR: Unknown workflow step type ID: ', master%steps(i)%step_type
            error stop 'Execution failed due to unknown step type.'
        end select

        write(output_unit, '(3A)') '<<< Finished Step: ', trim(master%steps(i)%step_name), ' [DONE]'
    end do

    write(output_unit, '(A,A)') 'Finished master workflow: ', trim(current_master)
end subroutine execute_master_tree

function resolve_master_step_path(master_filename, configured_path) result(resolved_path)
    implicit none

    character(len=*), intent(in) :: master_filename, configured_path
    character(len=1024) :: resolved_path, normalized_config, parent_dir, candidate
    logical :: file_exists
    integer :: slash_position

    normalized_config = trim(adjustl(configured_path))
    if (len_trim(normalized_config) == 0) error stop 'Empty workflow configuration path.'
    if (normalized_config(1:1) == '/') then
        resolved_path = normalize_workflow_path(trim(normalized_config))
        return
    end if

    slash_position = scan(trim(master_filename), '/', back=.true.)
    if (slash_position == 0) then
        parent_dir = '.'
    else if (slash_position == 1) then
        parent_dir = '/'
    else
        parent_dir = master_filename(:slash_position - 1)
    end if

    if (trim(parent_dir) == '.') then
        candidate = normalize_workflow_path(trim(normalized_config))
    else if (trim(parent_dir) == '/') then
        candidate = normalize_workflow_path('/' // trim(normalized_config))
    else
        candidate = normalize_workflow_path(trim(parent_dir) // '/' // trim(normalized_config))
    end if
    inquire(file=trim(candidate), exist=file_exists)
    if (file_exists) then
        resolved_path = candidate
        return
    end if

    candidate = normalize_workflow_path(trim(normalized_config))
    inquire(file=trim(candidate), exist=file_exists)
    if (file_exists) then
        resolved_path = candidate
    else
        resolved_path = normalize_workflow_path(trim(parent_dir) // '/' // trim(normalized_config))
    end if
end function resolve_master_step_path

function normalize_workflow_path(path) result(normalized_path)
    implicit none

    character(len=*), intent(in) :: path
    character(len=1024) :: normalized_path, source, components(256)
    integer :: source_length, position, component_end, component_count, i, output_length
    logical :: is_absolute

    source = trim(adjustl(path))
    source_length = len_trim(source)
    if (source_length == 0) error stop 'Cannot normalize an empty path.'
    is_absolute = source(1:1) == '/'
    component_count = 0
    position = 1

    do while (position <= source_length)
        if (source(position:position) == '/') then
            position = position + 1
            cycle
        end if
        component_end = position
        do while (component_end <= source_length)
            if (source(component_end:component_end) == '/') exit
            component_end = component_end + 1
        end do

        select case (source(position:component_end - 1))
        case ('.')
            continue
        case ('..')
            if (component_count > 0) then
                if (trim(components(component_count)) /= '..') then
                    component_count = component_count - 1
                else if (.not. is_absolute) then
                    if (component_count >= size(components)) error stop 'Workflow path has too many components.'
                    component_count = component_count + 1
                    components(component_count) = '..'
                end if
            else if (.not. is_absolute) then
                if (component_count >= size(components)) error stop 'Workflow path has too many components.'
                component_count = component_count + 1
                components(component_count) = '..'
            end if
        case default
            if (component_count >= size(components)) error stop 'Workflow path has too many components.'
            component_count = component_count + 1
            components(component_count) = source(position:component_end - 1)
        end select
        position = component_end + 1
    end do

    normalized_path = ''
    output_length = 0
    if (is_absolute) then
        normalized_path(1:1) = '/'
        output_length = 1
    end if
    do i = 1, component_count
        if (output_length > 0) then
            if (normalized_path(output_length:output_length) /= '/') then
                output_length = output_length + 1
                normalized_path(output_length:output_length) = '/'
            end if
        end if
        if (output_length + len_trim(components(i)) > len(normalized_path)) then
            error stop 'Normalized workflow path is too long.'
        end if
        normalized_path(output_length + 1:output_length + len_trim(components(i))) = trim(components(i))
        output_length = output_length + len_trim(components(i))
    end do
    if (output_length == 0) normalized_path = '.'
end function normalize_workflow_path


subroutine run_library_workflow(job, topologies, sc_lib, dom_lib, user_lib)
    use sidechains, only: sidechain_library_type, init_resnames_and_sidechain_head_groups, init_sidechain_head, compute_freq_and_resids, write_sidechain_library
    use domains_module, only: domain_library_type, init_domain_maps, write_domain_library
    use user_defined_cv_module, only: user_cv_library_type, read_user_cv_library, write_user_cv_library, clear_user_cv_library
    use sidechains, only: clear_sidechain_library
    use domains_module, only: clear_domain_library
    use topology_module, only: topology_set_type
    implicit none

    type(library_workflow_type), intent(in)       :: job
    type(topology_set_type),     intent(in)       :: topologies
    type(sidechain_library_type), intent(out)     :: sc_lib
    type(domain_library_type),    intent(out)     :: dom_lib
    type(user_cv_library_type),   intent(out)     :: user_lib
    type(time_window) :: window
    integer :: i, j, valid_count
    logical :: is_active
    logical, allocatable :: local_filters(:)
    integer :: filter_val
    logical, allocatable :: local_domain_filters(:)
integer :: domain_val
    write(*,*)
    write(*,*) '======================================'
    write(*,*) '      Executing Library Workflow      '
    write(*,*) '======================================'
    write(*,*)


    ! --------------------------------------------------
    ! 2. Sidechain-Bibliothek initialisieren & filtern
    ! --------------------------------------------------
    if (job%run_params%initialize_sidechain_cvs) then
        write(*,*) '-> Initializing and computing Sidechain library...'
        call init_resnames_and_sidechain_head_groups(sc_lib)
        call init_sidechain_head(sc_lib)
        call make_sc_res_maps(sc_lib,job%model_params%istate, topologies)
        call compute_freq_and_resids(sc_lib, job%model_params%istate, topologies)
        
        do i=1,sc_lib%ncv
        write(0,*) sc_lib%residue(i)%freq,'resids',sc_lib%residue(i)%resname
        end do
       
        
       write(0,*) 'Sidechain lib preinitialization done',job%io_params%read_sidechain_library
       if (job%io_params%read_sidechain_library) then
        call read_sidechain_library(trim(job%io_params%sidechain_library_file),sc_lib)
       end if
       
       write(0,*) shape(job%sc_type_filters)
       write(0,*) job%sc_type_filters
       write(0,*) sc_lib%ncv
       
      ! pause
! Workflow-Filter anwenden: Inaktive Typen komplett auf Frequenz 0 setzen
! Workflow-Filter anwenden: Inaktive Typen komplett auf Frequenz 0 setzen
if (allocated(job%sc_type_filters)) then
    if (size(job%sc_type_filters) == sc_lib%ncv) then
        allocate(local_filters(sc_lib%ncv))
        local_filters = job%sc_type_filters
        write(0,*) local_filters
        
        
  do i = 1, sc_lib%ncv
    ! Wir mappen den Filter direkt auf einen Integer
    if (local_filters(i)) then
        filter_val = 1
    else
        filter_val = 0
    end if

    !write(0,*) i, 'Filter-Int:', filter_val, 'out'

   
    if (filter_val == 0) then
    !    write(0,*) i, 'Filter-Int:', filter_val, 'in'
        sc_lib%residue(i)%freq = 0
        write(0,*) i, trim(sc_lib%residue(i)%resname), 'aussortiert'
        if (allocated(sc_lib%residue(i)%res_ids)) then
            deallocate(sc_lib%residue(i)%res_ids)
        end if
    end if
end do
        
        deallocate(local_filters)
    end if
end if
        
        !do i=1,sc_lib%ncv
        !write(0,*) sc_lib%residue(i)%freq,'resids',sc_lib%residue(i)%resname
        !end do
        

end if
    ! --------------------------------------------------
    ! 3. Domain-Bibliothek initialisieren & filtern (Nur aktive behalten)
    ! --------------------------------------------------
    if (job%run_params%initialize_domain_cvs) then
        write(*,*) '-> Initializing Domain library...'
        call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)
        
        if (job%io_params%read_domain_library) then
           call read_domain_library(trim(job%io_params%domain_library_file), dom_lib)
        end if

        ! Inaktive Domains herausfiltern und Array kompaktieren (mit robustem Integer-Check)
        if (allocated(job%domain_filters)) then
            if (size(job%domain_filters) == dom_lib%ndo) then
                
                ! Sicher lokal kopieren
                if (allocated(local_domain_filters)) deallocate(local_domain_filters)
                allocate(local_domain_filters(dom_lib%ndo))
                local_domain_filters = job%domain_filters

                ! Gültige Domains über robusten Integer-Check zählen
                valid_count = 0
                do i = 1, dom_lib%ndo
                    if (local_domain_filters(i)) then
                        domain_val = 1
                    else
                        domain_val = 0
                    end if

                    if (domain_val == 1) then
                        valid_count = valid_count + 1
                    end if
                end do

                if (valid_count < dom_lib%ndo) then
                    block
                        type(domain_library_type) :: temp_lib
                        if (valid_count > 0) then
                            allocate(temp_lib%domain(valid_count))
                            temp_lib%ndo = valid_count
                            
                            j = 0
                            do i = 1, dom_lib%ndo
                                if (local_domain_filters(i)) then
                                    domain_val = 1
                                else
                                    domain_val = 0
                                end if

                                if (domain_val == 1) then
                                    j = j + 1
                                    temp_lib%domain(j) = dom_lib%domain(i)
                                end if
                            end do

                            call clear_domain_library(dom_lib)
                            !call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)
                            dom_lib = temp_lib
                        else
                            call clear_domain_library(dom_lib)
                            dom_lib%ndo = 0
                        end if
                    end block
                end if

                deallocate(local_domain_filters)
            end if
        end if
    end if

!write(0,*) dom_lib%ndo, 'nach loop',job%domain_filters,shape(dom_lib%domain)
!pause
!--------------------------------------------------
! --------------------------------------------------
    ! 4. User-Defined CV Bibliothek einlesen & filtern (Direkt auf user_lib)
    ! --------------------------------------------------
    if (job%run_params%initialize_user_cvs) then
        write(*,'(3A)') '-> Reading User CV library from "', trim(job%io_params%user_cv_library_file), '"...'
        call read_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
        
        if (user_lib%n_cvs > 0 .and. allocated(user_lib%cvs)) then
            call init_user_cv_maps(topologies%state(job%model_params%istate), user_lib)

            ! Wenn Workflow-Filter vorliegen, wenden wir sie direkt an und schrumpfen das Array
            if (allocated(job%user_cv_filters) .and. size(job%user_cv_filters) == user_lib%n_cvs) then
                block
                    type(user_cv_type), allocatable :: temp_cvs(:)
                    integer :: valid_count, j

                    ! 1. Zählen, wie viele übrig bleiben
                    valid_count = 0
                    do i = 1, user_lib%n_cvs
                        if (job%user_cv_filters(i)) valid_count = valid_count + 1
                    end do

                    ! 2. Nur filtern, wenn auch wirklich welche ausgeschlossen wurden
                    if (valid_count < user_lib%n_cvs) then
                        if (valid_count > 0) then
                            allocate(temp_cvs(valid_count))
                            j = 0
                            do i = 1, user_lib%n_cvs
                                if (job%user_cv_filters(i)) then
                                    j = j + 1
                                    temp_cvs(j) = user_lib%cvs(i)
                                end if
                            end do
                            
                            ! Direkt in die echte user_lib verschieben (verhindert den Alloc-Crash)
                            call move_alloc(temp_cvs, user_lib%cvs)
                        else
                            ! Falls alle abgewählt wurden
                            if (allocated(user_lib%cvs)) deallocate(user_lib%cvs)
                        end if
                        
                        user_lib%n_cvs = valid_count
                    end if
                end block
            end if
        else
            user_lib%n_cvs = 0
        end if
    else
        user_lib%n_cvs = 0
    end if
   ! --------------------------------------------------
    ! 5. Optionale Referenzstruktursuche / Festlegung
    ! --------------------------------------------------
    if (job%run_params%search_reference_structures) then
        write(*,*) '-> Processing reference structures (Method: ', trim(job%model_params%reference_method), ')...'
        
        select case (trim(job%model_params%reference_method))
        case ('cv_orthogonality')
            window%istate          = job%reference_window%istate
            window%irep            = job%reference_window%irep
            window%start_time_step = job%reference_window%start_time_step
            window%end_time_step   = job%reference_window%end_time_step

            call ref_struc_search_user( &
                project, topologies, sc_lib, dom_lib, user_lib, &
                window, &
                job%model_params%n_candidates, &
                job%model_params%method_window)
                
        case ('rmsd')
            window%istate          = job%reference_window%istate
            window%irep            = job%reference_window%irep
            window%start_time_step = job%reference_window%start_time_step
            window%end_time_step   = job%reference_window%end_time_step

            call ref_struc_search_rmsd_user( &
                project, topologies, sc_lib, dom_lib, user_lib, &
                window, &
                job%model_params%n_candidates, &
                job%model_params%method_window)

        case ('fixed_frame')
            call ref_struc_fixed_user( &
                project, topologies, sc_lib, dom_lib, user_lib, &
                job%reference_window%istate, &
                job%reference_window%irep, &
                job%reference_window%start_time_step)
                
        case default
            error stop 'Invalid or unsupported reference method in workflow job'
        end select
        
        write(*,*) '-> Reference structure setup completed successfully.'
    end if

    ! --------------------------------------------------
    ! 6. Finales Abspeichern aller Bibliotheken (Allerletzter Schritt)
    ! --------------------------------------------------
    if (job%io_params%write_sidechain_library) then
        write(*,'(3A)') '-> Writing final Sidechain library to "', trim(job%io_params%sidechain_library_file), '"...'
        call write_sidechain_library(trim(job%io_params%sidechain_library_file), sc_lib)
    end if

    if (job%io_params%write_domain_library) then
        write(*,'(3A)') '-> Writing final Domain library to "', trim(job%io_params%domain_library_file), '"...'
        call write_domain_library(trim(job%io_params%domain_library_file), dom_lib)
    end if

    if (job%io_params%write_user_cv_library .and. user_lib%n_cvs > 0) then
        write(*,'(3A)') '-> Writing final User CV library to "', trim(job%io_params%user_cv_library_file), '"...'
        call write_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
    end if


    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Library Workflow Execution Done!   '
    write(*,*) '======================================'
    write(*,*)

end subroutine run_library_workflow








subroutine run_euclidean_lifting_workflow(job, project, topologies)
    use workflow_module           ! Passe den Modulnamen ggf. an
    use project_module
    use topology_module
    use replicate_selector_workflow_module, only: replica_selector_type, read_replica_selector_file
    use cv_plan_io_module, only: read_cv_plan
    ! use dein_cv_plan_module          ! Falls der cv_plan einen eigenen Typ hat
    
    implicit none

    type(euclidean_lifting_workflow_type), intent(in) :: job
    type(project_type), intent(inout) :: project
    type(topology_set_type), intent(in) :: topologies
    type(cv_plan) :: cv_p

   
    type (refstruc_type)         :: ref1      ! Referenzstruktur für PBC
    type(frame2plot) :: test_frame             ! Für den optionalen PBC-Test
    type(lifting_plan_type) :: lifting_plan
    type(replica_selector_type) :: selector
    integer :: selector_state, irep, active_replicates

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Executing Euclidean Lifting Pipeline'
    write(*,*) '======================================'

    if (.not. allocated(project%state)) error stop 'Project has no states for lifting'
    if (job%run_params%lift_istate < 1 .or. job%run_params%lift_istate > project%nstate) then
        error stop 'Selected lifting state is outside the project'
    end if
    if (size(project%state) < project%nstate) error stop 'Project state array is inconsistent'
    if (project%state(job%run_params%lift_istate)%nrep < 1 .or. &
        .not. allocated(project%state(job%run_params%lift_istate)%rep)) then
        error stop 'Selected lifting state has no replicates'
    end if
    if (size(project%state(job%run_params%lift_istate)%rep) /= &
        project%state(job%run_params%lift_istate)%nrep) then
        error stop 'Selected lifting state has an inconsistent replicate array'
    end if

    ! --------------------------------------------------
    ! 1. Optionale PBC-Korrektur vorbereiten und testen
    ! --------------------------------------------------
    if (job%run_params%use_pbc_correction) then
        write(*,*) 'Preparing PBC correction...'
        if (job%run_params%ref_irep > project%state(job%run_params%lift_istate)%nrep) then
            error stop 'PBC reference replicate is outside the selected state'
        end if
        if (job%run_params%ref_itime > &
            project%state(job%run_params%lift_istate)%rep(job%run_params%ref_irep)%n_frames_use_total) then
            error stop 'PBC reference time is outside the selected replicate'
        end if
        if (job%run_params%test_pbc) then
            if (job%run_params%pbc_irep > project%state(job%run_params%lift_istate)%nrep) then
                error stop 'PBC test replicate is outside the selected state'
            end if
            if (job%run_params%pbc_itime > &
                project%state(job%run_params%lift_istate)%rep(job%run_params%pbc_irep)%n_frames_use_total) then
                error stop 'PBC test time is outside the selected replicate'
            end if
        end if
        
        ! State ist fix über job%run_params%lift_istate definiert, 
        ! Replik und Timestep kommen aus den Parametern
        call prepare_pbc_correction( &
            job%run_params%lift_istate, &
            job%run_params%ref_irep, &
            job%run_params%ref_itime, &
            project, ref1)

        if (job%run_params%test_pbc) then
            write(*,*) 'Testing PBC correction...'
            
            ! Test-Frame zusammenbauen
            test_frame%istate = job%run_params%lift_istate
            test_frame%irep   = job%run_params%pbc_irep
            test_frame%itime  = job%run_params%pbc_itime

            call test_pbc_correction(project, topologies, test_frame, ref1)
        end if
    else
        write(*,*) 'PBC correction skipped (disabled in workflow).'
    end if

    ! --------------------------------------------------
    ! 2. Lifting Plan initialisieren
    ! --------------------------------------------------
    write(*,*) 'Initializing lifting plan...'
    call init_lifting_plan_state( job%run_params%lift_istate,&
        lifting_plan, &
        job%run_params%offset, &
        job%run_params%stride, &
        job%run_params%size_block, &
        project)

    ! --------------------------------------------------
    ! 3. Euclidean Lifting ausführen (mit 100 GB Trajektorien-Pfaden)
    ! --------------------------------------------------
    if (job%run_params%execute_lifting) then
        write(*,*) 'Executing lifting '

        call read_cv_plan(cv_p, trim(job%io_params%cv_plan_filename))
        call read_replica_selector_file(trim(job%io_params%selector_file), selector, selector_state)
        if (selector_state /= job%run_params%lift_istate) then
            error stop 'Replicate selector state does not match lifting state'
        end if
        if (selector%nrep_total /= project%state(job%run_params%lift_istate)%nrep) then
            error stop 'Replicate selector size does not match lifting state'
        end if
        if (selector%nrep_active < 1) error stop 'Replicate selector contains no active replicates'

        active_replicates = 0
        do irep = 1, lifting_plan%nreps(job%run_params%lift_istate)
            if (.not. selector%is_active(irep)) then
                lifting_plan%nblock(job%run_params%lift_istate, irep) = 0
            else
                active_replicates = active_replicates + 1
            end if
        end do
        write(*,'(A,I0,A,I0)') 'Lifting active replicates: ', active_replicates, ' of ', selector%nrep_total
        call write_lifting_plan(lifting_plan, trim(job%io_params%lifting_plan_filename))
        write(*,'(3A)') 'Lifting plan saved to "', trim(job%io_params%lifting_plan_filename), '"'
        
        call execute_lifting_plan_centroid( &
            cv_p                 = cv_p, &
            lifting_plan         = lifting_plan, &
            project              = project, &
            ref                  = ref1, & ! Wird ignoriert/muss übergeben werden je nach Implementierung
            path_lifted_tra      = trim(job%io_params%lifted_path)//'/', &
            prefix_lifted_tra    = trim(job%io_params%lifted_prefix), &
            path_centroid_tra    = trim(job%io_params%centroid_path)//'/', &
            prefix_centroid_tra  = trim(job%io_params%centroid_prefix), &
            output_policy        = trim(job%io_params%output_policy), &
            use_pbc_correction   = job%run_params%use_pbc_correction)

        write(*,*) 'Euclidean lifting finished successfully.'
    else
        write(*,*) 'Execution skipped (execute_lifting = .false.).'
        write(*,*) 'Lifting plan and setup completed successfully.'
    end if

end subroutine run_euclidean_lifting_workflow




subroutine project_setting(project)
implicit none

integer :: i,j,k,ndcd_dummy

type(project_type), intent(inout) :: project


do i=1,project%nstate

           
    do j=1,project%state(i)%nrep
    ndcd_dummy=project%state(i)%rep(j)%ndcd
    write(0,*) ndcd_dummy
     
    allocate(project%state(i)%rep(j)%ignore_first_per_dcd(ndcd_dummy))
    allocate(project%state(i)%rep(j)% ignore_last_per_dcd(ndcd_dummy))
    allocate(project%state(i)%rep(j)% nframes_per_dcd_raw(ndcd_dummy))
    allocate(project%state(i)%rep(j)% nframes_per_dcd_use(ndcd_dummy))

        do k=1,project%state(i)%rep(j)%ndcd
        
             project%state(i)%rep(j)%ignore_first_per_dcd(k) = 1
             project%state(i)%rep(j)%ignore_last_per_dcd(k)  = 1
             project%state(i)%rep(j)%offset_start=0!100010!100000
             project%state(i)%rep(j)%offset_end  =0
             project%state(i)%rep(j)%nlast_natom_ignore=0
             project%state(i)%rep(j)%nfirst_natom_to_read&
             = topologies%state(i)%topo%natom
        end do

    end do

end do

end subroutine





subroutine test_topology(topo)

    use topology_module
    implicit none

    type(topology_type), intent(in) :: topo

    integer :: i
    integer :: nshow

    write(0,*) ''
    write(0,*) '================================'
    write(0,*) 'TOPOLOGY DEBUG TEST'
    write(0,*) '================================'

    !====================================================
    ! BASIC INFO
    !====================================================
    write(0,*) 'natom =', topo%natom
    write(0,*) 'nres  =', topo%nres

    if (topo%natom <= 0) then
        write(0,*) 'ERROR: no atoms loaded'
        return
    end if

    !====================================================
    ! ATOM TABLE (first entries)
    !====================================================
    write(0,*) ''
    write(0,*) 'First atoms:'

    nshow = min(20, topo%natom)

    do i = 1, nshow
        write(0,'(I6,2X,A4,2X,A3,2X,I6,2X,3F10.3)') &
            i, &
            topo%atom_id(i), &
            topo%acid_id(i), &
            topo%acid_number(i), &
            topo%refx(i), topo%refy(i), topo%refz(i)
    end do

    !====================================================
    ! RESIDUE CHECK
    !====================================================
    write(0,*) ''
    write(0,*) 'Residue mapping check:'

    do i = 1, nshow
        write(0,*) 'atom', i, '-> residue', topo%residx(i)
    end do

    !====================================================
    ! CONSISTENCY CHECKS
    !====================================================
    write(0,*) ''
    write(0,*) 'Consistency checks:'

    if (.not. allocated(topo%atom_id)) then
        write(0,*) 'ERROR: atom_id not allocated'
    end if

    if (.not. allocated(topo%residx)) then
        write(0,*) 'ERROR: residx not allocated'
    end if

    if (size(topo%residx) /= topo%natom) then
        write(0,*) 'WARNING: residx size mismatch'
    end if

    write(0,*) ''
    write(0,*) '================================'
    write(0,*) 'DEBUG TEST DONE'
    write(0,*) '================================'

end subroutine test_topology





end program main