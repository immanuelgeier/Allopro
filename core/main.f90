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
!integer :: irep,istate

logical :: found
integer :: n, clock,i
logical, allocatable :: matrix_active_final(:,:)

type(random_seed_workflow_type) ::seed_job, initial_seed_workflow,prodcution_seed_workflow
logical :: run


!call configure_euclidean_lifting_workflow(euclidian_lifting_workflow)
!stop
!call run_workflow_driver()
call configure_cv_plan_workflow(cv_plan_workflow1)

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

call configure_cv_plan_workflow(cv_plan_workflow1)



stop


write(0,*) 'EUCLIDIAN LIFTING WORKFLWO CONFUGURATION:'
call configure_euclidean_lifting_workflow(euclidian_lifting_workflow)



if(run) then
call run_euclidean_lifting_workflow(euclidian_lifting_workflow, project, topologies, cv_plan1)	
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

subroutine configure_sub_workflows()
    use iso_fortran_env, only: output_unit, error_unit
    
    implicit none

    character(len=32) :: input
    integer :: counter
    logical :: continue_loop

    ! --------------------------------------------------
    ! 1. SEED WORKFLOWS SCHLEIFE
    ! --------------------------------------------------
    write(output_unit, *)
    write(output_unit, *) '======================================'
    write(output_unit, *) '       SEED WORKFLOW CONFIGURATION    '
    write(output_unit, *) '======================================'
    
    counter = 0
    continue_loop = .true.

    do while (continue_loop)
        counter = counter + 1
        write(output_unit, '(A,I0)') '--- Configuring Seed Workflow #', counter
        
        ! Ruft die Routine auf, die intern auch das Speichern übernimmt
        call configure_seed_workflow(initial_seed_workflow)

        write(output_unit, '(A)', advance='no') 'Configure another Seed workflow? (y/n) [n]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'n'
        input = trim(adjustl(input))
        continue_loop = (input == 'y' .or. input == 'Y' .or. input == 'yes')
    end do


    ! --------------------------------------------------
    ! 2. SYSTEM & PROJEKT INITIALISIERUNG
    ! --------------------------------------------------
    write(output_unit, *)
    write(output_unit, *) '======================================'
    write(output_unit, *) '     SYSTEM & PROJECT INITIALIZATION  '
    write(output_unit, *) '======================================'
    write(output_unit, *) '-> Initializing topologies for librarys'
    
     call initialize_system(project, topologies)
     call initialize_trajectory_settings(project, topologies)
     !call analyze_project_frames(project)
     !call print_project_summary(project)
     !
     !call init_resnames_and_sidechain_head_groups(sc_lib)
     !call init_sidechain_head(sc_lib)
     !
     !call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)

     
  write(output_unit, *)
    write(output_unit, *) '======================================'
    write(output_unit, *) '    USER CV LIBRARY CONFIGURATION     '
    write(output_unit, *) '======================================'

    ! 1. Erstmal fragen, ob überhaupt User-CV-Libraries erstellt werden sollen
    write(output_unit, '(A)', advance='no') 'Do you want to create user-defined CV libraries? (y/n) [n]: '
    read(*,'(A)') input
    if (len_trim(input) == 0) input = 'n'
    input = trim(adjustl(input))

    if (input == 'y' .or. input == 'Y' .or. input == 'yes') then
        counter = 0
        continue_loop = .true.

        do while (continue_loop)
            counter = counter + 1
            write(output_unit, '(A,I0)') '--- Creating User CV Library Set #', counter
            
            ! Aufruf deiner Routine
            call configure_user_cv_library(topologies)
            
            write(output_unit, '(A)', advance='no') 'Create another set of user CV libraries? (y/n) [n]: '
            read(*,'(A)') input
            if (len_trim(input) == 0) input = 'n'
            input = trim(adjustl(input))
            continue_loop = (input == 'y' .or. input == 'Y' .or. input == 'yes')
        end do
    else
        write(output_unit, *) '-> Skipping User CV libraries creation.'
    end if

    ! --------------------------------------------------
    ! 3. CV-LIBRARY WORKFLOWS SCHLEIFE
    ! --------------------------------------------------
    write(output_unit, *)
    write(output_unit, *) '======================================'
    write(output_unit, *) '    CV-LIBRARY WORKFLOW CONFIGURATION '
    write(output_unit, *) '======================================'
    
    counter = 0
    continue_loop = .true.

    do while (continue_loop)
        counter = counter + 1
        write(output_unit, '(A,I0)') '--- Configuring CV-Library Workflow #', counter
        
        call configure_library_workflow(cv_lib_workflow_use, topologies)

        write(output_unit, '(A)', advance='no') 'Configure another CV-Library workflow? (y/n) [n]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'n'
        input = trim(adjustl(input))
        continue_loop = (input == 'y' .or. input == 'Y' .or. input == 'yes')
    end do


    ! --------------------------------------------------
    ! 4. CV-PLAN WORKFLOWS SCHLEIFE
    ! --------------------------------------------------
    write(output_unit, *)
    write(output_unit, *) '======================================'
    write(output_unit, *) '    CV-PLAN WORKFLOW CONFIGURATION    '
    write(output_unit, *) '======================================'
    
    counter = 0
    continue_loop = .true.

    do while (continue_loop)
        counter = counter + 1
        write(output_unit, '(A,I0)') '--- Configuring CV-Plan Workflow #', counter
        
         call configure_cv_plan_workflow(cv_plan_workflow1)

        write(output_unit, '(A)', advance='no') 'Configure another CV-Plan workflow? (y/n) [n]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'n'
        input = trim(adjustl(input))
        continue_loop = (input == 'y' .or. input == 'Y' .or. input == 'yes')
    end do

    write(output_unit, *) '-> All interactive sub-workflow configurations completed.'
end subroutine configure_sub_workflows

subroutine run_workflow_driver()
    use iso_fortran_env, only: output_unit, error_unit
    use master_workflow_io_module
    use random_seed_workflow_module, only: random_seed_workflow_type, read_seed_workflow, run_seed_workflow
    use trajectory_workflow_io_module, only: trajectory_workflow_type, read_trajectory_workflow, run_trajectory_workflow
    use cv_lib_workflow_io_module ! Angenommen, hier liegt der Library-Workflow
    implicit none

    character(len=1024) :: input
    character(len=1024) :: master_cfg_filename
    integer :: main_choice, ios
    type(master_workflow_type) :: master
    integer :: i
    logical :: config_exists

    ! --- HAUPTMENÜ ANZEIGEN ---
    write(output_unit, *)
    write(output_unit, *) '======================================'
    write(output_unit, *) '    PROTEIN WORKFLOW DRIVER SYSTEM    '
    write(output_unit, *) '======================================'
    write(output_unit, *) '1 = Configure sub-workflows interactively'
    write(output_unit, *) '2 = Configure master workflow only (using existing configs)'
    write(output_unit, *) '3 = Validate, Dry-Run and Execute master workflow (Batch)'
    write(output_unit, *)
    write(output_unit, '(A)', advance='no') 'Select option [3]: '
    
    read(*,'(A)') input
    main_choice = 3
    if (len_trim(input) > 0) then
        read(input,*,iostat=ios) main_choice
        if (ios /= 0) error stop 'Invalid main option selected'
    end if

    ! --- AUSWAHL VERARBEITEN ---
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
        
        ! Master-Protokoll Dateinamen interaktiv abfragen
        write(output_unit, '(A)', advance='no') 'Master protocol filename to execute [master_workflow.cfg]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) then
            master_cfg_filename = 'master_workflow.cfg'
        else
            master_cfg_filename = trim(adjustl(input))
        end if

        ! 0. Projekt und Topologien einmalig initialisieren, damit nstate/nrep bekannt sind
        call initialize_system(project, topologies)

        ! 1. Master-Protokoll einlesen
        call read_master_workflow(master, trim(master_cfg_filename))

        ! 2. Dry-Run / Konsistenzprüfung (Prüfung aller Config-Dateien)
        write(output_unit, *) '-> Running Dry-Run consistency check...'
        do i = 1, master%n_steps
            inquire(file=trim(master%steps(i)%config_file), exist=config_exists)
            if (.not. config_exists) then
                write(error_unit, '(3A)') 'ERROR: Configuration file "', trim(master%steps(i)%config_file), '" was not found!'
                error stop 'Dry-run failed due to missing configuration files.'
            end if
        end do
        write(output_unit, *) '-> Dry-run passed successfully! All configuration files are present.'

        ! 3. Ausführungs-Loop über die numerischen IDs (ohne die Master-Datei zu überschreiben)
        write(output_unit, *) '-> Starting Execution Loop...'
        do i = 1, master%n_steps
            write(output_unit, '(3A)') '>>> Executing Step: ', trim(master%steps(i)%step_name), ' ...'

            ! Sauberer Dispatch über die eindeutige numerische ID
            select case (master%steps(i)%step_type)
            
            case (1) ! seed
                block
                    type(random_seed_workflow_type) :: seed_job
                    call read_seed_workflow(trim(master%steps(i)%config_file), seed_job)
                    call run_seed_workflow(seed_job)
                end block

            case (2) ! traj_settings
                block
                    type(trajectory_workflow_type) :: traj_job
                    
                    call read_trajectory_workflow(traj_job, trim(master%steps(i)%config_file))
                    call run_trajectory_workflow(traj_job, project,topologies)
                end block

            case (3) ! cv_lib
                block
                    type(library_workflow_type) :: cv_job
                    call read_library_workflow(cv_job, trim(master%steps(i)%config_file))
                    call validate_library_workflow(cv_job)
                    call run_library_workflow(cv_job, topologies, sc_lib, dom_lib, user_lib)
                end block

            case (4) ! cv_plan
                block
                    ! Hier den passenden Aufruf für cv_plan / lifting einsetzen
                end block

            case (5) ! lifting
                block
                    ! Entsprechender Workflow-Aufruf für lifting
                end block

            case default
                write(error_unit, '(A,I0)') 'ERROR: Unknown workflow step type ID: ', master%steps(i)%step_type
                error stop 'Execution failed due to unknown step type.'
            end select

            write(output_unit, '(3A)') '<<< Finished Step: ', trim(master%steps(i)%step_name), ' [DONE]'
        end do

        write(output_unit, *) 'All workflow steps completed successfully!'
        
        
    case default
        error stop 'Invalid option selected in workflow driver.'

    end select
end subroutine run_workflow_driver


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

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '      Executing Library Workflow      '
    write(*,*) '======================================'
    write(*,*)


    ! --------------------------------------------------
    ! 2. Sidechain-Bibliothek initialisieren & filtern
    ! --------------------------------------------------
    if (job%run_params%initialize_sidechain_maps) then
        write(*,*) '-> Initializing and computing Sidechain library...'
        call init_resnames_and_sidechain_head_groups(sc_lib)
        call init_sidechain_head(sc_lib)
        call  make_sc_res_maps(sc_lib,job%model_params%istate, topologies)
        call compute_freq_and_resids(sc_lib, job%model_params%istate, topologies)
        
        write(0,*) 'sidechainlib preinitialization done'
        pause
        
        
! Workflow-Filter anwenden: Inaktive Typen komplett auf Frequenz 0 setzen
        if (allocated(job%sc_type_filters)) then
            if (size(job%sc_type_filters) == sc_lib%ncv) then
                do i = 1, sc_lib%ncv
                    if (.not. job%sc_type_filters(i)) then
                        sc_lib%residue(i)%freq = 0
                        if (allocated(sc_lib%residue(i)%res_ids)) then
                            deallocate(sc_lib%residue(i)%res_ids)
                        end if
                    end if
                end do
            end if
        end if
        end if

    ! --------------------------------------------------
    ! 3. Domain-Bibliothek initialisieren & filtern (Nur aktive behalten)
    ! --------------------------------------------------
    if (job%run_params%initialize_domain_maps) then
        write(*,*) '-> Initializing Domain library...'
        call init_domain_maps(topologies%state(job%model_params%istate), dom_lib)

! Inaktive Domains herausfiltern und Array kompaktieren
        if (allocated(job%domain_filters)) then
            if (size(job%domain_filters) == dom_lib%ndo) then
                block
                    type(domain_library_type) :: temp_lib

                    valid_count = 0
                    do i = 1, dom_lib%ndo
                        if (job%domain_filters(i)) valid_count = valid_count + 1
                    end do

                    if (valid_count < dom_lib%ndo) then
                        allocate(temp_lib%domain(valid_count))
                        temp_lib%ndo = valid_count
                        
                        j = 0
                        do i = 1, dom_lib%ndo
                            if (job%domain_filters(i)) then
                                j = j + 1
                                temp_lib%domain(j) = dom_lib%domain(i)
                            end if
                        end do

                        call clear_domain_library(dom_lib)
                        dom_lib = temp_lib
                    end if
                end block
            end if
        end if
    end if


    ! --------------------------------------------------
    ! 4. User-Defined CV Bibliothek einlesen & filtern (Nur aktive behalten)
    ! --------------------------------------------------
    if (job%io_params%read_user_cv_library) then
        write(*,'(3A)') '-> Reading User CV library from "', trim(job%io_params%user_cv_library_file), '"...'
        call read_user_cv_library(trim(job%io_params%user_cv_library_file), user_lib)
        call init_user_cv_maps(topologies%state(job%model_params%istate), user_lib)

        if (allocated(job%user_cv_filters) .and. size(job%user_cv_filters) == user_lib%n_cvs) then
            block
                type(user_cv_library_type) :: temp_user_lib

                valid_count = 0
                do i = 1, user_lib%n_cvs
                    if (job%user_cv_filters(i)) valid_count = valid_count + 1
                end do

                if (valid_count < user_lib%n_cvs) then
                    allocate(temp_user_lib%cvs(valid_count))
                    temp_user_lib%n_cvs = valid_count
                    
                    j = 0
                    do i = 1, user_lib%n_cvs
                        if (job%user_cv_filters(i)) then
                            j = j + 1
                            temp_user_lib%cvs(j) = user_lib%cvs(i)
                        end if
                    end do

                    call clear_user_cv_library(user_lib)
                    user_lib = temp_user_lib
                end if
            end block
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





subroutine build_cv_plan_from_workflow(job, istate, topologies, plan)
    implicit none    

    type(cv_plan_workflow_type), intent(in) :: job
    integer, intent(in) :: istate
    type(topology_set_type), intent(in) :: topologies
    type(cv_plan), intent(out) :: plan

    type(sidechain_library_type) :: sc_lib
    type(domain_library_type)    :: dom_lib
    integer :: i

   ! ! 1. Sidechain-Bibliothek einlesen, falls im Workflow aktiviert
   ! if (job%run_params%include_sidechain_cvs) then
   !     !call read_sidechain_library(trim(job%io_params%sidechain_library_file), sc_lib)
   !     
   !     ! Schleife über Einträge/Frequenzen in sc_lib und Hinzufügen zum Plan
   !     ! do i = 1, ...
   !     !     call add_cv_pair_sidechain(...)
   !     ! end do
   ! end if

    ! 2. Domain-Bibliothek einlesen, falls im Workflow aktiviert
   ! if (job%run_params%include_domain_cvs) then
   !    ! call read_domain_library(trim(job%io_params%domain_library_file), dom_lib)
   !     
   !     ! do i = 1, ...
   !     !     call add_cv_pair_domaine(...)
   !     ! end do
   ! end if

end subroutine build_cv_plan_from_workflow



subroutine run_euclidean_lifting_workflow(job, project, topologies, cv_p)
    use workflow_module           ! Passe den Modulnamen ggf. an
    use project_module
    use topology_module
    ! use dein_cv_plan_module          ! Falls der cv_plan einen eigenen Typ hat
    
    implicit none

    type(euclidean_lifting_workflow_type), intent(in) :: job
    type(project_type), intent(inout) :: project
    type(topology_set_type), intent(in) :: topologies
    type(cv_plan), intent(in)     :: cv_p  ! Typ für deinen cv_plan1

   
    type (refstruc_type)         :: ref1      ! Referenzstruktur für PBC
    type(frame2plot) :: test_frame             ! Für den optionalen PBC-Test
    type(lifting_plan_type) :: lifting_plan

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Executing Euclidean Lifting Pipeline'
    write(*,*) '======================================'

    ! --------------------------------------------------
    ! 1. Optionale PBC-Korrektur vorbereiten und testen
    ! --------------------------------------------------
    if (job%run_params%use_pbc_correction) then
        write(*,*) 'Preparing PBC correction...'
        
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
        
        call execute_lifting_plan_centroid( &
            cv_p                 = cv_p, &
            lifting_plan         = lifting_plan, &
            project              = project, &
            ref                  = ref1, & ! Wird ignoriert/muss übergeben werden je nach Implementierung
            path_lifted_tra      = trim(job%io_params%lifted_path)//'/', &
            prefix_lifted_tra    = trim(job%io_params%lifted_prefix), &
            path_centroid_tra    = trim(job%io_params%centroid_path)//'/', &
            prefix_centroid_tra  = trim(job%io_params%centroid_prefix), &
            output_policy        = trim(job%io_params%output_policy))

        write(*,*) 'Euclidean lifting finished successfully.'
    else
        write(*,*) 'Execution skipped (execute_lifting = .false.).'
        write(*,*) 'Lifting plan and setup completed successfully.'
    end if

end subroutine run_euclidean_lifting_workflow




        
subroutine find_common_ref(istat,irep,sc_lib_opt,dom_lib_opt)

implicit none
integer,intent(in) :: istat,irep
type(time_window) :: window
integer :: start_search,end_search,stride_search
integer :: itest,nttest,ct,i,j1(1)
real(8), allocatable :: stats_b(:),stats_s(:)
real(8), allocatable :: all_stats_b(:,:),all_stats_s(:,:)
real(8) :: sum_res
type(sidechain_library_type) ,intent(inout):: sc_lib_opt
type(domain_library_type)    ,intent(inout):: dom_lib_opt
integer :: isc,idom,cc,ca,j
real(8), allocatable  :: min_sc(:),min_do(:)
integer, allocatable  :: mini_sc(:),mini_do(:)


start_search =1
end_search   =100000
stride_search=1000

nttest=(end_search-start_search+1)
nttest=nttest/stride_search


call ref_from_frames(1,1,1,sc_lib,topologies,project)


call init_domain_maps(topologies%state(1),dom_lib)
call ref_from_frames_domains(1,1,1,dom_lib,topologies,project)
ca=0
do isc=1,sc_lib%ncv
do i=1,min(sc_lib%residue(isc)%freq,1)
ca=ca+1
end do
end do
write(0,*)'number of occurrent res',ca
pause

call test_simple_vector_cv_statisics(topologies,sc_lib,dom_lib,project,stats_b,stats_s)

!write(0,*) nttest,size(stats_b),size(stats_s)
!pause
allocate(all_stats_b(nttest,size(stats_b)))
allocate(all_stats_s(nttest,size(stats_s)))
allocate(min_do(dom_lib%ndo))
allocate(min_sc(ca))
allocate(mini_do(size(min_do)))
allocate(mini_sc(size(min_sc)))

ct=0
min_sc=100000._8
min_do=100000._8
do itest=start_search,end_search,stride_search
     ct=ct+1
     call ref_from_frames(itest,1,1,sc_lib,topologies,project)
     
     !call test_get_ref_from_frames(sc_lib,topologies,project)
     
     
     !call ref_from_frames_domains(1,dom_lib,topologies,project)
     
     call ref_from_frames_domains(itest,1,1,dom_lib,topologies,project)
     
     call test_simple_vector_cv_statisics(topologies,sc_lib,dom_lib,project,stats_b,stats_s)
     all_stats_b(ct,:)=stats_b
     all_stats_s(ct,:)=stats_s
     
     cc=0
     ca=0
     do isc=1,sc_lib%ncv
        sum_res=0._8
        do i=1,sc_lib%residue(isc)%freq
            if(i==1) ca=ca+1
            j=sc_lib%residue(isc)%res_ids(i)
            cc=cc+1
            sum_res=sum_res+stats_s(cc)
            write(0,*)ca, cc,stats_s(cc),sc_lib%residue(isc)%resname!,sum_res,'sum',min_sc(ca)>sum_res,CA
        if(min_sc(ca)>sum_res) then
           sc_lib_opt%residue(isc)=sc_lib%residue(isc)
        min_sc(ca) =sum_res
        mini_sc(ca)=ct
        end if
	    end do
     end do
     
     do i=1,dom_lib%ndo
      if(stats_b(i)<min_do(i)) then
      	dom_lib_opt%domain(i)=dom_lib%domain(i)
      	mini_do(i)=ct
      end if
     end do
end do


ca=0
     do isc=1,sc_lib%ncv
        do i=1,min(sc_lib%residue(isc)%freq,1) 
            ca=ca+1   
                write(0,*)sc_lib_opt%residue(isc)%resname
                write(0,*)ca,   sc_lib_opt%residue(isc)%cv(1,:)
                write(0,*)ca,   sc_lib_opt%residue(isc)%cv(2,:)
                write(0,*) sqrt(min_sc(ca)/dble(sc_lib%residue(isc)%freq)),mini_sc(ca)
	            write(0,*)
         end do
     end do

     
     do i=1,dom_lib_opt%ndo
         write(0,*) all_stats_b(:,i),'domain instability',i,mini_do(i)
         write(0,*) dom_lib_opt%domain(i)%cv(1,:)
         write(0,*) dom_lib_opt%domain(i)%cv(2,:)
         write(0,*) 
     end do

     
     
     
   
end subroutine






subroutine test_window_loader
implicit none
integer ::i,ntstep,nt
type(time_window) :: window
type(window_map)  :: map
type(raw_data_block) :: block
logical,allocatable  :: window_mask(:) 


window%istate =1
window%irep   =1
window%start_time_step =1
window%end_time_step   =2

ntstep=window%end_time_step-window%start_time_step+1

allocate(window_mask(ntstep))
window_mask=.true.


call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)

    nt=count(window_mask)
    write(0,*) block%istate
    
    do i=1,nt
    write(0,*) block%xyz(:,1,i)
    end do
    
window%istate =1
window%irep   =1
window%start_time_step =2
window%end_time_step   =3

ntstep=window%end_time_step-window%start_time_step+1

!allocate(window_mask(ntstep))
window_mask=.true.

!call get_window_map(window, project, map)
call print_window_map(window, map)

call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)

    nt=count(window_mask)
    write(0,*) block%istate
    
    do i=1,nt
    write(0,*) block%xyz(:,1,i),i
    end do
    
    
end subroutine test_window_loader


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






subroutine make_plan() 
implicit none
integer :: i,j,k,fdcd,natomdcd

    write(0,*)
    write(0,*) '================================'
    write(0,*) 'PROJECT SUMMARY'
    write(0,*) '================================'

    write(0,*) 'states:',project%nstate

    do i=1,project%nstate
        write(0,*)
        write(0,*) 'state:',trim(project%state(i)%name)
        write(0,*) 'pdb  :',trim(project%state(i)%pdb_ref)
        write(0,*) 'replicates:',project%state(i)%nrep

        do j=1,project%state(i)%nrep
            write(0,*) '  rep',j
            write(0,*) '     dir :',trim(project%state(i)%rep(j)%dcd_dir)
            write(0,*) '     ndcd:',project%state(i)%rep(j)%ndcd
            do k=1,project%state(i)%rep(j)%ndcd
                write(0,*) 'file:',trim(project%state(i)%rep(j)%dcd_files(k))
                call read_dcd_header(trim(project%state(i)%rep(j)%dcd_files(k)),fdcd,natomdcd)
                write(0,*) 'frames in dcd',fdcd,'atoms in dcd', natomdcd
                pause
            end do
        end do
    end do

    write(0,*)

end subroutine make_plan 



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



subroutine test_simple_vector_cv_statisics(topologies,sc_lib,dom_lib,project,stats_b,stats_s)
implicit none

type(time_window) :: window
type(project_type), intent(in)           ::project
type(sidechain_library_type), intent(in) ::sc_lib
type(domain_library_type) , intent(in)   ::dom_lib
type(topology_set_type)   , intent(in)   ::topologies

logical,allocatable  ::activ_res(:)
logical,allocatable  ::activ_dom(:)
logical,allocatable  ::window_mask(:) 
integer :: ntstep
real(8),allocatable :: stats_b(:),stats_s(:)

!call sidechain_cv_inventur(topologies,1,sc_lib)


window%istate =1
window%irep   =1
window%start_time_step =1
window%end_time_step   =1000
ntstep=window%end_time_step-window%start_time_step+1

allocate(activ_dom(dom_lib%ndo))
allocate(window_mask(ntstep))
activ_dom=.true.
window_mask=.false.
window_mask=.true.
!window_mask(ntstep-10:ntstep)=.true.
call simple_vector_cv_statisics(topologies,window,window_mask,activ_res,activ_dom,sc_lib,dom_lib,project,stats_b,stats_s)

end subroutine




subroutine simple_vector_cv_statisics(topologies,window,window_mask,activ_res&
,activ_dom,sc_lib,dom_lib,project,stats_b,stats_s)
implicit none
type(topology_set_type), intent(in)      ::topologies
type(project_type), intent(in)           ::project
type(domain_library_type) , intent(in)   ::dom_lib
type(sidechain_library_type), intent(in) ::sc_lib
type(time_window),intent(in) :: window
type(raw_data_block) :: block
type(window_map)     :: map

logical,allocatable, intent(in) :: activ_res(:),activ_dom(:)
logical,allocatable, intent(in)  :: window_mask(:) 
integer :: time_step,i,iatom,idom,ntstep
real(8) :: b1(3),b2(3),s1(3),s2(3),lb1,lb2,b12,basis(3,3),t1(3),t2(3),t3(3),t4(3),t5(3)
real(8) :: b1_res(3),b2_res(3),s1_res(3),s2_res(3),lb1_res,lb2_res,b12_res,basis_res(3,3)
real(8) ::q1(3),q2(3)
real(8),allocatable :: l1diffsq(:),l2diffsq(:),b12diffsq(:),l1s,l2s,b12s
real(8),allocatable :: l1diffsq_res(:),l2diffsq_res(:),b12diffsq_res(:),l1s_res,l2s_res,b12s_res
real(8),allocatable :: stats_b(:),stats_s(:)
integer :: ndomain,nres,ncv,ncvtype,j, ires


!call get_window_map(window, project, map)
call load_raw_data_block( &
    project, &
    window, &
    window_mask, &
    block)
    
    ndomain=dom_lib%ndo
    
allocate( l1diffsq(ndomain))
allocate( l2diffsq(ndomain))
allocate(b12diffsq(ndomain))

nres=topologies%state(window%istate)%topo%nres
allocate( l1diffsq_res(nres))
allocate( l2diffsq_res(nres))
allocate(b12diffsq_res(nres))

ntstep=count(window_mask)

 l1diffsq=0._8
 l2diffsq=0._8
b12diffsq=0._8

l1diffsq_res = 0.d0
l2diffsq_res = 0.d0
b12diffsq_res = 0.d0

do time_step=1,ntstep




    do idom=1,dom_lib%ndo
    
    b1=0._8
    b2=0._8
    s1=0._8
    s2=0._8
    
        do i=1,dom_lib%domain(idom)%natom_keep
        
        iatom=dom_lib%domain(idom)%map(i)+dom_lib%domain(idom)%start_atom-1
        b1=b1+block%xyz(:,iatom,time_step)*dom_lib%domain(idom)%cv(1,i)
        b2=b2+block%xyz(:,iatom,time_step)*dom_lib%domain(idom)%cv(2,i)
        
        end do
        !write(0,*) block%xyz(:,iatom,time_step),'#############'
        !write(0,*)shape(block%xyz)
    lb1=dot_product(b1,b1)
    lb2=dot_product(b2,b2)
    lb1=sqrt(lb1)
    lb2=sqrt(lb2)
    b1=b1/lb1
    b2=b2/lb2
    b12=dot_product(b1,b2)
    !write(0,*) lb1,lb2,b12
    
    l1s =1._8-lb1
    l2s =1._8-lb1
    b12s=b12
             
     l1diffsq(idom)= l1diffsq(idom)+l1s *l1s 
     l2diffsq(idom)= l2diffsq(idom)+l2s *l2s 
    b12diffsq(idom)=b12diffsq(idom)+b12s*b12s
    
    basis(1,:)=b1
    basis(2,:)=b2
    basis(3,:)=cross_product(b1,b2)
    
    
    end do
    


do ncvtype = 1, size(sc_lib%residue)

    if (.not. allocated(sc_lib%residue(ncvtype)%cv)) cycle

    do j = 1, sc_lib%residue(ncvtype)%freq

        ires = sc_lib%residue(ncvtype)%res_ids(j)

        b1_res = 0.d0
        b2_res = 0.d0

        do i = 1, sc_lib%residue(ncvtype)%natom_head

            iatom = sc_lib%residue(ncvtype)%map(i) + &
                    topologies%state(window%istate)%topo%acid_start_list(ires) - 1

            b1_res = b1_res + block%xyz(:,iatom,time_step) * &
                               sc_lib%residue(ncvtype)%cv(1,i)

            b2_res = b2_res + block%xyz(:,iatom,time_step) * &
                               sc_lib%residue(ncvtype)%cv(2,i)

         !write(0,*) topologies%state(window%istate)%topo%atom_id(iatom),sc_lib%residue(ncvtype)%map(i)
        
         
       
         
         end do
         
         

!write(0,*) dot_product(b1_res,b2_res)       
!write(0,*) b1_res, b2_res,sc_lib%residue(ncvtype)%resname
        lb1_res = sqrt(dot_product(b1_res,b1_res))
        lb2_res = sqrt(dot_product(b2_res,b2_res))
!write(0,*)sc_lib%residue(ncvtype)%cv(1,:)
!write(0,*)sc_lib%residue(ncvtype)%cv(2,:)
!pause
        !if (lb1_res <= tiny(1.d0) .or. lb2_res <= tiny(1.d0)) cycle

        b1_res = b1_res / lb1_res
        b2_res = b2_res / lb2_res

        b12_res = dot_product(b1_res,b2_res)
        write(0,*) lb1_res,lb2_res,b12_res,'llc'&
        ,topologies%state(window%istate)%topo%resnames(ires)&
        ,ires,sc_lib%residue(ncvtype)%resname
        
        l1s_res  = 1.d0 - lb1_res
        l2s_res  = 1.d0 - lb2_res
        b12s_res = b12_res

        l1diffsq_res(ires)  = l1diffsq_res(ires)  + l1s_res*l1s_res
        l2diffsq_res(ires)  = l2diffsq_res(ires)  + l2s_res*l2s_res
        b12diffsq_res(ires) = b12diffsq_res(ires) + b12s_res*b12s_res

        basis_res(1,:) = b1_res
        basis_res(2,:) = b2_res
        basis_res(3,:) = cross_product(b1_res,b2_res)

    end do

end do
 
!pause
end do



 l1diffsq= l1diffsq/dble(ntstep)
 l2diffsq= l2diffsq/dble(ntstep)
b12diffsq=b12diffsq/dble(ntstep)





write(0,*)
write(0,*) 'local domain basis instability scores'
write(0,*)
do idom=1,ndomain
write(0,*) dom_lib%domain(idom)%domain_name,l1diffsq(idom),l2diffsq(idom),b12diffsq(idom)
write(0,*) dom_lib%domain(idom)%natom_keep/3,-log(l1diffsq(idom)),-log(l2diffsq(idom)),-log(b12diffsq(idom))
write(0,*) idom,-log(l1diffsq(idom))-log(l2diffsq(idom))-log(b12diffsq(idom))
end do


 l1diffsq_res= l1diffsq_res/dble(ntstep)
 l2diffsq_res= l2diffsq_res/dble(ntstep)
b12diffsq_res=b12diffsq_res/dble(ntstep)


i=0
write(0,*)
write(0,*) 'local sidechain basis instability scores'
write(0,*)
do idom=2,topologies%state(window%istate)%topo%nres-1
if(topologies%state(window%istate)%topo%resnames(idom)=='GLY'.or. topologies%state(window%istate)%topo%resnames(idom)=='ALA') cycle
i=i+1
write(0,*) topologies%state(window%istate)%topo%resnames(idom)
write(0,*) l1diffsq_res(idom),l2diffsq_res(idom),b12diffsq_res(idom)
!write(0,*) -log(l1diffsq_res(idom)),-log(l2diffsq_res(idom)),-log(b12diffsq_res(idom))
!write(0,*) idom,-log(l1diffsq_res(idom))-log(l2diffsq_res(idom))-log(b12diffsq_res(idom))
!write(0,*)
end do

!write(0,*) ndomain,'???',dom_lib%ndo
if(allocated(stats_b)) deallocate(stats_b)
allocate(stats_b(dom_lib%ndo))


if(allocated(stats_s)) deallocate(stats_s)
allocate(stats_s(i))
i=0
do ncvtype = 1, size(sc_lib%residue)

    if (.not. allocated(sc_lib%residue(ncvtype)%cv)) cycle

    do j = 1, sc_lib%residue(ncvtype)%freq
        i=i+1
        ires = sc_lib%residue(ncvtype)%res_ids(j)
        
        stats_s(i)=b12diffsq_res(ires)

  end do
end do
        
        
        
stats_b=b12diffsq

!pause
end subroutine



end program main