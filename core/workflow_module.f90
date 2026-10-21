module workflow_module
use path_module

implicit none

type :: project_workflow_type
    logical :: initialized = .false.
    logical :: load_checkpoint = .false.

    character(len=1024) :: checkpoint_file = &
        'project_checkpoint.dat'
end type project_workflow_type



!type :: random_seed_workflow_type
!    logical :: initialized = .false.
!    type(random_seed_config_type) :: config
!end type random_seed_workflow_type

	
	
	
    type :: workflow_window_type
        integer :: istate = 1
        integer :: irep = 1
        integer :: start_time_step = 1
        integer :: end_time_step = 1000
    end type workflow_window_type

    type :: library_model_params_type
        ! Zustand, fuer den die Bibliotheken aufgebaut werden
        integer :: istate = 1
        ! Methode zur Referenzstruktursuche:
        ! 'rmsd' oder 'cv_orthogonality'
        character(len=32) :: reference_method = 'cv_orthogonality'
        integer :: reference_method_index = 1
        ! Parameter der Referenzstruktursuche
        integer :: n_candidates = 100
        !(stride=1, random=2)
        integer :: method_window = 1
    end type library_model_params_type

    type :: library_run_params_type
        ! Initialisierung und Pruefung
        logical :: initialize_sidechain_maps = .true.
        logical :: compute_sidechain_frequencies = .true.
        logical :: test_sidechain_maps = .true.
        logical :: initialize_domain_maps = .true.
        
        ! User-defined CV workflow steps
        logical :: initialize_user_cvs = .true.

        ! Referenzstrukturen suchen
        logical :: search_reference_structures = .true.
    end type library_run_params_type

    type :: library_io_params_type
        ! Bibliotheken laden oder neu erzeugen
        logical :: read_sidechain_library = .false.
        logical :: read_domain_library = .false.
        logical :: read_user_cv_library = .false.

        ! Bibliotheken speichern
        logical :: write_sidechain_library = .true.
        logical :: write_domain_library = .true.
        logical :: write_user_cv_library = .true.

        character(len=1024) :: sidechain_library_file = 'sidechain_lib.txt'
        character(len=1024) :: domain_library_file = 'domain_lib.txt'
        character(len=1024) :: user_cv_library_file = 'user_cv_lib.txt'
    end type library_io_params_type

    type :: library_workflow_type
        type(library_model_params_type) :: model_params
        type(library_run_params_type)   :: run_params
        type(library_io_params_type)    :: io_params
        type(workflow_window_type)      :: reference_window
        integer :: sidechain_filter_mode = 3 ! 1=All true, 2=All false, 3=By type, 4=Individual
       logical, allocatable :: sc_type_filters(:)     ! Größe: sc_lib%ncv
       logical, allocatable :: sc_residue_filters(:)  ! Größe: Gesamtanzahl der Residuen / Instanzen
       logical, allocatable :: domain_filters(:)      ! Größe: dom_lib%ndo
       logical, allocatable :: user_cv_filters(:)     ! Größe: user_lib%n_cvs
    end type library_workflow_type

!EULICIAN LIFTING
type :: euclidean_lifting_model_params_type
    ! Parameter des geometrischen Modells
    integer :: basis_index = 9
end type euclidean_lifting_model_params_type


type :: euclidean_lifting_run_params_type
    ! 1. Auswahl des Referenz-Frames
    integer :: ref_istate = 1
    integer :: ref_irep   = 1
    integer :: ref_itime  = 1

    ! 2. Auswahl des zu liftenden States und Replikats
    integer :: lift_state_mode = 1  ! 1 = All states, 2 = Specific state
    integer :: lift_istate     = 1  ! Verwendet wenn mode == 2
    integer :: lift_rep_mode   = 1  ! 1 = All reps, 2 = Specific rep
    integer :: lift_irep       = 1  ! Verwendet wenn mode == 2

    ! 3. Auswahl des Test-Frames für PBC (nur wenn test_pbc == .true.)
    integer :: pbc_istate      = 1
    integer :: pbc_irep        = 1
    integer :: pbc_itime       = 1

    ! Parameter für init_lifting_plan_all
    integer :: offset      = 0
    integer :: stride      = 10
    integer :: size_block  = 1000000

    logical :: test_pbc        = .false.
    logical :: execute_lifting = .false.
    logical :: use_pbc_correction
end type euclidean_lifting_run_params_type


type :: euclidean_lifting_io_params_type
    character(len=1024) :: centroid_path   = '.'
    character(len=256)  :: centroid_prefix = 'centroid'

    character(len=1024) :: lifted_path     = '.'
    character(len=256)  :: lifted_prefix   = 'lifted_trajectory'

    character(len=16)   :: output_policy   = 'ask'
end type euclidean_lifting_io_params_type


type :: euclidean_lifting_workflow_type
    type(euclidean_lifting_model_params_type) :: model_params
    type(euclidean_lifting_run_params_type)   :: run_params
    type(euclidean_lifting_io_params_type)    :: io_params
end type euclidean_lifting_workflow_type


type :: lifting_model_params_type
    integer :: lag = 10
    integer :: basis_index = 9
    integer :: nlevel_target = 4

    real(8) :: tau_mem = 100.0_8
    real(8) :: dt = 1.0_8
end type

type :: lifting_run_params_type
    character(len=16) :: output_policy = 'ask'
end type

type :: lifting_io_params_type
    type(file_ref_type) :: input_trajectory
    type(file_ref_type) :: output_lifted_trajectory
end type

type :: s2_lifting_workflow_type
    type(lifting_model_params_type) :: model_params
    type(lifting_run_params_type)   :: run_params
    type(lifting_io_params_type)    :: io_params
end type

    type :: odf_model_params_type
        integer :: basis_index = 9
        integer :: nlevel = 6
        integer :: max_lag = 20
        integer :: nobj = 162

        real(8) :: tau_mem = 100.0_8
        real(8) :: dt = 1.0_8
    end type odf_model_params_type

    type :: odf_run_params_type
        real(8) :: max_ram_gb = 8.0_8
        integer :: max_chunk_size = 2*10**8
    end type odf_run_params_type

    type :: odf_io_params_type
        type(file_ref_type) :: input_trajectory
        type(file_ref_type) :: output_odf

        character(len=16) :: output_policy = 'ask'
    end type odf_io_params_type

    type :: odf_workflow_type
        type(odf_model_params_type) :: model_params
        type(odf_run_params_type)   :: run_params
        type(odf_io_params_type)    :: io_params
    end type odf_workflow_type
    
        type :: histogram_model_params_type
        integer :: nlevel_target = 4
        real(8) :: tau_mem = 100.0_8
        real(8) :: dt = 1.0_8
    end type histogram_model_params_type

    type :: histogram_run_params_type
        real(8) :: max_ram_gb = 8.0_8
        integer :: max_chunk_size = 2*10**8
    end type histogram_run_params_type

    type :: histogram_io_params_type
        type(file_ref_type) :: input_trajectory
        type(file_ref_type) :: output_histogram

        character(len=16) :: output_policy = 'ask'
    end type histogram_io_params_type

    type :: histogram_workflow_type
        type(histogram_model_params_type) :: model_params
        type(histogram_run_params_type)   :: run_params
        type(histogram_io_params_type)    :: io_params
    end type histogram_workflow_type
    
        type :: te_model_params_type
        integer :: lag = 10
        integer :: basis_index = 9
        integer :: nlevel_target = 4

        real(8) :: tau_mem = 100.0_8
        real(8) :: dt = 1.0_8

        logical :: block_shuffle = .false.
        integer :: block_size_shuffle = 10000
    end type te_model_params_type

    type :: te_run_params_type
        real(8) :: max_ram_gb = 8.0_8
        integer :: max_chunk_size = 2*10**8

        integer :: n_surrogates = 5
        integer :: random_seed = 12345
    end type te_run_params_type

    type :: te_io_params_type
        type(file_ref_type) :: input_histogram
        type(file_ref_type) :: input_odf
        type(file_ref_type) :: output_raw_te
        type(file_ref_type) :: output_surrogate_te

        character(len=16) :: output_policy = 'ask'
    end type te_io_params_type

    type :: te_workflow_type
        type(te_model_params_type) :: model_params
        type(te_run_params_type)   :: run_params
        type(te_io_params_type)    :: io_params
    end type te_workflow_type

    

!CV PLAN 

    ! Konfiguration für eine einzelne Sidechain-Bibliothek
    type :: sidechain_lib_config_type
        character(len=1024) :: filename = 'sidechain_lib.txt'
        integer :: selection_mode = 1 ! 1 = Alle, 2 = Spezifisch
        logical, allocatable :: filters(:) ! Filter-Array für die CVs DIESER Bibliothek
    end type sidechain_lib_config_type

    ! Konfiguration für eine einzelne Domain-Bibliothek
    type :: domain_lib_config_type
        character(len=1024) :: filename = 'domain_lib.txt'
        integer :: selection_mode = 1
        logical, allocatable :: filters(:)
    end type domain_lib_config_type

    ! Konfiguration für eine einzelne User-Bibliothek
    type :: user_lib_config_type
        character(len=1024) :: filename = 'user_lib.txt'
        integer :: selection_mode = 1
        logical, allocatable :: filters(:)
    end type user_lib_config_type

    ! Run- und I/O-Parameter zusammengefasst im Workflow-Typ
    type :: cv_plan_workflow_type
        logical :: initialized = .false.
        
        logical :: include_sidechain_cvs = .true.
        logical :: include_domain_cvs    = .true.
        logical :: include_user_cvs      = .true.

        ! Arrays von Bibliotheken (jede hat Dateiname, Modus und eigene Filter)
        type(sidechain_lib_config_type), allocatable :: sidechain_libs(:)
        type(domain_lib_config_type),    allocatable :: domain_libs(:)
        type(user_lib_config_type),      allocatable :: user_libs(:)
    end type cv_plan_workflow_type
contains




end module