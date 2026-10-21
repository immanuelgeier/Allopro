module path_module

    implicit none
    private

    public :: file_ref_type

    ! Eine Datei bzw. ein Datensatz
    type :: file_ref_type
        character(len=256) :: prefix = ''
    end type file_ref_type

    ! Zentrale Verzeichnisse
    type :: path_config_type
        character(len=1024) :: project_root = ''
        character(len=1024) :: project_output = ''

        character(len=1024) :: lifted_trajectories = ''
        character(len=1024) :: histograms = ''
        character(len=1024) :: odf = ''
        character(len=1024) :: te_matrices = ''
    end type path_config_type
    
end module path_module    