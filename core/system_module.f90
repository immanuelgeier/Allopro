module system_module

use project_module
use topology_module
use io_pdb_module

implicit none

contains

!=========================================================
! INIT ENVIRONMENT
!=========================================================
subroutine init_environment(project)

    type(project_type), intent(inout) :: project
    character(400) :: env

    call get_environment_variable("MDGRAPH_PROJECT", env)

    if (len_trim(env) == 0) then
        write(0,*) 'ERROR: MDGRAPH_PROJECT not set'
        stop
    end if

    project%rootdir = trim(env)

end subroutine init_environment


!=========================================================
! MAIN INITIALIZATION PIPELINE
!=========================================================
subroutine initialize_system(project, topologies)

    type(project_type),      intent(inout) :: project
    type(topology_set_type), intent(inout) :: topologies

    character(400) :: pdb_file
    character(400) :: domain_file
    integer :: istate,idom
    integer :: natom

    !-----------------------------------------------------
    ! 1. ENVIRONMENT SETUP
    !-----------------------------------------------------
    call init_environment(project)

    !-----------------------------------------------------
    ! 2. LOAD PROJECT STRUCTURE
    !-----------------------------------------------------
    call load_project_structure(project)

    write(0,*)
    write(0,*) 'project loaded'
    write(0,*) 'states:', project%nstate

    !-----------------------------------------------------
    ! 3. ALLOCATE STATE TOPOLOGIES
    !-----------------------------------------------------
    call allocate_topology_set( &
        topologies, &
        project%nstate)

    !-----------------------------------------------------
    ! 4. LOAD ALL STATE TOPOLOGIES
    !-----------------------------------------------------
    do istate = 1, project%nstate

        write(0,*)
        write(0,*) '-----------------------------------'
        write(0,*) 'loading topology for state', istate
        write(0,*) '-----------------------------------'

        pdb_file = project%state(istate)%pdb_ref
        domain_file = project%state(istate)%domains_ref
        
        topologies%state(istate)%pdb_file    = pdb_file
        topologies%state(istate)%domain_file = domain_file
        write(0,*) 'reference pdb:'
        write(0,*) trim(pdb_file)
        write(0,*) 'domaine file:'
        write(0,*) trim(domain_file)

        !-------------------------------------------------
        ! count atoms
        !-------------------------------------------------
        call count_atoms_pdb( &
            pdb_file, &
            natom)

        topologies%state(istate)%topo%natom = natom

        write(0,*) 'natom =', natom

        !-------------------------------------------------
        ! allocate topology
        !-------------------------------------------------
        call allocate_topology( &
            topologies%state(istate)%topo)

        !-------------------------------------------------
        ! read pdb
        !-------------------------------------------------
        call read_pdb( &
            pdb_file, &
            topologies%state(istate)%topo)
        
        !-------------------------------------------------
        ! read domains
        !-------------------------------------------------    
        call read_domains(domain_file,topologies%state(istate))
        !-------------------------------------------------
        ! build residue structure
        !-------------------------------------------------
        call finalize_topology( &
            topologies%state(istate)%topo)

        !-------------------------------------------------
        ! summary
        !-------------------------------------------------
        write(0,*) 'atoms    =', &
            topologies%state(istate)%topo%natom

        write(0,*) 'residues =', &
            topologies%state(istate)%topo%nres
         
        write(0,*)'ndomains', topologies%state(istate)%ndomain
        
        do idom=1,topologies%state(istate)%ndomain
             write(0,*)'idom', idom,topologies%state(istate)%domain(idom)%name
        end do
        

    end do

    write(0,*)
    write(0,*) '==================================='
    write(0,*) 'ALL TOPOLOGIES INITIALIZED'
    write(0,*) 'nstates =', topologies%nstate
    write(0,*) '==================================='
    !pause
end subroutine initialize_system


!=========================================================
! TOPOLOGY PIPELINE
!=========================================================
subroutine build_topology(topo)

    type(topology_type), intent(inout) :: topo

    call make_acid_size_list(topo)
    call make_acid_start_list(topo)
    call make_residx(topo)

end subroutine build_topology


end module system_module