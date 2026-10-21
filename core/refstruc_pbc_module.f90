module refstruc_pbc_module

use iso_fortran_env, only : real64
use atom_mst_module
use project_module
use time_window_loader
use time_window_file_mapping
use topology_module
use graphical_rep_module

implicit none

private

public :: refstruc_type
public :: refstruc_state_type
public :: prepare_pbc_correction
public :: init_refstruc
public :: unwrap_xyz_block
public :: test_pbc_correction
!=========================================================
! Reference structure
!=========================================================
type :: refstruc_type

    real(4), allocatable :: ref_xyz(:,:)

    integer, allocatable :: parent(:)

    real(4), allocatable :: edge_len(:)

end type refstruc_type


!=========================================================
! Reference structures for states/replicates
!=========================================================
type :: refstruc_state_type

    type(refstruc_type), allocatable :: refstruc(:)

end type refstruc_state_type


contains


subroutine test_pbc_correction( &
    project, &
    topologies, &
    frame, &
    ref)

    implicit none

    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(frame2plot), intent(in) :: frame

    ! Bereits vorher erzeugte Referenzstruktur + MST
    type(refstruc_type), intent(in) :: ref

    type(raw_data_block) :: block
    type(time_window) :: window

    logical, allocatable :: mask(:)

    real(4), allocatable :: pos_before(:,:)
    real(4), allocatable :: pos_after(:,:)

    real(4), allocatable :: occ(:)
    real(4), allocatable :: beta(:)

    logical, allocatable :: group(:)

    character(len=400) :: outfile_before
    character(len=400) :: outfile_reference
    character(len=400) :: outfile_after

    integer :: natom


    !=====================================================
    ! Set time window
    !=====================================================

    window%istate          = frame%istate
    window%irep            = frame%irep
    window%start_time_step = frame%itime
    window%end_time_step   = frame%itime

    allocate(mask(1))
    mask = .true.


    !=====================================================
    ! Load test trajectory frame
    !=====================================================

    call load_raw_data_block( &
        project, &
        window, &
        mask, &
        block)


    natom = size(block%xyz,2)


    !=====================================================
    ! Check reference structure
    !=====================================================

    if (size(ref%ref_xyz,1) /= 3) then
        error stop 'test_pbc_correction: invalid reference xyz'
    end if

    if (size(ref%ref_xyz,2) /= natom) then
        error stop &
            'test_pbc_correction: reference and test frame have different natom'
    end if

    if (size(ref%parent) /= natom) then
        error stop &
            'test_pbc_correction: reference MST has wrong size'
    end if


    !=====================================================
    ! Store original coordinates
    !=====================================================

    allocate(pos_before(3,natom))

    pos_before = block%xyz(:,:,1)


    !=====================================================
    ! Print diagnostics
    !=====================================================

    write(0,*) ' '
    write(0,*) '======================================'
    write(0,*) 'PBC TEST'
    write(0,*) 'state  = ', frame%istate
    write(0,*) 'rep    = ', frame%irep
    write(0,*) 'itime  = ', frame%itime
    write(0,*) 'natom  = ', natom
    write(0,*) 'box    = ', block%boxl(:,1)

    write(0,*) 'REFERENCE MST'
    write(0,*) 'MST max edge = ', &
        maxval(ref%edge_len)

    if (natom > 1) then

        write(0,*) 'MST mean edge = ', &
            sum(ref%edge_len(2:)) / &
            real(natom-1,8)

    end if

    write(0,*) '======================================'
    write(0,*) ' '


    !=====================================================
    ! Allocate PDB auxiliary data
    !=====================================================

    allocate(occ(natom))
    allocate(beta(natom))
    allocate(group(natom))

    occ   = 0.0
    beta  = 0.0
    group = .true.


    !=====================================================
    ! Write ORIGINAL test frame
    !=====================================================

    outfile_before = 'pbc_before.pdb'

    call write_subgroup_pdb( &
        outfile_before, &
        topologies%state(frame%istate)%topo, &
        pos_before, &
        occ, &
        beta, &
        group)


    !=====================================================
    ! Write REFERENCE structure
    !=====================================================

    outfile_reference = 'pbc_reference.pdb'

    call write_subgroup_pdb( &
        outfile_reference, &
        topologies%state(frame%istate)%topo, &
        real(ref%ref_xyz,4), &
        occ, &
        beta, &
        group)


    !=====================================================
    ! Apply PBC correction
    !
    ! IMPORTANT:
    ! ref is only read here.
    ! The reference structure is NOT modified.
    !=====================================================

    call unwrap_xyz_block( &
        ref, &
        block%xyz, &
        block%boxl)


    !=====================================================
    ! Extract corrected coordinates
    !=====================================================

    allocate(pos_after(3,natom))

    pos_after = block%xyz(:,:,1)


    !=====================================================
    ! Write CORRECTED test frame
    !=====================================================

    outfile_after = 'pbc_after.pdb'

    call write_subgroup_pdb( &
        outfile_after, &
        topologies%state(frame%istate)%topo, &
        pos_after, &
        occ, &
        beta, &
        group)


    !=====================================================
    ! Cleanup
    !=====================================================

    deallocate(mask)

    deallocate(pos_before)
    deallocate(pos_after)

    deallocate(occ)
    deallocate(beta)
    deallocate(group)

end subroutine test_pbc_correction

!=========================================================
! Prepare PBC correction
!
! Loads one reference frame and constructs the MST from
! its coordinates.
!
! istate : state
! irep   : replicate
! itime  : time step of reference frame
!
! ref    : resulting reference structure + MST
!=========================================================
subroutine prepare_pbc_correction( &
    istate, &
    irep, &
    itime, &
    project, &
    ref)

    implicit none

    integer, intent(in) :: istate
    integer, intent(in) :: irep
    integer, intent(in) :: itime

    type(project_type), intent(in) :: project

    type(refstruc_type), intent(out) :: ref

    type(raw_data_block) :: block
    type(time_window) :: window

    logical, allocatable :: mask(:)

    real, allocatable :: pos(:,:)


    !-----------------------------------------------------
    ! Define reference time window
    !-----------------------------------------------------
    window%istate           = istate
    window%irep             = irep
    window%start_time_step  = itime
    window%end_time_step    = itime


    !-----------------------------------------------------
    ! Only keep this one frame
    !-----------------------------------------------------
    allocate(mask(1))
    mask = .true.


    !-----------------------------------------------------
    ! Load trajectory frame
    !-----------------------------------------------------
    call load_raw_data_block( &
        project, &
        window, &
        mask, &
        block)


    !-----------------------------------------------------
    ! Extract coordinates of reference frame
    !-----------------------------------------------------
    allocate(pos(3,size(block%xyz,2)))

    pos = real(block%xyz(:,:,1),real64)


    !-----------------------------------------------------
    ! Construct reference structure and MST
    !-----------------------------------------------------
    call init_refstruc( &
        ref, &
        pos)


    !-----------------------------------------------------
    ! Clean up
    !-----------------------------------------------------
    deallocate(pos)
    deallocate(mask)

end subroutine prepare_pbc_correction


!=========================================================
! Initialize reference structure and calculate MST
!=========================================================
subroutine init_refstruc(ref, xyz)

    implicit none

    type(refstruc_type), intent(out) :: ref

    real(4), intent(in) :: xyz(:,:)

    integer :: natom


    natom = size(xyz,2)


    !-----------------------------------------------------
    ! Sanity check
    !-----------------------------------------------------
    if (size(xyz,1) /= 3) then

        error stop &
            'init_refstruc: xyz must have dimension (3,natom)'

    end if


    !-----------------------------------------------------
    ! Store reference coordinates in double precision
    !-----------------------------------------------------
    allocate(ref%ref_xyz(3,natom))

    ref%ref_xyz = real(xyz,real64)


    !-----------------------------------------------------
    ! Allocate MST
    !-----------------------------------------------------
    allocate(ref%parent(natom))
    allocate(ref%edge_len(natom))


    !-----------------------------------------------------
    ! Build minimum spanning tree
    !-----------------------------------------------------
    call build_atom_mst( &
        ref%ref_xyz, &
        ref%parent, &
        ref%edge_len)

end subroutine init_refstruc


!=========================================================
! PBC unwrap a trajectory block
!
! xyz(:,:,t) contains wrapped coordinates.
!
! xyz is modified in place.
!
! boxl(:,t) = (Lx,Ly,Lz)
!
! The MST calculated from the reference structure is used
! to reconstruct a continuous representation of the
! protein.
!=========================================================
subroutine unwrap_xyz_block(ref, xyz, boxl)

    implicit none

    type(refstruc_type), intent(in) :: ref

    real, intent(inout) :: xyz(:,:,:)
    real, intent(in)    :: boxl(:,:)

    integer :: natom
    integer :: nframe

    integer :: t
    integer :: i
    integer :: j
    integer :: p

    logical, allocatable :: visited(:)

    real(real64) :: Lx
    real(real64) :: Ly
    real(real64) :: Lz

    real(real64) :: dx
    real(real64) :: dy
    real(real64) :: dz

    real(real64), allocatable :: xyz_unwrapped(:,:)


    natom  = size(xyz,2)
    nframe = size(xyz,3)


    !-----------------------------------------------------
    ! Sanity checks
    !-----------------------------------------------------

    if (size(xyz,1) /= 3) then

        error stop &
            'unwrap_xyz_block: xyz must have dimension (3,natom,nframe)'

    end if


    if (size(boxl,1) /= 3) then

        error stop &
            'unwrap_xyz_block: boxl must have dimension (3,nframe)'

    end if


    if (size(boxl,2) /= nframe) then

        error stop &
            'unwrap_xyz_block: boxl and xyz have different nframe'

    end if


    if (size(ref%parent) /= natom) then

        error stop &
            'unwrap_xyz_block: reference and trajectory have different natom'

    end if


    !-----------------------------------------------------
    ! Temporary coordinates for one frame
    !-----------------------------------------------------

    allocate(xyz_unwrapped(3,natom))
    allocate(visited(natom))


    !=====================================================
    ! Loop over frames
    !=====================================================

    do t = 1, nframe


        !-------------------------------------------------
        ! Box lengths
        !-------------------------------------------------

        Lx = real(boxl(1,t),real64)
        Ly = real(boxl(2,t),real64)
        Lz = real(boxl(3,t),real64)


        if (Lx <= 0.0_real64 .or. &
            Ly <= 0.0_real64 .or. &
            Lz <= 0.0_real64) then

            error stop &
                'unwrap_xyz_block: invalid box length'

        end if


        !-------------------------------------------------
        ! Reset traversal
        !-------------------------------------------------

        visited = .false.


        !-------------------------------------------------
        ! Root atom
        !
        ! Atom 1 is the root because
        ! parent(1) = 0.
        !-------------------------------------------------

        xyz_unwrapped(:,1) = &
            real(xyz(:,1,t),real64)

        visited(1) = .true.


        !=================================================
        ! Traverse MST
        !
        ! We cannot simply use:
        !
        !   do i = 2,natom
        !
        ! because parent(i) is NOT necessarily smaller
        ! than i.
        !
        ! Instead, find an unvisited atom whose parent
        ! has already been visited.
        !=================================================

        do


            i = 0
            p = 0


            !---------------------------------------------
            ! Find next atom whose parent is already known
            !---------------------------------------------

            do j = 1, natom

                if (visited(j)) cycle

                if (ref%parent(j) <= 0) cycle

                if (visited(ref%parent(j))) then

                    i = j
                    p = ref%parent(j)

                    exit

                end if

            end do


            !---------------------------------------------
            ! No more atoms to process
            !---------------------------------------------

            if (i == 0) exit


            !---------------------------------------------
            ! Raw displacement between atom and parent
            !---------------------------------------------

            dx = real(xyz(1,i,t),real64) - &
                 xyz_unwrapped(1,p)

            dy = real(xyz(2,i,t),real64) - &
                 xyz_unwrapped(2,p)

            dz = real(xyz(3,i,t),real64) - &
                 xyz_unwrapped(3,p)


            !---------------------------------------------
            ! Minimum-image displacement
            !---------------------------------------------

            dx = dx - Lx*nint(dx/Lx)
            dy = dy - Ly*nint(dy/Ly)
            dz = dz - Lz*nint(dz/Lz)


            !---------------------------------------------
            ! Reconstruct atom position
            !---------------------------------------------

            xyz_unwrapped(1,i) = &
                xyz_unwrapped(1,p) + dx

            xyz_unwrapped(2,i) = &
                xyz_unwrapped(2,p) + dy

            xyz_unwrapped(3,i) = &
                xyz_unwrapped(3,p) + dz


            !---------------------------------------------
            ! Mark atom as processed
            !---------------------------------------------

            visited(i) = .true.


        end do


        !-------------------------------------------------
        ! Verify that the complete MST was traversed
        !-------------------------------------------------

        if (count(visited) /= natom) then

            write(0,*) 'unwrap_xyz_block: MST traversal failed'
            write(0,*) 'frame      = ', t
            write(0,*) 'visited    = ', count(visited)
            write(0,*) 'natom      = ', natom

            error stop &
                'unwrap_xyz_block: not all atoms reachable from root'

        end if


        !-------------------------------------------------
        ! Copy corrected coordinates back
        !-------------------------------------------------

        xyz(:,:,t) = real(xyz_unwrapped)


    end do


    !-----------------------------------------------------
    ! Cleanup
    !-----------------------------------------------------

    deallocate(xyz_unwrapped)
    deallocate(visited)

end subroutine unwrap_xyz_block

!=========================================================

end module refstruc_pbc_module