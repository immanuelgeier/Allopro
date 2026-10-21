module atom_mst_module

    use iso_fortran_env, only : real64

    implicit none

    private

    public :: build_atom_mst

contains

!=========================================================
! Build a minimum spanning tree from atom coordinates.
!
! xyz      : coordinates, dimension (3,natom)
! parent   : parent(i) in the spanning tree
!            parent(root) = 0
! edge_len : distance between atom i and parent(i)
!
! The graph is complete: every atom is connected to every
! other atom with the Euclidean distance as edge weight.
!
! Prim's algorithm is used without storing the complete
! distance matrix.
!=========================================================
subroutine build_atom_mst(xyz, parent, edge_len)

    real(4), intent(in)  :: xyz(:,:)

    integer,      intent(out) :: parent(:)
    real(4), intent(out) :: edge_len(:)

    integer :: natom
    integer :: i, j
    integer :: current
    real(4) :: d2
    real(4) :: best_d2

    logical, allocatable :: in_tree(:)
    real(4), allocatable :: min_d2(:)

    natom = size(xyz,2)

    if (size(xyz,1) /= 3) then
        error stop 'build_atom_mst: xyz must have dimension (3,natom)'
    end if

    if (size(parent) /= natom) then
        error stop 'build_atom_mst: wrong size of parent'
    end if

    if (size(edge_len) /= natom) then
        error stop 'build_atom_mst: wrong size of edge_len'
    end if

    allocate(in_tree(natom))
    allocate(min_d2(natom))

    parent  = 0
    edge_len = 0.

    in_tree = .false.
    min_d2  = huge(1.)

    !-----------------------------------------------------
    ! Choose atom 1 as root.
    !-----------------------------------------------------
    current = 1
    in_tree(current) = .true.

    !-----------------------------------------------------
    ! Initialize distances from root.
    !-----------------------------------------------------
    do j = 1, natom

        if (j == current) cycle

        min_d2(j) = distance_squared(xyz(:,current), &
                                     xyz(:,j))

        parent(j) = current

    end do

    !-----------------------------------------------------
    ! Prim's algorithm
    !-----------------------------------------------------
    do i = 2, natom

        ! Find closest atom not yet in tree.
        best_d2 = huge(1.)
        current = 0

        do j = 1, natom

            if (in_tree(j)) cycle

            if (min_d2(j) < best_d2) then
                best_d2 = min_d2(j)
                current = j
            end if

        end do

        if (current == 0) then
            error stop 'build_atom_mst: failed to construct MST'
        end if

        ! Add current atom to tree.
        in_tree(current) = .true.

        edge_len(current) = sqrt(best_d2)

        !-------------------------------------------------
        ! Update distances to the new tree vertex.
        !-------------------------------------------------
        do j = 1, natom

            if (in_tree(j)) cycle

            d2 = distance_squared(xyz(:,current), &
                                  xyz(:,j))

            if (d2 < min_d2(j)) then
                min_d2(j) = d2
                parent(j) = current
            end if

        end do

    end do

    ! Root has no parent.
    parent(1) = 0
    edge_len(1) = 0.

    deallocate(in_tree)
    deallocate(min_d2)

end subroutine build_atom_mst


!=========================================================
pure function distance_squared(a,b) result(d2)

    real(4), intent(in) :: a(3)
    real(4), intent(in) :: b(3)

    real(4) :: d2
    real(4) :: dx, dy, dz

    dx = a(1) - b(1)
    dy = a(2) - b(2)
    dz = a(3) - b(3)

    d2 = dx*dx + dy*dy + dz*dz

end function distance_squared

!=========================================================

end module atom_mst_module