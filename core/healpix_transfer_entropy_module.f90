module healpix_transfer_entropy_module
    use iso_fortran_env, only: output_unit
    use healpix_marginals_module
    use healpix_type_module
    
    implicit none
    private

    ! Öffentliche Schnittstellen
    public :: joint_entropy_storage
    public :: calculate_all_entropies
    public :: compute_te_from_storage
    public :: save_entropies_gnuplot
    public :: free_entropy_storage
    public :: calculate_te_matrix

    !! Platzhalter für deine bestehenden Strukturen (zur Kompilierbarkeit)
    !type, public :: refinement_lut
    !    integer :: sparse_size
    !end type refinement_lut

    !type, public :: healpix_marginals
    !    integer :: nobject
    !    integer :: npix
    !    integer(4), allocatable :: xy(:,:,:,:)
    !    integer(4), allocatable :: yz(:,:,:,:)
    !    ! Hinweis: Falls du x- und z-Marginale hast, nutzt die Routine unten 
    !    ! die exakt gleichen total_pairs-gewichteten Aggregationen aus xy und yz, 
    !    ! um statistische Konsistenz für das Paar (i,j) zu garantieren.
    !end type healpix_marginals

    ! Typ zur Speicherung ALLER Shannon-Entropie-Kombinationen
    type :: joint_entropy_storage
        integer :: nobject = 0
        ! Einzelentropien
        real(8), allocatable :: h_x(:,:)
        real(8), allocatable :: h_y(:,:)
        real(8), allocatable :: h_z(:,:)
        ! 2D Verbundentropien
        real(8), allocatable :: h_xy(:,:)
        real(8), allocatable :: h_yz(:,:)
        real(8), allocatable :: h_xz(:,:)
        ! 3D Verbundentropie
        real(8), allocatable :: h_xyz(:,:)
    end type joint_entropy_storage

contains
subroutine calculate_te_matrix(raw_counts, fine_to_sparse, lut, marg, te_matrix)
        integer(4), intent(in)                 :: raw_counts(:,:,:)
        integer, intent(in)                    :: fine_to_sparse(0:,0:,0:)
        type(refinement_lut), intent(in)       :: lut
        type(healpix_marginals), intent(in)    :: marg                    
        real(8), intent(out)                   :: te_matrix(:,:)         

        integer :: i, j, ix, iy, iz, idx_sparse
        real(8) :: total_pairs
        real(8) :: p_xyz, p_xy, p_yz, p_y
        real(8) :: h_xy, h_yz, h_y, h_xyz
        real(8) :: val
        real(8), parameter :: ln2 = 0.6931471805599453_8

        te_matrix = 0.0_8

        do j = 1, marg%nobject
            do i = 1, marg%nobject
                if (i == j) cycle 

                total_pairs = real(sum(marg%xy(:, :, i, j)), 8)
                if (total_pairs <= 0.0_8) cycle

                h_xy = 0.0_8; h_yz = 0.0_8; h_y  = 0.0_8

                ! H(Y)
                do iy = 0, marg%npix-1
                    val = real(sum(marg%xy(:, iy, i, j)), 8)
                    if (val > 0.0_8) then
                        p_y = val / total_pairs
                        h_y = h_y - p_y * (log(p_y) / ln2)
                    end if
                end do

                ! H(X,Y)
                do iy = 0, marg%npix-1
                    do ix = 0, marg%npix-1
                        val = real(marg%xy(ix, iy, i, j), 8)
                        if (val > 0.0_8) then
                            p_xy = val / total_pairs
                            h_xy = h_xy - p_xy * (log(p_xy) / ln2)
                        end if
                    end do
                end do

                ! H(Y,Z)
                do iz = 0, marg%npix-1
                    do iy = 0, marg%npix-1
                        val = real(marg%yz(iy, iz, i, j), 8)
                        if (val > 0.0_8) then
                            p_yz = val / total_pairs
                            h_yz = h_yz - p_yz * (log(p_yz) / ln2)
                        end if
                    end do
                end do

                ! H(X,Y,Z)
                h_xyz = 0.0_8
                do idx_sparse = 1, lut%sparse_size
                    val = real(raw_counts(idx_sparse, i, j), 8)
                    if (val <= 0.0_8) cycle
                    
                    p_xyz = val / total_pairs
                    h_xyz = h_xyz - p_xyz * (log(p_xyz) / ln2)
                end do

                ! TE = H(X,Y) + H(Y,Z) - H(Y) - H(X,Y,Z)
                te_matrix(i, j) = max(0.0_8, h_xy + h_yz - h_y - h_xyz)
            end do
        end do
    end subroutine calculate_te_matrix

    ! =========================================================================
    ! 1. Hauptroutine: Berechnet jede Kombination von Entropien für alle Paare
    ! =========================================================================
    subroutine calculate_all_entropies(raw_counts, fine_to_sparse, lut, marg, storage)
        integer(4), intent(in)                 :: raw_counts(:,:,:)
        integer, intent(in)                    :: fine_to_sparse(0:,0:,0:)
        type(refinement_lut), intent(in)       :: lut
        type(healpix_marginals), intent(in)    :: marg                    
        type(joint_entropy_storage), intent(inout) :: storage

        integer :: i, j, ix, iy, iz, idx_sparse
        real(8) :: total_pairs
        real(8) :: p_x, p_y, p_z, p_xy, p_yz, p_xz, p_xyz
        real(8) :: val
        real(8), parameter :: ln2 = 0.6931471805599453_8

        ! 1. Dynamische Speicherallokation
        storage%nobject = marg%nobject
        if (.not. allocated(storage%h_x))   allocate(storage%h_x(marg%nobject, marg%nobject))
        if (.not. allocated(storage%h_y))   allocate(storage%h_y(marg%nobject, marg%nobject))
        if (.not. allocated(storage%h_z))   allocate(storage%h_z(marg%nobject, marg%nobject))
        if (.not. allocated(storage%h_xy))  allocate(storage%h_xy(marg%nobject, marg%nobject))
        if (.not. allocated(storage%h_yz))  allocate(storage%h_yz(marg%nobject, marg%nobject))
        if (.not. allocated(storage%h_xz))  allocate(storage%h_xz(marg%nobject, marg%nobject))
        if (.not. allocated(storage%h_xyz)) allocate(storage%h_xyz(marg%nobject, marg%nobject))

        ! 2. Null-Initialisierung
        storage%h_x   = 0.0_8; storage%h_y   = 0.0_8; storage%h_z   = 0.0_8
        storage%h_xy  = 0.0_8; storage%h_yz  = 0.0_8; storage%h_xz  = 0.0_8
        storage%h_xyz = 0.0_8

        ! 3. Berechnungsschleifen
        do j = 1, marg%nobject
            do i = 1, marg%nobject
                !if (i == j) cycle 
                
                write(0,*) sum(raw_counts(:,i,j))
                write(0,*) sum(marg%xy(:,:,i,j))
                !pause

                ! Gesamtsumme der gültigen Paare für diese Kombination (i, j) bestimmen
                total_pairs = real(sum(marg%xy(:, :, i, j)), 8)
                if (total_pairs <= 0.0_8) cycle

                ! --- H(X), H(Y) und H(X,Y) via marg%xy ---
                do iy = 0, marg%npix-1
                    ! H(Y) Aggregation
                    val = real(sum(marg%xy(:, iy, i, j)), 8)
                    if (val > 0.0_8) then
                        p_y = val / total_pairs
                        storage%h_y(i, j) = storage%h_y(i, j) - p_y * (log(p_y) / ln2)
                    end if

                    do ix = 0, marg%npix-1
                        ! H(X,Y) Element
                        val = real(marg%xy(ix, iy, i, j), 8)
                        if (val > 0.0_8) then
                            p_xy = val / total_pairs
                            storage%h_xy(i, j) = storage%h_xy(i, j) - p_xy * (log(p_xy) / ln2)
                        end if
                    end do
                end do

                ! H(X) Aggregation (separat über Spalten von xy, um konsistent zu total_pairs zu bleiben)
                do ix = 0, marg%npix-1
                    val = real(sum(marg%xy(ix, :, i, j)), 8)
                    if (val > 0.0_8) then
                        p_x = val / total_pairs
                        storage%h_x(i, j) = storage%h_x(i, j) - p_x * (log(p_x) / ln2)
                    end if
                end do

                ! --- H(Z) und H(Y,Z) via marg%yz ---
                do iz = 0, marg%npix-1
                    ! H(Z) Aggregation
                    val = real(sum(marg%yz(:, iz, i, j)), 8)
                    if (val > 0.0_8) then
                        p_z = val / total_pairs
                        storage%h_z(i, j) = storage%h_z(i, j) - p_z * (log(p_z) / ln2)
                    end if

                    do iy = 0, marg%npix-1
                        ! H(Y,Z) Element
                        val = real(marg%yz(iy, iz, i, j), 8)
                        if (val > 0.0_8) then
                            p_yz = val / total_pairs
                            storage%h_yz(i, j) = storage%h_yz(i, j) - p_yz * (log(p_yz) / ln2)
                        end if
                    end do
                end do

                ! --- H(X,Z) und H(X,Y,Z) via raw_counts & LUT ---
                ! Da raw_counts die 3D Verteilung (X,Y,Z) enthält, können wir daraus
                ! sowohl H(X,Y,Z) als auch die fehlende Paar-Kombination H(X,Z) exakt berechnen.
                do idx_sparse = 1, lut%sparse_size
                    val = real(raw_counts(idx_sparse, i, j), 8)
                    if (val <= 0.0_8) cycle
                    
                    p_xyz = val / total_pairs
                    ! H(X,Y,Z)
                    storage%h_xyz(i, j) = storage%h_xyz(i, j) - p_xyz * (log(p_xyz) / ln2)
                end do

                ! H(X,Z) Rekonstruktion aus der vollen 3D-Dichte (Umrechnung über fine_to_sparse Koordinaten)
                ! Wir aggregieren Y aus der 3D-Dichte heraus, um das exakte P(X,Z) für H(X,Z) zu erhalten
                do iz = 0, marg%npix-1
                    do ix = 0, marg%npix-1
                        val = 0.0_8
                        do iy = 0, marg%npix-1
                            idx_sparse = fine_to_sparse(ix, iy, iz)
                            if (idx_sparse > 0) then
                                val = val + real(raw_counts(idx_sparse, i, j), 8)
                            end if
                        end do
                        
                        if (val > 0.0_8) then
                            p_xz = val / total_pairs
                            storage%h_xz(i, j) = storage%h_xz(i, j) - p_xz * (log(p_xz) / ln2)
                        end if
                    end do
                end do

            end do
        end do
        

    end subroutine calculate_all_entropies

    ! =========================================================================
    ! 2. Die TE-Berechnung: Reine Matrixarithmetik auf Basis des Speichertyps
    ! =========================================================================
    subroutine compute_te_from_storage(storage, te_matrix)
        type(joint_entropy_storage), intent(in) :: storage
        real(8), intent(out)                    :: te_matrix(:,:)

        ! TE = H(X,Y) + H(Y,Z) - H(Y) - H(X,Y,Z)
        te_matrix = storage%h_xy + storage%h_yz - storage%h_y - storage%h_xyz
        
        ! Numerisch bedingtes negatives Rauschen kappen
        te_matrix = max(0.0_8, te_matrix)
    end subroutine compute_te_from_storage

    ! =========================================================================
    ! 3. Speicher-Routine für Gnuplot (Alle Entropieformen + TE zeilenweise)
    ! =========================================================================
    subroutine save_entropies_gnuplot(filename, storage)
        character(len=*), intent(in)            :: filename
        type(joint_entropy_storage), intent(in) :: storage
        
        integer :: u, i, j
        real(8) :: current_te

        if (.not. allocated(storage%h_y)) then
            write(output_unit, '(A)') "Fehler: Die Entropie-Daten wurden noch nicht berechnet!"
            return
        end if

        open(newunit=u, file=filename, status='replace', action='write')

        ! Umfassender Header für Gnuplot Dokumentation
        write(u, '(A)') "# Gnuplot Comprehensive Information Theory Export"
        write(u, '(A)') "# columns: i j H(X) H(Y) H(Z) H(XY) H(YZ) H(XZ) H(XYZ) TE"

        do i = 1, storage%nobject
            do j = 1, storage%nobject
                
                if (i == j) then
                    current_te = 0.0_8
                else
                    current_te = max(0.0_8, storage%h_xy(i,j) + storage%h_yz(i,j) - storage%h_y(i,j) - storage%h_xyz(i,j))
                end if

                write(u, '(2(I4,1X),8(F12.6,1X))') &
                    i, j, &
                    storage%h_x(i, j),   storage%h_y(i, j),   storage%h_z(i, j), &
                    storage%h_xy(i, j),  storage%h_yz(i, j),  storage%h_xz(i, j), &
                    storage%h_xyz(i, j), current_te
            end do
            ! Leerzeile für korrekte Gnuplot 3D-Datenblockgrenzen (splot)
            write(u, *)
        end do

        close(u)
    end subroutine save_entropies_gnuplot

    ! =========================================================================
    ! 4. Destruktor zur sauberen Speicherfreigabe
    ! =========================================================================
    subroutine free_entropy_storage(storage)
        type(joint_entropy_storage), intent(inout) :: storage
        if (allocated(storage%h_x))   deallocate(storage%h_x)
        if (allocated(storage%h_y))   deallocate(storage%h_y)
        if (allocated(storage%h_z))   deallocate(storage%h_z)
        if (allocated(storage%h_xy))  deallocate(storage%h_xy)
        if (allocated(storage%h_yz))  deallocate(storage%h_yz)
        if (allocated(storage%h_xz))  deallocate(storage%h_xz)
        if (allocated(storage%h_xyz)) deallocate(storage%h_xyz)
        storage%nobject = 0
    end subroutine free_entropy_storage

end module healpix_transfer_entropy_module