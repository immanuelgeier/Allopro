module mod_te_io_reconstruction

    use iso_fortran_env, only: real32, int32, int64
    use sparse_lut_module, only: hitlist_type, pair_map_type, get_lut_index
    implicit none

    private
    public :: read_pair_histogram, make_te_chunk_filename

contains

    ! =========================================================================
    ! Ein EINZELNES Paar (i,j) aus den Chunks lesen & als 3D-Array rekonstruieren
    ! =========================================================================
! =========================================================================
! Ein EINZELNES Paar (i,j) aus den Chunks lesen & als 3D-Array rekonstruieren
!
! Ausgabe:
!
!   hist_3d
!       vollständiges npix x npix x npix Histogramm
!
!   hist_3d_reduced
!       reduziertes Histogramm, bei dem alle vollständig leeren
!       Zeilen/Spalten/Ebenen der drei Dimensionen entfernt wurden
!
!   red_idx1
!       ursprüngliche Indizes der in Dimension 1 behaltenen Bins
!
!   red_idx2
!       ursprüngliche Indizes der in Dimension 2 behaltenen Bins
!
!   red_idx3
!       ursprüngliche Indizes der in Dimension 3 behaltenen Bins
!
!   nred1, nred2, nred3
!       Größen der drei reduzierten Dimensionen
! =========================================================================
subroutine read_pair_histogram(input_prefix, pair_map, hits, npix, i, j, &
                               hist_3d, hist_3d_reduced, &
                               red_idx1, red_idx2, red_idx3, &
                               nred1, nred2, nred3, found)

    character(len=*),        intent(in)  :: input_prefix
    type(pair_map_type),     intent(in)  :: pair_map
    type(hitlist_type),      intent(in)  :: hits
    integer,                 intent(in)  :: npix
    integer,                 intent(in)  :: i, j

    real(real32), allocatable, intent(out) :: hist_3d(:,:,:)
    real(real32), allocatable, intent(out) :: hist_3d_reduced(:,:,:)

    integer, allocatable, intent(out) :: red_idx1(:)
    integer, allocatable, intent(out) :: red_idx2(:)
    integer, allocatable, intent(out) :: red_idx3(:)

    integer, intent(out) :: nred1
    integer, intent(out) :: nred2
    integer, intent(out) :: nred3

    logical, intent(out) :: found

    character(len=800) :: chunk_file

    integer            :: ichunk
    integer            :: u_in
    integer            :: pair_size

    integer            :: i1, i2, i3
    integer            :: k1, k2, k3

    integer(int64)     :: read_pos
    integer(int64)     :: global_idx
    integer(int64)     :: local_idx

    real(real32), allocatable :: sparse_pair_data(:)

    logical, allocatable :: active1(:)
    logical, allocatable :: active2(:)
    logical, allocatable :: active3(:)


    ! =====================================================================
    ! Initialisierung
    ! =====================================================================

    found = .false.

    nred1 = 0
    nred2 = 0
    nred3 = 0


    ! =====================================================================
    ! 1. Prüfen, ob für das Paar Einträge vorhanden sind
    ! =====================================================================

    pair_size = pair_map%size(i, j)

    if (pair_size <= 0) then
        return
    end if

    ichunk = pair_map%chunk(i, j)

    if (ichunk <= 0) then
        return
    end if


    ! =====================================================================
    ! 2. Chunk-Dateinamen ermitteln
    ! =====================================================================

    call make_te_chunk_filename(input_prefix, ichunk, chunk_file)


    ! =====================================================================
    ! 3. Vollständiges Histogramm allokieren
    ! =====================================================================

    allocate(hist_3d(npix, npix, npix))

    hist_3d = 0.0_real32


    ! =====================================================================
    ! 4. Aktive Bins initialisieren
    !
    ! active1(k) = .true., wenn Bin k in Dimension 1 mindestens
    !              einen Nicht-Null-Eintrag besitzt.
    !
    ! Entsprechend für active2 und active3.
    ! =====================================================================

    allocate(active1(npix))
    allocate(active2(npix))
    allocate(active3(npix))

    active1 = .false.
    active2 = .false.
    active3 = .false.


    ! =====================================================================
    ! 5. Komprimierte Binärdaten genau dieses Paares lesen
    ! =====================================================================

    allocate(sparse_pair_data(pair_size))


    ! Byte-Offset in der Stream-Datei
    read_pos = (pair_map%local_start(i, j) - 1_int64) * &
               4_int64 + 1_int64


    open(newunit=u_in, file=trim(chunk_file), access='stream', &
         form='unformatted', status='old', action='read')


    read(u_in, pos=read_pos) sparse_pair_data


    close(u_in)


    ! =====================================================================
    ! 6. Sparse Flat Array -> dichtes 3D-Histogramm
    !
    ! Gleichzeitig werden die aktiven Bins für das reduzierte
    ! Histogramm bestimmt.
    ! =====================================================================

    do i1 = 0, npix - 1

        if (.not. hits%hits_lagged(i1 + 1, i)) cycle

        do i2 = 0, npix - 1

            if (.not. hits%hits_present(i2 + 1, i)) cycle

            do i3 = 0, npix - 1

                if (.not. hits%hits_present(i3 + 1, j)) cycle


                ! ---------------------------------------------------------
                ! O(1) Sparse Index ermitteln
                ! ---------------------------------------------------------

                global_idx = get_lut_index( &
                    i, j, i1, i2, i3, npix, hits, pair_map)


                if (global_idx > 0_int64) then

                    local_idx = global_idx - &
                                pair_map%start(i, j) + 1_int64


                    ! -----------------------------------------------------
                    ! Wert ins dichte Histogramm schreiben
                    ! -----------------------------------------------------

                    hist_3d(i1 + 1, i2 + 1, i3 + 1) = &
                        sparse_pair_data(local_idx)


                    ! -----------------------------------------------------
                    ! Falls dieser Eintrag tatsächlich ungleich 0 ist,
                    ! sind alle drei entsprechenden Bins aktiv.
                    !
                    ! Dadurch müssen wir später nicht noch einmal das
                    ! komplette npix^3 Histogramm durchsuchen.
                    ! -----------------------------------------------------

                    if (sparse_pair_data(local_idx) /= 0.0_real32) then

                        active1(i1 + 1) = .true.
                        active2(i2 + 1) = .true.
                        active3(i3 + 1) = .true.

                    end if

                end if

            end do
        end do
    end do


    ! Sparse Daten werden nicht mehr benötigt
    deallocate(sparse_pair_data)


    ! =====================================================================
    ! 7. Anzahl der aktiven Bins bestimmen
    ! =====================================================================

    nred1 = count(active1)
    nred2 = count(active2)
    nred3 = count(active3)


    ! =====================================================================
    ! 8. Sonderfall:
    !    Das Histogramm ist vollständig leer.
    ! =====================================================================

    if (nred1 == 0 .or. nred2 == 0 .or. nred3 == 0) then

        allocate(hist_3d_reduced(0, 0, 0))

        allocate(red_idx1(0))
        allocate(red_idx2(0))
        allocate(red_idx3(0))

        deallocate(active1)
        deallocate(active2)
        deallocate(active3)

        found = .true.

        return

    end if


    ! =====================================================================
    ! 9. Arrays mit den ursprünglichen Bin-Indizes allokieren
    ! =====================================================================

    allocate(red_idx1(nred1))
    allocate(red_idx2(nred2))
    allocate(red_idx3(nred3))


    ! =====================================================================
    ! 10. Ursprüngliche Indizes der aktiven Bins bestimmen
    ! =====================================================================

    k1 = 0

    do i1 = 1, npix

        if (active1(i1)) then

            k1 = k1 + 1
            red_idx1(k1) = i1

        end if

    end do


    k2 = 0

    do i2 = 1, npix

        if (active2(i2)) then

            k2 = k2 + 1
            red_idx2(k2) = i2

        end if

    end do


    k3 = 0

    do i3 = 1, npix

        if (active3(i3)) then

            k3 = k3 + 1
            red_idx3(k3) = i3

        end if

    end do


    ! =====================================================================
    ! 11. Reduziertes Histogramm allokieren
    ! =====================================================================

    allocate(hist_3d_reduced(nred1, nred2, nred3))


    ! =====================================================================
    ! 12. Reduziertes Histogramm aus dem vollständigen Histogramm kopieren
    ! =====================================================================

    do k1 = 1, nred1

        i1 = red_idx1(k1)

        do k2 = 1, nred2

            i2 = red_idx2(k2)

            do k3 = 1, nred3

                i3 = red_idx3(k3)

                hist_3d_reduced(k1, k2, k3) = &
                    hist_3d(i1, i2, i3)

            end do

        end do

    end do


    ! =====================================================================
    ! 13. Temporäre Arrays freigeben
    ! =====================================================================

    deallocate(active1)
    deallocate(active2)
    deallocate(active3)


    ! =====================================================================
    ! Fertig
    ! =====================================================================

    found = .true.

end subroutine read_pair_histogram



    subroutine make_te_chunk_filename(prefix, ichunk, filename)
        character(len=*), intent(in)  :: prefix
        integer,          intent(in)  :: ichunk
        character(len=*), intent(out) :: filename
        write(filename, '(A, "_te_hist_chunk_", I4.4, ".real4.bin")') trim(prefix), ichunk
    end subroutine make_te_chunk_filename


    
end module mod_te_io_reconstruction