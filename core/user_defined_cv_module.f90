module user_defined_cv_module

    use topology_module, only: topology_set_type, state_topology_type, topology_type
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    private

    public :: user_cv_type
    public :: user_cv_library_type
    public :: edit_user_cv_library
    public :: read_user_cv_library
    public :: write_user_cv_library
    public :: clear_user_cv_library
    public :: create_multiple_user_cv_libraries
    public :: get_user_cv_ref,init_user_cv_maps
    public :: configure_user_cv_library
    !=============================================================
    ! Einzelner benutzerdefinierter CV-Typ
    !=============================================================
    type :: user_cv_type
        logical :: active=.true.
        character(len=64) :: cv_name = ''
        integer :: istate = 1
        character(len=32) :: reference_method = 'atom_pair' ! 'gyration' oder 'atom_pair'

        ! Atom- und Residuenauswahl für den CV
        integer :: natom = 0
        integer :: start_atom
        integer, allocatable :: atom_indices(:)

        ! Richtungsreferenzen (nur relevant, wenn reference_method == 'atom_pair')
        integer :: dir1_atom_1 = 0
        integer :: dir1_atom_2 = 0
        integer :: dir2_atom_1 = 0
        integer :: dir2_atom_2 = 0

        ! Kompatibilitätsfelder
        integer :: nres = 0
        integer, allocatable :: res_ids(:)
        character(len=3), allocatable :: res_names(:)
        character(len=4), allocatable :: atom_names(:)
        character(len=1), allocatable :: chain_ids(:)
        real(8), allocatable :: ref_xyz(:, :)
        real(8), allocatable :: cv_values(:, :)
        integer, allocatable :: map(:)
    end type user_cv_type

    !=============================================================
    ! Bibliothek für benutzerdefinierte CVs
    !=============================================================
    type :: user_cv_library_type
        integer :: n_cvs = 0
        type(user_cv_type), allocatable :: cvs(:)
    end type user_cv_library_type

contains


subroutine init_user_cv_maps(topo_state, user_lib)
    use topology_module, only: topology_type
    implicit none

    type(state_topology_type), intent(in) :: topo_state
    type(user_cv_library_type), intent(inout) :: user_lib

    integer :: icv, i, min_idx, n

    if (.not. allocated(user_lib%cvs)) return

    do icv = 1, user_lib%n_cvs
        if (.not. allocated(user_lib%cvs(icv)%atom_indices)) cycle

        n = size(user_lib%cvs(icv)%atom_indices)
        user_lib%cvs(icv)%natom = n

        if (n > 0) then
            ! Kleinsten Atomindex als Startpunkt definieren (analog zu Domain-Startatomen)
            min_idx = minval(user_lib%cvs(icv)%atom_indices)
            user_lib%cvs(icv)%start_atom = min_idx

            ! Map-Array für die relativen Atompositionen allokieren und befüllen
            if (allocated(user_lib%cvs(icv)%map)) deallocate(user_lib%cvs(icv)%map)
            allocate(user_lib%cvs(icv)%map(n))

            do i = 1, n
                user_lib%cvs(icv)%map(i) = user_lib%cvs(icv)%atom_indices(i) - min_idx + 1
            end do
        end if
    end do

end subroutine init_user_cv_maps


subroutine get_user_cv_ref(time,istate, icv, user_lib, topologies, block, xyz)
use time_window_loader
!use project_module
implicit none
    integer :: istate, irep, time, icv, i, iatom
    type(raw_data_block), intent(in) :: block
    type(user_cv_library_type), intent(inout) :: user_lib
    real(8), allocatable, intent(out) :: xyz(:,:)
    type(topology_set_type), intent(in) :: topologies

    ! Überprüfen, ob Atome für diesen User-CV definiert sind
    if (user_lib%cvs(icv)%natom < 1) then
        error stop 'User CV has no atoms defined for coordinate extraction'
    end if

    ! Speicher für die Koordinaten allokieren (3 Koordinaten x Anzahl Atome)
    allocate(xyz(3, user_lib%cvs(icv)%natom))

    ! Koordinaten über die im User-CV gespeicherten Atom-Indizes aus dem Block holen
    do i = 1, user_lib%cvs(icv)%natom
        iatom = user_lib%cvs(icv)%atom_indices(i)
        
        ! Falls die Indizes direkt globale Atomnummern im System sind:
        xyz(:, i) = block%xyz(:, iatom, time)
    end do

end subroutine get_user_cv_ref

    !=============================================================
    ! Bibliothek komplett aufräumen (Allocatables deallozieren)
    !=============================================================
    subroutine clear_user_cv_library(lib)
        type(user_cv_library_type), intent(inout) :: lib
        integer :: i

        if (allocated(lib%cvs)) then
            do i = 1, size(lib%cvs)
                if (allocated(lib%cvs(i)%atom_indices)) deallocate(lib%cvs(i)%atom_indices)
                if (allocated(lib%cvs(i)%res_ids)) deallocate(lib%cvs(i)%res_ids)
                if (allocated(lib%cvs(i)%res_names)) deallocate(lib%cvs(i)%res_names)
                if (allocated(lib%cvs(i)%atom_names)) deallocate(lib%cvs(i)%atom_names)
                if (allocated(lib%cvs(i)%chain_ids)) deallocate(lib%cvs(i)%chain_ids)
                if (allocated(lib%cvs(i)%ref_xyz)) deallocate(lib%cvs(i)%ref_xyz)
                if (allocated(lib%cvs(i)%cv_values)) deallocate(lib%cvs(i)%cv_values)
                if (allocated(lib%cvs(i)%map)) deallocate(lib%cvs(i)%map)
            end do
            deallocate(lib%cvs)
        end if
        lib%n_cvs = 0
    end subroutine clear_user_cv_library


   ! Interaktives Erstellen / Editieren einer einzelnen User-Defined CV Library
    ! =============================================================
    subroutine edit_user_cv_library(lib, topologies)
        type(user_cv_library_type), intent(inout) :: lib
        type(topology_set_type), intent(in)       :: topologies

        character(len=1024) :: input
        character(len=4096) :: line_buffer
        character(len=1024) :: filename
        integer :: ios, i, n_add, chosen_state, method_choice
        integer :: count, val, pos, a_idx

        write(*,*)
        write(*,*) '======================================'
        write(output_unit, *) '   Configure User-Defined CV Library  '
        write(*,*) '======================================'
        write(*,*)

        ! Vorherige Reste aufräumen
        call clear_user_cv_library(lib)

        ! 1. State-Index zuerst abfragen, um die passende Topologie zu wählen
        write(*,'(A)', advance='no') 'Enter state index [1]: '
        read(*,'(A)') input
        chosen_state = 1
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) chosen_state
            if (ios /= 0) error stop 'Invalid state index'
        end if

        ! 2. Anzahl der CVs abfragen (Jetzt direkt hier im Fluss!)
        write(*,'(A,I0,A)', advance='no') 'Number of custom CVs to define [', lib%n_cvs, ']: '
        read(*,'(A)') input
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) n_add
            if (ios /= 0) error stop 'Invalid number of CVs'
            lib%n_cvs = n_add
        end if

        ! Falls 0 CVs gewünscht sind, sauber abbrechen ohne Fehler
        if (lib%n_cvs <= 0) then
            write(*,*) '-> No custom CVs defined. Skipping library creation.'
            return
        end if

        if (lib%n_cvs > 0) allocate(lib%cvs(lib%n_cvs))

        ! Erweiterte Atom- und Residuenübersicht aus der gewählten Topologie ausgeben
        write(*,*)
        write(*,*) '----------------------------------------------------------------------'
        write(*,*) ' Topology atoms for State ', chosen_state, ':'
        write(*,*) '----------------------------------------------------------------------'
        write(*,*) ' Index    Atom     Residue    ResNum    Chain   '
        write(*,*) '----------------------------------------------------------------------'
        
        do a_idx = 1, topologies%state(chosen_state)%topo%natom
            write(*, '(I6, 4X, A4, 4X, A3, I8, 4X, A1)') &
                a_idx, &
                trim(topologies%state(chosen_state)%topo%atom_id(a_idx)), &
                trim(topologies%state(chosen_state)%topo%acid_id(a_idx)), &
                topologies%state(chosen_state)%topo%acid_number(a_idx), &
                topologies%state(chosen_state)%topo%strain_id(a_idx)
        end do
        write(*,*) '----------------------------------------------------------------------'

        do i = 1, lib%n_cvs
            write(*,*)
            write(*,*) '--- Configuring Custom CV #', i, ' ---'
            
            lib%cvs(i)%istate = chosen_state
            lib%cvs(i)%active = .true.

            ! Name wählen
            write(*,'(A)', advance='no') 'Enter CV Name: '
            read(*,'(A)') input
            lib%cvs(i)%cv_name = trim(adjustl(input))

            ! Atome direkt über Indizes eingeben
            write(*,*) 'Enter atom indices for this CV group (space-separated, e.g., 12 14 15 20):'
            write(*,'(A)', advance='no') 'Atom indices: '
            read(*,'(A)') line_buffer

            if (len_trim(line_buffer) > 0) then
                count = 0
                pos = 1
                do
                    read(line_buffer(pos:), *, iostat=ios) val
                    if (ios /= 0) exit
                    count = count + 1
                    pos = pos + index(line_buffer(pos:), ' ')
                    if (pos <= 1 .or. pos > len_trim(line_buffer)) exit
                end do

                if (count == 0) count = 50

                lib%cvs(i)%natom = count
                if (allocated(lib%cvs(i)%atom_indices)) deallocate(lib%cvs(i)%atom_indices)
                allocate(lib%cvs(i)%atom_indices(count))
                
                read(line_buffer, *, iostat=ios) lib%cvs(i)%atom_indices(1:count)
                if (ios /= 0) then
                    error stop 'Error parsing atom indices'
                end if
            end if

            ! Referenzrichtungsmethode wählen
            write(*,*) 'Select reference direction method:'
            write(*,*) '  1 = Gyration'
            write(*,*) '  2 = Atom pair'
            write(*,'(A)', advance='no') 'Select method [2]: '
            read(*,'(A)') input
            method_choice = 2
            if (len_trim(input) > 0) then
                read(input, *, iostat=ios) method_choice
                if (ios /= 0) error stop 'Invalid reference method option'
            end if

            select case (method_choice)
            case (1)
                lib%cvs(i)%reference_method = 'gyration'

            case (2)
                lib%cvs(i)%reference_method = 'atom_pair'

                write(*,'(A)', advance='no') '  Direction 1 - Atom 1 index [0]: '
                read(*,'(A)') input
                if (len_trim(input) > 0) read(input, *) lib%cvs(i)%dir1_atom_1

                write(*,'(A)', advance='no') '  Direction 1 - Atom 2 index [0]: '
                read(*,'(A)') input
                if (len_trim(input) > 0) read(input, *) lib%cvs(i)%dir1_atom_2

                write(*,'(A)', advance='no') '  Direction 2 - Atom 1 index [0]: '
                read(*,'(A)') input
                if (len_trim(input) > 0) read(input, *) lib%cvs(i)%dir2_atom_1

                write(*,'(A)', advance='no') '  Direction 2 - Atom 2 index [0]: '
                read(*,'(A)') input
                if (len_trim(input) > 0) read(input, *) lib%cvs(i)%dir2_atom_2

            case default
                error stop 'Invalid reference direction method choice'
            end select

        end do

        ! 3. Direkt am Ende das Speichern der Bibliothek anbieten (analog zu deinen anderen Workflows)
        write(*,*)
        write(*,'(A)', advance='no') 'Save user CV library? (y/n) [y]: '
        read(*,'(A)') input
        if (len_trim(input) == 0) input = 'y'

        if (trim(input) == 'y' .or. trim(input) == 'Y') then
            write(*,'(A,A,A)', advance='no') 'User CV library filename [user_cv_lib.txt]: '
            read(*,'(A)') filename
            if (len_trim(filename) == 0) filename = 'user_cv_lib.txt'

            call write_user_cv_library(trim(filename), lib)
            write(*,'(3A)') 'User CV library saved to "', trim(filename), '"'
        end if

        write(output_unit, *) 'User-defined CV library configuration complete.'

    end subroutine edit_user_cv_library

    !=============================================================
    ! Bibliothek in Datei schreiben
    !=============================================================
    subroutine write_user_cv_library(filename, lib)
        character(*), intent(in) :: filename
        type(user_cv_library_type), intent(in) :: lib

        integer :: i, j, u, n

        open(newunit=u, file=filename, status="replace", action="write")

        write(u, '(A,I0)') "NCVS ", lib%n_cvs

        do i = 1, lib%n_cvs

            write(u,*)
            write(u, '(A)') "USER_CV"

            write(u, '(A,1X,A)') "NAME", trim(lib%cvs(i)%cv_name)
            write(u, '(A,I0)') "ISTATE", lib%cvs(i)%istate
            write(u, '(A,1X,A)') "REFERENCE_METHOD", trim(lib%cvs(i)%reference_method)

            if (trim(lib%cvs(i)%reference_method) == 'atom_pair') then
                write(u, '(A,1X,I0,1X,I0)') "DIR1", lib%cvs(i)%dir1_atom_1, lib%cvs(i)%dir1_atom_2
                write(u, '(A,1X,I0,1X,I0)') "DIR2", lib%cvs(i)%dir2_atom_1, lib%cvs(i)%dir2_atom_2
            end if

            write(u, '(A)') "ATOM_INDICES"
            n = lib%cvs(i)%natom
            write(u,*) n
            if (n > 0 .and. allocated(lib%cvs(i)%atom_indices)) then
                write(u,*) (lib%cvs(i)%atom_indices(j), j=1, n)
            end if

            if (allocated(lib%cvs(i)%ref_xyz)) then
                write(u, '(A)') "REF_XYZ"
                write(u,*) size(lib%cvs(i)%ref_xyz, 2)
                do j = 1, size(lib%cvs(i)%ref_xyz, 2)
                    write(u, '(3F24.16)') lib%cvs(i)%ref_xyz(:, j)
                end do
            end if

            if (allocated(lib%cvs(i)%cv_values)) then
                write(u, '(A)') "CV_VALUES"
                write(u,*) size(lib%cvs(i)%cv_values, 2)
                do j = 1, size(lib%cvs(i)%cv_values, 2)
                    write(u, '(3F24.16)') lib%cvs(i)%cv_values(:, j)
                end do
            end if

            write(u, '(A)') "END_USER_CV"

        end do

        close(u)

    end subroutine write_user_cv_library


    !=============================================================
    ! Bibliothek aus Datei einlesen (inkl. automatischem Cleanup)
    !=============================================================
!=============================================================
    ! Bibliothek aus Datei einlesen (Fehlerfrei für Namen und Blöcke)
    !=============================================================
    subroutine read_user_cv_library(filename, lib)
        character(*), intent(in) :: filename
        type(user_cv_library_type), intent(inout) :: lib

        integer :: u, ios, p
        integer :: i, j, n, val
        character(256) :: line
        character(32) :: key

        ! 1. Zuerst sicherstellen, dass die Zielstruktur komplett sauber ist
        call clear_user_cv_library(lib)

        open(newunit=u, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            error stop 'Error opening user CV library file.'
        end if

        ! Anzahl CVs einlesen
        read(u, *, iostat=ios) key, lib%n_cvs
        if (ios /= 0) then
            lib%n_cvs = 0
            close(u)
            return
        end if

        if (lib%n_cvs > 0) then
            if (allocated(lib%cvs)) deallocate(lib%cvs)
            allocate(lib%cvs(lib%n_cvs))
        end if

        i = 0 ! Zähler für den aktuellen CV

        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle

            if (index(line, 'USER_CV') > 0) then
                i = i + 1
                if (i > lib%n_cvs) exit
                ! Standardwerte für diesen CV setzen
                lib%cvs(i)%active = .true.
                lib%cvs(i)%cv_name = 'unnamed_cv'
                lib%cvs(i)%istate = 1
                lib%cvs(i)%reference_method = 'atom_pair'
                lib%cvs(i)%natom = 0
                lib%cvs(i)%dir1_atom_1 = 0
                lib%cvs(i)%dir1_atom_2 = 0
                lib%cvs(i)%dir2_atom_1 = 0
                lib%cvs(i)%dir2_atom_2 = 0
                cycle
            end if

            if (i < 1 .or. i > lib%n_cvs) cycle

            if (index(line, 'NAME') > 0) then
                ! Alles nach dem Keyword NAME herausschneiden
                p = index(line, 'NAME') + 4
                lib%cvs(i)%cv_name = trim(adjustl(line(p:)))
                ! Anführungszeichen entfernen falls vorhanden
                if (len_trim(lib%cvs(i)%cv_name) >= 2) then
                    if (lib%cvs(i)%cv_name(1:1) == '"' .and. lib%cvs(i)%cv_name(len_trim(lib%cvs(i)%cv_name):len_trim(lib%cvs(i)%cv_name)) == '"') then
                        lib%cvs(i)%cv_name = lib%cvs(i)%cv_name(2:len_trim(lib%cvs(i)%cv_name)-1)
                    end if
                end if

            else if (index(line, 'ISTATE') > 0) then
                read(line, *, iostat=ios) key, lib%cvs(i)%istate

            else if (index(line, 'REFERENCE_METHOD') > 0) then
                p = index(line, 'REFERENCE_METHOD') + 16
                lib%cvs(i)%reference_method = trim(adjustl(line(p:)))

            else if (index(line, 'DIR1') > 0) then
                read(line, *, iostat=ios) key, lib%cvs(i)%dir1_atom_1, lib%cvs(i)%dir1_atom_2

            else if (index(line, 'DIR2') > 0) then
                read(line, *, iostat=ios) key, lib%cvs(i)%dir2_atom_1, lib%cvs(i)%dir2_atom_2

            else if (index(line, 'ATOM_INDICES') > 0) then
                read(u, *, iostat=ios) n
                if (ios /= 0) n = 0
                lib%cvs(i)%natom = n
                if (n > 0) then
                    if (allocated(lib%cvs(i)%atom_indices)) deallocate(lib%cvs(i)%atom_indices)
                    allocate(lib%cvs(i)%atom_indices(n))
                    lib%cvs(i)%atom_indices(:) = 0
                    do j = 1, n
                        read(u, *, iostat=ios) val
                        if (ios == 0) then
                            lib%cvs(i)%atom_indices(j) = val
                        else
                            exit
                        end if
                    end do
                end if

            else if (index(line, 'REF_XYZ') > 0) then
                read(u, *, iostat=ios) n
                if (n > 0) then
                    if (allocated(lib%cvs(i)%ref_xyz)) deallocate(lib%cvs(i)%ref_xyz)
                    allocate(lib%cvs(i)%ref_xyz(3, n))
                    lib%cvs(i)%ref_xyz(:, :) = 0.0d0
                    do j = 1, n
                        read(u, *, iostat=ios) lib%cvs(i)%ref_xyz(:, j)
                        if (ios /= 0) exit
                    end do
                end if

            else if (index(line, 'CV_VALUES') > 0) then
                read(u, *, iostat=ios) n
                if (n > 0) then
                    if (allocated(lib%cvs(i)%cv_values)) deallocate(lib%cvs(i)%cv_values)
                    allocate(lib%cvs(i)%cv_values(2, n))
                    lib%cvs(i)%cv_values(:, :) = 0.0d0
                    do j = 1, n
                        read(u, *, iostat=ios) lib%cvs(i)%cv_values(:, j)
                        if (ios /= 0) exit
                    end do
                end if
            end if

        end do

        close(u)

    end subroutine read_user_cv_library
    
    subroutine configure_user_cv_library(topologies)
    use topology_module, only: topology_set_type
    implicit none

    type(topology_set_type), intent(in) :: topologies

    type(user_cv_library_type) :: user_lib
    character(len=1024) :: input
    character(len=1024) :: filename

    write(*,*)
    write(*,*) '======================================'
    write(*,*) '   Configure User-Defined CV Library  '
    write(*,*) '======================================'
    write(*,*)

    ! 1. Aufräumen und eine einzelne Lib konfigurieren
    call clear_user_cv_library(user_lib)
    call edit_user_cv_library(user_lib, topologies)

    ! 2. Standard-Dateiname festlegen und abfragen
    filename = 'user_cv_library.cfg'

    write(*,'(A,A,A)', advance='no') 'Save this library to file [', trim(filename), ']: '
    read(*,'(A)') input
    if (len_trim(input) > 0) then
        filename = trim(input)
    end if

    ! 3. Abspeichern
    call write_user_cv_library(trim(filename), user_lib)
    write(*,'(3A)') '-> Successfully saved User CV library to "', trim(filename), '"'

end subroutine configure_user_cv_library
    ! =============================================================
    ! Routine zum interaktiven Erstellen und Abspeichern mehrerer User-CV-Bibliotheken
    ! =============================================================
    subroutine create_multiple_user_cv_libraries(topologies)
        use topology_module, only: topology_set_type
        implicit none

        type(topology_set_type), intent(in) :: topologies

        type(user_cv_library_type) :: user_lib
        character(len=1024) :: input
        character(len=1024) :: filename
        integer :: n_libs, i, ios

        write(*,*)
        write(*,*) '=================================================='
        write(*,*) '   Interactive Creator for Multiple User CV Libs  '
        write(*,*) '=================================================='
        write(*,*)

        write(*,'(A)', advance='no') 'How many separate User CV libraries do you want to create? [1]: '
        read(*,'(A)') input
        
        n_libs = 1
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) n_libs
            if (ios /= 0 .or. n_libs < 1) error stop 'Invalid number of libraries'
        end if

        do i = 1, n_libs
            write(*,*)
            write(*,*) '--------------------------------------------------'
            write(*,*) ' Configuring User CV Library #', i, ' of ', n_libs
            write(*,*) '--------------------------------------------------'

            call clear_user_cv_library(user_lib)

            call edit_user_cv_library(user_lib, topologies)

            if (i == 1) then
                filename = 'user_cv_library.cfg'
            else
                write(filename, '(A,I0,A)') 'user_cv_library_', i, '.cfg'
            end if

            write(*,'(A,A,A)', advance='no') 'Save this library to file [', trim(filename), ']: '
            read(*,'(A)') input
            if (len_trim(input) > 0) then
                filename = trim(input)
            end if

            call write_user_cv_library(trim(filename), user_lib)
            write(*,'(3A)') '-> Successfully saved User CV library to "', trim(filename), '"'
        end do

        write(*,*)
        write(*,*) '=================================================='
        write(*,*) ' All User CV libraries successfully created & saved!'
        write(*,*) '=================================================='

    end subroutine create_multiple_user_cv_libraries

end module user_defined_cv_module