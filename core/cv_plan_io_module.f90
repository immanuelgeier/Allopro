module cv_plan_io_module
    use direction_vectors_analyse_plan_module
    use iso_fortran_env, only: output_unit, error_unit
    implicit none

    private
    public :: write_cv_plan, read_cv_plan, print_cv_plan, write_cv_names_file

contains

    ! ==========================================================
    ! CV-PLAN IN DATEI SCHREIBEN
    ! ==========================================================
    subroutine write_cv_plan(plan, filename)
        type(cv_plan), intent(in) :: plan
        character(len=*), intent(in) :: filename
        character(len=2048) :: names_filename
        integer :: unit, ios
        integer :: separator

        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Cannot open CV plan file for writing: ', trim(filename)
            error stop
        end if

        call write_cv_plan_unit(plan, unit)
        separator = scan(trim(filename), '/', back=.true.)
        if (separator > 0) then
            names_filename = filename(:separator)//'names_'//trim(filename(separator+1:))
        else
            names_filename = 'names_'//trim(filename)
        end if
        call write_cv_names_file(plan, trim(names_filename))
        close(unit, iostat=ios)
    end subroutine write_cv_plan

    ! ==========================================================
    ! CV-PLAN AUF EINHEIT (z.B. stdout) AUSGEBEN
    ! ==========================================================
    subroutine print_cv_plan(plan)
        type(cv_plan), intent(in) :: plan
        call write_cv_plan_unit(plan, output_unit)
    end subroutine print_cv_plan

    ! Interne Methode zum Schreiben der Daten
    subroutine write_cv_plan_unit(plan, unit)
        type(cv_plan), intent(in) :: plan
        integer, intent(in) :: unit
        integer :: i

        write(unit, '(A)') '# CV plan data structure'
        write(unit, '(A)') 'format_version = 1'
        write(unit, '(A)') ''

        ! 1. Terms schreiben
        if (allocated(plan%term)) then
            write(unit, '(A,I0)') 'n_terms = ', size(plan%term)
            do i = 1, size(plan%term)
                write(unit, '(A,I0,A)') 'term_block_', i, ' {'
                write(unit, '(A,I0)')   '  atom   = ', plan%term(i)%atom
                write(unit, '(A,I0)')   '  target = ', plan%term(i)%target
                write(unit, '(A,I0)')   '  comp   = ', plan%term(i)%comp
                write(unit, '(A,ES24.16)') '  coeff  = ', plan%term(i)%coeff
                if (plan%term(i)%has_c2_sym) then
                    write(unit, '(A)')   '  has_c2_sym = .true.'
                else
                    write(unit, '(A)')   '  has_c2_sym = .false.'
                end if
                write(unit, '(A)') '}'
            end do
        else
            write(unit, '(A)') 'n_terms = 0'
        end if

        write(unit, '(A)') ''

        ! 2. Objects schreiben
        if (allocated(plan%object)) then
            write(unit, '(A,I0)') 'n_objects = ', size(plan%object)
            do i = 1, size(plan%object)
                write(unit, '(A,I0,A)') 'object_block_', i, ' {'
                write(unit, '(A,I0)')   '  first_term  = ', plan%object(i)%first_term
                write(unit, '(A,I0)')   '  last_term   = ', plan%object(i)%last_term
                write(unit, '(A,I0)')   '  object_type = ', plan%object(i)%object_type
                write(unit, '(A,I0)')   '  object_id   = ', plan%object(i)%object_id
                write(unit, '(A,I0)')   '  resid       = ', plan%object(i)%resid
                write(unit, '(A,I0)')   '  domid       = ', plan%object(i)%domid
                write(unit, '(A,A)')    '  point_group = ', trim(plan%object(i)%point_group)
                write(unit, '(A,A)')    '  name_1      = ', trim(plan%object(i)%name(1))
                write(unit, '(A,A)')    '  name_2      = ', trim(plan%object(i)%name(2))
                write(unit, '(A,I0)')   '  refid       = ', plan%object(i)%refid
                write(unit, '(A)') '}'
            end do
        else
            write(unit, '(A)') 'n_objects = 0'
        end if
    end subroutine write_cv_plan_unit


    ! ==========================================================
    ! NEU: NUR CV-INDIZES UND NAMEN SPEICHERN
    ! Ungerade Indizes = name_1, Gerade Indizes = name_2
    ! ==========================================================
    subroutine write_cv_names_file(plan, filename)
        type(cv_plan), intent(in) :: plan
        character(len=*), intent(in) :: filename
        integer :: unit, ios, i, idx_1, idx_2

        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Cannot open CV names file for writing: ', trim(filename)
            error stop
        end if

        write(unit, '(A)') '# CV Indices and Names mapping'
        write(unit, '(A)') '# Format: Index    CV_Name'
        
        if (allocated(plan%object)) then
            do i = 1, size(plan%object)
                idx_1 = 2 * i - 1
                idx_2 = 2 * i

                write(unit, '(I6, 2X, A)') idx_1, trim(plan%object(i)%name(1))
                write(unit, '(I6, 2X, A)') idx_2, trim(plan%object(i)%name(2))
            end do
        end if

        close(unit, iostat=ios)
    end subroutine write_cv_names_file


    ! ==========================================================
    ! CV-PLAN AUS DATEI EINLESEN
    ! ==========================================================
    subroutine read_cv_plan(plan, filename)
        type(cv_plan), intent(out) :: plan
        character(len=*), intent(in) :: filename

        character(len=2048) :: line
        character(len=128)  :: key
        character(len=1900) :: value

        integer :: unit, ios, line_number, eq_pos, n_items, current_idx
        logical :: in_term_block, in_object_block

        ! Initialisierung
        if (allocated(plan%term)) deallocate(plan%term)
        if (allocated(plan%object)) deallocate(plan%object)

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Cannot open CV plan file for reading: ', trim(filename)
            error stop
        end if

        line_number = 0
        in_term_block = .false.
        in_object_block = .false.
        current_idx = 0

        do
            read(unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) then
                close(unit)
                error stop 'Error reading CV plan file'
            end if

            line_number = line_number + 1
            line = adjustl(line)

            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle

            ! Block-Enden erkennen
            if (trim(line) == '}') then
                in_term_block = .false.
                in_object_block = .false.
                cycle
            end if

            ! Block-Starts erkennen
            if (index(trim(line), 'term_block_') == 1) then
                read(line(12:), *, iostat=ios) current_idx
                in_term_block = .true.
                cycle
            end if

            if (index(trim(line), 'object_block_') == 1) then
                read(line(14:), *, iostat=ios) current_idx
                in_object_block = .true.
                cycle
            end if

            eq_pos = index(line, '=')
            if (eq_pos == 0) then
                write(error_unit, '(A,I0,A)') 'Error at line ', line_number, ': Missing "="'
                error stop
            end if

            key = trim(adjustl(line(:eq_pos-1)))
            value = trim(adjustl(line(eq_pos+1:)))

            ! Globale Parameter oder Block-Inhalte verarbeiten
            if (.not. in_term_block .and. .not. in_object_block) then
                select case (trim(key))
                case ('n_terms')
                    read(value, *, iostat=ios) n_items
                    if (n_items > 0) allocate(plan%term(n_items))
                case ('n_objects')
                    read(value, *, iostat=ios) n_items
                    if (n_items > 0) allocate(plan%object(n_items))
                end select
            else if (in_term_block .and. current_idx > 0) then
                select case (trim(key))
                case ('atom')
                    read(value, *, iostat=ios) plan%term(current_idx)%atom
                case ('target')
                    read(value, *, iostat=ios) plan%term(current_idx)%target
                case ('comp')
                    read(value, *, iostat=ios) plan%term(current_idx)%comp
                case ('coeff')
                    read(value, *, iostat=ios) plan%term(current_idx)%coeff
                case ('has_c2_sym')
                    read(value, *, iostat=ios) plan%term(current_idx)%has_c2_sym
                end select
            else if (in_object_block .and. current_idx > 0) then
                select case (trim(key))
                case ('first_term')
                    read(value, *, iostat=ios) plan%object(current_idx)%first_term
                case ('last_term')
                    read(value, *, iostat=ios) plan%object(current_idx)%last_term
                case ('object_type')
                    read(value, *, iostat=ios) plan%object(current_idx)%object_type
                case ('object_id')
                    read(value, *, iostat=ios) plan%object(current_idx)%object_id
                case ('resid')
                    read(value, *, iostat=ios) plan%object(current_idx)%resid
                case ('domid')
                    read(value, *, iostat=ios) plan%object(current_idx)%domid
                case ('point_group')
                    plan%object(current_idx)%point_group = trim(value)
                case ('name_1')
                    plan%object(current_idx)%name(1) = trim(value)
                case ('name_2')
                    plan%object(current_idx)%name(2) = trim(value)
                case ('refid')
                    read(value, *, iostat=ios) plan%object(current_idx)%refid
                end select
            end if

        end do

        close(unit)
    end subroutine read_cv_plan

end module cv_plan_io_module