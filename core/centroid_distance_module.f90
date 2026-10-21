
module centroid_distance_histogram_module

    use position_lift_module
    use lifting_plan_module
    use lifted_trajectory_io_module
    use lifting_plan_execution_module

    implicit none

contains


    !===================================================================
    !
    ! Erzeugt ein 1D-Abstandshistogramm für JEDES ungeordnete Paar.
    !
    ! Alle Paare werden gemeinsam in EINE Datei geschrieben.
    !
    ! Ausgabe:
    !
    !   pair_i  pair_j  r_center  count  probability  density
    !
    ! Das Histogramm ist ausschließlich über den skalaren Abstand r.
    !
    ! Für jedes Paar (i,j) werden alle Zeitpunkte aller Blöcke
    ! gemeinsam akkumuliert.
    !
    !===================================================================

    subroutine run_centroid_distance_histogram(lifting_plan, prefix)

        implicit none

        type(lifting_plan_type), intent(in) :: lifting_plan
        character(len=*), intent(in) :: prefix

        integer, parameter :: nbins = 1000

        real(8), parameter :: rmin = 0.0d0
        real(8), parameter :: rmax = 20.0d0

        integer :: nobject
        integer :: npairs

        integer(8), allocatable :: histogram(:,:)
        integer(8), allocatable :: nsamples(:)

        character(len=800) :: filename


        !---------------------------------------------------------------
        ! nobject aus dem ersten Trajektorienblock bestimmen
        !---------------------------------------------------------------

        call get_first_trajectory_info( &
            lifting_plan, &
            prefix, &
            nobject)


        npairs = nobject * (nobject - 1) / 2


        write(0,*) 'number of objects = ', nobject
        write(0,*) 'number of pairs   = ', npairs


        !---------------------------------------------------------------
        ! histogram(pair,bin)
        !
        ! Jedes Paar besitzt sein eigenes 1D-Histogramm.
        !---------------------------------------------------------------

        allocate(histogram(npairs,nbins))
        allocate(nsamples(npairs))

        histogram = 0_8
        nsamples  = 0_8


        !---------------------------------------------------------------
        ! Alle Blöcke einlesen und Paar-Histogramme akkumulieren
        !---------------------------------------------------------------

        call sample_centroid_distance_distribution( &
            lifting_plan, &
            prefix, &
            rmin, &
            rmax, &
            nbins, &
            nobject, &
            histogram, &
            nsamples)


        !---------------------------------------------------------------
        ! Genau EINE Datei für ALLE Paare
        !---------------------------------------------------------------

        filename = trim(prefix)//'_distance_histogram.dat'


        call write_centroid_distance_histogram( &
            filename, &
            rmin, &
            rmax, &
            histogram, &
            nsamples, &
            nobject)


        write(0,*) 'pair distance histograms written to ', &
                   trim(filename)


        deallocate(histogram)
        deallocate(nsamples)


    end subroutine run_centroid_distance_histogram



    !===================================================================
    !
    ! Liest den ersten vorhandenen Trajektorienblock und bestimmt
    ! die Anzahl der Objekte.
    !
    !===================================================================

    subroutine get_first_trajectory_info( &
        lifting_plan, &
        prefix, &
        nobject)

        implicit none

        type(lifting_plan_type), intent(in) :: lifting_plan
        character(len=*), intent(in) :: prefix

        integer, intent(out) :: nobject

        type(lifted_trajectory) :: tra

        character(len=800) :: filename

        integer :: istate
        integer :: irep
        integer :: iblock


        do istate = 1, size(lifting_plan%istates)

            do irep = 1, lifting_plan%nreps(istate)

                do iblock = 1, lifting_plan%nblock(istate,irep)

                    call make_filename_lifted_tra( &
                        istate, &
                        irep, &
                        iblock, &
                        prefix, &
                        filename)


                    call read_trajectory(filename, tra)

                    nobject = tra%nobject


                    if (allocated(tra%vec)) then
                        deallocate(tra%vec)
                    endif

                    return

                enddo

            enddo

        enddo


        write(0,*) 'ERROR: no trajectory file found'
        error stop


    end subroutine get_first_trajectory_info



    !===================================================================
    !
    ! Liest ALLE Trajektorienblöcke und akkumuliert für jedes Paar
    ! ein eigenes 1D-Histogramm.
    !
    !===================================================================

    subroutine sample_centroid_distance_distribution( &
        lifting_plan, &
        prefix, &
        rmin, &
        rmax, &
        nbins, &
        nobject, &
        histogram, &
        nsamples)

        implicit none

        type(lifting_plan_type), intent(in) :: lifting_plan

        character(len=*), intent(in) :: prefix

        real(8), intent(in) :: rmin
        real(8), intent(in) :: rmax

        integer, intent(in) :: nbins
        integer, intent(in) :: nobject

        integer(8), intent(inout) :: histogram(:,:)
        integer(8), intent(inout) :: nsamples(:)

        integer :: istate
        integer :: irep
        integer :: iblock

        character(len=800) :: filename


        do istate = 1, size(lifting_plan%istates)

            do irep = 1, lifting_plan%nreps(istate)

                do iblock = 1, lifting_plan%nblock(istate,irep)


                    call make_filename_lifted_tra( &
                        istate, &
                        irep, &
                        iblock, &
                        prefix, &
                        filename)


                    write(0,*) 'reading ', trim(filename)


                    call accumulate_centroid_distances( &
                        filename, &
                        rmin, &
                        rmax, &
                        nbins, &
                        nobject, &
                        histogram, &
                        nsamples)


                enddo

            enddo

        enddo


    end subroutine sample_centroid_distance_distribution



    !===================================================================
    !
    ! Akkumuliert die skalaren Abstände aller Paare.
    !
    ! histogram(ipair,ibin)
    !
    ! Für jedes Paar:
    !
    !       i < j
    !
    ! und für jeden Zeitpunkt:
    !
    !       r = |x_i - x_j|
    !
    !===================================================================

    subroutine accumulate_centroid_distances( &
        filename, &
        rmin, &
        rmax, &
        nbins, &
        nobject, &
        histogram, &
        nsamples)

        implicit none

        character(len=*), intent(in) :: filename

        real(8), intent(in) :: rmin
        real(8), intent(in) :: rmax

        integer, intent(in) :: nbins
        integer, intent(in) :: nobject

        integer(8), intent(inout) :: histogram(:,:)
        integer(8), intent(inout) :: nsamples(:)

        type(lifted_trajectory) :: tra

        integer :: ntimes

        integer :: itime
        integer :: i
        integer :: j

        integer :: ii
        integer :: jj

        integer :: ibin
        integer :: ipair

        real(8) :: dr
        real(8) :: r


        call read_trajectory(filename, tra)

        ntimes = tra%ntimes

        dr = (rmax-rmin)/real(nbins,8)


        do itime = 1, ntimes

            ipair = 0

            do i = 1, nobject-1

                ii = 2*i - 1

                do j = i+1, nobject

                    jj = 2*j - 1

                    ipair = ipair + 1


                    !---------------------------------------------------
                    ! Skalarer Abstand des Paares (i,j)
                    !---------------------------------------------------

                    r = sqrt(sum( &
                        (tra%vec(:,ii,itime) - &
                         tra%vec(:,jj,itime))**2))


                    !---------------------------------------------------
                    ! In das Histogramm dieses Paares einsortieren
                    !---------------------------------------------------

                    if (r >= rmin .and. r < rmax) then

                        ibin = int((r-rmin)/dr) + 1

                        if (ibin >= 1 .and. ibin <= nbins) then

                            histogram(ipair,ibin) = &
                                histogram(ipair,ibin) + 1_8

                            nsamples(ipair) = &
                                nsamples(ipair) + 1_8

                        endif

                    endif

                enddo

            enddo

        enddo


        if (allocated(tra%vec)) then
            deallocate(tra%vec)
        endif


    end subroutine accumulate_centroid_distances



    !===================================================================
    !
    ! Schreibt ALLE Paarhistogramme in EINE Datei.
    !
    ! Jede Zeile beschreibt einen Bin eines bestimmten Paares.
    !
    ! Spalten:
    !
    !   pair_i
    !   pair_j
    !   r_center
    !   count
    !   probability
    !   density
    !
    !===================================================================

    subroutine write_centroid_distance_histogram( &
        filename, &
        rmin, &
        rmax, &
        histogram, &
        nsamples, &
        nobject)

        implicit none

        character(len=*), intent(in) :: filename

        real(8), intent(in) :: rmin
        real(8), intent(in) :: rmax

        integer(8), intent(in) :: histogram(:,:)
        integer(8), intent(in) :: nsamples(:)

        integer, intent(in) :: nobject

        integer :: unit

        integer :: nbins

        integer :: i
        integer :: j
        integer :: ipair
        integer :: ibin

        real(8) :: dr
        real(8) :: rcenter
        real(8) :: probability
        real(8) :: density


        nbins = size(histogram,2)

        dr = (rmax-rmin)/real(nbins,8)


        open(newunit=unit, &
             file=filename, &
             status='replace', &
             action='write')


        write(unit,'(A)') '# scalar pair distance histograms'
        write(unit,'(A)') '# one histogram for each pair (i,j)'
        write(unit,'(A)') '#'
        write(unit,'(A)') &
            '# pair_i  pair_j  r_center  count  probability  density'


        ipair = 0

        do i = 1, nobject-1

            do j = i+1, nobject

                ipair = ipair + 1


                do ibin = 1, nbins

                    rcenter = rmin + &
                              (real(ibin,8)-0.5d0)*dr


                    if (nsamples(ipair) > 0_8) then

                        probability = &
                            real(histogram(ipair,ibin),8) / &
                            real(nsamples(ipair),8)

                        density = probability / dr

                    else

                        probability = 0.0d0
                        density     = 0.0d0

                    endif


                    write(unit,'(2(I8,1X),ES24.16,1X,I16,1X,ES24.16,1X,ES24.16)') &
                        i, &
                        j, &
                        rcenter, &
                        histogram(ipair,ibin), &
                        probability, &
                        density

                enddo


            enddo

        enddo


        close(unit)


    end subroutine write_centroid_distance_histogram


end module centroid_distance_histogram_module

