module te_pair_statistics_module

    use iso_fortran_env, only : real32, real64, int64

    implicit none

    private

    public :: te_pair_statistics_type
    public :: calculate_te_pair_statistics
    public :: calculate_markov_spectrum
    public :: clear_te_pair_statistics


    ! =====================================================================
    ! LAPACK DGEEV
    !
    ! Berechnet die Eigenwerte einer allgemeinen reellen Matrix.
    !
    ! Wir benötigen hier nur die Eigenwerte, nicht die Eigenvektoren.
    ! =====================================================================

   interface

    subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, &
                     vr, ldvr, work, lwork, info)

        use iso_fortran_env, only : real64

        implicit none

        character(len=1), intent(in) :: jobvl
        character(len=1), intent(in) :: jobvr

        integer, intent(in) :: n
        integer, intent(in) :: lda
        integer, intent(in) :: ldvl
        integer, intent(in) :: ldvr
        integer, intent(in) :: lwork

        real(real64), intent(inout) :: a(lda,*)

        real(real64), intent(out) :: wr(*)
        real(real64), intent(out) :: wi(*)

        real(real64), intent(out) :: vl(ldvl,*)
        real(real64), intent(out) :: vr(ldvr,*)

        real(real64), intent(inout) :: work(*)

        integer, intent(out) :: info

    end subroutine dgeev

end interface



    ! =====================================================================
    ! Statistikobjekt für EIN Paar
    !
    ! Achsenkonvention:
    !
    !   hist(f,p,q)
    !
    !   f = future
    !   p = present des betrachteten Systems
    !   q = present des Quellsystems
    ! =====================================================================

    type :: te_pair_statistics_type

        ! ---------------------------------------------------------------
        ! Dimensionen des reduzierten Histogramms
        ! ---------------------------------------------------------------

        integer :: n_future  = 0
        integer :: n_present = 0
        integer :: n_source  = 0


        ! ---------------------------------------------------------------
        ! Gesamtzahl der Counts
        ! ---------------------------------------------------------------

        integer(int64) :: n_counts = 0_int64


        ! ---------------------------------------------------------------
        ! Reduziertes 3D-Histogramm
        !
        ! hist(f,p,q)
        ! ---------------------------------------------------------------

        real(real32), allocatable :: hist(:,:,:)


        ! ---------------------------------------------------------------
        ! 2D Count-Marginalen
        !
        ! counts_fp = sum_q H(f,p,q)
        ! counts_fq = sum_p H(f,p,q)
        ! counts_pq = sum_f H(f,p,q)
        ! ---------------------------------------------------------------

        real(real64), allocatable :: counts_fp(:,:)
        real(real64), allocatable :: counts_fq(:,:)
        real(real64), allocatable :: counts_pq(:,:)


        ! ---------------------------------------------------------------
        ! 1D Count-Marginalen
        ! ---------------------------------------------------------------

        real(real64), allocatable :: counts_f(:)
        real(real64), allocatable :: counts_p(:)
        real(real64), allocatable :: counts_q(:)


        ! ---------------------------------------------------------------
        ! 2D Wahrscheinlichkeiten
        ! ---------------------------------------------------------------

        real(real64), allocatable :: p_fp(:,:)
        real(real64), allocatable :: p_fq(:,:)
        real(real64), allocatable :: p_pq(:,:)


        ! ---------------------------------------------------------------
        ! 1D Wahrscheinlichkeiten
        ! ---------------------------------------------------------------

        real(real64), allocatable :: p_f(:)
        real(real64), allocatable :: p_p(:)
        real(real64), allocatable :: p_q(:)


        ! ---------------------------------------------------------------
        ! Markov-Kernel
        !
        ! K(f,p) = P(f|p)
        !
        ! Jede Spalte summiert auf 1.
        ! ---------------------------------------------------------------

        real(real64), allocatable :: markov_k(:,:)


        ! ---------------------------------------------------------------
        ! Korrelationsfaktor
        !
        ! corr(f,p)
        !   = P(f,p) / [P(f) P(p)]
        !
        ! log_corr = log(corr)
        ! ---------------------------------------------------------------

        real(real64), allocatable :: correlation_fp(:,:)
        real(real64), allocatable :: log_correlation_fp(:,:)


        ! ---------------------------------------------------------------
        ! Pointwise Mutual Information
        !
        ! identisch zu log_correlation_fp
        ! ---------------------------------------------------------------

        real(real64), allocatable :: pmi_fp(:,:)


        ! ---------------------------------------------------------------
        ! Pointwise Transfer Entropy
        !
        ! TE(f,p,q)
        !
        ! = log [
        !       P(f,p,q) P(p)
        !       ----------------
        !       P(f,p) P(p,q)
        !     ]
        !
        ! Nur für nichtverschwindende Wahrscheinlichkeiten.
        ! ---------------------------------------------------------------

        real(real64), allocatable :: pointwise_te(:,:,:)


        ! ---------------------------------------------------------------
        ! Entropien
        ! ---------------------------------------------------------------

        real(real64) :: entropy_f   = 0.0_real64
        real(real64) :: entropy_p   = 0.0_real64
        real(real64) :: entropy_q   = 0.0_real64

        real(real64) :: entropy_fp  = 0.0_real64
        real(real64) :: entropy_fq  = 0.0_real64
        real(real64) :: entropy_pq  = 0.0_real64

        real(real64) :: entropy_fpq = 0.0_real64


        ! ---------------------------------------------------------------
        ! Bedingte Entropien
        !
        ! H(f|p)
        ! H(f|q)
        ! H(p|q)
        ! H(f|p,q)
        ! ---------------------------------------------------------------

        real(real64) :: conditional_entropy_f_given_p   = 0.0_real64
        real(real64) :: conditional_entropy_f_given_q   = 0.0_real64
        real(real64) :: conditional_entropy_p_given_q   = 0.0_real64

        real(real64) :: conditional_entropy_f_given_pq  = 0.0_real64


        ! ---------------------------------------------------------------
        ! Mutual Informations
        ! ---------------------------------------------------------------

        real(real64) :: mutual_information_fp = 0.0_real64
        real(real64) :: mutual_information_fq = 0.0_real64
        real(real64) :: mutual_information_pq = 0.0_real64


        ! ---------------------------------------------------------------
        ! Transfer Entropy
        !
        ! TE(q -> f | p)
        ! ---------------------------------------------------------------

        real(real64) :: transfer_entropy = 0.0_real64


        ! ---------------------------------------------------------------
        ! Effektive Anzahl von Zuständen
        !
        ! exp(H)
        ! ---------------------------------------------------------------

        real(real64) :: effective_states_f   = 1.0_real64
        real(real64) :: effective_states_p   = 1.0_real64
        real(real64) :: effective_states_q   = 1.0_real64

        real(real64) :: effective_states_fp  = 1.0_real64
        real(real64) :: effective_states_fq  = 1.0_real64
        real(real64) :: effective_states_pq  = 1.0_real64

        real(real64) :: effective_states_fpq = 1.0_real64


        ! ---------------------------------------------------------------
        ! Spektrale Eigenschaften des Markov-Kernels
        ! ---------------------------------------------------------------

        logical :: spectrum_available = .false.

        complex(real64), allocatable :: eigenvalues(:)

        real(real64) :: spectral_gap     = 0.0_real64
        real(real64) :: relaxation_time  = 0.0_real64

        real(real64), allocatable :: stationary_distribution(:)


        ! ---------------------------------------------------------------
        ! Entropierate des Markov-Prozesses
        ! ---------------------------------------------------------------

        real(real64) :: entropy_rate = 0.0_real64


        ! ---------------------------------------------------------------
        ! Irreversibilität
        ! ---------------------------------------------------------------

        real(real64) :: detailed_balance_error = 0.0_real64

        real(real64) :: entropy_production = 0.0_real64


        ! ---------------------------------------------------------------
        ! Ursprüngliche Bin-Indizes des reduzierten Histogramms
        ! ---------------------------------------------------------------

        integer, allocatable :: red_idx_future(:)
        integer, allocatable :: red_idx_present(:)
        integer, allocatable :: red_idx_source(:)

    end type te_pair_statistics_type


contains


    ! =====================================================================
    ! Hauptfunktion
    ! =====================================================================

    subroutine calculate_te_pair_statistics(hist_3d_reduced, &
                                            red_idx_future, &
                                            red_idx_present, &
                                            red_idx_source, &
                                            stats)

        real(real32), intent(in) :: hist_3d_reduced(:,:,:)

        integer, intent(in) :: red_idx_future(:)
        integer, intent(in) :: red_idx_present(:)
        integer, intent(in) :: red_idx_source(:)

        type(te_pair_statistics_type), intent(inout) :: stats


        integer :: nf, np, nq
        integer :: f, p, q

        real(real64) :: total_counts
        real(real64) :: x


        ! =================================================================
        ! Altes Statistikobjekt löschen
        ! =================================================================

        call clear_te_pair_statistics(stats)


        ! =================================================================
        ! Dimensionen
        ! =================================================================

        nf = size(hist_3d_reduced, 1)
        np = size(hist_3d_reduced, 2)
        nq = size(hist_3d_reduced, 3)

        stats%n_future  = nf
        stats%n_present = np
        stats%n_source  = nq


        ! =================================================================
        ! Bin-Indizes übernehmen
        ! =================================================================

        allocate(stats%red_idx_future(size(red_idx_future)))
        allocate(stats%red_idx_present(size(red_idx_present)))
        allocate(stats%red_idx_source(size(red_idx_source)))

        stats%red_idx_future = red_idx_future
        stats%red_idx_present = red_idx_present
        stats%red_idx_source = red_idx_source


        ! =================================================================
        ! Histogramm übernehmen
        ! =================================================================

        allocate(stats%hist(nf,np,nq))

        stats%hist = hist_3d_reduced


        ! =================================================================
        ! 2D Counts
        ! =================================================================

        allocate(stats%counts_fp(nf,np))
        allocate(stats%counts_fq(nf,nq))
        allocate(stats%counts_pq(np,nq))

        stats%counts_fp = sum( &
            real(hist_3d_reduced,real64), dim=3)

        stats%counts_fq = sum( &
            real(hist_3d_reduced,real64), dim=2)

        stats%counts_pq = sum( &
            real(hist_3d_reduced,real64), dim=1)


        ! =================================================================
        ! 1D Counts
        ! =================================================================

        allocate(stats%counts_f(nf))
        allocate(stats%counts_p(np))
        allocate(stats%counts_q(nq))

        stats%counts_f = sum(stats%counts_fp, dim=2)
        stats%counts_p = sum(stats%counts_fp, dim=1)
        stats%counts_q = sum(stats%counts_fq, dim=1)


        ! =================================================================
        ! Gesamtzahl
        ! =================================================================

        total_counts = sum(stats%counts_fp)

        if (total_counts <= 0.0_real64) then

            stats%n_counts = 0_int64

            return

        end if

        stats%n_counts = int(nint(total_counts),int64)


        ! =================================================================
        ! 2D Wahrscheinlichkeiten
        ! =================================================================

        allocate(stats%p_fp(nf,np))
        allocate(stats%p_fq(nf,nq))
        allocate(stats%p_pq(np,nq))

        stats%p_fp = stats%counts_fp / total_counts
        stats%p_fq = stats%counts_fq / total_counts
        stats%p_pq = stats%counts_pq / total_counts


        ! =================================================================
        ! 1D Wahrscheinlichkeiten
        ! =================================================================

        allocate(stats%p_f(nf))
        allocate(stats%p_p(np))
        allocate(stats%p_q(nq))

        stats%p_f = stats%counts_f / total_counts
        stats%p_p = stats%counts_p / total_counts
        stats%p_q = stats%counts_q / total_counts


        ! =================================================================
        ! Markov-Kernel K(f|p)
        ! =================================================================

        allocate(stats%markov_k(nf,np))

        stats%markov_k = 0.0_real64

        do p = 1, np

            if (stats%counts_p(p) > 0.0_real64) then

                stats%markov_k(:,p) = &
                    stats%counts_fp(:,p) / stats%counts_p(p)

            end if

        end do


        ! =================================================================
        ! Entropien
        ! =================================================================

        stats%entropy_f = &
            entropy_1d(stats%p_f)

        stats%entropy_p = &
            entropy_1d(stats%p_p)

        stats%entropy_q = &
            entropy_1d(stats%p_q)

        stats%entropy_fp = &
            entropy_2d(stats%p_fp)

        stats%entropy_fq = &
            entropy_2d(stats%p_fq)

        stats%entropy_pq = &
            entropy_2d(stats%p_pq)


        ! =================================================================
        ! 3D Entropie
        !
        ! Kein p_3d-Array wird angelegt, um Speicher zu sparen.
        ! =================================================================

        stats%entropy_fpq = 0.0_real64

        do f = 1, nf
            do p = 1, np
                do q = 1, nq

                    if (hist_3d_reduced(f,p,q) > 0.0_real32) then

                        x = real(hist_3d_reduced(f,p,q),real64) / &
                            total_counts

                        stats%entropy_fpq = &
                            stats%entropy_fpq - x*log(x)

                    end if

                end do
            end do
        end do


        ! =================================================================
        ! Bedingte Entropien
        ! =================================================================

        stats%conditional_entropy_f_given_p = &
            stats%entropy_fp - stats%entropy_p

        stats%conditional_entropy_f_given_q = &
            stats%entropy_fq - stats%entropy_q

        stats%conditional_entropy_p_given_q = &
            stats%entropy_pq - stats%entropy_q

        stats%conditional_entropy_f_given_pq = &
            stats%entropy_fpq - stats%entropy_pq


        ! =================================================================
        ! Mutual Informations
        ! =================================================================

        stats%mutual_information_fp = &
            stats%entropy_f + stats%entropy_p - &
            stats%entropy_fp

        stats%mutual_information_fq = &
            stats%entropy_f + stats%entropy_q - &
            stats%entropy_fq

        stats%mutual_information_pq = &
            stats%entropy_p + stats%entropy_q - &
            stats%entropy_pq


        ! =================================================================
        ! Transfer Entropy
        !
        ! TE(q -> f | p)
        !
        ! = H(f|p) - H(f|p,q)
        ! =================================================================

        stats%transfer_entropy = &
            stats%conditional_entropy_f_given_p - &
            stats%conditional_entropy_f_given_pq


        ! =================================================================
        ! Korrelation / PMI
        ! =================================================================

        allocate(stats%correlation_fp(nf,np))
        allocate(stats%log_correlation_fp(nf,np))
        allocate(stats%pmi_fp(nf,np))

        stats%correlation_fp = 0.0_real64
        stats%log_correlation_fp = 0.0_real64
        stats%pmi_fp = 0.0_real64

        do f = 1, nf
            do p = 1, np

                if (stats%p_fp(f,p) > 0.0_real64 .and. &
                    stats%p_f(f) > 0.0_real64 .and. &
                    stats%p_p(p) > 0.0_real64) then

                    x = stats%p_fp(f,p) / &
                        (stats%p_f(f)*stats%p_p(p))

                    stats%correlation_fp(f,p) = x

                    stats%log_correlation_fp(f,p) = log(x)

                    stats%pmi_fp(f,p) = log(x)

                end if

            end do
        end do


        ! =================================================================
        ! Pointwise TE
        !
        ! Hier wird das Array nur angelegt, wenn es sinnvoll ist.
        !
        ! ACHTUNG:
        ! Dieses Array kann bei großen Zustandsräumen sehr groß werden.
        ! =================================================================

        allocate(stats%pointwise_te(nf,np,nq))
        stats%pointwise_te = 0.0_real64

        do f = 1, nf
            do p = 1, np
                do q = 1, nq

                    if (hist_3d_reduced(f,p,q) > 0.0_real32 .and. &
                        stats%p_fp(f,p) > 0.0_real64 .and. &
                        stats%p_pq(p,q) > 0.0_real64 .and. &
                        stats%p_p(p) > 0.0_real64) then

                        x = ( &
                            real(hist_3d_reduced(f,p,q),real64) / &
                            total_counts * &
                            stats%p_p(p) ) / &
                            (stats%p_fp(f,p)*stats%p_pq(p,q))

                        stats%pointwise_te(f,p,q) = log(x)

                    end if

                end do
            end do
        end do


        ! =================================================================
        ! Effektive Zustandszahlen
        ! =================================================================

        stats%effective_states_f = &
            exp(stats%entropy_f)

        stats%effective_states_p = &
            exp(stats%entropy_p)

        stats%effective_states_q = &
            exp(stats%entropy_q)

        stats%effective_states_fp = &
            exp(stats%entropy_fp)

        stats%effective_states_fq = &
            exp(stats%entropy_fq)

        stats%effective_states_pq = &
            exp(stats%entropy_pq)

        stats%effective_states_fpq = &
            exp(stats%entropy_fpq)


    end subroutine calculate_te_pair_statistics


    ! =====================================================================
    ! Entropie eines 1D-Vektors
    ! =====================================================================

    pure function entropy_1d(p) result(h)

        real(real64), intent(in) :: p(:)

        real(real64) :: h
        integer :: i

        h = 0.0_real64

        do i = 1, size(p)

            if (p(i) > 0.0_real64) then
                h = h - p(i)*log(p(i))
            end if

        end do

    end function entropy_1d


    ! =====================================================================
    ! Entropie einer 2D-Verteilung
    ! =====================================================================

    pure function entropy_2d(p) result(h)

        real(real64), intent(in) :: p(:,:)

        real(real64) :: h
        integer :: i, j

        h = 0.0_real64

        do j = 1, size(p,2)
            do i = 1, size(p,1)

                if (p(i,j) > 0.0_real64) then
                    h = h - p(i,j)*log(p(i,j))
                end if

            end do
        end do

    end function entropy_2d


    ! =====================================================================
    ! Spektrale Eigenschaften des Markov-Kernels
    !
    ! K(f,p) ist spaltenstochastisch:
    !
    !   sum_f K(f,p) = 1
    !
    ! Eine echte Spektralanalyse ist nur möglich, wenn
    !
    !   n_future == n_present
    !
    ! gilt.
    ! =====================================================================

subroutine calculate_markov_spectrum(stats, dt)

    use iso_fortran_env, only : real64

    implicit none

    type(te_pair_statistics_type), intent(inout) :: stats
    real(real64),                  intent(in)    :: dt

    integer :: n
    integer :: i
    integer :: j
    integer :: info
    integer :: lwork

    real(real64), allocatable :: a(:,:)
    real(real64), allocatable :: wr(:)
    real(real64), allocatable :: wi(:)
    real(real64), allocatable :: vl(:,:)
    real(real64), allocatable :: vr(:,:)
    real(real64), allocatable :: work(:)

    real(real64) :: lambda_abs
    real(real64) :: lambda_second
    real(real64) :: x
    real(real64) :: y


    !=====================================================================
    ! Initialisierung
    !=====================================================================

    stats%spectrum_available = .false.

    stats%spectral_gap           = 0.0_real64
    stats%relaxation_time        = 0.0_real64
    stats%entropy_rate           = 0.0_real64
    stats%detailed_balance_error = 0.0_real64
    stats%entropy_production     = 0.0_real64


    !=====================================================================
    ! Nur quadratische Übergangsmatrix
    !=====================================================================

    if (size(stats%markov_k,1) /= &
        size(stats%markov_k,2)) then

        return

    end if


    n = size(stats%markov_k,1)

    if (n <= 0) return


    !=====================================================================
    ! DGEEV benötigt eine Kopie der Matrix,
    ! da A überschrieben wird.
    !=====================================================================

    allocate(a(n,n))

    a = stats%markov_k


    allocate(wr(n))
    allocate(wi(n))

    ! Bei JOBVL='N' und JOBVR='N' werden diese Arrays nicht benutzt.
    ! LDVL/LDVR müssen trotzdem >= 1 sein.
    allocate(vl(1,1))
    allocate(vr(1,1))


    !=====================================================================
    ! Workspace Query
    !=====================================================================

    lwork = -1

    allocate(work(1))

    call dgeev( &
        'N', 'N', &
        n, &
        a, n, &
        wr, wi, &
        vl, 1, &
        vr, 1, &
        work, lwork, &
        info )


    if (info /= 0) then

        deallocate(a, wr, wi, vl, vr, work)

        return

    end if


    lwork = max(1, int(work(1)))

    deallocate(work)

    allocate(work(lwork))


    !=====================================================================
    ! Eigentliche Eigenwertberechnung
    !=====================================================================

    a = stats%markov_k

    call dgeev( &
        'N', 'N', &
        n, &
        a, n, &
        wr, wi, &
        vl, 1, &
        vr, 1, &
        work, lwork, &
        info )


    if (info /= 0) then

        deallocate(a, wr, wi, vl, vr, work)

        return

    end if


    !=====================================================================
    ! Eigenwerte speichern
    !=====================================================================

    if (allocated(stats%eigenvalues)) then
        deallocate(stats%eigenvalues)
    end if

    allocate(stats%eigenvalues(n))

    do i = 1, n

        stats%eigenvalues(i) = &
            cmplx(wr(i), wi(i), kind=real64)

    end do


    !=====================================================================
    ! Zweitgrößter Eigenwert nach Betrag
    !
    ! Der Eigenwert lambda = 1 gehört zur stationären Mode und
    ! wird explizit ausgeschlossen.
    !
    ! Wichtig:
    !
    ! Wir verwenden NICHT:
    !
    !     if (abs(abs(lambda)-1) > ...)
    !
    ! weil damit beispielsweise ein Eigenwert -1 ebenfalls
    ! fälschlicherweise als stationäre Mode entfernt würde.
    !
    ! Stattdessen wird nur der Eigenwert ausgeschlossen, der
    ! tatsächlich bei +1 liegt.
    !=====================================================================

    lambda_second = 0.0_real64

    do i = 1, n

        lambda_abs = abs(stats%eigenvalues(i))

        ! Eigenwert +1 erkennen
        if (abs(real(stats%eigenvalues(i),real64) - 1.0_real64) < &
            1.0e-8_real64 .and. &
            abs(aimag(stats%eigenvalues(i))) < 1.0e-8_real64) then

            cycle

        end if

        if (lambda_abs > lambda_second) then

            lambda_second = lambda_abs

        end if

    end do


    ! Numerische Rundungsfehler abfangen
    lambda_second = min(lambda_second, 1.0_real64)


    !=====================================================================
    ! Spectral Gap
    !
    ! gap = 1 - |lambda_2|
    !=====================================================================

    stats%spectral_gap = &
        1.0_real64 - lambda_second


    !=====================================================================
    ! Relaxationszeit
    !
    ! tau = -dt / log(|lambda_2|)
    !=====================================================================

    if (lambda_second > 0.0_real64 .and. &
        lambda_second < 1.0_real64) then

        stats%relaxation_time = &
            -dt / log(lambda_second)

    else

        stats%relaxation_time = 0.0_real64

    end if


    !=====================================================================
    ! Stationäre Verteilung
    !
    ! Für K(future,present):
    !
    !     pi_future = sum_present K(future,present) pi_present
    !
    ! also
    !
    !     pi_new = K * pi
    !=====================================================================

    if (allocated(stats%stationary_distribution)) then
        deallocate(stats%stationary_distribution)
    end if

    allocate(stats%stationary_distribution(n))

    call calculate_stationary_distribution(stats)


    !=====================================================================
    ! Entropierate
    !
    ! H_rate =
    !
    !   - sum_i pi_i sum_j K(j,i) log K(j,i)
    !
    ! weil die Spalten des Kernels K(future,present)
    ! die Übergangswahrscheinlichkeiten bilden.
    !=====================================================================

    stats%entropy_rate = 0.0_real64

    do i = 1, n

        if (stats%stationary_distribution(i) <= 0.0_real64) cycle

        do j = 1, n

            if (stats%markov_k(j,i) <= 0.0_real64) cycle

            stats%entropy_rate = &
                stats%entropy_rate - &
                stats%stationary_distribution(i) * &
                stats%markov_k(j,i) * &
                log(stats%markov_k(j,i))

        end do

    end do


    !=====================================================================
    ! Detailed Balance und Entropieproduktion
    !
    ! Vorwärtsfluss:
    !
    !     J_ij = pi_i K(j,i)
    !
    ! Rückwärtsfluss:
    !
    !     J_ji = pi_j K(i,j)
    !
    ! Detailed Balance:
    !
    !     J_ij = J_ji
    !
    !=====================================================================

    stats%detailed_balance_error = 0.0_real64
    stats%entropy_production     = 0.0_real64


    do i = 1, n

        if (stats%stationary_distribution(i) < 0.0_real64) cycle

        do j = 1, n

            if (stats%stationary_distribution(j) < 0.0_real64) cycle

            !-------------------------------------------------------------
            ! Vorwärtsfluss i -> j
            !-------------------------------------------------------------

            x = &
                stats%stationary_distribution(i) * &
                stats%markov_k(j,i)


            !-------------------------------------------------------------
            ! Rückwärtsfluss j -> i
            !-------------------------------------------------------------

            y = &
                stats%stationary_distribution(j) * &
                stats%markov_k(i,j)


            !-------------------------------------------------------------
            ! Detailed-Balance-Fehler
            !-------------------------------------------------------------

            stats%detailed_balance_error = &
                stats%detailed_balance_error + &
                (x-y)**2


            !-------------------------------------------------------------
            ! Entropieproduktion
            !
            ! Nur wenn beide Flüsse positiv sind.
            !-------------------------------------------------------------

            if (x > 0.0_real64 .and. &
                y > 0.0_real64) then

                stats%entropy_production = &
                    stats%entropy_production + &
                    x * log(x/y)

            end if

        end do

    end do


    !=====================================================================
    ! Norm des Detailed-Balance-Fehlers
    !=====================================================================

    stats%detailed_balance_error = &
        sqrt(stats%detailed_balance_error)


    !=====================================================================
    ! Erfolgreich
    !=====================================================================

    stats%spectrum_available = .true.


    !=====================================================================
    ! Speicher freigeben
    !=====================================================================

    deallocate(a)
    deallocate(wr)
    deallocate(wi)
    deallocate(vl)
    deallocate(vr)
    deallocate(work)


end subroutine calculate_markov_spectrum

    ! =====================================================================
    ! Stationäre Verteilung bestimmen
    !
    ! K pi = pi
    ! =====================================================================

    subroutine calculate_stationary_distribution(stats)

        type(te_pair_statistics_type), intent(inout) :: stats

        integer :: n
        integer :: i
        integer :: iter

        real(real64), allocatable :: old_pi(:)
        real(real64), allocatable :: new_pi(:)

        real(real64) :: diff
        real(real64) :: norm


        n = size(stats%markov_k,1)

        allocate(old_pi(n))
        allocate(new_pi(n))


        old_pi = 1.0_real64 / real(n,real64)


        do iter = 1, 100000

            new_pi = matmul(stats%markov_k,old_pi)


            ! -------------------------------------------------------------
            ! Numerische Normalisierung
            ! -------------------------------------------------------------

            norm = sum(new_pi)

            if (norm > 0.0_real64) then
                new_pi = new_pi / norm
            end if


            diff = maxval(abs(new_pi-old_pi))


            old_pi = new_pi


            if (diff < 1.0e-12_real64) exit

        end do


        stats%stationary_distribution = new_pi


        deallocate(old_pi)
        deallocate(new_pi)

    end subroutine calculate_stationary_distribution


    ! =====================================================================
    ! Statistikobjekt vollständig freigeben
    ! =====================================================================

    subroutine clear_te_pair_statistics(stats)

        type(te_pair_statistics_type), intent(inout) :: stats


        if (allocated(stats%hist)) &
            deallocate(stats%hist)

        if (allocated(stats%counts_fp)) &
            deallocate(stats%counts_fp)

        if (allocated(stats%counts_fq)) &
            deallocate(stats%counts_fq)

        if (allocated(stats%counts_pq)) &
            deallocate(stats%counts_pq)

        if (allocated(stats%counts_f)) &
            deallocate(stats%counts_f)

        if (allocated(stats%counts_p)) &
            deallocate(stats%counts_p)

        if (allocated(stats%counts_q)) &
            deallocate(stats%counts_q)


        if (allocated(stats%p_fp)) &
            deallocate(stats%p_fp)

        if (allocated(stats%p_fq)) &
            deallocate(stats%p_fq)

        if (allocated(stats%p_pq)) &
            deallocate(stats%p_pq)

        if (allocated(stats%p_f)) &
            deallocate(stats%p_f)

        if (allocated(stats%p_p)) &
            deallocate(stats%p_p)

        if (allocated(stats%p_q)) &
            deallocate(stats%p_q)


        if (allocated(stats%markov_k)) &
            deallocate(stats%markov_k)

        if (allocated(stats%correlation_fp)) &
            deallocate(stats%correlation_fp)

        if (allocated(stats%log_correlation_fp)) &
            deallocate(stats%log_correlation_fp)

        if (allocated(stats%pmi_fp)) &
            deallocate(stats%pmi_fp)

        if (allocated(stats%pointwise_te)) &
            deallocate(stats%pointwise_te)


        if (allocated(stats%eigenvalues)) &
            deallocate(stats%eigenvalues)

        if (allocated(stats%stationary_distribution)) &
            deallocate(stats%stationary_distribution)


        if (allocated(stats%red_idx_future)) &
            deallocate(stats%red_idx_future)

        if (allocated(stats%red_idx_present)) &
            deallocate(stats%red_idx_present)

        if (allocated(stats%red_idx_source)) &
            deallocate(stats%red_idx_source)


        stats%n_future = 0
        stats%n_present = 0
        stats%n_source = 0

        stats%n_counts = 0_int64


        stats%entropy_f = 0.0_real64
        stats%entropy_p = 0.0_real64
        stats%entropy_q = 0.0_real64

        stats%entropy_fp = 0.0_real64
        stats%entropy_fq = 0.0_real64
        stats%entropy_pq = 0.0_real64
        stats%entropy_fpq = 0.0_real64


        stats%conditional_entropy_f_given_p = 0.0_real64
        stats%conditional_entropy_f_given_q = 0.0_real64
        stats%conditional_entropy_p_given_q = 0.0_real64
        stats%conditional_entropy_f_given_pq = 0.0_real64


        stats%mutual_information_fp = 0.0_real64
        stats%mutual_information_fq = 0.0_real64
        stats%mutual_information_pq = 0.0_real64

        stats%transfer_entropy = 0.0_real64


        stats%effective_states_f = 1.0_real64
        stats%effective_states_p = 1.0_real64
        stats%effective_states_q = 1.0_real64

        stats%effective_states_fp = 1.0_real64
        stats%effective_states_fq = 1.0_real64
        stats%effective_states_pq = 1.0_real64
        stats%effective_states_fpq = 1.0_real64


        stats%spectrum_available = .false.

        stats%spectral_gap = 0.0_real64
        stats%relaxation_time = 0.0_real64
        stats%entropy_rate = 0.0_real64

        stats%detailed_balance_error = 0.0_real64
        stats%entropy_production = 0.0_real64

    end subroutine clear_te_pair_statistics


end module te_pair_statistics_module