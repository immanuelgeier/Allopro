module sparse_lut_module

    use iso_fortran_env, only : int32, int64
    use healpix_state_space_exploration
    
    implicit none
    private

    public :: build_idx_chunks
    public :: pair_map_type
    public :: chunk_map_type
    public :: allocate_pair_map
    public :: read_chunk
    public :: test_mapping
    public :: prep_lut
    public :: test_mapping_real
    public :: hitlist_type
    public :: get_lut_index
    public :: make_hitlist
    public :: compute_hit_counts
    public :: build_pair_map_odf,build_pair_map_odf_red
    
    integer, parameter :: i8 = selected_int_kind(18)


type :: pair_map_type
    integer, allocatable     :: local_start(:,:)
    integer(i8), allocatable :: start(:,:)
    integer, allocatable     :: size(:,:)
    integer, allocatable     :: chunk(:,:)
    integer, allocatable     :: active_i(:)
    integer, allocatable     :: active_j(:)
    integer :: n_active 

end type pair_map_type



type :: chunk_map_type

    integer :: n_chunks

    integer(i8), allocatable :: global_start(:)
    integer(i8), allocatable :: chunk_size(:)

end type chunk_map_type


type :: hitlist_type
logical, allocatable :: hits_lagged(:,:)
logical, allocatable :: hits_present(:,:)
integer, allocatable :: rank_lagged(:,:)
integer, allocatable :: rank_present(:,:)
integer, allocatable :: n_lagged(:)
integer, allocatable :: n_present(:)
end type

contains



function get_lut_index( &
    i, j, i1, i2, i3, &
    npix, &
    hits, &
    pair_map) result(idx)

    use iso_fortran_env, only : i8 => int64

    implicit none

    integer, intent(in)             :: i, j
    integer, intent(in)             :: i1, i2, i3
    integer, intent(in)             :: npix
    type(hitlist_type), intent(in)  :: hits
    type(pair_map_type), intent(in) :: pair_map

    integer(i8) :: idx              ! <-- WICHTIG: 64-Bit Rückgabe!
    integer     :: k1, k2, k3
    integer(i8) :: local_idx

    ! ---------------------------------------------------------
    ! 1. Rang innerhalb der Support-Listen
    ! ---------------------------------------------------------
    k1 = hits%rank_lagged(i1, i)
    k2 = hits%rank_present(i2, i)
    k3 = hits%rank_present(i3, j)

    ! ---------------------------------------------------------
    ! 2. Triplett muss im Support liegen
    ! ---------------------------------------------------------
    if (k1 <= 0 .or. k2 <= 0 .or. k3 <= 0) then
        idx = -1_i8
        return
    end if

    ! ---------------------------------------------------------
    ! 3. Lokaler Index innerhalb des Objektpaares (in 64-Bit)
    ! ---------------------------------------------------------
    local_idx = int(k1, i8) &
              + int(hits%n_lagged(i), i8) * int(k2 - 1, i8) &
              + int(hits%n_lagged(i), i8) * int(hits%n_present(i), i8) * int(k3 - 1, i8)

    ! ---------------------------------------------------------
    ! 4. Globaler LUT-Index
    ! ---------------------------------------------------------
    idx = pair_map%start(i, j) + local_idx - 1_i8

end function get_lut_index


subroutine compute_hit_counts(hits, nobj)
    implicit none
    type(hitlist_type), intent(inout) :: hits
    integer, intent(in)               :: nobj
    integer                           :: i

    if (.not. allocated(hits%n_lagged))  allocate(hits%n_lagged(nobj))
    if (.not. allocated(hits%n_present)) allocate(hits%n_present(nobj))

    do i = 1, nobj
        hits%n_lagged(i)  = count(hits%hits_lagged(:, i))
        hits%n_present(i) = count(hits%hits_present(:, i))
    end do
end subroutine compute_hit_counts

subroutine test_mapping_real(odf_lib, prefix_odf,path_odf)

    use iso_fortran_env, only : int32, int64

    implicit none

    integer, parameter :: i8 = selected_int_kind(18)

    type(odf_lib_type), intent(inout) :: odf_lib
    type(hitlist_type)                 :: hits
    type(pair_map_type)                :: pair_map

    integer :: nlevel_target
    integer :: nobj, iobj
    integer :: npix

    integer :: ri, rj
    integer :: i1, i2, i3
    integer :: idx_binary
    integer(i8) :: idx_direct
    integer :: idx_in_chunk

    integer :: block_start
    integer :: block_end

    integer :: ntest

    integer(4) :: max_chunk_size
    integer(int64) :: n_total

    integer :: n_chunks,key

    real(8) :: rn

    character(len=800) :: prefix
    character(len=800) :: filename
    character(len=*) ::prefix_odf, path_odf
    integer(int32), allocatable :: vec(:)

    real(8) :: time_start, time_end
    real(8) :: time_search_start, time_search_end
    real(8) :: total_time_binary, total_time_direct
    integer(int32) :: idx_direct_in_chunk

    !=========================================================
    ! Einstellungen
    !=========================================================
    nlevel_target = 3
    prefix = 'test_idx_real'
    max_chunk_size = 10**9
    


    !=========================================================
    ! ODF Library einlesen
    !=========================================================
    write(0,*) '----------------------------------------'
    write(0,*) 'Reading ODF library'
    write(0,*) '----------------------------------------'

    call read_odf_lib(odf_lib, prefix_odf,path_odf)


    !=========================================================
    ! Parameter
    !=========================================================
    nobj = size(odf_lib%odf(1)%hist, 2)

    npix = 2**(nlevel_target - 1)
    npix = 12 * npix * npix

    pair_map%n_active = nobj
    allocate(pair_map%active_i(nobj))
    allocate(pair_map%active_j(nobj))
    do iobj = 1, nobj
        pair_map%active_i(iobj) = iobj
        pair_map%active_j(iobj) = iobj
    end do

    write(0,*) 'nobj          = ', nobj
    write(0,*) 'nlevel source = ', odf_lib%nlevel
    write(0,*) 'nlevel target = ', nlevel_target
    write(0,*) 'npix          = ', npix
    write(0,*) 'n_odf         = ', odf_lib%n_odf


    !=========================================================
    ! Hitlists erzeugen & Counts vorberechnen
    !=========================================================
    write(0,*) '----------------------------------------'
    write(0,*) 'Building hitlists & computing counts'
    write(0,*) '----------------------------------------'

    call make_hitlist(hits, odf_lib, nlevel_target)
    
    ! NEU: Hit-Counts vorab im Typ speichern!
    call compute_hit_counts(hits, nobj)

    write(0,*) 'total present hits = ', count(hits%hits_present)
    write(0,*) 'total lagged hits  = ', count(hits%hits_lagged)


    !=========================================================
    ! Pair map erzeugen
    !=========================================================
    write(0,*) '----------------------------------------'
    write(0,*) 'Building LUT chunks'
    write(0,*) '----------------------------------------'

    call cpu_time(time_start)

    
! 1. Schnell die Pair-Map berechnen (Dauert < 0.1 Sekunden)
call build_pair_map_odf( &
        pair_map, nobj, npix, hits, &
        max_chunk_size, n_total, n_chunks )

! === SCHNELL-CHECK GEGEN NEGATIVE WERTE ===
write(0,*) 'Check Global Start (1,1) = ', pair_map%start(1,1)
write(0,*) 'Check Total LUT Keys     = ', n_total
if (pair_map%start(1,1) < 0 .or. n_total < 0) then
    write(0,*) 'FEHLER: Immer noch ein negatives Feld vorhanden!'
    stop
endif

!! 2. Chunk-Dateien erzeugen (Das dauert je nach Datengröße etwas)
!call write_idx_chunks_odf( &
!        prefix, pair_map, nobj, npix, hits, max_chunk_size )
!
   
    call cpu_time(time_end)

    write(0,*) 'time for build_idx_chunks = ', time_end - time_start
    write(0,*) 'n_total  = ', n_total
    write(0,*) 'n_chunks = ', n_chunks


    !=========================================================
    ! Zufälliges Objektpaar
    !=========================================================
    call random_number(rn)
    ri = floor(dble(nobj) * rn) + 1

    call random_number(rn)
    rj = floor(dble(nobj) * rn) + 1

    write(0,*) '----------------------------------------'
    write(0,*) 'Random object pair'
    write(0,*) '----------------------------------------'
    write(0,*) 'ri = ', ri
    write(0,*) 'rj = ', rj
    write(0,*) 'global start = ', pair_map%start(ri, rj)
    write(0,*) 'local start  = ', pair_map%local_start(ri, rj)
    write(0,*) 'size         = ', pair_map%size(ri, rj)
    write(0,*) 'chunk        = ', pair_map%chunk(ri, rj)


    !=========================================================
    ! Chunk laden
    !=========================================================
    call make_chunk_filename( &
            prefix, &
            pair_map%chunk(ri, rj), &
            filename )

    write(0,*) 'reading chunk: ', trim(filename)
    call read_chunk(filename, vec)


    !=========================================================
    ! Block-Grenzen innerhalb des Chunks
    !=========================================================
    block_start = pair_map%local_start(ri, rj)
    block_end   = pair_map%local_start(ri, rj) + pair_map%size(ri, rj) - 1

    write(0,*) 'block start = ', block_start
    write(0,*) 'block end   = ', block_end


    !=========================================================
    ! Test: Vergleiche Binärsuche vs. Direkte Berechnung
    !=========================================================
  !=========================================================
    ! Test: Vergleiche Binärsuche vs. Direkte Berechnung
    !=========================================================
    total_time_binary = 0.0d0
    total_time_direct = 0.0d0

    do ntest = 1, 10000

        ! 1. Zufällige Kombination holen
        !    - (i1, i2, i3): Die einzelnen Pixelindizes
        !    - key: Der gesuchte Pixel-Kombinations-Schlüssel (z. B. 9336)
        call get_random_pixel_kombi_from_hitlist( &
                npix, &
                ri, rj, &
                hits, &
                i1, i2, i3, &  
                key )          

        ! ----------------------------------------------------
        ! Methode A: Binärsuche im Chunk
        ! ----------------------------------------------------
        call cpu_time(time_search_start)

        idx_in_chunk = binary_search( &
                vec, &
                block_start, &
                block_end, &
                key )

        call cpu_time(time_search_end)
        total_time_binary = total_time_binary + (time_search_end - time_search_start)


! ----------------------------------------------------
        ! Methode B: Direkte O(1) Index-Berechnung (Global)
        ! ----------------------------------------------------
        call cpu_time(time_search_start)

        idx_direct = get_lut_index( &
                ri, rj, &
                i1, i2, i3, &
                npix, &
                hits, &
                pair_map )

        call cpu_time(time_search_end)
        total_time_direct = total_time_direct + (time_search_end - time_search_start)

        ! ----------------------------------------------------
        ! Umrechnung: Globaler Index -> Relativer Index im geladenen Chunk
        ! ----------------------------------------------------
        ! idx_direct ist der globale Index. Um ihn im Chunk 'vec' abzufragen,
        ! berechnen wir die relative Position innerhalb des Objektpaars:
        idx_direct_in_chunk = int(idx_direct - pair_map%start(ri, rj) &
                                  + int(pair_map%local_start(ri, rj), i8), int32)

        ! ----------------------------------------------------
        ! Korrekte Prüfungen:
        ! ----------------------------------------------------
        
        ! 1. Wurde der Schlüssel mit der Binärsuche überhaupt gefunden?
        if (idx_in_chunk < block_start .or. idx_in_chunk > block_end) then
            write(0,*) 'ERROR: Key not found via binary search'
            write(0,*) '  Key = ', key
            stop
        end if

        ! 2. Stimmt der berechnete Chunk-Index mit dem der Binärsuche überein?
        if (idx_direct_in_chunk /= idx_in_chunk) then
            write(0,*) 'ERROR: Index mismatch!'
            write(0,*) '  Key (searched)          = ', key
            write(0,*) '  Index via Binary Search = ', idx_in_chunk
            write(0,*) '  Index via Direct Calc   = ', idx_direct_in_chunk
            write(0,*) '  Global Direct Index     = ', idx_direct
            stop
        end if

        ! 3. Steht an der direkt berechneten Position im Chunk auch der Schlüssel?
        if (vec(idx_direct_in_chunk) /= key) then
            write(0,*) 'ERROR: Direct index points to wrong array content!'
            write(0,*) '  Expected Key = ', key
            write(0,*) '  Found Key    = ', vec(idx_direct_in_chunk)
            write(0,*) '  At Index     = ', idx_direct_in_chunk
            stop
        end if
    
    
    end do

    !=========================================================
    ! Ergebnis
    !=========================================================
    write(0,*) '----------------------------------------'
    write(0,*) 'Mapping test successful'
    write(0,*) '----------------------------------------'
    write(0,*) 'Total time (Binary Search) = ', total_time_binary
    write(0,*) 'Total time (Direct O(1))   = ', total_time_direct
    write(0,*) 'Speedup                    = ', total_time_binary / max(total_time_direct, 1.0d-12)
    write(0,*) '----------------------------------------'

    deallocate(vec)
    
    !call sample_te_histograms_single_pass( &
    !        lifting_plan, hits, pair_map, &
    !        npix, 1, output_prefix)

end subroutine test_mapping_real

subroutine test_mapping_real_old(odf_lib, prefix_odf,path_odf)

use iso_fortran_env, only : int32,int64

implicit none

integer, parameter :: i8 = selected_int_kind(18)

type(odf_lib_type),intent(inout) :: odf_lib
type(hitlist_type) :: hits
type(pair_map_type) :: pair_map
!character(len=*), intent(in) :: prefix_odf,path_odf
integer :: nlevel_target
integer :: nobj,iobj
integer :: npix

integer :: ri
integer :: rj
integer :: idx
integer :: idx_in_chunk

integer :: block_start
integer :: block_end

integer :: ntest

integer(4) :: max_chunk_size
integer(int64) :: n_total

integer :: n_chunks

real(8) :: rn

character(len=800) :: prefix
character(len=800) :: filename

integer(int32), allocatable :: vec(:)

real(8) :: time_start
real(8) :: time_end

real(8) :: time_search_start
real(8) :: time_search_end

real(8) :: total_time
character(len=*),intent(in) :: prefix_odf,path_odf

!=========================================================
! Einstellungen
!=========================================================

nlevel_target = 3!odf_lib%nlevel

prefix = 'test_idx_real'

max_chunk_size = 10**9


ntest = 1000


!=========================================================
! ODF Library einlesen
!=========================================================

write(0,*) '----------------------------------------'
write(0,*) 'Reading ODF library'
write(0,*) '----------------------------------------'

call read_odf_lib(odf_lib, prefix_odf,path_odf)

!=========================================================
! Parameter
!=========================================================

nobj = size(odf_lib%odf(1)%hist,2)

npix = 2**(nlevel_target-1)
npix = 12*npix*npix


pair_map%n_active=nobj
allocate(pair_map%active_i(nobj))
allocate(pair_map%active_j(nobj))


do iobj=1,nobj
pair_map%active_i(iobj)=iobj
pair_map%active_j(iobj)=iobj
end do


write(0,*) 'nobj          = ',nobj
write(0,*) 'nlevel source = ',odf_lib%nlevel
write(0,*) 'nlevel target = ',nlevel_target
write(0,*) 'npix          = ',npix
write(0,*) 'n_odf         = ',odf_lib%n_odf


!=========================================================
! Hitlists erzeugen
!=========================================================

write(0,*) '----------------------------------------'
write(0,*) 'Building hitlists'
write(0,*) '----------------------------------------'

call make_hitlist(hits,odf_lib,nlevel_target)


write(0,*) 'total present hits = ', &
           count(hits%hits_present)

write(0,*) 'total lagged hits  = ', &
           count(hits%hits_lagged)


!=========================================================
! Pair map erzeugen
!=========================================================

write(0,*) '----------------------------------------'
write(0,*) 'Building LUT chunks'
write(0,*) '----------------------------------------'


call cpu_time(time_start)

call build_idx_chunks_odf( &
        prefix, &
        pair_map, &
        nobj, &
        npix, &
        hits, &
        max_chunk_size, &
        n_total, &
        n_chunks )

call cpu_time(time_end)


write(0,*) 'time for build_idx_chunks = ', &
           time_end-time_start

write(0,*) 'n_total  = ',n_total
write(0,*) 'n_chunks = ',n_chunks


!=========================================================
! Zufälliges Objektpaar
!=========================================================

call random_number(rn)

ri = floor(dble(nobj)*rn)+1

call random_number(rn)

rj = floor(dble(nobj)*rn)+1


write(0,*) '----------------------------------------'
write(0,*) 'Random object pair'
write(0,*) '----------------------------------------'

write(0,*) 'ri = ',ri
write(0,*) 'rj = ',rj

write(0,*) 'global start = ', &
           pair_map%start(ri,rj)

write(0,*) 'local start  = ', &
           pair_map%local_start(ri,rj)

write(0,*) 'size         = ', &
           pair_map%size(ri,rj)

write(0,*) 'chunk        = ', &
           pair_map%chunk(ri,rj)


!=========================================================
! Chunk laden
!=========================================================

call make_chunk_filename( &
        prefix, &
        pair_map%chunk(ri,rj), &
        filename )

write(0,*) 'reading chunk: ',trim(filename)

call read_chunk(filename,vec)


!=========================================================
! Block-Grenzen innerhalb des Chunks
!=========================================================

block_start = pair_map%local_start(ri,rj)

block_end = pair_map%local_start(ri,rj) + &
            pair_map%size(ri,rj) - 1


write(0,*) 'block start = ',block_start
write(0,*) 'block end   = ',block_end


!=========================================================
! Test
!=========================================================

total_time = 0.0d0


call cpu_time(time_start)


do ntest=1,1000


    !-----------------------------------------------------
    ! Zufällige gültige Kombination aus dem Support
    !
    ! idx entspricht:
    !
    ! X_ri(t+tau)
    ! X_ri(t)
    ! X_rj(t)
    !-----------------------------------------------------

    call get_random_pixel_kombi_from_hitlist_old( &
            npix, &
            ri, &
            rj, &
            hits, &
            idx )


    !-----------------------------------------------------
    ! Im LUT-Block suchen
    !-----------------------------------------------------

    call cpu_time(time_search_start)

    idx_in_chunk = binary_search( &
            vec, &
            block_start, &
            block_end, &
            idx )

    call cpu_time(time_search_end)


    total_time = total_time + &
                 (time_search_end-time_search_start)


    !-----------------------------------------------------
    ! Prüfen
    !-----------------------------------------------------

    if(idx_in_chunk < block_start .or. &
       idx_in_chunk > block_end) then

        write(0,*) 'ERROR: key not found'
        write(0,*) '  key         = ',idx
        write(0,*) '  result      = ',idx_in_chunk
        write(0,*) '  block start = ',block_start
        write(0,*) '  block end   = ',block_end

        stop

    endif


    if(idx /= vec(idx_in_chunk)) then

        write(0,*) 'ERROR: wrong LUT entry'
        write(0,*) '  key          = ',idx
        write(0,*) '  found        = ',vec(idx_in_chunk)
        write(0,*) '  position     = ',idx_in_chunk

        stop

    endif


enddo


call cpu_time(time_end)


!=========================================================
! Ergebnis
!=========================================================

write(0,*) '----------------------------------------'
write(0,*) 'Mapping test successful'
write(0,*) '----------------------------------------'

write(0,*) 'total loop time   = ', &
           time_end-time_start

write(0,*) 'total search time = ', &
           total_time

write(0,*) 'average search    = ', &
           total_time/dble(ntest)

write(0,*) '----------------------------------------'


deallocate(vec)

end subroutine test_mapping_real_old




subroutine prep_lut(odf_lib,nlevel_target,prefix_odf,path_odf)
implicit none

integer :: nlevel_target
type(odf_lib_type), intent(inout):: odf_lib
type(hitlist_type) :: hits
character*800 :: prefix,filename
integer :: max_chunk_size,n_chunks,no,npix
integer(int64) ::n_total
type(pair_map_type) :: paire_map
character(len=*),intent(in) :: prefix_odf,path_odf

call read_odf_lib(odf_lib,prefix_odf,path_odf)


write(0,*) odf_lib%nlevel

if(nlevel_target> odf_lib%nlevel) then
write(0,*) 'Error: Resolution of stored odf to low'	
stop
end if

write(0,*) odf_lib%n_odf
write(0,*) odf_lib%max_lag
no=size(odf_lib%odf(1)%hist,2)

call make_hitlist(hits,odf_lib,nlevel_target)

call  allocate_pair_map(paire_map,no)

prefix='test_idx_lut'
npix=2**(nlevel_target-1)
npix=npix*npix*12

max_chunk_size=1e8
call build_idx_chunks_odf(prefix,paire_map,no,npix,hits, &
                            max_chunk_size,n_total,n_chunks)

end subroutine


subroutine make_hitlist(hits,odf_lib,nlevel_target)

implicit none

type(odf_lib_type), intent(in) :: odf_lib
type(hitlist_type), intent(out) :: hits

integer :: iobj
integer :: iodf
integer :: ipix
integer :: ipix_target

integer :: nobj
integer :: npix_source
integer :: npix_target

integer :: nlevel_source
integer :: nlevel_target

integer :: level_diff
integer :: factor
integer :: i,k

!=========================================================
! Parameter
!=========================================================

nobj = size(odf_lib%odf(1)%hist,2)

nlevel_source = odf_lib%nlevel

if (nlevel_target > nlevel_source) then

    write(0,*) 'ERROR in make_hitlist:'
    write(0,*) 'nlevel_target > nlevel_source'
    write(0,*) '  source = ',nlevel_source
    write(0,*) '  target = ',nlevel_target

    stop

endif


!=========================================================
! Anzahl Pixel des Zielgitters
!=========================================================

npix_target = 2**(nlevel_target-1)
npix_target = 12*npix_target*npix_target


!=========================================================
! Speicher
!=========================================================

allocate(hits%hits_present(npix_target,nobj))
allocate(hits%hits_lagged (npix_target,nobj))

hits%hits_present = .false.
hits%hits_lagged  = .false.


!=========================================================
! Kein Downsampling
!=========================================================

if (nlevel_target == nlevel_source) then

    npix_source = odf_lib%odf(1)%npix


    !=====================================================
    ! Present ODF
    !=====================================================

    do iobj = 1,nobj

        do ipix = 0,npix_source-1

            if (odf_lib%odf(1)%hist(ipix+1,iobj) > 0) then

                hits%hits_present(ipix+1,iobj) = .true.

            endif

        enddo

    enddo


    !=====================================================
    ! Lagged ODFs
    !
    ! odf(2) ... odf(n_odf)
    !
    ! Die Supports werden vereinigt.
    !=====================================================

    do iodf = 2,odf_lib%n_odf

        npix_source = odf_lib%odf(iodf)%npix

        do iobj = 1,nobj

            do ipix = 0,npix_source-1

                if (odf_lib%odf(iodf)%hist(ipix+1,iobj) > 0) then

                    hits%hits_lagged(ipix+1,iobj) = .true.

                endif

            enddo

        enddo

    enddo


!=========================================================
! Nested Downsampling
!=========================================================

else

    level_diff = nlevel_source - nlevel_target

    factor = 4**level_diff


    write(0,*) 'HEALPix NESTED downsampling'
    write(0,*) '  source level = ',nlevel_source
    write(0,*) '  target level = ',nlevel_target
    write(0,*) '  factor       = ',factor


    !=====================================================
    ! Present ODF
    !=====================================================

    npix_source = odf_lib%odf(1)%npix

    do iobj = 1,nobj

        do ipix = 0,npix_source-1

            if (odf_lib%odf(1)%hist(ipix+1,iobj) > 0) then

                ipix_target = ipix / factor

                hits%hits_present(ipix_target+1,iobj) = .true.

            endif

        enddo

    enddo


    !=====================================================
    ! Lagged ODFs
    !
    ! Alle ODFs 2 ... n_odf werden auf das Target-Level
    ! abgebildet und logisch vereinigt.
    !=====================================================

    do iodf = 2,odf_lib%n_odf

        npix_source = odf_lib%odf(iodf)%npix

        do iobj = 1,nobj

            do ipix = 0,npix_source-1

                if (odf_lib%odf(iodf)%hist(ipix+1,iobj) > 0) then

                    ipix_target = ipix / factor

                    hits%hits_lagged(ipix_target+1,iobj) = .true.

                endif

            enddo

        enddo

    enddo

endif


!=========================================================
! Ausgabe zur Kontrolle
!=========================================================

do iobj = 1,nobj

    write(0,*) 'object ',iobj, &
               ' present support = ', &
               count(hits%hits_present(:,iobj))

    write(0,*) 'object ',iobj, &
               ' lagged support  = ', &
               count(hits%hits_lagged(:,iobj))

enddo

allocate(hits%rank_lagged(0:npix_target-1,nobj))
allocate(hits%rank_present(0:npix_target-1,nobj))

hits%rank_lagged  = 0
hits%rank_present = 0

do i=1,nobj

    k = 0

    do ipix=0,npix_target-1

        if(hits%hits_lagged(ipix+1,i)) then

            k = k + 1
            hits%rank_lagged(ipix,i) = k

        endif

    enddo


    k = 0

    do ipix=0,npix_target-1

        if(hits%hits_present(ipix+1,i)) then

            k = k + 1
            hits%rank_present(ipix,i) = k

        endif

    enddo

enddo

end subroutine make_hitlist



subroutine allocate_pair_map(map,no)

implicit none

type(pair_map_type),intent(out) :: map
integer,intent(in) :: no


allocate(map%start(no,no))
allocate(map%local_start(no,no))
allocate(map%size(no,no))
allocate(map%chunk(no,no))


end subroutine allocate_pair_map


subroutine test_mapping()
implicit none
integer, allocatable :: odf1(:,:),odf2(:,:),odf3(:,:),odf4(:,:)
logical, allocatable :: hit1(:,:),hit2(:,:),hit3(:,:),hit4(:,:)
integer ::nobj
integer ::nside
integer(i8) ::max_chunk_size,n_total
integer :: npix
integer :: nrand_try,ri,rj
character*800 :: prefix,filename
integer :: n_chunks,idx,idx_in_chunk,ntest
real(8) :: rn
integer, allocatable :: vec(:)
type(pair_map_type) :: pair_map
type(chunk_map_type):: chunk_map
integer :: block_start,block_end

real(8) :: time_start, time_end
real(8) :: time_search_start, time_search_end
real(8) :: total_time


nrand_try=20
nobj=100

nside=3
npix=(nside-1)
npix=4**npix
npix=12*npix

allocate(odf1(0:npix-1,nobj))
allocate(odf2(0:npix-1,nobj))
allocate(odf3(0:npix-1,nobj))
allocate(odf4(0:npix-1,nobj))

allocate(hit1(0:npix-1,nobj))
allocate(hit2(0:npix-1,nobj))
allocate(hit3(0:npix-1,nobj))
allocate(hit4(0:npix-1,nobj))

write(0,*) 'initialisation random test odf'
call  init_sparse_random_odf(nobj,npix,nrand_try,odf1)
call  init_sparse_random_odf(nobj,npix,nrand_try,odf2)
call  init_sparse_random_odf(nobj,npix,nrand_try,odf3)
call  init_sparse_random_odf(nobj,npix,nrand_try,odf4)

write(0,*) 'odf to logical array'
call get_hits(nobj,npix,odf1,hit1)
call get_hits(nobj,npix,odf2,hit2)
call get_hits(nobj,npix,odf3,hit3)
call get_hits(nobj,npix,odf4,hit4)

prefix='test_idx'

max_chunk_size=10**8
write(0,*) 'allocate pair map'
call allocate_pair_map(pair_map,nobj)
write(0,*) 'build idx chunks'

call cpu_time(time_start)
call build_idx_chunks(prefix,pair_map,nobj,npix,hit1,hit2,hit3,hit4, &
                            max_chunk_size,n_total,n_chunks)

call cpu_time(time_end)
write(0,*) 'time for build_idx_chunk',time_end-time_start
call random_number(rn)
ri=floor(dble(nobj)*rn)+1
call random_number(rn)
rj=floor(dble(nobj)*rn)+1
 
write(0,*)'random objek pair',ri,rj
write(0,*)'gobal start',pair_map%start(ri,rj)
write(0,*)'local start',pair_map%local_start(ri,rj)
write(0,*)'size',pair_map%size (ri,rj)
write(0,*)'chunk',pair_map%chunk(ri,rj)
write(0,*)'n_total',n_total 
write(0,*)'n_chunk',n_chunks
!write(0,*)  chunk_map%global_start(pair_map%chunk(ri,rj))
!write(0,*)  chunk_map%chunk_size  (pair_map%chunk(ri,rj))

call make_chunk_filename(prefix,pair_map%chunk(ri,rj),filename)
write(0,*) 'read chunk'
call read_chunk(filename,vec)

total_time = 0.0d0

call cpu_time(time_start)

do ntest=1,10**3
call get_random_pixel_kombi_from_hit_space(nobj,npix,ri,rj,hit1,hit2,hit3,hit4,idx)
 
!write(0,*) 'index',idx

!write(0,*) 'local start, local end',pair_map%local_start(ri,rj),pair_map%local_start(ri,rj)+pair_map%size(ri,rj)-1


block_start=pair_map%local_start(ri,rj)
block_end=pair_map%local_start(ri,rj)+pair_map%size(ri,rj)-1
call cpu_time(time_search_start)
idx_in_chunk= binary_search(vec,block_start,block_end,idx)

call cpu_time(time_search_end)
total_time = total_time + (time_search_end-time_search_start)
!write(0,*)'pos in chunk', idx_in_chunk
!write(0,*)'index from chunk', vec(idx_in_chunk)
if(idx/=vec(idx_in_chunk)) stop
end do

call cpu_time(time_end)

write(0,*) '--------------------------------'
write(0,*) 'total loop time:',time_end-time_start
write(0,*) 'total search time:',total_time
write(0,*) 'average search time:',total_time/dble(ntest)
write(0,*) '--------------------------------'

end subroutine




subroutine get_random_pixel_kombi_from_hit_space(nobj,npix,o1,o2,hit1,hit2,hit3,hit4,idx)
implicit none
integer,intent(in)  :: o1,o2,npix,nobj
integer,intent(out) :: idx
logical,intent(in) :: hit1(0:npix-1,nobj),hit2(0:npix-1,nobj),hit3(0:npix-1,nobj),hit4(0:npix-1,nobj)
logical :: not_found_combi
integer :: nmax_search,nsearch,i1,i2,i3,i4
real(8) :: rn

!write(0,*) 'Searching random kombi'
nmax_search=10**9

not_found_combi=.true.
nsearch=0

do while(not_found_combi .and. nsearch<nmax_search+1)
       nsearch=nsearch+1
      call random_number(rn)
      i1=floor(rn*dble(npix))
      call random_number(rn)
      i2=floor(rn*dble(npix))
      call random_number(rn)
      i3=floor(rn*dble(npix))


        if((hit1(i1,o1) .and. hit2(i2,o1) ).and. &
          ( hit3(i3,o2) .or.  hit4(i3,o2)) ) then


            idx=lut_key(i1,i2,i3,npix)
            not_found_combi=.false.
             
        endif
enddo

if(nsearch==nmax_search .and. not_found_combi) then
	write(0,*) 'no kombi'
	stop
end if	

end subroutine get_random_pixel_kombi_from_hit_space

subroutine get_hits(nobj,npix,odf,hit)
implicit none
integer,intent(in) :: odf(0:npix-1,nobj)
logical,intent(inout) :: hit(0:npix-1,nobj)
integer :: i,j,npix,nobj

hit=.false.

do i=1,size(odf,2)
do j=0,npix-1
if(odf(j,i)>0) hit(j,i)=.true.
end do
end do

end subroutine



subroutine init_sparse_random_odf(nobj,npix,n,odf)
implicit none
integer, intent(in) :: nobj,npix,n
integer :: i,j,k,ipix,t1,t2,t3
integer,intent(inout) :: odf(0:npix-1,nobj)
real(8) :: rn
odf=0



do j=1,n
  do i=1,nobj
    call random_number(rn)
    ipix=floor(rn*dble(npix))
    !write(0,*) ipix
    odf(ipix,i)=odf(ipix,i)+1
  end do
end do




end subroutine init_sparse_random_odf

subroutine build_pair_map_odf_red(pair_map, no, npix, hits, max_chunk_size, n_total, n_chunks)

    use iso_fortran_env, only : int64

    implicit none

    type(pair_map_type), intent(inout) :: pair_map
    type(hitlist_type),  intent(in)    :: hits

    integer, intent(in) :: no
    integer, intent(in) :: npix

    integer(int64), intent(out) :: n_total
    integer,        intent(out) :: n_chunks

    integer(4), intent(in) :: max_chunk_size

    ! Lokale Variablen
    integer :: i, j
    integer :: idx_pair
    integer(int64) :: pair_size
    integer(int64) :: global_idx
    integer(int64) :: chunk_entries
    integer :: chunk
    integer :: n_lagged(no), n_present(no)

    !=========================================================
    ! Konsistenzchecks
    !=========================================================
    if (size(hits%hits_lagged, 1) /= npix .or. size(hits%hits_present, 1) /= npix) then
        write(*,*) 'ERROR: hitlists have wrong npix'
        stop
    endif
    if (size(hits%hits_lagged, 2) /= no .or. size(hits%hits_present, 2) /= no) then
        write(*,*) 'ERROR: hitlists have wrong number of objects'
        stop
    endif

    !=========================================================
    ! Counts pro Objekt bestimmen
    !=========================================================
    do i = 1, no
        n_lagged(i)  = count(hits%hits_lagged(:, i))
        n_present(i) = count(hits%hits_present(:, i))
    end do

    !=========================================================
    ! Pair map Arrays allozieren & zurücksetzen (falls noch nicht geschehen)
    !=========================================================
    if (.not. allocated(pair_map%start)) then
        allocate(pair_map%start(no, no))
        allocate(pair_map%local_start(no, no))
        allocate(pair_map%size(no, no))
        allocate(pair_map%chunk(no, no))
    endif

    pair_map%start        = 0_int64
    pair_map%local_start = 0
    pair_map%size        = 0
    pair_map%chunk       = 0

    !=========================================================
    ! Chunk- & Startindizes direkt für die aktiven Paare berechnen
    ! (Nutzt die bereits vorgegebenen active_i / active_j Listen)
    !=========================================================
    chunk         = 1
    chunk_entries = 0_int64
    global_idx    = 1_int64

    do idx_pair = 1, pair_map%n_active

        i = pair_map%active_i(idx_pair)
        j = pair_map%active_j(idx_pair)

        ! Optionaler Sicherheitscheck, ob das aktive Paar auch Treffer besitzt
        if (n_lagged(i) == 0 .or. n_present(i) == 0 .or. n_present(j) == 0) then
            write(*,*) 'ERROR: Active pair (', i, ',', j, ') has zero hits!'
            stop
        endif

        pair_size = int(n_lagged(i), int64) * &
                    int(n_present(i), int64) * &
                    int(n_present(j), int64)

        if (pair_size <= 0_int64 .or. pair_size > max_chunk_size) then
            write(*,*) 'ERROR: Invalid pair size for (', i, ',', j, '):', pair_size
            stop
        endif

        ! Prüfen, ob neues Chunk gestartet werden muss
        if (chunk_entries > 0_int64 .and. &
            chunk_entries + pair_size > max_chunk_size) then
            chunk = chunk + 1
            chunk_entries = 0_int64
        endif

        pair_map%chunk(i, j)       = chunk
        pair_map%local_start(i, j) = int(chunk_entries + 1_int64)
        pair_map%start(i, j)       = global_idx  
        pair_map%size(i, j)        = int(pair_size)

        chunk_entries = chunk_entries + pair_size
        global_idx    = global_idx + pair_size

    end do

    n_chunks = chunk
    n_total  = global_idx - 1_int64

    write(0,*) '----------------------------------------'
    write(0,*) 'Pair Map built successfully (Reduced)'
    write(0,*) 'Active pairs = ', pair_map%n_active
    write(0,*) 'Total keys   = ', n_total
    write(0,*) 'Total chunks = ', n_chunks
    write(0,*) '----------------------------------------'

end subroutine build_pair_map_odf_red

subroutine build_pair_map_odf(pair_map, no, npix, hits, max_chunk_size, n_total, n_chunks)

    use iso_fortran_env, only : int64

    implicit none

    type(pair_map_type), intent(inout) :: pair_map
    type(hitlist_type),  intent(in)    :: hits

    integer, intent(in) :: no
    integer, intent(in) :: npix

    integer(int64), intent(out) :: n_total
    integer,        intent(out) :: n_chunks

    integer(4), intent(in) :: max_chunk_size

    ! Lokale Variablen
    integer :: i, j
    integer(int64) :: pair_size
    integer(int64) :: global_idx
    integer(int64) :: chunk_entries
    integer :: chunk
    integer :: n_active, active_count
    integer :: n_lagged(no), n_present(no)

    !=========================================================
    ! Konsistenzchecks
    !=========================================================
    if (size(hits%hits_lagged, 1) /= npix .or. size(hits%hits_present, 1) /= npix) then
        write(*,*) 'ERROR: hitlists have wrong npix'
        stop
    endif
    if (size(hits%hits_lagged, 2) /= no .or. size(hits%hits_present, 2) /= no) then
        write(*,*) 'ERROR: hitlists have wrong number of objects'
        stop
    endif

    !=========================================================
    ! Counts pro Objekt bestimmen
    !=========================================================
    do i = 1, no
        n_lagged(i)  = count(hits%hits_lagged(:, i))
        n_present(i) = count(hits%hits_present(:, i))
    end do

    !=========================================================
    ! Pair map allozieren & zurücksetzen
    !=========================================================
    if (.not. allocated(pair_map%start)) then
        allocate(pair_map%start(no, no))
        allocate(pair_map%local_start(no, no))
        allocate(pair_map%size(no, no))
        allocate(pair_map%chunk(no, no))
    endif

    pair_map%start       = 0_int64
    pair_map%local_start = 0
    pair_map%size        = 0
    pair_map%chunk       = 0

    !=========================================================
    ! Aktive Objektpaare ermitteln
    !=========================================================
    n_active = 0
    do i = 1, no
        if (n_lagged(i) == 0 .or. n_present(i) == 0) cycle
        do j = 1, no
            if (n_present(j) == 0) cycle
            n_active = n_active + 1
        end do
    end do

    if (allocated(pair_map%active_i)) deallocate(pair_map%active_i)
    if (allocated(pair_map%active_j)) deallocate(pair_map%active_j)

    allocate(pair_map%active_i(n_active))
    allocate(pair_map%active_j(n_active))
    pair_map%n_active = n_active

    active_count = 0
    do i = 1, no
        if (n_lagged(i) == 0 .or. n_present(i) == 0) cycle
        do j = 1, no
            if (n_present(j) == 0) cycle
            active_count = active_count + 1
            pair_map%active_i(active_count) = i
            pair_map%active_j(active_count) = j
        end do
    end do

    !=========================================================
    ! Chunk- & Startindizes berechnen (O(1) Pro Paar)
    !=========================================================
    chunk         = 1
    chunk_entries = 0_int64
    global_idx    = 1_int64

    do active_count = 1, pair_map%n_active

        i = pair_map%active_i(active_count)
        j = pair_map%active_j(active_count)

        pair_size = int(n_lagged(i), int64) * &
                    int(n_present(i), int64) * &
                    int(n_present(j), int64)

        if (pair_size <= 0_int64 .or. pair_size > max_chunk_size) then
            write(*,*) 'ERROR: Invalid pair size for (', i, ',', j, '):', pair_size
            stop
        endif

        ! Prüfen, ob neues Chunk gestartet werden muss
        if (chunk_entries > 0_int64 .and. &
            chunk_entries + pair_size > max_chunk_size) then
            chunk = chunk + 1
            chunk_entries = 0_int64
        endif

        pair_map%chunk(i, j)       = chunk
        pair_map%local_start(i, j) = int(chunk_entries + 1_int64)
        pair_map%start(i, j)       = global_idx  ! saubere 64-Bit Zuweisung
        pair_map%size(i, j)        = int(pair_size)

        chunk_entries = chunk_entries + pair_size
        global_idx    = global_idx + pair_size

    end do

    n_chunks = chunk
    n_total  = global_idx - 1_int64

    write(0,*) '----------------------------------------'
    write(0,*) 'Pair Map built successfully'
    write(0,*) 'Active pairs = ', pair_map%n_active
    write(0,*) 'Total keys   = ', n_total
    write(0,*) 'Total chunks = ', n_chunks
    write(0,*) '----------------------------------------'

end subroutine build_pair_map_odf

subroutine write_idx_chunks_odf(prefix, pair_map, no, npix, hits, max_chunk_size)

    use iso_fortran_env, only : int32, int64

    implicit none

    character(len=*),   intent(in) :: prefix
    type(pair_map_type), intent(in) :: pair_map
    type(hitlist_type),  intent(in) :: hits

    integer,    intent(in) :: no
    integer,    intent(in) :: npix
    integer(4), intent(in) :: max_chunk_size

    ! Lokale Variablen
    integer :: i, j, ip
    integer :: k1, k2, k3
    integer :: current_chunk
    integer :: active_count
    integer(int64) :: chunk_entries

    integer(int32) :: key
    integer(int32), allocatable :: chunk_vec(:)

    integer, allocatable :: pix_lagged(:,:)
    integer, allocatable :: pix_present(:,:)
    integer :: n_lagged(no), n_present(no)

    character(len=200) :: filename
    integer :: progress_step, next_progress, progress_percent

    !=========================================================
    ! Support-Listen aufbauen
    !=========================================================
    allocate(pix_lagged(npix, no))
    allocate(pix_present(npix, no))

    do i = 1, no
        n_lagged(i)  = 0
        n_present(i) = 0
        do ip = 0, npix - 1
            if (hits%hits_lagged(ip + 1, i)) then
                n_lagged(i) = n_lagged(i) + 1
                pix_lagged(n_lagged(i), i) = ip
            endif
            if (hits%hits_present(ip + 1, i)) then
                n_present(i) = n_present(i) + 1
                pix_present(n_present(i), i) = ip
            endif
        end do
    end do

    allocate(chunk_vec(int(max_chunk_size)))

    current_chunk = 1
    chunk_entries = 0_int64

    progress_step = max(1, pair_map%n_active / 10)
    next_progress = progress_step

    !=========================================================
    ! Chunks befüllen & speichern
    !=========================================================
    do active_count = 1, pair_map%n_active

        i = pair_map%active_i(active_count)
        j = pair_map%active_j(active_count)

        ! Wenn das Paar in ein neues Chunk gehört, das alte vorher wegschreiben
        if (pair_map%chunk(i, j) /= current_chunk) then
            call make_chunk_filename(prefix, current_chunk, filename)
            call write_chunk(filename, chunk_vec(1:int(chunk_entries)))

            current_chunk = pair_map%chunk(i, j)
            chunk_entries = 0_int64
        endif

        ! Keys berechnen & ablegen
        do k3 = 1, n_present(j)
            do k2 = 1, n_present(i)
                do k1 = 1, n_lagged(i)

                    key = lut_key( &
                            pix_lagged(k1, i), &
                            pix_present(k2, i), &
                            pix_present(k3, j), &
                            npix )

                    chunk_entries = chunk_entries + 1_int64
                    chunk_vec(int(chunk_entries)) = key

                end do
            end do
        end do

        ! Fortschritts-Anzeige
        if (active_count >= next_progress .or. active_count == pair_map%n_active) then
            progress_percent = int(100.0d0 * dble(active_count) / dble(pair_map%n_active))
            write(0,'(A,I3,A,I8,A,I8)') 'Writing LUT Chunks: ', progress_percent, ' %  (', &
                                        active_count, ' / ', pair_map%n_active, ')'
            next_progress = next_progress + progress_step
        endif

    end do

    ! Letzten Chunk schreiben
    if (chunk_entries > 0_int64) then
        call make_chunk_filename(prefix, current_chunk, filename)
        call write_chunk(filename, chunk_vec(1:int(chunk_entries)))
    endif

    deallocate(chunk_vec)
    deallocate(pix_lagged)
    deallocate(pix_present)

    write(0,*) '----------------------------------------'
    write(0,*) 'All Chunk Files written successfully'
    write(0,*) '----------------------------------------'

end subroutine write_idx_chunks_odf

subroutine build_idx_chunks_odf(prefix,pair_map,no,npix,hits, &
                                max_chunk_size,n_total,n_chunks)

use iso_fortran_env, only : int32,int64

implicit none

character(len=*),intent(in) :: prefix

type(pair_map_type),intent(inout) :: pair_map
type(hitlist_type),intent(in) :: hits

integer,intent(in) :: no
integer,intent(in) :: npix

integer(int64),intent(out) :: n_total
integer,intent(out) :: n_chunks

integer(4),intent(in) :: max_chunk_size


!=========================================================
! lokale Variablen
!=========================================================

integer :: i,j
integer :: k1,k2,k3
integer :: chunk

integer(int64) :: pair_size
integer(int64) :: global_idx
integer(int64) :: chunk_entries

integer(int32) :: key

integer(int32),allocatable :: chunk_vec(:)

integer,allocatable :: pix_lagged(:,:)
integer,allocatable :: pix_present(:,:)

integer :: n_lagged(no)
integer :: n_present(no)

integer :: ip

integer :: n_active
integer :: active_count

character(len=200) :: filename
integer :: progress_step
integer :: next_progress
integer :: progress_percent

!=========================================================
! Konsistenzcheck
!=========================================================

if(size(hits%hits_lagged,1) /= npix) then

    write(*,*) 'ERROR: lagged hitlist has wrong npix'
    write(*,*) '  expected = ',npix
    write(*,*) '  actual   = ',size(hits%hits_lagged,1)

    stop

endif


if(size(hits%hits_present,1) /= npix) then

    write(*,*) 'ERROR: present hitlist has wrong npix'
    write(*,*) '  expected = ',npix
    write(*,*) '  actual   = ',size(hits%hits_present,1)

    stop

endif


if(size(hits%hits_lagged,2) /= no) then

    write(*,*) 'ERROR: lagged hitlist has wrong number of objects'
    write(*,*) '  expected = ',no
    write(*,*) '  actual   = ',size(hits%hits_lagged,2)

    stop

endif


if(size(hits%hits_present,2) /= no) then

    write(*,*) 'ERROR: present hitlist has wrong number of objects'
    write(*,*) '  expected = ',no
    write(*,*) '  actual   = ',size(hits%hits_present,2)

    stop

endif


!=========================================================
! Pair map
!=========================================================

if(.not.allocated(pair_map%start)) then

    allocate(pair_map%start(no,no))
    allocate(pair_map%local_start(no,no))
    allocate(pair_map%size(no,no))
    allocate(pair_map%chunk(no,no))

endif


pair_map%start       = 0
pair_map%local_start = 0
pair_map%size        = 0
pair_map%chunk       = 0


!=========================================================
! Support-Listen erzeugen
!=========================================================

allocate(pix_lagged(npix,no))
allocate(pix_present(npix,no))


do i=1,no

    n_lagged(i)  = 0
    n_present(i) = 0


    do ip=0,npix-1

        !-------------------------------------------------
        ! X_i(t+tau)
        !-------------------------------------------------

        if(hits%hits_lagged(ip+1,i)) then

            n_lagged(i) = n_lagged(i)+1

            pix_lagged(n_lagged(i),i) = ip

        endif


        !-------------------------------------------------
        ! X_i(t)
        !-------------------------------------------------

        if(hits%hits_present(ip+1,i)) then

            n_present(i) = n_present(i)+1

            pix_present(n_present(i),i) = ip

        endif

    enddo

enddo


!=========================================================
! Anzahl aktiver Objektpaare bestimmen
!
! Ein Paar (i,j) ist aktiv, wenn
!
! n_lagged(i)  > 0
! n_present(i) > 0
! n_present(j) > 0
!
!=========================================================

n_active = 0


do i=1,no

    if(n_lagged(i) == 0) cycle
    if(n_present(i) == 0) cycle

    do j=1,no

        if(n_present(j) == 0) cycle

        n_active = n_active + 1

    enddo

enddo


!=========================================================
! Aktive-Paar-Liste allozieren
!=========================================================

if(allocated(pair_map%active_i)) then
    deallocate(pair_map%active_i)
endif

if(allocated(pair_map%active_j)) then
    deallocate(pair_map%active_j)
endif


if(n_active > 0) then

    allocate(pair_map%active_i(n_active))
    allocate(pair_map%active_j(n_active))

else

    allocate(pair_map%active_i(0))
    allocate(pair_map%active_j(0))

endif


pair_map%n_active = n_active


!=========================================================
! Aktive-Paar-Liste füllen
!=========================================================

active_count = 0


do i=1,no

    if(n_lagged(i) == 0) cycle
    if(n_present(i) == 0) cycle

    do j=1,no

        if(n_present(j) == 0) cycle

        active_count = active_count + 1

        pair_map%active_i(active_count) = i
        pair_map%active_j(active_count) = j

    enddo

enddo


!=========================================================
! Information
!=========================================================

write(0,*) '----------------------------------------'
write(0,*) 'LUT active object pairs'
write(0,*) '----------------------------------------'
write(0,*) 'number of objects = ',no
write(0,*) 'active pairs      = ',pair_map%n_active
write(0,*) 'possible pairs    = ',no*no
write(0,*) '----------------------------------------'


!=========================================================
! Speicher für Chunk
!=========================================================

allocate(chunk_vec(int(max_chunk_size)))


!=========================================================
! Initialisierung
!=========================================================

chunk         = 1
chunk_entries = 0_int64
global_idx    = 1_int64


!=========================================================
! Nur aktive Objektpaare bearbeiten
!=========================================================

do active_count=1,pair_map%n_active

    i = pair_map%active_i(active_count)
    j = pair_map%active_j(active_count)


    !-----------------------------------------------------
    ! Größe des Triplett-Produkts
    !
    ! X_i(t+tau) X_i(t) X_j(t)
    !-----------------------------------------------------

    pair_size = int(n_lagged(i),int64) * &
                int(n_present(i),int64) * &
                int(n_present(j),int64)


    !-----------------------------------------------------
    ! Sicherheitscheck
    !-----------------------------------------------------

    if(pair_size <= 0_int64) then

        write(*,*) 'ERROR: active pair has zero size'
        write(*,*) '  i = ',i
        write(*,*) '  j = ',j
        stop

    endif


    if(pair_size > max_chunk_size) then

        write(*,*) 'ERROR: pair larger than chunk'
        write(*,*) '  object i = ',i
        write(*,*) '  object j = ',j
        write(*,*) '  lagged   = ',n_lagged(i)
        write(*,*) '  present i= ',n_present(i)
        write(*,*) '  present j= ',n_present(j)
        write(*,*) '  pair size= ',pair_size
        write(*,*) '  max chunk= ',max_chunk_size

        stop

    endif


    !-----------------------------------------------------
    ! Neuen Chunk beginnen, falls nötig
    !-----------------------------------------------------

    if(chunk_entries > 0_int64 .and. &
       chunk_entries + pair_size > max_chunk_size) then

        call make_chunk_filename( &
            prefix, &
            chunk, &
            filename )

        call write_chunk( &
            filename, &
            chunk_vec(1:int(chunk_entries)) )


        chunk = chunk + 1

        chunk_entries = 0_int64

    endif


    !-----------------------------------------------------
    ! Pair Map
    !-----------------------------------------------------

    pair_map%chunk(i,j) = chunk

    pair_map%local_start(i,j) = &
        int(chunk_entries + 1)

    pair_map%start(i,j) = &
         global_idx

    pair_map%size(i,j) = &
        int(pair_size)


    !=====================================================
    ! Triplett erzeugen
    !
    ! k1 = X_i(t+tau)
    ! k2 = X_i(t)
    ! k3 = X_j(t)
    !=====================================================

    do k3=1,n_present(j)

        do k2=1,n_present(i)

            do k1=1,n_lagged(i)

                key = lut_key( &
                    pix_lagged(k1,i), &
                    pix_present(k2,i), &
                    pix_present(k3,j), &
                    npix )


                chunk_entries = chunk_entries + 1_int64

                chunk_vec(int(chunk_entries)) = key

            enddo

        enddo

    enddo


    !-----------------------------------------------------
    ! Globalen Index weitersetzen
    !-----------------------------------------------------

    global_idx = global_idx + pair_size
    
        !=====================================================
    ! Fortschritt
    !=====================================================

    if(active_count >= next_progress .or. &
       active_count == pair_map%n_active) then

        progress_percent = int( &
            100.0d0 * dble(active_count) / &
            dble(pair_map%n_active) )

        write(0,'(A,I3,A,I8,A,I8)') &
            'LUT: ',progress_percent,' %  (', &
            active_count,' / ',pair_map%n_active,')'

        next_progress = next_progress + progress_step

    endif

enddo


!=========================================================
! Letzten Chunk schreiben
!=========================================================

if(chunk_entries > 0_int64) then

    call make_chunk_filename( &
        prefix, &
        chunk, &
        filename )

    call write_chunk( &
        filename, &
        chunk_vec(1:int(chunk_entries)) )

endif


!=========================================================
! Ergebnis
!=========================================================

n_chunks = chunk

n_total = global_idx - 1_int64


!=========================================================
! Speicher freigeben
!=========================================================

deallocate(chunk_vec)
deallocate(pix_lagged)
deallocate(pix_present)


write(0,*) '----------------------------------------'
write(0,*) 'LUT construction finished'
write(0,*) 'active pairs = ',pair_map%n_active
write(0,*) 'total keys   = ',n_total
write(0,*) 'chunks       = ',n_chunks
write(0,*) '----------------------------------------'


end subroutine build_idx_chunks_odf




subroutine build_idx_chunks(prefix,pair_map,no,npix,hits1,hits2,hits3,hits_hyp, &
                            max_chunk_size,n_total,n_chunks)

use iso_fortran_env, only : int32,int64

implicit none

character(len=*),intent(in) :: prefix
type(pair_map_type),intent(inout) :: pair_map

integer,intent(in) :: no,npix
integer(int64),intent(out) :: n_total
integer,intent(out) :: n_chunks
integer(int64),intent(in) :: max_chunk_size

logical,intent(in) :: hits1(0:npix-1,no)
logical,intent(in) :: hits2(0:npix-1,no)
logical,intent(in) :: hits3(0:npix-1,no)
logical,intent(in) :: hits_hyp(0:npix-1,no)

integer :: i,j
integer :: k1,k2,k3
integer :: chunk
integer(int64) :: pair_size
integer(int64) :: global_idx
integer(int64) :: chunk_entries

integer(int32) :: key

integer(int32),allocatable :: chunk_vec(:)

integer,allocatable :: pix1(:,:),pix2(:,:),pix3(:,:)
integer :: n1(no),n2(no),n3(no)

integer :: ip

character(len=200) :: filename

!---------------------------------------------------------
! pair map
!---------------------------------------------------------

if(.not.allocated(pair_map%start)) then

    allocate(pair_map%start(no,no))
    allocate(pair_map%local_start(no,no))
    allocate(pair_map%size(no,no))
    allocate(pair_map%chunk(no,no))

endif

pair_map%start       = 0
pair_map%local_start = 0
pair_map%size        = 0
pair_map%chunk       = 0

!---------------------------------------------------------
! Speicher
!---------------------------------------------------------

allocate(chunk_vec(int(max_chunk_size)))

allocate(pix1(npix,no))
allocate(pix2(npix,no))
allocate(pix3(npix,no))

!---------------------------------------------------------
! Trefferlisten einmal erzeugen
!---------------------------------------------------------

do i=1,no

    n1(i)=0
    n2(i)=0
    n3(i)=0

    do ip=0,npix-1

        if(hits1(ip,i)) then
            n1(i)=n1(i)+1
            pix1(n1(i),i)=ip
        endif

        if(hits2(ip,i)) then
            n2(i)=n2(i)+1
            pix2(n2(i),i)=ip
        endif

        if(hits3(ip,i) .or. hits_hyp(ip,i)) then
            n3(i)=n3(i)+1
            pix3(n3(i),i)=ip
        endif

    enddo

enddo

!---------------------------------------------------------
! Initialisierung
!---------------------------------------------------------

chunk=1
chunk_entries=0_int64
global_idx=1_int64

!---------------------------------------------------------
! Objektpaare
!---------------------------------------------------------

do i=1,no
do j=1,no

    pair_size=int(n1(i),int64) * &
              int(n2(i),int64) * &
              int(n3(j),int64)

    if(pair_size>max_chunk_size) then
        write(*,*) "pair larger than chunk"
        stop
    endif

    if(chunk_entries>0_int64 .and. &
       chunk_entries+pair_size>max_chunk_size) then

        call make_chunk_filename(prefix,chunk,filename)

        call write_chunk(filename,chunk_vec(1:int(chunk_entries)))

        chunk=chunk+1
        chunk_entries=0_int64

    endif

    pair_map%chunk(i,j)=chunk
    pair_map%local_start(i,j)=chunk_entries+1
    pair_map%start(i,j)=global_idx
    pair_map%size(i,j)=pair_size

    do k3=1,n3(j)
    do k2=1,n2(i)
    do k1=1,n1(i)

        key=lut_key( &
             pix1(k1,i), &
             pix2(k2,i), &
             pix3(k3,j), &
             npix)

        chunk_entries=chunk_entries+1

        chunk_vec(chunk_entries)=key

    enddo
    enddo
    enddo

    global_idx=global_idx+pair_size

enddo
enddo

!---------------------------------------------------------
! letzten Chunk schreiben
!---------------------------------------------------------

if(chunk_entries>0) then

    call make_chunk_filename(prefix,chunk,filename)

    call write_chunk(filename,chunk_vec(1:int(chunk_entries)))

endif

n_chunks=chunk
n_total=global_idx-1

deallocate(chunk_vec)
deallocate(pix1,pix2,pix3)

end subroutine


subroutine list_pix(no,npix,hits1,hits2,hits3,hits_hyp)
implicit none
logical, intent(in), dimension(:,:) ::hits1,hits2,hits3,hits_hyp
integer :: no1,no2,no,npix
integer ::i1,i2,i3,i,j
integer :: idx
integer, allocatable :: idx_vec(:)

idx=0
do i=1,no
do j=1,no
do i1=0,npix-1
do i2=0,npix-1
do i3=0,npix-1

if((hits1(i1,i) .and. (hits2(i2,i)) .and. (hits3(i3,j) &
.or. hits_hyp(i3,j) ))) then
idx=idx+1
end if

end do
end do
end do
end do
end do
write(0,*) idx,'oben'

allocate(idx_vec(idx))

idx=1
do i=1,no
do j=1,no
!pair_kombi_start(i,j)=idx

do i3=0,npix-1
do i2=0,npix-1
do i1=0,npix-1

if(hits1(i1,i) .and. hits2(i2,i) .and. hits3(i3,j) ) then

idx_vec(idx)=i1+npix*i2+npix**2*i3

idx=idx+1
end if

end do
end do
end do

!pair_kombi_size(i,j)=idx-pair_kombi_start(i,j)

end do
end do

write(0,*) idx,'unten'

end subroutine


subroutine write_chunk(filename,vec)

use iso_fortran_env, only : int32

implicit none

character(len=*), intent(in) :: filename
integer(int32), intent(in)   :: vec(:)

integer :: iu
integer(int32) :: n


n=size(vec)


open(newunit=iu, &
     file=filename, &
     access="stream", &
     form="unformatted", &
     status="replace", &
     action="write")


write(iu) n
write(iu) vec


close(iu)


end subroutine write_chunk


subroutine read_chunk(filename,vec)

use iso_fortran_env, only : int32

implicit none

character(len=*), intent(in) :: filename
integer(int32), allocatable, intent(out) :: vec(:)

integer :: iu
integer(int32) :: n


open(newunit=iu, &
     file=filename, &
     access="stream", &
     form="unformatted", &
     status="old", &
     action="read")


read(iu) n


allocate(vec(n))


read(iu) vec


close(iu)


end subroutine read_chunk



pure integer function binary_search(keys,first,last,key) result(pos)

implicit none

integer,intent(in) :: keys(:)
integer,intent(in) :: first,last,key

integer :: left,right,mid,value


left=first
right=last

pos=-1


do while(left<=right)

    mid=left+(right-left)/2

    value=keys(mid)


    if(value < key) then

        left=mid+1


    else if(value > key) then

        right=mid-1


    else

        pos=mid
        return

    endif


enddo


end function binary_search


pure integer function lut_key(i1,i2,i3,npix)

implicit none

integer,intent(in) :: i1,i2,i3,npix


lut_key=i1+npix*i2+npix*npix*i3


end function lut_key

subroutine make_chunk_filename(prefix,chunk,filename)

implicit none

character(len=*), intent(in)  :: prefix
integer, intent(in)            :: chunk
character(len=*), intent(out)  :: filename


write(filename,'(A,"_chunk_",I6.6,".bin")') trim(prefix),chunk


end subroutine make_chunk_filename


subroutine get_random_pixel_kombi_from_hitlist_old( &
            npix,ri,rj,hits,idx)

use iso_fortran_env, only : int32

implicit none

integer,intent(in) :: npix
integer,intent(in) :: ri
integer,intent(in) :: rj

type(hitlist_type),intent(in) :: hits

integer(int32),intent(out) :: idx

integer :: ipix1
integer :: ipix2
integer :: ipix3

integer :: n1
integer :: n2
integer :: n3

integer :: k

real(8) :: rn


!=========================================================
! Anzahl gültiger Pixel
!=========================================================

n1 = count(hits%hits_lagged(:,ri))

n2 = count(hits%hits_present(:,ri))

n3 = count(hits%hits_present(:,rj))


if(n1 == 0) then
    write(0,*) 'ERROR: no lagged support for object ',ri
    stop
endif

if(n2 == 0) then
    write(0,*) 'ERROR: no present support for object ',ri
    stop
endif

if(n3 == 0) then
    write(0,*) 'ERROR: no present support for object ',rj
    stop
endif


!=========================================================
! Zufälliges Pixel aus lagged support
!=========================================================
call random_number(rn)
k = 1 + int(n1*rn)

ipix1 = 0

do while(k > 0)

    if(hits%hits_lagged(ipix1+1,ri)) then
        k = k-1
    endif

    if(k > 0) ipix1 = ipix1+1

enddo


!=========================================================
! Zufälliges Pixel aus present support von ri
!=========================================================
call random_number(rn)
k = 1 + int(n2*rn)

ipix2 = 0

do while(k > 0)

    if(hits%hits_present(ipix2+1,ri)) then
        k = k-1
    endif

    if(k > 0) ipix2 = ipix2+1

enddo


!=========================================================
! Zufälliges Pixel aus present support von rj
!=========================================================
call random_number(rn)
k = 1 + int(n3*rn)

ipix3 = 0

do while(k > 0)

    if(hits%hits_present(ipix3+1,rj)) then
        k = k-1
    endif

    if(k > 0) ipix3 = ipix3+1

enddo


!=========================================================
! LUT Key
!
!   ipix1 = X_ri(t+tau)
!   ipix2 = X_ri(t)
!   ipix3 = X_rj(t)
!=========================================================

idx = lut_key( &
        ipix1, &
        ipix2, &
        ipix3, &
        npix )


end subroutine get_random_pixel_kombi_from_hitlist_old

subroutine get_random_pixel_kombi_from_hitlist( &
            npix,ri,rj,hits,ipix1,ipix2,ipix3,idx)

use iso_fortran_env, only : int32

implicit none

integer,intent(in) :: npix
integer,intent(in) :: ri
integer,intent(in) :: rj

type(hitlist_type),intent(in) :: hits

integer(int32),intent(out) :: idx

integer,intent(out) :: ipix1
integer,intent(out) :: ipix2
integer :: ipix3

integer :: n1
integer :: n2
integer :: n3

integer :: k

real(8) :: rn


!=========================================================
! Anzahl gültiger Pixel
!=========================================================

n1 = count(hits%hits_lagged(:,ri))

n2 = count(hits%hits_present(:,ri))

n3 = count(hits%hits_present(:,rj))


if(n1 == 0) then
    write(0,*) 'ERROR: no lagged support for object ',ri
    stop
endif

if(n2 == 0) then
    write(0,*) 'ERROR: no present support for object ',ri
    stop
endif

if(n3 == 0) then
    write(0,*) 'ERROR: no present support for object ',rj
    stop
endif


!=========================================================
! Zufälliges Pixel aus lagged support
!=========================================================
call random_number(rn)
k = 1 + int(n1*rn)

ipix1 = 0

do while(k > 0)

    if(hits%hits_lagged(ipix1+1,ri)) then
        k = k-1
    endif

    if(k > 0) ipix1 = ipix1+1

enddo


!=========================================================
! Zufälliges Pixel aus present support von ri
!=========================================================
call random_number(rn)
k = 1 + int(n2*rn)

ipix2 = 0

do while(k > 0)

    if(hits%hits_present(ipix2+1,ri)) then
        k = k-1
    endif

    if(k > 0) ipix2 = ipix2+1

enddo


!=========================================================
! Zufälliges Pixel aus present support von rj
!=========================================================
call random_number(rn)
k = 1 + int(n3*rn)

ipix3 = 0

do while(k > 0)

    if(hits%hits_present(ipix3+1,rj)) then
        k = k-1
    endif

    if(k > 0) ipix3 = ipix3+1

enddo


!=========================================================
! LUT Key
!
!   ipix1 = X_ri(t+tau)
!   ipix2 = X_ri(t)
!   ipix3 = X_rj(t)
!=========================================================

idx = lut_key( &
        ipix1, &
        ipix2, &
        ipix3, &
        npix )


end subroutine get_random_pixel_kombi_from_hitlist

end module sparse_lut_module