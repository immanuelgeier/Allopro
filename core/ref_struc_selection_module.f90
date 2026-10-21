module ref_struc_selection_module

use project_module
use topology_module
use direction_vectors_analyse_plan_module
use position_lift_module
use time_window_file_mapping
use time_window_loader
use domains_module
use sidechainhead_refstruc
use domain_refstruc
use  kabsch_module
use user_defined_cv_module

implicit none

contains

subroutine ref_struc_fixed_user(project, topologies, sc_lib, dom_lib, user_lib, istate, irep, itime)
    use user_defined_cv_module, only: user_cv_library_type
    implicit none

    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(domain_library_type), intent(inout) :: dom_lib
    type(sidechain_library_type), intent(inout) :: sc_lib
    type(user_cv_library_type), intent(inout) :: user_lib
    
    integer, intent(in) :: istate
    integer, intent(in) :: irep
    integer, intent(in) :: itime

    type(time_window) :: window1
    logical, allocatable :: window_mask1(:)
    type(raw_data_block) :: block
    real(8), allocatable :: xyz(:, :), local_cv(:, :), cv(:, :)
    real(8) :: d1_orig(3), d2_orig(3), d1_ref(3), d2_ref(3)
    real(8) :: v11(3), v12(3), v21(3), v22(3)
    integer :: i, icv, jres
    integer :: iatom, natom_user, res_j
    character(len=4) :: atom_name
    character(len=32) :: ref_method
    logical :: f11, f12, f21, f22

    ! ---------------------------------------------------------
    ! 1. Time Window und Maske für den spezifischen Frame bauen
    ! ---------------------------------------------------------
    window1%istate = istate
    window1%irep = irep
    window1%start_time_step = itime
    window1%end_time_step = itime

    allocate(window_mask1(1))
    window_mask1(1) = .true.

    call load_raw_data_block(project, window1, window_mask1, block)

    write(0,*) '-> Setting fixed reference from state=', istate, ' rep=', irep, ' time=', itime


    ! ==========================================================
    ! 2. DOMAINS
    ! ==========================================================
    do icv = 1, dom_lib%ndo
        ! Lokaler Frame-Index im Block ist 1
        call get_dom_ref(1, irep, istate, icv, dom_lib, topologies, project, block, xyz)
        call direction_vecs_ref_inv_gyration(size(xyz, 2), xyz, .true., d1_orig, d2_orig, cv)

        dom_lib%domain(icv)%ref_xyz = xyz
        dom_lib%domain(icv)%d1_ref = d1_orig
        dom_lib%domain(icv)%d2_ref = d2_orig
        dom_lib%domain(icv)%cv = cv
        
        deallocate(xyz, cv)
    end do


    ! ==========================================================
    ! 3. SIDECHAINS (Standardmäßig Instanz jres = 1)
    jres = 1 
    ! ==========================================================
    do icv = 1, sc_lib%ncv
        if (sc_lib%residue(icv)%freq == 0) cycle

        call get_res_ref(1, istate, icv, jres, sc_lib, topologies, block, xyz)
        res_j = sc_lib%residue(icv)%res_ids(jres)

        v11 = 0.0d0; v12 = 0.0d0; v21 = 0.0d0; v22 = 0.0d0
        f11 = .false.; f12 = .false.; f21 = .false.; f22 = .false.

        if (.not. sc_lib%residue(icv)%has_sym) then
            do i = 1, sc_lib%residue(icv)%natom_head
                iatom = sc_lib%residue(icv)%map(i) + topologies%state(istate)%topo%acid_start_list(res_j) - 1
                atom_name = topologies%state(istate)%topo%atom_id(iatom)

                if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then; v11 = xyz(:, i); f11 = .true.; end if
                if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then; v12 = xyz(:, i); f12 = .true.; end if
                if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then; v21 = xyz(:, i); f21 = .true.; end if
                if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then; v22 = xyz(:, i); f22 = .true.; end if
            end do

            if (.not. f11 .or. .not. f12 .or. .not. f21 .or. .not. f22) then
                error stop 'Reference direction atom not found for sidechain'
            end if

            d1_ref = v12 - v11
            d2_ref = v22 - v21
            call direction_vecs_ref_inv(sc_lib%residue(icv)%natom_head, xyz, d1_ref, d2_ref, local_cv)
        else
            call direction_vecs_sym(sc_lib%residue(icv)%natom_head, xyz, sc_lib%residue(icv)%cv_id, local_cv, sc_lib%residue(icv)%order)
            d1_ref = sc_lib%residue(icv)%d1_ref
            d2_ref = sc_lib%residue(icv)%d2_ref
        end if

        sc_lib%residue(icv)%ref_xyz = xyz
        sc_lib%residue(icv)%d1_ref = d1_ref
        sc_lib%residue(icv)%d2_ref = d2_ref
        sc_lib%residue(icv)%cv = local_cv

        deallocate(xyz, local_cv)
    end do


    ! ==========================================================
    ! 4. USER-DEFINED CVs
    ! ==========================================================
    do icv = 1, user_lib%n_cvs
        if (.not. user_lib%cvs(icv)%active) cycle

        natom_user = user_lib%cvs(icv)%natom
        if (natom_user < 1) cycle

        allocate(xyz(3, natom_user))
        call get_user_cv_ref(1, istate, icv, user_lib, topologies, block, xyz)

        ref_method = trim(user_lib%cvs(icv)%reference_method)

        select case (ref_method)
        case ('atom_pair')
            v11 = xyz(:, user_lib%cvs(icv)%dir1_atom_1)
            v12 = xyz(:, user_lib%cvs(icv)%dir1_atom_2)
            v21 = xyz(:, user_lib%cvs(icv)%dir2_atom_1)
            v22 = xyz(:, user_lib%cvs(icv)%dir2_atom_2)
            
            d1_ref = v12 - v11
            d2_ref = v22 - v21

            allocate(local_cv(2, natom_user))
            call direction_vecs_ref_inv(natom_user, xyz, d1_ref, d2_ref, local_cv)
            
            user_lib%cvs(icv)%ref_xyz = xyz
            if (allocated(user_lib%cvs(icv)%cv_values)) deallocate(user_lib%cvs(icv)%cv_values)
            allocate(user_lib%cvs(icv)%cv_values(2, natom_user))
            user_lib%cvs(icv)%cv_values = local_cv
            deallocate(local_cv)

        case default ! 'gyration'
            call direction_vecs_ref_inv_gyration(natom_user, xyz, .true., d1_orig, d2_orig, cv)
            user_lib%cvs(icv)%ref_xyz = xyz
            if (allocated(user_lib%cvs(icv)%cv_values)) deallocate(user_lib%cvs(icv)%cv_values)
            allocate(user_lib%cvs(icv)%cv_values(size(cv, 1), size(cv, 2)))
            user_lib%cvs(icv)%cv_values = cv
            deallocate(cv)
        end select

        deallocate(xyz)
    end do

end subroutine ref_struc_fixed_user

subroutine ref_struc_search_rmsd_user( &
    project, topologies, sc_lib, dom_lib, user_lib, window1, &
    nsample, method_window)

    use user_defined_cv_module, only: user_cv_library_type
    implicit none

    type(time_window), intent(in) :: window1
    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(domain_library_type), intent(inout) :: dom_lib
    type(sidechain_library_type), intent(inout) :: sc_lib
    type(user_cv_library_type), intent(inout) :: user_lib

    integer, intent(in) :: nsample
    integer, intent(in) :: method_window

    logical, allocatable :: window_mask1(:)
    type(raw_data_block) :: block

    if (nsample < 1) then
        error stop 'nsample must be >= 1'
    end if

    call fill_window_mask(window1, window_mask1, method_window, nsample)
    call load_raw_data_block(project, window1, window_mask1, block)

    ! Modularer Aufruf der einzelnen Sektionen
    call search_domain_references(project, topologies, dom_lib, window1, nsample, block)
    call search_sidechain_references(sc_lib, topologies, window1, nsample, block)
    call search_user_cv_references(user_lib, topologies, window1, nsample, block)

end subroutine ref_struc_search_rmsd_user

subroutine search_domain_references(project, topologies, dom_lib, window1, nsample, block)
    use domains_module, only: domain_library_type
    implicit none
    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(domain_library_type), intent(inout) :: dom_lib
    type(time_window), intent(in) :: window1
    integer, intent(in) :: nsample
    type(raw_data_block), intent(in) :: block

    real(8), allocatable :: xyz(:,:), all_dom_xyz(:,:,:), cv(:,:)
    real(8) :: d1_orig(3), d2_orig(3), rot(3,3), rmsd
    real(8) :: sum_score, mean_score, min_score
    integer :: icv, itimes, itime2, natom_dom, best_j, irep

    do icv = 1, dom_lib%ndo
        write(0,*) 'RMSD search for domain: ', trim(dom_lib%domain(icv)%domain_name)

        call get_dom_ref(1, irep, window1%istate, icv, dom_lib, topologies, project, block, xyz)
        natom_dom = size(xyz, 2)

        allocate(all_dom_xyz(3, natom_dom, nsample))
        all_dom_xyz(:,:,1) = xyz

        do itimes = 2, nsample
            call get_dom_ref(itimes, irep, window1%istate, icv, dom_lib, topologies, project, block, xyz)
            all_dom_xyz(:,:,itimes) = xyz
        end do

        min_score = huge(1.0d0)
        best_j = 0

        do itimes = 1, nsample
            sum_score = 0.0d0
            do itime2 = 1, nsample
                if (itime2 == itimes) cycle
                call kabsch_rotation(natom_dom, all_dom_xyz(:,:,itimes), all_dom_xyz(:,:,itime2), rot, rmsd)
                sum_score = sum_score + rmsd
            end do
            mean_score = sum_score / dble(nsample - 1)
            if (mean_score < min_score) then
                min_score = mean_score
                best_j = itimes
            end if
        end do

        xyz = all_dom_xyz(:,:,best_j)
        call direction_vecs_ref_inv_gyration(size(xyz,2), xyz, .true., d1_orig, d2_orig, cv)

        dom_lib%domain(icv)%ref_xyz = xyz
        dom_lib%domain(icv)%d1_ref = d1_orig
        dom_lib%domain(icv)%d2_ref = d2_orig
        dom_lib%domain(icv)%cv = cv

        deallocate(all_dom_xyz)
        deallocate(cv)
    end do
end subroutine search_domain_references


subroutine search_sidechain_references(sc_lib, topologies, window1, nsample, block)
    use sidechains, only: sidechain_library_type
    implicit none
    type(sidechain_library_type), intent(inout) :: sc_lib
    type(topology_set_type), intent(in) :: topologies
    type(time_window), intent(in) :: window1
    integer, intent(in) :: nsample
    type(raw_data_block), intent(in) :: block

    real(8), allocatable :: xyz(:,:), all_res_xyz(:,:,:,:), local_cv(:,:)
    real(8) :: d1_ref(3), d2_ref(3), v11(3), v12(3), v21(3), v22(3)
    real(8) :: rot(3,3), rmsd, sum_score, mean_score, min_score
    integer :: icv, ires, jres, itimes, itime2, natom_res, freq, best_j, best_jres, res_j, j, iatom
    character(len=4) :: atom_name
    logical :: f11, f12, f21, f22

    do icv = 1, sc_lib%ncv
        freq = sc_lib%residue(icv)%freq
        if (freq == 0) cycle

        write(0,*) 'RMSD search for residue type: ', trim(sc_lib%residue(icv)%resname)

        call get_res_ref(1, window1%istate, icv, 1, sc_lib, topologies, block, xyz)
        natom_res = size(xyz, 2)

        allocate(all_res_xyz(3, natom_res, freq, nsample))
        
        do itimes = 1, nsample
            do jres = 1, freq
                call get_res_ref(itimes, window1%istate, icv, jres, sc_lib, topologies, block, xyz)
                all_res_xyz(:,:,jres,itimes) = xyz
            end do
        end do

        min_score = huge(1.0d0)
        best_j = 0
        best_jres = 0

        do itimes = 1, nsample
            do jres = 1, freq
                sum_score = 0.0d0
                do itime2 = 1, nsample
                    do ires = 1, freq
                        if (itime2 == itimes .and. ires == jres) cycle
                        call kabsch_rotation(natom_res, all_res_xyz(:,:,jres,itimes), all_res_xyz(:,:,ires,itime2), rot, rmsd)
                        sum_score = sum_score + rmsd
                    end do
                end do
                mean_score = sum_score / dble(nsample * freq - 1)
                if (mean_score < min_score) then
                    min_score = mean_score
                    best_j = itimes
                    best_jres = jres
                end if
            end do
        end do

        res_j = sc_lib%residue(icv)%res_ids(best_jres)
        xyz = all_res_xyz(:,:,best_jres,best_j)

        v11 = 0.0d0; v12 = 0.0d0; v21 = 0.0d0; v22 = 0.0d0
        f11 = .false.; f12 = .false.; f21 = .false.; f22 = .false.

        if (sc_lib%residue(icv)%has_sym) then
            call direction_vecs_sym(sc_lib%residue(icv)%natom_head, xyz, sc_lib%residue(icv)%cv_id, local_cv, sc_lib%residue(icv)%order)
            d1_ref = sc_lib%residue(icv)%d1_ref
            d2_ref = sc_lib%residue(icv)%d2_ref
        else
            do j = 1, sc_lib%residue(icv)%natom_head
                iatom = sc_lib%residue(icv)%map(j) + topologies%state(window1%istate)%topo%acid_start_list(res_j) - 1
                atom_name = topologies%state(window1%istate)%topo%atom_id(iatom)
                if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then; v11 = xyz(:,j); f11 = .true.; end if
                if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then; v12 = xyz(:,j); f12 = .true.; end if
                if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then; v21 = xyz(:,j); f21 = .true.; end if
                if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then; v22 = xyz(:,j); f22 = .true.; end if
            end do
            d1_ref = v12 - v11
            d2_ref = v22 - v21
            call direction_vecs_ref_inv(sc_lib%residue(icv)%natom_head, xyz, d1_ref, d2_ref, local_cv)
        end if

        sc_lib%residue(icv)%ref_xyz = xyz
        sc_lib%residue(icv)%d1_ref = d1_ref
        sc_lib%residue(icv)%d2_ref = d2_ref
        sc_lib%residue(icv)%cv = local_cv

        deallocate(all_res_xyz)
    end do
end subroutine search_sidechain_references


subroutine search_user_cv_references(user_lib, topologies, window1, nsample, block)
    use user_defined_cv_module, only: user_cv_library_type
    implicit none
    type(user_cv_library_type), intent(inout) :: user_lib
    type(topology_set_type), intent(in) :: topologies
    type(time_window), intent(in) :: window1
    integer, intent(in) :: nsample
    type(raw_data_block), intent(in) :: block

    real(8), allocatable :: xyz(:,:), all_user_xyz(:,:,:), local_cv(:,:)
    real(8) :: d1_orig(3), d2_orig(3), rot(3,3), rmsd, sum_score, mean_score, min_score
    integer :: icv, natom_cv, itimes, itime2, best_j

    do icv = 1, user_lib%n_cvs
        if (.not. user_lib%cvs(icv)%active) cycle
        natom_cv = user_lib%cvs(icv)%natom
        if (natom_cv < 1) cycle

        write(0,*) 'RMSD reference search for User CV: ', trim(user_lib%cvs(icv)%cv_name)

        allocate(all_user_xyz(3, natom_cv, nsample))
        do itimes = 1, nsample
            call get_user_cv_ref(itimes, window1%istate, icv, user_lib, topologies, block, xyz)
            all_user_xyz(:,:,itimes) = xyz
        end do

        min_score = huge(1.0d0)
        best_j = 0

        do itimes = 1, nsample
            sum_score = 0.0d0
            do itime2 = 1, nsample
                if (itime2 == itimes) cycle
                call kabsch_rotation(natom_cv, all_user_xyz(:,:,itimes), all_user_xyz(:,:,itime2), rot, rmsd)
                sum_score = sum_score + rmsd
            end do
            mean_score = sum_score / dble(nsample - 1)
            if (mean_score < min_score) then
                min_score = mean_score
                best_j = itimes
            end if
        end do

        if (best_j == 0) best_j = 1
        xyz = all_user_xyz(:,:,best_j)

        ! Berechnung der Vektoren / Koeffizienten
        call direction_vecs_ref_inv_gyration(natom_cv, xyz, .true., d1_orig, d2_orig, local_cv)

        ! WICHTIG: Hier werden nun ref_xyz UND cv_values korrekt in die Library geschrieben!
        user_lib%cvs(icv)%ref_xyz = xyz
        if (allocated(user_lib%cvs(icv)%cv_values)) deallocate(user_lib%cvs(icv)%cv_values)
        allocate(user_lib%cvs(icv)%cv_values(size(local_cv, 1), size(local_cv, 2)))
        user_lib%cvs(icv)%cv_values = local_cv

        deallocate(all_user_xyz)
        deallocate(local_cv)
    end do
end subroutine search_user_cv_references

subroutine ref_struc_search_rmsd( &
    project,topologies,sc_lib,dom_lib,window1, &
    nsample,method_window)

implicit none

type(time_window), intent(in) :: window1
type(project_type), intent(in) :: project
type(topology_set_type), intent(in) :: topologies
type(domain_library_type), intent(inout) :: dom_lib
type(sidechain_library_type), intent(inout) :: sc_lib

integer, intent(in) :: nsample
integer, intent(in) :: method_window

logical, allocatable :: window_mask1(:)

type(raw_data_block) :: block

real(8), allocatable :: xyz(:,:)
real(8), allocatable :: all_dom_xyz(:,:,:)
real(8), allocatable :: all_res_xyz(:,:,:,:)
real(8), allocatable :: local_cv(:,:),cv(:,:)

real(8) :: d1_orig(3),d2_orig(3)
real(8) :: d1_ref(3),d2_ref(3)

real(8) :: v11(3),v12(3)
real(8) :: v21(3),v22(3)

real(8) :: rot(3,3),rmsd
real(8) :: sum_score,mean_score,min_score

integer :: icv,ires,jres
integer :: itimes,itime2
integer :: i,j,iatom
integer :: natom_res,natom_dom
integer :: freq,ncount
integer :: best_j,best_jres
integer :: res_j,irep

character(len=4) :: atom_name

logical :: f11,f12,f21,f22


!=========================================================
! Daten laden
!=========================================================

if (nsample < 1) then
    error stop 'nsample must be >= 1'
end if

call fill_window_mask( &
    window1,window_mask1,method_window,nsample)

call load_raw_data_block( &
    project,window1,window_mask1,block)


!=========================================================
! 1. DOMAINS
!
! Kandidaten: alle nsample Frames einer Domain.
! Score: mittlere Kabsch-RMSD zu allen anderen Frames.
!=========================================================

do icv = 1, dom_lib%ndo

    write(0,*) 'RMSD search for domain: ', &
        trim(dom_lib%domain(icv)%domain_name)

    !-----------------------------------------------------
    ! Referenzstrukturen aller Frames laden
    !-----------------------------------------------------

    call get_dom_ref( &
        1,irep,window1%istate,icv, &
        dom_lib,topologies,project,block,xyz)

    natom_dom = size(xyz,2)

    allocate(all_dom_xyz(3,natom_dom,nsample))

    all_dom_xyz(:,:,1) = xyz

    do itimes = 2, nsample

        call get_dom_ref( &
            itimes,irep,window1%istate,icv, &
            dom_lib,topologies,project,block,xyz)

        if (size(xyz,2) /= natom_dom) then
            error stop 'Inconsistent domain atom count'
        end if

        all_dom_xyz(:,:,itimes) = xyz

    end do

    !-----------------------------------------------------
    ! Besten Frame bestimmen
    !-----------------------------------------------------

    min_score = huge(1.0d0)
    best_j = 0

    do itimes = 1, nsample

        sum_score = 0.0d0
        ncount = 0

        do itime2 = 1, nsample

            ! Selbstvergleich auslassen
            if (itime2 == itimes) cycle

            call kabsch_rotation( &
                natom_dom, &
                all_dom_xyz(:,:,itimes), &
                all_dom_xyz(:,:,itime2), &
                rot,rmsd)

            sum_score = sum_score + rmsd
            ncount = ncount + 1

        end do

        if (ncount > 0) then
            mean_score = sum_score / dble(ncount)
        else
            mean_score = 0.0d0
        end if

        if (mean_score < min_score) then
            min_score = mean_score
            best_j = itimes
        end if

    end do

    if (best_j == 0) then
        error stop 'No domain reference selected'
    end if

    !-----------------------------------------------------
    ! Ausgewählte Domain-Referenz speichern
    !-----------------------------------------------------

    xyz = all_dom_xyz(:,:,best_j)

    call direction_vecs_ref_inv_gyration( &
        size(xyz,2),xyz,.true.,d1_orig,d2_orig,cv)
    
    !dom_lib%domain(icv)%ref_frame = best_j
    !dom_lib%domain(icv)%ref_score = min_score

    dom_lib%domain(icv)%ref_xyz = xyz
    dom_lib%domain(icv)%d1_ref = d1_orig
    dom_lib%domain(icv)%d2_ref = d2_orig
    dom_lib%domain(icv)%cv = cv

    write(0,*) 'Selected domain reference: ', &
        trim(dom_lib%domain(icv)%domain_name), &
        ' frame=',best_j, &
        ' mean RMSD=',min_score

    deallocate(all_dom_xyz)

end do


!=========================================================
! 2. AMINOSÄURETYPEN
!
! Kandidaten: alle Kombinationen aus Frame und Residuum.
! Score: mittlere Kabsch-RMSD zu allen anderen Strukturen
!        desselben Aminosäuretyps.
!=========================================================

do icv = 1, sc_lib%ncv

    freq = sc_lib%residue(icv)%freq

    if (freq == 0) cycle

    write(0,*) 'RMSD search for residue type: ', &
        trim(sc_lib%residue(icv)%resname), &
        ' freq=',freq

    !-----------------------------------------------------
    ! Alle Strukturen dieses Typs laden
    !-----------------------------------------------------

    call get_res_ref( &
        1,window1%istate,icv,1, &
        sc_lib,topologies,block,xyz)

    natom_res = size(xyz,2)

    allocate(all_res_xyz(3,natom_res,freq,nsample))

    all_res_xyz(:,:,1,1) = xyz

    do itimes = 1, nsample

        do jres = 1, freq

            if (itimes == 1 .and. jres == 1) cycle

            call get_res_ref( &
                itimes,window1%istate,icv,jres, &
                sc_lib,topologies,block,xyz)

            if (size(xyz,1) /= 3 .or. &
                size(xyz,2) /= natom_res) then
                error stop 'Inconsistent sidechain atom count'
            end if

            all_res_xyz(:,:,jres,itimes) = xyz

        end do

    end do

    !-----------------------------------------------------
    ! Besten Kandidaten aus allen Frames und Residuen
    ! bestimmen
    !-----------------------------------------------------

    min_score = huge(1.0d0)
    best_j = 0
    best_jres = 0

    do itimes = 1, nsample

        do jres = 1, freq

            sum_score = 0.0d0
            ncount = 0

            do itime2 = 1, nsample

                do ires = 1, freq

                    ! Selbstvergleich auslassen
                    if (itime2 == itimes .and. &
                        ires == jres) cycle

                    call kabsch_rotation( &
                        natom_res, &
                        all_res_xyz(:,:,jres,itimes), &
                        all_res_xyz(:,:,ires,itime2), &
                        rot,rmsd)

                    sum_score = sum_score + rmsd
                    ncount = ncount + 1

                end do

            end do

            if (ncount > 0) then
                mean_score = sum_score / dble(ncount)
            else
                mean_score = 0.0d0
            end if

            if (mean_score < min_score) then
                min_score = mean_score
                best_j = itimes
                best_jres = jres
            end if

        end do

    end do

    if (best_j == 0 .or. best_jres == 0) then
        error stop 'No sidechain reference selected'
    end if

    !-----------------------------------------------------
    ! Gewählte Referenzstruktur übernehmen
    !-----------------------------------------------------

    res_j = sc_lib%residue(icv)%res_ids(best_jres)

    xyz = all_res_xyz(:,:,best_jres,best_j)

    ! Richtungsvektoren initialisieren
    v11 = 0.0d0
    v12 = 0.0d0
    v21 = 0.0d0
    v22 = 0.0d0

    f11 = .false.
    f12 = .false.
    f21 = .false.
    f22 = .false.

    !-----------------------------------------------------
    ! Symmetriebehandlung wie in der bisherigen Routine
    !-----------------------------------------------------

    if (sc_lib%residue(icv)%has_sym) then

        call direction_vecs_sym( &
            sc_lib%residue(icv)%natom_head, &
            xyz, &
            sc_lib%residue(icv)%cv_id, &
            local_cv, &
            sc_lib%residue(icv)%order)

        d1_ref = sc_lib%residue(icv)%d1_ref
        d2_ref = sc_lib%residue(icv)%d2_ref

    else

        ! Richtungsvektoren aus den Referenzatomen bestimmen
        do j = 1, sc_lib%residue(icv)%natom_head

            iatom = sc_lib%residue(icv)%map(j) + &
                topologies%state(window1%istate)%topo% &
                acid_start_list(res_j) - 1

            atom_name = &
                topologies%state(window1%istate)%topo%atom_id(iatom)

            if (trim(sc_lib%residue(icv)%dir1_atom_1) == &
                trim(atom_name)) then
                v11 = xyz(:,j)
                f11 = .true.
            end if

            if (trim(sc_lib%residue(icv)%dir1_atom_2) == &
                trim(atom_name)) then
                v12 = xyz(:,j)
                f12 = .true.
            end if

            if (trim(sc_lib%residue(icv)%dir2_atom_1) == &
                trim(atom_name)) then
                v21 = xyz(:,j)
                f21 = .true.
            end if

            if (trim(sc_lib%residue(icv)%dir2_atom_2) == &
                trim(atom_name)) then
                v22 = xyz(:,j)
                f22 = .true.
            end if

        end do

        if (.not.f11 .or. .not.f12 .or. &
            .not.f21 .or. .not.f22) then

            write(0,*) 'Missing reference atom: ', &
                trim(sc_lib%residue(icv)%resname), &
                ' frame=',best_j, &
                ' residue=',res_j

            error stop 'Reference direction atom missing'

        end if

        d1_ref = v12 - v11
        d2_ref = v22 - v21

        call direction_vecs_ref_inv( &
            sc_lib%residue(icv)%natom_head, &
            xyz,d1_ref,d2_ref,local_cv)

    end if

    !-----------------------------------------------------
    ! Referenzdaten und Score speichern
    !-----------------------------------------------------

    !sc_lib%residue(icv)%ref_frame = best_j
    !sc_lib%residue(icv)%ref_resid = res_j
    !sc_lib%residue(icv)%ref_score = min_score

    sc_lib%residue(icv)%ref_xyz = xyz
    sc_lib%residue(icv)%d1_ref = d1_ref
    sc_lib%residue(icv)%d2_ref = d2_ref
    sc_lib%residue(icv)%cv = local_cv

    write(0,*) 'Selected residue reference: ', &
        trim(sc_lib%residue(icv)%resname), &
        ' frame=',best_j, &
        ' residue=',res_j, &
        ' mean RMSD=',min_score

    deallocate(all_res_xyz)

end do

end subroutine ref_struc_search_rmsd

subroutine ref_struc_search_rmsd_user_old( &
    project, topologies, sc_lib, dom_lib, user_lib, window1, &
    nsample, method_window)

    use user_defined_cv_module, only: user_cv_library_type
    implicit none

    type(time_window), intent(in) :: window1
    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(domain_library_type), intent(inout) :: dom_lib
    type(sidechain_library_type), intent(inout) :: sc_lib
    type(user_cv_library_type), intent(inout) :: user_lib

    integer, intent(in) :: nsample
    integer, intent(in) :: method_window

    logical, allocatable :: window_mask1(:)

    type(raw_data_block) :: block

    real(8), allocatable :: xyz(:,:)
    real(8), allocatable :: all_dom_xyz(:,:,:)
    real(8), allocatable :: all_res_xyz(:,:,:,:)
    real(8), allocatable :: all_user_xyz(:,:,:) ! Für User-CVs
    real(8), allocatable :: local_cv(:,:),cv(:,:)

    real(8) :: d1_orig(3),d2_orig(3)
    real(8) :: d1_ref(3),d2_ref(3)

    real(8) :: v11(3),v12(3)
    real(8) :: v21(3),v22(3)

    real(8) :: rot(3,3),rmsd
    real(8) :: sum_score,mean_score,min_score

    integer :: icv,ires,jres
    integer :: itimes,itime2
    integer :: i,j,iatom,natom_cv
    integer :: natom_res,natom_dom
    integer :: freq,ncount
    integer :: best_j,best_jres
    integer :: res_j,irep

    character(len=4) :: atom_name
    logical :: f11,f12,f21,f22


    !=========================================================
    ! Daten laden
    !=========================================================

    if (nsample < 1) then
        error stop 'nsample must be >= 1'
    end if

    call fill_window_mask( &
        window1,window_mask1,method_window,nsample)

    call load_raw_data_block( &
        project,window1,window_mask1,block)


    !=========================================================
    ! 1. DOMAINS
    !=========================================================

    do icv = 1, dom_lib%ndo

        write(0,*) 'RMSD search for domain: ', &
            trim(dom_lib%domain(icv)%domain_name)

        call get_dom_ref( &
            1,irep,window1%istate,icv, &
            dom_lib,topologies,project,block,xyz)

        natom_dom = size(xyz,2)

        allocate(all_dom_xyz(3,natom_dom,nsample))
        all_dom_xyz(:,:,1) = xyz

        do itimes = 2, nsample
            call get_dom_ref( &
                itimes,irep,window1%istate,icv, &
                dom_lib,topologies,project,block,xyz)

            if (size(xyz,2) /= natom_dom) then
                error stop 'Inconsistent domain atom count'
            end if

            all_dom_xyz(:,:,itimes) = xyz
        end do

        min_score = huge(1.0d0)
        best_j = 0

        do itimes = 1, nsample
            sum_score = 0.0d0
            ncount = 0

            do itime2 = 1, nsample
                if (itime2 == itimes) cycle

                call kabsch_rotation( &
                    natom_dom, &
                    all_dom_xyz(:,:,itimes), &
                    all_dom_xyz(:,:,itime2), &
                    rot,rmsd)

                sum_score = sum_score + rmsd
                ncount = ncount + 1
            end do

            if (ncount > 0) then
                mean_score = sum_score / dble(ncount)
            else
                mean_score = 0.0d0
            end if

            if (mean_score < min_score) then
                min_score = mean_score
                best_j = itimes
            end if
        end do

        if (best_j == 0) then
            error stop 'No domain reference selected'
        end if

        xyz = all_dom_xyz(:,:,best_j)

        call direction_vecs_ref_inv_gyration( &
            size(xyz,2),xyz,.true.,d1_orig,d2_orig,cv)
        
        dom_lib%domain(icv)%ref_xyz = xyz
        dom_lib%domain(icv)%d1_ref = d1_orig
        dom_lib%domain(icv)%d2_ref = d2_orig
        dom_lib%domain(icv)%cv = cv

        write(0,*) 'Selected domain reference: ', &
            trim(dom_lib%domain(icv)%domain_name), &
            ' frame=',best_j, &
            ' mean RMSD=',min_score

        deallocate(all_dom_xyz)

    end do


    !=========================================================
    ! 2. AMINOSÄURETYPEN (Sidechains)
    !=========================================================

    do icv = 1, sc_lib%ncv

        freq = sc_lib%residue(icv)%freq
        if (freq == 0) cycle

        write(0,*) 'RMSD search for residue type: ', &
            trim(sc_lib%residue(icv)%resname), &
            ' freq=',freq

        call get_res_ref( &
            1,window1%istate,icv,1, &
            sc_lib,topologies,block,xyz)

        natom_res = size(xyz,2)

        allocate(all_res_xyz(3,natom_res,freq,nsample))
        all_res_xyz(:,:,1,1) = xyz

        do itimes = 1, nsample
            do jres = 1, freq
                if (itimes == 1 .and. jres == 1) cycle

                call get_res_ref( &
                    itimes,window1%istate,icv,jres, &
                    sc_lib,topologies,block,xyz)

                if (size(xyz,1) /= 3 .or. size(xyz,2) /= natom_res) then
                    error stop 'Inconsistent sidechain atom count'
                end if

                all_res_xyz(:,:,jres,itimes) = xyz
            end do
        end do

        min_score = huge(1.0d0)
        best_j = 0
        best_jres = 0

        do itimes = 1, nsample
            do jres = 1, freq
                sum_score = 0.0d0
                ncount = 0

                do itime2 = 1, nsample
                    do ires = 1, freq
                        if (itime2 == itimes .and. ires == jres) cycle

                        call kabsch_rotation( &
                            natom_res, &
                            all_res_xyz(:,:,jres,itimes), &
                            all_res_xyz(:,:,ires,itime2), &
                            rot,rmsd)

                        sum_score = sum_score + rmsd
                        ncount = ncount + 1
                    end do
                end do

                if (ncount > 0) then
                    mean_score = sum_score / dble(ncount)
                else
                    mean_score = 0.0d0
                end if

                if (mean_score < min_score) then
                    min_score = mean_score
                    best_j = itimes
                    best_jres = jres
                end if
            end do
        end do

        if (best_j == 0 .or. best_jres == 0) then
            error stop 'No sidechain reference selected'
        end if

        res_j = sc_lib%residue(icv)%res_ids(best_jres)
        xyz = all_res_xyz(:,:,best_jres,best_j)

        v11 = 0.0d0
        v12 = 0.0d0
        v21 = 0.0d0
        v22 = 0.0d0

        f11 = .false.
        f12 = .false.
        f21 = .false.
        f22 = .false.

        if (sc_lib%residue(icv)%has_sym) then
            call direction_vecs_sym( &
                sc_lib%residue(icv)%natom_head, &
                xyz, &
                sc_lib%residue(icv)%cv_id, &
                local_cv, &
                sc_lib%residue(icv)%order)

            d1_ref = sc_lib%residue(icv)%d1_ref
            d2_ref = sc_lib%residue(icv)%d2_ref
        else
            do j = 1, sc_lib%residue(icv)%natom_head
                iatom = sc_lib%residue(icv)%map(j) + &
                    topologies%state(window1%istate)%topo% &
                    acid_start_list(res_j) - 1

                atom_name = &
                    topologies%state(window1%istate)%topo%atom_id(iatom)

                if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then
                    v11 = xyz(:,j)
                    f11 = .true.
                end if
                if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then
                    v12 = xyz(:,j)
                    f12 = .true.
                end if
                if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then
                    v21 = xyz(:,j)
                    f21 = .true.
                end if
                if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then
                    v22 = xyz(:,j)
                    f22 = .true.
                end if
            end do

            if (.not.f11 .or. .not.f12 .or. .not.f21 .or. .not.f22) then
                write(0,*) 'Missing reference atom for sidechain: ', trim(sc_lib%residue(icv)%resname)
                error stop 'Reference direction atom missing'
            end if

            d1_ref = v12 - v11
            d2_ref = v22 - v21

            call direction_vecs_ref_inv( &
                sc_lib%residue(icv)%natom_head, &
                xyz,d1_ref,d2_ref,local_cv)
        end if

        sc_lib%residue(icv)%ref_xyz = xyz
        sc_lib%residue(icv)%d1_ref = d1_ref
        sc_lib%residue(icv)%d2_ref = d2_ref
        sc_lib%residue(icv)%cv = local_cv

        write(0,*) 'Selected residue reference: ', &
            trim(sc_lib%residue(icv)%resname), &
            ' frame=',best_j, &
            ' residue=',res_j, &
            ' mean RMSD=',min_score

        deallocate(all_res_xyz)

    end do


    !=========================================================
    ! 3. USER-DEFINED CVs (Strukturennähe via Kabsch-RMSD)
    !=========================================================

    do icv = 1, user_lib%n_cvs

        if (.not. user_lib%cvs(icv)%active) cycle

        natom_cv = user_lib%cvs(icv)%natom

        write(0,*) 'RMSD reference search for User CV: ', &
            trim(user_lib%cvs(icv)%cv_name)

        if (natom_cv < 1) cycle

        ! Koordinaten aller Frames für die Atome dieses User-CVs laden
        allocate(all_user_xyz(3, natom_cv, nsample))

        do itimes = 1, nsample
            call get_user_cv_ref(itimes, window1%istate, icv, user_lib, topologies, block, xyz)
            
            if (size(xyz, 2) /= natom_cv) then
                error stop 'Inconsistent user-defined CV atom count'
            end if
            all_user_xyz(:,:,itimes) = xyz
        end do

        ! RMSD- bzw. Kabsch-basierte Auswahl über die Frames
        min_score = huge(1.0d0)
        best_j = 0

        do itimes = 1, nsample
            sum_score = 0.0d0
            ncount = 0

            do itime2 = 1, nsample
                if (itime2 == itimes) cycle

                call kabsch_rotation( &
                    natom_cv, &
                    all_user_xyz(:,:,itimes), &
                    all_user_xyz(:,:,itime2), &
                    rot, rmsd)

                sum_score = sum_score + rmsd
                ncount = ncount + 1
            end do

            if (ncount > 0) then
                mean_score = sum_score / dble(ncount)
            else
                mean_score = 0.0d0
            end if

            if (mean_score < min_score) then
                min_score = mean_score
                best_j = itimes
            end if
        end do

        if (best_j == 0) best_j = 1

        xyz = all_user_xyz(:,:,best_j)

        call direction_vecs_ref_inv_gyration( &
            size(xyz,2), xyz, .true., d1_orig, d2_orig, local_cv)

        user_lib%cvs(icv)%ref_xyz = xyz

        write(0,*) 'Selected User CV reference: ', &
            trim(user_lib%cvs(icv)%cv_name), ' frame=', best_j, &
            ' mean RMSD=', min_score

        deallocate(all_user_xyz)

    end do

end subroutine ref_struc_search_rmsd_user_old

subroutine ref_struc_search(project,topologies,sc_lib,dom_lib,window1,nsample,method_window)
implicit none


type(time_window), intent(in) :: window1
type(project_type), intent(in) :: project
type(topology_set_type)  , intent(in)     :: topologies
type(domain_library_type), intent(inout) ::dom_lib
type(sidechain_library_type), intent(inout) :: sc_lib
integer, intent(in) :: nsample
integer, intent(in)  :: method_window

logical, allocatable :: window_mask1(:)
real(8), allocatable  :: refstrucs(:,:,:)
integer :: nref,iref,natom_all,refid,resid
real(8), allocatable :: xyz(:,:),score_vec(:),norm_score_vec(:,:)
integer, allocatable :: atoms(:)
real(8) ::d1_orig(3),d2_orig(3),d1_ref(3),d2_ref(3)
type(raw_data_block) :: block
real(8), allocatable :: local_cv(:,:),cv(:,:)
type(cv_plan)   :: plan
integer :: i,icv,irep,ires,jres,istate,itimes
type(lifted_trajectory) :: lifted
real(8) :: v11(3), v12(3)
real(8) :: v21(3), v22(3)
integer :: iatom
integer :: nhead,j,natom_res,k,res_start_atom,l
character(len=4) :: atom_name
logical :: f11,f12,f21,f22
integer :: res_i,res_j,best_j,jres2,ires2
real(8), allocatable :: scores(:,:,:)
real(8) ::min_score,cos_score
integer :: icv2, ncount
real(8) :: sum_score, mean_score
integer :: block_start,freq,idx,best_ires,best_jres


!stop

call fill_window_mask(window1,window_mask1,method_window,nsample) 
!call fill_window_mask(window2,window_mask2,method_window,nsample) 

call load_raw_data_block( &
    project, &
    window1, &
    window_mask1, &
    block)
nref=0
natom_all=0

do icv=1,sc_lib%ncv
 !if (sc_lib%residue(icv)%has_sym) cycle
        
   nref=nref+sc_lib%residue(icv)%freq*sc_lib%residue(icv)%freq
   do ires=1,sc_lib%residue(icv)%freq
      natom_all=natom_all+sc_lib%residue(icv)%natom_head
   end do
end do

nref=nref+dom_lib%ndo
do i=1,dom_lib%ndo
natom_all=natom_all+dom_lib%domain(i)%natom_keep
end do

allocate(scores(3,nsample,nref))
iref=0
write(0,*)'nref', nref
!pause

do itimes=1,nsample

nref=0

 do icv=1,dom_lib%ndo
 nref=nref+1
     
     call get_dom_ref(itimes,irep,window1%istate,icv,dom_lib,topologies,project,block,xyz)
     call get_atoms_domaine(istate,dom_lib,topologies,icv,atoms)
     call direction_vecs_ref_inv_gyration(size(xyz,2),xyz,.true.,d1_orig,d2_orig,cv)
     call add_domain_cv_pair(plan, window1%istate, dom_lib%domain(icv)%domain_name, &
                              'none', atoms, cv)

     !call ref_from_frame_single_domain(time,icv,irep,istate,lib,topologies,project,block)    
 end do
 
    do icv = 1, sc_lib%ncv
        !if (sc_lib%residue(icv)%has_sym) cycle
        
  
        
        do ires = 1, sc_lib%residue(icv)%freq
           
            resid = sc_lib%residue(icv)%res_ids(ires)
            res_i=sc_lib%residue(icv)%res_ids(ires)
            
            ! Zielstruktur / Atome durch ires bestimmen
            call get_atoms_sidechain( &
                window1%istate, sc_lib, topologies, icv, ires, atoms)
            
            do jres = 1, sc_lib%residue(icv)%freq
                nref=nref+1
                refid = sc_lib%residue(icv)%res_ids(jres)
                res_j=sc_lib%residue(icv)%res_ids(jres)
                
               
                ! Referenz für jres bestimmen
                call get_res_ref( &
                    itimes, window1%istate, icv, jres, sc_lib, topologies, &
                    block, xyz)
                write(0,*) shape(xyz),'initial shape'
               ! stop
                write(0,*) sc_lib%residue(icv)%resname, ires, jres

   ! Vektoren vor der Suche zurücksetzen
                v11 = 0.0_8
                v12 = 0.0_8
                v21 = 0.0_8
                v22 = 0.0_8
                
                f11=.false.
                f12=.false.
                f21=.false.
                f22=.false.
                
                do j = 1, sc_lib%residue(icv)%natom_head
                    iatom = sc_lib%residue(icv)%map(j) + &
                            topologies%state(window1%istate)%topo%acid_start_list(res_j) - 1

                   ! xyz(:,j) = block%xyz(:,iatom,itimes)
                    atom_name = topologies%state(window1%istate)%topo%atom_id(iatom)

              
                   !write(0,*) atom_name,'atomname',ires,jres,res_i,res_j
                   !write(0,*) sc_lib%residue(icv)%dir1_atom_1
                   !write(0,*) sc_lib%residue(icv)%dir1_atom_2
                   !write(0,*) sc_lib%residue(icv)%dir2_atom_1
                   !write(0,*) sc_lib%residue(icv)%dir2_atom_2
                    if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then
                        v11 = xyz(:,j)
                        f11=.true.
                    end if
                    if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then
                        v12 = xyz(:,j)
                        f12=.true.
                    end if
                    if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then
                        v21 = xyz(:,j)
                        f21=.true.
                    end if
                    if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then
                        v22 = xyz(:,j)
                        f22=.true.
                    end if
                end do
                
                         
         if (.not. f11 .or. .not. f12 .or. &
             .not. f21 .or. .not. f22) then
             write(0,*) sc_lib%residue(icv)%res_ids(:),'res_ids'
             
             write(0,*) 'Missing reference atom: ', &
                 trim(sc_lib%residue(icv)%resname), &
                 ' icv=', icv, ' jres=', jres
         
             write(0,*) 'dir1: ', &
                 trim(sc_lib%residue(icv)%dir1_atom_1), &
                 trim(sc_lib%residue(icv)%dir1_atom_2)
         
             write(0,*) 'dir2: ', &
                 trim(sc_lib%residue(icv)%dir2_atom_1), &
                 trim(sc_lib%residue(icv)%dir2_atom_2)
         
             write(0,*) 'found: ', f11, f12, f21, f22
         
             error stop 'Reference atom missing'
         end if


                d1_ref = v12 - v11
                d2_ref = v22 - v21
                write(0,*) d1_ref,d2_ref
                
                !stop
                if(sc_lib%residue(icv)%has_sym) then	
               call direction_vecs_sym(sc_lib%residue(icv)%natom_head&
               ,xyz,sc_lib%residue(icv)%cv_id,local_cv,sc_lib%residue(icv)%order)
                else
                ! CV direkt in die lokale Matrix schreiben lassen (2 x natom_head)
                call direction_vecs_ref_inv( &
                    sc_lib%residue(icv)%natom_head, &
                    xyz, &
                    d1_ref, &
                    d2_ref, &
                    local_cv)  ! <-- Direktes Übergeben der lokalen Matrix
              end if
                ! Die lokale Matrix direkt an add_cv_pair_sidechain übergeben
                call add_cv_pair_sidechain( &
                    plan, &
                    object_type=2, &
                    resid=resid, &
                    refid=refid, &
                    name=sc_lib%residue(icv)%resname, &
                    atoms=atoms, &
                    cv=local_cv, &  ! <-- Hier wird die korrekte, aktuelle Matrix übergeben
                    point_group='none', &
                    has_c2_sym=.false.)

            end do
        end do
        
    end do


  call lift_cv_plan(plan, block, lifted)

  
  call calculate_orthogonality_score(lifted,score_vec)
  call calculate_norm_score(lifted, norm_score_vec)
  

  write(0,*)  shape(score_vec),'shape score_vec'
  write(0,*)  shape(lifted%vec),'shape lifted%vec'
  write(0,*)  shape(scores),'shape scores'
 
  do i=1,size(score_vec)
  scores(1,itimes,i)=score_vec(i)
  scores(2:3,itimes,i)=norm_score_vec(:,i)
  write(0,*) i,score_vec(i),norm_score_vec(:,i),itimes
  end do
  
  ! pause
call destroy_cv_plan(plan)

end do


nref = 0

do i = 1, dom_lib%ndo
    nref = nref + 1

    min_score = huge(1.0d0)
    best_j = 0

    do j = 1, nsample
        cos_score = scores(1,j,nref)

        if (cos_score < min_score) then
            min_score = cos_score
            best_j = j
        end if
    end do

    if (best_j == 0) error stop 'No domain reference selected'

    ! Optimale Referenz erneut laden
    call get_dom_ref(best_j,irep,window1%istate,i, &
                     dom_lib,topologies,project,block,xyz)

    ! Referenzvektoren und CV berechnen
    call direction_vecs_ref_inv_gyration( &
         size(xyz,2),xyz,.true.,d1_orig,d2_orig,cv)

    ! Ausgewählte Referenz in dom_lib speichern
    !dom_lib%domain(i)%ref_frame = best_j
    !dom_lib%domain(i)%ref_score = min_score

    dom_lib%domain(i)%ref_xyz = xyz
    dom_lib%domain(i)%d1_ref = d1_orig
    dom_lib%domain(i)%d2_ref = d2_orig
    dom_lib%domain(i)%cv = cv

    write(0,*) 'Best domain reference: ', &
               trim(dom_lib%domain(i)%domain_name), &
               ' frame=',best_j, &
               ' score=',min_score
end do

nref = dom_lib%ndo

do icv = 1, sc_lib%ncv

    if (  sc_lib%residue(icv)%freq == 0) cycle

    freq = sc_lib%residue(icv)%freq
    block_start = nref + 1

    min_score = huge(1.0d0)
    best_j = 0
    best_jres = 0

    ! Optimales Referenzresiduum und optimalen Frame finden
    do jres = 1, freq
        do j = 1, nsample

            sum_score = 0.0d0

            do ires = 1, freq
                idx = block_start + (ires-1)*freq + (jres-1)
                sum_score = sum_score + scores(1,j,idx)
            end do

            if (sum_score < min_score) then
                min_score = sum_score
                best_j = j
                best_jres = jres
            end if

        end do
    end do

    if (best_j == 0 .or. best_jres == 0) then
        error stop 'No sidechain reference selected'
    end if

    ! Identität des gewählten Referenzresiduums
    res_j = sc_lib%residue(icv)%res_ids(best_jres)

    ! Ausgewählte Referenzstruktur erneut laden
    call get_res_ref(best_j,window1%istate,icv,best_jres, &
                     sc_lib,topologies,block,xyz)

    ! Richtungsvektoren aus genau dieser Referenz rekonstruieren
    v11 = 0.0d0
    v12 = 0.0d0
    v21 = 0.0d0
    v22 = 0.0d0

    f11 = .false.
    f12 = .false.
    f21 = .false.
    f22 = .false.

    do j = 1, sc_lib%residue(icv)%natom_head

        iatom = sc_lib%residue(icv)%map(j) + &
            topologies%state(window1%istate)%topo%acid_start_list(res_j) - 1

        atom_name = topologies%state(window1%istate)%topo%atom_id(iatom)

        if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then
            v11 = xyz(:,j)
            f11 = .true.
        end if

        if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then
            v12 = xyz(:,j)
            f12 = .true.
        end if

        if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then
            v21 = xyz(:,j)
            f21 = .true.
        end if

        if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then
            v22 = xyz(:,j)
            f22 = .true.
        end if

    end do

    if (.not.f11 .or. .not.f12 .or. .not.f21 .or. .not.f22) then
        error stop 'Reference direction atom not found'
    end if

    d1_ref = v12 - v11
    d2_ref = v22 - v21

    call direction_vecs_ref_inv( &
         sc_lib%residue(icv)%natom_head,xyz,d1_ref,d2_ref,local_cv)

    ! Ausgewählte Referenz und Score speichern
    !sc_lib%residue(icv)%ref_frame = best_j
    !sc_lib%residue(icv)%ref_resid = res_j
    !sc_lib%residue(icv)%ref_score = min_score

    sc_lib%residue(icv)%ref_xyz = xyz
    sc_lib%residue(icv)%d1_ref = d1_ref
    sc_lib%residue(icv)%d2_ref = d2_ref
    sc_lib%residue(icv)%cv = local_cv

    write(0,*) 'Best sidechain reference: ', &
               trim(sc_lib%residue(icv)%resname), &
               ' frame=',best_j, &
               ' residue=',res_j, &
               ' score=',min_score

    nref = nref + freq*freq

end do

end subroutine ref_struc_search

subroutine ref_struc_search_user(project, topologies, sc_lib, dom_lib, user_lib, window1, nsample, method_window)
    use user_defined_cv_module, only: user_cv_library_type
    implicit none

    type(time_window), intent(in) :: window1
    type(project_type), intent(in) :: project
    type(topology_set_type), intent(in) :: topologies
    type(domain_library_type), intent(inout) :: dom_lib
    type(sidechain_library_type), intent(inout) :: sc_lib
    type(user_cv_library_type), intent(inout) :: user_lib
    integer, intent(in) :: nsample
    integer, intent(in) :: method_window

    logical, allocatable :: window_mask1(:)
    real(8), allocatable :: refstrucs(:,:,:)
    integer :: nref, iref, natom_all, refid, resid
    real(8), allocatable :: xyz(:,:), score_vec(:), norm_score_vec(:,:)
    integer, allocatable :: atoms(:)
    real(8) :: d1_orig(3), d2_orig(3), d1_ref(3), d2_ref(3)
    type(raw_data_block) :: block
    real(8), allocatable :: local_cv(:,:), cv(:,:)
    type(cv_plan) :: plan
    integer :: i, icv, irep, ires, jres, istate, itimes
    type(lifted_trajectory) :: lifted
    real(8) :: v11(3), v12(3)
    real(8) :: v21(3), v22(3)
    integer :: iatom
    integer :: nhead, j, natom_res, k, res_start_atom, l
    character(len=4) :: atom_name
    logical :: f11, f12, f21, f22
    integer :: res_i, res_j, best_j, jres2, ires2
    real(8), allocatable :: scores(:,:,:)
    real(8) :: min_score, cos_score
    integer :: icv2, ncount
    real(8) :: sum_score, mean_score
    integer :: block_start, freq, idx, best_ires, best_jres
    integer :: natom_user
    character(len=32) :: ref_method

    call fill_window_mask(window1, window_mask1, method_window, nsample) 

    call load_raw_data_block( &
        project, &
        window1, &
        window_mask1, &
        block)

    nref = 0
    natom_all = 0

    ! Sidechain-Beitrag zu nref und natom_all
    do icv = 1, sc_lib%ncv
        nref = nref + sc_lib%residue(icv)%freq * sc_lib%residue(icv)%freq
        do ires = 1, sc_lib%residue(icv)%freq
            natom_all = natom_all + sc_lib%residue(icv)%natom_head
        end do
    end do

    ! Domain-Beitrag zu nref und natom_all
    nref = nref + dom_lib%ndo
    do i = 1, dom_lib%ndo
        natom_all = natom_all + dom_lib%domain(i)%natom_keep
    end do

    ! User-Defined CV Beitrag zu nref und natom_all
    do icv = 1, user_lib%n_cvs
        if (user_lib%cvs(icv)%active) then
            nref = nref + 1
            natom_all = natom_all + user_lib%cvs(icv)%natom
        end if
    end do

    allocate(scores(3, nsample, nref))
    iref = 0
    write(0,*) 'Total reference slots (nref):', nref

    do itimes = 1, nsample

        nref = 0

        ! 1. Domains zum Plan hinzufügen
        do icv = 1, dom_lib%ndo
            nref = nref + 1
            call get_dom_ref(itimes, irep, window1%istate, icv, dom_lib, topologies, project, block, xyz)
            call get_atoms_domaine(istate, dom_lib, topologies, icv, atoms)
            call direction_vecs_ref_inv_gyration(size(xyz,2), xyz, .true., d1_orig, d2_orig, cv)
            call add_domain_cv_pair(plan, window1%istate, dom_lib%domain(icv)%domain_name, &
                                    'none', atoms, cv)
        end do

        ! 2. Sidechains zum Plan hinzufügen
        do icv = 1, sc_lib%ncv
            do ires = 1, sc_lib%residue(icv)%freq
                resid = sc_lib%residue(icv)%res_ids(ires)
                res_i = sc_lib%residue(icv)%res_ids(ires)
                
                call get_atoms_sidechain( &
                    window1%istate, sc_lib, topologies, icv, ires, atoms)
                
                do jres = 1, sc_lib%residue(icv)%freq
                    nref = nref + 1
                    refid = sc_lib%residue(icv)%res_ids(jres)
                    res_j = sc_lib%residue(icv)%res_ids(jres)
                    
                    call get_res_ref( &
                        itimes, window1%istate, icv, jres, sc_lib, topologies, &
                        block, xyz)

                    v11 = 0.0_8; v12 = 0.0_8; v21 = 0.0_8; v22 = 0.0_8
                    f11 = .false.; f12 = .false.; f21 = .false.; f22 = .false.
                    
                    do j = 1, sc_lib%residue(icv)%natom_head
                        iatom = sc_lib%residue(icv)%map(j) + &
                                topologies%state(window1%istate)%topo%acid_start_list(res_j) - 1

                        atom_name = topologies%state(window1%istate)%topo%atom_id(iatom)

                        if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then
                            v11 = xyz(:,j); f11 = .true.
                        end if
                        if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then
                            v12 = xyz(:,j); f12 = .true.
                        end if
                        if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then
                            v21 = xyz(:,j); f21 = .true.
                        end if
                        if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then
                            v22 = xyz(:,j); f22 = .true.
                        end if
                    end do
                    
                    if (.not. f11 .or. .not. f12 .or. .not. f21 .or. .not. f22) then
                        error stop 'Reference atom missing for sidechain'
                    end if

                    d1_ref = v12 - v11
                    d2_ref = v22 - v21
                    
                    if (sc_lib%residue(icv)%has_sym) then    
                        call direction_vecs_sym(sc_lib%residue(icv)%natom_head, &
                                                xyz, sc_lib%residue(icv)%cv_id, local_cv, sc_lib%residue(icv)%order)
                    else
                        call direction_vecs_ref_inv( &
                            sc_lib%residue(icv)%natom_head, &
                            xyz, d1_ref, d2_ref, local_cv)
                    end if

                    call add_cv_pair_sidechain( &
                        plan, object_type=2, resid=resid, refid=refid, &
                        name=sc_lib%residue(icv)%resname, atoms=atoms, cv=local_cv, &
                        point_group='none', has_c2_sym=.false.)
                end do
            end do
        end do

        ! 3. User-Defined CVs zum Plan hinzufügen (mit Methoden-Fallunterscheidung)
        do icv = 1, user_lib%n_cvs
            if (.not. user_lib%cvs(icv)%active) cycle
            nref = nref + 1
            natom_user = user_lib%cvs(icv)%natom

            allocate(xyz(3, natom_user))
            call get_user_cv_ref(itimes, window1%istate, icv, user_lib, topologies, block, xyz)
            
            ref_method = trim(user_lib%cvs(icv)%reference_method)

            select case (ref_method)
            case ('atom_pair')
                v11 = xyz(:, user_lib%cvs(icv)%dir1_atom_1)
                v12 = xyz(:, user_lib%cvs(icv)%dir1_atom_2)
                v21 = xyz(:, user_lib%cvs(icv)%dir2_atom_1)
                v22 = xyz(:, user_lib%cvs(icv)%dir2_atom_2)
                
                d1_ref = v12 - v11
                d2_ref = v22 - v21

                allocate(local_cv(2, natom_user))
                call direction_vecs_ref_inv(natom_user, xyz, d1_ref, d2_ref, local_cv)
                
                call add_user_cv_pair(plan, window1%istate, user_lib%cvs(icv)%cv_name, user_lib%cvs(icv)%atom_indices, local_cv)
                deallocate(local_cv)

            case default ! 'gyration'
                call direction_vecs_ref_inv_gyration(natom_user, xyz, .true., d1_orig, d2_orig, cv)
                call add_user_cv_pair(plan, window1%istate, user_lib%cvs(icv)%cv_name, user_lib%cvs(icv)%atom_indices, cv)
                deallocate(cv)
            end select

            deallocate(xyz)
        end do

        ! LIFTING & SCORING
        call lift_cv_plan(plan, block, lifted)
        call calculate_orthogonality_score(lifted, score_vec)
        call calculate_norm_score(lifted, norm_score_vec)

        do i = 1, size(score_vec)
            scores(1, itimes, i) = score_vec(i)
            scores(2:3, itimes, i) = norm_score_vec(:, i)
        end do
        
        call destroy_cv_plan(plan)

    end do


    ! ==========================================================
    ! BESTE REFERENZEN AUSWÄHLEN UND SPEICHERN
    ! ==========================================================

    ! 1. Domains
    nref = 0
    do i = 1, dom_lib%ndo
        nref = nref + 1
        min_score = huge(1.0d0)
        best_j = 0

        do j = 1, nsample
            cos_score = scores(1, j, nref)
            if (cos_score < min_score) then
                min_score = cos_score
                best_j = j
            end if
        end do

        if (best_j == 0) error stop 'No domain reference selected'

        call get_dom_ref(best_j, irep, window1%istate, i, dom_lib, topologies, project, block, xyz)
        call direction_vecs_ref_inv_gyration(size(xyz,2), xyz, .true., d1_orig, d2_orig, cv)

        dom_lib%domain(i)%ref_xyz = xyz
        dom_lib%domain(i)%d1_ref = d1_orig
        dom_lib%domain(i)%d2_ref = d2_orig
        dom_lib%domain(i)%cv = cv
    end do

    ! 2. Sidechains
    nref = dom_lib%ndo
    do icv = 1, sc_lib%ncv
        if (sc_lib%residue(icv)%freq == 0) cycle

        freq = sc_lib%residue(icv)%freq
        block_start = nref + 1
        min_score = huge(1.0d0)
        best_j = 0
        best_jres = 0

        do jres = 1, freq
            do j = 1, nsample
                sum_score = 0.0d0
                do ires = 1, freq
                    idx = block_start + (ires-1)*freq + (jres-1)
                    sum_score = sum_score + scores(1, j, idx)
                end do

                if (sum_score < min_score) then
                    min_score = sum_score
                    best_j = j
                    best_jres = jres
                end if
            end do
        end do

        if (best_j == 0 .or. best_jres == 0) error stop 'No sidechain reference selected'

        res_j = sc_lib%residue(icv)%res_ids(best_jres)
        call get_res_ref(best_j, window1%istate, icv, best_jres, sc_lib, topologies, block, xyz)

        v11 = 0.0d0; v12 = 0.0d0; v21 = 0.0d0; v22 = 0.0d0
        f11 = .false.; f12 = .false.; f21 = .false.; f22 = .false.

        do j = 1, sc_lib%residue(icv)%natom_head
            iatom = sc_lib%residue(icv)%map(j) + topologies%state(window1%istate)%topo%acid_start_list(res_j) - 1
            atom_name = topologies%state(window1%istate)%topo%atom_id(iatom)

            if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then; v11 = xyz(:,j); f11 = .true.; end if
            if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then; v12 = xyz(:,j); f12 = .true.; end if
            if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then; v21 = xyz(:,j); f21 = .true.; end if
            if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then; v22 = xyz(:,j); f22 = .true.; end if
        end do

        if (.not. f11 .or. .not. f12 .or. .not. f21 .or. .not. f22) error stop 'Reference direction atom not found'

        d1_ref = v12 - v11
        d2_ref = v22 - v21
        call direction_vecs_ref_inv(sc_lib%residue(icv)%natom_head, xyz, d1_ref, d2_ref, local_cv)

        sc_lib%residue(icv)%ref_xyz = xyz
        sc_lib%residue(icv)%d1_ref = d1_ref
        sc_lib%residue(icv)%d2_ref = d2_ref
        sc_lib%residue(icv)%cv = local_cv

        nref = nref + freq*freq
    end do

    ! 3. User-Defined CVs (mit Methoden-Fallunterscheidung beim Speichern)
   do icv = 1, user_lib%n_cvs
        if (.not. user_lib%cvs(icv)%active) cycle
        nref = nref + 1
        min_score = huge(1.0d0)
        best_j = 0

        do j = 1, nsample
            cos_score = scores(1, j, nref)
            if (cos_score < min_score) then
                min_score = cos_score
                best_j = j
            end if
        end do

        if (best_j == 0) error stop 'No reference selected for user-defined CV'

        natom_user = user_lib%cvs(icv)%natom
        allocate(xyz(3, natom_user))
        call get_user_cv_ref(best_j, window1%istate, icv, user_lib, topologies, block, xyz)

        ref_method = trim(user_lib%cvs(icv)%reference_method)

        select case (ref_method)
        case ('atom_pair')
            v11 = xyz(:, user_lib%cvs(icv)%dir1_atom_1)
            v12 = xyz(:, user_lib%cvs(icv)%dir1_atom_2)
            v21 = xyz(:, user_lib%cvs(icv)%dir2_atom_1)
            v22 = xyz(:, user_lib%cvs(icv)%dir2_atom_2)
            
            d1_ref = v12 - v11
            d2_ref = v22 - v21

            allocate(local_cv(2, natom_user))
            call direction_vecs_ref_inv(natom_user, xyz, d1_ref, d2_ref, local_cv)
            
            user_lib%cvs(icv)%ref_xyz = xyz
            if (allocated(user_lib%cvs(icv)%cv_values)) deallocate(user_lib%cvs(icv)%cv_values)
            allocate(user_lib%cvs(icv)%cv_values(2, natom_user))
            user_lib%cvs(icv)%cv_values = local_cv
            deallocate(local_cv)

        case default ! 'gyration'
            call direction_vecs_ref_inv_gyration(natom_user, xyz, .true., d1_orig, d2_orig, cv)
            user_lib%cvs(icv)%ref_xyz = xyz
            if (allocated(user_lib%cvs(icv)%cv_values)) deallocate(user_lib%cvs(icv)%cv_values)
            allocate(user_lib%cvs(icv)%cv_values(size(cv, 1), size(cv, 2)))
            user_lib%cvs(icv)%cv_values = cv
            deallocate(cv)
        end select

        deallocate(xyz)
    end do

end subroutine ref_struc_search_user


subroutine calculate_orthogonality_score(lifted, score_vec)

    type(lifted_trajectory), intent(in) :: lifted
    real(8), allocatable, intent(out) :: score_vec(:)

    real(8) :: v1(3), v2(3)
    real(8) :: norm1, norm2, cosine
    integer :: i, j, npairs, nframes

    npairs  = size(lifted%vec, 2) / 2
    nframes = size(lifted%vec, 3)

    if (mod(size(lifted%vec, 2), 2) /= 0) then
        error stop 'Odd number of vectors'
    end if

    allocate(score_vec(npairs))
    score_vec = 0.0_8

    do i = 1, npairs
        do j = 1, nframes

            v1 = lifted%vec(:, 2*i-1, j)
            v2 = lifted%vec(:, 2*i,   j)

            norm1 = sqrt(dot_product(v1, v1))
            norm2 = sqrt(dot_product(v2, v2))

            if (norm1 > 1.0e-12_8 .and. &
                norm2 > 1.0e-12_8) then

                cosine = dot_product(v1, v2) / (norm1 * norm2)

                score_vec(i) = score_vec(i) + cosine**2

            end if

        end do

        score_vec(i) = score_vec(i) / real(nframes, 8)

    end do

end subroutine calculate_orthogonality_score

subroutine calculate_norm_score(lifted, score_vec)

    type(lifted_trajectory), intent(in) :: lifted
    real(8), allocatable, intent(out) :: score_vec(:,:)

    real(8) :: v(3)
    real(8) :: norm, mean_norm
    integer :: i, j, npairs, nframes, k

    npairs  = size(lifted%vec, 2) / 2
    nframes = size(lifted%vec, 3)

    if (mod(size(lifted%vec, 2), 2) /= 0) then
        error stop 'Odd number of vectors'
    end if

    if (nframes < 1) then
        error stop 'No frames in lifted trajectory'
    end if

    ! Zwei Varianzwerte pro Vektorpaar
    allocate(score_vec(2, npairs))
    score_vec = 0.0_8

    do i = 1, npairs
        do k = 1, 2

            mean_norm = 0.0_8

            ! Mittelwert der Vektorlänge
            do j = 1, nframes
                v = lifted%vec(:, 2*i-2+k, j)
                mean_norm = mean_norm + sqrt(dot_product(v, v))
            end do

            mean_norm = mean_norm / real(nframes, 8)

            ! Varianz der Vektorlänge
            do j = 1, nframes
                v = lifted%vec(:, 2*i-2+k, j)
                norm = sqrt(dot_product(v, v))

                score_vec(k,i) = score_vec(k,i) + &
                    (norm - mean_norm)**2
            end do

            score_vec(k,i) = score_vec(k,i) / real(nframes, 8)

        end do
    end do

end subroutine calculate_norm_score

function cross_product(x1,x2) result(c)

implicit none

real(8), intent(in) :: x1(3), x2(3)
real(8) :: c(3)

c(1) = x1(2)*x2(3) - x1(3)*x2(2)
c(2) = x1(3)*x2(1) - x1(1)*x2(3)
c(3) = x1(1)*x2(2) - x1(2)*x2(1)

end function cross_product



subroutine debug_prolin_mapping(project,topologies,sc_lib,dom_lib,window1,window2,nsample,method_window,method_score)
implicit none


type(time_window), intent(in) :: window1,window2
type(project_type), intent(in) :: project
type(topology_set_type)  , intent(in)     :: topologies
type(domain_library_type), intent(inout) ::dom_lib
type(sidechain_library_type), intent(inout) :: sc_lib
integer, intent(in) :: nsample
integer, intent(in)  :: method_window,method_score 

logical, allocatable :: window_mask1(:),window_mask2(:)
real(8), allocatable  :: refstrucs(:,:,:)
integer :: nref,iref,natom_all,refid,resid
real(8), allocatable :: xyz(:,:)
integer, allocatable :: atoms(:)
real(8) ::d1_orig(3),d2_orig(3),d1_ref(3),d2_ref(3)
type(raw_data_block) :: block
real(8), allocatable :: local_cv(:,:),cv(:,:)
type(cv_plan)   :: plan
integer :: i,icv,irep,ires,jres,istate,itimes
type(lifted_trajectory) :: lifted
real(8) :: v11(3), v12(3)
real(8) :: v21(3), v22(3)
integer :: iatom
integer :: nhead,j,natom_res,k,res_start_atom,l
character(len=4) :: atom_name
logical :: f11,f12,f21,f22
integer :: res_i,res_j


icv=17
        !if (sc_lib%residue(icv)%has_sym) cycle
        
        ! Lokale 2 x natom_head Matrix für die CVs definieren, 
        ! damit sie pro jres-Schleife nicht überschrieben werden:
        ! (Falls moderneres Fortran / allocatable:)
        ! real(8), allocatable :: local_cv(:,:)
        ! allocate(local_cv(2, sc_lib%residue(icv)%natom_head))
        
        do ires = 1, sc_lib%residue(icv)%freq
            resid = sc_lib%residue(icv)%res_ids(ires)
            res_i=sc_lib%residue(icv)%res_ids(ires)
            
            ! Zielstruktur / Atome durch ires bestimmen
            call get_atoms_sidechain( &
                window1%istate, sc_lib, topologies, icv, ires, atoms)
            
            do jres = 1, sc_lib%residue(icv)%freq
                refid = sc_lib%residue(icv)%res_ids(jres)
                res_j=sc_lib%residue(icv)%res_ids(jres)
                
                ! Referenz für jres bestimmen
               !call get_res_ref( &
               !    itimes, window1%istate, icv, jres, sc_lib, topologies, &
               !    block, xyz)
               !
                write(0,*) sc_lib%residue(icv)%resname, ires, jres

   ! Vektoren vor der Suche zurücksetzen
                v11 = 0.0_8
                v12 = 0.0_8
                v21 = 0.0_8
                v22 = 0.0_8
                
                f11=.false.
                f12=.false.
                f21=.false.
                f22=.false.
                
                do j = 1, sc_lib%residue(icv)%natom_head
                    iatom = sc_lib%residue(icv)%map(j) + &
                            topologies%state(window1%istate)%topo%acid_start_list(res_j) - 1

                    !sc_lib%residue(icv)%ref_xyz(:,j) = block%xyz(:,iatom,itimes)
                    atom_name = topologies%state(window1%istate)%topo%atom_id(iatom)

                   write(0,*) 'map',sc_lib%residue(icv)%map
                   write(0,*) atom_name,'atomname',ires,jres,res_i,res_j,iatom
                   write(0,*) sc_lib%residue(icv)%dir1_atom_1
                   write(0,*) sc_lib%residue(icv)%dir1_atom_2
                   write(0,*) sc_lib%residue(icv)%dir2_atom_1
                   write(0,*) sc_lib%residue(icv)%dir2_atom_2
                    if (trim(sc_lib%residue(icv)%dir1_atom_1) == trim(atom_name)) then
                       ! v11 = xyz(:,j)
                        f11=.true.
                    end if
                    if (trim(sc_lib%residue(icv)%dir1_atom_2) == trim(atom_name)) then
                        !v12 = xyz(:,j)
                        f12=.true.
                    end if
                    if (trim(sc_lib%residue(icv)%dir2_atom_1) == trim(atom_name)) then
                        !v21 = xyz(:,j)
                        f21=.true.
                    end if
                    if (trim(sc_lib%residue(icv)%dir2_atom_2) == trim(atom_name)) then
                        !v22 = xyz(:,j)
                        f22=.true.
                    end if
                end do
                
                
      if (.not. f11 .or. .not. f12 .or. &
    .not. f21 .or. .not. f22) then
    write(0,*) sc_lib%residue(icv)%res_ids(:),'res_ids'
    
    write(0,*) 'Missing reference atom: ', &
        trim(sc_lib%residue(icv)%resname), &
        ' icv=', icv, ' jres=', jres

    write(0,*) 'dir1: ', &
        trim(sc_lib%residue(icv)%dir1_atom_1), &
        trim(sc_lib%residue(icv)%dir1_atom_2)

    write(0,*) 'dir2: ', &
        trim(sc_lib%residue(icv)%dir2_atom_1), &
        trim(sc_lib%residue(icv)%dir2_atom_2)

    write(0,*) 'found: ', f11, f12, f21, f22

    error stop 'Reference atom missing'
end if

write(0,*) 'state=', window1%istate, &
           ' res_j=', res_j, &
           ' start=', topologies%state(window1%istate)%topo%acid_start_list(res_j)

write(0,*) 'map=', sc_lib%residue(icv)%map(:)

write(0,*) 'topology atoms=', &
    topologies%state(window1%istate)%topo%atom_id( &
        topologies%state(window1%istate)%topo%acid_start_list(res_j): &
        topologies%state(window1%istate)%topo%acid_start_list(res_j)+ &
        topologies%state(window1%istate)%topo%acid_size(res_j)-1)

                d1_ref = v12 - v11
                d2_ref = v22 - v21
                write(0,*) d1_ref,d2_ref
                ! CV direkt in die lokale Matrix schreiben lassen (2 x natom_head)
                call direction_vecs_ref_inv( &
                    sc_lib%residue(icv)%natom_head, &
                    xyz, &
                    d1_ref, &
                    d2_ref, &
                    local_cv)  ! <-- Direktes Übergeben der lokalen Matrix
              
                ! Die lokale Matrix direkt an add_cv_pair_sidechain übergeben
                call add_cv_pair_sidechain( &
                    plan, &
                    object_type=2, &
                    resid=resid, &
                    refid=refid, &
                    name=sc_lib%residue(icv)%resname, &
                    atoms=atoms, &
                    cv=local_cv, &  ! <-- Hier wird die korrekte, aktuelle Matrix übergeben
                    point_group='none', &
                    has_c2_sym=.false.)

            end do
        end do
        
        if (trim(sc_lib%residue(icv)%resname) == 'PRO') then
    write(0,*) 'PRO target residue index = ', resid
    write(0,*) 'PRO target residue start = ', &
         topologies%state(window1%istate)%topo%acid_start_list(resid)

    do j = 1, sc_lib%residue(icv)%natom_head
        iatom = sc_lib%residue(icv)%map(j) + &
                topologies%state(window1%istate)%topo%acid_start_list(resid) - 1

        write(0,*) trim(sc_lib%residue(icv)%atom_names(j)), &
                   ' map = ', sc_lib%residue(icv)%map(j), &
                   ' target atom = ', trim( topologies%state(window1%istate)%topo%atom_id(iatom))
    end do
end if
        
        ! deallocate(local_cv) ! falls oben allokiert





end subroutine debug_prolin_mapping



subroutine fill_window_mask(window,window_mask,method,nsample)
implicit none
integer, intent(in) ::method,nsample
type(time_window), intent(in) :: window
integer :: ntimes,stride,nfound,i
real    :: rn
logical, allocatable, intent(out) :: window_mask(:)



ntimes=window%end_time_step-window%start_time_step+1
write(0,*) window%start_time_step,window%end_time_step
allocate(window_mask(ntimes))

window_mask=.false.
if(method==1) then

stride=ntimes/nsample
  if(stride<1 .or. stride>ntimes) then
  	  write(0,*) 'error in fill_window_mask',ntimes,nsample,stride
  	  stop
  end if
  do i=1,ntimes,stride
    window_mask(i)=.true.
  end do
return  
end if	

if(method==2) then

nfound=0
do while(nfound<nsample)
call random_number(rn)
i=floor(rn*dble(ntimes))+1
if(window_mask(i)) then
	cycle
else
	nfound=nfound+1
end if

end do
return
end if	

write(0,*) 'unknowen method to fill window mask',method
stop

end subroutine



end module ref_struc_selection_module