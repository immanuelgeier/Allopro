module te_plot_module
    use iso_fortran_env, only: error_unit, real64
    use te_io_module, only: te_plot_data_type
    use direction_vectors_analyse_plan_module, only: cv_plan
    implicit none
    private

    public :: write_te_summary_table
    public :: fill_pdb_occ_and_beta_sidechains

contains

subroutine write_te_summary_table(filename, data)
    character(len=*), intent(in) :: filename
    type(te_plot_data_type), intent(in) :: data
    integer :: unit, ios, i

    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'Cannot write TE summary table: ',trim(filename)
        error stop 'TE summary output failed'
    end if
    write(unit,'(A)') '# cv row_sum column_sum antisymmetric_row_sum antisymmetric_column_sum'
    do i = 1, data%matrix_size
        write(unit,'(I0,4(1X,ES24.16E3))',iostat=ios) i,data%row_sum(i),data%column_sum(i), &
            data%antisymmetric_row_sum(i),data%antisymmetric_column_sum(i)
        if (ios /= 0) then
            close(unit)
            error stop 'Error writing TE summary table'
        end if
    end do
    close(unit,iostat=ios)
    if (ios /= 0) error stop 'Error closing TE summary table'
end subroutine write_te_summary_table

subroutine fill_pdb_occ_and_beta_sidechains(cv_p, cv_values, cv_offset, scale, occupancy, beta)
    type(cv_plan), intent(in) :: cv_p
    real(real64), intent(in) :: cv_values(:)
    integer, intent(in) :: cv_offset
    real(real64), intent(in) :: scale
    real(4), intent(out) :: occupancy(:), beta(:)
    integer :: i, cv_index, atom_index

    occupancy = 0.0_4
    beta = 0.0_4
    do i = 1, size(cv_p%term)
        if (cv_p%term(i)%target <= cv_offset) cycle
        cv_index = cv_p%term(i)%target
        atom_index = cv_p%term(i)%atom
        if (2*cv_index > size(cv_values)) error stop 'CV plan target exceeds the TE matrix size'
        if (atom_index < 1 .or. atom_index > size(occupancy)) &
            error stop 'CV plan atom index exceeds the selected topology'
        occupancy(atom_index) = real(cv_values(2*cv_index-1)*scale,4)
        beta(atom_index) = real(cv_values(2*cv_index)*scale,4)
    end do
end subroutine fill_pdb_occ_and_beta_sidechains

end module te_plot_module
