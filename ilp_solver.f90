module ilp_solver
    implicit none
contains

    subroutine solve_ilp(n, atomic_masses, target_mass, max_atoms, solutions, deviations, num_solutions, tolerance)
        integer, intent(in) :: n, max_atoms(n)
        real(8), intent(in) :: atomic_masses(n), target_mass, tolerance
        integer, intent(out) :: solutions(1000, n)
        real(8), intent(out) :: deviations(1000)
        integer, intent(out) :: num_solutions
        integer :: temp_solution(n)

        num_solutions = 0
        call branch_and_bound(n, atomic_masses, target_mass, max_atoms, temp_solution, 1, solutions, deviations, num_solutions, tolerance)
    end subroutine solve_ilp

    recursive subroutine branch_and_bound(n, atomic_masses, target_mass, max_atoms, temp_solution, level, solutions, deviations, num_solutions, tolerance)
        integer, intent(in) :: n, max_atoms(n), level
        real(8), intent(in) :: atomic_masses(n), target_mass, tolerance
        integer, intent(inout) :: temp_solution(n), solutions(1000, n)
        real(8), intent(inout) :: deviations(1000)
        integer, intent(inout) :: num_solutions
        integer :: i, j
        real(8) :: current_mass, current_deviation

        if (level > n) then
            ! Compute the total molecular mass
            current_mass = 0.0
            do j = 1, n
                current_mass = current_mass + temp_solution(j) * atomic_masses(j)
            end do
            current_deviation = abs(current_mass - target_mass)

            ! ✅ Store only solutions within the tolerance
            if (current_deviation <= tolerance .and. num_solutions < 1000) then
                num_solutions = num_solutions + 1
                solutions(num_solutions, :) = temp_solution
                deviations(num_solutions) = current_deviation
            end if
            return
        end if

        ! Iterate over possible atom counts
        do i = 0, max_atoms(level)
            temp_solution(level) = i
            call branch_and_bound(n, atomic_masses, target_mass, max_atoms, temp_solution, level + 1, solutions, deviations, num_solutions, tolerance)
        end do
    end subroutine branch_and_bound

end module ilp_solver
