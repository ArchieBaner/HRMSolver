program hrms_solver
    use ilp_solver
    implicit none

    integer, parameter :: max_elements = 10, max_solutions = 1000
    integer :: num_elements, max_atoms(max_elements), num_solutions, i, j, db_size
    real(8) :: target_mass, atomic_masses(max_elements), deviations(max_solutions), tolerance
    character(len=2) :: elements(max_elements), db_elements(118)
    real(8) :: db_masses(118)
    integer :: solutions(max_solutions, max_elements)

    ! Read element data from elements.dat
    call load_element_database("elements.dat", db_elements, db_masses, db_size)

    ! Get user input
    print*, "Enter the target molar mass:"
    read*, target_mass

    print*, "Enter number of elements in the molecule (max 10):"
    read*, num_elements

    print*, "Enter element symbols:"
    do i = 1, num_elements
        read*, elements(i)
        atomic_masses(i) = find_atomic_mass(elements(i), db_elements, db_masses, db_size)
    end do

    print*, "Enter max count of each element in the molecule:"
    do i = 1, num_elements
        read*, max_atoms(i)
    end do

    print*, "Enter acceptable deviation tolerance:"
    read*, tolerance

    ! Solve the ILP problem and find all valid formulas
    call solve_ilp(num_elements, atomic_masses, target_mass, max_atoms, solutions, deviations, num_solutions, tolerance)

    ! Print results
    if (num_solutions == 0) then
        print*, "No valid formulas found!"
    else
        print*, "Valid molecular formulas:"
        do i = 1, num_solutions
            do j = 1, num_elements
                if (solutions(i, j) > 0) then
                    write(*, '(A, I2)', advance='no') trim(elements(j)), solutions(i, j)
                end if
            end do
            print*, "  Deviation:", deviations(i)
        end do
    end if

contains

    ! Load element data from file
    subroutine load_element_database(filename, elements, masses, size)
        character(len=*), intent(in) :: filename
        character(len=2), intent(out) :: elements(118)
        real(8), intent(out) :: masses(118)
        integer, intent(out) :: size
        integer :: i, io_status
        open(10, file=filename, status="old", action="read", iostat=io_status)
        if (io_status /= 0) then
            print*, "Error: Unable to open elements.dat"
            stop
        end if

        size = 0
        do i = 1, 118
            read(10, *, iostat=io_status) elements(i), masses(i)
            if (io_status /= 0) exit
            size = size + 1
        end do
        close(10)
    end subroutine load_element_database

    ! Find atomic mass for a given element
    function find_atomic_mass(element, elements, masses, size) result(mass)
        character(len=2), intent(in) :: element, elements(118)
        real(8), intent(in) :: masses(118)
        integer, intent(in) :: size
        real(8) :: mass
        integer :: i

        mass = -1.0
        do i = 1, size
            if (trim(elements(i)) == trim(element)) then
                mass = masses(i)
                return
            end if
        end do
        print*, "Error: Element ", element, " not found in database."
        stop
    end function find_atomic_mass

end program hrms_solver
