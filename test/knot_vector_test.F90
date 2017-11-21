module knot_vector_test
    use fruit
    implicit none

contains

    subroutine test_knot1()
        use knot_vector, only : FillOpenKnot
        use test_utils, only : compare_matrices
        implicit none
        integer(kind = 4) :: n, p
        real (kind = 8), allocatable, dimension(:) :: U
        real (kind = 8), allocatable, dimension(:) :: res

        n = 0
        p = 0
        allocate(res(n + p + 2))
        call FillOpenKnot(res, n, p)
        U = (/ 1.d0, 1.d0 /)
        call compare_matrices(U, res)
        deallocate(res)

    end subroutine test_knot1



    subroutine test_knot2()
        use knot_vector, only: CountSpans
        implicit none
        integer :: res
        real (kind = 8), DIMENSION(8) :: A

        A = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
        res = CountSpans(1, 1, A)
        call assert_equals(1, res)
    end subroutine test_knot2


    subroutine test_knot3()
        use knot_vector, only: PrepareKnot
        implicit none
        integer(kind = 4) :: n, p
        real (kind = 8), allocatable, dimension(:) :: U
        integer(kind = 4) :: nelem


    end subroutine test_knot3




end module knot_vector_test
