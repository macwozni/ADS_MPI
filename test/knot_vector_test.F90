module knot_vector_test
    use fruit
    implicit none

contains

    subroutine test_knot1()
        use knot_vector, only: CountSpans
        integer :: result

        real (kind = 8), DIMENSION(8) :: A = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
        result = CountSpans(1, 1, A)
        call assert_equals(-1, result)
    end subroutine test_knot1

end module knot_vector_test
