module test_utils
    use fruit
    implicit none

contains

    subroutine compare_matrices(U, R)
        use fruit
        implicit none
        integer (kind = 4) :: len1, len2
        real (kind = 8), allocatable, dimension(:) :: U
        real (kind = 8), allocatable, dimension(:) :: R

        len1 = size(U)
        len2 = size(R)
        call assert_equals(len1, len2)
    end subroutine compare_matrices
end module test_utils