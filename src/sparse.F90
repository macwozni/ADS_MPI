module sparse
    
implicit none

type sparse_matrix_entry
    integer (kind=4) :: x
    integer (kind=4) :: y
    real    (kind=8) :: val
    type (sparse_matrix_entry), pointer :: next
end type sparse_matrix_entry
    
type sparse_matrix_line
    integer (kind=4) :: x
    type (sparse_matrix_line), pointer :: next
    type (sparse_matrix_entry), pointer :: first
end type sparse_matrix_line

type sparse_matrix
    integer (kind=4) :: x
    integer (kind=4) :: y
    integer (kind=8) :: total_entries
    type (sparse_matrix_line), pointer :: first
end type sparse_matrix



contains



subroutine initialize_sparse(x,y,matrix)
    implicit none
    integer (kind=4), intent(in) :: x
    integer (kind=4), intent(in) :: y
    type (sparse_matrix), pointer, intent(out) :: matrix
    
    allocate(matrix)
    matrix%x = x
    matrix%y = y
    matrix%total_entries = 0
    matrix%first => NULL()
end subroutine initialize_sparse


subroutine find(matrix,x,y,entr)
    implicit none
    type(sparse_matrix), pointer, intent(inout) :: matrix
    integer (kind=4), intent(in) :: x
    integer (kind=4), intent(in) :: y
    type(sparse_matrix_entry), pointer, intent(out)    :: entr
    type(sparse_matrix_line), pointer  :: line
    type(sparse_matrix_entry), pointer :: tmp
    type(sparse_matrix_entry), pointer :: tmp2
    
    call find_line(matrix,x,line)
    
    if(.NOT. associated(line%first)) then
        allocate(line%first)
        line%first%x = x
        line%first%y = y
        line%first%next => NULL()
        line%first%val = 0.d0
        entr => line%first
        return
    end if
    
    tmp => line%first
    
    do while (associated(tmp%next))
        if (tmp%next%y .EQ. y) then
            entr => tmp%next
            return
        end if
        if (tmp%next%y .GT. y) then
            tmp2 => tmp%next
            allocate(tmp%next)
            tmp%next%x = x
            tmp%next%y = y
            tmp%next%next => tmp2
            tmp%val = 0.d0
            entr => tmp%next
            return
        end if
        tmp => tmp%next
    end do
end subroutine find


subroutine find_line(matrix,x,line)
    implicit none
    type(sparse_matrix), intent(inout) :: matrix
    integer (kind=4), intent(in) :: x
    type(sparse_matrix_line), pointer, intent(out) :: line
    type(sparse_matrix_line), pointer :: tmp
    type(sparse_matrix_line), pointer :: tmp2
    
    if (.NOT.associated(matrix%first)) then
        allocate(matrix%first)
        matrix%first%x = x
        matrix%first%next => NULL()
        matrix%first%first => NULL()
        line => matrix%first
        return
    endif
    
    tmp => matrix%first
    
    do while (associated(tmp%next))
        if (tmp%next%x .EQ. x) then
            line => tmp%next
            return
        end if
        if (tmp%next%x .GT. x) then
            tmp2 => tmp%next
            allocate(tmp%next)
            tmp%next%x = x
            tmp%next%next => tmp2
            tmp%next%first => NULL()
            line => tmp%next
            return
        end if
        tmp => tmp%next
    end do
    
    allocate(tmp%next)
    tmp%next%x = x
    tmp%next%next => NULL()
    tmp%next%first => NULL()
    line => tmp%next
end subroutine find_line



subroutine add(matrix,x,y,val)
    implicit none
    type(sparse_matrix), pointer, intent(inout) :: matrix
    integer (kind=4), intent(in) :: x
    integer (kind=4), intent(in) :: y
    real (kind=8), intent(in) :: val
    type(sparse_matrix_entry), pointer :: entr
    
    call find(matrix,x,y,entr)
    entr%val = entr%val + val
end subroutine add


subroutine clear_matrix(matrix)
    implicit none
    type(sparse_matrix), pointer, intent(inout) :: matrix
    type(sparse_matrix_line), pointer :: line
    
    if(.NOT. associated(matrix)) return
    line => matrix%first
    deallocate(matrix)
    call clear_line(line)
end subroutine clear_matrix



recursive subroutine clear_line(line)
    implicit none
    type(sparse_matrix_line), pointer, intent(inout) :: line
    type(sparse_matrix_line), pointer :: tmp
    
    if (.NOT. associated(line)) return
    tmp => line%next
    call clear_entry(line%first)
    deallocate(line)
    call clear_line(tmp)
end subroutine clear_line



recursive subroutine clear_entry(entr)
    implicit none
    type(sparse_matrix_entry), pointer, intent(inout) :: entr
    type(sparse_matrix_entry), pointer :: tmp
    
    if(.NOT. associated(entr)) return
    tmp => entr%next
    deallocate(entr)
    call clear_entry(tmp)
end subroutine clear_entry


subroutine to_dense_matrix(matrix, dmatrix)
    implicit none
    type(sparse_matrix), intent(in) :: matrix
    real (kind=8), dimension(matrix%x,matrix%y), intent(out) :: dmatrix
    type(sparse_matrix_line), pointer :: line
    type(sparse_matrix_entry), pointer :: entr
    
    dmatrix = 0.d0
    line => matrix%first
    
    do while (associated(line))
        entr => line%first
        do while(associated(entr))
            dmatrix(entr%x, entr%y) = entr%val
            entr => entr%next
        end do
        line => line%next
    end do
end subroutine to_dense_matrix

end module sparse
    
    
