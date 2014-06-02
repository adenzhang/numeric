integer function ludcmp(a, indx, d)
    ! given MxN matrix a, replace it by LU decomposition of a rowwise permutation of itself.
    ! indx is an vector of length N, recording the row permutation effected by partial pivating.
    ! d is output as +-1 depending on whether the number of row interchanges was even or odd.
    ! return 0 if success
    use nrtype
    use nrutil, only: assert_eq, imaxloc, nrerror, outerprod, swap
    implicit none
    real(SP), dimension(:,:), intent(inout) :: a
    integer(i4b), dimension(:), intent(out) :: indx
    real(sp), intent(out) :: d

    real(sp), dimension(size(a,1)) :: vv     ! vv stores the implicit scaling of each row
    real(sp), parameter :: tiny = 1.02-10_sp ! a small number
    INTEGER(I4B) :: j,n,imax

    ludcmp = 0
    write(*,*) 'size a:', size(a,1), 'x', size(a,2)
    n = assert_eq(size(a,1), size(a,2), size(indx), 'ludcmp')
    d = 1.0                                  ! no row interchanges yet.
    vv = maxval(abs(a),dim=2)                ! loop over rows to get implicit scaling info.
    if( any(vv == 0.0) ) then
        call nrerror('singular matrix in ludcmp')
        ludcmp = -1
    end if
    vv = 1.0_sp/vv                           ! save the scaling.
    do j=1,n
        imax = j-1 + imaxloc(vv(j:n)*abs(a(j:N, j))) ! find the pivot row
        if( j /= imax) then                  ! change rows
             call swap(a(imax,:), a(j,:))
             d = -d                          ! change the parity of d
             vv(imax) = vv(j)                ! change scale factor
        end if
        indx(j) = imax
        if( a(j,j) == 0.0 ) then
            a(j,j) = tiny    ! singular matrix
            ludcmp = -2
        end if
        a(j+1:n, j) = a(j+1:n,n)/a(j,j)      ! divide by privot element.
        a(j+1:n, j+1:n) = a(j+1:n, j+1:n) - outerprod(a(j+1:n, j), a(j, j+1:n)) ! reduce remaining sub matrix
    end do
end function

integer function lubksb(a,indx,b)
    !Solves the set of N linear equations A * X = B.Herethe N × N matrix a is input, not
    !as the original matrix A, but rather as its LU decomposition, determined by the routine
    !ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
    !input as the right-hand-side vector B,alsooflength N, and returns with the solution vector
    !X. a and indx are not modiﬁed by this routine and can be left in place for successive calls
    !with diﬀerent right-hand sides b. This routine takes into account the possibility that b will
    !begin with many zero elements, so it is eﬃcient for use in matrix inversion.
    USE nrtype; USE nrutil, ONLY : assert_eq
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER(I4B) :: i,n,ii,ll
    REAL(SP) :: summ

    lubksb = 0
    n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii=0 !When ii is set to a positive value, it will become the in-
!    dex of the ﬁrst nonvanishing element of b.Wenowdo
!    the forward substitution, equation (2.3.6). The only new
!    wrinkle is to unscramble the permutation as we go.
    do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii /= 0) then
            summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ /= 0.0) then
            ii=i !A nonzero element was encountered, so from now on we will
        !    have to do the dot product above.
        end if
        b(i)=summ
    end do
    do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
END function lubksb


