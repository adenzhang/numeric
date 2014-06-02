program integrationarray
implicit none

integer  :: i, j
integer, parameter :: n=10
real(8), dimension(n)  :: x, y
real(8), dimension(n,n) :: val
real(8)  :: integ

integ = 0.0D0

do j=1,10
  y(j) = dble(j)*1.0D0
  x(j) = dble(j) * 2.5D-1
enddo
call subr(x,y,val,n)
do j=1,n
  do i=1,n
     integ = integ + val(i,j)
   enddo
enddo

write (*,'(a,e20.10e3)')' Integration value = ',integ

end program integrationarray

subroutine subr(a,b,c,n)
  implicit none
  integer :: n, i, j
  real(8) :: a(n),b(n),c(n,n)
  write(*,*) ' size a:', size(a,1)
  do j=1,n
     do i=1,n
     c(i,j) = sin(a(i)+b(j))
     enddo
  enddo
  return
end subroutine subr