program app
    use nrtype
    use nrutil, only: get_time
    integer, parameter :: N = 3000
    real :: A(N,N), alu(N,N), b(N), x(N)
    integer :: iptv(N), info, d(N)

      real start_time, finish_time
     integer count_0, count_1, count_rate, count_max

    call random_number(A)
    call random_number(b)

!   test lapack
    write(*,*) '-- using lapack'
    alu = a;
    x = b;
    start_time = get_time()
    write (6, '(8x, 1a, 1f16.6)') 'begin time:  ', start_time

    call sgetrf(N,N,alu,N,iptv,info)
    if( info/=0 ) print *, "Error! sgetrf"
    call sgetrs('N', N, 1, alu, N, iptv, x, N, info)
    if( info/=0 ) print *, "Error! sgetrs"

    finish_time =  get_time()
    write (6, '(8x, 2(1a, 1f16.6))') 'end time:    ', finish_time, ' timing: ',finish_time - start_time

    call sgemv('N',N,N, 1.0, A,N,x,1,-1.0,b,1)
    write(*,*) "error:", snrm2(N, 1, b)

!   test nr
    write(*,*) '-- using nr'
    alu = a;
    x = b;
    start_time = get_time()
    write (6, '(8x, 1a, 1f16.6)') 'begin time:  ', start_time

    info = ludcmp(alu, iptv, d)
    if( info/=0 ) print *, "Error! ludcmp"
    info = lubksb(alu, iptv, x)
    if( info/=0 ) print *, "Error! lubksb"

    finish_time =  get_time()
    write (6, '(8x, 2(1a, 1f16.6))') 'end time:    ', finish_time, ' timing: ',finish_time - start_time

    call sgemv('N',N,N, 1.0, A,N,x,1,-1.0,b,1)
    write(*,*) "error:", snrm2(N, 1, b)

end program app




