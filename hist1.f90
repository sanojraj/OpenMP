program hist
implicit none 
integer(8)           :: N, N1, step, i, j, 
real*8,dimension(0:600000)::rx,ry,ry1, coun1                !Number of coordinates within the file

N1=249   ! No of particle in a cluster
step=2   ! bining size
    
    open (10, file="test222", status='OLD')
    do i = 1,N1
        read(10,*)rx(i),ry(i)
    end do

    do i = 0,1000
        coun1(i) = 0.0
    end do

    do i =0,180,step
        do j = 1,N1
            if ((int(ry(j)) .gt. i) .and. (int(ry(j)) .le. i+step)) then
                coun1(i) = coun1(i)+1.0
            end if
      
        end do
        print*,i,coun1(i)
        coun1 = 0.0
     end do

end program

