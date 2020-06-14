program omp_par_do
 use omp_lib
 implicit none
 character(len=40)    :: filename           !Storing the file name
 character(len=20)    :: nam                !Collecting the file name within the file
 integer(8)           :: N
 real*8,dimension(1:700000)::rx,ry,rz,rx1,rx2,ry1,ry2,rz1,rz2,result1                !Number of coordinates within  the file
 !real*8,dimension(30000:30000)::coa1
 real(8), allocatable :: coa(:,:),coa1(:,:),coa2(:,:)         !Array containing xyz coordinates
 character(len=6), allocatable :: atom(:,:) !Array containing the atomic make up
 integer(8)          ::l, i,j,k,start,count1,count2,count3,count4,count5,count6,count7,count8,start1,start2,N1,     a          !Do loop parameters
 real*8::summ,summ1,sum1,potential,charge,q1,dist,vale
   integer :: num_threads = 16
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 N1=76800
 charge = 130 !! interms of electron
 q1 = real(charge) / real(N1)
     open (10, file="52161.xyz", status='OLD')
     open (1, file="52x16-1.xyz", status='unknown')
     open (2, file="52x16-2.xyz", status='unknown')
 
     do i = 1,N1
         read(10,*)rx(i),ry(i),rz(i)
     end do
 
     do i = 1,N1
         write(1,*) rx(i), ry(i), rz(i)
         write(2,*) rx(i)+160.0D0, ry(i), rz(i)
     end do
 
     close (1)
     close (2)
 
     open (1, file="52x16-1.xyz", status='unknown')
     open (2, file="52x16-2.xyz", status='unknown')
 
     do i = 1,N1
         read(1,*) rx1(i), ry1(i), rz1(i)
         read(2,*) rx2(i), ry2(i), rz2(i)
     end do
 
     do k = 1,100
         open (101, file="charge130.txt",access = 'append', status='old')
         summ = 0.0
         dist = 0.0
         potential = 0.0
         !$ call omp_set_num_threads(num_threads)
         !$omp parallel default(none) &
         !$omp private(i,j,summ,potential,dist)     &     
         ! !$omp private (Eall,Ms2,NotNeigh0) &
         !$omp shared(N1,q1,result1,rx1,ry1,rz1,rx2,ry2,rz2)
         !$omp do
         do i =1,N1
             do j = 1,N1
                 !  summ = summ + (rx1(i)+rx2(j))
                 dist = dsqrt ((rx1(i)-rx2(j))**2.0D0 + (ry1(i)-ry2(j))**2.0D0 + (rz1(i)-rz2(j))**2.0D0)
                 potential = (((q1*1.6d0)**2)*9.0d0)/(dist*80.0D0)
                 summ = summ +potential*(1.44D0)*10.0D0
                 ! potential = (((-23.0d0*2.3d0)*(1.78D0)))*((1.7d0/(dist))**6.0D0)
                 ! potential = (((23.0d0*2.3d0)*(1.0d0/(3.14)**2.0)))/(dist**6.0D0)
                 ! summ = summ +potential
             end do
         result1(i) = summ
         summ = 0.0
         end do
         !$omp end do
         !$omp end parallel
 
         summ = 0.0
         ! do i =1,N1
             summ = sum(result1(1:N1))
         ! end do
         result1=0.0
         write(101,*)rx2(1)-rx1(1)-160.0d0,summ!,result1(1:2)
         do l =  1, N1
             rx2(l) = rx2(l) + 2.0d0
         end do
         close(101)
     end do
     end program omp_par_do
