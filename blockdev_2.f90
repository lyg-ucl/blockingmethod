!==========================================================================
! Program: BLOCKDEV                                              
! Version: BLOCKDEV_2.0
! Copyright (C) Yunguo Li @ UCL                              
! Email :: yunguo.li@ucl.ac.uk                               
! Function :: calc. statistical uncertainty by using the blocking method
!             Ref.: H. Flyvbjerg, J. Chem. Phys. 91, 461 (1989)
! Usage ::  blockdev.x datafile num_of_column                             
! Input ::  column-wise datafile                          
!==========================================================================
PROGRAM BLOCKDEV
!==========================================================================
IMPLICIT  NONE

INTEGER, PARAMETER :: ap = selected_real_kind(15,300)
INTEGER :: n,m,ntype,nprim,i,j,k,IOstatus=0
REAL(ap) :: mxddev
REAL(ap), DIMENSION(:,:), ALLOCATABLE :: array1,array2,array3,array4,delta
REAL(ap),DIMENSION(:),ALLOCATABLE :: c,average,stddev,ddelta
CHARACTER(len=255) :: arg1,arg2

 CALL get_command_argument(1,arg1)
 CALL get_command_argument(2,arg2)
 arg2=trim(arg2)
 READ(arg2,'(I3)') ntype
 OPEN(unit = 15, file = trim(arg1), status='old')

 n=0
 DO WHILE (IOstatus .eq. 0)
   READ(15,*,IOSTAT=IOstatus)   
   n = n + 1
 END DO
 n = n -1
 write(*,'(A16,I10)') "number of type: ",ntype
 write(*,'(A16,I10)') "number of data: ",n 

 ALLOCATE(array1(ntype,n))
 ALLOCATE(array2(ntype,n/2))
 ALLOCATE(average(ntype))
 ALLOCATE(c(ntype))
 ALLOCATE(delta(ntype,2))
 ALLOCATE(ddelta(ntype))
 ALLOCATE(stddev(ntype))

 average=0.0
 REWIND(15)
 DO i=1,n
   READ(15,*) array1(1:ntype,i)
   average(1:ntype) = average(1:ntype) + array1(1:ntype,i)
 END DO
 average = average / n
 WRITE(*,*) "Average :             ", average 

 DO i=1,n
   stddev(1:ntype)=stddev(1:ntype)+(array1(1:ntype,i)-average(1:ntype))**2
 END DO
 stddev(1:ntype)=SQRT(stddev(1:ntype)/n)
 WRITE(*,*) "Standard deviation:   ", stddev(1:ntype)

 OPEN(16, file ='temp.dat', status='replace')
 DO i=1, n / 2
   array2(1:ntype,i) = 0.5*(array1(1:ntype,2*i-1) + array1(1:ntype,2*i))
   c(1:ntype) = c(1:ntype) + (array2(1:ntype,i) - average(1:ntype))**2
   WRITE(16,*) array2(1:ntype,i)
 END DO
 c = 2*c / n
 CLOSE(16)
 delta(:,1) = SQRT(c(:)/(n/2-1))
! deldev(:)=   c/(n/2-1) * (1/SQRT(2*(REAL(n/2)-1))) )
 delta(:,2) = delta(:,1)
 j=1
 DO WHILE (m.GT.4)
   OPEN(16, file ='temp.dat', status='old')
   m=n/2
   ALLOCATE(array3(ntype,m))
   ALLOCATE(array4(ntype,m/2))
   DO i=1, m
     READ(16,*) array3(1:ntype,i)
   ENDDO
   CLOSE(16)
   OPEN(16, file ='temp.dat', status='replace')
   DO i=1, m / 2
     array4(1:ntype,i) = 0.5*(array3(1:ntype,2*i-1) + array3(1:ntype,2*i))
     c(1:ntype) = c(1:ntype) + (array4(1:ntype,i) - average(1:ntype))**2
     WRITE(16,*) array4(1:ntype,i) 
   ENDDO
   c = 2*c / m
   CLOSE(16)
   delta(:,2) = SQRT(c(:)/(m-1))
   ddelta(:) = (delta(:,2)-delta(:,1))/delta(:,1)
   mxddev = MAXVAL(ddelta)
   WRITE(*,*) "Block size :         ", j
!   WRITE(*,*) "Stastical error :    ", delta(:,2),ddelta(:)
!   WRITE(*,*) "Stastical error (%): ", delta(:,2)/average
   DEALLOCATE(array3)
   DEALLOCATE(array4)
   n=m ; j=j+1; delta(:,1) = delta(:,2)
 END DO

 WRITE(*,*) "Stastical error :    ", delta(:,2)
 WRITE(*,*) "Stastical error (%): ", delta(:,2)/average 

 DEALLOCATE(array1)
 DEALLOCATE(array2)
 DEALLOCATE(average)
 DEALLOCATE(c)
 DEALLOCATE(delta)
 DEALLOCATE(ddelta)
 DEALLOCATE(stddev)

END PROGRAM BLOCKDEV


