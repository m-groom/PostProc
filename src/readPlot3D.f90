! To compile: gfortran -O3 -ffree-line-length-none -cpp -fstack-arrays -fno-protect-parens -shared -fPIC readPlot3D.f90 -o libreadplot3d.so
! With debugging: gfortran -O3 -ffree-line-length-none -cpp -fstack-arrays -fno-protect-parens -shared -fPIC -g -fbacktrace readPlot3D.f90 -o libreadplot3d.so
! cp libreadplot3d.so ../lib

! Subroutine for reading per processor binary Plot 3D output for a 3D domain
subroutine readPerProcPlot3D(sTimeStep, iNVar, iNX, iNY, iNZ, rX, rY, rZ, rQ)
   implicit none
   ! Name of input file
   character(32), intent(in) :: sTimeStep
   ! Number of variables
   integer, intent(in) :: iNVar
   ! Number of cells in x, y, z
   integer, intent(in) :: iNX, iNY, iNZ
   ! Coordinates
   real(4), allocatable, dimension(:,:,:), intent(out) :: rX, rY, rZ
   ! Solution array
   real(4), allocatable, dimension(:,:,:,:), intent(out) :: rQ
   ! Number of domains
   integer :: iNPR, iProc, iBlocks
   ! Start and stop of each domain
   integer, allocatable, dimension(:) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
   integer :: iE, jE, kE
   ! File reference
   integer :: iFileReference, data_unit
   ! Integers for loops
   integer :: i, j, k, n
   ! File name
   character(64) :: sFileName, sPath, sProc

   sPath="./data"

   ! Open prempi.dat file
   iFileReference = 300
   open(iFileReference,file=trim(adjustl(sPath))//"/"//"prempi.dat",action="read")
   read(iFileReference,*) iNPR
   allocate(iStart(iNPR),iEnd(iNPR),jStart(iNPR),jEnd(iNPR),kStart(iNPR),kEnd(iNPR))
   do n = 1, iNPR
      read(iFileReference,*) iProc
      read(iFileReference,*) kStart(iProc+1), kEnd(iProc+1)
      read(iFileReference,*) jStart(iProc+1), jEnd(iProc+1)
      read(iFileReference,*) iStart(iProc+1), iEnd(iProc+1)
   end do
   close(iFileReference)

   ! Allocate arrays
   allocate(rX(iNX+1,iNY+1,iNZ+1), rY(iNX+1,iNY+1,iNZ+1), rZ(iNX+1,iNY+1,iNZ+1))
   allocate(rQ(iNX+1,iNY+1,iNZ+1,iNVar))

   ! Read grid file
   do n = 1, iNPR
      ! Get file name for current domain
      write(sProc,fmt="(i8)") n-1
      sFileName = trim(adjustl(sPath))//"/"//"out_0.00000000/0.00000000"//"."//trim(adjustl(sProc))//".g"
      ! Open file
      open(newunit=data_unit,form="unformatted",file=sFileName,action="read")
      read(data_unit) iBlocks
      read(data_unit) iE, jE, kE
      if (iE.ne.(iEnd(n)-1 - iStart(n)).or.jE.ne.(jEnd(n)-1 - jStart(n)).or.kE.ne.(kEnd(n)-1 - kStart(n))) then
         write(*,*) "Error: Grid size mismatch."
         write(*,*) "Expected: ", iEnd(n)-1 - iStart(n), jEnd(n)-1 - jStart(n), kEnd(n)-1 - kStart(n)
         write(*,*) "Found: ", iE, jE, kE
         stop
      end if
      read(data_unit)(((rX(i,j,k),i=iStart(n),iEnd(n)-2),j=jStart(n),jEnd(n)-2),k=kStart(n),kEnd(n)-2),&
                     (((rY(i,j,k),i=iStart(n),iEnd(n)-2),j=jStart(n),jEnd(n)-2),k=kStart(n),kEnd(n)-2),&
                     (((rZ(i,j,k),i=iStart(n),iEnd(n)-2),j=jStart(n),jEnd(n)-2),k=kStart(n),kEnd(n)-2)
      close(data_unit)
   end do

   rQ = 0.0

end subroutine readPerProcPlot3D


