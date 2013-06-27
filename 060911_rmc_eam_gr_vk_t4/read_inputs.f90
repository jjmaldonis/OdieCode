!************************************************************
!******A module is defined here to load experiment data******
!*******This subroutine is written by FY on 12/18/2008*******
!************************************************************
MODULE ReadInputs

CONTAINS

!It reads input file names, temperature, maximum atom movement distance
!cutoff distance, used_data_sets, weight factors calculating chi square
!structure factor pairs derived from electron scattering, neutron scattering
!and X-ray scattering
!FEM data V(k), k and V_err triplets
!Rotation number and fem algorithm
!Pixel distance
!status

!   change log
!   gr_e_err is added - jwh 04/25/2009
!   modified to print what's read from the fem file right away - JWH 05/08/2009
!   The order of reading FEM angles are changed to nphi, npsi, ntheta = JWH 09/04/09




!If this subroutine used a module or modules, do not include implicit none line

SUBROUTINE read_inputs(param_filename,temperature, max_move, cutoff_r, used_data_sets, weights, gr_e, r_e, gr_e_err, &
gr_n, r_n, gr_x, r_x, V, k, V_err, V_background, ntheta, nphi, npsi, scale_fac, Q, fem_algorithm, pixel_distance, total_steps, &
rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

!param_filename=input file name containing initilizing parameters,such as temperature
!temperature=beginning temperature for RMC
!max_move=maximum allowable movement distance
!cutoff_r=cut off distance for different pairs
!used_data_sets=A logical array determine whether electron, neutron, X-ray or FEM data are available or not
!gr_e and r_e=electron G(r) data
!gr_n and r_n=neutron G(r) data
!gr_x and r_x=X-ray G(r) data
!V=FEM intensity variance
!k=scattering vector
!V_err=FEM measurement variance error
!ntheta, nphi, npsi=rotation number for calculating simulated V
!fem_algorithm=which algorithm is used to do variance calculation
!Q=resolution for FEM
!pixel_distance=pixel distance
!status=whether param_filename is opened with success or not

IMPLICIT NONE
!INTEGER, INTENT(IN) :: name_length
CHARACTER (LEN=*), INTENT(IN) :: param_filename  !Assume size array, be careful

REAL, INTENT(OUT) :: temperature
REAL, INTENT(OUT) :: max_move
REAL, INTENT(OUT), DIMENSION(:,:) :: cutoff_r !
LOGICAL, INTENT(OUT), DIMENSION(4) :: used_data_sets
REAL, INTENT(OUT), DIMENSION(4) :: weights
REAL, POINTER, DIMENSION(:) :: gr_e,r_e,gr_e_err
REAL, POINTER, DIMENSION(:) :: gr_n,r_n
REAL, POINTER, DIMENSION(:) :: gr_x,r_x
REAL, POINTER, DIMENSION(:) :: V, k, V_err, V_background
REAL, INTENT(OUT) :: rmin_e, rmax_e
REAL, INTENT(OUT) :: rmin_n, rmax_n
REAL, INTENT(OUT) :: rmin_x, rmax_x
INTEGER, INTENT(OUT) :: ntheta, nphi, npsi
REAL, INTENT(OUT) :: Q
REAL, INTENT(OUT) :: scale_fac
INTEGER, INTENT(OUT) ::fem_algorithm
REAL, INTENT(OUT) :: pixel_distance
INTEGER, INTENT(OUT) ::total_steps
INTEGER, INTENT(OUT) :: status2

!Variable declared in this subroutine, local variables
!comment1=comment in param_filename
CHARACTER (LEN=80) comment1  
CHARACTER (LEN=20) ScatteringFile ! The scattering file name, at most 20 characters
    ! Note: Electron, neutron, and X-ray data all use this name, but in order.
CHARACTER (LEN=20) FemFile ! The FEM data file name, at most 20 characters
INTEGER FileNameLength !The file name length in ScatteringFile or FemFile
INTEGER i, j
real Indicator_End !Indicator_End=-1 means the end of file reaches
INTEGER status1 !Indicate the status of opening file in this subroutine
    !Files for electron, neutron, X-ray scattering and FEM data
LOGICAL File_End !Indicate FileEnd has reached
INTEGER Num_Line !Number of lines in each data file except comment line
REAL, POINTER, DIMENSION(:) :: TempData !Temperature data
INTEGER Stat_Allocate1, Stat_Allocate2, Stat_Allocate3,Stat_Allocate4 !Allocate status, 0 means success

OPEN(20, FILE=param_filename,IOSTAT=status2, STATUS='OLD')
IF(status2 .NE. 0) THEN !Open fails
  PRINT *, 'Cannot open file with name: ', param_filename
  RETURN
ENDIF
READ(20, '(A80)') comment1 !Read comment and it is neglected in the input
!PRINT *, comment1
read(20, *) total_steps
read(20, *) temperature
read(20, *) max_move
read(20, * )
read(20, *) cutoff_r !This is a nelements by nelements matrix

!Fortran read each number in the file in column order

i=0

!********************************************************************
!*******************READ electron scattering file**********************
!********************************************************************
i=i+1
Num_Line=0

READ(20, '(A)') ScatteringFile
ScatteringFile = ADJUSTL(ScatteringFile)
!write(*,*)ScatteringFile
READ(20, *) weights(i)
READ(20, *) rmin_e, rmax_e



IF(ScatteringFile(1:1) .NE. '#') THEN !Electron scattering file exists
  used_data_sets(i) = .TRUE.

  FileNameLength=LEN_TRIM(ScatteringFile)
  File_End=.FALSE.

  !Read r_e and gr_e data
  !read r_e first, then gr_e
  !First line is comment
  !At first, count how many data pairs are in the file
  !-1 is chosen as the last point
  !PRINT *, 'ScatteringFile(1:FileNameLength) 1 is:  ', ScatteringFile(1:FileNameLength) !Debug
  OPEN(30,FILE=ScatteringFile(1:FileNameLength),IOSTAT=status1,STATUS='OLD') 
   IF(status1 .EQ. 0) THEN !Open succeeds
     READ(30, '(A80)') comment1
     DO WHILE( .NOT. File_End)
        READ(30, *) Indicator_End
        IF(ABS(Indicator_End+1.0) .LE. 1E-6) THEN !File end is reached
          !PRINT *, 'Indicator_End is: ', Indicator_End
          !PRINT *, 'Electron scattering file end is reached'
          EXIT !Exit this loop
        ELSE
          Num_Line = Num_Line + 1
          !READ(30, *) 
        ENDIF
     ENDDO

     REWIND(30) !Go to the beginning of the file

     READ(30, '(A80)') comment1
     ALLOCATE(TempData(3*Num_Line),STAT=Stat_Allocate1)
     IF(Stat_Allocate1 .EQ. 0) THEN
       READ(30, *) TempData
       ALLOCATE(r_e(Num_Line), STAT=Stat_Allocate2)
       ALLOCATE(gr_e(Num_Line), STAT=Stat_Allocate3)
          allocate(gr_e_err(num_line), stat=stat_allocate4)   !JWH 04/25/2009
       IF ((Stat_Allocate2 .EQ. 0) .AND. (Stat_Allocate3 .EQ. 0) .and. (Stat_Allocate4 .EQ. 0)  ) THEN
          r_e=TempData(1:3*Num_Line:3)
          gr_e=TempData(2:3*Num_Line:3)
          gr_e_err=tempdata(3:3*num_line:3)
       ELSE
         PRINT *, 'Electron scattering part 2 or 3 fails!'
         RETURN
       ENDIF !Allocate2 and 3
       
     ELSE
       PRINT *, 'Electron scattering part allocation fail'
       RETURN
     ENDIF

     DEALLOCATE(TempData)
    
   ELSE
    PRINT *, 'Open electron scattering file fails'
   ENDIF
  CLOSE(30)
ELSE 
  used_data_sets(i) = .FALSE.
ENDIF !Electron scattering file exists


!*******************READ electron scattering file**********************

!*******************************************************************
!*****************Read neutron scattering file**********************
!*******************************************************************

i=i+1
Num_Line=0

READ(20, '(A)') ScatteringFile
ScatteringFile = ADJUSTL(ScatteringFile)
READ(20, *) weights(i)
READ(20, *) rmin_n, rmax_n

IF(ScatteringFile(1:1) .NE. '#') THEN !Electron scattering file exists
  used_data_sets(i) = .TRUE.
  FileNameLength=LEN_TRIM(ScatteringFile)
  File_End=.FALSE.

  !Read r_n and gr_n data
  !read r_n first, then gr_n
  !First line is comment
  !At first, count how many data pairs are in the file
  !-1 is chosen as the last point
  OPEN(30,FILE=ScatteringFile(1:FileNameLength),IOSTAT=status1,STATUS='OLD') 
   IF(status1 .EQ. 0) THEN !Open succeeds
     READ(30, '(A80)') comment1
     DO WHILE( .NOT. File_End)
        READ(30, *) Indicator_End
        IF(ABS(Indicator_End+1.0) .LE. 1E-6) THEN !File end is reached
          !PRINT *, 'Neutron scattering file end is reached!'
          EXIT !Exit this loop
        ELSE
          Num_Line = Num_Line + 1
          !READ(30, *) 
        ENDIF
     ENDDO

     REWIND(30) !Go to the beginning of the file

     READ(30, '(A80)') comment1
     ALLOCATE(TempData(2*Num_Line),STAT=Stat_Allocate1)
     IF(Stat_Allocate1 .EQ. 0) THEN
       READ(30, *) TempData
       ALLOCATE(r_n(Num_Line), STAT=Stat_Allocate2)
       ALLOCATE(gr_n(Num_Line), STAT=Stat_Allocate3)

       IF ((Stat_Allocate2 .EQ. 0) .AND. (Stat_Allocate3 .EQ. 0)) THEN
          r_n=TempData(1:2*Num_Line:2)
          gr_n=TempData(2:2*Num_Line:2)
       ELSE
         PRINT *, 'Neutron scattering part 2 or 3 fails!'
         RETURN
       ENDIF !Allocate2 and 3
       
     ELSE
       PRINT *, 'Neutron scattering part allocation fail'
       RETURN
     ENDIF

     DEALLOCATE(TempData)
    
   ELSE
    PRINT *, 'Open Neutron scattering file fails'
   ENDIF
  CLOSE(30)
ELSE 
  used_data_sets(i) = .FALSE.
ENDIF !Neutron scattering file exists

!*******************READ neutron scattering file**********************

!*******************************************************************
!*****************Read X-ray scattering file**********************
!*******************************************************************

i=i+1
Num_Line=0

READ(20, '(A)') ScatteringFile
ScatteringFile = ADJUSTL(ScatteringFile)
READ(20, *) weights(i)
READ(20, *) rmin_x, rmax_x


IF(ScatteringFile(1:1) .NE. '#') THEN !X-ray scattering file exists 

  used_data_sets(i) = .TRUE.
  FileNameLength=LEN_TRIM(ScatteringFile)
  File_End=.FALSE.
  !Read r_x and gr_x data
  !read r_x first, then gr_x
  !First line is comment
  !At first, count how many data pairs are in the file
  !-1 is chosen as the last point
  !OPEN(30,FILE=ScatteringFile(1:FileNameLength),IOSTAT=status1,STATUS='OLD') - This sucks. It never works. Total junk. God damn wasting time. - JWH 052210
  open(30, file="xrd_zr50_kramer_reduced.txt", IOSTAT=status1,STATUS='OLD')
   IF(status1 .EQ. 0) THEN !Open succeeds
     READ(30, '(A80)') comment1
     DO WHILE( .NOT. File_End)
        READ(30, *) Indicator_End
        IF(ABS(Indicator_End+1.0) .LE. 1E-6) THEN !File end is reached
          !PRINT *, 'X-ray scattering file end is reached!'
          EXIT !Exit this loop
        ELSE
          Num_Line = Num_Line + 1
          !READ(30, *) 
        ENDIF
     ENDDO                          
     REWIND(30) !Go to the beginning of the file

     READ(30, '(A80)') comment1
     ALLOCATE(TempData(2*Num_Line),STAT=Stat_Allocate1)
     IF(Stat_Allocate1 .EQ. 0) THEN
       READ(30, *) TempData
       ALLOCATE(r_x(Num_Line), STAT=Stat_Allocate2)
       ALLOCATE(gr_x(Num_Line), STAT=Stat_Allocate3)

       IF ((Stat_Allocate2 .EQ. 0) .AND. (Stat_Allocate3 .EQ. 0)) THEN
          r_x=TempData(1:2*Num_Line:2)

          gr_x=TempData(2:2*Num_Line:2)
       ELSE
         PRINT *, 'X-ray scattering part 2 or 3 fails!'
         RETURN
       ENDIF !Allocate2 and 3
       
     ELSE
       PRINT *, 'X-ray scattering part allocation fail'
       RETURN
     ENDIF

     DEALLOCATE(TempData)
    
   ELSE
    PRINT *, 'Open X-ray scattering file fails'
   ENDIF
  CLOSE(30)
ELSE 
  used_data_sets(i) = .FALSE.
ENDIF !X-ray scattering file exists

!*******************READ X-ray scattering file**********************


!*******************************************************************
!******************Read FEM file************************************
!*******************************************************************

i=i+1
Num_Line=0

READ(20, '(A)') FemFile
FemFile = ADJUSTL(FemFile)
!write(*,*)FemFile
READ(20, *) weights(i)
READ(20, *) scale_fac
READ(20, *) nphi, npsi, ntheta
READ(20, *) Q
!write(*,*)Q
READ(20, *) fem_algorithm
!write(*,*)fem_algorithm
READ(20, *) pixel_distance
!write(*,*)pixel_distance

IF(FemFile(1:1) .NE. '#') THEN !FEM file exists
  used_data_sets(i) = .TRUE.
  FileNameLength=LEN_TRIM(FemFile)
  File_End=.FALSE.

  !Read k ,V, and V_err data
  !read k first, then V, last V_err
  !First line is comment
  !At first, count how many data pairs are in the file
  !-1 is chosen as the last point
  OPEN(30,FILE=FemFile(1:FileNameLength),IOSTAT=status1,STATUS='OLD') 
   IF(status1 .EQ. 0) THEN !Open succeeds
     READ(30, '(A80)') comment1
     DO WHILE( .NOT. File_End)
        READ(30, *) Indicator_End
        IF(ABS(Indicator_End+1.0) .LE. 1E-6) THEN !File end is reached
          !PRINT *, 'FEM file end is reached'
          EXIT !Exit this loop
        ELSE
          Num_Line = Num_Line + 1
          !READ(30, *) 
        ENDIF
     ENDDO

     REWIND(30) !Go to the beginning of the file

     READ(30, '(A80)') comment1
     ALLOCATE(TempData(4*Num_Line),STAT=Stat_Allocate1)
     IF(Stat_Allocate1 .EQ. 0) THEN
       READ(30, *) TempData
       ALLOCATE(k(Num_Line), STAT=Stat_Allocate2)
       ALLOCATE(V(Num_Line), STAT=Stat_Allocate3)
       ALLOCATE(V_err(Num_Line), STAT=Stat_Allocate4)
       ALLOCATE(V_background(Num_Line))

       IF ((Stat_Allocate2 .EQ. 0) .AND. (Stat_Allocate3 .EQ. 0) .AND. (Stat_Allocate4 .EQ. 0)) THEN
           k=TempData(1:4*Num_Line:4)
           V=TempData(2:4*Num_Line:4)
           V_err=TempData(3:4*Num_Line:4)
           V_background=TempData(4:4*Num_Line:4)
       ELSE
         PRINT *, 'FEM part 2 or 3, or 4 fails!'
         RETURN
       ENDIF !Allocate2, 3 and 4
       
     ELSE
       PRINT *, 'FEM part allocation fail'
       RETURN
     ENDIF

     DEALLOCATE(TempData)
    
   ELSE
    PRINT *, 'Open FEM file fails'
   ENDIF
  CLOSE(30)
ELSE 
  used_data_sets(i) = .FALSE.
ENDIF !FEM file exists

!*******************READ FEM file**********************

CLOSE(20)





!CHARACTER, ALLOCATABLE :: filename


END SUBROUTINE read_inputs

END MODULE ReadInputs

