! femsim, completely written from scratch in Fortran, based
! on modules developed for rmc.
!
! usage is:
!
!     femsim <model file>
!
! parameters can be fed in by redirecting a text file or by hand.
!
! Changelog
!   first version, 03-06-09 pmv
!
!   Define a logical variable use_femsim by Feng Yi on 03/19/2009
!   which is an optional argument in subroutine fem
!
!   Add call subroutine I_average to calculate intensity average as a function of k by FY on 03/30/2009


program femsim

  use model_mod
  use fem_mod

  implicit none

  type(model) :: m
  real, allocatable, dimension(:) :: k, vk, i_k
  real res
  real, dimension(:,:), pointer :: scatfac_e
  character(len=256) :: c, model_filename, outbase, outfile
  integer :: len, istat, nk, i
  real :: kstart, kstep
  integer :: nphi, ntheta, npsi
  LOGICAL use_femsim

  ! parse command line for model name
  if (command_argument_count() /= 1) then
     write (*,*) 'femsim invoked with the wrong number of arguments.  Proper usage is:'
     write (*,*) ' '
     write (*,*) '    femsim <input model file>'
     write (*,*) ' '
     stop
  else
     call get_command_argument(1, c, len, istat)
     if (istat == 0) then
        model_filename = trim(c)
     else
        write (*,*) 'Cannot read model file name.  Exiting.'
        stop
     end if
  endif

  write (*,*) 'Fortran femsim v 1.0; 03-06-09'
  write (*,*) ' '

  call read_model(model_filename, c, m, istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open model file: ', model_filename
     write (*,*) 'Exiting.'
     stop
  endif
  
  ! FEM parameters
  write (*,*) ' '
  write (*,*) 'Please enter the real-space pixel diameter in Angstroms:'
  read (*,*) res
  write (*,*) 'Please enter the first k point in 1/Angstroms:'
  read (*,*) kstart
  write (*,*) 'Please enter the stepsize in k in 1/Angstroms:'
  read (*,*) kstep
  write (*,*) 'Please enter the number of points in k:'
  read (*,*) nk

  ! Simulation parameters
  write (*,*) 'Please enter the number of rotations in theta:'
  read (*,*) ntheta
  write (*,*) 'Please enter the number of rotations in phi:'
  read (*,*) nphi
  write (*,*) 'Please enter the number of rotations in psi:'
  read (*,*) npsi

  write (*,*) 'Please enter the output filename:'
  read (*,*) outbase

  ! Set-up parameters and arrays for fem_initialize
  res = 0.5*res  ! fem_initialize wants pixel radius

  ALLOCATE(k(0))  
  DEALLOCATE(k)
  WRITE(*,*) 'DEALLCATE k succeeds!'
  allocate(k(nk), vk(nk), i_k(nk),stat=istat)
  if( istat /= 0) then
     write (*,*) 'Cannot allcoate memory for V(k) in top level program.'
     write (*,*) 'Exiting.'
     stop
  endif

  k = (/ ( (kstart + kstep*real(i)), i=0, nk-1 ) /)

  ! Open the V(k) output file
  outfile = trim(outbase) // '_vk.out'
  open(unit=300, file=outfile, form='formatted', status='replace', iostat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open output file: ',outfile
     write (*,*) 'Exiting.'
     stop
  endif

  ! Check to make sure we can open the intensities output file before doing
  ! the full calculation
  outfile = trim(outbase) // '_int.out'
  open(unit=301, file=outfile, form='formatted', status='replace', iostat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open output file: ',outfile
     write (*,*) 'Exiting.'
     stop
  endif
  close(301)

  ! fem initialize
  write (*,*) ' '
  write (*,*) 'Initializing FEM calculation.'
  call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfac_e, istat)
  if (istat /= 0) then
     write (*,*) 'Failure in fem_initialize.'
     write (*,*) 'Exiting.'
     stop
  endif

  ! fem calculate with hutch
  use_femsim = .TRUE.
  write (*,*) ' '
  write (*,*) 'Calculating V(k)'
  call fem(m, res, k, vk, scatfac_e, istat, use_femsim)
  if (istat /=0 ) then
     write (*,*) 'Failure in fem subroutine.'
     write (*,*) 'Exiting.'
     stop
  endif

  ! Write the output files
  write (*,*) ' '
  write (*,*) 'Writing output to files.'
  
  write (300, *) 'k   V(k)'
  do i=1, nk
     write (300, *) k(i), vk(i)
  enddo
  close(unit=300)

  call write_intensities(outfile, k, istat)
  if( istat /= 0) then
     write (*,*) 'Intensities file output failed to file: ',outfile
  endif

  
  close(301)
 



  
  

  CALL I_average(i_k)
  
  outfile = trim(outbase) // '_average_i.out'
  open(unit=302, file=outfile, form='formatted', status='replace', iostat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open output file: ',outfile
     write (*,*) 'Exiting.'
     stop
  endif
  WRITE(302,*) 'k I(k)'
  DO i=1, nk
	WRITE(*,*) k(i), i_k(i) !debug
	WRITE(302,*) k(i), i_k(i)
  ENDDO

  CLOSE(302)


end program femsim
