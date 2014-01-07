program fft_test2

  use MKL_DFTI

  implicit none

  integer, parameter :: PIX = 1024

  complex, pointer, dimension(:) :: c
  integer :: status, s2d(2), s3d(3)

  type(DFTI_DESCRIPTOR), pointer :: c1d_handle, c2d_handle, c3d_handle

  write (*,*) 'Program fft_test.'
  write (*,*)

  write (*,*) 'Computing 1D transform.'
  allocate(c(PIX), stat=status)
  c = cmplx(0.0, 0.0)
  status = DftiCreateDescriptor(c1d_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, 128)
  status = DftiCommitDescriptor(c1d_handle)
  status = DftiComputeForward(c1d_handle, c)
  status = DftiFreeDescriptor(c1d_handle)
  deallocate(c)
  write (*,*) '1D transform complete.'
  write (*,*)

  write (*,*) 'Computing 2D transform.'
  allocate(c(PIX*PIX), stat=status)
  c = cmplx(0.0, 0.0)
  s2d = PIX
  status = DftiCreateDescriptor(c2d_handle, DFTI_SINGLE, DFTI_COMPLEX, 2, s2d)
  status = DftiCommitDescriptor(c2d_handle)
  status = DftiComputeForward(c2d_handle, c)
  status = DftiFreeDescriptor(c2d_handle)
  deallocate(c)
  write (*,*) '2D transform complete.'
  write (*,*)

  write (*,*) 'Computing 3D transform.'
  allocate(c(PIX*PIX*PIX), stat=status)
  c = cmplx(0.0, 0.0)
  s3d = PIX
  status = DftiCreateDescriptor(c3d_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, s3d)
  status = DftiCommitDescriptor(c3d_handle)
  status = DftiComputeForward(c3d_handle, c)
  status = DftiFreeDescriptor(c3d_handle)
  write (*,*) '3D transform complete.'

end program fft_test2
