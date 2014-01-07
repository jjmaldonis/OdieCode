program fft_test

  use MKL_DFTI

  implicit none

  integer, parameter :: PIX = 512

  complex :: c1d(PIX), c2d(PIX,PIX), c3d(PIX,PIX,PIX)
  complex :: c2d_1(PIX*PIX), c3d_1(PIX*PIX*PIX)
  integer :: status, s2d(2), s3d(3)

  equivalence(c2d, c2d_1)
  equivalence(c3d, c3d_1)

  type(DFTI_DESCRIPTOR), pointer :: c1d_handle, c2d_handle, c3d_handle

  write (*,*) 'Program fft_test.'
  write (*,*)

  write (*,*) 'Initializing data.'
  c1d = cmplx(0.0, 0.0)
  c2d = cmplx(0.0, 0.0)
  c3d = cmplx(0.0, 0.0)
  write (*,*)

  write (*,*) 'Computing 1D transform.'
  status = DftiCreateDescriptor(c1d_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, 128)
  status = DftiCommitDescriptor(c1d_handle)
  status = DftiComputeForward(c1d_handle, c1d)
  status = DftiFreeDescriptor(c1d_handle)
  write (*,*) '1D transform complete.'
  write (*,*)

  write (*,*) 'Computing 2D transform.'
  s2d = PIX
  status = DftiCreateDescriptor(c2d_handle, DFTI_SINGLE, DFTI_COMPLEX, 2, s2d)
  status = DftiCommitDescriptor(c2d_handle)
  status = DftiComputeForward(c2d_handle, c2d_1)
  status = DftiFreeDescriptor(c2d_handle)
  write (*,*) '2D transform complete.'
  write (*,*)

  write (*,*) 'Computing 3D transform.'
  s3d = PIX
  status = DftiCreateDescriptor(c3d_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, s3d)
  status = DftiCommitDescriptor(c3d_handle)
  status = DftiComputeForward(c3d_handle, c3d_1)
  status = DftiFreeDescriptor(c3d_handle)
  write (*,*) '3D transform complete.'


end program fft_test
