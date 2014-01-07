program fft_test3

  use MKL_DFTI

  implicit none

  real, dimension(:,:,:), pointer :: dat
  type(DFTI_DESCRIPTOR), pointer :: ftd

  integer d1, d2(2), d3(3)
  integer istat


  d2(1) = 5
  d2(2) = 4

  d3(1) = 5
  d3(2) = 4
  d3(3) = 3

  istat = DftiCreateDescriptor(ftd, DFTI_DOUBLE, DFTI_REAL, 3, d3)
  write (*,*) 'After CreateDescriptor, istat =',istat
  write (*,*) 'error message: ',trim(DftiErrorMessage(istat))

  istat = DftiSetValue(ftd, DFTI_PACKED_FORMAT, DFTI_PERM_FORMAT)
  write (*,*) 'After SetValue, istat = ',istat
  write (*,*) 'error message: ',trim(DftiErrorMessage(istat))

  istat = DftiCommitDescriptor(ftd)
  write (*,*) 'After CommitDescriptor, istat=',istat
  write (*,*) 'error mesasge: ',trim(DftiErrorMessage(istat))

  istat = DftiFreeDescriptor(ftd)

end program fft_test3
  
