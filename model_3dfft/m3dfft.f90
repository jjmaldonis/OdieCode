program m3dfft

  use mkl_dfti
  use model_mod
  use atompot
  use gfx
  use ft_filters_3d

  implicit none

  ! parameters
  integer, parameter :: TOO_BIG = 1024  ! how many pixels is too many pixels?

  ! input parameters from the user
  character (len=256) :: modelname, outname, potname, ftname
  logical :: pot_out_yn, ft_out_yn
  integer :: npix, pot_out_pix, ft_out_pix

  ! internal variables for the FT
  real, pointer, dimension(:,:,:) :: modelpot
  complex, pointer, dimension(:) :: modelpot_ft
  real, pointer, dimension(:,:,:) :: outblock
  type(model) :: m
  type(DFTI_DESCRIPTOR), pointer :: ft_desc

  ! real space and reciprocal space size and filtering variables
  real :: dxx, dyy, dzz, dkx, dky, dkz, kxmax, kymax, kzmax


  ! counters, input buffers, etc.
  character (len=512) :: comment, line
  integer :: istat, i, j, k, l, lout, hout, mid, s(3)
  complex :: t

  interface 
     function t2o(ii, jj, kk, nx, ny)
       integer :: t2o
       integer, intent(in) :: ii, jj, kk, nx, ny
     end function t2o
  end interface


  write (*,*) 'Model 3D FFT: compute the 3D Fourier transform of a model atomic potential.'
  write (*,*) 

  ! Get user input parameters
  write (*,*) 'Please enter the model file name:'
  read (*,*) modelname

  call read_model(modelname, comment, m, istat)
  if( istat /= 0) then
     write (*,*) 'Unable to open model file ',trim(modelname),'.  Exiting.'
     stop
  endif

  write (*,*) 'How many pixels in along one dimension of the FFT?'
  read (*,*) npix
  if (npix > TOO_BIG) then
     write (*,*) 'Memory required is large.  This may run VERY slowly.'
  endif
  if (mod(npix, 2) /= 0) then
     write (*,*) 'Number of pixels must be divisible by 2.'
     stop
  endif

  write (*,*) 'Output the calculated potential (y/n)?'
  read (*,*) line
  if (line == 'y') then 
     pot_out_yn = .TRUE.
  elseif (line == 'n') then
     pot_out_yn = .FALSE.
  else
     write (*,*) 'Unrecognized input: ',line
     stop
  endif

  if (pot_out_yn) then
     write (*,*) 'How many pixels along one dimension in the potential output?'
     read (*,*) pot_out_pix
     if (pot_out_pix > npix) then
        write (*,*) 'Can only output the total number of pixels, ',npix,'.'
        pot_out_pix = npix
     endif
     if (mod(pot_out_pix, 2) /= 0) then
        write (*,*) 'Number of pixels must be divisible by 2.'
        stop
     endif
  endif

  write (*,*) 'Output the FT modulus (y/n)?'
  read (*,*) line
  if (line == 'y') then
     ft_out_yn = .TRUE.
  elseif (line == 'n') then
     ft_out_yn = .FALSE.
  else
     write (*,*) 'Unrecognized input: ',line
  endif

  if(ft_out_yn) then
     write (*,*) 'How many pixels along one dimension in the FT output?'
     read (*,*) ft_out_pix
     if (ft_out_pix > npix) then
        write (*,*) 'Can only onput the total number of pixels, ',npix,'.'
        ft_out_pix = npix
     endif
     if (mod(ft_out_pix, 2) /= 0) then
        write (*,*) 'Number of pixels must be divisible by 2.'
        stop
     endif
  endif

  if (.NOT. (ft_out_yn .OR. pot_out_yn) ) then
     write (*,*) 'You really should output something,or the calculation isnt worth doing.  Exiting.'
     stop
  endif

  write (*,*) 'Please enter the output filename base: '
  read (*,*) outname
  potname = trim(outname) // '_pot.gfx'
  ftname = trim(outname) // '_ft.gfx'

  ! step sizes
  dxx = m%lx/npix
  dyy = m%ly/npix
  dzz = m%lz/npix
  dkx = 1.0 / m%lx
  dky = 1.0 / m%ly
  dkz = 1.0 / m%lz
  kxmax = 1.0/(2.0*dxx)
  kymax = 1.0/(2.0*dyy)
  kzmax = 1.0/(2.0*dzz)

  ! Echo the input parameters back
  write (*,*)
  write (*,*) 'Calculating 3D FFT with the following parameters:'
  write (*,*) 'Input model: ',trim(modelname)
  write (*,*) 'FFT pixels: ',npix
  write (*,*) 'Real-space sampling in Angstroms is: '
  write (*,*) '   x: ',m%lx/npix
  write (*,*) '   y: ',m%ly/npix
  write (*,*) '   z: ',m%lz/npix
  write (*,*) 'Reciprocal space sampling in 1/Angstroms is: '
  write (*,*) '  kx: start: ',-1.0*kxmax,'  step:',dkx
  write (*,*) '  ky: start: ',-1.0*kymax,'  step:',dky
  write (*,*) '  kz: start: ',-1.0*kzmax,'  step:',dkz
  if(pot_out_yn) then
     write (*,*) 'Potential output to file: ',trim(potname)
     write (*,*) 'Using the center ',pot_out_pix,' pixels.'
     write (*,*) 'Real-space sampling in Angstroms is: '
     write (*,*) '   x: ',m%lx/npix
     write (*,*) '   y: ',m%ly/npix
     write (*,*) '   z: ',m%lz/npix
  endif
  if(ft_out_yn) then
     write (*,*) 'FT output to file: ',trim(ftname)
     write (*,*) 'Using the center ',ft_out_pix,' pixels.'
     write (*,*) 'Reciprocal space sampling in 1/Angstroms is: '
     write (*,*) '  kx: start: ',-0.5*dkx*ft_out_pix,'  step:',dkx
     write (*,*) '  ky: start: ',-0.5*dky*ft_out_pix,'  step:',dky
     write (*,*) '  kz: start: ',-0.5*dkz*ft_out_pix,'  step:',dkz
  endif

  

  ! Do the calculations and output the results
  write (*,*)
  write (*,*)
  write (*,*) 'Calculating model potential.'
  call int_model_pot3D(m, npix, npix, npix, modelpot)
  write (*,*) 'Model potential calculated.'
  
  if(pot_out_yn) then
     write (*,*)
     write (*,*) 'Writing output potential.'
     lout = (npix/2) - (pot_out_pix/2)
     hout = (npix/2) + (pot_out_pix/2)
     call Write3DGFX(potname, modelpot(lout:hout, lout:hout, lout:hout), istat)
  endif
  
  write (*,*) 'Calculating the Fourier Transform.'

  ! copy the real modelpot to the complex, 1D modelpot_ft, deallocate modelpot
  allocate(modelpot_ft(npix*npix*npix), stat=istat)
  if( istat /= 0) then
     write (*,*) 'Cannot allocate memory for modelpot_ft.'
  endif

  l = 1
  do k=1,npix
     do j=1,npix
        do i=1,npix
           modelpot_ft(l) = cmplx(modelpot(i,j,k), 0.0)
           l=l+1
        enddo
     enddo
  enddo
  deallocate(modelpot)

  ! compute the FT
  s = npix
  istat = DftiCreateDescriptor(ft_desc, DFTI_SINGLE, DFTI_COMPLEX, 3, s)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT CreateDescriptor: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiCommitDescriptor(ft_desc)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT CommitDescriptor: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiComputeForward(ft_desc, modelpot_ft)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT ComputeForward: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiFreeDescriptor(ft_desc)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT FreeDescriptor: ',trim(DftiErrorMessage(istat))
  endif
  write (*,*) 'Done.'

  ! swap the octants around so that zero frequency is in the middle
  write (*,*) 'Octant swap.'
  mid = npix/2
  do k=1,npix
     do j=1,npix
        do i=1,mid
           t = modelpot_ft(t2o(i, j, k, npix, npix))
           modelpot_ft(t2o(i,j,k, npix, npix)) = modelpot_ft(t2o(i+mid,j,k,npix,npix))
           modelpot_ft(t2o(i+mid,j,k,npix,npix)) = t
        enddo
     enddo
  enddo
  do k=1,npix
     do j=1,mid
        do i=1,npix
           t = modelpot_ft(t2o(i, j, k, npix, npix))
           modelpot_ft(t2o(i,j,k, npix, npix)) = modelpot_ft(t2o(i,j+mid,k,npix,npix))
           modelpot_ft(t2o(i,j+mid,k,npix,npix)) = t
        enddo
     enddo
  enddo
  do k=1,mid
     do j=1,npix
        do i=1,npix
           t = modelpot_ft(t2o(i, j, k, npix, npix))
           modelpot_ft(t2o(i,j,k, npix, npix)) = modelpot_ft(t2o(i,j,k+mid,npix,npix))
           modelpot_ft(t2o(i,j,k+mid,npix,npix)) = t
        enddo
     enddo
  enddo
  write (*,*) 'Done.'
     
  ! write the FT to a file
  if(ft_out_yn) then

     ! Calculate the modulus of the desired sub-block and write it to a file.
     write (*,*) 'Calculating FT modulus and writing to file.'
     allocate(outblock(ft_out_pix, ft_out_pix, ft_out_pix), stat=istat)
     if (istat /= 0) then
        write (*,*) 'Cannot allocate memory for FT output.'
        stop
     endif
     lout = (npix/2) - (ft_out_pix)/2
     hout = (npix/2) + (ft_out_pix)/2
     l=1
     do k=1,npix
        do j=1,npix
           do i=1, npix
              if( (k < hout) .AND. (k > lout) ) then
                 if( (j < hout) .AND. (j > lout) ) then
                    if( (i < hout) .AND. (i > lout) ) then
                       outblock( (i-lout),(j-lout),(k-lout)) = cabs(modelpot_ft(l))
                    endif
                 endif
              endif
              l = l+1
           enddo
        enddo
     enddo
     
     call Write3DGFX(ftname, outblock, istat)
     deallocate(outblock)
  endif
  
  ! filter and back-transform
  ! for now do a simple bandwidth limit as a test
  write (*,*) 'Bandwidth limiting to 1.0 1/Angstroms.'
  do k=1,npix
     do j=1,npix
        do i=1,npix
           if ( sqrt( (real(i)*dkx - kxmax)**2 + (real(j)*dky-kymax)**2 + (real(k)*dkz-kzmax)**2 ) > 1.0) then
              modelpot_ft(t2o(i, j, k, npix, npix)) = 0.0
           endif
        enddo
     enddo
  enddo

  ! compute the IFT
  write (*,*) 'Computing the inverse FT.'
  s = npix
  istat = DftiCreateDescriptor(ft_desc, DFTI_SINGLE, DFTI_COMPLEX, 3, s)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT CreateDescriptor: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiSetValue(ft_desc, DFTI_BACKWARD_SCALE, 1.0/(npix**3))
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT DftiSetValue: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiCommitDescriptor(ft_desc)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT CommitDescriptor: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiComputeBackward(ft_desc, modelpot_ft)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT ComputeForward: ',trim(DftiErrorMessage(istat))
  endif
  istat = DftiFreeDescriptor(ft_desc)
  if (istat /= 0) then
     write (*,*) 'Error from MKL FFT FreeDescriptor: ',trim(DftiErrorMessage(istat))
  endif
  write (*,*) 'Done.'

  ! Calculate the modulus of the desired sub-block and write it to a file.
  write (*,*) 'Writing back transform to a file.'
  allocate(outblock(pot_out_pix, pot_out_pix, pot_out_pix), stat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot allocate memory for IFT output.'
     stop
  endif
  lout = (npix/2) - (pot_out_pix)/2
  hout = (npix/2) + (pot_out_pix)/2
  l=1
  do k=1,npix
     do j=1,npix
        do i=1, npix
           if( (k < hout) .AND. (k > lout) ) then
              if( (j < hout) .AND. (j > lout) ) then
                 if( (i < hout) .AND. (i > lout) ) then
                    outblock( (i-lout),(j-lout),(k-lout)) = cabs(modelpot_ft(l))
                 endif
              endif
           endif
           l = l+1
        enddo
     enddo
  enddo
  
  call Write3DGFX(potname, outblock, istat)
  deallocate(outblock)
  
  write (*,*) 'Program model 3D FFT finished.'

end program m3dfft


function t2o(ii, jj, kk, nx, ny)
  integer :: t2o
  integer, intent(in) :: ii, jj, kk, nx, ny
  
  t2o = ii + (jj-1)*nx + (kk-1)*(nx*ny)
  
end function t2o

