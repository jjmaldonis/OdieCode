! functions for applying 3D Fourier-space filters to complex waves
! intended for use with m3dfft.f90
!
! begun 05-11-09 pmv

module ft_filters_3D

  !interface 
  !   function t2o(ii, jj, kk, nx, ny)
  !     integer :: t2o
  !     integer, intent(in) :: ii, jj, kk, nx, ny
  !   end function t2o
  !end interface


contains

  subroutine LowPass(modelpot_ft, klimit, krolloff, npix_x, npix_y, npix_z, dkx, dky, dkz, kxmax, kymax, kzmax)
    complex, pointer, dimension(:), intent(inout) :: modelpot_ft
    real, intent(in) :: klimit, krolloff
    integer, intent(in) :: npix_x, npix_y, npix_z
    real, intent(in) :: dkx, dky, dkz, kxmax, kymax, kzmax

    integer :: i, j, k
    real kx, ky, kz, klimit2, kr, kr2
    klimit2 = klimit*klimit

    do k=1,npix_z
       kz = real(k)*dkz - kzmax
       do j=1,npix_y
          ky = real(j)*dky - kymax
          do i=1,npix_x
             kx = real(i)*dkx - kxmax
             kr2 = kx*kx + ky*ky + kz*kz
             if (kr2 > klimit2) then 
                kr = sqrt(kr2)
                modelpot_ft(t2o(i, j, k, npix_x, npix_y)) = &
                     modelpot_ft(t2o(i, j, k, npix_x, npix_y))*exp( -1.0 * (kr-klimit)**2 / krolloff )
             endif
          enddo
       enddo
    enddo

  end subroutine LowPass

  
  subroutine HighPass(modelpot_ft, klimit, krolloff, npix_x, npix_y, npix_z, dkx, dky, dkz, kxmax, kymax, kzmax)
    complex, pointer, dimension(:), intent(inout) :: modelpot_ft
    real, intent(in) :: klimit, krolloff
    integer, intent(in) :: npix_x, npix_y, npix_z
    real, intent(in) :: dkx, dky, dkz, kxmax, kymax, kzmax

    integer :: i, j, k
    real kx, ky, kz, klimit2, kr, kr2
    klimit2 = klimit*klimit

    do k=1,npix_z
       kz = real(k)*dkz - kzmax
       do j=1,npix_y
          ky = real(j)*dky - kymax
          do i=1,npix_x
             kx = real(i)*dkx - kxmax
             kr2 = kx*kx + ky*ky + kz*kz
             if (kr2 < klimit2) then 
                kr = sqrt(kr2)
                modelpot_ft(t2o(i, j, k, npix_x, npix_y)) = &
                     modelpot_ft(t2o(i, j, k, npix_x, npix_y))*exp( -1.0 * (kr-klimit)**2 / krolloff )
             endif
          enddo
       enddo
    enddo

  end subroutine HighPass


  subroutine BandPass(modelpot_ft, klow, khigh, krolloff, npix_x, npix_y, npix_z, dkx, dky, dkz, kxmax, kymax, kzmax)
    complex, pointer, dimension(:), intent(inout) :: modelpot_ft
    real, intent(in) :: klow, khigh, krolloff
    integer, intent(in) :: npix_x, npix_y, npix_z
    real, intent(in) :: dkx, dky, dkz, kxmax, kymax, kzmax

    integer :: i, j, k
    real kx, ky, kz, klow2, khigh2, kr, kr2
    klow2 = klow*klow
    khigh2 = khigh*khigh

    do k=1,npix_z
       kz = real(k)*dkz - kzmax
       do j=1,npix_y
          ky = real(j)*dky - kymax
          do i=1,npix_x
             kx = real(i)*dkx - kxmax
             kr2 = kx*kx + ky*ky + kz*kz
             if (kr2 < klow2) then 
                kr = sqrt(kr2)
                modelpot_ft(t2o(i, j, k, npix_x, npix_y)) = &
                     modelpot_ft(t2o(i, j, k, npix_x, npix_y))*exp( -1.0 * (kr-klow)**2 / krolloff )
             elseif (kr2 > khigh2) then
                kr = sqrt(kr2)
                modelpot_ft(t2o(i, j, k, npix_x, npix_y)) = &
                     modelpot_ft(t2o(i, j, k, npix_x, npix_y))*exp( -1.0 * (kr-khigh)**2 / krolloff )
             endif
          enddo
       enddo
    enddo 

  end subroutine BandPass


  subroutine SpheresPass(modelpot_ft, kcent, krad, krolloff, npix_x, npix_y, npix_z, dkx, dky, dkz, kxmax, kymax, kzmax)
    complex, pointer, dimension(:), intent(inout) :: modelpot_ft
    real, dimension(:,:), intent(in) :: kcent
    real, intent(in) :: krad, krolloff
    integer, intent(in) :: npix_x, npix_y, npix_z
    real, intent(in) :: dkx, dky, dkz, kxmax, kymax, kzmax

    integer :: i, j, k, nsphere
    real kx, ky, kz, klow2, khigh2, kr, kr2

    nsphere = size(kcent, 1)

    klow2 = klow*klow
    khigh2 = khigh*khigh

    write(*,*) "WARNING! THE SPHERESPASS SUBROUTINE PROBABLY DOESNT WORK!  -Jason"
    write(*,*) "klow and khigh are undefined but used."

    do k=1,npix_z
       kz = real(k)*dkz - kzmax
       do j=1,npix_y
          ky = real(j)*dky - kymax
          do i=1,npix_x
             kx = real(i)*dkx - kxmax
             kr2 = kx*kx + ky*ky + kz*kz
             if (kr2 < klow2) then 
                kr = sqrt(kr2)
                modelpot_ft(t2o(i, j, k, npix_x, npix_y)) = &
                     modelpot_ft(t2o(i, j, k, npix_x, npix_y))*exp( -1.0 * (kr-klow)**2 / krolloff )
             elseif (kr2 > khigh2) then
                kr = sqrt(kr2)
                modelpot_ft(t2o(i, j, k, npix_x, npix_y)) = &
                     modelpot_ft(t2o(i, j, k, npix_x, npix_y))*exp( -1.0 * (kr-khigh)**2 / krolloff )
             endif
          enddo
       enddo
    enddo 

  end subroutine SpheresPass

end module ft_filters_3D

