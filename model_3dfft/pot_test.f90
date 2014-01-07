! Test vatom in atompot module by generating data for Fig 5.4 from
! Kirkland's book.  As of 02-05-09, everything looks fine.  pmv
program pot_test

  use atompot
  implicit none

  real, dimension(200) :: r, C, Si, Au, Cu, U
  integer i
  real ipot

  r = (/ ( (0.5/200.0)*real(i), i=1, 200 ) /)

  do i=1, 200
     C(i) = vatom(6, r(i))
     Si(i) = vatom(14, r(i))
     Cu(i) = vatom(29, r(i))
     Au(i) = vatom(79, r(i))
     U(i) = vatom(92, r(i))
  enddo

  write (*,*) 'r   C   Si   Cu   Au    U'
!  write (*,*) ( (r(i), C(i), Si(i), Cu(i), Au(i)), i=1,200 )

  ipot = int_pot(14, 1.0, 1.0, 1.0, 0.5, 1.5, 0.5, 1.5, 0.5, 1.5)
  write (*,*) 'ipot = ',ipot

end program pot_test
