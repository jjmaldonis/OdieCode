!*********************************************************
!**********This module is written by FY on 12/19/2008****
!**********Include three functions***********************
!*********************************************************

!    function chi_square, check_cutoffs and generate_move are included
!

! Changelog
!    unused arguments (r_e etc) removed from chi_sqaure, 12/20/08 pmv
!    subroutine random_move updated. 01/12/08 Jinwoo Hwang
!    function check_cufoff is completed by Feng Yi on 01/16/2009, may need a test program
!    check_cutoffs is modified by Feng Yi on 01/21/2009
!    removed "use ReadInputs" which doesn't seem to be needed; pmv 1/23/09
!    added interface block for function ran2 to enable clean compilation pmv 03/12/09
!    In check_cutoffs, pbc is added to the pair calculation after the hutch is calling. jwh 04/14/2009
!    Added deallocate(atoms) in check_cutoffs, pmv 04/17/09
!    change if condition when deallocate(atoms) in check_cutoffs, FY 04/17/2009
!    gr_e_err is added - jwh 04/25/2009

MODULE rmc_functions

  use model_mod
  implicit none

  interface
     function ran2(idum)
       real :: ran2
       integer :: idum
     end function ran2
  end interface

CONTAINS

!***************************************************************
!This function is used to calculate chi square
  FUNCTION chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, V, V_err,&
       gr_e_sim, gr_n_sim, gr_x_sim, V_sim, scale_fac,&
	rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)

    !used_data_sets=A logical array determine whether electron, neutron, X-ray or FEM data are available
    !gr_e electron G(r) data
    !gr_n neutron G(r) data
    !gr_x X-ray G(r) data
    !V=FEM intensity variance
    !V_err=FEM measurement variance error
    !gr_e_sim=simulated g(r) for electron scattering
    !gr_n_sim=simulated g(r) for neutron scattering
    !gr_x_sim=simulated g(r) for X-ray scattering
    !V_sim=simulated V(k) for variance data
    !All these all read from main program

    LOGICAL, DIMENSION(4) :: used_data_sets
    REAL,  DIMENSION(4) :: weights
    REAL, POINTER, DIMENSION(:) :: gr_e, gr_e_err
    REAL, POINTER, DIMENSION(:) :: gr_n
    REAL, POINTER, DIMENSION(:) :: gr_x
    REAL, POINTER, DIMENSION(:) :: V, V_err
    REAL, POINTER, DIMENSION(:) :: gr_e_sim,gr_n_sim,gr_x_sim
    REAL, POINTER, DIMENSION(:) :: V_sim
    real, intent(in) :: scale_fac, rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x 
    real, intent (in) :: del_r_e, del_r_n, del_r_x	
    integer, intent(in) ::nk
    real, intent(out) :: chi2_gr, chi2_vk
    
    !********The following are local variables******
    INTEGER i, j
    integer nf   !normalization factor   - JWH 04/25/2009
    REAL chi_square
    REAL,DIMENSION(4) :: sum1 !record summation of each contribution from diffraction data
    !and FEM data
    
    chi_square=1.0
    DO i=1, 4
       sum1(i)=0.0
    ENDDO
    
    !Electron diffraction data
    i=1
    IF(used_data_sets(i)) THEN
       !sum1(i)=sum((gr_e_sim-gr_e)**2)*weights(i)
	nf = 0
	do j=int(rmin_e/del_r_e)+1, int(rmax_e/del_r_e)+1    !jwh -032409
		nf = nf + 1
		sum1(i) = sum1(i)+weights(i)*((gr_e_sim(j)-gr_e(j))/gr_e_err(j))**2
	enddo
	sum1(i) = sum1(i)/nf

    ENDIF
    
    !Neutron diffraction data
    i=i+1
    IF(used_data_sets(i)) THEN
       !sum1(i)=sum((gr_n_sim-gr_n)**2)*weights(i)
	nf = 0
	do j=int(rmin_n/del_r_n)+1, int(rmax_n/del_r_n)+1    !jwh -032409
		nf = nf +1
		sum1(i) = sum1(i)+weights(i)*(gr_n_sim(j)-gr_n(j))**2
	enddo
	sum1(i) = sum1(i)/nf
    ENDIF
    
    !X-ray diffraction data
    i=i+1
    IF(used_data_sets(i)) THEN
       !sum1(i)=sum((gr_x_sim-gr_x)**2)*weights(i)
	nf = 0
	do j=int(rmin_x/del_r_x)+1, int(rmax_x/del_r_x)+1    !jwh -032409
		nf = nf + 1
		sum1(i) = sum1(i)+weights(i)*(gr_x_sim(j)-gr_x(j))**2
	enddo
	sum1(i) = sum1(i)/nf
    ENDIF
    
    !FEM data
    i=i+1
    IF(used_data_sets(i)) THEN
       !sum1(i)=sum(((V_sim-V*scale_fac)/(scale_fac*V_err))**2)*weights(i)
	nf = 0
	do j=1,nk
		nf = nf + 1
		sum1(i) = sum1(i)+weights(4)*((v(j)*scale_fac-v_sim(j))/(scale_fac*V_err(j)))**2
	enddo
	sum1(i) = sum1(i)/nf
    ENDIF

    chi2_gr = sum1(3)
    chi2_vk = sum1(4)

    chi_square=sum(sum1)
    
  END FUNCTION chi_square


!*******************************************************
!*********A new function************************
!*******************************************************

!SUBROUTINE Calc_Type_Order(m, Type_Order) must be called first, then 
!this function can be called
!The check_cufoffs function return .TRUE. if 
!****no any two atoms are within cutoff distance
!Otherwise return .FALSE.

  FUNCTION check_cutoffs(m,cutoff_r,moved_atom)

    LOGICAL check_cutoffs
    REAL,DIMENSION(:,:) ::cutoff_r
    INTEGER moved_atom
    TYPE(model) m

    INTEGER, DIMENSION(:), POINTER :: atoms
    INTEGER  nlist
    INTEGER istat
    REAL radius, temp_x, temp_y, temp_z
    INTEGER i,j
    INTEGER num1, num2
    REAL dist_pair

    !Find the maximum cut-off distance
    radius=MAXVAL(MAXVAL(cutoff_r, 1),1)

    !PRINT *, 'radius is: ', radius

    ! Makes list (in atoms) of the indices of atoms that are in a cube of side 
    ! length radius, centered on the hutch containing  the point (px, py, pz),
    ! using the hutch array.  Useful for calculating G(r).  Returns 1 in istat
    ! if memory allocation fails and -1 if no atoms are found.

    CALL hutch_list_3D(m, m%xx(moved_atom),m%yy(moved_atom),m%zz(moved_atom), radius, atoms, istat, nlist)

    IF (istat .EQ. 1) THEN
       PRINT *, 'memory allocation fails!'
       RETURN
    ENDIF
    
    IF (istat .EQ. -1) THEN
       PRINT *, 'No atom is found'
       check_cutoffs = .TRUE. 
       RETURN
    ENDIF

  !Begin to calculate pair distance
  !Notice, there are totally (nlist-1) atoms

  !First, determine the type of moved_atom 
    DO i=1, m%nelements
       IF(m%znum(moved_atom) .EQ. m%atom_type(i)) THEN
          num1 = i
          EXIT
       ENDIF
    ENDDO
    
    DO i=1, (nlist-1)
       !Do not count the pair distance to itself
       IF (atoms(i) .NE. moved_atom) THEN
          
          DO j=1, m%nelements
             IF(m%atom_type(j) .EQ. m%znum(atoms(i))) THEN
                num2 = j
                EXIT
             ENDIF
	  ENDDO !j=1, m%nelements
          
       !ENDIF  !032409 - jwh
       
       !Calculate the atomic distance
       !Compare with cutoff_r

	temp_x = abs(m%xx(moved_atom) - m%xx(atoms(i)))  !pbc added - JWH 04/14/2009
	temp_y = abs(m%yy(moved_atom) - m%yy(atoms(i)))
	temp_z = abs(m%zz(moved_atom) - m%zz(atoms(i)))

	!write(*,*)"abs=", temp_x, temp_y, temp_z
	
	temp_x = temp_x - m%lx*anint(temp_x/m%lx)
	temp_y = temp_y - m%ly*anint(temp_y/m%ly)
	temp_z = temp_z - m%lz*anint(temp_z/m%lz)
       
       dist_pair = temp_x**2 + temp_y**2 + temp_z**2
             
       dist_pair = SQRT(dist_pair)
       
       IF (dist_pair  .LT. cutoff_r(num1, num2)) THEN
	!write(*,*)dist_pair  , cutoff_r(num1, num2) !debug - jwh 032409
          check_cutoffs=.FALSE.
	   EXIT
	else
		check_cutoffs = .TRUE.  !jwh 032409
       ENDIF

	endif !032409 - jwh
       
    ENDDO !i=1, (nlist-1)

    !if (associated(atoms)) deallocate(atoms) ! pmv 4/17/09
    if (nlist .GT. 1) deallocate(atoms) ! FY 4/17/09

    
  END FUNCTION check_cutoffs
  


  !subroutine random_move---------------------------------------------------------------
  !Generates a random move of one atom.
  !Written for testing purpose here.
  !The mc move will be done some place else later, but it should still update the variables xx_cur, and xx_new (and etc).
  subroutine random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, alpha)
    use mpi  
    type(model), intent(inout) :: m
    integer iseed, w
    real alpha, aa, bb, cc
    real, intent(out) :: xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new    !Cur and new positions of atoms
    real :: rand1, rand2, rand3, rand4


    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    rand3 = ran2(iseed)
    rand4 = ran2(iseed)

    w = int(m%natoms*rand1)+1

    !write(*,*)myid, rand1, rand2, rand3, rand4

    xx_cur = m%xx(w)            !Cur positions of the atom before random move
    yy_cur = m%yy(w)
    zz_cur = m%zz(w)
    
    aa = alpha*(rand2 - 0.5)
    bb = alpha*(rand3 - 0.5)
    cc = alpha*(rand4 - 0.5)
	
    m%xx(w) = m%xx(w) + aa 
    m%yy(w) = m%yy(w) + bb 
    m%zz(w) = m%zz(w) + cc
	
    if(m%xx(w)>m%lx*0.5) m%xx(w)=m%xx(w)-m%lx       !pbc 
    if(m%yy(w)>m%ly*0.5) m%yy(w)=m%yy(w)-m%ly
    if(m%zz(w)>m%lz*0.5) m%zz(w)=m%zz(w)-m%lz
    if(m%xx(w)<-m%lx*0.5) m%xx(w)=m%xx(w)+m%lx
    if(m%yy(w)<-m%ly*0.5) m%yy(w)=m%yy(w)+m%ly
    if(m%zz(w)<-m%lz*0.5) m%zz(w)=m%zz(w)+m%lz
    
    xx_new=m%xx(w)              !new positions of the atom after random move
    yy_new=m%yy(w)
    zz_new=m%zz(w)
  end subroutine random_move


END MODULE rmc_functions
