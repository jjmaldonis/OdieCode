! Data and subroutines to calculate FEM V(k,R) from a model.
!
! WARNING: FEM RESOLUTION R RADIUS VS DIAMETER ISSUE IS NOT RESOLVED!
! For objective aperture radius Q, 0.61/Q is the diameter of the FEM
! pixel.  What we usually call FEM "resolution" is a diameter, not a 
! radius.
!
! WARNING: ACCEPT_MOVE OR REJECT_MOVE MUST BE CALLED AFTER FEM_UPDATE.  OTHERWISE,
! THE NEXT CALL TO FEM_UPDATE WILL CRASH.
! 
! need to write dump_intensities or dump_image subroutine to output to a file either 
! the entire intensities array or just one image at a particular k and rotation
!
! Module data, all private, all persistent
!   nk: total number of k points
!   npix: total number of pixels
!   nrot: total number of rotations
!   pix: (npix x 2) list of pixel (x,y) positions
!   rot: (nrot x 3) list of rotations (theta, phi, psi) Euler angles
!   int: calculated intensities (nk, npix, nrot)
!   int_sq: point by point square of the intensities array
!   int_sum: (nk) list of the sum of all the intensities at a particular k
!   int_sq_sum: (nk) list of the sum of all the intensities**2 at a particular k
!   mrot: (nrot) array of rotated model derived types
!   old_index: (nrot) a list of atoms which a particular rmc atom move has changed in each rotated model
!   old_pos: (nrot) a list of the old positions of the atoms which a particular rmc atom move has
!            changed in each rotated model
!
! Public subroutines:
!
! subroutine fem_initialize(model, resolution, k, nk, npix_x, npix_y, ntheta, nphi, npsi, besselfunc1, besselfunc2, scat_fact, istat)
!   sets up the arrays for a fem calculation, and calculates the lists of pixel positions
!   and rotation angles.  must be called before fem or fem_update. all the numbers of
!   things (nk, ntheta, etc.) must be at least one.
!
! subroutine fem(model, resolution, k, vk, besselfunc1, besselfunc2, scat_fact, istat, use_femsim)
!   calculates the fem v(k) for a model from scratch.  sets up the rotated model array,
!   rotates all the models, calculates all the intensities, and calculates V.  k must
!   contain the list of k points, and vk will contain V(k).  returns non-zero in istat
!   if some memory cannot be allocated.  must be called before the first call of fem_update.
!   fem_initialize must be called before fem.
!   use_femsim is an optional argument determining whether algorithm for femsim calculation
!   or rmc is used, this is added by Feng Yi on 03/19/2009
!
! subroutine fem_update(model, moved_atom_index, resolution, k, vk, besselfunc1, besselfunc2, scat_fact, istat)
!   re-calcualtes fem for the model in which the atom at moved_atom_index has moved.  that atom
!   MUST already be moved in the model!  updates the position of the moved atom in each stored,
!   rotated model, then only recalculates the intensities from pixels which either contained the
!   original atom position or that now contain the new atom position.  fem() must be called
!   before this is called the first time.  fem_accept_move or fem_reject_move must be called before
!   fem_update is called again.
!
! subroutine fem_accept_move()
!   accepts the rmc move tested in fem_update.  what it really does it reset the old_pos and
!   old_index arrays for the next call of fem_update.
!
! subroutine fem_reject_move(model)
!   rejects the rmc move tested by fem_update.  resets all the positions of the moved atoms in all
!   the rotated models.  model MUST have the moved atom put back to its original position.
! 
!! subroutine print_image(m, res, k, besselfunc1, besselfunc2, scatfact_e, istat)
!   Generates all images and prints them at particular k and rotation.
!   It can be called from main program while rmc running.

! subroutine print_image2(m, res, k, besselfunc1, besselfunc2, scatfact_e, istat)
!   Generates all images and prints them
!   It stops the program after printing, so it should be called independantly.
!
!
! Private subroutines:
!
! subroutine init_pix(m, npix_x, npix_y, istat)
!   calculates the list of pixel positions in x and y on the model, based on the model size, the
!   resolution, and the number of pixels.  returns non-zero in istat if memory allocation fails.
!
! subroutine init_rot(ntheta, nphi, npsi, istat)
!   calculates the list of rotation Euler angles based on ntheta, nphi, and npsi.  Returns non-
!   zero in istat if memory allocation fails.
!
! subroutine intensity(m, res, px, py, k, int, j0, a1, scatfact_e, istat)
!   calculates the fem intensity of the model m, at pixel position (px, py) and resolution res
!   for all the k in the list k, and puts them in the 1-D array int.  this is where multiple
!   algorithms for the intensity calculation could be implemented.
!
! subroutine add_pos(p, xx, yy, zz, istat)
!   adds the position (xx, yy, zz) to the end of the position list pos.
!
! subroutine destroy_pos_list(p, np)
!   deallocates all the pointers, etc. in the position_list array p, which in np position_lists
!   long.
!
! subroutine fem_reset_old()
!   resets the old_index and old_pos arrays so fem_update can resuse them
!
!  
!
!
! Changelog
!   first verion, 01-06-09, pmv
!   init_rot is modified by Feng Yi on 01/15/2009
!
!	Modified by Jinwoo Hwang 01/27/2009	
!	subroutines fem_initialize, init_pix were modified. 
!	- fem_initilize assumes 'resolution=pix diameter'.
!	- Now pixel spacing (dr) and position are determined by the resolution and the model size.
!	subroutine fem_initialize reads bessel function and calls electron scattering factors. 
!	- These functions are passed to subroutines fem, fem_update, and intensity.
!	subroutine intensity was modified to calculate I(k) and V(K). 
!	- The code is based on Raymond-Biswas rmc_Vk.
!	subroutine print_image2 is written. It calculates I(k) and prints them out. 
!	-It stops the program after printing, so it should be called separately from the RMC process.  
!
!   Modified by Jinwoo Hwang 01/30/2009
!   gr_i allocation corrected. Now the gr_i is only allocated up to the maximum possible bin size.
!   print_image1 was written for printing I(k) during RMC.
!   For consistency for now, program was set to be "resolution = radius of pixel".
!	Miscellaneous changes made for parameter max_r, npix_x, npix_y
!	Bessel function J0 and aperture function A1 is now calculated in fem_initialize. 
!	
!	Modified by Jinwoo Hwang 02/18/2009
!	Array bound error fixed in const4.
!	All pixels were in positive x-y plane. Now they are corrected to cover the area centered at (0,0).
!
!	Modified by Jinwoo Hwang 02/24/2009
!	Debug done - Now the V(k) caculation yields same result as the old RMC_JWH.
!
!	Debug by Jinwoo Hwang 03/04/09
!	For testing purpose the pixels were fixed with 8x8 specific positions, and rotations were fixed to 3.
!	In fem_update, the double loop (do m=1,npix, do n=1,rot_atom%natoms)didn't have PBC. It's now added.
!	In fem_update, in the same double loop, some typos were fixed (m was written as j). 
!	
!	fem accept and reject process tested and debug was done. The fem_reject_move now includes the return of the
!	int_i, int_sq, m%positions to their old values. - Jinwoo Hwang 03/05/09
!
!       added subroutines write_intensities to output the entire int_i array to a single file.
!       it's also added to the list of public subroutines - pmv 03-06-09
!
!  in intensity, gr_i no longer allocated with 'save' attribute, 3/18/09 pmv
!  in intesntiy, znum_r no longer recalculated from scratch 3/18/09 pmv
!  scatfact_e no longer initialized to zero, then calculated in fem_initialize 3/18/09 pmv
!  unncessary moved_atom initialization in fem_updtae removed 3/18/09 pmv
!  old useless temporary code remove from init_pix 3/18/09 pmv
!  added explicit deallocate(scratch) to add_pos 3/18/09 pmv 
!  added deallocates to fem_reset_old 3/18/09 pmv 
!
!  add a variable femsim in fem subroutine, set it for .TRUE. for only femsim calculation
!  set it for .false. for rmc calculation. Did by Feng Yi on 03/19/2009
!  also in fem subroutine, for femsim calculation, only size 1 of mrot is needed	
!
!  add an optional argument use_femsim to fem subroutine. If use_femsim is present, it will 
!  be assigned to femsim. The default value of femsim is .FALSE.  By Feng Yi on 03/19/2009
!
!  define variable original of type pointer to be allocatable in intensity, and deallocate them
!  depending on whether they are allocated or size of array is 0 in intensity subroutine by Feng Yi on 03/20/2009
!
!    in fem_reset_old, the /=was changed to > -JWH 03/26/2009
!
!  Add two optional arguments rot_begin, rot_end to subroutine fem by FY on 03/26/2009
!  debug - in fem_initialize, a1(0)=1.0 was added. - jwh 03/26/2009

!  add subroutine I_average to calculate average intensity as a function of k
!  change variable type of vk in fem to dummy type by Feng Yi on 03/30/2009
!
!  add output rotation angle in fem subroutine by Feng Yi on 03/31/2009
!  
!  In fem_update, the condition for the "cycle" was corrected. - JWH 04/14/2009
!
!  Deallocate ind in rot_atom%rot_i om fem_update subroutine by FY on 04/17/2009
!  Change  if (associated(p%pos)) to if(p%nat .GT. 0) in add_pos subroutine by FY on 04/17/2009
!
!    in fem_initialize, the size(k) in the const4 is changed to k(size(k)). - JWH 04/23/2009
!
!  in fem_update, the conditions for deallocating rot_atom%rot_i(n)%ind is changed - JWH 04/30/2009
!
!  in fem_initialize, the Al was calculated as Al(r). Now it is fixed to calculate Al(2*pi*Q*r). by JWH 05/01/2009
!  The change on 05/01/2009 is wrong. changed back by JWH on 05/02/2009
!
!  Change const4 in fem_initialize to make sure enough space is allocated
!  Change 2. *t1*t2 to t1*t2 when i .NE. j in calculating gr_i to avoid double counting by FY on 05/04/2009
!
!  Change back to 2*t1*t2 in intensity subroutine by FY on 05/04/2009
!  in fem_initialize, when calculate j0 and j1, i begins from 1 instead of 0 by FY on 05/04/2009
!
!  in the subroutine "intensity", the data tyoe of the scatfact_e was modified
!  in the subroutine "bessel_func", the input data was re-formatted.   - JWH /05/08/2009
!*********************************************************************

module fem_mod

  use  model_mod
  use  RMC_Global
  use  scattering_factors 

  implicit none
  private
  public :: fem_initialize, fem, fem_update, fem_accept_move, fem_reject_move, I_average
  public :: write_intensities, print_image1, print_image2

  type pos_list
     integer :: nat
     real, pointer, dimension(:,:) :: pos
  end type pos_list

  integer, save :: nk, npix, nrot  ! number of k points, pixels, and rotations

  real, save, dimension(:,:), pointer :: pix ! npix x 2 list of pixel positions
  real, save, dimension(:,:), pointer :: rot ! nrot x 3 list of (theta, phi, psi) rotation angles

  real, save, dimension(:,:,:), pointer :: int_i, int_sq  ! nk x npix x nrot
  real, save, dimension(:,:,:), pointer :: old_int, old_int_sq
  real, save, dimension(:), pointer :: int_sum, int_sq_sum  ! nk long sums of int and int_sq arrays for calculating V(k)
  real, save, allocatable, dimension(:) :: j0, A1                                               

  type(model), save, dimension(:), pointer :: mrot  ! array of rotated models
  type(index_list), save, dimension(:), pointer :: old_index
  type(pos_list), save, dimension(:), pointer :: old_pos  


contains

!
! public subroutines
!
  subroutine fem_initialize(m, res, k, nki, ntheta, nphi, npsi, scatfact_e, istat)
    type(model), intent(in) :: m
    real, intent(in) :: res
    real, dimension(:), intent(in) :: k
    integer, intent(in) :: nki, ntheta, nphi, npsi
    real, dimension(:,:), pointer :: scatfact_e
    integer, intent(out) :: istat

    real dr      !Distance between pixels
    real r_max, const1, const2, const3
    integer npix_x, npix_y , bin_max
    integer i, j
    integer const4 
    double precision b_x, b_j0 , b_j1

    r_max = 2*res     !assuming resolution=radius

    bin_max = int(r_max/fem_bin_width)+1

    const1 = twopi*(0.61/res)/fem_bin_width  !(0.61/res = Q) 
    const2 = 1/fem_bin_width
    const3 = (const1/(0.61/res))/const2
    !const4 = int(bin_max*const3*k(size(k)))+1
    const4 = int(bin_max*const3*CEILING(k(SIZE(k))))+1

    allocate(j0(0:const4),a1(0:const4), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Failed to allocate memory for Bessel and Airy functions.'
    endif

    !calculating bessel function
    j0=0.0
    do i=1, const4
       b_x = i*fem_bin_width
       call bessel_func(b_x, b_j0,b_j1)
       j0(i) = sngl(b_j0)
       a1(i) = 2*sngl(b_j1)/b_x
    enddo
    a1(0)=1.0
    j0(0)=1.0

    !OPEN(101,FILE='j0_file.xyz', STATUS='UNKNOWN')
    !DO i=0, const4
    !WRITE(101, *) i*fem_bin_width, j0(i)
    !ENDDO
    !CLOSE(101)
    
    !dr=2*res/1.414214   !resolution=pixel radius. dr=pixel spacing
    !Temporary for check point spread function PSF
    !dr = 2*0.32/1.41424
    dr = res   !052509  JWH - see 050409 note     	

    npix_x=int(m%lx/dr)
    npix_y=int(m%ly/dr)

    write(*,*)"pixels", npix_x, "by", npix_y

    !temporary - 022409
    !dr=13.0
    !npix_x=8
    !npix_y=8
    !temporary- end

    nk = nki
    npix = npix_x*npix_y
    nrot = ntheta*nphi*npsi

    allocate(int_i(nk, npix, nrot), old_int(nk, npix, nrot), old_int_sq(nk, npix, nrot), &
     int_sq(nk, npix, nrot), int_sum(nk), int_sq_sum(nk), stat=istat)
    nullify(old_index, old_pos)

    if (istat /= 0) then
       write (*,*) 'Cannot allocate memory in fem_initialize.'
       return
    endif

    call init_pix(m,npix_x, npix_y, dr, istat)
    if (istat /= 0) return
    call init_rot(ntheta, nphi, npsi, istat)
    if (istat /= 0) return

    call read_f_e 
    allocate(scatfact_e(m%nelements,nk), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Allocation of electron scattering factors table failed.'
       return
    endif

    !do j=1,m%nelements	   
    !   do i=1, nk
    !      scatfact_e(j,i)=0.0    !loop is required to initialize with keeping k intact.  !tr RE-ok-jwh
    !   enddo
    !enddo

    do j=1,m%nelements
       do i=1, nk
          scatfact_e(j,i)=f_e(m%atom_type(j),k(i)) 
       enddo
    enddo
    
  end subroutine fem_initialize


  
  subroutine fem(m, res, k, vk, scatfact_e, comm,  istat, use_femsim, rot_begin, rot_end)
    use mpi
    IMPLICIT NONE

    type(model), intent(in) :: m
    real, intent(in) :: res
    real, dimension(:), intent(in) :: k
    real, dimension(:), INTENT(OUT) :: Vk
    real, dimension(:,:), pointer :: scatfact_e
    integer, intent(out) :: istat
    LOGICAL, OPTIONAL, INTENT(IN) :: use_femsim
    INTEGER, OPTIONAL, INTENT(IN) :: rot_begin, rot_end
    LOGICAL femsim !added by Feng Yi on 03/19/2009
    real, dimension(:), allocatable :: psum_int, psum_int_sq, sum_int, sum_int_sq   !for mpi
    integer :: comm

    integer :: i, j, i1, j1
    INTEGER begin_rot, end_rot
    REAL test1
    
   IF( PRESENT(use_femsim)) THEN
    femsim = use_femsim
   ELSE
    femsim = .FALSE.
   ENDIF

   IF( PRESENT(rot_begin)) THEN
    begin_rot = rot_begin
   ELSE
    begin_rot  = 1
   ENDIF

   IF( PRESENT(rot_end)) THEN
    end_rot = rot_end
   ELSE
    end_rot = nrot
   ENDIF
   IF(femsim) THEN   !This is addeb by Feng Yi on 03/19/2009 only for femsim calculation

    !debug for memory leak
    !rot = 0.0
    !i= 1
    !rot(i,1) =0.0
    !rot(i,2) =0.0
    !rot(i,3) = 0.0 
    !ALLOCATE(mrot(1), STAT=istat)
    !call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m, mrot(1), istat)

    do i=begin_rot, end_rot
      ! initialize the rotated models
       allocate(mrot(1), stat=istat) !debug memory leak
       if (istat /= 0) then
          write(*,*) 'Cannot allocate rotated model array.'
          return
       endif
    !WRITE(*,*) i, 'Before rotate_model' !debug
    !rot(i,1) =0.0
    !rot(i,2) =0.0
    !rot(i,3) = 0.0 !for PSF and resolution test
       call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m, mrot(1), istat) !memory leak
    !WRITE(*,*) 'Atom positions!'
    !WRITE(*,*) 1, mrot(1)%xx(1), mrot(1)%yy(1), mrot(1)%zz(1)
    !WRITE(*,*) 2, mrot(1)%xx(2), mrot(1)%yy(2), mrot(1)%zz(2)

!	WRITE(*,*) i, rot(i,1), rot(i,2), rot(i,3)
       if (istat /= 0) then
          write (*,*) 'Failed to rotate model ',i
          return
       endif
    ! calculate intensities
       do j=1, npix
          call intensity(mrot(1), res, pix(j, 1), pix(j, 2), k, int_i(1:nk, j, i), scatfact_e, istat)
      !test1=sin(exp(j*0.2345)) !test memory leak
       enddo
    !WRITE(*,*) i, 'After intensity cal!' !debug
       !DEALLOCATE mrot(1)
       CALL destroy_model(mrot(1)) !memory leak
       DEALLOCATE(mrot) !memory leak
       !WRITE (*,*) 'rotated ', i !debug
    enddo !end i=1, nrot

      
    int_sq = int_i*int_i
    do i=1, nk
       Vk(i) = (sum(int_sq(i,1:npix,1:nrot)/(npix*nrot))/(sum(int_i(i,1:npix,1:nrot)/(npix*nrot))**2) ) - 1.0
    end do
    !***********************************************************
    
    ELSE

    !*************************************

    allocate(psum_int(size(k)), psum_int_sq(size(k)), sum_int(size(k)), sum_int_sq(size(k)), stat=istat)

    sum_int = 0.0
    sum_int_sq = 0.0
    psum_int = 0.0
    psum_int_sq = 0.0

    ! initialize the rotated models

    allocate(mrot(nrot), stat=istat)
    if (istat /= 0) then
       write(*,*) 'Cannot allocate rotated model array.'
       return
    endif

    do i= myid+1, nrot, numprocs
       call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m, mrot(i), istat)
       if (istat /= 0) then
          write (*,*) 'Failed to rotate model ',i
          return
       endif
    enddo

    allocate(old_index(nrot), old_pos(nrot), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Cannot allocate memory for old indices and positions in fem_initialize.'
       return
    endif

    ! initialize old_index and old_pos arrays
    do i=myid+1, nrot, numprocs
       old_index(i)%nat = 0
       nullify(old_index(i)%ind)
       old_pos(i)%nat = 0
       nullify(old_pos(i)%pos)
    enddo

    ! calculate intensities
    do i=myid+1, nrot, numprocs
       do j=1, npix
          call intensity(mrot(i), res, pix(j, 1), pix(j, 2), k, int_i(1:nk, j, i), scatfact_e, istat)
       int_sq(1:nk, j, i) = int_i(1:nk, j,i)**2

          psum_int(1:nk) = psum_int(1:nk) + int_i(1:nk, j, i)
          psum_int_sq(1:nk) = psum_int_sq(1:nk) + int_sq(1:nk, j, i)
    enddo

    enddo

    call mpi_reduce(psum_int, sum_int , size(k), mpi_real, mpi_sum, 0, comm, mpierr)
    call mpi_reduce(psum_int_sq, sum_int_sq, size(k), mpi_real, mpi_sum, 0, comm, mpierr)
    
    ! calculate the variance
    if(myid.eq.0)then   
    do i=1, nk
        !Vk(i) = ( sum(int_sq(i,1:npix,1:nrot)/(npix*nrot))/(sum(int_i(i,1:npix,1:nrot)/(npix*nrot))**2) ) - 1.0
    Vk(i) = ( sum_int_sq(i)/(npix*nrot) ) / ( ( sum_int(i)/(npix*nrot) )**2 ) - 1.0
    end do
    endif   

    deallocate(psum_int, psum_int_sq, sum_int, sum_int_sq)

   ENDIF



  end subroutine fem

!****************************************

   SUBROUTINE I_average(i_k)
    IMPLICIT NONE

	REAL ,DIMENSION(:), INTENT(OUT) :: i_k
    INTEGER i

    i_k = 0.0
    
    DO i=1, nk
      i_k(i) = SUM(int_i(i,1:npix,1:nrot))/(npix * nrot)
    ENDDO

   END SUBROUTINE I_average

!*******************************
  subroutine fem_update(m_in, atom, res, k, vk, scatfact_e, comm, istat)
    use mpi
    !include 'mpif.h'
    type(model), intent(inout) :: m_in
    integer, intent(in) :: atom
    real, intent(in) :: res
    real, dimension(:), intent(in) :: k
    real, dimension(:), intent(out) :: vk
    real, dimension(:,:), pointer :: scatfact_e
    integer, intent(out) :: istat
    real, dimension(:), allocatable :: psum_int, psum_int_sq, sum_int, sum_int_sq   !for mpi
    integer :: comm



    type(model) :: moved_atom, rot_atom
    integer :: i, j, m, n
    real :: res2, rot_dist_sq, orig_dist_sq, temp1, temp2

    ! res is diameter of a pixel, but later we need the square of the radius
    res2 = (res)**2    !temporary - radius = resolution assumed - JWH 02/25/09

    ! copy the original position of the moved atom in the original orientation into the model
    ! moved atom
    allocate(moved_atom%xx(1), moved_atom%yy(1), moved_atom%zz(1), moved_atom%znum(1), &
         moved_atom%atom_type(1), moved_atom%znum_r(1), moved_atom%composition(1), stat=istat)

    allocate(psum_int(size(k)), psum_int_sq(size(k)), sum_int(size(k)), sum_int_sq(size(k)), stat=istat)

    if(istat /= 0) return
    
    !moved_atom%xx(1) = 0.0   !tr RE-jwh
    !moved_atom%yy(1) = 0.0
    !moved_atom%zz(1) = 0.0

    moved_atom%natoms = 1
    moved_atom%xx(1) = m_in%xx(atom)
    moved_atom%yy(1) = m_in%yy(atom)
    moved_atom%zz(1) = m_in%zz(atom)
    moved_atom%znum(1) = m_in%znum(atom)
    moved_atom%znum_r(1) = m_in%znum_r(atom)
    moved_atom%lx = m_in%lx
    moved_atom%ly = m_in%ly
    moved_atom%lz = m_in%lz
    moved_atom%nelements = 1
    moved_atom%atom_type(1) = m_in%znum(atom)
    moved_atom%composition(1) = 1.0    
    moved_atom%rotated = .FALSE.

    istat = 0

    ! update each rotated model with the position of the moved atom, then recalculate the intensities
    ! of only the pixels involving the moved atom
    
    old_int=0.0
    old_int_sq=0.0

    sum_int = 0.0
    sum_int_sq = 0.0
    psum_int = 0.0
    psum_int_sq = 0.0

    rotations: do i = myid+1, nrot, numprocs
       do m=1, npix   
          old_int(1:nk, m, i) = int_i(1:nk, m, i)
          old_int_sq(1:nk, m, i) = int_sq(1:nk, m, i)
       enddo

!       write (*,*) 'Working on rotation ',i

       ! rotate that moved atom
       call rotate_model(rot(i,1), rot(i, 2), rot(i, 3), moved_atom, rot_atom, istat)
       if(istat /= 0) return

!       if( (rot_atom%natoms == 0) .and. (mrot(i)%rot_i(atom)%nat == 0) ) &
!            write(*,*)"cycle", mrot(i)%rot_i(atom)%nat, rot_atom%natoms

       !if( (rot_atom%natoms == 0) .and. (mrot(i)%rot_i(atom)%nat == 0) )  cycle   !JWH - 04/14/2009 
       if( (rot_atom%natoms == 0) .and. (mrot(i)%rot_i(atom)%nat == 0) )  go to 100

       ! the rotated atom has left the model, so there is no structural
       ! change - skip the rest of the loop
       ! if moving the atom changes the number of times it appears in the rotated model,
       ! just give up: re-rotation the entire model, and recalculate all the intensities
       ! (might not be necessary to recalculate all the intensities.  have to think about that
   
       if(mrot(i)%rot_i(atom)%nat /= rot_atom%natoms) then

!          write(*,*)"scratch", mrot(i)%rot_i(atom)%nat, rot_atom%natoms

          !write (*,*) 'Rerotating the model from scratch.'
          call destroy_model(mrot(i))
          call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m_in, mrot(i), istat)
          old_index(i)%nat = -1   ! -1 is not a valid # of atoms, so use this as a flag in reject_move
                                  ! that the input model must be rerotated from scratch
    
          do m=1, npix
             call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e, istat)
      int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
          enddo

       else

!          write(*,*)"reg", mrot(i)%rot_i(atom)%nat, rot_atom%natoms
          
          ! replace the original, rotated atom with its new position(s).  if it appears more than once
          ! in the new model, the rotated model atoms should appear in the same order as they do in the 
          ! full-size, rotated model
          !write (*,*) 'replacing rotated atoms with new positions.'

          do j=1,mrot(i)%rot_i(atom)%nat

             ! store the original index and position in old_index and old_pos
             call add_index(old_index(i), mrot(i)%rot_i(atom)%ind(j))
             call add_pos(old_pos(i), mrot(i)%xx(mrot(i)%rot_i(atom)%ind(j)), & 
                  mrot(i)%yy(mrot(i)%rot_i(atom)%ind(j)), mrot(i)%zz(mrot(i)%rot_i(atom)%ind(j)), istat)

             ! change the model atom position to the new, rotated position
             mrot(i)%xx(mrot(i)%rot_i(atom)%ind(j)) = rot_atom%xx(j)
             mrot(i)%yy(mrot(i)%rot_i(atom)%ind(j)) = rot_atom%yy(j)
             mrot(i)%zz(mrot(i)%rot_i(atom)%ind(j)) = rot_atom%zz(j)
            
          enddo

          ! now check if either the original position of the moved atom or the rotated position of the moved
          ! atom are inside each pixel.  If so, that pixel intensity must be recalculated.
          ! replace this with one loop over the rotated atom positions and another over the original atom
          ! positions, and if either is true, calculate the intensities.  that way it works for either 
          ! case above.

 
          do m=1, npix

             do n=1, rot_atom%natoms
                !This didn't have PBC. Now it's added below - 030409 - JWH
                !rot_dist_sq = (mrot(i)%xx(n) - pix(m,1))**2 + (mrot(i)%yy(n) - pix(m,2))**2  

                temp1 = rot_atom%xx(n) - pix(m,1)
                temp2 = rot_atom%yy(n) - pix(m,2)
                temp1 = temp1 - mrot(i)%lx*anint(temp1/mrot(i)%lx)
                temp2 = temp2 - mrot(i)%ly*anint(temp2/mrot(i)%ly)

                rot_dist_sq = (temp1)**2 + (temp2)**2

                temp1 = old_pos(i)%pos(n,1) - pix(m,1)
                temp2 = old_pos(i)%pos(n,2) - pix(m,2)
                temp1 = temp1 - mrot(i)%lx*anint(temp1/mrot(i)%lx)
                temp2 = temp2 - mrot(i)%ly*anint(temp2/mrot(i)%ly)

                orig_dist_sq = (temp1)**2 + (temp2)**2

                if ( (rot_dist_sq <= res2) .OR. (orig_dist_sq <= res2) ) then
                   call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e,istat)
            int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                endif
             enddo
          enddo
       endif

100 continue
    do m=1, npix
        psum_int(1:nk) = psum_int(1:nk) + int_i(1:nk, m, i)
        psum_int_sq(1:nk) = psum_int_sq(1:nk) + int_sq(1:nk, m, i)
    enddo
    
    !Deallocate ind in rot_atom%rot_i
    
    DO n=1, SIZE(rot_atom%rot_i)
!	IF (rot_atom%rot_i(j)%nat .GT. 0) THEN
    !write(*,*)"rot_atom%rot_i(j)%nat", rot_atom%rot_i(j)%nat
    if(associated(rot_atom%rot_i(n)%ind))then
      DEALLOCATE(rot_atom%rot_i(n)%ind)
    else
!	write(*,*)"not associated"
    ENDIF
!	endif
    ENDDO
    
    deallocate(rot_atom%xx, rot_atom%yy, rot_atom%zz, rot_atom%znum,  rot_atom%rot_i, rot_atom%znum_r, stat=istat)
    
    enddo rotations
    
    call mpi_reduce(psum_int, sum_int , size(k), mpi_real, mpi_sum, 0, comm, mpierr)
    call mpi_reduce(psum_int_sq, sum_int_sq, size(k), mpi_real, mpi_sum, 0, comm, mpierr)
    
    !write(*,*)myid, sum_int 

    ! recalculate the variance
    if(myid.eq.0)then
    do i=1, nk
        !Vk(i) = ( sum(int_sq(i,1:npix,1:nrot)/(npix*nrot))/(sum(int_i(i,1:npix,1:nrot)/(npix*nrot))**2) ) - 1.0
    Vk(i) = ( sum_int_sq(i)/(npix*nrot) ) / ( ( sum_int(i)/(npix*nrot) )**2 ) - 1.0
    end do
    endif   
    
    deallocate(moved_atom%xx, moved_atom%yy, moved_atom%zz, moved_atom%znum, &
         moved_atom%atom_type, moved_atom%znum_r, moved_atom%composition, stat=istat)
    deallocate(psum_int, psum_int_sq, sum_int, sum_int_sq)
    

  end subroutine fem_update








  ! accept the move.  don't need to change any of the atom positions in any of the rotated
  ! models, but we do need to clear the old_index and old_pos arrays for reuse.
  subroutine fem_accept_move(comm)

    use mpi
    integer :: comm
    call fem_reset_old(comm)

  end subroutine fem_accept_move



  ! reject the move.  replace all the moved atoms with their original positions.  requires
  ! that the model argument have the original (unmoved) positions in it.
  subroutine fem_reject_move(m, comm)

    use mpi
    type(model), intent(inout) :: m
    integer :: i, j, istat
    integer :: comm


    do i=myid+1, nrot, numprocs

       ! the rotated atom wasn't in the model, so this model doesn't need to be changed
       if(old_index(i)%nat == 0) cycle  

       ! the move changed the number of atoms in the model, so the model must be re-rotated
       ! from scratch
       if(old_index(i)%nat == -1) then
          call destroy_model(mrot(i))
          call rotate_model(rot(i,1), rot(i,2), rot(i,3), m, mrot(i), istat)
          cycle
       endif

       ! otherwise, copy the old positions back into the model at the correct indices

       do j=1,old_index(i)%nat
          mrot(i)%xx(old_index(i)%ind(j)) = old_pos(i)%pos(j,1)
          mrot(i)%yy(old_index(i)%ind(j)) = old_pos(i)%pos(j,2)
          mrot(i)%zz(old_index(i)%ind(j)) = old_pos(i)%pos(j,3)
       enddo

       !The saved intensity values must return to their old values - JWH 03/05/09
       do j=1, npix  
          int_i(1:nk, j, i) = old_int(1:nk, j, i)
          int_sq(1:nk, j, i) = old_int_sq(1:nk, j, i)
       enddo
    enddo

    call fem_reset_old(comm)

  end subroutine fem_reject_move



  subroutine write_intensities(outfile, k, istat)
    character(len=*), intent(in) :: outfile
    real, dimension(:), intent(in) :: k
    integer, intent(out) :: istat

    integer :: ik, irot, ipix

    istat = 0

    open(unit=314, file=outfile, form='formatted',status='replace',iostat=istat)
    if( istat /= 0) then
       write (*,*) 'Cannot open intensities output file: ',istat
       return
    endif


    do ik=1, nk
       if(ik /= 1) write (314,*) ' '
       write (314,*) 'k =',k(ik)
       write (314,*) 'theta  phi   psi  pix_x   pix_y  intensity'
       do irot=1, nrot
          do ipix=1, npix
             write (314,'(G16.8,G16.8,G16.8,G16.8,G16.8,G16.8)') rot(irot, 1), rot(irot, 2), & 
                  rot(irot, 3), pix(ipix, 1), pix(ipix, 2), int_i(ik, ipix, irot)
          enddo
       enddo
    enddo

    close(314)

  end subroutine write_intensities



  subroutine print_image1
    integer :: i, j, kk

    open(unit=903,file="test_img_intensity1.txt",form='formatted',status='unknown')
    write(903,*)'rotation'," ",'pixel'," ",'intensity'
    
    kk = 6.0   !nth k point
    i = 2.0       !nt rotation
 
    do j=1, npix   !pixel loop
       write(903,*)i, j, int_i(kk, j, i)
    enddo
    write(903,*)
    close(903)
 
  end subroutine print_image1



  subroutine print_image2(m, res, k, scatfact_e, istat)
    type(model), intent(in) :: m
    real, intent(in) :: res
    real, dimension(:), intent(in) :: k
	real, dimension(:,:), pointer :: scatfact_e
    integer, intent(out) :: istat
    integer :: i, j

    ! initialize the rotated models
    allocate(mrot(nrot), stat=istat)
    if (istat /= 0) then
       write(*,*) 'Cannot allocate rotated model array.'
       return
    endif

    do i=1, nrot
       call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m, mrot(i), istat)
       if (istat /= 0) then
          write (*,*) 'Failed to rotate model ',i
          return
       endif
    enddo

    allocate(old_index(nrot), old_pos(nrot), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Cannot allocate memory for old indices and positions in fem_initialize.'
       return
    endif

    ! initialize old_index and old_pos arrays
    do i=1, nrot
       old_index(i)%nat = 0
       nullify(old_index(i)%ind)
       old_pos(i)%nat = 0
       nullify(old_pos(i)%pos)
    enddo

    open(unit=902,file="test_img_intensity_initial.txt",form='formatted',status='unknown')
    write(902,*)'rotation', 'pixel', 'intensity'
    ! calculate intensities
    do i=1, nrot
       do j=1, npix
          call intensity(mrot(i), res, pix(j, 1), pix(j, 2), k, int_i(1:nk, j, i), scatfact_e,istat)
          write(902,*)i, j,pix(j,1),pix(j,2),int_i(6, j, i)
       enddo
    enddo


    stop
  end subroutine print_image2

!
! private subroutines
!

  subroutine init_pix(m, npix_x, npix_y, dr, istat)
    type(model), intent(in) :: m
    integer, intent(in) :: npix_x, npix_y
    real dr
    integer, intent(out) :: istat
    integer i, j, k

    allocate(pix(npix, 2), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Cannot allocate pixel position array.'
       return
    endif

    k=1
    do i=1, npix_x
       do j=1, npix_y
          pix(k,1) = (i-1)*dr+(dr/2.0)-(m%lx/2.0)
          pix(k,2) = (j-1)*dr+(dr/2.0)-(m%ly/2.0)
!          write(*,*)pix(k,1),pix(k,2)
          k = k+1
       enddo
    enddo

    !temporary - to compare with old t11
    !k=1
    !do i=1, npix_x
     !  do j=1, npix_y
      !    pix(k,1) = -45.5 + 13.0*(i-1)
       !   pix(k,2) = -45.5 + 13.0*(j-1)
        !  k=k+1
       !enddo
    !enddo

  end subroutine init_pix


  subroutine init_rot(ntheta, nphi, npsi, istat)
    integer, intent(in) :: ntheta, nphi, npsi
    integer, intent(out) :: istat

    integer :: i,j, k, l
	real,dimension(3) :: step_size

    allocate(rot(nrot, 3), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Cannot allocate rotations array.'
       return
    endif

    !phi runs from 0 to 2 PI
    !theta runs from 0 to PI
    !psi runs from 0 to 2 PI
    
    !step_size(1) for theta step
    !step_size(2) for phi step
    !step_size(3) for psi step
    
    step_size(1) = PI / (ntheta + 1.0)
    step_size(2) = TWOPI / nphi
    step_size(3) = TWOPI / npsi
    
    l = 1
    do i=1,ntheta
       do j=1, nphi
          do k=1, npsi
             rot(l, 1) = step_size(1) * i      !real(ntheta)
             rot(l, 2) = step_size(2) * (j-1)  !real(nphi)
             rot(l, 3) = step_size(3) * (k-1)  !real(npsi)
             l = l+1
          enddo
       enddo
    enddo

    !temporary 02/19/09 - Generates same rotations as the old RMC (for comparison)
    !rot(1,1)=0.0
    !rot(1,2)=0.0
    !rot(1,3)=0.0
    !rot(2,1)=pi/2.0
    !rot(2,2)=0.0
    !rot(2,3)=0.0
    !rot(3,1)=pi/2.0
    !rot(3,2)=pi/2.0
    !rot(3,3)=0.0
    !temporary 02/19/09 - end
   
  end subroutine init_rot




  subroutine intensity(m_int, res, px, py, k, int_i, scatfact_e, istat)
    type(model), intent(in) :: m_int
    real, intent(in) :: res, px, py
    real, dimension(nk), intent(in) :: k
    real, dimension(nk), intent(out) :: int_i
    !real, dimension(:,:), pointer, intent(in) :: scatfact_e
    real, dimension(:,:), pointer :: scatfact_e
    integer, intent(out) :: istat

    real, dimension(:,:,:), allocatable :: gr_i   ! unneeded 'save' keyword removed pmv 03/18/09  !tr RE-ok -jwh  
    real, dimension(:), allocatable ::x1, y1, rr_a
    real, dimension(:,:), allocatable :: sum1
    real :: x2, y2, rr, t1, t2, const1, const2, const3, pp, r_max 
    integer, pointer, dimension(:) :: pix_atoms, znum_r
    integer :: i,j,ii,jj,kk
    integer :: bin_max

    r_max = 2*res     !assuming resolution=radius
    bin_max = int(r_max/fem_bin_width)+1

    call hutch_list_pixel(m_int, px, py, res, pix_atoms, istat)

    allocate(gr_i(m_int%nelements,m_int%nelements, 0:bin_max), stat=istat)
    allocate(x1(size(pix_atoms)),y1(size(pix_atoms)),rr_a(size(pix_atoms)), stat=istat)
    allocate(sum1(m_int%nelements,size(pix_atoms)), stat=istat)
    allocate(znum_r(size(pix_atoms)), stat=istat)


    ! replaced code recalculating znum_r with code copying it from previous calculations 3/18/09 pmv  !tr RE-ok-jwh
    !WRITE(*,*) 'before using znum_r'
    do i=1, size(pix_atoms)
       znum_r(i) = m_int%znum_r(pix_atoms(i))
    enddo
    !WRITE(*,*) 'After assigning znum_r' !debug
    !do i=1, size(pix_atoms)
    !   do j=1, m_int%nelements
    !      if( m_int%znum(pix_atoms(i)) .eq. m_int%atom_type(j) )then
    !         znum_r(i) = j
    !      endif
    !   enddo
    !enddo

    gr_i = 0.0
    int_i = 0.0
    x1 = 0.0
    y1 = 0.0
    rr_a = 0.0
    x2 = 0.0
    y2 = 0.0
    
    const1 = twopi*(0.61/res)/fem_bin_width  !(0.61/res = Q) 
    const2 = 1/fem_bin_width
    const3 = (const1/(0.61/res))/const2
    
    do i=1,size(pix_atoms)
       x2=m_int%xx(pix_atoms(i))-px           
       y2=m_int%yy(pix_atoms(i))-py
       x2=x2-m_int%lx*anint(x2/m_int%lx)       
       y2=y2-m_int%ly*anint(y2/m_int%ly)
       rr_a(i)=sqrt(x2*x2 + y2*y2)
       if(rr_a(i).le.res)then
          x1(i)=x2
          y1(i)=y2
          j=int(const1*rr_a(i))
          sum1(znum_r(i),i)=A1(j)

!          if(a1(j).ne.a1(j))then !debug - jwh
!             write(*,*)j,a1(j)
!             pause
!          endif

       endif
    enddo
    !WRITE(*,*) 'Finish sum1 in intensity calc' !debug 
    do i=1,size(pix_atoms)
       if(rr_a(i).le.res)then
          do j=i,size(pix_atoms)
             if(rr_a(j).le.res)then
                x2=x1(i)-x1(j)
                y2=y1(i)-y1(j)
                rr=sqrt(x2*x2 + y2*y2)
                kk=int(const2*rr)
        IF(kk .GE. bin_max) THEN !debug
          WRITE(*,*) 'kk is larger than bin_max', kk, bin_max
        ENDIF
                if(i == j)then
                   t1=sum1(znum_r(i),i)
                   gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+t1*t1 
                else
                   t1=sum1(znum_r(i),i)
                   t2=sum1(znum_r(j),j)
                   gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+2.0*t1*t2  !changed by FY on 05/04/2009
                end if
                
!                if(gr_i(znum_r(i),znum_r(j),kk).ne.gr_i(znum_r(i),znum_r(j),kk))then
!                  write(*,*)t1,t2,gr_i(znum_r(i),znum_r(j),kk)
!                  pause
!                endif
             endif
          end do
       endif
    end do
    
    !open(111, file='nan_test_print.txt', status='unknown')
    !WRITE(*,*) 'Finish gr_i calc in intensity!'  !debug
    do i=1,nk
       do j=0,bin_max
          do ii=1,m_int%nelements
             do jj=1,m_int%nelements
                !write(111,*)gr_i(ii,jj,j)
                pp=const3*j*k(i)
                int_i(i)=int_i(i)+scatfact_e(ii,i)*scatfact_e(jj,i)*J0(INT(pp))*gr_i(ii,jj,j)
             enddo
          enddo
       end do
    end do
   !WRITE(*,*) 'Finish int_i calcualtion!'
   IF(ALLOCATED(gr_i)) THEN
     DEALLOCATE(gr_i)
   ENDIF
   !WRITE(*,*) 'gr_i deall!'
   IF(ALLOCATED(x1)) THEN
     DEALLOCATE(x1,y1, rr_a, znum_r)
   ENDIF
   !WRITE(*,*) 'x1, y1, rr_a, znum_r'
   IF(SIZE(pix_atoms) .GT. 0) THEN
    DEALLOCATE(pix_atoms)
   ENDIF
   !WRITE(*,*) 'pix_atoms deal', 'sum1 size: ', SIZE(sum1, 1), SIZE(sum1, 2)
   IF(ALLOCATED(sum1)) THEN
     DEALLOCATE(sum1)
   ENDIF
   !WRITE(*,*) 'sum1 deal'
   
    !deallocate(gr_i, x1, y1, sum1, znum_r, rr_a, pix_atoms)
  end subroutine intensity



  subroutine add_pos(p, xx, yy, zz, istat)
    type(pos_list), intent(inout) :: p
    real, intent(in) :: xx, yy,  zz
    integer, intent(out) :: istat

    
    real, dimension(:,:), allocatable :: scratch
    
    !if (associated(p%pos)) then
    if (p%nat .GT. 0) then
       allocate(scratch(p%nat+1,3), stat=istat)
       if (istat /= 0) continue
       scratch(1:p%nat, 1:3) = p%pos
       p%nat = p%nat+1
       scratch(p%nat,1) = xx
       scratch(p%nat,2) = yy
       scratch(p%nat,3) = zz
       deallocate(p%pos)
       allocate(p%pos(p%nat,3), stat=istat)
       if (istat /= 0) continue
       p%pos = scratch
    else
       p%nat = 1
       allocate(p%pos(1,3), stat=istat)
       if (istat /= 0) continue
       p%pos(1,1) = xx
       p%pos(1,2) = yy
       p%pos(1,3) = zz
    endif

    if (istat /= 0) then
       write (*,*) 'Error allocating memory in add_pos.'
    endif
 
    if (allocated(scratch)) then
    deallocate(scratch)  ! added 3/18/09 pmv 
    endif

  end subroutine add_pos


  subroutine destroy_pos_list(p, np)
    type(pos_list), pointer, dimension(:) :: p
    integer, intent(in) :: np

    integer i
    do i=1,np
       deallocate(p(i)%pos)
    enddo

    deallocate(p)

  end subroutine destroy_pos_list


  
  subroutine fem_reset_old(comm)
   
    use mpi
    integer :: comm
    integer :: i

    do i=myid+1, nrot, numprocs
       !if(old_index(i)%nat /= 0) deallocate(old_index(i)%ind)  ! deallocate lists added 3/18/09, pmv  !tr-jwh
    if(old_index(i)%nat > 0) deallocate(old_index(i)%ind)
       nullify(old_index(i)%ind)
       old_index(i)%nat = 0
       !if(old_pos(i)%nat /= 0) deallocate(old_pos(i)%pos) ! deallocate added 3/18/09 pmv  !tr-jwh
    if(old_pos(i)%nat > 0) deallocate(old_pos(i)%pos)
       nullify(old_pos(i)%pos)
       old_pos(i)%nat = 0
    enddo

  end subroutine fem_reset_old



  subroutine bessel_func(x,bj0,bj1)
 
    IMPLICIT none
    doubleprecision A,B,A1,B1,BJ0,BJ1,BY0,BY1,DY0,DY1,X,X2,RP2,DJ0,DJ1,R,DABS
    doubleprecision EC,CS0,W0,R0,CS1,W1,R1,T1,T2,P0,P1,Q0,Q1,CU,DCOS,DSIN
    !integer i,j,k,l,m,n,k0
    integer k,k0
    DIMENSION A(12),B(12),A1(12),B1(12)
   
        RP2=0.63661977236758D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN
           BJ0=1.0D0
           R=1.0D0
           DO 5 K=1,30
              R=-0.25D0*R*X2/(K*K)
              BJ0=BJ0+R
              IF (DABS(R).LT.DABS(BJ0)*1.0D-15) GO TO 10
5          CONTINUE
10         BJ1=1.0D0
           R=1.0D0
           DO 15 K=1,30
              R=-0.25D0*R*X2/(K*(K+1.0D0))
              BJ1=BJ1+R
              IF (DABS(R).LT.DABS(BJ1)*1.0D-15) GO TO 20
15         CONTINUE
20         BJ1=0.5D0*X*BJ1
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           CS0=0.0D0
           W0=0.0D0
           R0=1.0D0
           DO 25 K=1,30
              W0=W0+1.0D0/K
              R0=-0.25D0*R0/(K*K)*X2
              R=R0*W0
              CS0=CS0+R
              IF (DABS(R).LT.DABS(CS0)*1.0D-15) GO TO 30
25         CONTINUE
30         BY0=RP2*(EC*BJ0-CS0)
           CS1=1.0D0
           W1=0.0D0
           R1=1.0D0
           DO 35 K=1,30
              W1=W1+1.0D0/K
              R1=-0.25D0*R1/(K*(K+1))*X2
              R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS1=CS1+R
              IF (DABS(R).LT.DABS(CS1)*1.0D-15) GO TO 40
35         CONTINUE
40         BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00, &
            -.5725014209747314D+00,.6074042001273483D+01, &
            -.1100171402692467D+03,.3038090510922384D+04, &
            -.1188384262567832D+06,.6252951493434797D+07, &
            -.4259392165047669D+09,.3646840080706556D+11, &
            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00, &
            .1727727502584457D+01,-.2438052969955606D+02, &
            .5513358961220206D+03,-.1825775547429318D+05, &
            .8328593040162893D+06,-.5006958953198893D+08, &
            .3836255180230433D+10,-.3649010818849833D+12, &
            .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00, &
            .6765925884246826D+00,-.6883914268109947D+01, &
            .1215978918765359D+03,-.3302272294480852D+04, &
            .1276412726461746D+06,-.6656367718817688D+07, &
            .4502786003050393D+09,-.3833857520742790D+11, &
            .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00, &
            -.1993531733751297D+01,.2724882731126854D+02, &
            -.6038440767050702D+03,.1971837591223663D+05, &
            -.8902978767070678D+06,.5310411010968522D+08, &
            -.4043620325107754D+10,.3827011346598605D+12, &
            -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO 45 K=1,K0
              P0=P0+A(K)*X**(-2*K)
45            Q0=Q0+B(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO 50 K=1,K0
              P1=P1+A1(K)*X**(-2*K)
50            Q1=Q1+B1(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN

  end subroutine bessel_func

end module fem_mod
