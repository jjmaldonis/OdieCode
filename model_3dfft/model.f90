! Model module: data types and functions for manipulating model structures
!
! Parameters that must be defined elsewhere:
!   ATOMS_PER_HUTCH: the average # of atoms in each hutch
!
! Defined data types:
!
!  model, consisting of:
!   xx, yy, zz: coordinate lists in Angstroms
!   znum: atomic number list
!   lx, ly, lz: box size in Angstroms
!   nelements: # of elements
!   atom_type: list of atomic numbers present
!   composition: list of fractional compositions, in the order of atom_type
!   ha: hutch array data structure
!   logical rotated: TRUE if model has been rotated, FALSE otherwise
!   unrot_atoms: number of atoms in the model before it was rotated
!   rot_i: array of index_lists; lists which atom(s) atom i has become in a rotated model
!   
! hutch:
!   at: list of atom indices inside the hutch
!   nat: # of atoms inside the hutch
!
! hutch_array:
!   h: 3D array of hutches, covering the model
!   nhutch_x, nhutch_y, nhutch_z: # of hutches in x, y, and z
!   hutch_size: physical size of one (cubic) hutch in Angstroms
!   atom_hutch: list of which hutch (i,j,k) a particular atom belongs to
!
! index_list:
!   ind: list of atom indices
!   nat: # of atoms in the list
! INDEX_LIST IS STUPID, SINCE IT IS THE SAME AS HUTCH.  SHOULD CONSOLIDATE.
!
! hutch_2D_list:
!    stores relative position (only in x and y dimension)
!    with respect to a square, this is for 2D calculation
!
! hutch_2D_array:
!    2-D array of hutch_2D_list type
!    first index refers to relative position of a point in a square
!    second index refers to ratio of radius to square size
!    position=relative position of a point in a square 
!    position=j*(number_x-1) + i, where j refers to which y position the point is in
!    i refers to which x position the point is in
!    size_position = position dimension size
!    ratio_radius_square = ratio array
!    size_ratio = ratio dimension size
!
! hutch_3D_list:
!    store relative position (only in x ,y, and z dimension)
!    with respect to a hutch, this is for 3D calculation
!
! hutch_3D_array stores 
!    relative position in x and y dimension hutch_size as the unit for counting
!    list_3D is a 3-D array of hutch_3D_list type
!    first index refers to relative position of a point in a hutch
!    second index refers to ratio of radius to hutch size
!    position=relative position of a point in a hutch 
!    position=k*(number_y-1)*(number_x-1)+j*(number_x-1) + i, 
!    where k refers to which z position the point is in
!    j refers to which y position the point is in
!    i refers to which x position the point is in
!    size_position = position dimension size
!    ratio_radius_square = ratio array
!    size_ratio = ratio dimension size
!
!
! Subroutines:
!
! read_model(input_filename, comment, model, status)
!    Reads the file input_filename into the model type model.  The comment line from the model
!    ends up in comment, and status is nonzero if the file won't open or the memory can't be allocated.
!
! write_model(output_filename, comment, model, status)
!    Dumps the model in model to the file output_filename, with first comment line comment. 
!    Returns nonzero status if the file won't open.
!
! composition_model(model)
!  Calculates the number of different elements in the model, sets nelements, and fills in 
!  the lists atom_type and composition.
!
! recenter_model(xc, yc, zc, model)
!    Shifts the atom positions in model so that the mid-point between the maximum and 
!    and minimum atom positions in each dimensions sits at the position (xc, yc, zc),
!    measured in units of the model supercell.
!
! check_model(model, status)
!    Simple error checking on the model.  Currently makes sure that there are no atoms outside 
!    the box and that all the atomic numbers are within the range 1 to 103.  Returns nonzero in status
!    if there's a problem.
!
! periodic_continue_model(xp, yp, zp, model_in, model_out, initialize_hutches, status)
!    Uses periodic boundary conditions to make (xp, yp, zp) copies of the input model in the
!    output model.  If initialize_hutches is .TRUE., the hutch array is initialized.  If it's .FALSE.
!    it's not initialized.  Returns nonzero in status if there's a memory allocation problem.
!    
! rotate_model(theta, phi, psi, model_in, model_out, istat)
!    Takes the input model, expands it by 3x3x3, rotates it through the Euler angles (theta, phi, 
!    psi), then cuts it back down to size and puts the results in model_out.  The output model
!    no longer obeys periodic  boundary conditions.  istat is nonzero for memory allocation failure
!    
! copy_model(model_in, model_out, status)
!    Copies the input model to the output model.  Returns nonzero in istat if memory can't be
!    allocated.
!
! destroy_model(model)
!    deallocates all the dynamically allocated / pointer etc bits of a model
!
! add_index(index_list, index)
!    adds the integer index to the end of the ind array in index_list and updates the # of 
!    indices.
! THIS FUNCTIONS IS STUPID - IT DUPLICATE HUTCH_ADD_ATOM BELOW.
!
! destroy_rot_indices(unrot_atoms, ri)
!    Deallocates all the of the index lists in ri.  probably duplicates destroy_hutch.
!
! model_init_hutch(model, status)
!    Populates the hutch array with appropriate values based on current model the global atom 
!    positions, cell parameters, etc.  Does NOT check to see if hutch_array has already been 
!    initialized, so it should be called only once for each hutch array.  Returns 0 in status 
!    if all is well, and 1 if one of the memory allocations failed.
!
! hutch_move_atom(model, atom, new_x, new_y, new_z)
!    moves the atom with index atom from it's current position in hutch_array to the hutch 
!    corresponding to the position (new_x, new_y, new_z) in Angstroms.  Does no error checking,
!    so it will crash if the new position is outside the supercell, for example.  Useful
!    for updating the hutch array to reflect a Monte Carlo atom move.
!
! hutch_list_3D(model, px, py, pz, radius, atom_indices, status)
!    creates a  list of atom indices in atom_indices of all the atoms that are within a distance
!    radius of the position (px, py, pz) in hutch_array, plus a few more.  Useful for
!    calculated G(r) out to some radius r.  What it actually returns is the list of
!    atoms in a cube in which a sphere of radius radius could be inscribed.  This is
!    to avoid calculating all the distances between atoms, which the G(r) routine is
!    going to do anyway.  Returns 1 in status if atom_indices can't be allocated.
!
! hutch_list_pixel(model, px, py, diameter, atom_indices, status)
!    creates a list of atom indices in atom_indices of all the atoms in a column along z within
!    a diameter diameter centered at the (x,y) position (px,py) in the hutch array, 
!    plus a few more.  Useful for calculating the FEM intensity from a pixel at (px,
!    py).  What it actually returns is the list of atoms in a rectangular prism in which
!    a cylinder of diameter diameter and length lz could be inscribed.  Returns 1 in
!    status if atom_indices can't be allocated.
!
! hutch_position(model, xx, yy, zz, hx, hy, hz)
!    calculates the hutch indices hx, hy, hz that correspond to the 3D position (xx, yy, zz)
!
! hutch_add_atom(model, atom_index, hx, hy, hz)
!    adds an atom atom_index to the list of the (hx, hy, hz) hutch.  called by model_init_hutch
!    and hutch_move_atom
!
! hutch_remove_atom(model, atom_index)
!    removes the atom atom_index the from whatever hutch it currently occupies.  called by
!    hutch_move_atom
!
! destroy_hutch(hutch_array): deallocates all the dynamically allocated / pointer etc
!    arrays in hutch_array
!
! dump_hutches(model): writes a formatted list of which atoms are in each hutch
!    and a list of the hutch for each atom.  Useful for debugging.
!
! pre_calc_2D_hutch:
!    pre-calculates the hutch_2D_array for the improved hutch algorithm in hutch_list_pixel.  
!    Each 2D array each element stores the relative position of each hutch (only in x and y dimension)
!
! destroy_hutch_2D_array:
!    deallocates memory for hutch_2D_array
!
! pre_calc_3D_hutch:
!    pre-calculates the hutch_3D_array for the improved hutch algorithm in hutch_list_3D.
!    each element stores the relative position of each hutch (x, y and z dimension
!
! destroy_hutch_3D_array:
!    deallocates the hutch_3D_array memory created by pre_calc_3D_hutch
!
! hutch_position_eff(m, xx, yy, zz, hx, hy, hz,p_relative_3D, p_relative_2D)
!    Gives the relative position of point (xx, yy,zz) in hutch (hx,hy,hz)
!    in both 2D and 3D dimension
!
! hutch_list_3D_eff:
!    similar to hutch_list_3D, but should list fewer atoms in atom list
!
! hutch_list_pixel_eff:
!    similar to hutch_list_pixel, but should list fewer atoms in atom list
!
!
! Changelog
!
! from old hutch_mod module:
!    12-19-08 first version, pmv
!    atom_hutch in type hutch changed from allocatable to pointer, 12/20/08, fy
!    hutch_list routines modified to make only one loop through the hutches,
!       and to do it in row-major order.  should be faster.   12/29/08, pmv
!
! model module:
!   Written by Jinwoo Hwang 12/19/2008
!   turned into a module by pmv 12-20-08
!   added intent to subroutine arguments, 12/20/08, pmv
!   rewrote to accomodate model defined type and added everything below read & write model, 1/1/09, pmv
!   combined with hutch module, rewrote the hutch functions to use the model dervied type 1/3/09, pmv
!   added rot_i to model dervied type, rewrote rotate_model to use it, for fem_update 1/5/09 pmv
!   fixed pointer assignment bug in hutch_list_pixel and hutch_list_3D /1/6/09 pmv
!   minor change in data types to make it compatible with visual fortran. 1/7/09 Jinwoo Hwang
!   add subroutine Calc_Type_Order and a varialbe Type_Order, by Feng Yi on 01/16/2009 
!      Puts the atomic number in increasing order, totally nelements  elements
!      in an 1-D array of size nelements
!   Remove Calc_Type_Order(m, Type_Order) and move its function to composition_model.  nelements
!      and composition are now sorted by znum, low to high; by Feng Yi on 01/21/2009
!   data types hutch_2d_list, hutch_2D_array, hutch_3D_list and hutch_3D_array,
!      subroutines pre_calc_2D_hutch, destroy_hutch_2D_array, pre_calc_3D_hutch,
!      destroy_hutch_3D_array, hutch_position_eff, hutch_list_3D_eff, hutch_list_pixel_eff, all
!      added by Feng Yi on 01/22/2009
!

module model_mod

  !use RMC_Global  ! Global variables

  implicit none

  ! derived data type for the hutches.  at is a pointer to a list of the indices
  ! of the atoms in this hutch.  nat is the number of atoms in this hutch
  type hutch
     integer, dimension(:), pointer :: at
     integer :: nat
  end type hutch
  
  ! derived data type for the hutch array: contains an array of hutches, plus
  ! supporting information
  type hutch_array
     type(hutch), dimension(:,:,:), pointer :: h     ! array of hutch objects
     integer :: nhutch_x, nhutch_y, nhutch_z             ! number of hutches in x, y, and z
     real :: hutch_size                                  ! physical size of a hutch in Angstroms
     integer, pointer, dimension(:,:) :: atom_hutch  ! list of the hutch indices for every atom
  end type hutch_array

  ! adjustable-size list of atom indices
  type index_list
     integer :: nat
     integer, pointer, dimension(:) :: ind
  end type index_list

  ! list for a rotated model that is the original model's natoms long.  each item in the
  ! list is itself a list of the indices in the new, rotated model of the atoms corresponding
  ! to the the original atom in the unrotated model.
  ! for a rotated model, a list of the indices in the rotated model corresponding to each atom
  ! index in the original, unrotated model.  the original, unrotated model's index long.

  ! Defined type for a structural model with atoms positions and a bunch of metadata
  type model
     integer :: natoms                              ! number of atoms in the model
     real, pointer, dimension(:) :: xx, yy, zz      ! atom positions in Angstroms
     integer, pointer, dimension(:) :: znum         ! atom atomic numbers
     real :: lx, ly, lz                             ! box size, in Angstroms
     integer :: nelements                           ! # of elements in the model
     integer, pointer, dimension(:) :: atom_type    ! array listing atomic numbers present
     real, pointer, dimension(:) :: composition     ! fractional composition in the order of atom_type
     type(hutch_array) :: ha                        ! hutch data structure
     logical :: rotated                             ! TRUE if model has been rotated, FALSE otherwise
     integer :: unrot_natoms
     type(index_list), dimension(:), pointer :: rot_i ! list of which atoms in the rotated model correspond
                                                      ! to the index i in the unrotated model
  end type model


  ! list the relative x, and y position of squre to the center square
  ! for the atom at a certain position in a squre
  ! with certain ratio of radius to the square size
  TYPE hutch_2D_list
     INTEGER, DIMENSION(:), POINTER :: list_x
     INTEGER, DIMENSION(:), POINTER :: list_y
     INTEGER size_d !number of elements in x, y dimension
  END TYPE hutch_2D_list


  ! list the relative x, and y position of square to the center square
  ! divide the square into equal 4 parts
  ! the ratio of radius to the square size range from 0 to 10.5, can be changed
  ! depending on calculation needed
  TYPE hutch_2D_array
     TYPE(hutch_2D_list), DIMENSION(:,:), POINTER :: list_2D
     !first index refers to relative position of a point in a square
     !second index refers to ratio of radius to square size
     !position=relative position of a point in a square 
     !position=j*(number_x-1) + i, where j refers to which y position the point is in
     !i refers to which x position the point is in
     INTEGER size_position
     REAL, DIMENSION(:), POINTER :: ratio_radius_square
     INTEGER size_ratio
  END TYPE hutch_2D_array
  
  
  ! list the relative x, y and Z position of squre to the center hutch
  ! for the atom at a certain position in a hutch
  ! with certain ratio of radius to the hutch size
  TYPE hutch_3D_list
     INTEGER, DIMENSION(:), POINTER :: list_x
     INTEGER, DIMENSION(:), POINTER :: list_y
     INTEGER, DIMENSION(:), POINTER :: list_z
     INTEGER size_d !number of elements in x, y and z dimension
  END TYPE hutch_3D_list
  
  ! list the relative x, and y position of hutch to the center square
  ! divide the square into equal 8 parts
  ! the ratio of radius to the hutch size range from 0 to 10.5, can be changed
  ! depending on calculation needed
  TYPE hutch_3D_array
     TYPE(hutch_3D_list), DIMENSION(:, :), POINTER :: list_3D
     !first index refers to relative position of a point in a hutch
     !second index refers to ratio of radius to hutch size
     !position=relative position of a point in a hutch 
     !position=k*(number_y-1)*(number_x-1)+j*(number_x-1) + i, 
     !where k refers to which z position the point is in
     !j refers to which y position the point is in
     !i refers to which x position the point is in
     INTEGER size_position
     REAL, DIMENSION(:), POINTER :: ratio_radius_hutch
     INTEGER size_ratio
  END TYPE hutch_3D_array

  integer, parameter :: ATOMS_PER_HUTCH=1
    
contains

  ! Reads a model in the Kirkland .xyz file format from the file model_filename.
  ! Puts the first-line comment in "comment", the model in m, and returns 0 in
  ! istat if the file can't be opened or the memory allocation fails.
  subroutine read_model(model_filename, comment, m, istat)

    implicit none

    character (len=*),intent(in) :: model_filename 
    character (LEN=*),intent(out) :: comment
    type(model), intent(out) :: m
    integer, intent(out) :: istat      !0 for successful open, others for failure.
    
    integer i
	
    integer atom_count, nat
    
    nat=0
    
    open(1,file=model_filename,iostat=istat,status='old')
    if(istat.ne.0) then !Open fails
       write(*,*)"Error in opening flie,"," ",model_filename," status=",istat
       return
    endif
    
    atom_count = 0
    read(1,*)
    read(1,*)
    do while(atom_count.ne.-1)
       read(1,*)atom_count
       nat=nat+1.0
    enddo
    nat=nat-1.0

    rewind(1)

    m%natoms = nat
    allocate(m%xx(nat), m%yy(nat), m%zz(nat), m%znum(nat), stat=istat)
    
    if(istat /= 0) then
       write (*,*) 'Unable to allocate memory for the model being read.'
       return
    endif

    read(1,'(a80)') comment
    read(1,*) m%lx,m%ly,m%lz
    do i=1,nat
       read(1,*) m%znum(i),m%xx(i),m%yy(i),m%zz(i)
    enddo
    close(1)

    call composition_model(m)
    m%rotated = 0

    call recenter_model(0.0, 0.0, 0.0, m)
    call model_init_hutch(m, istat)

    write (*,*) 'Read ',nat,' atoms from the model described as:'
    write (*,*) trim(comment)

  end subroutine read_model


  
  
  ! Writes the model in the data structure m to the file output_filename with
  ! the first-line comment "comment".  istat is 0 if the file can't be written.
  subroutine write_model(output_filename, comment, m, istat)
    character (LEN=*),intent(in):: output_filename
    character (LEN=*), intent(in):: comment
    type(model), intent(in) :: m
    integer, intent(out) :: istat

    integer i

    open(2,file=output_filename,iostat=istat,status='unknown')
    
    if(istat.ne.0) then !Open fails
       write(*,*)"Error in opening flie,"," ",output_filename,"status3=",istat 
    endif
    
    write(2,'(a80)') trim(comment)
    write(2,*) m%lx, m%ly, m%lz
    do i=1,m%natoms
       write(2,*) m%znum(i),m%xx(i),m%yy(i),m%zz(i)
    enddo
    write(2,*)"-1"
    close(2)
    
  end subroutine write_model


  ! Calculates the composition of the model and fills in nelements, atom_type, and composition
  subroutine composition_model(m)
    type(model), intent(inout) :: m

    integer, dimension(103) :: znum_list
    integer :: i, j, isnew
    INTEGER temp1
    REAL temp2 ! for temporary storage
    
    m%nelements=1
    znum_list(1) = m%znum(1)
    
    do i=1,m%natoms
       isnew = 1
       do j=1,m%nelements
          if(m%znum(i) == znum_list(j)) isnew = 0
       enddo
       if (isnew /= 0) then
          m%nelements=m%nelements+1
          znum_list(m%nelements) = m%znum(i)
       endif
    enddo

    allocate(m%atom_type(m%nelements), m%composition(m%nelements))
    m%atom_type = znum_list(1:m%nelements)
    m%composition = 0.0

    do i = 1, m%natoms
       do j=1,m%nelements
          if(m%atom_type(j) == m%znum(i)) then
             m%composition(j) = m%composition(j) + 1.0
             cycle
          endif
       enddo
    enddo

    m%composition = m%composition / real(m%natoms)

    !Put the elements in atom_type and composition in atomic
    !number increasing order. Added by Feng Yi on 01/21/2009
    !Floating bubble sort method
    DO i=1, (SIZE(m%atom_type)-1)
       DO j=1, (SIZE(m%atom_type)-i)
          IF(m%atom_type(j) .GT. m%atom_type(j+1)) THEN
             temp1 = m%atom_type(j)
             m%atom_type(j) = m%atom_type(j+1) 
             m%atom_type(j+1) = temp1
             
             temp2 = m%composition(j)
             m%composition(j) = m%composition(j+1)
             m%composition(j+1) = temp2
          ENDIF
       ENDDO
    ENDDO
    
  end subroutine composition_model


  ! Shifts the atom positions in model so that the mid-point between the maximum and 
  ! and minimum atom positions in each dimensions sits at the position (xc, yc, zc),
  ! measured in units of the model supercell.
  subroutine recenter_model(xc, yc, zc, m)
    real, intent(in) :: xc, yc, zc
    type(model), intent(inout) :: m

    real :: xshift, yshift, zshift

    xshift = xc*m%lx - (maxval(m%xx) + minval(m%xx))/2.0
    yshift = yc*m%ly - (maxval(m%yy) + minval(m%yy))/2.0
    zshift = zc*m%lz - (maxval(m%zz) + minval(m%zz))/2.0

    m%xx = m%xx+xshift
    m%yy = m%yy+yshift
    m%zz = m%zz+zshift

  end subroutine recenter_model


  ! simple error checking on a model.  Currently checks: are all the 
  ! atoms in the box?  Are all the atomic numbers between 1 and 103
  ! (the range for which Kirkland calculated electron scattering factors)
  ! More should be added as we think of it.
  subroutine check_model(m, istat)
    type(model), intent(in) :: m
    integer, intent(out) :: istat

    real xlen, ylen, zlen

    istat = 0

    xlen = maxval(m%xx) - minval(m%xx)
    ylen = maxval(m%yy) - minval(m%yy)
    zlen = maxval(m%zz) - minval(m%zz)

    if ( xlen > m%lx ) then
       write (*,*) 'Maximum x distance of ',xlen,' Ang exceeds box size ',m%lx,' Ang.'
       istat = 1
    end if

    if ( ylen > m%ly ) then
       write (*,*) 'Maximum y distance of ',ylen,' Ang exceeds box size ',m%ly,' Ang.'
       istat = 1
    end if

    if ( zlen > m%lz ) then
       write (*,*) 'Maximum z distance of ',zlen,' Ang exceeds box size ',m%lz,' Ang.'
       istat = 1
    end if

    if (minval(m%znum) < 1) then
       write (*,*) 'Minimum atomic number of ', minval(m%znum, 1), 'is less than zero.'
       istat = 1
    end if

    if (maxval(m%znum) > 103) then
       write (*,*) 'Maximum atomic number of ', maxval(m%znum, 1), 'is greater than 103.'
       istat = 1
    end if

  end subroutine check_model


  ! Makes (xp, yp, zp) copies of the input model min and puts them in the output model 
  ! mout.  Returns non-zero in istat is the memory can't be allocated.
  subroutine periodic_continue_model(xp, yp, zp, min, mout, init_hutch, istat)
    integer, intent(in):: xp, yp, zp
    type(model), intent(in) :: min
    type(model), intent(out) :: mout
    logical, intent(in) :: init_hutch
    integer, intent(out) :: istat

    integer :: i, j, k, c
    real :: shift_x, shift_y, shift_z

    mout%natoms = min%natoms*xp*yp*zp

    allocate(mout%xx(mout%natoms), mout%yy(mout%natoms), mout%zz(mout%natoms), &
         mout%znum(mout%natoms), stat=istat)
    if(istat /= 0) then
       write (*,*) 'Error allocating memory for the periodic continued model.'
       return
    end if

    mout%lx = min%lx*real(xp)
    mout%ly = min%ly*real(yp)
    mout%lz = min%lz*real(zp)

    c=0
    do i = 0, xp-1
       shift_x = real(i)*min%lx
       do j = 0, yp-1
          shift_y = real(j)*min%ly
          do k = 0, zp-1
             shift_z = real(k)*min%lz
             mout%xx(c*min%natoms+1:(c+1)*min%natoms) = min%xx + shift_x
             mout%yy(c*min%natoms+1:(c+1)*min%natoms) = min%yy + shift_y
             mout%zz(c*min%natoms+1:(c+1)*min%natoms) = min%zz + shift_z
             mout%znum(c*min%natoms+1:(c+1)*min%natoms) = min%znum
             c = c+1
          end do
       end do
    end do

    call recenter_model(0., 0., 0., mout)

    ! continuation doesn't change the composition, so just copy it
    mout%nelements = min%nelements
    allocate(mout%atom_type(mout%nelements), mout%composition(mout%nelements), stat=istat)
    if(istat /= 0) then
       write (*,*) 'Problem allocating memory in periodic_continue_model.'
       return
    endif
    mout%atom_type = min%atom_type
    mout%composition = mout%composition

    if(init_hutch) then
       call model_init_hutch(mout, istat)
       if(istat /= 0) then
          write (*,*) 'Cannot allocate memeory for the new hutch_array.'
          return
       endif
    endif

  end subroutine periodic_continue_model


  ! rotates model min by angles theta, phi, psi, and puts the results in mrot
  subroutine rotate_model(theta, phi, psi, min, mrot, istat)
    real, intent(in) :: theta, phi, psi
    type(model), intent(in) :: min
    type(model), intent(out) :: mrot
    integer, intent(out) :: istat

    real, dimension(3,3) :: r                         ! rotation matrix
    real :: cpsi, cphi, ctheta, sphi, spsi, stheta    ! sines and cosines of the angles
    integer :: i, j                                   ! loop counters
    real :: x, y, z                                   ! temporary positions
    real :: lx2, ly2, lz2                             ! half box sizes
    type(model) :: mt                                 ! temporary oversize model
    integer, dimension(:), allocatable :: orig_indices

    ! periodic continue mt to 3x3x3 of the original model
    call copy_model(min, mt, istat)
    if (istat /= 0) return
    call periodic_continue_model(3, 3, 3, min, mt, .FALSE., istat)
    if (istat /= 0) return

    allocate(orig_indices(mt%natoms), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Memory allocation failure in rotate_model.'
       return
    endif
    orig_indices = (/ (mod(i,min%natoms)+1, i=1,mt%natoms) /)

    ! generate the members of a 3x3 rotation matrix.  Use the Goldstein "x-convention" 
    ! and Euler angles phi theta, psi.
    cpsi = cos(psi)
    cphi = cos(phi)
    ctheta = cos(theta)
    sphi = sin(phi)
    spsi = sin(psi)
    stheta = sin(theta)

    r(1,1) = cpsi*cphi - ctheta*sphi*spsi 
    r(1,2) = cpsi*sphi + ctheta*cphi*spsi 
    r(1,3) = spsi*stheta 
    r(2,1) = -1.0*spsi*cphi - ctheta*sphi*cpsi 
    r(2,2) = -1.0*spsi*sphi + ctheta*cphi*cpsi 
    r(2,3) = cpsi*stheta 
    r(3,1) = stheta*sphi 
    r(3,2) = -1.0*stheta*cphi 
    r(3,3) = ctheta 

    ! Rotate the position vectors in mt
    do i=1,mt%natoms
       x = mt%xx(i)*r(1,1) + mt%yy(i)*r(1,2) + mt%zz(i)*r(1,3) 
       y = mt%xx(i)*r(2,1) + mt%yy(i)*r(2,2) + mt%zz(i)*r(2,3) 
       z = mt%xx(i)*r(3,1) + mt%yy(i)*r(3,2) + mt%zz(i)*r(3,3) 
       mt%xx(i) = x
       mt%yy(i) = y
       mt%zz(i) = z
    end do


    ! cut the temporary model back to the original box size
    ! first count the atoms in the box
    mrot%natoms = 0
    lx2 = min%lx / 2.0
    ly2 = min%ly / 2.0
    lz2 = min%lz / 2.0
    do i=1, mt%natoms
       if (mt%xx(i) <= lx2 .AND. mt%xx(i) >= -1.0*lx2) then
          if (mt%yy(i) <= ly2 .AND. mt%yy(i) >= -1.0*ly2) then
             if (mt%zz(i) <= lz2 .AND. mt%zz(i) >= -1.0*lz2) then
                mrot%natoms = mrot%natoms + 1
             endif
          endif
       endif
    enddo

    if (abs( real(mrot%natoms - min%natoms) / real(min%natoms) ) >= 0.02 & 
         .AND. (mrot%natoms - min%natoms) > 2 ) then
       write (*,*) 'Potential problem: # of atoms in the rotated model differs from the original'
       write (*,*) 'by more than 2%.  Original number is ',min%natoms,'.  New number is',mrot%natoms,'.'
    endif

    ! allocate memory for the new atoms
    mrot%unrot_natoms = min%natoms
    allocate(mrot%xx(mrot%natoms), mrot%yy(mrot%natoms), mrot%zz(mrot%natoms), &
         mrot%znum(mrot%natoms),  mrot%rot_i(mrot%unrot_natoms), stat=istat)
    if( istat /= 0) then
       write (*,*) 'Problem allocating memory in rotate_model.'
       return
    endif

    do i=1,mrot%unrot_natoms
       mrot%rot_i(i)%nat = 0
       nullify(mrot%rot_i(i)%ind)
    enddo

    ! now copy just the atoms inside the box from the temp model to the rotated one
    j=1
    do i=1, mt%natoms
       if (mt%xx(i) <= lx2 .AND. mt%xx(i) >= -1.0*lx2) then
          if (mt%yy(i) <= ly2 .AND. mt%yy(i) >= -1.0*ly2) then
             if (mt%zz(i) <= lz2 .AND. mt%zz(i) >= -1.0*lz2) then
                mrot%xx(j) = mt%xx(i)
                mrot%yy(j) = mt%yy(i)
                mrot%zz(j) = mt%zz(i)
                mrot%znum(j) = mt%znum(i)
                call add_index(mrot%rot_i(orig_indices(i)), j)
                j = j+1
             endif
          endif
       endif
    enddo

    ! set the rest of of the rotated model paramters
    mrot%lx = min%lx
    mrot%ly = min%ly
    mrot%lz = min%lz
    mrot%rotated = .TRUE.
    mrot%unrot_natoms = min%natoms
    call composition_model(mrot) ! have to recalculate this because the # of atoms
                                 ! may have changed a little

    call model_init_hutch(mrot, istat)

  end subroutine rotate_model

  ! copies the contents of model min to model mout
  subroutine copy_model(min, mout, istat)
    type(model), intent(in) :: min
    type(model), intent(out) :: mout
    integer, intent(out) :: istat

    integer i

    mout%lx = min%lx
    mout%ly = min%ly
    mout%lz = min%lz
    mout%natoms = min%natoms
    mout%rotated = min%rotated
    mout%nelements = min%nelements

    allocate(mout%xx(mout%natoms), mout%yy(mout%natoms), mout%zz(mout%natoms), &
         mout%znum(mout%natoms), mout%atom_type(mout%nelements), mout%composition(mout%nelements), &
         stat=istat)

    if (istat /= 0) then
       write (*,*) 'Problem allocating memory in copy_model.'
       return
    end if

    mout%xx = min%xx
    mout%yy = min%yy
    mout%zz = min%zz
    mout%znum = min%znum
    mout%atom_type = min%atom_type
    mout%composition = min%composition
    mout%rotated = min%rotated

    if(min%rotated) then
       allocate(mout%rot_i(min%unrot_natoms), stat=istat)
       if(istat /= 0) then
          write (*,*) 'Problem allocating memory in copy_model.'
          return
       endif
       do i=1,min%unrot_natoms
          allocate(mout%rot_i(i)%ind(min%rot_i(i)%nat))
          mout%rot_i(i)%nat = min%rot_i(i)%nat
          mout%rot_i(i)%ind = min%rot_i(i)%ind
       enddo
    endif
    
    call model_init_hutch(mout, istat)

  end subroutine copy_model


  ! Deallocates all the various allocatable arrays and sub-arrays in a model.
  subroutine destroy_model(m)
    type(model), intent(inout) :: m

    deallocate(m%xx, m%yy, m%zz, m%znum, m%atom_type, m%composition)
    call destroy_hutch(m%ha)
    call destroy_rot_indices(m%unrot_natoms, m%rot_i)
  end subroutine destroy_model

  ! add index i to the index list inside il
  subroutine add_index(il, i)
    type(index_list), intent(inout) :: il
    integer, intent(in) :: i

    integer, dimension(:), allocatable :: scratch

    if( il%nat >= 1 ) then
       allocate(scratch(il%nat))
       scratch = il%ind
       il%nat = il%nat+1
       deallocate(il%ind)
       allocate(il%ind(il%nat))
       il%ind(1:il%nat-1) = scratch
       il%ind(il%nat) = i
    else
       il%nat = 1
       allocate(il%ind(1))
       il%ind(1) = i
    endif

  end subroutine add_index

  ! deallocates all of the allocatable arrays and sub-arrays in an index list.  used by destroy_model
  subroutine destroy_rot_indices(unrot_natoms, ri)
    integer, intent(in) :: unrot_natoms
    !type(index_list), pointer, dimension(:), intent(inout) :: ri
   type(index_list), pointer, dimension(:) :: ri
    integer i
    
    do i=1, unrot_natoms
       deallocate(ri(i)%ind)
    enddo

    deallocate(ri)
  end subroutine destroy_rot_indices


  ! initializes the hutch_array ha within the model.  It calcualtes the hutch_size based on the
  ! model box size and the parameter ATOMS_PER_HUTCH, the number of hutches in 
  ! the array nhutch_x, nhutch_y, and hhutch_z, then assigns all the atoms in
  ! the current model atom position arrays xa, ya, and za to the appropriate
  ! hutches.  It does NOT check whether ha has already been initialized, so this
  ! routine should NEVER be called more than once for the same hutch_array.
  subroutine model_init_hutch(m, status)
    type(model), intent(inout) :: m
    integer, intent(out) :: status

    integer :: istat, i, hx, hy, hz


    write (*,*) 'Initializing hutch array'

    status = 0

    m%ha%hutch_size = ((m%lx*m%ly*m%lz)/m%natoms*ATOMS_PER_HUTCH)**(1./3.)
    write (*,*) 'Hutch size is ',m%ha%hutch_size,' Angstroms.'

    m%ha%nhutch_x = ceiling( m%lx / m%ha%hutch_size)
    m%ha%nhutch_y = ceiling( m%ly / m%ha%hutch_size)
    m%ha%nhutch_z = ceiling( m%lz / m%ha%hutch_size)
	!write (*,*)m%ha%nhutch_x,m%ha%nhutch_y,m%ha%nhutch_z

    allocate(m%ha%h(m%ha%nhutch_x, m%ha%nhutch_y, m%ha%nhutch_z), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Cannot allocate memory for hutch algorithm.  Exiting.'
       status = 1
       return
    end if

    allocate(m%ha%atom_hutch(m%natoms, 3), stat=istat)
    if (istat /= 0) then
       write(*,*) 'Cannot allocate memory for hutch algorithm.  Exiting.'
       status = 1
       return
    end if
    m%ha%atom_hutch = 0

    do hx = 1, m%ha%nhutch_x
       do hy = 1, m%ha%nhutch_y
          do hz = 1, m%ha%nhutch_z
             nullify(m%ha%h(hx, hy, hz)%at)
             m%ha%h(hx, hy, hz)%nat = 0
          end do
       end do
    end do

    do i=1, m%natoms
       call hutch_position(m, m%xx(i), m%yy(i), m%zz(i), hx, hy, hz)
       if(hz < 1) then
          write (*,*) 'Atom i=',i,' position is ',m%xx(i),', ',m%yy(i),', ',m%zz(i)
       endif
       call hutch_add_atom(m, i, hx, hy, hz)
    end do

  end subroutine model_init_hutch
  

  ! Moves the atom with index atom from its current hutch in hutch_array to
  ! to the hutch that encompasses position (xx, yy, zz).  Used to update the
  ! hutch_array for a Monte Carlo move of one atom.
  subroutine hutch_move_atom(m, atom, xx, yy, zz)
    type(model), target, intent(inout) :: m
    integer, intent(in) :: atom
    real, intent(in) :: xx, yy, zz

    integer :: hx, hy, hz
    type(hutch_array), pointer :: ha
    ha => m%ha
    
    !write (*,*) 'Inside move_atom.'

    call hutch_remove_atom(m, atom)

    !write(*,*) 'hutch_size = ',hutch_size
    call hutch_position(m, xx, yy, zz, hx, hy, hz)

    !write (*,*) 'just before add_atom, hx,hy,hz = ',hz,hy,hz
    call hutch_add_atom(m, atom, hx, hy, hz)
    ha%atom_hutch(atom, 1) = hx
    ha%atom_hutch(atom, 2) = hy
    ha%atom_hutch(atom, 3) = hz

  end subroutine hutch_move_atom


  ! Makes list (in atoms) of the indices of atoms that are in a cube of side 
  ! length radius, centered on the hutch containing  the point (px, py, pz),
  ! using the hutch array.  Useful for calculating G(r).  Returns 1 in istat
  ! if memory allocation fails and -1 if no atoms are found.
  subroutine hutch_list_3D(m, px, py, pz, radius, atoms, istat, nlist)

	type(model), target, intent(in) :: m
    real, intent(in) :: px, py, pz
    real, intent(in) :: radius
   !integer,  dimension(:), pointer, intent(out) :: atoms
	integer,  dimension(:), pointer :: atoms
    integer, intent(out) :: nlist        ! number of atoms in list
   !integer :: nlist        ! number of atoms in list
    integer, intent(out) :: istat
    integer :: hx, hy, hz   ! hutch of position (px, py, pz)
    integer :: nh           ! number of hutches corresponding to radius
    integer :: i, j, k      ! counting variables
    integer :: hi, hj, hk   ! counting variables with periodic boundary conditions
    integer, allocatable, dimension(:), target :: temp_atoms  ! temporary atoms list
    type(hutch_array), pointer :: ha

	
!	write(*,*)"testline1",px,py,pz
    ha => m%ha
    allocate(temp_atoms(m%natoms), stat=istat)
!	allocate(temp_atoms(m%natoms*6), stat=istat) !for debug

    if (istat /= 0) then
       write (*,*) 'Cannot allocate index list in hutch_list_3D.'
       return
    endif

    call hutch_position(m, px, py, pz, hx, hy, hz)
!	write(*,*)"testline2",radius
    nh = ceiling(radius / ha%hutch_size)

    ! accumulate the atoms we want to keep into temp_atoms
    nlist = 1
    do k = (hz-nh), (hz+nh)
      if (k > ha%nhutch_z) then
         hk = k - ha%nhutch_z
       else if (k < 1) then
          hk = k + ha%nhutch_z
       else
          hk = k
       end if
       do j = (hy-nh), (hy+nh)
          if (j > ha%nhutch_y) then
             hj = j - ha%nhutch_y
          else if (j < 1) then
             hj = j+ ha%nhutch_y
          else
             hj = j
          end if
          do i = (hx-nh), (hx+nh)
             if (i > ha%nhutch_x) then !Periodic boundary condition
                hi = i - ha%nhutch_x
            else if (i < 1) then
                hi = i + ha%nhutch_x
             else
                hi = i
             end if
             if (ha%h(hi, hj, hk)%nat > 0) then
			    !PRINT *, 'nlist+ha%h(hi,hj,hk)%nat-1 is: ', nlist+ha%h(hi,hj,hk)%nat-1 !for debug
                temp_atoms(nlist:nlist+ha%h(hi,hj,hk)%nat-1) = ha%h(hi,hj,hk)%at
                nlist = nlist+ha%h(hi,hj,hk)%nat
             endif
          end do
       end do
    end do
!	write(*,*)"testline3",nlist
    allocate(atoms(nlist-1), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
       return
    endif

	!assign atoms to the subset of temp_atoms that was filled in
    if (nlist > 1) then
       atoms = temp_atoms(1:nlist-1)
    else
       nullify(atoms)
       istat = -1
    endif

  end subroutine hutch_list_3D


  ! Makes a list of atom indices (in atoms) of the atoms in a rectangular
  ! prism with side length diameter in x and y, through the model thickness
  ! in z, centered on the hutch containing the point (px, py).  Useful for
  ! calculating the FEM intensity at (px, py).  Returns 1 in istat if the
  ! atoms array cannot be allocated.
  subroutine hutch_list_pixel(m, px, py, diameter, atoms, istat)
    type(model), target, intent(in) :: m
    real, intent(in) :: px, py, diameter
    integer, pointer, dimension(:) :: atoms
    integer, intent(out) :: istat

    integer :: hx, hy, hz   ! hutch of position (px, py, pz)
    integer :: nh           ! number of hutches corresponding to diameter
    integer :: nlist        ! number of atoms in list
    integer :: i, j         ! counting variables
    integer :: hi, hj, hk   ! counting variables with periodic boundary conditions

    integer, dimension(:), allocatable, target :: temp_atoms
    type(hutch_array), pointer :: ha

    ha => m%ha

    allocate(temp_atoms(m%natoms), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel'
       return
    end if

    call hutch_position(m, px, py, 0.0, hx, hy, hz)
    nh = ceiling( (diameter/2.) / ha%hutch_size)

    !write (*,*) 'hutch size is :', ha%hutch_size
    !write (*,*) 'center top hutch is ', hx, hy, hz
    !write (*,*) 'hutch radius is: ', nh

    nlist = 1
    do hk = 1, ha%nhutch_z
       do j = (hy-nh), (hy+nh)
          if (j > ha%nhutch_y) then
             hj = j - ha%nhutch_y
          else if (j < 1) then
             hj = j + ha%nhutch_y
          else
             hj = j
          end if
          do i = (hx-nh), (hx+nh)
             if (i > ha%nhutch_x) then
                hi = i - ha%nhutch_x
             else if (i < 1) then
                hi = i + ha%nhutch_x
             else
                hi = i
             end if
              if (ha%h(hi, hj, hk)%nat > 0) then
                temp_atoms(nlist:nlist+ha%h(hi,hj,hk)%nat-1) = ha%h(hi,hj,hk)%at
                nlist = nlist+ha%h(hi,hj,hk)%nat
              endif
          end do
       end do
    end do
    
    allocate(atoms(nlist-1), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
       return
    endif

    ! assign atoms to the subset of temp_atoms that was filled in
    if (nlist > 1) then
       atoms = temp_atoms(1:nlist-1)

    else
       nullify(atoms)
       istat = -1
    endif

  end subroutine hutch_list_pixel


  ! returns the indices of the hutch that encompasses position (xx, yy, zz) in
  ! the hutch_array in the integers (hx, hy, hz).  It assumes that the model 
  ! extends from -lx/2 to lx/2, -ly/2 to ly/2 and -lz/2 to lz/2 and does no error
  ! checking.
  subroutine hutch_position(m, xx, yy, zz, hx, hy, hz)
    type(model), intent(in) :: m
    real, intent(in) :: xx, yy, zz
    integer, intent(out) :: hx, hy, hz

    hx = ceiling( (xx + 0.5*m%lx) / m%ha%hutch_size)
    hy = ceiling( (yy + 0.5*m%ly) / m%ha%hutch_size)
    hz = ceiling( (zz + 0.5*m%lz) / m%ha%hutch_size)

    if (hx == 0) hx = 1
    if (hy == 0) hy = 1
    if (hz == 0) hz = 1

  end subroutine hutch_position


  ! Adds the atom with index atom to the hutch_array in hutch hx, hy, hz.
  subroutine hutch_add_atom(m, atom, hx, hy, hz)
    type(model), target, intent(inout) :: m
    integer, intent(in) :: atom, hx, hy, hz

    integer :: nat
    integer, dimension(m%ha%h(hx, hy, hz)%nat+1) :: scratch_atoms
    type(hutch_array), pointer :: ha

    ha => m%ha

    !write (*,*) 'Inside add_atom atom = ',atom,' and (hx,hy,hz) = ',hx,',',hy,',',hz

    nat = ha%h(hx,hy,hz)%nat

    if(nat > 0) then

       scratch_atoms(1:nat) = ha%h(hx, hy, hz)%at
       scratch_atoms(nat+1) = atom

       deallocate(ha%h(hx,hy,hz)%at)
       allocate(ha%h(hx,hy,hz)%at(1:nat+1)) ! +1 for fencepost 

       ha%h(hx,hy,hz)%at = scratch_atoms

    else
       allocate(ha%h(hx,hy,hz)%at(1:1))
       ha%h(hx,hy,hz)%at(1) = atom
    end if

    ha%h(hx,hy,hz)%nat = nat+1
    ha%atom_hutch(atom, 1) = hx
    ha%atom_hutch(atom, 2) = hy
    ha%atom_hutch(atom, 3) = hz

  end subroutine hutch_add_atom


  ! remove atom atom from its current hutch.  This reduces the number of atoms in hutch
  ! array by one, so it should only be used in conjunction with hutch_add_atom.
  subroutine hutch_remove_atom(m, atom)
    type(model), target, intent(inout) :: m
    integer, intent(in) :: atom

    integer, dimension(m%ha%h(m%ha%atom_hutch(atom,1),m%ha%atom_hutch(atom,2),m%ha%atom_hutch(atom,3))%nat) :: scratch_atoms
    integer :: hx, hy, hz, i, j
    type(hutch_array), pointer :: ha

    ha => m%ha

    hx = ha%atom_hutch(atom,1)
    hy = ha%atom_hutch(atom,2)
    hz = ha%atom_hutch(atom,3)

    scratch_atoms = ha%h(hx,hy,hz)%at

    deallocate(ha%h(hx,hy,hz)%at)
    allocate(ha%h(hx,hy,hz)%at(1:ha%h(hx,hy,hz)%nat-1))

    j=1
    do i=1, ha%h(hx,hy,hz)%nat
       if (scratch_atoms(i) /= atom) then
          ha%h(hx,hy,hz)%at(j) = scratch_atoms(i)
          j=j+1
       end if
    end do

    ha%h(hx,hy,hz)%nat = ha%h(hx,hy,hz)%nat-1
    ha%atom_hutch(atom,1) = 0
    ha%atom_hutch(atom,2) = 0
    ha%atom_hutch(atom,3) = 0


  end subroutine hutch_remove_atom


  ! deallocates the hutch_array ha and the atom lists inside it.  used as part of destroy_model.
  subroutine destroy_hutch(ha)
    type(hutch_array), intent(inout) :: ha
    
    integer i, j, k
    do i=1, ha%nhutch_x
       do j=1, ha%nhutch_y
          do k=1, ha%nhutch_z
             deallocate(ha%h(i,j,k)%at)
          enddo
       enddo
    enddo

    deallocate(ha%h, ha%atom_hutch)

  end subroutine destroy_hutch


  ! diagnostic routine dump a hutch data structure to the stand output
  subroutine dump_hutches(m)
    type(model), target, intent(in) :: m

    type(hutch_array), pointer :: ha
    integer :: i, j, k

    ha => m%ha

    do i=1, ha%nhutch_x
       do j=1, ha%nhutch_y
          do k=1, ha%nhutch_z
             if(ha%h(i,j,k)%nat > 0) then
                write (*, "('Hutch ',i3,', ',i3,', ',i3,'  contains atoms: ',10i4)"), i, j, k, ha%h(i,j,k)%at
             else
                write (*,"( 'Hutch ',i3,', ',i3,', ',i3,'  contains no atoms.')"), i, j, k
             end if
          end do
       end do
    end do

    do i=1, m%natoms
       write (*,"('atom ',i5,' is in hutch ', 3i4)"), i, ha%atom_hutch(i,1), ha%atom_hutch(i,2), ha%atom_hutch(i,3)
    end do

  end subroutine dump_hutches



  !pre-calculate the relative hutch position with respect to the center hutch
  !consider the worst condition for 2D case
  SUBROUTINE pre_calc_2D_hutch(list1,ratio_list,size_ref)
    TYPE(hutch_2D_array) list1
    REAL, DIMENSION(:) :: ratio_list
    !ratio_list from 0 to 0.5, 0.5 to 1.5, then step size with 1
    !ratio_list(1)=0.5, ratio_list(2)=1.5, ratio_list(3)=2.5, etc. 
    REAL size_ref !square size
    
    
    INTEGER num_x, num_y
    INTEGER i1,j1,k1 
    INTEGER n1,n2
    INTEGER nh
    INTEGER num_temp1
    !INTEGER,ALLOCATABLE, DIMENSION(:, :) :: num_temp2
    INTEGER num_temp2
    INTEGER,ALLOCATABLE, DIMENSION(:, :) :: temp_list
    REAL r_dist, ref_dist
    REAL p_x, p_y
    
    REAL x_square, y_square
    
    
    num_x=2
    num_y=2
    
    list1%size_ratio = SIZE(ratio_list,1)
    list1%size_position =  num_x*num_y
    
    !POINTER allocation
    ALLOCATE(list1%list_2D(list1%size_position ,list1%size_ratio ))
    
    ALLOCATE(list1%ratio_radius_square(list1%size_ratio))
    list1%ratio_radius_square = ratio_list
    
    !IF (.NOT. ALLOCATED(num_temp2)) THEN
    ! ALLOCATE(num_temp2(list1%size_position ,list1%size_ratio))
    !ENDIF
    
    IF(.NOT. ALLOCATED(temp_list)) THEN
       ALLOCATE(temp_list(2, 30*30))  !this size can be modified
       !first dimension refers to x and y, 1 is x and 2 is y
       !second dimension refers to x or y relative position
    ENDIF
    
    
    
    
    DO i1=1, num_y
       DO j1=1, num_x
          num_temp1 = (i1-1) * num_x + j1
          
          
          DO k1=1, list1%size_ratio
             num_temp2 = 0
             !Calculate which hutch is within distance r to the point
             !consider the worst case
             ref_dist = size_ref * ratio_list(k1)
             nh = CEILING(ratio_list(k1))
             DO n1=-nh, nh
                !***********x dimension***********
                IF(n1 .GT. 0) THEN
                   !pick rightmost corner of square
                   p_x=(n1-1.0)*0.5
                ELSEIF (n1 .LT. 0) THEN
                   !pick leftmost corner of square
                   p_x=(n1-2.0)*0.5
                ELSE
                   p_x=0
                ENDIF
                DO n2=-nh, nh
                   !***********y dimension***********
                   IF(n2 .GT. 0) THEN
                      !pick rightmost corner of square
                      p_y=(n2-1.0)*0.5
                   ELSEIF (n2 .LT. 0) THEN
                      !pick leftmost corner of square
                      p_y=(n2-2.0)*0.5
                   ELSE
                      p_y=0
                   ENDIF
                   
                   
                   !closest point in a hutch to the origin
                   !pick the reference hutch center as the origin
                   IF(n1 .NE. 0) THEN
                      x_square = (n1-n1/ABS(n1)*0.5)*size_ref
                   ELSE
                      x_square = 0.0
                   ENDIF
                   
                   IF(n2 .NE. 0) THEN
                      y_square = (n2-n2/ABS(n2)*0.5)*size_ref
                   ELSE
                      y_square = 0.0
                   ENDIF
                   
                   !Calculate the closest distance between one part of the reference square
                   !and a square with relative number position n1 and n2
                   r_dist = sqrt((x_square-p_x)**2+(y_square-p_y)**2)
                   
                   !compare r_dist with ref_dist
                   IF(r_dist .LE. ref_dist) THEN
                      num_temp2 = num_temp2 + 1
                      temp_list(1,num_temp2) = n1
                      temp_list(2,num_temp2) = n2
                   ENDIF
                   
                   
                ENDDO !n2
             ENDDO !n1
             
             
             list1%list_2D(num_temp1,k1)%size_d = num_temp2
             ALLOCATE(list1%list_2D(num_temp1,k1)%list_x(num_temp2))
             ALLOCATE(list1%list_2D(num_temp1,k1)%list_y(num_temp2))
             
             list1%list_2D(num_temp1,k1)%list_x = temp_list(1, 1:num_temp2)
             list1%list_2D(num_temp1,k1)%list_y = temp_list(2, 1:num_temp2)
             
             
             
          ENDDO ! k1
          
          
          
       ENDDO !j1
    ENDDO  !i1
    
    IF (ALLOCATED(temp_list)) THEN
       DEALLOCATE(temp_list)
    ENDIF
    
  END SUBROUTINE pre_calc_2D_hutch
  


  !destroy hutch_2D_array
  SUBROUTINE destroy_hutch_2D_array(list1)
    TYPE(hutch_2D_array) list1
    
    INTEGER i1, j1
    
    DO i1=1, list1%size_position
       DO j1=1,list1%size_ratio
          DEALLOCATE(list1%list_2D(i1,j1)%list_x)
          DEALLOCATE(list1%list_2D(i1,j1)%list_y)
       ENDDO
    ENDDO
    
    DEALLOCATE(list1%list_2D)
    DEALLOCATE(list1%ratio_radius_square)
  END SUBROUTINE destroy_hutch_2D_array
  


  !pre-calculate the relative hutch position with respect to the center hutch
  !consider the worst condition for 3D case

  SUBROUTINE pre_calc_3D_hutch(list1,ratio_list,size_ref)
    TYPE(hutch_3D_array) list1
    REAL, DIMENSION(:) :: ratio_list
    !ratio_list from 0 to 0.5, 0.5 to 1.5, then step size with 1
    !ratio_list(1)=0.5, ratio_list(2)=1.5, ratio_list(3)=2.5, etc. 
    REAL size_ref !square size
	   

    INTEGER num_x, num_y, num_z
    INTEGER i1,j1,k1, k2 
    INTEGER n1,n2, n3
    INTEGER nh
    INTEGER num_temp1
    !INTEGER,ALLOCATABLE, DIMENSION(:, :) :: num_temp2
    INTEGER num_temp2
    INTEGER,ALLOCATABLE, DIMENSION(:, :) :: temp_list
    REAL r_dist, ref_dist
    REAL p_x, p_y, p_z
    
    REAL x_hutch, y_hutch, z_hutch
    
    
    num_x=2
    num_y=2
    num_z=2
    
    list1%size_ratio = SIZE(ratio_list,1)
    list1%size_position =  num_x*num_y*num_z
    
    !POINTER allocation
    ALLOCATE(list1%list_3D(list1%size_position ,list1%size_ratio ))
    
    ALLOCATE(list1%ratio_radius_hutch(list1%size_ratio))
    list1%ratio_radius_hutch = ratio_list
    
    !IF (.NOT. ALLOCATED(num_temp2)) THEN
    ! ALLOCATE(num_temp2(list1%size_position ,list1%size_ratio))
    !ENDIF
    
    IF(.NOT. ALLOCATED(temp_list)) THEN
       ALLOCATE(temp_list(3, 30*30*30))  !this size can be modified
       !first dimension refers to x and y, 1 is x and 2 is y
       !second dimension refers to x or y relative position
    ENDIF
    
    
    DO i1=1, num_z
       DO j1=1, num_y
          DO k1=1,num_x

             num_temp1 = (i1-1) * num_y * num_x + (j1-1) * num_x + k1
             
             
             DO k2=1, list1%size_ratio
                num_temp2 = 0
                !Calculate which hutch is within distance r to the point
                !consider the worst case
                ref_dist = size_ref * ratio_list(k2)
                nh = CEILING(ratio_list(k2))
                DO n1=-nh, nh
                   !***********x dimension***********
                   IF(n1 .GT. 0) THEN
                      !pick rightmost corner of square
                      p_x=(n1-1.0)*0.5
                   ELSEIF (n1 .LT. 0) THEN
                      !pick leftmost corner of square
                      p_x=(n1-2.0)*0.5
                   ELSE
                      p_x=0
                   ENDIF
                   DO n2=-nh, nh
                      !***********y dimension***********
                      IF(n2 .GT. 0) THEN
                         !pick rightmost corner of square
                         p_y=(n2-1.0)*0.5
                      ELSEIF (n2 .LT. 0) THEN
                         !pick leftmost corner of square
                         p_y=(n2-2.0)*0.5
                      ELSE
                         p_y=0
                      ENDIF
                      
                      DO n3=-nh, nh
                         !***********y dimension***********
                         IF(n3 .GT. 0) THEN
                            !pick rightmost corner of square
                            p_z=(n3-1.0)*0.5
                         ELSEIF (n3 .LT. 0) THEN
                            !pick leftmost corner of square
                            p_z=(n3-2.0)*0.5
                         ELSE
                            p_z=0
                         ENDIF
                         
                         
                         !closest point in a hutch to the origin
                         !pick the reference hutch center as the origin
                         IF(n1 .NE. 0) THEN
                            x_hutch = (n1-n1/ABS(n1)*0.5)*size_ref
                         ELSE
                            x_hutch = 0.0
                         ENDIF
                         
                         IF(n2 .NE. 0) THEN
                            y_hutch = (n2-n2/ABS(n2)*0.5)*size_ref
                         ELSE
                            y_hutch = 0.0
                         ENDIF
                         
                         IF(n3 .NE. 0) THEN
                            z_hutch = (n3-n3/ABS(n3)*0.5)*size_ref
                         ELSE
                            z_hutch = 0.0
                         ENDIF
                         
                         !Calculate the closest distance between one part of the reference square
                         !and a square with relative number position n1 and n2
                         r_dist = sqrt((x_hutch-p_x)**2+(y_hutch-p_y)**2+(z_hutch-p_z)**2)
                         
                         !compare r_dist with ref_dist
                         IF(r_dist .LE. ref_dist) THEN
                            num_temp2 = num_temp2 + 1
                            temp_list(1,num_temp2) = n1
                            temp_list(2,num_temp2) = n2
                            temp_list(3,num_temp2) = n3
                         ENDIF
                         
                      ENDDO !n3
                      
                   ENDDO !n2
                   
                ENDDO !n1
                
                
                list1%list_3D(num_temp1,k2)%size_d = num_temp2
                ALLOCATE(list1%list_3D(num_temp1,k2)%list_x(num_temp2))
                ALLOCATE(list1%list_3D(num_temp1,k2)%list_y(num_temp2))
                
                list1%list_3D(num_temp1,k2)%list_x = temp_list(1, 1:num_temp2)
                list1%list_3D(num_temp1,k2)%list_y = temp_list(2, 1:num_temp2)
                list1%list_3D(num_temp1,k2)%list_y = temp_list(3, 1:num_temp2)
                
                
                
             ENDDO ! k2
             
             
          ENDDO !k1
       ENDDO !j1
    ENDDO  !i1
    
    IF (ALLOCATED(temp_list)) THEN
       DEALLOCATE(temp_list)
    ENDIF
    
  END SUBROUTINE pre_calc_3D_hutch
  
  

  !destroy hutch_3D_array
  SUBROUTINE destroy_hutch_3D_array(list1)
    TYPE(hutch_3D_array) list1
    
    INTEGER i1, j1
	  
    DO i1=1, list1%size_position
       DO j1=1,list1%size_ratio
          DEALLOCATE(list1%list_3D(i1,j1)%list_x)
          DEALLOCATE(list1%list_3D(i1,j1)%list_y)
          DEALLOCATE(list1%list_3D(i1,j1)%list_z)
       ENDDO
    ENDDO
    
    DEALLOCATE(list1%list_3D)
    DEALLOCATE(list1%ratio_radius_hutch)
  END SUBROUTINE destroy_hutch_3D_array
  
  

  ! Makes list (in atoms) of the indices of atoms that are in a space of 
  ! length radius, centered on the hutch containing  the point (px, py, pz),
  ! using the hutch array.  Useful for calculating G(r).  Returns 1 in istat
  ! if memory allocation fails and -1 if no atoms are found.
  ! more efficient than old algorithm
  subroutine hutch_list_3D_eff(m, px, py, pz, radius, atoms, istat, nlist,list1)

    type(model), target, intent(in) :: m
    TYPE(hutch_3D_array), INTENT(IN) :: list1 !store relative position of each hutch
	                                      !which are need for calculation
    real, intent(in) :: px, py, pz
    real, intent(in) :: radius
    !integer,  dimension(:), pointer, intent(out) :: atoms
    integer,  dimension(:), pointer :: atoms
    integer, intent(out) :: nlist        ! number of atoms in list
    !integer :: nlist        ! number of atoms in list
    integer, intent(out) :: istat
    integer :: hx, hy, hz   ! hutch of position (px, py, pz)
    integer :: nh           ! number of hutches corresponding to radius
    integer :: i, j, k, k1    ! counting variables
    integer :: hi, hj, hk   ! counting variables with periodic boundary conditions
    integer, allocatable, dimension(:), target :: temp_atoms  ! temporary atoms list
    type(hutch_array), pointer :: ha

    INTEGER p_relative_3D, p_relative_2D
    INTEGER ratio_position
    REAL ratio1
    LOGICAL New_Algorithm
!	write(*,*)"testline1",px,py,pz
    ha => m%ha
    allocate(temp_atoms(m%natoms), stat=istat)
    !	allocate(temp_atoms(m%natoms*6), stat=istat) !for debug
    
    if (istat /= 0) then
       write (*,*) 'Cannot allocate index list in hutch_list_3D.'
       return
    endif
    
    call hutch_position_eff(m, px, py, pz, hx, hy, hz,p_relative_3D, p_relative_2D)
!	write(*,*)"testline2",radius
    ratio1 = radius / ha%hutch_size
    nh = ceiling(ratio1)

    New_Algorithm = .FALSE.

    !first, find which ratio range 
    DO i=1, list1%size_ratio
       IF(i .EQ. 1) THEN
          IF((ratio1 .GE. 0) .AND. (ratio1 .LE. list1%ratio_radius_hutch(i))) THEN
             ratio_position = i
             New_Algorithm = .TRUE. 
             EXIT
          ENDIF
       ELSE
          IF((ratio1 .GE. list1%ratio_radius_hutch(i)) .AND. (ratio1 .LE. list1%ratio_radius_hutch(i+1))) THEN
             ratio_position = i
             New_Algorithm = .TRUE. 
             EXIT
          ENDIF
       ENDIF
    ENDDO
    
    
    IF (New_Algorithm) THEN
       nlist = 1
       DO k1=1, list1%list_3D(p_relative_3D, ratio_position)%size_d
	  k = list1%list_3D(p_relative_3D, ratio_position)%list_z(k1)
	  k = hz + k
	  if (k > ha%nhutch_z) then
             hk = k - ha%nhutch_z
          else if (k < 1) then
             hk = k + ha%nhutch_z
          else
             hk = k
          end if

	  j = list1%list_3D(p_relative_3D, ratio_position)%list_y(k1)
	  j = hy + j
	  if (j > ha%nhutch_y) then
             hj = j - ha%nhutch_y
          else if (j < 1) then
             hj = j+ ha%nhutch_y
          else
             hj = j
          end if
          
	  i = list1%list_3D(p_relative_3D, ratio_position)%list_x(k1)
	  i = hx + i
	  if (i > ha%nhutch_x) then !Periodic boundary condition
             hi = i - ha%nhutch_x
          else if (i < 1) then
             hi = i + ha%nhutch_x
          else
             hi = i
          end if
          
	  if (ha%h(hi, hj, hk)%nat > 0) then
             !PRINT *, 'nlist+ha%h(hi,hj,hk)%nat-1 is: ', nlist+ha%h(hi,hj,hk)%nat-1 !for debug
             temp_atoms(nlist:nlist+ha%h(hi,hj,hk)%nat-1) = ha%h(hi,hj,hk)%at
             nlist = nlist+ha%h(hi,hj,hk)%nat
          endif


       ENDDO !k1
    ELSE
   !Old algorithm
   ! accumulate the atoms we want to keep into temp_atoms
       nlist = 1
       do k = (hz-nh), (hz+nh)
          if (k > ha%nhutch_z) then
             hk = k - ha%nhutch_z
          else if (k < 1) then
             hk = k + ha%nhutch_z
          else
             hk = k
          end if
          do j = (hy-nh), (hy+nh)
             if (j > ha%nhutch_y) then
                hj = j - ha%nhutch_y
             else if (j < 1) then
                hj = j+ ha%nhutch_y
             else
                hj = j
             end if
             do i = (hx-nh), (hx+nh)
                if (i > ha%nhutch_x) then !Periodic boundary condition
                   hi = i - ha%nhutch_x
                else if (i < 1) then
                   hi = i + ha%nhutch_x
                else
                   hi = i
                end if
                if (ha%h(hi, hj, hk)%nat > 0) then
                   !PRINT *, 'nlist+ha%h(hi,hj,hk)%nat-1 is: ', nlist+ha%h(hi,hj,hk)%nat-1 !for debug
                   temp_atoms(nlist:nlist+ha%h(hi,hj,hk)%nat-1) = ha%h(hi,hj,hk)%at
                   nlist = nlist+ha%h(hi,hj,hk)%nat
                endif
             end do
          end do
       end do
       
    ENDIF !which algorithm to be picked up
    !	write(*,*)"testline3",nlist
    allocate(atoms(nlist-1), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
       return
    endif
    
    !assign atoms to the subset of temp_atoms that was filled in
    if (nlist > 1) then
       atoms = temp_atoms(1:nlist-1)
    else
       nullify(atoms)
       istat = -1
    endif
    
  end subroutine hutch_list_3D_eff


  ! returns the indices of the hutch that encompasses position (xx, yy, zz) in
  ! the hutch_array in the integers (hx, hy, hz).  It assumes that the model 
  ! extends from -lx/2 to lx/2, -ly/2 to ly/2 and -lz/2 to lz/2 and does no error
  ! checking.
  ! also return relative position of point (xx, yy, zz) in the hutch
  ! in 2D and 3D cases
  ! the hutch is divided into equal 8 parts 2*2*2
  subroutine hutch_position_eff(m, xx, yy, zz, hx, hy, hz,p_relative_3D, p_relative_2D)
    type(model), intent(in) :: m
    real, intent(in) :: xx, yy, zz
    integer, intent(out) :: hx, hy, hz
    INTEGER, INTENT(OUT) :: p_relative_2D, p_relative_3D

    REAL,DIMENSION(3) :: hutch_center !hutch center 
    INTEGER i
    INTEGER n1, n2, n3
    INTEGER num_x, num_y, num_z

    num_x=2
    num_y=2
    num_z=2
    
    hx = ceiling( (xx + 0.5*m%lx) / m%ha%hutch_size)
    hy = ceiling( (yy + 0.5*m%ly) / m%ha%hutch_size)
    hz = ceiling( (zz + 0.5*m%lz) / m%ha%hutch_size)
    
    if (hx == 0) hx = 1
    if (hy == 0) hy = 1
    if (hz == 0) hz = 1
    
    hutch_center(1) = (hx-0.5) * m%ha%hutch_size
    hutch_center(2) = (hy-0.5) * m%ha%hutch_size
    hutch_center(3) = (hz-0.5) * m%ha%hutch_size
    
    IF((xx + 0.5*m%lx) .GE. hutch_center(1)) THEN
       n3 = 2
    ELSE
       n3=1
    ENDIF
    
    IF((yy + 0.5*m%ly) .GE. hutch_center(2)) THEN
       n2 = 2
    ELSE
       n2=1
    ENDIF
    
    IF((zz + 0.5*m%lz) .GE. hutch_center(1)) THEN
       n1 = 2
    ELSE
       n1=1
    ENDIF
    
    p_relative_2D = (n2-1) * num_x + n3
    p_relative_3D = (n1-1) * num_y * num_x + (n2-1) * num_x + n3
    
    
  end subroutine hutch_position_eff
  

  ! Makes a list of atom indices (in atoms) of the atoms in a rectangular
  ! prism with side length diameter in x and y, through the model thickness
  ! in z, centered on the hutch containing the point (px, py).  Useful for
  ! calculating the FEM intensity at (px, py).  Returns 1 in istat if the
  ! atoms array cannot be allocated.
  subroutine hutch_list_pixel_eff(m, px, py, diameter, atoms, istat, list1)
    type(model), target, intent(in) :: m
    TYPE(hutch_2D_array), INTENT(IN) :: list1
    real, intent(in) :: px, py, diameter
    integer, pointer, dimension(:) :: atoms
    integer, intent(out) :: istat

    integer :: hx, hy, hz   ! hutch of position (px, py, pz)
    integer :: nh           ! number of hutches corresponding to diameter
    integer :: nlist        ! number of atoms in list
    integer :: i, j, k1        ! counting variables
    integer :: hi, hj, hk   ! counting variables with periodic boundary conditions

    integer, dimension(:), allocatable, target :: temp_atoms
    type(hutch_array), pointer :: ha

    INTEGER p_relative_3D, p_relative_2D
    INTEGER ratio_position
    REAL ratio1
    LOGICAL New_Algorithm
    
    ha => m%ha

    allocate(temp_atoms(m%natoms), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel'
       return
    end if
    
    call hutch_position_eff(m, px, py, 0.0, hx, hy, hz, p_relative_3D, p_relative_2D)
    
    ratio1 = (diameter/2.) / ha%hutch_size
    nh = ceiling( ratio1)
    
    !write (*,*) 'hutch size is :', ha%hutch_size
    !write (*,*) 'center top hutch is ', hx, hy, hz
    !write (*,*) 'hutch radius is: ', nh
    !first, find which ratio range 
    DO i=1, list1%size_ratio
       IF(i .EQ. 1) THEN
          IF((ratio1 .GE. 0) .AND. (ratio1 .LE. list1%ratio_radius_square(i))) THEN
             ratio_position = i
             New_Algorithm = .TRUE. 
             EXIT
          ENDIF
       ELSE
          IF((ratio1 .GE. list1%ratio_radius_square(i)) .AND. (ratio1 .LE. list1%ratio_radius_square(i+1))) THEN
             ratio_position = i
             New_Algorithm = .TRUE. 
             EXIT
          ENDIF
       ENDIF
    ENDDO
    
    IF(New_Algorithm) THEN
       nlist = 1
       DO hk = 1, ha%nhutch_z
	  DO k1=1, list1%list_2D(p_relative_2D, ratio_position)%size_d
             j = list1%list_2D(p_relative_2D, ratio_position)%list_y(k1)
             j = hy + j
             if (j > ha%nhutch_y) then
                hj = j - ha%nhutch_y
             else if (j < 1) then
                hj = j+ ha%nhutch_y
             else
                hj = j
             end if
             
             i = list1%list_2D(p_relative_2D, ratio_position)%list_x(k1)
             i = hx + i
             if (i > ha%nhutch_x) then !Periodic boundary condition
                hi = i - ha%nhutch_x
             else if (i < 1) then
                hi = i + ha%nhutch_x
             else
                hi = i
             end if
             
             if (ha%h(hi, hj, hk)%nat > 0) then
                temp_atoms(nlist:nlist+ha%h(hi,hj,hk)%nat-1) = ha%h(hi,hj,hk)%at
                nlist = nlist+ha%h(hi,hj,hk)%nat
             endif
             
	  ENDDO !end k1
       ENDDO !hk
    ELSE
       nlist = 1
       do hk = 1, ha%nhutch_z
          do j = (hy-nh), (hy+nh)
             if (j > ha%nhutch_y) then
                hj = j - ha%nhutch_y
             else if (j < 1) then
                hj = j + ha%nhutch_y
             else
                hj = j
             end if
             do i = (hx-nh), (hx+nh)
                if (i > ha%nhutch_x) then
                   hi = i - ha%nhutch_x
                else if (i < 1) then
                   hi = i + ha%nhutch_x
                else
                   hi = i
                end if
                if (ha%h(hi, hj, hk)%nat > 0) then
                   temp_atoms(nlist:nlist+ha%h(hi,hj,hk)%nat-1) = ha%h(hi,hj,hk)%at
                   nlist = nlist+ha%h(hi,hj,hk)%nat
                endif
             end do
          end do
       end do
       
    ENDIF !which algorithm to be picked up
    
    allocate(atoms(nlist-1), stat=istat)
    if (istat /= 0) then
       write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
       return
    endif
    
    ! assign atoms to the subset of temp_atoms that was filled in
    if (nlist > 1) then
       atoms = temp_atoms(1:nlist-1)
       
    else
       nullify(atoms)
       istat = -1
    endif
    
  end subroutine hutch_list_pixel_eff
  
end module model_mod

