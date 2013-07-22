!
!This module contains the main functions for the model.
!It was written by Jason Maldonis from the old model.f90 file on 06/28/13.
!

! Reads a model in the Kirkland .xyz file format from the file
! model_filename.
! Puts the first-line comment in "comment", the model in m, and
! returns 0 in
! istat if the file cant be opened or the memory allocation
! fails.

module model_mod
    use RMC_Global  ! Global variables
    implicit none
    ! derived data type for the hutches.  at is a pointer to a list
    ! of the indices
    ! of the atoms in this hutch.  nat is the number of atoms in
    ! this hutch
    type hutch
        integer, dimension(:), pointer :: at
        integer :: nat
    end type hutch

    ! derived data type for the hutch array: contains an array of hutches,
    ! plus supporting information
    type hutch_array
        ! array of hutch objects
        type(hutch), dimension(:,:,:), pointer :: h
        ! number of hutches in x, y, and z
        integer :: nhutch_x, nhutch_y, nhutch_z
        ! physical size of a hutch in Angstroms
        real :: hutch_size
        ! list of the hutch indices for every atom
        integer, pointer, dimension(:,:) :: atom_hutch
    end type hutch_array

    ! adjustable-size list of atom indices
    type index_list
        integer :: nat
        integer, pointer, dimension(:) :: ind
    end type index_list


    ! list for a rotated model that is the original models natoms long.  each atom in the list
    ! list is itself a list of the indices in the new, rotated model of the atoms corresponding
    ! to the the original atom in the unrotated model.
    ! for a rotated model, a list of the indices in the rotated model corresponding to each atom
    ! index in the original, unrotated model.  the original, unrotated model's index long.
    ! Defined type for a structural model with atoms positions and a bunch of metadata
    type model
        integer :: natoms                              ! number of atoms in the model
        real, pointer, dimension(:) :: xx, yy, zz      ! atom positions in Angstroms
        integer, pointer, dimension(:) :: znum, znum_r         ! atom atomic numbers, and reduced z numbners
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

    TYPE hutch_2D_list
    ! list the relative x, and y position of squre to the center square
    ! for the atom at a certain position in a squre
    ! with certain ratio of radius to the square size
        INTEGER(SELECTED_INT_KIND(3)), DIMENSION(:), POINTER :: list_x
        INTEGER(SELECTED_INT_KIND(3)), DIMENSION(:), POINTER :: list_y
        INTEGER size_d !number of elements in x, y dimension
    END TYPE hutch_2D_list

    TYPE hutch_2D_array
    ! list the relative x, and y position of square to the center square
    ! divide the square into equal 4 parts
    ! the ratio of radius to the square size range from 0 to 10.5, can be changed
    ! depending on calculation needed
        TYPE(hutch_2D_list), DIMENSION(:,:), POINTER :: list_2D
        !first index refers to relative position of a point in a square
        !second index refers to ratio of radius to square size
        !position=relative position of a point in a square
        !position=j*(number_x-1) + i, where j refers to which y position the
        !point is in
        !i refers to which x position the point is in
        INTEGER size_position
        REAL, DIMENSION(:), POINTER :: ratio_radius_square
        INTEGER size_ratio
    END TYPE hutch_2D_array

    TYPE hutch_3D_list
    ! list the relative x, y and Z position of squre to the center hutch
    ! for the atom at a certain position in a hutch
    ! with certain ratio of radius to the hutch size
        INTEGER(SELECTED_INT_KIND(3)), DIMENSION(:), POINTER :: list_x
        INTEGER(SELECTED_INT_KIND(3)), DIMENSION(:), POINTER :: list_y
        INTEGER(SELECTED_INT_KIND(3)), DIMENSION(:), POINTER :: list_z
        INTEGER size_d !number of elements in x, y and z dimension
    END TYPE hutch_3D_list

    TYPE hutch_3D_array
    ! list the relative x, and y position of hutch to the center square
    ! divide the square into equal 8 parts
    ! the ratio of radius to the hutch size range from 0 to 10.5, can be changed
    ! depending on calculation needed
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

    ! Global variables
    TYPE(hutch_2D_array),PRIVATE :: list_1_2D
    TYPE(hutch_3D_array),PRIVATE :: list_1_3D
    LOGICAL, PRIVATE, SAVE :: hlist_2D_calc = .FALSE.
    LOGICAL, PRIVATE, SAVE :: hlist_3D_calc = .FALSE.


contains

    subroutine read_model(model_filename, comment, m, istat)
    ! Reads a model in the Kirkland .xyz file format from the file model_filename.
    ! Puts the first-line comment in "comment", the model in m, and returns 0 in
    ! istat if the file cant be opened or the memory allocation fails.
    ! The subroutine composition_model has been incorporated into this
    ! subroutine by Jason on 06/28/13.
        implicit none
        character (len=*),intent(in) :: model_filename
        character (LEN=*),intent(out) :: comment
        type(model), intent(out) :: m
        integer, intent(out) :: istat      !0 for successful open, others for failure.
        integer :: i, j, atom_count=0, nat=0, atom_temp
        integer, dimension(103) :: elements=0
        real :: comp_temp

        ! Open file that contains the model information.
        open(1,file=model_filename,iostat=istat,status='old')
        if(istat.ne.0) then !Open fails
            write(*,*)"Error in opening flie,"," ",model_filename," status=",istat
            return
        endif

        read(1,*) ! Comment line.
        read(1,*) ! This line contains the box size (lx, ly, lz).
        ! Count how many atoms there are in the model.
        do while( atom_count .ne. -1)
            read(1,*) atom_count ! Read line.
            nat=nat+1.0
        enddo
        nat=nat-1.0

        rewind(1)

        ! Set the number of atoms in the model m and allocate space for each
        ! coordinate.
        m%natoms = nat
        allocate(m%xx(nat), m%yy(nat), m%zz(nat), m%znum(nat), stat=istat)
        ! Allocate should return 0 if successful.
        if(istat /= 0) then
            write (*,*) 'Unable to allocate memory for the model being read.'
            return
        endif

        ! Read in the first 80 characters of the comment
        read(1,'(a80)') comment
        ! Read in the box size.
        read(1,*) m%lx,m%ly,m%lz
        ! If the model is not a perfect cube then the rest of the calculations
        ! wont work, so we really should check that.
        if((m%lx /= m%ly) .or. (m%lx /= m%lz)) then
            write(*,*) "The model is not a cube and will work correctly. Exiting."
            return
        endif
        ! Read the atomic numbers and atom positions directly into the model.
        do i=1,nat
            read(1,*) m%znum(i),m%xx(i),m%yy(i),m%zz(i)
            ! If this atom has atomic number z, then increment the z position in
            ! the array elements. This counts the number of each atom type we have.
            elements(m%znum(i)) = elements(m%znum(i)) + 1
        enddo
        close(1)

        ! Count the number of elements we have in our model.
        m%nelements=0
        do i=1, 103
            if(elements(i) /= 0) then
                m%nelements = m%nelements + 1
            end if
        end do

        ! Note: nelements is usually between 1 and 5 so these loops are tiny.
        ! Set m%atom_type to contain the atom types; set m%composition to
        ! contain the percent composition (as a number between 0 and 1).
        allocate(m%atom_type(m%nelements), m%composition(m%nelements), stat=istat)
        if(istat /= 0) then
            write (*,*) 'Unable to allocate memory for m%atom_type and m%composition.'
            return
        endif
        ! Initialize the composition array to 0.0
        m%composition = 0.0
        ! i corresponds to the atomic number.
        ! j is the next open position in composition and atom_type.
        j = 1
        do i=1, 103
            if(elements(i) /= 0) then
                ! If we reach a non-zero element in elements then there are
                ! atoms with atomic number i in the model. Append this atomc
                ! number to atom_types and calculate the fractional composition
                ! of this element in the model, storing it in m%composition.
                ! Increment j to move to the next open spot.
                m%atom_type(j) = i
                m%composition(j) = real(elements(i)) / real(m%natoms)
                j = j + 1
            end if
        end do
        ! Note that m%atom_type and m%composition are now linked. You can look
        ! at m%composition, but you still won't know which element has this
        ! fractional composition. You need to look at the same index in
        ! m%atom_type in order to figure this out.

        ! Sort atom_type by increasing atomic order. Re-order composition in the
        ! same way so the indices stay in sync. (Insertion sort)
        do i=1, m%nelements
            do j=1, i
                if( m%atom_type(i) < m%atom_type(j) ) then
                    atom_temp = m%atom_type(i)
                    comp_temp = m%composition(i)
                    m%atom_type(i) = m%atom_type(j)
                    m%composition(i) = m%composition(j)
                    m%atom_type(j) = atom_temp
                    m%composition(j) = comp_temp
                end if
            end do
        end do

        ! For each atom i, add a parameter znum_r(i) that corresponds to
        ! m%atom_type and m%composition for fast lookup.
        allocate(m%znum_r(m%natoms), stat=istat)
        if(istat /= 0) then
            write (*,*) 'Unable to allocate memory for m%znum_r.'
            return
        endif
        m%znum_r = 0.0
        do i=1, m%natoms
            do j=1, m%nelements
                if(m%znum(i) .eq. m%atom_type(j)) then
                    m%znum_r(i) = j
                end if
            end do
        end do

        m%rotated = .FALSE.

        call recenter_model(0.0, 0.0, 0.0, m)

        ! Calls hutch_position and hutch_add_atom in loops.
        ! It does some allocation too.
        call model_init_hutches(m, istat)
    end subroutine read_model

    subroutine recenter_model(xc, yc, zc, m)
    ! Shifts the atom positions in model so that the mid-point between the maximum and
    ! and minimum atom positions in each dimensions sits at the position (xc, yc, zc),
    ! measured in units of the model supercell.
    ! BIGO(natoms)*6
        real, intent(in) :: xc, yc, zc
        type(model), intent(inout) :: m
        real :: xshift, yshift, zshift
        ! maxval calculates the maximum value in the array. There are some
        ! nice parameters for it described online by the way.
        xshift = xc*m%lx - (maxval(m%xx) + minval(m%xx))/2.0
        yshift = yc*m%ly - (maxval(m%yy) + minval(m%yy))/2.0
        zshift = zc*m%lz - (maxval(m%zz) + minval(m%zz))/2.0

        m%xx = m%xx+xshift
        m%yy = m%yy+yshift
        m%zz = m%zz+zshift
    end subroutine recenter_model

    subroutine model_init_hutches(m, status)
    ! Initializes the hutch_array ha within the model m. It calcualtes the
    ! hutch_size (based on the model box size and the parameter
    ! ATOMS_PER_HUTCH) and the number of hutches in the array (nhutch_x, nhutch_y,
    ! and hhutch_z). It then assigns all the atoms in the current model atom
    ! position arrays xa, ya, and za to the appropriate hutches.  It does NOT
    ! check whether ha has already been initialized, so this routine should
    ! NEVER be called more than once for the same hutch_array.
        type(model), intent(inout) :: m
        integer, intent(out) :: status
        integer :: istat, numhutches, hx, hy, hz, i
        ! Note: numhutches is not the total number of hutches, it is the number
        ! of hutches in each dimension. So numhutches^3 is the total.

        status = 0

        ! Jason rewrote this because it was compilcated. I made
        ! sure the current and the previous are mathematically equivalent.
        numhutches = anint( (m%natoms/ATOMS_PER_HUTCH)**(1./3.) )
        m%ha%hutch_size = m%lx / numhutches 
        !write (*,*) 'Hutch size is ',m%ha%hutch_size,' Angstroms.'
        !write (*,*) 'Number of hutch in each dimension is: ', numhutches

        m%ha%nhutch_x = numhutches 
        m%ha%nhutch_y = numhutches 
        m%ha%nhutch_z = numhutches 

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

        ! These hutch atom arrays are allocated and initialized in
        ! hutch_add_atom. We just need to initialize them to empty 
        ! and nat to 0 so that we can add atoms to them correctly.
        do hx = 1, m%ha%nhutch_x
            do hy = 1, m%ha%nhutch_y
                do hz = 1, m%ha%nhutch_z
                    ! I don't think we need the line below.
                    nullify(m%ha%h(hx, hy, hz)%at)
                    m%ha%h(hx, hy, hz)%nat = 0
                end do
            end do
        end do

        ! Calculate which hutch each atom should be in and add it to that hutch.
        do i=1, m%natoms
            call hutch_position(m, m%xx(i), m%yy(i), m%zz(i), hx, hy, hz)
            call hutch_add_atom(m, i, hx, hy, hz)
        end do
    end subroutine model_init_hutches

    subroutine hutch_position(m, xx, yy, zz, hx, hy, hz)
    ! returns the indices of the hutch that encompasses position (xx, yy, zz) in
    ! the hutch_array in the integers (hx, hy, hz).  It assumes that the model 
    ! extends from -lx/2 to lx/2, -ly/2 to ly/2 and -lz/2 to lz/2 and does no
    ! error checking.
        type(model), intent(in) :: m 
        real, intent(in) :: xx, yy, zz
        integer, intent(out) :: hx, hy, hz

        ! This makes the range of hx, hy, and hz from 0 to nhutch_i, however
        ! the only time one of them will be 0 is if the position is exactly on
        ! the left edge. Thats what the next set of if statements is for: If 
        ! they are on an edge just move them a little bit over. Technically you 
        ! can get statistically more atoms in hutches on the 3 "left" edges but 
        ! it will happen extremely rarely so it wont matter. By the time we 
        ! are done hx, hy, and hz are restrained from 1 to nhutch_i so we can 
        ! convienently put them in an array.
        hx = ceiling( (xx + 0.5*m%lx) / m%ha%hutch_size )
        hy = ceiling( (yy + 0.5*m%ly) / m%ha%hutch_size )
        hz = ceiling( (zz + 0.5*m%lz) / m%ha%hutch_size )

        if (hx == 0) hx = 1
        if (hy == 0) hy = 1
        if (hz == 0) hz = 1

        ! Jason 20130722 I commented this out because I think I was wrong above.
        ! I dont think these are necessary. I want to find out why someone
        ! thought they were too. I could definitely be wrong here.
        ! Maybe due to rounding errors if an atom is on the far right edge.
        ! But if thats the case, then hx = m%ha%nhutch_x not hx = 1, etc.
        !if(xx .ge. m%lx/2.0) hx = 1
        !if(yy .ge. m%ly/2.0) hy = 1
        !if(zz .ge. m%lz/2.0) hz = 1
    end subroutine hutch_position


    subroutine hutch_add_atom(m, atom, hx, hy, hz)
    ! Adds the atom with index atom to the hutch_array in hutch hx, hy, hz.
        type(model), target, intent(inout) :: m
        integer, intent(in) :: atom, hx, hy, hz
        integer :: nat
        integer, dimension(m%ha%h(hx, hy, hz)%nat+1) :: scratch_atoms
        type(hutch_array), pointer :: ha
        ha => m%ha
        ! ha%h(hx,hy,hz)%nat is set to 0 in a do loop in model_init_hutches,
        ! slightly before this function is called for each atom.
        nat = ha%h(hx,hy,hz)%nat
        if(nat > 0) then
            scratch_atoms(1:nat) = ha%h(hx, hy, hz)%at
            scratch_atoms(nat+1) = atom
            ! Reallocate with new size
            deallocate(ha%h(hx,hy,hz)%at)
            allocate(ha%h(hx,hy,hz)%at(1:nat+1)) ! +1 for extra atom
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

    subroutine check_model(m, istat)
    ! simple error checking on a model.  Currently checks: are all the  
    ! atoms in the box?  Are all the atomic numbers between 1 and 103
    ! (the range for which Kirkland calculated electron scattering factors)
    ! More should be added as we think of it.
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

    subroutine rotate_model(phi, psi, theta, min, mrot, istat)
        ! rotates model min by angles phi, psi, theta and puts the results in mrot. min is unchanged.
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
        istat = 0
        call periodic_continue_model(3, 3, 3, min, mt, .FALSE., istat)
        if (istat /= 0) return

        allocate(orig_indices(mt%natoms), stat=istat)
        if (istat /= 0) then 
           write (*,*) 'Memory allocation failure in rotate_model.'
           return
        endif

        do i=1,mt%natoms   !Array loop was temporarily changed due to error in visual fortran - Jinwoo Hwang
            if(mod(i,min%natoms) .eq. 0)then
                orig_indices(i) = min%natoms
            else
                orig_indices(i) = mod(i,min%natoms)
            endif
            if((orig_indices(i) .gt. min%natoms) .or. (orig_indices(i) .lt. 1)) then
                write(*,*) 'wrong here', i, orig_indices(i)
            endif
        enddo
        orig_indices = (/ (mod(i,min%natoms)+1, i=1,mt%natoms) /)

        ! generate the members of a 3x3 rotation matrix.  Use the Goldstein "x-convention"
        ! and Euler angles phi theta, psi.
        cpsi = cos(psi)
        cphi = cos(phi)
        ctheta = cos(theta)
        sphi = sin(phi)
        spsi = sin(psi)
        stheta = sin(theta)

        !phi ignored - JWH 09/02/09
        r(1,1) = cpsi
        r(1,2) = spsi
        r(1,3) = 0.0
        r(2,1) = -ctheta*spsi
        r(2,2) = ctheta*cpsi
        r(2,3) = stheta
        r(3,1) = stheta*spsi
        r(3,2) = -stheta*cpsi
        r(3,3) = ctheta

        ! Rotate the position vectors in mt (the temporary 3x3x3 model).
        do i=1,mt%natoms
            if(abs(mt%xx(i)).le.1.2*sqrt(2.0)*min%lx/2)then
                if(abs(mt%yy(i)).le.1.2*sqrt(2.0)*min%ly/2)then
                    if(abs(mt%zz(i)).le.1.2*sqrt(2.0)*min%lz/2)then
                        x = mt%xx(i)*r(1,1) + mt%yy(i)*r(1,2) + mt%zz(i)*r(1,3)
                        y = mt%xx(i)*r(2,1) + mt%yy(i)*r(2,2) + mt%zz(i)*r(2,3)
                        z = mt%xx(i)*r(3,1) + mt%yy(i)*r(3,2) + mt%zz(i)*r(3,3)
                        mt%xx(i) = x
                        mt%yy(i) = y
                        mt%zz(i) = z
                        !write(1008,*)i, mt%znum_r(i), mt%xx(i), mt%yy(i), mt%zz(i)
                    endif
                endif
            endif
        end do

        ! Cut the temporary model back to the original box size.
        ! First count the atoms in the box.
        mrot%natoms = 0
        lx2 = min%lx / 2.0
        ly2 = min%ly / 2.0
        lz2 = min%lz / 2.0
        do i=1, mt%natoms
            if((mt%xx(i) <= lx2 .AND. mt%xx(i) >= -1.0*lx2) .and. &
               (mt%yy(i) <= ly2 .AND. mt%yy(i) >= -1.0*ly2) .and. &
               (mt%zz(i) <= lz2 .AND. mt%zz(i) >= -1.0*lz2)) then
                mrot%natoms = mrot%natoms + 1
            !else
                !if(min%natoms /= 1425) write(*,*) "Atom outside world."
            endif
        enddo
        ! Allocate memory for the new atoms.
        mrot%unrot_natoms = min%natoms
        allocate(mrot%xx(mrot%natoms), mrot%yy(mrot%natoms), mrot%zz(mrot%natoms), &
            mrot%znum(mrot%natoms),  mrot%rot_i(mrot%unrot_natoms), mrot%znum_r(mrot%natoms), stat=istat) !add mrot%znum_r here by Feng Yi on 03/19/2009
        if( istat /= 0) then
           write (*,*) 'Problem allocating memory in rotate_model.'
           return
        endif

        do i=1,mrot%unrot_natoms
           mrot%rot_i(i)%nat = 0
           nullify(mrot%rot_i(i)%ind)
        enddo

        ! now copy just the atoms inside the original box size 
        ! from the temp model to the rotated one.
        j=1
        do i=1, mt%natoms
            if (mt%xx(i) <= lx2 .AND. mt%xx(i) >= -1.0*lx2) then
                if (mt%yy(i) <= ly2 .AND. mt%yy(i) >= -1.0*ly2) then
                    if (mt%zz(i) <= lz2 .AND. mt%zz(i) >= -1.0*lz2) then
                        mrot%xx(j) = mt%xx(i)
                        mrot%yy(j) = mt%yy(i)
                        mrot%zz(j) = mt%zz(i)
                        mrot%znum(j) = mt%znum(i)
                        mrot%znum_r(j) = mt%znum_r(i) !Added by Feng Yi on 03/19/2009   !Bug fixed : j to i -JWH 09/03/09
                        !write(*,*)"here",i, mt%znum_r(i), mt%xx(i),mt%yy(i),mt%zz(i)
                        ! add_index is basically just the general 
                        ! "append(list, element)" function except 
                        ! it takes a type object containing the list 
                        ! and an int equal to its size.
                        call add_index(mrot%rot_i(orig_indices(i)), j)
                        j = j+1
                    endif
                endif
            endif
        enddo

        !release the memory allocated to mt
        !call destroy_model(mt)   ! pmv 04/17/09
        deallocate(mt%atom_type, mt%composition)
        deallocate(mt%znum,mt%znum_r, mt%xx, mt%yy, mt%zz)
        ! TODO I think we need to destroy mt like pmv did above?

        ! set the rest of of the rotated model paramters
        mrot%lx = min%lx
        mrot%ly = min%ly
        mrot%lz = min%lz
        mrot%rotated = .TRUE.
        mrot%unrot_natoms = min%natoms

        if(mrot%natoms .ne. 0) then !added by jwh 03/26/2009
            call composition_model(mrot) ! have to recalculate this because the # of atoms may have changed a little
        endif

        call model_init_hutches(mrot, istat)

        if(allocated(orig_indices)) then !added by feng yi on 3/14/2009
            deallocate(orig_indices)
        endif
    end subroutine rotate_model

    subroutine destroy_model(m)
    ! Deallocates all the various allocatable arrays and sub-arrays in a model.
        type(model), intent(inout) :: m 
        deallocate(m%xx, m%yy, m%zz, m%znum, m%atom_type, m%znum_r, m%composition)
        call destroy_hutch(m%ha)
        call destroy_rot_indices(m%unrot_natoms, m%rot_i)
    end subroutine destroy_model

    subroutine destroy_hutch(ha)
    ! deallocates the hutch_array ha and the atom lists inside it.  used as part
    ! of destroy_model.
        type(hutch_array), intent(inout) :: ha
        integer i, j, k
        if(associated(ha%h)) then
            do i=1, ha%nhutch_x
                do j=1, ha%nhutch_y
                    do k=1, ha%nhutch_z
                        if (ha%h(i,j,k)%nat .gt. 0) then !added by feng yi on 03/14/2009
                            deallocate(ha%h(i,j,k)%at)
                        endif
                    enddo
                enddo
            enddo
        deallocate(ha%h, ha%atom_hutch)
        endif !if associated(ha%h)
    end subroutine destroy_hutch

    subroutine destroy_rot_indices(unrot_natoms, ri)
    ! deallocates all of the allocatable arrays and sub-arrays in an index list.
    ! used by destroy_model
        integer, intent(in) :: unrot_natoms
        type(index_list), pointer, dimension(:) :: ri
        integer i
        do i=1, unrot_natoms
            if(ri(i)%nat .gt. 0) then  !added by feng yi
                deallocate(ri(i)%ind)
            endif
        enddo
        if(associated(ri))then   !JWH - 042109
          deallocate(ri)
        endif
    end subroutine destroy_rot_indices

    subroutine hutch_list_pixel(m, px, py, diameter, atoms, istat)
        type(model), target, intent(in) :: m
        real, intent(in) :: px, py, diameter
        integer, pointer, dimension(:) :: atoms
        integer, intent(out) :: istat
        integer :: hx, hy, hz   ! hutch of position (px, py, pz)
        integer :: nh           ! number of hutches corresponding to diameter
        integer :: nlist        ! number of atoms in list
        integer :: i, j, k  ! counting variables
        integer, dimension(:), allocatable, target :: temp_atoms
        real, dimension(3) :: hcenter
        real :: dist
        integer :: i_start, i_end, j_start, j_end

        !write(*,*) "Number of hutches in the x, y, and z directions:", ha%nhutch_x, ha%nhutch_y, ha%nhutch_z

        if (.not. hlist_2d_calc) then
            call pre_calc_2d_hutch(m)
            hlist_2d_calc = .TRUE.
        endif

        ! Allocatae temp_atoms with the max number of atoms so that no matter
        ! how many we find, there will always be enough room.
        allocate(temp_atoms(m%natoms), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel'
            return
        end if

        ! I am going to do a slight approximation in this function, but it will
        ! be dang close. Considering the hutches are currently so small and
        ! contain only an atom or two, the additional hutches that will be
        ! included are not detrimental.
        ! The idea is to iterate through each hutch, calculate its center,
        ! and compute the distance from its center to (px,py) in the x-y plane;
        ! if this distance is <= diameter/2 + m%ha%hutch_size/sqrt(2.0) then we
        ! include that hutchs atoms. The above sum is the sum of the radius of
        ! the area we want to include + the "radius" (half diagonal) of the
        ! hutch. The half diagonal of the hutch may be a bit of an
        ! overapprximation, but it isnt much of one.

        ! I would like to update this function to pre-calculate the bounds on
        ! the do loops, but I dont know how due to the curvature of the circle.
        ! I could do the exact same thing as in hutch_list_pixel_sq with the
        ! bounds, but then also check in the if statement. This would reduce the
        ! iterations in the do loop, but the if statement would still be
        ! necessary. Doing this (which I have done) gives us a circle incribed
        ! in a square, and we need to exclude the hutches not in the circle but
        ! still in the square.

        ! Jason 20130712
        ! Precalcuate the hutches that are in range of a square with side length
        ! diameter. This is mathematically correct, but not intiutive, 
        ! unfortunately for the reader. This is also identical to the bounds in
        ! hutch_list_pixel_sq.
        i_start = ceiling( ( px - diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        i_end =   ceiling( ( px + diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        j_start = ceiling( ( py - diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        j_end =   ceiling( ( py + diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        !write(*,*) "i_start, i_end=", i_start, i_end
        !write(*,*) "j_start, j_end=", j_start, j_end
        !nh = (i_end-i_start+1)*(j_end-j_start+1)*(m%ha%nhutch_z) ! This wont work in the function.

        nh = 0
        nlist = 1
        do i=i_start, i_end
            do j=j_start, j_end
                do k=1, m%ha%nhutch_z
                    ! Calculate hutch centers.
                    hcenter(1) = -m%lx/2.0 + m%ha%hutch_size/2.0 + i*m%ha%hutch_size
                    hcenter(2) = -m%ly/2.0 + m%ha%hutch_size/2.0 + j*m%ha%hutch_size
                    hcenter(3) = -m%lz/2.0 + m%ha%hutch_size/2.0 + k*m%ha%hutch_size
                    ! Calculate distance.
                    dist = sqrt( (px-hcenter(1))**2 + (py-hcenter(2))**2 )
                    if( dist < diameter/2.0 + m%ha%hutch_size/sqrt(2.0) ) then
                        call hutch_position(m, hcenter(1), hcenter(2), hcenter(3), hx, hy, hz)
                        if(m%ha%h(hx, hy, hz)%nat /= 0) then
                            temp_atoms(nlist:nlist+m%ha%h(hx, hy, hz)%nat-1) = m%ha%h(hx, hy, hz)%at(1:m%ha%h(hx, hy, hz)%nat)
                            nlist = nlist + m%ha%h(hx, hy, hz)%nat
                        endif
                        nh = nh + 1
                    endif
                enddo
            enddo
        enddo


        ! Copy all the atoms we found in the previous loop into atoms.
        if (nlist > 1) then
            allocate(atoms(nlist-1), stat=istat)
            if (istat /= 0) then
                write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
                return
            endif
            atoms = temp_atoms(1:nlist-1)
        else
            nullify(atoms)
            istat = -1
        endif

        deallocate(temp_atoms)

        write(*,*) "pixel (", px,py, ") has diameter", diameter, "and contains", nlist, "atoms and ", nh, &
            "hutches !<= ", ( (ceiling(diameter/m%ha%hutch_size)+1) * (ceiling(diameter/m%ha%hutch_size)+1) * 11 ) ! debug

    end subroutine hutch_list_pixel

    subroutine hutch_position_eff(m, xx, yy, zz, hx, hy, hz,p_relative_3D, p_relative_2D)
    ! This must be called after the subroutine pre_calc_2D_hutch is called.
    ! It returns the indices of the hutch that encompasses position (xx, yy,
    ! zz) in the hutch_array in the integers (hx, hy, hz).  It assumes that the
    ! model extends from -lx/2 to lx/2, -ly/2 to ly/2 and -lz/2 to lz/2 and does no
    ! error checking. It also return relative position of point (xx, yy, zz) in the hutch
    ! in 2D and 3D cases the hutch is divided into equal 8 parts 2*2*2
        type(model), intent(in) :: m
        real, intent(in) :: xx, yy, zz
        integer, intent(out) :: hx, hy, hz
        integer, intent(out) :: p_relative_2d, p_relative_3d
        real,dimension(3) :: hutch_center !hutch center
        integer n1, n2, n3
        integer num_x, num_y, num_z

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

        if((xx + 0.5*m%lx) .ge. hutch_center(1)) then
            n3 = 2
        else
            n3=1
        endif
        if((yy + 0.5*m%ly) .ge. hutch_center(2)) then
            n2 = 2
        else
            n2=1
        endif
        if((zz + 0.5*m%lz) .ge. hutch_center(3)) then
            n1 = 2
        else
            n1=1
        endif
        p_relative_2d = (n2-1) * num_x + n3
        p_relative_3D = (n1-1) * num_y * num_x + (n2-1) * num_x + n3
    end subroutine hutch_position_eff

    subroutine pre_calc_2d_hutch(m)
    !pre-calculate the relative hutch position with respect to the center hutch
    !consider the worst condition for 2d case
    !use private variable list_1_2d
    !it is called only once. once, notice!!!
        type(model), target, intent(in) :: m
        real, allocatable,dimension(:) :: ratio_list
        !ratio_list from 0 to 0.5, 0.5 to 1.5, then step size with 1
        !ratio_list(1)=0.5, ratio_list(2)=1.5, ratio_list(3)=2.5, etc.
        real size_ref !square size
        integer num_x, num_y
        integer i1,j1,k1
        integer(selected_int_kind(3)) n1,n2
        integer(selected_int_kind(3)) nh
        integer num_temp1
        integer num_temp2
        integer(selected_int_kind(3)),allocatable, dimension(:, :) :: temp_list
        real r_dist, ref_dist
        real p_x, p_y
        real x_square, y_square

        size_ref = m%ha%hutch_size
        num_temp1 = ceiling(m%lx/2/size_ref)  !changed by feng yi on 03/25/2005

        if( .not. allocated(ratio_list)) then
            allocate(ratio_list(num_temp1))
        endif
        ratio_list(1) = 0.5
        do i1=2, num_temp1
            ratio_list(i1) = ratio_list(i1-1) + 1
        enddo
        num_x=2
        num_y=2

        list_1_2d%size_ratio = size(ratio_list,1)
        list_1_2d%size_position =  num_x*num_y

        allocate(list_1_2d%list_2d(list_1_2d%size_position ,list_1_2d%size_ratio))
        allocate(list_1_2d%ratio_radius_square(list_1_2d%size_ratio))
        list_1_2d%ratio_radius_square = ratio_list

        num_temp2 = ceiling(ratio_list(num_temp1)) * 2 + 1

        if(.not. allocated(temp_list)) then
            allocate(temp_list(2, num_temp2 * num_temp2))  !this size can be modified, revised by feng yi on 2/18/2009
            !first dimension refers to x and y, 1 is x and 2 is y
            !second dimension refers to x or y relative position
        endif

        do i1=1, num_y
            do j1=1, num_x
            num_temp1 = (i1-1) * num_x + j1
                do k1=1, list_1_2d%size_ratio
                    num_temp2 = 0
                    !calculate which hutch is within distance r to the point
                    !consider the worst case
                    ref_dist = size_ref * ratio_list(k1)
                    nh = ceiling(ratio_list(k1))
                    do n1=-nh, nh
                        !***********x dimension***********
                        if(n1 .gt. 0) then
                            !pick rightmost corner of square
                            p_x=(j1-1.0)*0.5
                        elseif (n1 .lt. 0) then
                            !pick leftmost corner of square
                            p_x=(j1-2.0)*0.5
                        else
                            p_x=0
                        endif
                        do n2=-nh, nh
                            !***********y dimension***********
                            if(n2 .gt. 0) then
                                !pick rightmost corner of square
                                p_y=(i1-1.0)*0.5
                            elseif (n2 .lt. 0) then
                                !pick leftmost corner of square
                                p_y=(i1-2.0)*0.5
                            else
                                p_y=0
                            endif

                            !closest point in a hutch to the origin
                            !pick the reference hutch center as the origin
                            if(n1 .ne. 0) then
                                x_square = (n1-n1/abs(n1)*0.5)*size_ref
                            else
                                x_square = 0.0
                            endif

                            if(n2 .ne. 0) then
                                y_square = (n2-n2/abs(n2)*0.5)*size_ref
                            else
                                y_square = 0.0
                            endif

                            !calculate the closest distance between one part of the
                            !reference square
                            !and a square with relative number position n1 and n2
                            r_dist = sqrt((x_square-p_x)**2+(y_square-p_y)**2)

                            !compare r_dist with ref_dist
                            if(r_dist .le. ref_dist) then
                                num_temp2 = num_temp2 + 1
                                temp_list(1,num_temp2) = n1
                                temp_list(2,num_temp2) = n2
                            endif
                        enddo !n2
                    enddo !n1

                    list_1_2d%list_2d(num_temp1,k1)%size_d = num_temp2
                    allocate(list_1_2d%list_2d(num_temp1,k1)%list_x(num_temp2))
                    allocate(list_1_2d%list_2d(num_temp1,k1)%list_y(num_temp2))
                    list_1_2d%list_2d(num_temp1,k1)%list_x = temp_list(1, 1:num_temp2)
                    list_1_2d%list_2d(num_temp1,k1)%list_y = temp_list(2, 1:num_temp2)
                enddo ! k1
            enddo !j1
        enddo  !i1

        if (allocated(temp_list)) then
            deallocate(temp_list)
        endif

        hlist_2d_calc = .true.

        if(allocated(ratio_list)) then
            deallocate(ratio_list)
        endif
    end subroutine pre_calc_2d_hutch

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
        if(allocated(scratch)) then
        deallocate(scratch)
        endif
    end subroutine add_index


    subroutine remove_element(il, elem)
        type(index_list), intent(inout) :: il
        integer, intent(in) :: elem
        integer, dimension(:), allocatable :: scratch
        integer :: i
        allocate(scratch(il%nat-1))
        do i=1,il%nat
            if(il%ind(i) == elem) then
                scratch(1:i-1) = il%ind(1:i-1)
                scratch(i:il%nat-1) = il%ind(i+1:il%nat)
                exit
            endif
        enddo
        deallocate(il%ind)
        allocate(il%ind( il%nat-1 ))
        il%ind = scratch
        deallocate(scratch)
        il%nat = il%nat - 1
    end subroutine remove_element

    subroutine composition_model(m)
    ! Calculates the composition of the model and fills in nelements, atom_type,
    ! and composition.
        type(model), intent(inout) :: m
        integer, dimension(103) :: znum_list
        integer :: i, j, isnew
        integer temp1
        real temp2 ! for temporary storage

        m%nelements=1
        znum_list(1) = m%znum(1)
        do i=1,m%natoms
            isnew = 1
            do j=1,m%nelements
                ! If atom i's atomic number is already in the list, don't add
                ! its atomic number to the list again.
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

        ! Floating bubble method
        do i=1, (size(m%atom_type)-1)
            do j=1, (size(m%atom_type)-i)
                if(m%atom_type(j) .gt. m%atom_type(j+1)) then
                    temp1 = m%atom_type(j)
                    m%atom_type(j) = m%atom_type(j+1)
                    m%atom_type(j+1) = temp1

                    temp2 = m%composition(j)
                    m%composition(j) = m%composition(j+1)
                    m%composition(j+1) = temp2
                endif
            enddo
        enddo
    end subroutine composition_model

    subroutine periodic_continue_model(xp, yp, zp, min, mout, init_hutch, istat)
    ! Makes (xp, yp, zp) copies of the input model min and puts them in the output model
    ! mout.  Returns non-zero in istat is the memory can't be allocated.
        integer, intent(in):: xp, yp, zp
        type(model), intent(in) :: min
        type(model), intent(out) :: mout
        logical, intent(in) :: init_hutch
        integer, intent(out) :: istat
        integer :: i, j, k, c
        real :: shift_x, shift_y, shift_z

        mout%natoms = min%natoms*xp*yp*zp
        allocate(mout%xx(mout%natoms), mout%yy(mout%natoms), mout%zz(mout%natoms), &
             mout%znum(mout%natoms), mout%znum_r(mout%natoms),stat=istat) !modified by Feng Yi on 03/19/2009
        if(istat /= 0) then
            write (*,*) 'Error allocating memory for the periodic continued model.'
            return
        end if

        mout%lx = min%lx*real(xp)
        mout%ly = min%ly*real(yp)
        mout%lz = min%lz*real(zp)

        c=0
        do i = -(xp-1)/2, (xp-1)/2     !jwh fyi 040809
            shift_x = real(i)*min%lx
            do j = -(yp-1)/2, (yp-1)/2
                shift_y = real(j)*min%ly
                do k = -(zp-1)/2, (zp-1)/2
                    shift_z = real(k)*min%lz
                    mout%xx(c*min%natoms+1:(c+1)*min%natoms) = min%xx + shift_x
                    mout%yy(c*min%natoms+1:(c+1)*min%natoms) = min%yy + shift_y
                    mout%zz(c*min%natoms+1:(c+1)*min%natoms) = min%zz + shift_z
                    mout%znum(c*min%natoms+1:(c+1)*min%natoms) = min%znum
                    mout%znum_r(c*min%natoms+1:(c+1)*min%natoms) = min%znum_r  !added by Feng Yi on 03/19/2009
                    c = c+1
                end do
            end do
        end do

        mout%nelements = min%nelements
        allocate(mout%atom_type(mout%nelements), mout%composition(mout%nelements), stat=istat)
        if(istat /= 0) then
            write (*,*) 'Problem allocating memory in periodic_continue_model.'
            return
        endif
        mout%atom_type = min%atom_type
        mout%composition = mout%composition

        if(init_hutch) then
            call model_init_hutches(mout, istat)
            if(istat /= 0) then
                write (*,*) 'Cannot allocate memeory for the new hutch_array.'
                return
            endif
        endif
    end subroutine periodic_continue_model


    subroutine hutch_list_3D(m, px, py, pz, radius, atoms, istat, nlist)
    ! Makes list (in atoms) of the indices of atoms that are in a space of
    ! length radius, centered on the hutch containing  the point (px, py, pz),
    ! using the hutch array.  Useful for calculating G(r).  Returns 1 in istat
    ! if memory allocation fails and -1 if no atoms are found.
    ! more efficient than old algorithm
    ! use private variable list_1_3D

    ! Returns the atoms (aka atom indices) that are within radius radius of 
    ! position (px,py,pz) in the list 'atoms'. Stores in nlist the number of
    ! atoms in this radius (i.e. size(atoms)+1 because of the original atom).

    ! TODO Check this to make it doesnt return too many like the others.

        type(model), target, intent(in) :: m
        real, intent(in) :: px, py, pz
        real, intent(in) :: radius
        !integer,  dimension(:), pointer, intent(out) :: atoms
        integer,  dimension(:), pointer :: atoms
        integer, intent(out) :: nlist        ! number of atoms in list
        integer, intent(out) :: istat
        integer :: hx, hy, hz   ! hutch of position (px, py, pz)
        integer :: nh           ! number of hutches corresponding to radius
        integer :: i, j, k, k1    ! counting variables
        integer :: hi, hj, hk   ! counting variables with periodic boundary conditions
        integer, allocatable, dimension(:), target :: temp_atoms  ! temporary atoms list
        type(hutch_array), pointer :: ha
        integer p_relative_3d, p_relative_2d
        integer ratio_position
        real ratio1
        logical :: use_new_alg
        !integer i2, j2, k2 !consider the hutch across the box boundary

        ha => m%ha
        allocate(temp_atoms(m%natoms), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Cannot allocate index list in hutch_list_3D.'
            return
        endif

        call hutch_position_eff(m, px, py, pz, hx, hy, hz,p_relative_3D, p_relative_2D)

        ratio1 = radius / ha%hutch_size
        nh = ceiling(ratio1)

        if(.not. hlist_3d_calc) then
            call pre_calc_3d_hutch(m)
        endif

        use_new_alg = .FALSE.
        do i=1, list_1_3d%size_ratio
            if(i .eq. 1) then
                if((ratio1 .ge. 0) .and. (ratio1 .le.  list_1_3d%ratio_radius_hutch(i))) then
                    ratio_position = i
                    use_new_alg = .TRUE.
                    exit
                endif
            else
                if((ratio1 .ge. list_1_3d%ratio_radius_hutch(i-1)) .and. (ratio1 .le.  list_1_3d%ratio_radius_hutch(i))) then
                    ratio_position = i
                    use_new_alg = .TRUE.
                    exit
                endif
            endif
        enddo

        if(use_new_alg) then
            ! Using new algorithm.
            nlist = 1
            do k1=1, list_1_3d%list_3d(p_relative_3d, ratio_position)%size_d
                k = list_1_3D%list_3D(p_relative_3D, ratio_position)%list_z(k1)
                k = hz + k
                if (k > ha%nhutch_z) then
                    hk = k - ha%nhutch_z
                else if (k < 1) then
                    hk = k + ha%nhutch_z
                else
                    hk = k
                end if
                j = list_1_3D%list_3D(p_relative_3D, ratio_position)%list_y(k1)
                j = hy + j
                if (j > ha%nhutch_y) then
                    hj = j - ha%nhutch_y
                else if (j < 1) then
                    hj = j+ ha%nhutch_y
                else
                    hj = j
                end if
                i = list_1_3D%list_3D(p_relative_3D, ratio_position)%list_x(k1)
                i = hx + i
                !Periodic boundary condition
                if (i > ha%nhutch_x) then
                    hi = i - ha%nhutch_x
                else if (i < 1) then
                    hi = i + ha%nhutch_x
                else
                    hi = i
                end if
                if (ha%h(hi, hj, hk)%nat > 0) then
                    temp_atoms(nlist:(nlist+ha%h(hi,hj,hk)%nat-1)) = ha%h(hi,hj,hk)%at
                    nlist = nlist+ha%h(hi,hj,hk)%nat
                endif
            enddo !k1
        else
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
                            !PRINT *, 'nlist+ha%h(hi,hj,hk)%nat-1 is: ', !nlist+ha%h(hi,hj,hk)%nat-1 !for debug
                            temp_atoms(nlist:(nlist+ha%h(hi,hj,hk)%nat-1)) = ha%h(hi,hj,hk)%at
                            nlist = nlist+ha%h(hi,hj,hk)%nat
                        endif
                    end do
                end do
            end do
        endif ! use_new_alg

        allocate(atoms(nlist-1), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Unable to allocate memory for atom indices in hutch_list_3D.'
            return
        endif

        !assign atoms to the subset of temp_atoms that was filled in
        if (nlist > 1) then
            atoms = temp_atoms(1:nlist-1)
        else
            nullify(atoms)
            istat = -1
        endif

        if(associated(ha)) then
            nullify(ha)
        endif
        deallocate(temp_atoms) !added by feng yi on 03/03/2009
    end subroutine hutch_list_3d

    subroutine hutch_list_pixel_sq(m, px, py, diameter, atoms, istat)
    ! Makes a list of atom indices (in atoms) of the atoms in a rectangular
    ! prism with side length diameter in x and y, through the model thickness
    ! in z, centered on the hutch containing the point (px, py). Useful for
    ! calculating the FEM intensity at (px, py).  Returns 1 in istat if the
    ! atoms array cannot be allocated.
    ! Technically, if a hutch is partly within the region described above then
    ! all its atoms are put in 'atoms'. Is this what we want???
        type(model), target, intent(in) :: m
        real, intent(in) :: px, py, diameter
        integer, pointer, dimension(:) :: atoms !output of atom indices
        integer, intent(out) :: istat
        integer :: hi, hj, hk, hx, hy, hz
        integer :: nh           ! number of hutches corresponding to diameter
        integer :: nlist        ! number of atoms in list
        integer :: i, j, k      ! counting variables
        integer, dimension(:), allocatable, target :: temp_atoms
        integer :: i_start, i_end, j_start, j_end

        !write(*,*) "Number of hutches in the x, y, and z directions:", m%ha%nhutch_x, m%ha%nhutch_y, m%ha%nhutch_z
        allocate(temp_atoms(m%natoms), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel'
            return
        end if

call hutch_position(m, px-diameter/2.0, py-diameter/2.0, 0.0, i_start, j_start, hz)
call hutch_position(m, px+diameter/2.0, py+diameter/2.0, 0.0, i_end, j_end, hz)
        ! Jason 20130712
        ! Precalcuate the hutches that are in range. This is mathematically
        ! correct, but not intiutive, unfortunately for the reader.
        !i_start = ceiling( ( px - diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size + 1 )
        !i_end =   ceiling( ( px + diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size + 1 )
        !j_start = ceiling( ( py - diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size + 1 )
        !j_end =   ceiling( ( py + diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size + 1 )
        !i_start = ( ( px - diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        !i_end =   ( ( px + diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        !j_start = ( ( py - diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        !j_end =   ( ( py + diameter/2.0 + m%lx/2.0 ) / m%ha%hutch_size )
        write(*,*) "i_start, i_end=", i_start, i_end
        write(*,*) "j_start, j_end=", j_start, j_end
        nh = (i_end-i_start+1)*(j_end-j_start+1)*(m%ha%nhutch_z)
        
        ! Fill in the list.
        nlist = 1
        do i = i_start, i_end
            do j = j_start, j_end
                do k = 1, m%ha%nhutch_z
                    if(m%ha%h(i, j, k)%nat /= 0) then
                        temp_atoms(nlist:nlist+m%ha%h(i, j, k)%nat-1) = m%ha%h(i, j, k)%at(1:m%ha%h(i, j, k)%nat)
                        nlist = nlist + m%ha%h(i, j, k)%nat
                    endif
                enddo
            enddo
        enddo

        ! Assign atoms to the subset of temp_atoms that was filled in.
        if( nlist > 1 ) then
            allocate(atoms(nlist-1), stat=istat)
            if (istat /= 0) then
                write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
                return
            endif
            atoms = temp_atoms
        else
            nullify(atoms)
            istat = -1
        endif

        write(*,*) "pixel (", px,py, ") has diameter", diameter, "and contains", nlist, "atoms and ", nh, &
            "hutches !<= ", ( (ceiling(diameter/m%ha%hutch_size)+1) * (ceiling(diameter/m%ha%hutch_size)+1) * 11 ) ! debug

        if(allocated(temp_atoms)) deallocate(temp_atoms)

    end subroutine hutch_list_pixel_sq


    subroutine pre_calc_3d_hutch(m)
    !pre-calculate the relative hutch position with respect to the center hutch
    !consider the worst condition for 3d case
    !assign values for list_1_3d
    !it is only called once. once,notice!!!!!
    !type(hutch_3d_array) list1
        type(model), target, intent(in) :: m
        real,allocatable, dimension(:) :: ratio_list
        !ratio_list from 0 to 0.5, 0.5 to 1.5, then step size with 1
        !ratio_list(1)=0.5, ratio_list(2)=1.5, ratio_list(3)=2.5, etc.
        real size_ref !square size
        integer num_x, num_y, num_z
        integer i1,j1,k1, k2
        integer(selected_int_kind (3)) n1,n2, n3
        integer(selected_int_kind(3)) nh
        integer num_temp1
        !integer,allocatable, dimension(:, :) :: num_temp2
        integer num_temp2
        integer(selected_int_kind (3)),allocatable, dimension(:, :) :: temp_list
        real r_dist, ref_dist
        real p_x, p_y, p_z
        real x_hutch, y_hutch, z_hutch

        size_ref = m%ha%hutch_size
        num_temp1 = ceiling(m%lx/size_ref) !changed by feng yi on 03/25/2009
        if( .not. allocated(ratio_list)) then
            allocate(ratio_list(num_temp1))
        endif
        ratio_list(1) = 0.5
        do i1=2, num_temp1
            ratio_list(i1) = ratio_list(i1-1) + 1
        enddo
        num_x=2
        num_y=2
        num_z=2

        list_1_3d%size_ratio = size(ratio_list,1)
        list_1_3d%size_position =  num_x*num_y*num_z
        allocate(list_1_3d%list_3d(list_1_3d%size_position ,list_1_3d%size_ratio))
        allocate(list_1_3d%ratio_radius_hutch(list_1_3d%size_ratio))
        list_1_3d%ratio_radius_hutch = ratio_list

        i1 = 2*ceiling(ratio_list(list_1_3d%size_ratio)) + 1
        j1 = i1
        k1 = i1

        if(.not. allocated(temp_list)) then
            allocate(temp_list(3, i1*j1*k1))  !this size can be modified
            !first dimension refers to x and y, 1 is x and 2 is y
            !second dimension refers to x or y relative position
        endif

        ! The variables over these do loops are small (0-10) as far as I can tell.
        do i1=1, num_z
            do j1=1, num_y
                do k1=1,num_x
                    num_temp1 = (i1-1) * num_y * num_x + (j1-1) * num_x + k1
                    do k2=1, list_1_3d%size_ratio
                        num_temp2 = 0
                        !calculate which hutch is within distance r to the point
                        !consider the worst case
                        ref_dist = size_ref * ratio_list(k2)
                        nh = ceiling(ratio_list(k2))
                        do n1=-nh, nh
                        !***********x dimension***********
                            if(n1 .gt. 0) then
                                !pick rightmost corner of square
                                p_x=(k1-1.0)*0.5
                            elseif (n1 .lt. 0) then
                                !pick leftmost corner of square
                                p_x=(k1-2.0)*0.5
                            else
                                p_x=0
                            endif
                            do n2=-nh, nh
                                !***********y dimension***********
                                if(n2 .gt. 0) then
                                    !pick rightmost corner of square
                                    p_y=(j1-1.0)*0.5
                                elseif (n2 .lt. 0) then
                                    !pick leftmost corner of square
                                    p_y=(j1-2.0)*0.5
                                else
                                    p_y=0
                                endif
                                do n3=-nh, nh
                                    !***********z dimension***********
                                    if(n3 .gt. 0) then
                                        !pick rightmost corner of square
                                        p_z=(i1-1.0)*0.5
                                    elseif (n3 .lt. 0) then
                                        !pick leftmost corner of square
                                        p_z=(i1-2.0)*0.5
                                    else
                                        p_z=0
                                    endif

                                    !closest point in a hutch to the origin
                                    !pick the reference hutch center as the
                                    !origin
                                    if(n1 .ne. 0) then
                                        x_hutch = (n1-n1/abs(n1)*0.5)*size_ref
                                    else
                                        x_hutch = 0.0
                                    endif

                                    if(n2 .ne. 0) then
                                        y_hutch = (n2-n2/abs(n2)*0.5)*size_ref
                                    else
                                        y_hutch = 0.0
                                    endif

                                    if(n3 .ne. 0) then
                                        z_hutch = (n3-n3/abs(n3)*0.5)*size_ref
                                    else
                                        z_hutch = 0.0
                                    endif

                                    !calculate the closest distance between one
                                    !part of the
                                    !reference square
                                    !and a square with relative number position
                                    !n1 and n2
                                    r_dist = sqrt((x_hutch-p_x)**2+(y_hutch-p_y)**2+(z_hutch-p_z)**2)

                                    !compare r_dist with ref_dist
                                    if(r_dist .le. ref_dist) then
                                        num_temp2 = num_temp2 + 1
                                        temp_list(1,num_temp2) = n1
                                        temp_list(2,num_temp2) = n2
                                        temp_list(3,num_temp2) = n3
                                    endif
                                enddo !n3
                            enddo !n2
                        enddo !n1

                        list_1_3d%list_3d(num_temp1,k2)%size_d = num_temp2
                        allocate(list_1_3d%list_3d(num_temp1,k2)%list_x(num_temp2))
                        allocate(list_1_3d%list_3d(num_temp1,k2)%list_y(num_temp2))
                        allocate(list_1_3d%list_3d(num_temp1,k2)%list_z(num_temp2))
                        list_1_3d%list_3d(num_temp1,k2)%list_x = temp_list(1, 1:num_temp2)
                        list_1_3d%list_3d(num_temp1,k2)%list_y = temp_list(2, 1:num_temp2)
                        list_1_3d%list_3d(num_temp1,k2)%list_z = temp_list(3, 1:num_temp2)
                    enddo ! k2
                enddo !k1
            enddo !j1
        enddo  !i1

        hlist_3d_calc = .true.

        if (allocated(temp_list)) then
            deallocate(temp_list)
        endif

        if(allocated(ratio_list)) then
            deallocate(ratio_list)
        endif

    end subroutine pre_calc_3d_hutch

    subroutine hutch_move_atom(m, atom, xx, yy, zz)
    ! Moves the atom with index atom from its current hutch in hutch_array to
    ! to the hutch that encompasses position (xx, yy, zz).  Used to update the
    ! hutch_array for a Monte Carlo move of one atom.
        type(model), target, intent(inout) :: m
        integer, intent(in) :: atom
        real, intent(in) :: xx, yy, zz
        integer :: hx, hy, hz
        type(hutch_array), pointer :: ha
        ha => m%ha
        call hutch_remove_atom(m, atom)
        call hutch_position(m, xx, yy, zz, hx, hy, hz)
        call hutch_add_atom(m, atom, hx, hy, hz)
        ha%atom_hutch(atom, 1) = hx
        ha%atom_hutch(atom, 2) = hy
        ha%atom_hutch(atom, 3) = hz
    end subroutine hutch_move_atom


    subroutine reject_position(m, atom, xx_cur, yy_cur, zz_cur)
        type(model), intent(inout) :: m
        integer, intent(in) :: atom
        real, intent(in) :: xx_cur, yy_cur, zz_cur
        !The moved atom in the original model, m, should return to their old position
        !when the random move is rejected - JWH 03/05/09
        m%xx(atom) = xx_cur
        m%yy(atom) = yy_cur
        m%zz(atom) = zz_cur
    end subroutine reject_position


    subroutine hutch_remove_atom(m, atom)
    ! remove atom atom from its current hutch.  This reduces the number of atoms in hutch
    ! array by one, so it should only be used in conjunction with hutch_add_atom.
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

        if(ha%h(hx, hy, hz)%nat .gt. 1) then  !added by feng yi on 03/19/2009
            allocate(ha%h(hx,hy,hz)%at(ha%h(hx,hy,hz)%nat-1))
            j=1
            do i=1, ha%h(hx,hy,hz)%nat
                if (scratch_atoms(i) /= atom) then
                    ha%h(hx,hy,hz)%at(j) = scratch_atoms(i)
                    j=j+1
                end if
            enddo

            ha%h(hx,hy,hz)%nat = ha%h(hx,hy,hz)%nat-1
            ha%atom_hutch(atom,1) = 0
            ha%atom_hutch(atom,2) = 0
            ha%atom_hutch(atom,3) = 0
        else
            ha%h(hx,hy, hz)%nat = 0
        endif
    end subroutine hutch_remove_atom


end module model_mod
