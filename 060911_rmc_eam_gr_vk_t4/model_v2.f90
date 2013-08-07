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
        integer, dimension(:), allocatable :: at
        integer :: nat
    end type hutch

    ! derived data type for the hutch array: contains an array of hutches,
    ! plus supporting information
    type hutch_array
        ! array of hutch objects
        type(hutch), dimension(:,:,:), allocatable :: h
        ! number of hutches in x, y, and z
        integer :: nhutch_x, nhutch_y, nhutch_z
        ! physical size of a hutch in Angstroms
        real :: hutch_size
        ! list of the hutch indices for every atom. we don't use it so im
        ! getting rid of it. Jason 20130729
        integer, allocatable, dimension(:,:) :: atom_hutch
    end type hutch_array

    ! adjustable-size list of atom indices
    type index_list
        integer :: nat
        integer, allocatable, dimension(:) :: ind
    end type index_list

    type real_index_list
        integer :: nat
        real, allocatable, dimension(:) :: ind
    end type real_index_list


    ! list for a rotated model that is the original models natoms long.  each atom in the list
    ! list is itself a list of the indices in the new, rotated model of the atoms corresponding
    ! to the the original atom in the unrotated model.
    ! for a rotated model, a list of the indices in the rotated model corresponding to each atom
    ! index in the original, unrotated model.  the original, unrotated model's index long.
    ! Defined type for a structural model with atoms positions and a bunch of metadata
    type model
        integer :: natoms                              ! number of atoms in the model
        !real, allocatable, dimension(:) :: xx, yy, zz      ! atom positions in Angstroms
        type(real_index_list) :: xx, yy, zz      ! atom positions in Angstroms
        !integer, allocatable, dimension(:) :: znum, znum_r         ! atom atomic numbers, and reduced z numbners
        type(index_list) :: znum, znum_r         ! atom atomic numbers, and reduced z numbners
        real :: lx, ly, lz                             ! box size, in Angstroms
        integer :: nelements                           ! # of elements in the model
        integer, allocatable, dimension(:) :: atom_type    ! array listing atomic numbers present
        real, allocatable, dimension(:) :: composition     ! fractional composition in the order of atom_type
        type(hutch_array) :: ha                        ! hutch data structure
        logical :: rotated                             ! TRUE if model has been rotated, FALSE otherwise
        integer :: unrot_natoms
        type(index_list), dimension(:), allocatable :: rot_i ! list of which atoms in the rotated model correspond
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
        call check_allocation(istat, "Error in opening flie, "//model_filename)

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
        !TODO change nat to nat + ceiling(nat*0.02) in the following line.
        allocate(m%xx%ind(nat), m%yy%ind(nat), m%zz%ind(nat), m%znum%ind(nat), stat=istat)
        m%xx%nat = nat
        m%yy%nat = nat
        m%zz%nat = nat
        m%znum%nat = nat
        m%znum_r%nat = nat
        call check_allocation(istat, 'Unable to allocate memory for the model being read.')

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
            read(1,*) m%znum%ind(i),m%xx%ind(i),m%yy%ind(i),m%zz%ind(i)
            ! If this atom has atomic number z, then increment the z position in
            ! the array elements. This counts the number of each atom type we have.
            elements(m%znum%ind(i)) = elements(m%znum%ind(i)) + 1
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
        call check_allocation(istat, 'Unable to allocate memory for m%atom_type and m%composition.')
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
        allocate(m%znum_r%ind(m%natoms), stat=istat)
        call check_allocation(istat, 'Unable to allocate memory for m%znum_r.')
        m%znum_r%ind = 0.0
        do i=1, m%natoms
            do j=1, m%nelements
                if(m%znum%ind(i) .eq. m%atom_type(j)) then
                    m%znum_r%ind(i) = j
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
        xshift = xc*m%lx - (maxval(m%xx%ind) + minval(m%xx%ind))/2.0
        yshift = yc*m%ly - (maxval(m%yy%ind) + minval(m%yy%ind))/2.0
        zshift = zc*m%lz - (maxval(m%zz%ind) + minval(m%zz%ind))/2.0

        m%xx%ind = m%xx%ind+xshift
        m%yy%ind = m%yy%ind+yshift
        m%zz%ind = m%zz%ind+zshift
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
                    if(allocated(m%ha%h(hx,hy,hz)%at)) deallocate(m%ha%h(hx,hy,hz)%at)
                    m%ha%h(hx, hy, hz)%nat = 0
                end do
            end do
        end do

        ! Calculate which hutch each atom should be in and add it to that hutch.
        do i=1, m%natoms
            call hutch_position(m, m%xx%ind(i), m%yy%ind(i), m%zz%ind(i), hx, hy, hz)
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
        integer, dimension(:,:), allocatable :: temp_atom_hutch
        type(hutch_array), pointer :: ha
        ha => m%ha
        ! ha%h(hx,hy,hz)%nat is set to 0 in a do loop in model_init_hutches,
        ! slightly before this function is called for each atom.
        ! TODO I should be able to use add_index and remove_index in functions
        ! like this. That would be ideal.
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
        ! Create space if there isnt already.
        ! I am lucky that this array allocation works.
        if( size(ha%atom_hutch) / 3 < atom ) then
            allocate(temp_atom_hutch( size(ha%atom_hutch) / 3, 3))
            temp_atom_hutch = ha%atom_hutch
            deallocate(ha%atom_hutch)
            allocate(ha%atom_hutch(m%natoms, 3))
            ha%atom_hutch = temp_atom_hutch
            deallocate(temp_atom_hutch)
        endif
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

        xlen = maxval(m%xx%ind) - minval(m%xx%ind)
        ylen = maxval(m%yy%ind) - minval(m%yy%ind)
        zlen = maxval(m%zz%ind) - minval(m%zz%ind)

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

        if (minval(m%znum%ind) < 1) then 
            write (*,*) 'Minimum atomic number of ', minval(m%znum%ind, 1), 'is less than zero.'
            istat = 1
        end if

        if (maxval(m%znum%ind) > 103) then 
            write (*,*) 'Maximum atomic number of ', maxval(m%znum%ind, 1), 'is greater than 103.'
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
            if(abs(mt%xx%ind(i)).le.1.2*sqrt(2.0)*min%lx/2)then
                if(abs(mt%yy%ind(i)).le.1.2*sqrt(2.0)*min%ly/2)then
                    if(abs(mt%zz%ind(i)).le.1.2*sqrt(2.0)*min%lz/2)then
                        x = mt%xx%ind(i)*r(1,1) + mt%yy%ind(i)*r(1,2) + mt%zz%ind(i)*r(1,3)
                        y = mt%xx%ind(i)*r(2,1) + mt%yy%ind(i)*r(2,2) + mt%zz%ind(i)*r(2,3)
                        z = mt%xx%ind(i)*r(3,1) + mt%yy%ind(i)*r(3,2) + mt%zz%ind(i)*r(3,3)
                        mt%xx%ind(i) = x
                        mt%yy%ind(i) = y
                        mt%zz%ind(i) = z
                        !write(1008,*)i, mt%znum_r%ind(i), mt%xx%ind(i), mt%yy%ind(i), mt%zz%ind(i)
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
            if((mt%xx%ind(i) <= lx2 .AND. mt%xx%ind(i) >= -1.0*lx2) .and. &
               (mt%yy%ind(i) <= ly2 .AND. mt%yy%ind(i) >= -1.0*ly2) .and. &
               (mt%zz%ind(i) <= lz2 .AND. mt%zz%ind(i) >= -1.0*lz2)) then
                mrot%natoms = mrot%natoms + 1
            !else
                !if(min%natoms /= 1425) write(*,*) "Atom outside world."
            endif
        enddo
        ! Allocate memory for the new atoms.
        mrot%unrot_natoms = min%natoms
        allocate(mrot%xx%ind(mrot%natoms), mrot%yy%ind(mrot%natoms), mrot%zz%ind(mrot%natoms), &
            mrot%znum%ind(mrot%natoms),  mrot%rot_i(mrot%unrot_natoms), mrot%znum_r%ind(mrot%natoms), stat=istat) !add mrot%znum_r here by Feng Yi on 03/19/2009
        mrot%xx%nat = mrot%natoms
        mrot%yy%nat = mrot%natoms
        mrot%zz%nat = mrot%natoms
        mrot%znum%nat = mrot%natoms
        mrot%znum_r%nat = mrot%natoms
        call check_allocation(istat, 'Problem allocating memory in rotate_model.')

        do i=1,mrot%unrot_natoms
           mrot%rot_i(i)%nat = 0
           if(allocated(mrot%rot_i(i)%ind)) deallocate(mrot%rot_i(i)%ind)
        enddo

        ! now copy just the atoms inside the original box size 
        ! from the temp model to the rotated one.
        j=1
        do i=1, mt%natoms
            if (mt%xx%ind(i) <= lx2 .AND. mt%xx%ind(i) >= -1.0*lx2) then
                if (mt%yy%ind(i) <= ly2 .AND. mt%yy%ind(i) >= -1.0*ly2) then
                    if (mt%zz%ind(i) <= lz2 .AND. mt%zz%ind(i) >= -1.0*lz2) then
                        mrot%xx%ind(j) = mt%xx%ind(i)
                        mrot%yy%ind(j) = mt%yy%ind(i)
                        mrot%zz%ind(j) = mt%zz%ind(i)
                        mrot%znum%ind(j) = mt%znum%ind(i)
                        mrot%znum_r%ind(j) = mt%znum_r%ind(i) !Added by Feng Yi on 03/19/2009   !Bug fixed : j to i -JWH 09/03/09
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
        ! I am assuming he commented out the above line and replaced it with the
        ! two below because he never created the hutch array or rot_i for mt.
        deallocate(mt%atom_type, mt%composition)
        deallocate(mt%znum%ind,mt%znum_r%ind, mt%xx%ind, mt%yy%ind, mt%zz%ind)

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
        deallocate(m%xx%ind, m%yy%ind, m%zz%ind, m%znum%ind, m%znum_r%ind)
        if(allocated(m%atom_type)) deallocate(m%atom_type)
        if(allocated(m%composition)) deallocate(m%composition)
        call destroy_hutch(m%ha)
        call destroy_rot_indices(m%unrot_natoms, m%rot_i)
    end subroutine destroy_model

    subroutine destroy_hutch(ha)
    ! Deallocates the hutch_array ha and the atom lists inside it.
    ! Used by destroy_model.
        type(hutch_array), intent(inout) :: ha
        integer i, j, k
        if(allocated(ha%h)) then
            do i=1,ha%nhutch_x
                do j=1,ha%nhutch_y
                    do k=1,ha%nhutch_z
                        if(ha%h(i,j,k)%nat .gt. 0) then
                            deallocate(ha%h(i,j,k)%at)
                        endif
                    enddo
                enddo
            enddo
            deallocate(ha%h, ha%atom_hutch)
        endif !if allocated(ha%h)
        ! I wonder if there is a memory leak because we dont actually delete the
        ! rest of the ha variables and ha itself. TODO
        ! This would occur for all the rot_atom models in fem_update.
    end subroutine destroy_hutch

    subroutine destroy_rot_indices(unrot_natoms, ri)
    ! Deallocates all of the allocatable arrays and sub-arrays in an index list.
    ! Used by destroy_model.
        integer, intent(in) :: unrot_natoms
        type(index_list), allocatable, dimension(:) :: ri
        integer :: i
        do i=1, unrot_natoms
            if(ri(i)%nat .gt. 0) then  !added by feng yi
                deallocate(ri(i)%ind)
            endif
        enddo
        if(allocated(ri)) deallocate(ri)
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
        integer :: i_start, i_end, j_start, j_end, trash

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

        ! Jason 20130722 Made this part much better.
        call hutch_position(m, px-diameter/2.0, py-diameter/2.0, 0.0, i_start, j_start, trash)
        call hutch_position(m, px+diameter/2.0, py+diameter/2.0, 0.0, i_end, j_end, trash)
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
            deallocate(atoms)
            istat = -1
        endif

        deallocate(temp_atoms)

        !write(*,*) "pixel (", px,py, ") has diameter", diameter, "and contains", nlist, "atoms and ", nh, &
            !"hutches !<= ", ( (ceiling(diameter/m%ha%hutch_size)+1) * (ceiling(diameter/m%ha%hutch_size)+1) * 11 ) ! debug

    end subroutine hutch_list_pixel

    subroutine hutch_position_eff(m, xx, yy, zz, hx, hy, hz,p_relative_3D, p_relative_2D)
    ! This must be called after the subroutine pre_calc_2D_hutch is called.
    ! It returns the indices of the hutch that encompasses position (xx, yy,
    ! zz) in the hutch_array in the integers (hx, hy, hz).  It assumes that the
    ! model extends from -lx/2 to lx/2, -ly/2 to ly/2 and -lz/2 to lz/2 and does no
    ! error checking. It also return relative position of point (xx, yy, zz) in the hutch
    ! in 2D and 3D cases the hutch is divided into equal 8 parts 2*2*2
    ! TODO Check if this function is necessary.
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
    ! TODO Figure out what the heck this function actually does.
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
        ! TODO Consider making this into a function that adds 10% of the size of
        ! il to il if we need to reallocate. This would reduce future needs to
        ! reallocate.
        ! Search for 'size(' to make sure nat is always used.
        type(index_list), intent(inout) :: il
        integer, intent(in) :: i
        integer, dimension(:), allocatable :: scratch
        if( il%nat >= 1 ) then
            ! If there is space no need to reallocate. If not, reallocate.
            if(size(il%ind) .ge. il%nat+1) then
                if(il%nat == -1) il%nat = 0 ! We set old_index(i) to -1 sometimes
                il%nat = il%nat + 1
                il%ind(il%nat) = i
            else
                allocate(scratch(il%nat))
                scratch = il%ind
                ! Increase the size by 2%. For the rot_i's this won't matter, but for
                ! the model lists it will increase them by a few atoms so that when
                ! we rotate in and out we don't need to reallocate every single time.
                ! 2% is barely a waste of space for the speed increase we will get.
                ! But only increment nat by 1. This is the size for the algorithm.
                il%nat = il%nat + 1
                deallocate(il%ind)
                !allocate(il%ind( il%nat + ceiling(il%nat*0.02) ))
                allocate(il%ind( il%nat + 1 )) !temp for debugging jason 20130729
                il%ind(1:il%nat-1) = scratch
                il%ind(il%nat) = i
            endif
        else
            il%nat = 1
            allocate(il%ind(1))
            il%ind(1) = i
        endif
        if(allocated(scratch)) then
            deallocate(scratch)
        endif
    end subroutine add_index

    subroutine add_index_real(il, i)
        ! TODO Consider making this into a function that adds 10% of the size of
        ! il to il if we need to reallocate. This would reduce future needs to
        ! reallocate.
        ! Search for 'size(' to make sure nat is always used.
        type(real_index_list), intent(inout) :: il
        real, intent(in) :: i
        integer, dimension(:), allocatable :: scratch
        if( il%nat >= 1 ) then
            ! If there is space no need to reallocate. If not, reallocate.
            if(size(il%ind) .ge. il%nat+1) then
                il%nat = il%nat + 1
                il%ind(il%nat) = i
            else
                allocate(scratch(il%nat))
                scratch = il%ind
                ! Increase the size by 2%. For the rot_i's this won't matter, but for
                ! the model lists it will increase them by a few atoms so that when
                ! we rotate in and out we don't need to reallocate every single time.
                ! 2% is barely a waste of space for the speed increase we will get.
                ! But only increment nat by 1. This is the size for the algorithm.
                il%nat = il%nat + 1
                deallocate(il%ind)
                !allocate(il%ind( il%nat + ceiling(il%nat*0.02) ))
                allocate(il%ind( il%nat + 1 )) !temp for debugging jason 20130729
                il%ind(1:il%nat-1) = scratch
                il%ind(il%nat) = i
            endif
        else
            il%nat = 1
            allocate(il%ind(1))
            il%ind(1) = i
        endif
        if(allocated(scratch)) then
            deallocate(scratch)
        endif
    end subroutine add_index_real


    subroutine remove_index(il, ind)
    ! TODO Consider not reallocating here unless there is a significant amount
    ! not being used. This would save time reallocating constantly.
    ! However, I would still need to do array deletion.
        type(index_list), intent(inout) :: il
        integer, intent(in) :: ind
        integer, dimension(:), allocatable :: scratch
        if( il%nat .gt. 1) then
            allocate(scratch(il%nat-1))
            ! First half
            scratch( 1:ind-1 ) = il%ind( 1:ind-1 )
            ! Second half
            scratch( ind:il%nat-1 ) = il%ind( ind+1:il%nat )
            deallocate(il%ind)
            allocate(il%ind( il%nat-1 ))
            il%ind = scratch
            deallocate(scratch)
            il%nat = il%nat - 1
        else
            deallocate(il%ind)
            il%nat = 0
        endif
    end subroutine remove_index

    subroutine remove_index_real(il, ind)
    ! TODO Consider not reallocating here unless there is a significant amount
    ! not being used. This would save time reallocating constantly.
    ! However, I would still need to do array deletion.
        type(real_index_list), intent(inout) :: il
        integer, intent(in) :: ind
        integer, dimension(:), allocatable :: scratch
        if( il%nat .gt. 1) then
            allocate(scratch(il%nat-1))
            ! First half
            scratch( 1:ind-1 ) = il%ind( 1:ind-1 )
            ! Second half
            scratch( ind:il%nat-1 ) = il%ind( ind+1:il%nat )
            deallocate(il%ind)
            allocate(il%ind( il%nat-1 ))
            il%ind = scratch
            deallocate(scratch)
            il%nat = il%nat - 1
        else
            deallocate(il%ind)
            il%nat = 0
        endif
    end subroutine remove_index_real


    subroutine composition_model(m)
    ! Calculates the composition of the model and fills in nelements, atom_type,
    ! and composition.
        type(model), intent(inout) :: m
        integer, dimension(103) :: znum_list
        integer :: i, j, isnew
        integer temp1
        real temp2 ! for temporary storage

        m%nelements=1
        znum_list(1) = m%znum%ind(1)
        do i=1,m%natoms
            isnew = 1
            do j=1,m%nelements
                ! If atom i's atomic number is already in the list, don't add
                ! its atomic number to the list again.
                if(m%znum%ind(i) == znum_list(j)) isnew = 0
            enddo
            if (isnew /= 0) then
                m%nelements=m%nelements+1
                znum_list(m%nelements) = m%znum%ind(i)
            endif
        enddo

        if( .not. allocated(m%atom_type) ) then
            allocate(m%atom_type(m%nelements))
        endif
        if( .not. allocated(m%composition) ) then
            allocate(m%composition(m%nelements))
        endif
        m%atom_type = znum_list(1:m%nelements)
        m%composition = 0.0

        do i = 1, m%natoms
            do j=1,m%nelements
                if(m%atom_type(j) == m%znum%ind(i)) then
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
        allocate(mout%xx%ind(mout%natoms), mout%yy%ind(mout%natoms), mout%zz%ind(mout%natoms), &
             mout%znum%ind(mout%natoms), mout%znum_r%ind(mout%natoms),stat=istat) !modified by Feng Yi on 03/19/2009
        mout%xx%nat = mout%natoms
        mout%yy%nat = mout%natoms
        mout%zz%nat = mout%natoms
        mout%znum%nat = mout%natoms
        mout%znum_r%nat = mout%natoms
        call check_allocation(istat, 'Error allocating memory for the periodic continued model.')

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
                    mout%xx%ind(c*min%natoms+1:(c+1)*min%natoms) = min%xx%ind + shift_x
                    mout%yy%ind(c*min%natoms+1:(c+1)*min%natoms) = min%yy%ind + shift_y
                    mout%zz%ind(c*min%natoms+1:(c+1)*min%natoms) = min%zz%ind + shift_z
                    mout%znum%ind(c*min%natoms+1:(c+1)*min%natoms) = min%znum%ind
                    mout%znum_r%ind(c*min%natoms+1:(c+1)*min%natoms) = min%znum_r%ind  !added by Feng Yi on 03/19/2009
                    c = c+1
                end do
            end do
        end do

        mout%nelements = min%nelements
        allocate(mout%atom_type(mout%nelements), mout%composition(mout%nelements), stat=istat)
        call check_allocation(istat, 'Problem allocating memory in periodic_continue_model.')
        mout%atom_type = min%atom_type
        mout%composition = mout%composition

        if(init_hutch) then
            call model_init_hutches(mout, istat)
            call check_allocation(istat, 'Cannot allocate memeory for the new hutch_array.')
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
            deallocate(atoms)
            istat = -1
        endif

        if(associated(ha)) then
            nullify(ha)
        endif
        deallocate(temp_atoms) !added by feng yi on 03/03/2009
    end subroutine hutch_list_3d

    subroutine hutch_list_pixel_sq(m, px, py, diameter, atoms, istat)
        type(model), target, intent(in) :: m
        real, intent(in) :: px, py, diameter
        integer, pointer, dimension(:) :: atoms !output of atom indices
        integer, intent(out) :: istat
        integer :: nh           ! number of hutches corresponding to diameter
        integer :: nlist        ! number of atoms in list
        integer :: i, j, k      ! counting variables
        integer, dimension(:), allocatable, target :: temp_atoms
        integer :: i_start, i_end, j_start, j_end, trash
        !logical :: found

        !write(*,*) "Number of hutches in the x, y, and z directions:", m%ha%nhutch_x, m%ha%nhutch_y, m%ha%nhutch_z
        allocate(temp_atoms(m%natoms), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel'
            return
        end if

        ! Jason 20130722 Made this part much better.
        call hutch_position(m, px-diameter/2.0, py-diameter/2.0, 0.0, i_start, j_start, trash)
        call hutch_position(m, px+diameter/2.0, py+diameter/2.0, 0.0, i_end, j_end, trash)
        !write(*,*) "i_start, i_end=", i_start, i_end
        !write(*,*) "j_start, j_end=", j_start, j_end
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
            deallocate(atoms)
            istat = -1
        endif

        !write(*,*) "pixel (", px,py, ") has diameter", diameter, "and contains", nlist, "atoms and ", nh, &
            !"hutches !<= ", ( (ceiling(diameter/m%ha%hutch_size)+1) * (ceiling(diameter/m%ha%hutch_size)+1) * 11 ) ! debug

        if(allocated(temp_atoms)) deallocate(temp_atoms)
!do i=1, m%natoms
!    found = .false.
!    do j=1, size(atoms)
!        if( atoms(j) == i ) then
!            found = .true.
!        endif
!    enddo
!    if( .not. found ) then
!        write(*,*) "HERE", i, m%xx(i), m%yy(i), m%zz(i)
!    endif
!enddo
!write(*,*) "size(atoms) = ", size(atoms)
!call sleep(60)
    end subroutine hutch_list_pixel_sq


    subroutine pre_calc_3d_hutch(m)
    !pre-calculate the relative hutch position with respect to the center hutch
    !consider the worst condition for 3d case
    !assign values for list_1_3d
    !it is only called once. once,notice!!!!!
    !type(hutch_3d_array) list1
    ! TODO Figure out what this does.
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
    ! to the hutch that encompasses position (xx, yy, zz). Used to update the
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
    end subroutine hutch_move_atom


    subroutine reject_position(m, atom, xx_cur, yy_cur, zz_cur)
        type(model), intent(inout) :: m
        integer, intent(in) :: atom
        real, intent(in) :: xx_cur, yy_cur, zz_cur
        !The moved atom in the original model, m, should return to their old position
        !when the random move is rejected - JWH 03/05/09
        m%xx%ind(atom) = xx_cur
        m%yy%ind(atom) = yy_cur
        m%zz%ind(atom) = zz_cur
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

        !call hutch_position(m, m%xx%ind(atom), m%yy%ind(atom), m%zz%ind(atom), hx, hy, hz)
        !write(*,*) "hx,hy,hz=", hx, hy, hz
        ! I wanted to get rid of atom_hutch, but if I accidently move the atom
        ! first before removing it from the hutch then I will royaly screw
        ! things up. So atom_hutch is more of a saftey feature right now, but it
        ! is a worthy one too. I may get rid of it later once I have all the
        ! bugs figured out. TODO

        hx = ha%atom_hutch(atom,1)
        hy = ha%atom_hutch(atom,2)
        hz = ha%atom_hutch(atom,3)
        !write(*,*) "hx,hy,hz=", hx, hy, hz

        scratch_atoms = ha%h(hx,hy,hz)%at
        deallocate(ha%h(hx,hy,hz)%at)

        if(ha%h(hx, hy, hz)%nat .gt. 1) then  !added by feng yi on 03/19/2009
            allocate(ha%h(hx,hy,hz)%at( ha%h(hx,hy,hz)%nat-1 ))
            j=1
            do i=1, ha%h(hx,hy,hz)%nat
                if (scratch_atoms(i) /= atom) then
                    ha%h(hx,hy,hz)%at(j) = scratch_atoms(i)
                    j=j+1
                end if
            enddo

            ha%h(hx,hy,hz)%nat = ha%h(hx,hy,hz)%nat-1
            ! I technically should reallocate here but I am going to do it in
            ! remove_atom because move_atom calls this, and if thats the case I
            ! dont need to reallocate.
            ha%atom_hutch(atom,1) = 0
            ha%atom_hutch(atom,2) = 0
            ha%atom_hutch(atom,3) = 0
        else
            ha%h(hx,hy, hz)%nat = 0
        endif
    end subroutine hutch_remove_atom


    subroutine move_atom(m, atom, xx, yy, zz)
        type(model), intent(inout) :: m
        integer, intent(in) :: atom
        real, intent(in) :: xx, yy, zz
        m%xx%ind(atom) = xx
        m%yy%ind(atom) = yy
        m%zz%ind(atom) = zz
        call hutch_move_atom(m, atom, xx, yy, zz)
    end subroutine move_atom

    subroutine add_atom(m, atom, xx, yy, zz, znum, znum_r)
        type(model), intent(inout) :: m
        integer, intent(in) :: atom ! index of atom in unroated model.
        real, intent(in) :: xx, yy, zz ! position of new atom
        integer, intent(in) :: znum, znum_r ! znum and znum_r of new atom
        integer :: hx, hy, hz
        ! We need to add an atom to xx, yy, zz, znum, znum_r, and the hutches.
        ! We need to increment natoms.
        ! We need to add a spot to rot_i with the correct index we used in xx, etc.
        ! We need to update the model's composition.
        ! We should check nelements and atom_type as well, but I am going to
        ! leave this out because it is so rare that we will remove the last atom
        ! of an atom type, and then re-add it later.
        !write(*,*) "A wild atom appeared!"

        ! We place the extra atom at the end of the above arrays, and therefore
        ! the new atom has index m%natoms+1.
        ! TODO Make sure if we remove an atom that every index gets updated
        ! correctly.
        ! Reallocate xx, yy, zz, znum, and znum_r bigger (leaving the
        ! end empty for the new atom to fit into).

        ! Add the atom to the model.
        call add_index(m%rot_i(atom), m%natoms + 1)
        call add_index_real(m%xx, xx)
        call add_index_real(m%yy, yy)
        call add_index_real(m%zz, zz)
        call add_index(m%znum, znum)
        call add_index(m%znum_r, znum_r)
        m%natoms = m%natoms + 1
        call hutch_position(m, xx, yy, zz, hx, hy, hz)
        ! This should give an out of bounds error on atom_hutch TODO
        ! I need to modify those functions to correctly increase the size of
        ! atom_hutch I think. But then I should be careful not to reallocate
        ! when the atom simply moves if I can help it... I don't want to do
        ! unnecessary reallocation.
        call hutch_add_atom(m, m%natoms, hx, hy, hz)

        ! Recalculate composition of model because it may have changed.
        call composition_model(m)
    end subroutine add_atom

    subroutine remove_atom(m, atom, ind)
        ! We need to remove atom ind from xx, yy, zz, znum, znum_r, and the hutches.
        ! We need to decrement natoms.
        ! We need to remove the spot from rot_i with the correct index we used in xx, etc.
        ! We need to update the model's composition.
        ! We should check nelements and atom_type as well, but I am going to
        ! leave this out because it is so rare that we will remove the last atom
        ! of an atom type.
        type(model), intent(inout) :: m
        integer, intent(in) :: atom ! index of atom in unroated model.
        integer, intent(in) :: ind ! index of atom to remove from m
        integer :: i, j, temp, hx, hy, hz
        integer, dimension(:,:), allocatable :: temp_atom_hutch
        !write(*,*) "An atom ran away!"

        temp = ind ! After I call remove_index on m%rot_i, ind is changed to 0
        ! since it is in an array. This seems like a fault of Fortran, but
        ! nevertheless I need to get around it by creating a new var temp and just
        ! setting it to ind before it gets changed. Then I will use temp after.
        ! I am assuming this happens because arrays are passed by reference, and
        ! since ind is part of an array, it is also passed by ref. This still
        ! shouldnt happen since I have it defined as an integer, intent(in).
        
        ! Remove ind from xx, yy, zz, znum, znum_r, and rot_i(atom).
        ! These calls decrement the index of all atoms with a higher index than
        ! ind. That is an inconvienence that we do need to deal with. It should
        ! only matter for rot_i hereafter, however, which we fix in the
        ! next forall loop.
        call remove_index_real(m%xx, ind)
        call remove_index_real(m%yy, ind)
        call remove_index_real(m%zz, ind)
        call remove_index(m%znum, ind)
        call remove_index(m%znum_r, ind)
        do i=1,m%rot_i(atom)%nat
            if(m%rot_i(atom)%ind(i) .eq. ind) then
                call remove_index(m%rot_i(atom), i)
                exit
            endif
        enddo
        ! Hereafter ind is not what it was before!

        ! Decrement every index in each rot_i that is higher than ind.
        ! We need to do this because we removed an element from each of the
        ! above arrays and therefore we need to correct the atom indices we
        ! are pointing to in rot_i.
        do i=1,m%unrot_natoms
            do j=1,m%rot_i(i)%nat
                if(m%rot_i(i)%ind(j) .gt. temp) then
                    m%rot_i(i)%ind(j) = m%rot_i(i)%ind(j) - 1
                endif
            enddo
        enddo

        ! Remove ind from its hutch.
        call hutch_remove_atom(m, temp)

        ! I also need to go through the hutches and decrement every index that
        ! is higher than ind for the same reason.
        ! Also, reallocate m%ha%atom_hutch here to one smaller.
        do i=temp+1, m%natoms
            hx = m%ha%atom_hutch(i,1)
            hy = m%ha%atom_hutch(i,2)
            hz = m%ha%atom_hutch(i,3)
            do j=1, m%ha%h(hx, hy, hz)%nat
                if( m%ha%h(hx, hy, hz)%at(j) .eq. i) then
                    m%ha%h(hx, hy, hz)%at(j) = m%ha%h(hx, hy, hz)%at(j) - 1
                endif
            enddo
        end do
        ! Reallocate m%ha%atom_hutch to one smaller.
        allocate(temp_atom_hutch( m%natoms, 3))
        temp_atom_hutch = m%ha%atom_hutch
        deallocate(m%ha%atom_hutch)
        allocate(m%ha%atom_hutch(m%natoms-1, 3))
        j = 1
        do i=1,m%natoms
            if( i /= temp) then
                m%ha%atom_hutch(j,1) = temp_atom_hutch(i,1)
                m%ha%atom_hutch(j,2) = temp_atom_hutch(i,2)
                m%ha%atom_hutch(j,3) = temp_atom_hutch(i,3)
                j = j + 1
            endif
        enddo
        deallocate(temp_atom_hutch)

        m%natoms = m%natoms - 1

        ! Recalculate composition of model because it may have changed.
        call composition_model(m)

    end subroutine remove_atom


    subroutine sort(il)
        ! Insertion sort on index list of ints
        type(index_list), intent(inout) :: il
        integer :: i, j, temp
        do i=1, il%nat
            do j=1, i
                if( il%ind(i) < il%ind(j) ) then
                    temp = il%ind(i)
                    il%ind(i) = il%ind(j)
                    il%ind(j) = temp
                end if
            end do
        end do
    end subroutine sort


    subroutine copy_model(m, mout)
        ! Copies model m to mout. If mout already contains information, it is
        ! deallocated and reallocated.
        ! Hopefully this function is faster than destroying m and re-rotating it
        type(model), intent(in) :: m
        type(model), intent(out) :: mout
        integer :: i, j, k
        
        mout%lx = m%lx
        mout%ly = m%ly
        mout%lz = m%lz

        mout%natoms = m%natoms
        if(allocated(mout%xx%ind)) deallocate(mout%xx%ind)
        if(allocated(mout%yy%ind)) deallocate(mout%xx%ind)
        if(allocated(mout%zz%ind)) deallocate(mout%xx%ind)
        if(allocated(mout%znum%ind)) deallocate(mout%znum%ind)
        if(allocated(mout%znum_r%ind)) deallocate(mout%znum_r%ind)
        allocate(mout%xx%ind(m%natoms), mout%yy%ind(m%natoms), mout%zz%ind(m%natoms), mout%znum%ind(m%natoms), mout%znum_r%ind(m%natoms))
        mout%xx%nat = m%xx%nat
        mout%yy%nat = m%yy%nat
        mout%zz%nat = m%zz%nat
        mout%znum%nat = m%znum%nat
        mout%znum_r%nat = m%znum_r%nat
        mout%xx%ind = m%xx%ind
        mout%yy%ind = m%yy%ind
        mout%zz%ind = m%zz%ind
        mout%znum%ind = m%znum%ind
        mout%znum_r%ind = m%znum_r%ind

        mout%nelements = m%nelements
        if(allocated(mout%atom_type)) deallocate(mout%atom_type)
        if(allocated(mout%composition)) deallocate(mout%composition)
        allocate(mout%atom_type(mout%nelements), mout%composition(mout%nelements))
        mout%atom_type = m%atom_type
        mout%composition = m%composition

        mout%rotated = m%rotated
        mout%unrot_natoms = m%unrot_natoms
        if(allocated(mout%rot_i)) then
            do i=1,m%unrot_natoms
                if(allocated(mout%rot_i(i)%ind)) then
                    deallocate(mout%rot_i(i)%ind)
                endif
            enddo
        endif
        if(allocated(mout%rot_i)) deallocate(mout%rot_i)
        allocate(mout%rot_i(mout%unrot_natoms))
        do i=1,m%unrot_natoms
            mout%rot_i(i)%nat = m%rot_i(i)%nat
            if(allocated(m%rot_i(i)%ind)) then
                allocate(mout%rot_i(i)%ind(mout%rot_i(i)%nat))
                mout%rot_i(i)%ind = m%rot_i(i)%ind
            endif
        enddo

        ! Destroy the old hutch, if necessary.
        call destroy_hutch(mout%ha)
       
        mout%ha%nhutch_x = m%ha%nhutch_x
        mout%ha%nhutch_y = m%ha%nhutch_y
        mout%ha%nhutch_z = m%ha%nhutch_z
        allocate(mout%ha%h(mout%ha%nhutch_x, mout%ha%nhutch_y,mout%ha%nhutch_y))
        do i=1,mout%ha%nhutch_x
            do j=1,mout%ha%nhutch_y
                do k=1,mout%ha%nhutch_z
                    mout%ha%h(i,j,k)%nat = m%ha%h(i,j,k)%nat
                    if(allocated(m%ha%h(i,j,k)%at)) then
                        if(allocated(mout%ha%h(i,j,k)%at)) then
                            deallocate(mout%ha%h(i,j,k)%at)
                        endif
                        allocate(mout%ha%h(i,j,k)%at(mout%ha%h(i,j,k)%nat))
                        mout%ha%h(i,j,k)%at = m%ha%h(i,j,k)%at
                    endif
                enddo
            enddo
        enddo
        mout%ha%hutch_size = m%ha%hutch_size
        if(allocated(m%ha%atom_hutch)) then
            if(allocated(mout%ha%atom_hutch)) then
                deallocate(mout%ha%atom_hutch)
            endif
            allocate(mout%ha%atom_hutch(mout%natoms, 3))
            mout%ha%atom_hutch = m%ha%atom_hutch
        else
            write(*,*) "There might be a problem. Atom_hutch should be associated for every input model."
        endif
    end subroutine copy_model


    subroutine check_allocation(istat, message)
        integer, intent(in) :: istat
        character(len=*), intent(in) :: message
        if (istat /= 0) then
            write (*,*) message
            return
        endif
    end subroutine check_allocation

end module model_mod
