!
!This model was written by Jason Maldonis from the old fem1.f90 file on
!06/29/13.
!


module fem_mod
    use  model_mod
    use  RMC_Global
    use  scattering_factors 

    implicit none
    private
    public :: fem_initialize, fem, I_average !femsim
    public:: fem_update, fem_accept_move, fem_reject_move !rmc
    public :: write_intensities
    !public :: print_image1, print_image2
    type pos_list
    integer :: nat
    real, pointer, dimension(:,:) :: pos
    end type pos_list
    integer, save :: nk, npix, nrot  ! number of k points, pixels, and rotations
    real, save, dimension(:,:), pointer :: pix ! npix x 2 list of pixel positions
    real, save, dimension(:,:), pointer :: rot ! nrot x 3 list of (phi, psi, theta) rotation angles
    real, save, dimension(:,:,:), pointer :: int_i, int_sq  ! nk x npix x nrot
    real, save, dimension(:,:,:), pointer :: old_int, old_int_sq
    real, save, dimension(:), pointer :: int_sum, int_sq_sum  ! nk long sums of int and int_sq arrays for calculating V(k)
    real, save, allocatable, dimension(:) :: j0, A1                                               
    type(model), save, dimension(:), pointer :: mrot  ! array of rotated models
    type(index_list), save, dimension(:), pointer :: old_index
    type(pos_list), save, dimension(:), pointer :: old_pos 

contains

    subroutine fem_initialize(m, res, k, nki, ntheta, nphi, npsi, scatfact_e, istat, square_pixel)
        type(model), intent(in) :: m 
        real, intent(in) :: res
        real, dimension(:), intent(in) :: k 
        integer, intent(in) :: nki, ntheta, nphi, npsi 
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        LOGICAL, OPTIONAL, INTENT(IN) :: square_pixel
        real dr      !Distance between pixels
        real r_max, const1, const2, const3
        integer npix_x, npix_y , bin_max
        integer i, j 
        integer const4 
        double precision b_x, b_j0 , b_j1 
        logical pixel_square

        if(present(square_pixel)) then 
            pixel_square = square_pixel
        else 
            pixel_square = .false.
        endif

        if(pixel_square) then 
            !r_max = SQRT(8.0) * res !diagonal in a square
            r_max = 2 * res !small pixel inscribed in Airy circle
        else 
            r_max = 2*res     !assuming resolution=radius
        endif
        !r_max = 6*res  !Check how cut-off affect v, 05/11/2009

        bin_max = int(r_max/fem_bin_width)+1

        const1 = twopi*(0.61/res)/fem_bin_width  !(0.61/res = Q) 
        const2 = 1/fem_bin_width
        const3 = (const1/(0.61/res))/const2
        const4 = int(bin_max*const3*CEILING(k(SIZE(k))))+1

        allocate(j0(0:const4),a1(0:const4), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Failed to allocate memory for Bessel and Airy functions.'
            return
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

        IF(pixel_square) THEN
            !dr = res ! condition 1, shift by res
            !npix_x=CEILING(m%lx/dr) 
            !npix_y=CEILING(m%ly/dr)  !integer of pixel number, shifted by dr
            !dr = m%lx/npix_x !changed by Feng Yi, fractional number of pixel
            !**********************
            !dr = res * 2.0  !no overlap between pixels
            !npix_x = ANINT(m%lx/dr)
            !npix_y = ANINT(m%ly/dr)
            !***********************
            dr = SQRT(2.0) * res !small pixel, inscribed in the Airy circle
            npix_x = ANINT(m%lx/dr)
            npix_y = ANINT(m%ly/dr)
            !dr = SQRT(2.0) * res/2.0 !small pixel, inscribed in the Airy
            !circle,shifted by half pixel size
            !npix_x = ANINT(m%lx/dr)
            !npix_y = ANINT(m%ly/dr)
        ELSE
            !dr=2*res/1.414214   !resolution=pixel radius. dr=pixel spacing
            !npix_x=int(m%lx/dr)
            !npix_y=int(m%ly/dr)
            !separated just by dr = res
            !dr = res
            dr = 2.0 *res
            npix_x=aint(m%lx/dr)
            npix_y=aint(m%ly/dr)
            !******************************
            !npix_x=CEILING(m%lx/dr) + 1 !changed by Feng Yi on 06/02/2009,
            !condition 1
            !npix_y=CEILING(m%ly/dr) + 1
            !dr = m%lx/(npix_x - 1.0) 
            !***************************************************************
            !npix_x=CEILING(m%lx/dr) ! changed by Feng Yi on 06/02/2009,
            !condition 2
            !npix_y=CEILING(m%ly/dr)
            !dr = m%lx/npix_x 
        ENDIF


        nk = nki
        npix = npix_x*npix_y
        !nrot = ntheta*nphi*npsi

        call init_rot(ntheta, nphi, npsi, nrot, istat)
        !if (istat /= 0) return

        allocate(int_i(nk, npix, nrot), old_int(nk, npix, nrot), old_int_sq(nk, npix, nrot), &
        int_sq(nk, npix, nrot), int_sum(nk), int_sq_sum(nk), stat=istat)
        nullify(old_index, old_pos)
        if (istat /= 0) then
            write (*,*) 'Cannot allocate memory in fem_initialize.'
            return
        endif

        call init_pix(m,npix_x, npix_y, dr, istat, pixel_square)
        !if (istat /= 0) return

        call read_f_e
        allocate(scatfact_e(m%nelements,nk), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Allocation of electron scattering factors table failed.'
            return
        endif

        do j=1,m%nelements
            do i=1, nk
                scatfact_e(j,i)=f_e(m%atom_type(j),k(i))
            enddo
        enddo

    end subroutine fem_initialize

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
5           CONTINUE
10          BJ1=1.0D0
            R=1.0D0
            DO 15 K=1,30
                R=-0.25D0*R*X2/(K*(K+1.0D0))
                BJ1=BJ1+R
                IF (DABS(R).LT.DABS(BJ1)*1.0D-15) GO TO 20
15          CONTINUE
20          BJ1=0.5D0*X*BJ1
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
25          CONTINUE
30          BY0=RP2*(EC*BJ0-CS0)
            CS1=1.0D0
            W1=0.0D0
            R1=1.0D0
            DO 35 K=1,30
                W1=W1+1.0D0/K
                R1=-0.25D0*R1/(K*(K+1))*X2
                R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
                CS1=CS1+R
                IF (DABS(R).LT.DABS(CS1)*1.0D-15) GO TO 40
35          CONTINUE
40          BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
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
45              Q0=Q0+B(K)*X**(-2*K-1)
            CU=DSQRT(RP2/X)
            BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
            BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
            T2=X-0.75D0*PI
            P1=1.0D0
            Q1=0.375D0/X
            DO 50 K=1,K0
                P1=P1+A1(K)*X**(-2*K)
50              Q1=Q1+B1(K)*X**(-2*K-1)
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

    subroutine init_rot(ntheta, nphi, npsi, num_rot, istat)
    ! Calculates the rotation angles and initializes them into the global
    ! rotation array rot. The rot_temp variable is probably unnecessary.
        integer, intent(in) :: ntheta, nphi, npsi
        integer, intent(out) :: istat
        integer :: i,j, k, jj, ii
        real,dimension(3) :: step_size
        integer :: ntheta_w(nphi*npsi)
        integer, intent(out) :: num_rot
        real, dimension(:,:), allocatable :: rot_temp
        real :: psi_temp
        integer :: pp

        allocate(rot_temp(ntheta*nphi*npsi,3), stat=istat)
        if (istat /= 0) then
           write (*,*) 'Cannot allocate temporary rotations array.'
           return
        endif

        !phi runs from 0 to 2 PI
        !psi runs from 0 to 2 PI
        !theta runs from 0 to PI   !not sure any more after weighting by psi angle - JWH 09/03/09
        !step_size(1) for phi step    
        !step_size(2) for psi step
        !step_size(3) for theta step
        step_size(1) = TWOPI / nphi
        step_size(2) = TWOPI / npsi
        step_size(3) = PI / ntheta  !not used any more after weighting by psi angle - JWH 09/03/09

        jj = 1
        do i=1, nphi
            do j=1, npsi/2
                psi_temp = (j-1)*step_size(2)
                ntheta_w(j) = int(sin(psi_temp)*ntheta)
                if(ntheta_w(j).ge.0)then
                    if(ntheta_w(j).gt.0)then
                        pp = 2*(ntheta_w(j)-1)
                    endif
                    if(ntheta_w(j).eq.0)then
                        pp = 1
                    endif
                    do k=1, pp
                        if(k*(pi/(ntheta_w(j)-1)).lt.pi)then
                            rot_temp(jj,1) = (i-1)*step_size(1)
                            rot_temp(jj,2) = (j-1)*step_size(2)
                            rot_temp(jj,3) = k*(pi/(ntheta_w(j)-1))
                            jj = jj + 1
                        endif
                    enddo
                endif
            enddo
        enddo

        num_rot = jj - 1

        allocate(rot(num_rot, 3), stat=istat)
        if (istat /= 0) then
           write (*,*) 'Cannot allocate rotations array.'
           return
        endif

        do i=1, num_rot
            rot(i,1) = rot_temp(i,1)
            rot(i,2) = rot_temp(i,2)
            rot(i,3) = rot_temp(i,3)
        enddo

        deallocate(rot_temp)

    end subroutine init_rot

    subroutine init_pix(m, npix_x, npix_y, dr, istat, square_pixel)
        type(model), intent(in) :: m
        integer, intent(in) :: npix_x, npix_y
        real :: dr
        integer, intent(out) :: istat
        logical, optional, intent(in) :: square_pixel
        integer :: i, j, k
        logical :: pixel_square

        allocate(pix(npix, 2), stat=istat)
        if (istat /= 0) then
           write (*,*) 'Cannot allocate pixel position array.'
           return
        endif

        IF(PRESENT(square_pixel)) THEN
            pixel_square = square_pixel
        ELSE
            pixel_square = .FALSE.
        ENDIF

        k=1
        IF(pixel_square) THEN
            !write(*,*)"pixels1"
            DO i=1, npix_x
                DO j=1, npix_y
                    pix(k,1) = (i-1)*dr+(dr/2.0)-(m%lx/2.0)
                    pix(k,2) = (j-1)*dr+(dr/2.0)-(m%ly/2.0)
                    !WRITE(*,*) k, pix(k,1), pix(k,2) !debug
                    k = k + 1
                ENDDO
            ENDDO
        ELSE
            !write(*,*)"pixels2"
            do i=1, npix_x
                do j=1, npix_y
                    pix(k,1) = (i-1)*dr+(dr/2.0)-(m%lx/2.0) !condition 2
                    pix(k,2) = (j-1)*dr+(dr/2.0)-(m%ly/2.0)
                    !*************************************
                    !pix(k,1) = (i-1)*dr-(m%lx/2.0)
                    !pix(k,2) = (j-1)*dr-(m%ly/2.0)  !changed by FY on 06/02/2009, condition 1
                    !write(*,*)pix(k,1),pix(k,2)
                    k = k+1
                enddo
            enddo
        ENDIF

        if(myid.eq.0)then
            write(*,*)"pixels=", npix_x, "by", npix_y
            k=1
            do i=1, npix_x
                do j=1, npix_y
                    write(*,*)pix(k,1), pix(k,2)
                    k=k+1
                enddo
            enddo
        endif
    end subroutine init_pix

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

    subroutine i_average(i_k)
        implicit none
        real ,dimension(:), intent(out) :: i_k
        integer i

        i_k = 0.0
        do i=1, nk
            i_k(i) = sum(int_i(i,1:npix,1:nrot))/(npix * nrot)
        enddo
    end subroutine i_average

    subroutine fem(m, res, k, vk, v_background, scatfact_e, comm, istat, square_pixel, use_femsim, rot_begin, rot_end)
        use mpi
        implicit none
        type(model), intent(in) :: m
        real, intent(in) :: res
        real, dimension(:), intent(in) :: k
        real, dimension(:), INTENT(OUT) :: Vk
        real, dimension(:), intent(in) :: v_background
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        logical, optional, intent(in) :: use_femsim
        logical, optional, intent(in) :: square_pixel
        integer, optional, intent(in) :: rot_begin, rot_end
        logical femsim !added by Feng Yi on 03/19/2009
        logical pixel_square !added by FY on 06/04/2009
        real, dimension(:), allocatable :: psum_int, psum_int_sq, sum_int, sum_int_sq  !mpi
        integer :: comm
        integer :: i, j, i1, j1
        integer begin_rot, end_rot
        real test1

        if( present(square_pixel)) then
            pixel_square = square_pixel
        else
            pixel_square = .FALSE.
        endif

        if( present(use_femsim)) then
            femsim = use_femsim
        else
            femsim = .FALSE.
        endif

        if( present(rot_begin)) then
            begin_rot = rot_begin
        else
            begin_rot  = 1
        endif

        if( present(rot_end)) then
            end_rot = rot_end
        else
            end_rot = nrot ! This is set in fem_initialize->init_rot
        endif

        if(femsim) then   !this is added by feng yi on 03/19/2009 only for femsim calculation

            write(*,*) "Rotating ", end_rot, " models."
            do i=begin_rot, end_rot
            ! initialize the rotated models
                allocate(mrot(1), stat=istat) !debug memory leak
                if (istat /= 0) then
                    write(*,*) 'Cannot allocate rotated model array.'
                    return
                endif
                call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m, mrot(1), istat) !memory leak
                if (istat /= 0) then
                    write (*,*) 'Failed to rotate model ',i
                    return
                endif

                ! Calculate intensities and store them in int_i(1:nk, j, i).
                do j=1, npix
                    call intensity(mrot(1), res, pix(j, 1), pix(j, 2), k, int_i(1:nk, j, i), scatfact_e, istat, pixel_square)
                enddo
                call destroy_model(mrot(1)) !memory leak
                deallocate(mrot) !memory leak
            enddo !end i=1, nrot
            write(*,*) "Model rotation completed."

            int_sq = int_i*int_i
            do i=1, nk
                Vk(i) = (sum(int_sq(i,1:npix,1:nrot)/(npix*nrot))/(sum(int_i(i,1:npix,1:nrot)/(npix*nrot))**2) ) - 1.0 !tempporally added by Feng Yi just for femsim
                !write(*,*) k(i), Vk(i), sum(int_i(i,1:npix,1:nrot)/(npix*nrot))
            enddo

        !***********************************************************
        ELSE        !RMC
        !*************************************

            allocate (psum_int(size(k)), psum_int_sq(size(k)), sum_int(size(k)), sum_int_sq(size(k)), stat=istat)
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

            ! I don't completely understand this loop. I thought what we wanted
            ! to do here was to calculate all the rotated models now. Are we???
            do i=myid+1, nrot, numprocs
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
            ! ??? I don't know what these are.
            do i=myid+1, nrot, numprocs
                old_index(i)%nat = 0
                nullify(old_index(i)%ind)
                old_pos(i)%nat = 0
                nullify(old_pos(i)%pos)
            enddo

            ! Calculate intensities. This is pretty expensive.
            do i=myid+1, nrot, numprocs
                do j=1, npix
                write(*,*) i, j
                    call intensity(mrot(i), res, pix(j, 1), pix(j, 2), k, int_i(1:nk, j, i), scatfact_e, istat, pixel_square)
                    int_sq(1:nk, j, i) = int_i(1:nk, j, i)**2
                    psum_int(1:nk) = psum_int(1:nk) + int_i(1:nk, j, i)
                    psum_int_sq(1:nk) = psum_int_sq(1:nk) + int_sq(1:nk, j, i)
                enddo
            enddo

            call mpi_reduce (psum_int, sum_int, size(k), mpi_real, mpi_sum, 0, comm, mpierr)
            call mpi_reduce (psum_int_sq, sum_int_sq, size(k), mpi_real, mpi_sum, 0, comm, mpierr)

            if(myid.eq.0)then
                do i=1, nk
                    Vk(i) = (sum_int_sq(i)/(npix*nrot))/((sum_int(i)/(npix*nrot))**2)-1.0
                    Vk(i) = Vk(i) - v_background(i)  ! background subtraction   052210 JWH
                end do
            endif

            deallocate(psum_int, psum_int_sq, sum_int, sum_int_sq)
        ENDIF !Femsim or RMC
    end subroutine fem

    subroutine intensity(m_int, res, px, py, k, int_i, scatfact_e, istat, square_pixel)
    ! Calculates int_i for output.
        type(model), intent(in) :: m_int
        real, intent(in) :: res, px, py
        real, dimension(nk), intent(in) :: k
        real, dimension(nk), intent(out) :: int_i
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        logical,optional, intent(in) :: square_pixel
        real, dimension(:,:,:), allocatable :: gr_i   ! unneeded 'save' keyword removed pmv 03/18/09  !tr re-ok -jwh
        real, dimension(:), allocatable ::x1, y1, rr_a
        real, dimension(:,:), allocatable :: sum1
        real :: x2, y2, rr, t1, t2, const1, const2, const3, pp, r_max
        integer, pointer, dimension(:) :: pix_atoms, znum_r
        integer :: i,j,ii,jj,kk
        integer :: bin_max
        logical pixel_square
        real, allocatable, dimension(:) :: rr_x, rr_y
        real sqrt1_2_res
        real k_1

        if(present(square_pixel)) then
            pixel_square = square_pixel
        else
            pixel_square = .FALSE.
        endif

        if(pixel_square) then
            sqrt1_2_res = SQRT(0.5) * res
        else
            sqrt1_2_res = res
        endif

        if(pixel_square) then
            !r_max = sqrt(8.0) * res
            r_max = 2 * res !small pixel inscribed in airy circle
        else
            r_max = 2*res     !assuming resolution=radius
        endif

        !r_max = 6*res    !check cut-off effect on 05/11/2009
        bin_max = int(r_max/fem_bin_width)+1

        if(pixel_square) then
            call hutch_list_pixel_sq(m_int, px, py, res, pix_atoms, istat) !small pixel inscribed in airy circle
        else
            call hutch_list_pixel(m_int, px, py, res, pix_atoms, istat)
        endif

        allocate(gr_i(m_int%nelements,m_int%nelements, 0:bin_max), stat=istat)
        allocate(x1(size(pix_atoms)),y1(size(pix_atoms)),rr_a(size(pix_atoms)), stat=istat)
        allocate(sum1(m_int%nelements,size(pix_atoms)), stat=istat)
        allocate(znum_r(size(pix_atoms)), stat=istat)

        if(pixel_square) then
        allocate( rr_x(size(pix_atoms)),rr_y(size(pix_atoms)), stat=istat)
        endif
        ! replaced code recalculating znum_r with code copying it from previous
        ! calculations 3/18/09 pmv  !tr RE-ok-jwh
        do i=1, size(pix_atoms)
            znum_r(i) = m_int%znum_r(pix_atoms(i))
        enddo

        gr_i = 0.0
        int_i = 0.0
        x1 = 0.0
        y1 = 0.0
        rr_a = 0.0
        x2 = 0.0
        y2 = 0.0

        const1 = twopi*(0.61/res)/fem_bin_width  !(0.61/res = Q)
        const2 = 1/fem_bin_width
        const3 = TWOPI

        ! Calculate sum1 for gr_i calculation in next loop.
        if(pixel_square) then
            do i=1,size(pix_atoms)
                x2=m_int%xx(pix_atoms(i))-px
                y2=m_int%yy(pix_atoms(i))-py
                x2=x2-m_int%lx*anint(x2/m_int%lx)
                y2=y2-m_int%ly*anint(y2/m_int%ly)
                rr_x(i) = ABS(x2)
                rr_y(i) = ABS(y2)
                rr_a(i)=sqrt(x2*x2 + y2*y2)
                !if((rr_x(i).le.res) .AND. (rr_y(i) .le. res))then
                if((rr_x(i) .le. sqrt1_2_res) .AND. (rr_y(i) .le.  sqrt1_2_res))then !small pixel inscribed in Airy circle
                    k_1=0.82333
                    x1(i)=x2
                    y1(i)=y2
                    j=int(const1*rr_a(i))
                    sum1(znum_r(i),i)=A1(j)
                endif
            enddo
        else
            do i=1,size(pix_atoms)
                x2=m_int%xx(pix_atoms(i))-px
                y2=m_int%yy(pix_atoms(i))-py
                x2=x2-m_int%lx*anint(x2/m_int%lx)
                y2=y2-m_int%ly*anint(y2/m_int%ly)
                rr_a(i)=sqrt(x2*x2 + y2*y2)
                if(rr_a(i).le.res)then
                    !if(rr_a(i) .le. res*3.0)then !check cut-off effect
                    x1(i)=x2
                    y1(i)=y2
                    j=int(const1*rr_a(i))
                    sum1(znum_r(i),i)=A1(j)
                endif
            enddo
        endif

        ! Calculate gr_i for int_i in next loop.
        if(pixel_square) then
            do i=1,size(pix_atoms)
                if((rr_x(i).le.sqrt1_2_res) .and. (rr_y(i) .le.  sqrt1_2_res))then
                    do j=i,size(pix_atoms)
                        if((rr_x(j).le.sqrt1_2_res) .and. (rr_y(j) .le.  sqrt1_2_res))then
                            x2=x1(i)-x1(j)
                            y2=y1(i)-y1(j)
                            rr=sqrt(x2*x2 + y2*y2)
                            kk=int(const2*rr)
                            if(i == j)then
                                t1=sum1(znum_r(i),i)
                                gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+t1*t1
                            else
                                t1=sum1(znum_r(i),i)
                                t2=sum1(znum_r(j),j)
                                gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+2.0*t1*t2 !changed by FY on 05/04/2009
                            endif
                        endif
                    enddo
                endif
            enddo
        else
            do i=1,size(pix_atoms)
                if(rr_a(i).le.res)then
                    !if(rr_a(i) .le. res*3.0)then  !check cut-off effect
                    do j=i,size(pix_atoms)
                        if(rr_a(j).le.res)then
                            !if(rr_a(j) .le. res*3.0)then  !check cut-off effect
                            x2=x1(i)-x1(j)
                            y2=y1(i)-y1(j)
                            rr=sqrt(x2*x2 + y2*y2)
                            kk=int(const2*rr)
                            if(i == j)then
                                t1=sum1(znum_r(i),i)
                                gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+t1*t1
                            else
                                t1=sum1(znum_r(i),i)
                                t2=sum1(znum_r(j),j)
                                gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+2.0*t1*t2 !changed by FY on 05/04/2009
                            endif
                        endif
                    enddo
                endif
            enddo
        endif !pixel_square

        ! Calculate int_i for output.
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

        if(allocated(gr_i)) then
            deallocate(gr_i)
        endif
        if(allocated(x1)) then
            deallocate(x1,y1, rr_a, znum_r)
        endif
        if(size(pix_atoms) .gt. 0) then
            deallocate(pix_atoms)
        endif
        if(allocated(sum1)) then
            deallocate(sum1)
        endif
        if(allocated(rr_x)) then
            deallocate(rr_x, rr_y)
        endif
    end subroutine intensity


    subroutine fem_update(m_in, atom, res, k, vk, v_background, scatfact_e, comm, istat,square_pixel)
        use mpi
        type(model), intent(in) :: m_in
        integer, intent(in) :: atom
        real, intent(in) :: res
        real, dimension(:), intent(in) :: k, v_background
        real, dimension(:), intent(out) :: vk
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        LOGICAL, OPTIONAL, INTENT(IN) :: square_pixel
        real, dimension(:), allocatable :: psum_int, psum_int_sq, sum_int, sum_int_sq    !mpi
        integer :: comm
        type(model) :: moved_atom, rot_atom
        integer :: i, j, m, n
        real :: res2, rot_dist_sq, orig_dist_sq, temp1, temp2
        REAL rr_x_old, rr_y_old, rr_x_new, rr_y_new
        LOGICAL pixel_square
        real :: sqrt1_2_res
        integer :: old_mrot_roti_nat
        logical :: no_int_recal
        istat = 0

        if( present(square_pixel)) then
            pixel_square = square_pixel
        else
            pixel_square = .FALSE.
        endif
        ! res is diameter of a pixel, but later we need the square of the radius
        res2 = (res)**2    !temporary - radius = resolution assumed - JWH 02/25/09

        ! Create a new model (moved_atom) with only one atom in it and put the
        ! position of the moved atom into it.
        allocate(moved_atom%xx(1), moved_atom%yy(1), moved_atom%zz(1), moved_atom%znum(1), &
        moved_atom%atom_type(1), moved_atom%znum_r(1), moved_atom%composition(1), stat=istat)
        allocate(psum_int(size(k)), psum_int_sq(size(k)), sum_int(size(k)), sum_int_sq(size(k)), stat=istat)
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

        ! update each rotated model with the position of the moved atom, then
        ! recalculate the intensities of only the pixels involving the moved atom

        old_int=0.0
        old_int_sq=0.0

        sum_int = 0.0
        sum_int_sq = 0.0

        psum_int = 0.0
        psum_int_sq = 0.0

        rotations: do i=myid+1, nrot, numprocs
        write(*,*) 'rotation ', i
            do m=1, npix
                old_int(1:nk, m, i) = int_i(1:nk, m, i)
                old_int_sq(1:nk, m, i) = int_sq(1:nk, m, i)
            enddo

            ! Rotate that moved atom
            call rotate_model(rot(i,1), rot(i, 2), rot(i, 3), moved_atom, rot_atom, istat)

            if( (rot_atom%natoms == 0) .and. (mrot(i)%rot_i(atom)%nat == 0) ) goto 100   !JWH - 04/14/2009
            ! the rotated atom has left the model, so there is no structural
            ! change - skip the rest of the loop
            ! if moving the atom changes the number of times it appears in the
            ! rotated model,
            ! just give up: re-rotate the entire model, and recalculate all
            ! the intensities
            ! (might not be necessary to recalculate all the intensities.  have
            ! to think about that
            ! TODO

            if(mrot(i)%rot_i(atom)%nat /= rot_atom%natoms) then
                !write(*,*)"scratch", mrot(i)%rot_i(atom)%nat, rot_atom%natoms
                ! store the original index and position in old_index and old_pos
                if(mrot(i)%rot_i(atom)%nat /= 0) then
                    do j=1, mrot(i)%rot_i(atom)%nat
                        call add_index(old_index(i), mrot(i)%rot_i(atom)%ind(j))
                        call add_pos(old_pos(i), mrot(i)%xx(mrot(i)%rot_i(atom)%ind(j)), &
                        mrot(i)%yy(mrot(i)%rot_i(atom)%ind(j)), mrot(i)%zz(mrot(i)%rot_i(atom)%ind(j)), istat)
                    enddo
                endif

                !write (*,*) 'Rerotating the model from scratch.'

                old_mrot_roti_nat = mrot(i)%rot_i(atom)%nat
                call destroy_model(mrot(i))
                call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m_in, mrot(i), istat)

                do m=1, npix
                    no_int_recal = .true.
                    if(rot_atom%natoms /= 0) then
                        do n=1, rot_atom%natoms
                            temp1 = rot_atom%xx(n) - pix(m,1)
                            temp2 = rot_atom%yy(n) - pix(m,2)
                            temp1 = temp1 - mrot(i)%lx*anint(temp1/mrot(i)%lx)
                            temp2 = temp2 - mrot(i)%ly*anint(temp2/mrot(i)%ly)
                            rr_x_new = ABS(temp1)
                            rr_y_new = ABS(temp2)
                            rot_dist_sq = (temp1)**2 + (temp2)**2

                            if(pixel_square) then
                                sqrt1_2_res = SQRT(0.5) * res !JASON added this line - taken from fem1.f90 in ~/Jinwoo/2011/060911_rmc_eam_gr_t0/ ???
                                if((rr_x_new .LE. sqrt1_2_res) .AND. (rr_y_new .LE.  sqrt1_2_res)) THEN
                                    call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e, istat, pixel_square)
                                    int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                                    no_int_recal = .false.
                                endif
                            else
                                sqrt1_2_res = res !JASON added this line - taken from fem1.f90 in ~/Jinwoo/2011/060911_rmc_eam_gr_t0/ ???
                                if (rot_dist_sq <= res2) then
                                    call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e, istat, pixel_square)
                                    int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                                    no_int_recal = .false.
                                endif
                            endif
                        enddo
                    endif

                    if(no_int_recal)then
                        if(old_mrot_roti_nat /= 0) then
                            do n=1, old_mrot_roti_nat
                                temp1 = old_pos(i)%pos(n,1) - pix(m,1)
                                temp2 = old_pos(i)%pos(n,2) - pix(m,2)
                                temp1 = temp1 - mrot(i)%lx*anint(temp1/mrot(i)%lx)
                                temp2 = temp2 - mrot(i)%ly*anint(temp2/mrot(i)%ly)
                                rr_x_old = ABS(temp1)
                                rr_y_old = ABS(temp2)
                                orig_dist_sq = (temp1)**2 + (temp2)**2

                                if(pixel_square) then
                                    sqrt1_2_res = SQRT(0.5) * res !JASON added this line - taken from fem1.f90 in ~/Jinwoo/2011/060911_rmc_eam_gr_t0/ ???
                                    if((rr_x_old .LE. sqrt1_2_res) .AND.  (rr_y_old .LE.  sqrt1_2_res)) then
                                        call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e, istat, pixel_square)
                                        int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                                    endif
                                else
                                    if (orig_dist_sq <= res2) then
                                        call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e, istat, pixel_square)
                                        int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                                    endif
                                endif
                            enddo
                        endif
                    endif
                enddo

                old_index(i)%nat = -1

            else
                !write(*,*)"reg", mrot(i)%rot_i(atom)%nat, rot_atom%natoms

                ! replace the original, rotated atom with its new position(s).
                ! if it appears more than once
                ! in the new model, the rotated model atoms should appear in the
                ! same order as they do in the full-size, rotated model
                !write (*,*) 'replacing rotated atoms with new positions.'

                do j=1,mrot(i)%rot_i(atom)%nat
                    ! store the original index and position in old_index and
                    ! old_pos
                    call add_index(old_index(i), mrot(i)%rot_i(atom)%ind(j))
                    call add_pos(old_pos(i), mrot(i)%xx(mrot(i)%rot_i(atom)%ind(j)), &
                    mrot(i)%yy(mrot(i)%rot_i(atom)%ind(j)), mrot(i)%zz(mrot(i)%rot_i(atom)%ind(j)), istat)

                    ! change the model atom position to the new, rotated
                    ! position
                    mrot(i)%xx(mrot(i)%rot_i(atom)%ind(j)) = rot_atom%xx(j)
                    mrot(i)%yy(mrot(i)%rot_i(atom)%ind(j)) = rot_atom%yy(j)
                    mrot(i)%zz(mrot(i)%rot_i(atom)%ind(j)) = rot_atom%zz(j)
                enddo

                ! now check if either the original position of the moved atom or
                ! the rotated position of the moved
                ! atom are inside each pixel.  If so, that pixel intensity must
                ! be recalculated.
                ! replace this with one loop over the rotated atom positions and
                ! another over the original atom
                ! positions, and if either is true, calculate the intensities.
                ! that way it works for either case above.

                do m=1, npix
                    do n=1, rot_atom%natoms
                        !This didn't have PBC. Now it's added below - 030409 -
                        !JWH
                        !rot_dist_sq = (mrot(i)%xx(n) - pix(m,1))**2 +
                        !(mrot(i)%yy(n) -
                        !pix(m,2))**2

                        temp1 = rot_atom%xx(n) - pix(m,1)
                        temp2 = rot_atom%yy(n) - pix(m,2)
                        temp1 = temp1 - mrot(i)%lx*anint(temp1/mrot(i)%lx)
                        temp2 = temp2 - mrot(i)%ly*anint(temp2/mrot(i)%ly)
                        rr_x_new = ABS(temp1)
                        rr_y_new = ABS(temp2)

                        rot_dist_sq = (temp1)**2 + (temp2)**2

                        temp1 = old_pos(i)%pos(n,1) - pix(m,1)
                        temp2 = old_pos(i)%pos(n,2) - pix(m,2)
                        temp1 = temp1 - mrot(i)%lx*anint(temp1/mrot(i)%lx)
                        temp2 = temp2 - mrot(i)%ly*anint(temp2/mrot(i)%ly)

                        rr_x_old = ABS(temp1)
                        rr_y_old = ABS(temp2)

                        orig_dist_sq = (temp1)**2 + (temp2)**2

                        if(pixel_square) then
                            if(((rr_x_old .LE. res) .AND. (rr_y_old .LE. res)) .OR. ((rr_x_new .LE. res) .AND. (rr_y_new .LE. res))) THEN
                                call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e,istat,pixel_square)
                                int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                            endif
                        else
                            if ( (rot_dist_sq <= res2) .OR. (orig_dist_sq <= res2) ) then
                                call intensity(mrot(i), res, pix(m, 1), pix(m, 2), k, int_i(1:nk, m, i), scatfact_e,istat,pixel_square)
                                int_sq(1:nk, m, i) = int_i(1:nk, m,i)**2
                            endif
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

            do n=1, size(rot_atom%rot_i,1)
                if(associated(rot_atom%rot_i(n)%ind))then
                    deallocate(rot_atom%rot_i(n)%ind)
                endif
            enddo

            deallocate(rot_atom%xx, rot_atom%yy, rot_atom%zz, rot_atom%znum, rot_atom%rot_i, rot_atom%znum_r, stat=istat)

        enddo rotations

        call mpi_reduce (psum_int, sum_int, size(k), mpi_real, mpi_sum, 0, comm, mpierr)
        call mpi_reduce (psum_int_sq, sum_int_sq, size(k), mpi_real, mpi_sum, 0, comm, mpierr)

        ! recalculate the variance
        if(myid.eq.0)then
            do i=1, nk
                Vk(i) = (sum_int_sq(i)/(npix*nrot))/((sum_int(i)/(npix*nrot))**2)-1.0
                Vk(i) = Vk(i) - v_background(i)   !background subtraction 052210 JWH
            end do
        endif

        deallocate(moved_atom%xx, moved_atom%yy, moved_atom%zz, moved_atom%znum, moved_atom%atom_type, moved_atom%znum_r, moved_atom%composition, stat=istat)
        deallocate(psum_int, psum_int_sq, sum_int, sum_int_sq)
    end subroutine fem_update

    subroutine fem_accept_move(comm)
    ! accept the move.  don't need to change any of the atom positions in any of the rotated
    ! models, but we do need to clear the old_index and old_pos arrays for reuse.
        use mpi
        integer :: comm
        call fem_reset_old(comm)
    end subroutine fem_accept_move

    subroutine fem_reset_old(comm)
    use mpi
    integer :: i, comm
    do i=myid+1, nrot, numprocs
        if(old_index(i)%nat > 0) deallocate(old_index(i)%ind)
            nullify(old_index(i)%ind)
        old_index(i)%nat = 0
        if(old_pos(i)%nat > 0) deallocate(old_pos(i)%pos)
            nullify(old_pos(i)%pos)
        old_pos(i)%nat = 0
    enddo
    end subroutine fem_reset_old

    subroutine fem_reject_move(m, comm)
    ! reject the move.  replace all the moved atoms with their original positions. 
    ! requires that the model argument have the original (unmoved) positions in it.
        use mpi
        type(model), intent(inout) :: m
        integer :: i, j, istat
        integer :: comm

        do i=myid+1, nrot, numprocs
            ! the rotated atom wasn't in the model, so this model doesn't need
            ! to be changed
            if(old_index(i)%nat == 0) cycle
            ! the move changed the number of atoms in the model, so the model
            ! must be re-rotated
            ! from scratch
            if(old_index(i)%nat == -1) then
                call destroy_model(mrot(i))
                call rotate_model(rot(i,1), rot(i,2), rot(i,3), m, mrot(i), istat)
                cycle
            endif

            ! otherwise, copy the old positions back into the model at the
            ! correct indices
            do j=1,old_index(i)%nat
                mrot(i)%xx(old_index(i)%ind(j)) = old_pos(i)%pos(j,1)
                mrot(i)%yy(old_index(i)%ind(j)) = old_pos(i)%pos(j,2)
                mrot(i)%zz(old_index(i)%ind(j)) = old_pos(i)%pos(j,3)
            enddo

            !The saved intensity values must return to their old values - JWH
            !03/05/09
            do j=1, npix
                int_i(1:nk, j, i) = old_int(1:nk, j, i)
                int_sq(1:nk, j, i) = old_int_sq(1:nk, j, i)
            enddo
        enddo
        call fem_reset_old(comm)
    end subroutine fem_reject_move

    subroutine add_pos(p, xx, yy, zz, istat)
        type(pos_list), intent(inout) :: p
        real, intent(in) :: xx, yy,  zz
        integer, intent(out) :: istat
        real, dimension(:,:), allocatable :: scratch
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


end module fem_mod

