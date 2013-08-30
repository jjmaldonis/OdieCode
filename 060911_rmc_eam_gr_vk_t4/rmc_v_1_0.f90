
! Reverse Monte Carlo structural simulation
!
! Performs RMC refinement against reduced density
! function data from electron, x-ray, or neutron
! diffraction data and against fluctuation electron
! microscopy V(k) data.
!
! Voyles research group: Jinwoo Hwang, Feng, Yi, and
! Paul Voyles, begun 12/16/08
! 
! merged with version of fy, 12/20/08, pmv
! more program flow added 12/20/08 pmv
! add new variable: electron scattering coefficients for calculating
! which are a_e, b_e, c_e and d_e
! add new module scattering_factors, add USE without :: following them for
! visual fortran debug by fy on 12/23
! need to not define gr_sim etc. arrays
! Re-written by Jinwoo Hwang - 03/20/2009
! Hard sphere cutoff applied.
! Debugged and tested by Jinwoo Hwang - 03/25/2009
! Vk allocated before initial Vk - JWH 04/02/2009
! hutch_move_atom added in reject part - JWH 04/15/2009.
! gr_e_err is added - jwh 04/25/2009
! See commit history for changes by jjm on Github: (check all branches)
!   https://github.com/refreshx2/OdieCode/tree/master/060911_rmc_eam_gr_vk_t4

program rmc

    use omp_lib
    use rmc_global
    use readinputs
    use model_mod
    use gr_mod
    use fem_mod
    use rmc_functions
    use eam_mod
    implicit none
    include 'mpif.h'
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: param_filename  
    character (len=256) :: outbase
    character (len=256) :: jobID, c
    character (len=512) :: comment
    character (len=256) :: time_elapsed, vki_fn, vku_fn, vkf_fn, output_model_fn, energy_fn, final_model_fn, gri_fn
    logical, dimension(4) :: used_data_sets
    real :: temperature
    real :: max_move
    real :: Q, res
    real :: pixel_distance
    real, dimension(4) :: weights
    real, pointer, dimension(:) :: gr_e,r_e,gr_e_err
    real, pointer, dimension(:) :: gr_n,r_n
    real, pointer, dimension(:) :: gr_x,r_x
    real, pointer, dimension(:) :: vk, vk_exp, k, vk_exp_err, v_background, vk_as !vk_as is for autoslice/multislice (whatever it's called)
    real, pointer, dimension(:,:) :: cutoff_r 
    real, pointer, dimension(:,:) :: scatfact_e
    real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: chi2_old, chi2_new, scale_fac, scale_fac_initial, del_chi, beta, chi2_gr, chi2_vk
    real :: rmin_e, rmax_e
    real :: rmin_n, rmax_n
    real :: rmin_x, rmax_x
    real :: R
    integer :: i, j
    integer :: w
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: fem_algorithm
    integer :: istat, status2, length
    integer :: total_steps
    integer :: iseed2
    real :: randnum
    real :: te1, te2
    logical :: square_pixel, accepted, use_rmc, use_multislice
    integer :: ipvd, nthr
    doubleprecision :: t0, t1, t2 !timers

    if(myid.eq.0)then
        write(*,*)
        write(*,*) "This is the dev version of rmc!"
        write(*,*)
    endif

    !call mpi_init(mpierr)
    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr) !http://www.open-mpi.org/doc/v1.5/man3/MPI_Init_thread.3.php
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)
    
    nthr = omp_get_max_threads()
    if(myid.eq.0) then
        write(*,*) "OMP found a max number of threads of", nthr
        write(*,*)
    endif
    
    t0 = omp_get_wtime()

    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        write (*,*) 'No jobID given.'
        write(*,*)
        jobID = ''
    end if

    ! Set input filenames.
    model_filename = 'model_040511c_t2_final.xyz'
    !model_filename = 'double_model.xyz'
    !model_filename = 'al50k_paul.xyz'
    param_filename = 'param_file.in'

    ! Set output filenames.
    outbase = ""
    write(time_elapsed, "(A12)") "time_elapsed"
    time_elapsed = trim(trim(time_elapsed)//jobID)//".txt"
    write(vki_fn, "(A10)") "vk_initial"
    vki_fn = trim(trim(vki_fn)//jobID)//".txt"
    write(vku_fn, "(A9)") "vk_update"
    vku_fn = trim(trim(vku_fn)//jobID)//".txt"
    write(vkf_fn, "(A8)") "vk_final"
    vkf_fn = trim(trim(vkf_fn)//jobID)//".txt"
    write(output_model_fn, "(A12)") "model_update"
    output_model_fn = trim(trim(output_model_fn)//jobID)//".txt"
    write(final_model_fn, "(A11)") "model_final"
    final_model_fn = trim(trim(final_model_fn)//jobID)//".txt"
    write(energy_fn, "(A6)") "energy"
    energy_fn = trim(trim(energy_fn)//jobID)//".txt"
    write(gri_fn, "(A10)") "gr_initial"
    gri_fn = trim(trim(gri_fn)//jobID)//".txt"

    ! Read input model
    call read_model(model_filename, comment, m, istat)
    call check_model(m, istat)
    call recenter_model(0.0, 0.0, 0.0, m)

    ! Read input parameters
    allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
    call read_inputs(param_filename,temperature, max_move, cutoff_r, &
        used_data_sets, weights, gr_e, r_e, gr_e_err, gr_n, r_n, gr_x, &
        r_x, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, &
        scale_fac_initial, Q, fem_algorithm, pixel_distance, total_steps, &
        rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

    scale_fac = scale_fac_initial
    res = 0.61/Q
    nk = size(k)
    beta=1./((8.6171e-05)*temperature)

    iseed2 = 104756

    square_pixel = .TRUE. ! RMC uses square pixels, not round.
    use_rmc = .TRUE.
    use_multislice = .FALSE.

    call read_eam(m)
    call eam_initial(m,te1)

    call scatt_power(m,used_data_sets,istat)
    call gr_initialize(m,r_e,gr_e,r_n,gr_n,r_x,gr_x,used_data_sets,istat)
    call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat,  square_pixel)
    allocate(vk(size(vk_exp)))
    allocate(vk_as(size(vk_exp)))
    if(myid.eq.0) call print_sampled_map(m, res, square_pixel)
    if(myid.eq.0) then
        if(pa%npix /= 1) then
            if(numprocs > 3*nrot) write(*,*) "WARNING: You are using too many cores!"
        else
            if(numprocs > nrot) write(*,*) "WARNING: You are using too many cores!"
        endif
    endif
    ! Fem updates vk based on the intensity calculations and v_background.
    call fem(m, res, k, vk, vk_as, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)

    t1 = omp_get_wtime()
    write(*,*) "Femsim took", t1-t0, "seconds on processor", myid

    if(myid.eq.0)then
        ! Write initial gr
        !open(unit=51,file=trim(gri_fn),form='formatted',status='unknown')
        !    do i=1, mbin_x
        !        R = del_r_x*(i)-del_r_x
        !        write(51,*)R, gr_x_sim_cur(i)
        !    enddo
        !close(51)
        ! Write initial vk 
        open(unit=52,file=trim(vki_fn),form='formatted',status='unknown')
            write(*,*) "Writing V(k)."
            do i=1, nk
                write (*,*) k(i),vk(i)
                write(52,*) k(i),vk(i)
            enddo
        close(52)
    endif

    if(use_rmc) then ! End here if we only want femsim. Set the variable above.

        ! Initial chi2
        chi2_old = chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, vk_exp, vk_exp_err, &
            gr_e_sim_cur, gr_n_sim_cur, gr_x_sim_cur, vk, scale_fac,&
            rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)
        chi2_old = chi2_old + te1 

        e2 = e1

        ! RMC step begins
        t0 = omp_get_wtime()
        i=1
        if(myid.eq.0)then
            write(*,*)
            write(*,*) "Initialization complete. Starting Monte Carlo."
            write(*,*) "Initial Conditions:"
            write(*,*) "   Step =      ", i
            write(*,*) "   Energy =     ", te1
            write(*,*) "   Temperature =", temperature
            write(*,*)
        endif

        ! Reset time_elapsed and energy_function
        open(35,file=trim(time_elapsed),form='formatted',status='unknown')
            write(35,*) numprocs, "processors being used. OpenMP + MPI."
            t1 = omp_get_wtime()
            write (35,*) i, t1-t0
        close(35)
        open(34,file=trim(energy_fn),form='formatted',status='unknown')
            write(34,*) i, te1
        close(34)


        ! The loops stops if the temperature goes below a certain temp (30.0).
        do while (i > 0)
            i=i+1
            t2 = omp_get_wtime()

            !if( i > 100) then
            !    write(*,*) "STOPPING MC AFTER 100 STEPS"
            !    call mpi_finalize(mpierr)
            !    stop ! Stop after 100 steps for timing runs.
            !endif

            call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            ! check_curoffs returns false if the new atom placement is too close to
            ! another atom. Returns true if the move is okay. (hard shere cutoff)
            do while( .not. check_cutoffs(m,cutoff_r,w) )
                ! Check_cutoffs returned false so reset positions and try again.
                m%xx%ind(w) = xx_cur
                m%yy%ind(w) = yy_cur
                m%zz%ind(w) = zz_cur
                call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            end do

            ! Update hutches, data for chi2, and chi2/del_chi
            call hutch_move_atom(m,w,xx_new, yy_new, zz_new)
            call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
            ! Use multislice every 10k steps if specified.
            if(use_multislice .and. mod(i,10000) .eq. 0) then
                call fem_update(m, w, res, k, vk, vk_as, v_background, scatfact_e, mpi_comm_world, istat, square_pixel, .true.)
                call update_scale_factor(scale_fac, scale_fac_initial, vk, vk_as)
            else
                call fem_update(m, w, res, k, vk, vk_as, v_background, scatfact_e, mpi_comm_world, istat, square_pixel, .false.)
            endif
            !write(*,*) "Finished updating eam, gr, and fem data."
            
            chi2_new = chi_square(used_data_sets,weights,gr_e, gr_e_err, &
                gr_n, gr_x, vk_exp, vk_exp_err, gr_e_sim_new, gr_n_sim_new, &
                gr_x_sim_new, vk, scale_fac, rmin_e, rmax_e, rmin_n, rmax_n, &
                rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)

            chi2_new = chi2_new + te2
            del_chi = chi2_new - chi2_old

            call mpi_bcast(del_chi, 1, mpi_real, 0, mpi_comm_world, mpierr)
               
            randnum = ran2(iseed2)
            ! Test if the move should be accepted or rejected based on del_chi
            !if(del_chi <0.0)then
            if(.true.)then !For timing purposes, always accept the move.
                ! Accept the move
                e1 = e2
                call fem_accept_move(mpi_comm_world)
                chi2_old = chi2_new
                accepted = .true.
                write(*,*) "MC move accepted outright."
            else
                ! Based on the random number above, even if del_chi is negative, decide
                ! whether to move or not (statistically).
                if(log(1.-randnum)<-del_chi*beta)then
                    ! Accept move
                    e1 = e2
                    call fem_accept_move(mpi_comm_world)
                    chi2_old = chi2_new
                    accepted = .true.
                    if(myid.eq.0)then
                        write(*,*) "MC move accepted due to probability. del_chi*beta = ", del_chi*beta
                    endif
                else
                    ! Reject move
                    e2 = e1
                    call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                    call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
                    call fem_reject_move(m, mpi_comm_world)
                    accepted = .false.
                    write(*,*) "MC move rejected."
                endif
            endif

            ! Every 150,000 steps lower the temp, max_move, and reset beta.
            if(mod(i,150000)==0)then
                temperature = temperature * sqrt(0.7)
                if(myid.eq.0)then
                    write(*,*) "Lowering temp to", temperature, "at step", i
                endif
                max_move = max_move * sqrt(0.94)
                beta=1./((8.6171e-05)*temperature)
                if(myid.eq.0)then
                    if(temperature.lt.30.0)then
                        if(myid.eq.0)then
                            write(*,*) "STOPPING MC DUE TO TEMP"
                        endif
                        stop
                    endif
                endif
            endif

            ! Periodically save data.
            if(mod(i,1)==0)then
                if(myid.eq.0)then
                    ! Write to vk_update
                    open(32,file=trim(vku_fn),form='formatted',status='unknown')
                        do j=1, nk
                            write(32,*)k(j),vk(j)
                        enddo
                    close(32)
                    ! Write to model_update
                    !open(33,file=trim(output_model_fn),form='formatted',status='unknown')
                        !write(33,*)"updated model"
                        !write(33,*)m%lx,m%ly,m%lz
                        !do j=1,m%natoms
                        !    write(33,*)m%znum%ind(j), m%xx%ind(j), m%yy%ind(j), m%zz%ind(j)
                        !enddo
                        !write(33,*)"-1"
                    !close(33)
                    open(34,file=trim(energy_fn),form='formatted',status='unknown',access='append')
                        if(accepted) then
                            ! Write to energy_function
                            write(*,*) i, te2
                            write(34,*) i, te2
                        endif
                    close(34)
                    ! Write to time_elapsed
                    open(35,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
                        t1 = omp_get_wtime()
                        write (*,*) "Step, time elapsed, temp:", i, t1-t0, temperature
                        write (35,*) i, t1-t0
                    close(35)
                    write(*,*) "Time per step = ", (t1-t0)/(i)
                    write(*,*) "This step's time = ", t1-t2
                    ! Time per step * num steps before decreasing temp * num
                    ! drops in temp necessary to get to temp=30.
                    write(*,*) "Approximate time remaining in seconds:", (t1-t0)/i * (150000 * log(30/temperature)/log(sqrt(0.7)) - i)! / nthr
                    !write(*,*) "Approximate time remaining in seconds (for 100 steps):", (t1-t0)/i * (100-i)! / nthr
                endif
            endif
        enddo

        write(*,*) "Monte Carlo Finished!"

        ! The rmc loop finished. Write final data.
        if(myid.eq.0)then
            ! Write final vk
            open(unit=54,file=trim(vkf_fn),form='formatted',status='unknown')
            do i=1, nk
                write(54,*)k(i),vk(i)
            enddo
            close(54)
            ! Write final model
            open(unit=55,file=trim(final_model_fn),form='formatted',status='unknown')
            write(55,*)"final model"
            write(55,*)m%lx,m%ly,m%lz
            do i=1,m%natoms
                write(55,*)m%znum%ind(i), m%xx%ind(i), m%yy%ind(i), m%zz%ind(i)
            enddo
            write(55,*)"-1"; close(55)
            ! Write final energy.
            open(56,file=trim(energy_fn),form='formatted', status='unknown',access='append')
            write(56,*) i, te2
            close(56)
            ! Write final time spent.
            open(57,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
            t1 = omp_get_wtime()
            write (57,*) i, t1-t0
            write(57,*) "Finshed.", numprocs, "processors."
            close(57)
        endif
    endif ! Use RMC
    call mpi_finalize(mpierr)

end program rmc
