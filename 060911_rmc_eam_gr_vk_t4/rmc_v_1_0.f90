
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


program rmc

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

character (len=256) model_filename
character (len=256):: param_filename  
character (len=512) comment

logical, dimension(4) :: used_data_sets
logical :: not_too_short
real :: temperature
real :: max_move
real :: Q, res
real :: pixel_distance
real, dimension(4) :: weights
real, pointer, dimension(:) :: gr_e,r_e,gr_e_err
real, pointer, dimension(:) :: gr_n,r_n
real, pointer, dimension(:) :: gr_x,r_x
real, pointer, dimension(:) :: vk, vk_exp, k, vk_exp_err, v_background
real, pointer, dimension(:,:) :: cutoff_r 
real, pointer, dimension(:,:) :: scatfact_e
real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
real :: chi2_old, chi2_new, scale_fac, del_chi, beta, chi2_gr, chi2_vk
real :: rmin_e, rmax_e
real :: rmin_n, rmax_n
real :: rmin_x, rmax_x
real :: R
integer :: i, j
integer :: w
integer :: nk
integer :: ntheta, nphi, npsi
integer :: fem_algorithm
integer :: istat, status2
integer :: total_steps
integer :: iseed2
real :: randnum
real :: te1, te2
logical :: square_pixel, use_femsim
doubleprecision :: t0, t1



call mpi_init(mpierr)
call mpi_comm_rank(mpi_comm_world, myid, mpierr)
call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

model_filename = 'model_1.xyz'
param_filename = 'param_file.in'


!read input model
call read_model(model_filename, comment, m, istat)
call check_model(m, istat)


if(istat.eq.0) then
!	write (*,*) 'Input model checks out OK.'
endif
call recenter_model(0.0, 0.0, 0.0, m)
!write (*,*) 'model recentered'

!read input parameters
allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
call read_inputs(param_filename,temperature, max_move, cutoff_r, used_data_sets, weights, gr_e, r_e, gr_e_err, gr_n, r_n, &
gr_x, r_x, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, scale_fac, Q, fem_algorithm, pixel_distance, total_steps, &
rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

!write(*,*) 'parameter file read'

!if(myid.eq.0)then
!do i=1, size(k)
!	write(*,*)k(i), vk_exp(i), vk_exp_err(i), V_background(i)
!enddo
!endif

res = 0.61/Q
nk = size(k)
beta=1./((8.6171e-05)*temperature)
!beta=1./(temperature)

iseed2 = 104756

square_pixel = .TRUE.
use_femsim = .FALSE.


call read_eam(m)   
call eam_initial(m, te1)


write(*,*)"te1=", te1

!initialize and calculate initial gr
call scatt_power(m,used_data_sets,istat)
!write (*,*) 'Znum reduced, scatt pow initialized.'
call gr_initialize(m,r_e,gr_e,r_n,gr_n,r_x,gr_x,used_data_sets,istat)


call gr_no_hutch(m,used_data_sets)

!initialize and calculate initial vk
call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat,  square_pixel)


allocate(vk(size(vk_exp)))



call fem(m, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)


if(myid.eq.0)then
!write initial gr

open(unit=51,file="test_gr_initial.txt",form='formatted',status='unknown')
do i=1, mbin_x
	R = del_r_x*(i)-del_r_x
	write(51,*)R, gr_x_sim_cur(i)
enddo
close(51)


!write initial vk 
open(unit=52,file="test_vk_initial.txt",form='formatted',status='unknown')
do i=1, nk
   !write (*,*) k(i),vk(i)
   write(52,*)k(i),vk(i)

enddo
close(52)
endif

!initial chi2
chi2_old = chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, vk_exp, vk_exp_err, &
       gr_e_sim_cur, gr_n_sim_cur, gr_x_sim_cur, vk, scale_fac,&
	rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)


chi2_old = chi2_old + te1 



e2 = e1

if(myid.eq.0)then
   write(*,*)"initial"
   write(*,*)chi2_old, chi2_gr, chi2_vk, te1, temperature
   write(*,*)"i, chi2_gr, chi2_vk, te1, temperature"
endif

!write(*,*)"MC start"
t0 = mpi_wtime()
!RMC step begins


!stop



!DO i=1, 2500000

!*********update here*******************
i=1315708
!***************************************

DO WHILE (i >0)
	i=i+1

	100 continue	

	call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)	!random move

	not_too_short = check_cutoffs(m,cutoff_r,w)
	if(not_too_short) then	!hard sphere cutoff
		continue
	else
		m%xx(w) = xx_cur
		m%yy(w) = yy_cur
		m%zz(w) = zz_cur
		goto 100
	endif

	call hutch_move_atom(m,w,xx_new, yy_new, zz_new)	  !update hutches.

	call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
       
	call gr_hutch_mc(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new,used_data_sets,istat)
	
   
       call fem_update(m, w, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)
	
	chi2_new = chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, vk_exp, vk_exp_err,&
        gr_e_sim_new, gr_n_sim_new, gr_x_sim_new, vk, scale_fac,&
	rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)
	

	chi2_new = chi2_new + te2

	del_chi = chi2_new - chi2_old

	call mpi_bcast(del_chi, 1, mpi_real, 0, mpi_comm_world, mpierr)
       
	randnum = ran2(iseed2)

	if(del_chi <0.0)then
		e1 = e2
       
		call accept_gr(m, used_data_sets)
		call fem_accept_move(mpi_comm_world)
                !if(mod(i,100)==0)then     
                   if(myid.eq.0)then
                      write(*,*)i, chi2_gr, chi2_vk, te2, temperature
                   endif
                !endif
		chi2_old = chi2_new

	else

		if(log(1.-randnum)<-del_chi*beta)then
			e1 = e2

			call accept_gr(m, used_data_sets)
			call fem_accept_move(mpi_comm_world)
                        
                        !if(mod(i,100)==0)then
                           if(myid.eq.0)then
                              write(*,*)i, chi2_gr, chi2_vk, te2, temperature
                           endif
                        !endif   
			chi2_old = chi2_new
               
		else	

			e2 = e1

			call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                     call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)	  !update hutches.
			call reject_gr(m,used_data_sets)
			call fem_reject_move(m, mpi_comm_world)
	              	
!                        if(mod(i,100)==0)then
!                           if(myid.eq.0)then
!                              write(*,*)i, chi2_new, chi2_vk, te2, "reject"
!                           endif
!                        endif   
		endif
	
      endif

	if(mod(i,50000)==0)then

		temperature = temperature * sqrt(0.7)
		max_move = max_move * sqrt(0.94)
		
		beta=1./((8.6171e-05)*temperature)
		!if(myid.eq.0)then
		!	write(*,*)"temp=", temperature, "max_move=", max_move
		!endif

		if(myid.eq.0)then
			if(temperature.lt.30.0)then
				stop
			endif
		endif
	endif	



	!periodically saving data
	if(mod(i,1000)==0)then
           if(myid.eq.0)then
		open(31,file='test_gr_update.txt',form='formatted',status='unknown')
		open(32,file='test_vk_update.txt',form='formatted',status='unknown')		
		open(33,file='test_model_update.txt',form='formatted',status='unknown')		
	
              do j=1, mbin_x
                   R = del_r_x*(j)-del_r_x
                   write(31,*)R, gr_x_sim_new(j)
		enddo
		
		do j=1, nk
			write(32,*)k(j),vk(j)
		enddo
		
		write(33,*)"updated model"
		write(33,*)m%lx,m%ly,m%lz
		
              do j=1,m%natoms
                   write(33,*)m%znum(j), m%xx(j), m%yy(j), m%zz(j)
		enddo
		write(33,*)"-1"
	
		close(31)
		close(32)
		close(33)
	endif
     endif
ENDDO
!close(56)


if(myid.eq.0)then
	t1 = mpi_wtime()
	write(*,*)"time=", t1-t0, "sec"
	write(*,*)t1, t0



!	!write final gr
!	open(unit=53,file="test_gr_update_final.txt",form='formatted',status='unknown')
!	do i=1, mbin_e
!		R = del_r_e*(i)-del_r_e
!		write(53,*)R, gr_e_sim_new(i)
!	enddo
!	close(53)
	
!	!write final vk
!	open(unit=54,file="test_vk_update_final.txt",form='formatted',status='unknown')
!	do i=1, nk
!		write(54,*)k(i),vk(i)
!	enddo
!	close(54)
	
	!write final model
	open(unit=55,file="test_model_update_final.txt",form='formatted',status='unknown')
	write(55,*)"updated model"
	write(55,*)m%lx,m%ly,m%lz
	do i=1,m%natoms
		write(55,*)m%znum(i), m%xx(i), m%yy(i), m%zz(i)
	enddo
	write(55,*)"-1"
	close(55)
endif
call mpi_finalize(mpierr)

end program rmc