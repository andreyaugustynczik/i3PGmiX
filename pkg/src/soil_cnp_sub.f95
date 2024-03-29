module soil_cnp_subroutines


IMPLICIT NONE
public :: compute_leaching, soil_temperature_swat, soil_temperature_toy, moisture_rate_modifier, &
			temperature_rate_modifier, dz_soil_cs, soilcarb_age, decay_withacc, nitrogen_gas, &
			calc_psi_soil, compute_fertility

contains

	!Subroutines for leaching------- ----------------------------------------------------
	subroutine compute_leaching( soil_class, leaching_org, leaching_min, excessSW, ASW )

		integer, intent(in) :: soil_class
		real (kind=8), intent(in) :: excessSW, ASW
		real (kind=8), intent(inout) :: leaching_org, leaching_min
		real (kind=8) :: sand_prop, leach_org, leach_min 

		if (soil_class == 1) then
			sand_prop = 0.85d0
		else if (soil_class == 2) then
			sand_prop = 0.725d0
		else
			sand_prop = 0.5d0
		end if

		leaching_org = 	0.000183d0 * excessSW * (0.01d0 + 0.04d0 * sand_prop)
		leaching_min = 	excessSW / ASW
		
	end subroutine compute_leaching


	! Subroutine to calculate nitrogen gas losses----------------------------------------------------
	subroutine nitrogen_gas(n_sp, n_inorg_T, n_gas, nminl_gas, minl_n_tot, immob_n_tot )
		! Input parameters
		integer, intent(in) :: n_sp
		! In/Out parameters
		real(kind=8), dimension(n_sp), intent(inout) :: n_inorg_T, n_gas

		real(kind=8), intent(in) :: nminl_gas
		real(kind=8), dimension(n_sp), intent(in) :: minl_n_tot, immob_n_tot



		! Local variables
		integer :: i
		real(kind=8), dimension(n_sp)::  n_minl_tot_pr

		do i=1, n_sp
			n_minl_tot_pr(i) = MAX(minl_n_tot(i) - immob_n_tot(i), 0.0)  ! calculate total mineralization of N
			n_gas(i) = nminl_gas * n_minl_tot_pr(i)  ! calculate fraction of net mineralization of N that is lost as gas
			n_inorg_T(i) = n_inorg_T(i) - n_gas(i)  ! subtract N lost as gas from inorganic N pool
		end do

	end subroutine nitrogen_gas


	!Subroutines for soil temperature ----------------------------------------------------
	subroutine soil_temperature_swat(t_soil, t_avg, t_year_avg, asw, snow, lai, m, sc)
	!Calculates the soil temperature based on the SWAT model implementation,
	!assumes temperature at 50 cm depth 

		!Input
		integer, intent(in) :: m, sc
		real(kind=8), intent(in) :: t_avg, t_year_avg, asw, snow, lai
		real(kind=8), intent(out) :: t_soil
		real(kind=8) :: t_soil_upper, t_soil_lower
	
		!Parameters
		real(kind=8) :: lambda, df, zd, z, dd, dd_max, z_tot, phi, t_bare, theta, bdod, t_soil_surf !, bcv, bcv_snow, bcv_veg
		!lag parameter
		lambda = 0.0d0
		
		!soil depth and profile depth (mm)
		z_tot = 1000.d0
		
		!get bare ground temperature
		call soil_temperature_toy(t_soil_surf, t_avg, m, lai, snow)
		

		! soil water correction factor
		if (sc == 1) then
			theta = (50.d0 + asw)
			bdod = 1.5d0
		else if (sc == 2) then
			theta = (80.d0 + asw)
			bdod = 1.45d0
		else if (sc == 3) then
			theta = (150.d0 + asw)
			bdod = 1.42d0
		else 
			theta = (200.d0 + asw)
			bdod = 1.4d0
		end if
		
		!maximum damping depth (mm)
		dd_max = 1000.d0 + 2500.d0 * bdod / ( bdod + 686.d0 * ( exp( -5.63d0 * bdod ) ) )

		phi = theta / ( ( 0.356d0 - 0.144d0 * bdod ) * z_tot )
		
		!damping depth
		dd = dd_max * exp(  log( 500.d0 / dd_max ) * ( ( 1.d0 - phi ) / ( 1.d0 + phi ) )**2.d0 )
		
		!layer depth
		z = 200.d0
		
		!ratio of depth to damping depth
		zd = z / dd
		
		!depth factor
		df = zd / ( zd + exp(-0.867d0  - 2.078d0 * zd ))
		
		!soil temperature update
		t_soil_upper = lambda * t_soil + ( 1.d0 - lambda ) * ( df * ( t_year_avg - t_soil_surf ) + t_soil_surf )
		
		!layer depth
		z = 700.d0
		
		!ratio of depth to damping depth
		zd = z / dd
		
		!depth factor
		df = zd / ( zd + exp(-0.867d0  - 2.078d0 * zd ))
		
		!soil temperature update
		t_soil_lower = lambda * t_soil + ( 1.d0 - lambda ) * ( df * ( t_year_avg - t_soil_surf ) + t_soil_surf )
		
		t_soil = 0.7d0 * t_soil_upper + 0.3d0 * t_soil_lower

	end subroutine soil_temperature_swat
	
	subroutine soil_temperature_toy(t_soil, t_avg, m, lai, snow)

		!Input
		integer, intent(in) :: m
		real(kind=8), intent(in) :: t_avg, lai, snow
		real(kind=8), intent(out) :: t_soil
		real(kind=8) :: bcv_snow, corr, t_soil_out
		

		! Calculates soil surface temperature based on Toy et al. (1979)
		if (m == 12 .or. m == 1 .or. m == 2) then
			t_soil_out = t_avg*0.656d0 + 2.3967d0
		else if (m == 3 .or. m == 4 .or. m == 5) then
			t_soil_out = t_avg*1.052d0 + 1.0239d0
		else if (m == 6 .or. m == 7 .or. m == 8) then
			t_soil_out = t_avg*0.856d0 + 6.3928d0
		else
			t_soil_out = t_avg*1.023d0 + 1.2856d0
		end if
		
		!snow correction
		bcv_snow = snow / ( snow + exp( 6.055d0 * snow - 0.3002d0 * snow ) )
		
		corr = min( (1.d0 - bcv_snow) ,Exp( -0.06d0 * lai ) )
		
		!correct for LAI/snow
		if (t_soil_out > t_avg) then
			t_soil = t_soil_out * corr
		end if

	end subroutine soil_temperature_toy

	!Soil moisture rate modifier calculation 
	subroutine moisture_rate_modifier(ASW, asw_max, alpha_moist, peat_io)
	
		REAL(8), INTENT(IN) :: ASW, asw_max
		INTEGER, INTENT(IN) :: peat_io
		REAL(8), INTENT(OUT) :: alpha_moist
		REAL(8) :: sth_wilt, sth_opt, sth_optl, sth_resp_min, fsth_dry, fsth_wet 
		! NOTE:  'sat': 0.4456, 'crit': 0.2677,'wilt' : 0.119 is for typical Temperate forest
		! and need to be later mapped from the site-specific paramter
		sth_wilt = 0.119d0 / 0.445d0
		sth_opt  = 0.5d0 * (1.d0 + sth_wilt)
		sth_optl = MIN( 0.2677 / 0.445d0 , sth_opt )
		sth_resp_min = sth_wilt
		alpha_moist = 0.2d0

		! additional variables for peat
		fsth_dry = 0.0 
		fsth_wet = 0.2
		
		! No peatland inclusion
		if(peat_io .eq. int(1)) then
			
			if((ASW / asw_max) > sth_resp_min .AND. (ASW / asw_max) <= sth_opt) then
				alpha_moist = 0.2d0 + (0.8d0 * (((ASW / asw_max) - sth_resp_min) / (sth_opt - sth_resp_min)))
			else if ((ASW / asw_max) > sth_opt) then
				alpha_moist = 1.d0 - (0.8d0 * ((ASW / asw_max) - sth_opt))
			end if 

		!Peatland model 
		else 
			if ((ASW / asw_max) <= sth_optl) then
				alpha_moist = ((1 - fsth_dry) / sth_optl) * (ASW / asw_max)  + fsth_dry
			else if ((ASW / asw_max)  > sth_opt) then
				alpha_moist = ((1 - fsth_wet) / (1 - sth_opt)) * (1 - (ASW / asw_max) ) +  &
								fsth_wet
			else
				alpha_moist = 1.0
			end if 

		end if
		
	end subroutine moisture_rate_modifier


	!Temparature rate modifier calculation 
	subroutine temperature_rate_modifier(t_soil, tmp_ave, tmp_year, ASW, snow_water, lai, month, soil_class, alpha_temp, Q10_soil)
		REAL(8), INTENT(IN) :: tmp_ave, tmp_year, ASW, snow_water, Q10_soil
		INTEGER, INTENT(IN) :: month, soil_class
		REAL(8), INTENT(INOUT) :: t_soil
		REAL(8), INTENT(IN), DIMENSION(:) :: lai
		REAL(8), INTENT(OUT) :: alpha_temp
		
		call soil_temperature_swat(t_soil, tmp_ave, tmp_year, ASW, snow_water, SUM(lai(:)), month, soil_class)
		alpha_temp = (Q10_soil ** ((t_soil - 10.d0) / 10.d0))
	end subroutine temperature_rate_modifier

	! Calculate soil matric potential
	subroutine calc_psi_soil(asw, maxasw, sc, psi)
		
		real(kind=8), intent(in) :: asw, maxasw
		integer, intent(in) :: sc
			
		real(kind=8) :: psi, a, n
		! Parameters from Schaap, M. G., & Van Genuchten, M. T. (2006). 
			
		if (sc .eq. int(1)) then
			a = 0.02630268d0	
			n = 2.233572223d0	
		else if (sc .eq. int(2)) then
			a = 0.040738028d0	
			n = 1.191242008d0	
		else if (sc .eq. int(3)) then
			a = 0.012022644d0		
			n = 1.377209469d0	
		else
			a = 0.011220185d0	
			n = 1.300169578d0	
		endif
			
		! Parameters from Lu, N., Godt, J. W., & Wu, D. T. (2010). 
		! if (sc .eq. int(1)) then
			! a = 0.3d0
			! n = 3.d0
		! else if (sc .eq. int(2) .or. sc .eq. int(3)) then
			! a = 0.05d0
			! n = 2.5d0
		! else
			! a = 0.01d0
			! n = 1.8d0
		! endif
		
		psi = -(((1.d0/ (asw/maxasw))**(n/(n-1)) - 1.d0)**(1.d0/n) / a) / 1000.d0
		
		if (psi < -10.d0) then
			psi = -10.d0
		else if (psi > -0.033d0) then
			psi = -0.033d0
		end if
			
	end subroutine calc_psi_soil


	! Nutrient competition

	subroutine compute_fertility(Nav, Pav, Un, Up, root_dist, root_area, n_dist, fertility, n_sp)
		
		integer, intent(in) :: n_sp
		integer :: i
		real(kind=8), intent(in) :: Nav, Pav, n_dist
		real(kind=8), dimension(n_sp), intent(in) :: Un, Up, root_dist, root_area 
		real(kind=8), dimension(n_sp), intent(out) :: fertility
		real(kind=8), dimension(n_sp) :: Nav_sp, Pav_sp, Nav_ratio, Pav_ratio, fr_n_sp, fr_p_sp
		
		!Partition available nitrogen according to  root distribution, N distribution and specific root area/length (1 m of soil)
		!Root distribution is computed based on a power law (1-gamma^z) in accordance with 4C model implementation and nitrogen distribution
		!is computed based on an exponential function (1-exp(-kz)) based on the ICP soil survey shares of nitrogen at 30 cm and 100 cm soil depth (k=0.031).
		!Both the product of the function are integrated up to 1m soil depth to calculate the ratio of available nutrients for each species

		!Nav_ratio(:) = root_area(:) * ((root_dist(:)**100.d0 * n_dist**100.d0 + 100.d0 * (log(root_dist(:)) + log(n_dist)) - 1.d0) / &
		!			 (log(root_dist(:)) + log(n_dist)) + (1.d0 - root_dist(:)**100.d0)/log(root_dist(:)) + &
		!			 (1.d0 - n_dist**100.d0)/log(n_dist))

		Nav_ratio(:) = root_area(:) * ( 100.d0 + (-1.d0 + exp(-100.d0*n_dist))/n_dist + (1-exp(-100*n_dist)*root_dist(:)**100.d0)/ &
						(n_dist - log(root_dist(:))) + (1.d0 - root_dist(:)**100.d0)/log(root_dist(:)) )
		
		Nav_ratio(:) = Nav_ratio(:) / sum(Nav_ratio(:))
		
		!Define available N for each species
		Nav_sp(:) = Nav * Nav_ratio(:)
		
		!Check if P cycle is included
		if (Pav >= 0.d0) then
			!Define available P for each species, for now assumes same distribution as N
			Pav_sp(:) = Pav * Nav_ratio(:)
		end if
		
		!Get fertility rating
		do i=1, n_sp	
			!Nitrogen driven fertility ratio
			if (Un(i) == 0.d0) then
				fr_n_sp(i) = 1.d0
			else
				fr_n_sp(i) = Nav_sp(i) / Un(i)
			end if
			
			if (fr_n_sp(i) > 1.d0) then
				fr_n_sp(i) = 1.d0
			else if (fr_n_sp(i) < 0.d0) then
				fr_n_sp(i) = 0.d0
			else if (isnan(fr_n_sp(i))) then
				fr_n_sp(i) = 1.d0
			end if	
			
			if (Pav >= 0.d0) then
				!Phosphorus driven fertility ratio
				if (Up(i) == 0.d0) then
					fr_p_sp(i) = 1.d0
				else
					fr_p_sp(i) = Pav_sp(i)/Up(i)
				end if
				
				if (fr_p_sp(i) > 1.d0) then
					fr_p_sp(i) = 1.d0
				else if (fr_p_sp(i) < 0.d0) then
					fr_p_sp(i) = 0.d0
				else if (isnan(fr_p_sp(i))) then
					fr_p_sp(i) = 1.d0
				end if

				!Define fertility rating
				fertility(i) = min(fr_n_sp(i),fr_p_sp(i))
			else
				fertility(i) = fr_n_sp(i)
			end if
		end do
		
	end subroutine compute_fertility


	! soil layer discretization and mixing based on PFT types (to be modified with species level)
	subroutine dz_soil_cs(n_sp, dim_cslyr, zz_l, dzsoil, cs_dpm, cs_rpm, cs_bio, cs_hum, cs_dpm_dz, &
		cs_rpm_dz, cs_bio_dz, cs_hum_dz, t_soil, mix_term, mix_s, month, dz_extra)
		INTEGER, INTENT(IN) :: n_sp, dim_cslyr, month
		REAL(8), INTENT(IN), DIMENSION(n_sp) :: cs_dpm, cs_rpm, cs_bio, cs_hum, dz_extra
		INTEGER :: i, j, z 
		REAL(8), DIMENSION(dim_cslyr) :: zd, zs, acc_rates, dz_tmp
		REAL(8), INTENT(IN), DIMENSION(dim_cslyr) :: dzsoil
		REAL(8), INTENT(OUT), DIMENSION(dim_cslyr) :: zz_l
		REAL(8), INTENT(OUT), DIMENSION(n_sp,dim_cslyr) :: cs_dpm_dz, cs_rpm_dz, cs_bio_dz, cs_hum_dz
		REAL(8) :: z0, zsoil, mixparam, botlinear, botuniform, bioturb_mix, cryoturb_mix
		REAL(8), DIMENSION(n_sp,dim_cslyr,4) :: soil_bgc
		REAL(8), INTENT(INOUT), DIMENSION(n_sp,dim_cslyr,4) :: mix_term, mix_s
		REAL(8), INTENT(INOUT) :: t_soil
		
		botuniform   = 1.0
		! Bottom of the uniform part of mixing in the soil (m).
		! Initialisation of mixing 
		botlinear    = 3.0
		! Bottom of the linearly reducing part of mixing in the soil (m)
		bioturb_mix  = 0.0001*month
		! BIOTURBATION mixing (convert units to month) 
		cryoturb_mix  = 0.0005*month
		! CRYOTURBATION mixing (convert units to month) 

		mix_s(:,:,:)      = 0.0
		mix_term(:,:,:)   = 0.0
		acc_rates         = (/1.0,1.0,1.0,1.0 /)
		dz_tmp(:)         = 0.0

		do j = 1, dim_cslyr
			if (j == 1) then
				zd(j) = dzsoil(j) / 2
			else
				zd(j) = sum(dzsoil(1:j-1)) + dzsoil(j) / 2
			end if
		end do

		! THIS SHOULD CHANGE. FOR NOW IT IS THE SINGLE VALUE FOR TEMPERATE DECIDUOUS FOREST
		! IT SHOULD BE BASED ON SPECIES (PFT-estimated from Nakhavali et al., 2018)
		z0 = 0.725914
	
		! Now calculate the depth decay factor and then weighted 
		
		zz_l(:) = (exp(-zd(:)/z0)) * dzsoil(:) 
		
		zs = sum(zz_l(:))
		zz_l(:) = zz_l(:) / zs
	
		! Get layered SOC pools
		do i = 1, n_sp
			do j = 1, dim_cslyr ! soil layers
				cs_dpm_dz(i,j) = zz_l(j) * cs_dpm(i)
				cs_rpm_dz(i,j) = zz_l(j) * cs_rpm(i)
				cs_bio_dz(i,j) = zz_l(j) * cs_bio(i)
				cs_hum_dz(i,j) = zz_l(j) * cs_hum(i)

				soil_bgc(i,j,1) = cs_dpm_dz(i,j)
				soil_bgc(i,j,2) = cs_rpm_dz(i,j)
				soil_bgc(i,j,3) = cs_bio_dz(i,j)
				soil_bgc(i,j,4) = cs_hum_dz(i,j)
			end do
		end do

		! Mixing calculation 	
		! Calculate the mixing parameters
		do i =1, n_sp
				dz_tmp = dzsoil
				dz_tmp(dim_cslyr)=dz_tmp(dim_cslyr)+dz_extra(i)

				if(t_soil < 0.d0) then
					mixparam = cryoturb_mix
				else
					mixparam = bioturb_mix 
				end if 
				zsoil = 0.d0

				do j = 1,dim_cslyr-1
					zsoil = zsoil + dz_tmp(j)
					if (zsoil < botuniform*1.0) then
						mix_s(i,j,1) = mixparam * acc_rates(1)
						mix_s(i,j,2) = mixparam * acc_rates(2)
						mix_s(i,j,3) = mixparam * acc_rates(3)
						mix_s(i,j,4) = mixparam * acc_rates(4)
					else if (zsoil < botlinear*1.0) then
						mix_s(i,j,1) = mixparam  * acc_rates(1) *                               &
										(1.0 - (zsoil - botuniform) / (botlinear - botuniform))
						mix_s(i,j,2) = mixparam * acc_rates(2) *                                &
										(1.0 - (zsoil - botuniform) / (botlinear - botuniform))
						mix_s(i,j,3) = mixparam * acc_rates(3) *                                &
										(1.0 - (zsoil - botuniform) / (botlinear - botuniform))
						mix_s(i,j,4) = mixparam * acc_rates(4) *                                &
										(1.0 - (zsoil - botuniform) / (botlinear - botuniform))
					else
						mix_s(i,j,1) = 0.00000000001
						mix_s(i,j,2) = 0.00000000001
						mix_s(i,j,3) = 0.00000000001
						mix_s(i,j,4) = 0.00000000001
					end if
				end do
		end do
	

			do i = 1,n_sp
				do z = 1,4  !soil carbon pools
				  	mix_term(i,1,z) = mix_s(i,1,z) * ( (soil_bgc(i,2,z) / dz_tmp(2))        &
									- (soil_bgc(i,1,z) / dz_tmp(1)) )                       &
									/ ( 0.5 * (dz_tmp(2) + dz_tmp(1)) )
					do j = 2,dim_cslyr-1
						mix_term(i,j,z) = mix_s(i,j,z) * ( (soil_bgc(i,j+1,z)                   &
										/ dz_tmp(j+1)) - (soil_bgc(i,j,z) / dz_tmp(j)) )      &
										/ ( 0.5 * (dz_tmp(j+1) + dz_tmp(j)) )                 &
										- mix_s(i,j-1,z) *( (soil_bgc(i,j,z)                  &
										/  dz_tmp(j)) - (soil_bgc(i,j-1,z) / dz_tmp(j-1)) )   &
										/ ( 0.5 * (dz_tmp(j) + dz_tmp(j-1)) )
					end do
				  	mix_term(i,dim_cslyr,z) = mix_s(i,dim_cslyr-1,z) *                    &
											  ( (soil_bgc(i,dim_cslyr-1,z)                &
											  / dz_tmp(dim_cslyr-1))                      &
											  - (soil_bgc(i,dim_cslyr,z)                  &
											  / dz_tmp(dim_cslyr)) )                      &
											  / ( 0.5 * (dz_tmp(dim_cslyr)                &
											  + dz_tmp(dim_cslyr-1)) )
				  ! The mixing has to be done in terms of density but
				  ! converted back to original unit for the increment, hence no overall factor
				  ! of 1/dzsoil(j) as it cancels out.
				end do
			end do 

			do i = 1,n_sp
				do j = 1, dim_cslyr 
					cs_dpm_dz(i,j) = cs_dpm_dz(i,j) + mix_term(i,j,1)
					cs_rpm_dz(i,j) = cs_rpm_dz(i,j) + mix_term(i,j,2)
					cs_bio_dz(i,j) = cs_bio_dz(i,j) + mix_term(i,j,3)
					cs_hum_dz(i,j) = cs_hum_dz(i,j) + mix_term(i,j,4)
				end do
			end do 


	end subroutine dz_soil_cs
	

	subroutine soilcarb_age(n_sp, dim_cslyr, dz_extra, mix_s, resp_s, resp_frac, &
		soilcarb_age_pool, dcs, Yr_C, Yl_C, lit_frac, dpm_ratio, dzsoil, soilage_max, month)

			INTEGER, INTENT(IN) :: n_sp, dim_cslyr, month
			REAL(8), INTENT(IN), DIMENSION(n_sp) :: dz_extra , Yr_C, Yl_C, dpm_ratio
			REAL(8), INTENT(INOUT) :: soilage_max
			REAL(8), INTENT(IN), DIMENSION(n_sp,dim_cslyr,4) :: mix_s
			REAL(8), INTENT(IN), DIMENSION(n_sp,dim_cslyr,4) :: resp_s, dcs, resp_frac
			REAL(8), INTENT(INOUT), DIMENSION(n_sp,dim_cslyr,4) :: soilcarb_age_pool
			REAL(8), INTENT(IN), DIMENSION(dim_cslyr) :: dzsoil, lit_frac
			INTEGER :: i, j, z 
			REAL(8), DIMENSION(dim_cslyr,4) :: bal_b4_in, fracsXin, total_in
			REAL(8), DIMENSION(dim_cslyr,2) ::litterin
			REAL(8) :: invdz2top, invdz2bot,tmp_resp
			REAL(8), DIMENSION(dim_cslyr) :: dz_tmp	
			
			!----------------------------------------------------------------------
			! Update "old" soil  carbon tracer
			!----------------------------------------------------------------------
			fracsXin(:,:) = 0.d0
			bal_b4_in(:,:) = 0.d0
			total_in(:,:) = 0.d0

			do i = 1,n_sp
			
				! Update dz
				dz_tmp(:) = dzsoil(:)
				dz_tmp(dim_cslyr) = dz_tmp(dim_cslyr) + dz_extra(i)
				
				invdz2top = 1.0 / (dz_tmp(2) * (dz_tmp(1) + dz_tmp(2)))
				invdz2bot = 1.0 / ( dz_tmp(dim_cslyr-1) * (dz_tmp(dim_cslyr) + &
								dz_tmp(dim_cslyr-1)) )
				
				! Update age by timestep (units of month)
				do j = 1,dim_cslyr
					do z = 1,4  
					soilcarb_age_pool(i,j,z) = soilcarb_age_pool(i,j,z) + month !number of month 
					end do
				end do
				soilage_max = soilage_max + month
				
				! Define the litter inputs to pools 1 (dpm) and 2 (rpm)
				do j = 1,dim_cslyr
					litterin(j,1) = dpm_ratio(i)  * Yl_C(i) * lit_frac(j)
					litterin(j,2) = (1.0 - dpm_ratio(i)) * Yr_C(i) * lit_frac(j)
				end do
				
				!---------------------------------------------
				! Boundary values (top and bottom layers)
				!---------------------------------------------
				
				! How much of the "original" carbon is left in the pool after the
				! timestep
				! before any other carbon is input
				do z = 1,4
					bal_b4_in(1,z) =  dcs(i,1,z) - (  resp_s(i,1,z) + &
									(2.0 * mix_s(i,1,z) * &
									dcs(i,1,z)) / (dz_tmp(1) * (dz_tmp(1) + dz_tmp(2))) ) 
									
					bal_b4_in(dim_cslyr,z) = dcs(i,dim_cslyr,z) - &
											( resp_s(i,dim_cslyr,z) + &
											(2.0 * mix_s(i,dim_cslyr-1,z) * &
											dcs(i,dim_cslyr,z)) / &
											(dz_tmp(dim_cslyr) * (dz_tmp(dim_cslyr)    &
											+ dz_tmp(dim_cslyr-1))) ) 
				end do
				
				! Amount of other carbon that goes into the pool 
				! multiplied and not multiplied by its oldc fraction
				! Pools 1 and 2 no input from other pools
				do  z = 1,2
					fracsXin(1,z) = ( soilcarb_age_pool(i,2,z) * 2.0 *   &
									mix_s(i,1,z) * dcs(i,2,z) ) /    &
									( dz_tmp(2) * (dz_tmp(1) + dz_tmp(2))  )
					total_in(1,z) = ( 2.0 * mix_s(i,1,z) * dcs(i,2,z) ) /   &
									( dz_tmp(2) * (dz_tmp(1) + dz_tmp(2))  )
					fracsXin(dim_cslyr,z) = ( soilcarb_age_pool(i,dim_cslyr-1,z) * 2.0  &
											* mix_s(i,dim_cslyr-1,z) * &
											dcs(i,dim_cslyr-1,z) ) / &
											( dz_tmp(dim_cslyr-1) * (dz_tmp(dim_cslyr) +&
											dz_tmp(dim_cslyr-1))  )
					total_in(dim_cslyr,z) = ( 2.0 * mix_s(i,dim_cslyr-1,z) * &
											dcs(i,dim_cslyr-1,z) ) / &
											( dz_tmp(dim_cslyr-1) * (dz_tmp(dim_cslyr) +&
											dz_tmp(dim_cslyr-1))  )
				end do
				
				! Pools 3 and 4 need input from pools 1 and 2
				
				total_in(1,3) = ( (2.0 * mix_s(i,1,3) * dcs(i,2,3) * invdz2top) + &
								0.46 * SUM(resp_s(i,1,1:4) * resp_frac(i,1,:)) ) 
				total_in(1,4) = ( (2.0 * mix_s(i,1,4) * dcs(i,2,4) * invdz2top)  + &
								0.54 * SUM(resp_s(i,1,1:4) * resp_frac(i,1,:)) ) 
				
				tmp_resp = ( resp_frac(i,1,1) * soilcarb_age_pool(i,1,1) * resp_s(i,1,1) + &
							resp_frac(i,1,2) * soilcarb_age_pool(i,1,2) * resp_s(i,1,2) + &
							resp_frac(i,1,3) * soilcarb_age_pool(i,1,3) * resp_s(i,1,3) + &
							resp_frac(i,1,4) * soilcarb_age_pool(i,1,4) * resp_s(i,1,4) )
				fracsXin(1,3) = ( (soilcarb_age_pool(i,2,3) * 2.0 * mix_s(i,1,3) * &
								dcs(i,2,3) * invdz2top) + 0.46 * tmp_resp ) 
				fracsXin(1,4) = ( (soilcarb_age_pool(i,2,4) * 2.0 * mix_s(i,1,4) * &
								dcs(i,2,4) * invdz2top) + 0.54 * tmp_resp ) 
				
				total_in(dim_cslyr,3) = ( (2.0 * mix_s(i,dim_cslyr-1,3) * &
											dcs(i,dim_cslyr-1,3) * invdz2bot) + &
											0.46 * SUM(resp_frac(i,dim_cslyr,:) * resp_s(i,dim_cslyr,1:4)) ) 
											
				total_in(dim_cslyr,4) = ( (2.0 * mix_s(i,dim_cslyr-1,4) * &
											dcs(i,dim_cslyr-1,4) * invdz2bot) + &
											0.54 * SUM(resp_frac(i,dim_cslyr,:) * resp_s(i,dim_cslyr,1:4)) ) 
											
				
				tmp_resp = (resp_frac(i,dim_cslyr,1) * soilcarb_age_pool(i,dim_cslyr,1) * &
							resp_s(i,dim_cslyr,1) + resp_frac(i,dim_cslyr,2) * soilcarb_age_pool(i,dim_cslyr,2) *&
							resp_s(i,dim_cslyr,2) + resp_frac(i,dim_cslyr,3) * soilcarb_age_pool(i,dim_cslyr,3) *&
							resp_s(i,dim_cslyr,3) + resp_frac(i,dim_cslyr,4) * soilcarb_age_pool(i,dim_cslyr,4) *&
							resp_s(i,dim_cslyr,4))
				fracsXin(dim_cslyr,3) = ( (soilcarb_age_pool(i,dim_cslyr-1,3) * 2.0 &
											* mix_s(i,dim_cslyr-1,3) * &
											dcs(i,dim_cslyr-1,3) * invdz2bot) + &
											0.46 * tmp_resp )  
				fracsXin(dim_cslyr,4) = ( (soilcarb_age_pool(i,dim_cslyr-1,4) * 2.0 &
											* mix_s(i,dim_cslyr-1,4) * &
											dcs(i,dim_cslyr-1,4) * invdz2bot) + &
											0.54 * tmp_resp ) 
				
				!----------------------------------------------
				! Non-boundary values
				!----------------------------------------------
				! How much of the "original" carbon is left in the pool after the
				! timestep
				! before any other carbon is input
				do j = 2,dim_cslyr-1
					do z = 1,4
						bal_b4_in(j,z) = dcs(i,j,z) - ( resp_s(i,j,z) + 2.0 * dcs(i,j,z) *&
									(mix_s(i,j,z) / (dz_tmp(j) + dz_tmp(j+1)) + &
									mix_s(i,j-1,z) / (dz_tmp(j-1) + dz_tmp(j))) &
									/ dz_tmp(j) ) 
					end do
				
					! Amount of other carbon that goes into the pool 
					! multiplied and not multiplied by its oldc fraction
					! Pools 1 and 2 no input from other pools
					do z = 1,2
						total_in(j,z) = ( 2.0 * mix_s(i,j,z) * dcs(i,j+1,z) / (dz_tmp(j+1) *  &
									(dz_tmp(j) + dz_tmp(j+1))) &
									+ 2.0 * mix_s(i,j-1,z) * dcs(i,j-1,z) / (dz_tmp(j-1) *  &
									(dz_tmp(j) + dz_tmp(j-1))) ) 
									
						fracsXin(j,z) = ( soilcarb_age_pool(i,j+1,z) * 2.0 * mix_s(i,j,z) *   &
									dcs(i,j+1,z) / (dz_tmp(j+1) * (dz_tmp(j) + &
									dz_tmp(j+1))) + soilcarb_age_pool(i,j-1,z) * 2.0*    &
									mix_s(i,j-1,z) * dcs(i,j-1,z) / &
									(dz_tmp(j-1) * (dz_tmp(j) + dz_tmp(j-1))) ) 
									
					end do
				
				
					! Pools 3 and 4 need input from pools 1 and 2
					total_in(j,3) = ( 2.0 * mix_s(i,j,3) * dcs(i,j+1,3) / &
									(dz_tmp(j+1) * (dz_tmp(j) + dz_tmp(j+1))) + &
									2.0 * mix_s(i,j-1,3) * dcs(i,j-1,3) / &
									(dz_tmp(j-1) * (dz_tmp(j) + dz_tmp(j-1))) &
									+ 0.46 * SUM(resp_frac(i,j,:) * resp_s(i,j,1:4)) )  
					total_in(j,4) = ( 2.0 * mix_s(i,j,4) * dcs(i,j+1,4) / &
									(dz_tmp(j+1) * (dz_tmp(j) + dz_tmp(j+1))) + &
									2.0 * mix_s(i,j-1,4) * dcs(i,j-1,4) / &
									(dz_tmp(j-1) * (dz_tmp(j) + dz_tmp(j-1))) &
									+ 0.54 * SUM(resp_frac(i,j,:) * resp_s(i,j,1:4)) ) 
				
					tmp_resp = (resp_frac(i,j,1) * soilcarb_age_pool(i,j,1) * resp_s(i,j,1) + &
							resp_frac(i,j,2) * soilcarb_age_pool(i,j,2) * resp_s(i,j,2) + &
							resp_frac(i,j,3) * soilcarb_age_pool(i,j,3) * resp_s(i,j,3) + &
							resp_frac(i,j,4) * soilcarb_age_pool(i,j,4) * resp_s(i,j,4) )
					fracsXin(j,3) = ( soilcarb_age_pool(i,j+1,3) * 2.0 * mix_s(i,j,3) *     &
									dcs(i,j+1,3) / (dz_tmp(j+1) * (dz_tmp(j) + dz_tmp(j+1)))  &
									+ soilcarb_age_pool(i,j-1,3) * 2.0 * mix_s(i,j-1,3) *   &
									dcs(i,j-1,3) / (dz_tmp(j-1) * (dz_tmp(j) + dz_tmp(j-1)))  &
									+ 0.46 * tmp_resp ) 
					fracsXin(j,4) = ( soilcarb_age_pool(i,j+1,4) * 2.0 * mix_s(i,j,4) *     &
									dcs(i,j+1,4) / (dz_tmp(j+1) * (dz_tmp(j) + dz_tmp(j+1)))  &
									+ soilcarb_age_pool(i,j-1,4) * 2.0 * mix_s(i,j-1,4) *   &
									dcs(i,j-1,4) / (dz_tmp(j-1) * (dz_tmp(j) + dz_tmp(j-1)))  &
									+ 0.54 * tmp_resp ) 
				end do
				
				! Add litter to the carbon going in OR out 
				
				do j = 1,dim_cslyr
					do z = 1,2
						if (litterin(j,z) < 0.0) THEN
							bal_b4_in(j,z) = bal_b4_in(j,z) + litterin(j,z) 
						else
							total_in(j,z) = total_in(j,z) + litterin(j,z) 
						end if
					end do
				end do
				
				!-------------------------------------------------------------------
				! Calculate old carbon fractions
				!-------------------------------------------------------------------
				do j = 1,dim_cslyr
					do z = 1,4
					if (bal_b4_in(j,z) > 0.0) THEN
						if( (total_in(j,z) + bal_b4_in(j,z)) > 1.0e-8 ) THEN
							soilcarb_age_pool(i,j,z) = (soilcarb_age_pool(i,j,z) *        &
														bal_b4_in(j,z) + fracsXin(j,z)) /   &
														(total_in(j,z) + bal_b4_in(j,z))
						end if
					else if (total_in(j,z) > 1.0e-8) THEN
						soilcarb_age_pool(i,j,z) = fracsXin(j,z) / total_in(j,z)
					end if
					if (soilcarb_age_pool(i,j,z) < 0.0) THEN
						soilcarb_age_pool(i,j,z) = 0.d0   
					end if
					end do                                                    
				end do
			
			end do 
	
	end subroutine soilcarb_age



	subroutine decay_withacc(n_sp, dim_cslyr, dz_extra, pc, cs, soilcarb_age_pool, &
		ns_dpm_rpm, soilage_max, dzsoil, compress_dr, compress_bh)

		INTEGER, INTENT(IN) :: n_sp, dim_cslyr
		REAL(8), INTENT(IN), DIMENSION(n_sp,dim_cslyr,4) :: pc
		REAL(8), INTENT(INOUT), DIMENSION(n_sp) :: dz_extra 
		REAL(8), INTENT(INOUT) :: soilage_max
		REAL(8), INTENT(IN), DIMENSION(dim_cslyr) :: dzsoil
		REAL(8), INTENT(INOUT), DIMENSION(n_sp,dim_cslyr,4) :: soilcarb_age_pool, cs
		REAL(8), INTENT(INOUT), DIMENSION(n_sp,dim_cslyr,2) :: ns_dpm_rpm
		REAL(8), INTENT(OUT), DIMENSION(dim_cslyr) :: compress_dr, compress_bh
		INTEGER :: i, j, z 
		REAL(8), DIMENSION(dim_cslyr) :: dz_eff, z_eff, dzsoil_ext, zsoil_ext, &
		compress_eff_dr, compress_eff_bh
		REAL(8), DIMENSION(dim_cslyr,4) :: cdens_eff, cdens, cs_tmp, soilage_tmp
		REAL(8), DIMENSION(dim_cslyr,2) :: ndens_eff, ndens, ns_tmp
		REAL(8) :: totvol, dz_extra_extra, cs_min, rho_dpm_rpm, rho_bio_hum
		REAL(8), DIMENSION(4) ::soilage_avg
				
		dz_eff(:)     = 0.d0
		z_eff(:)      = 0.d0
		dzsoil_ext(:) = 0.d0
		zsoil_ext(:)  = 0.d0
		cdens_eff(:,:)= 0.d0
		cdens(:,:)    = 0.d0
		cs_tmp(:,:)   = 0.d0
		soilage_tmp(:,:) = 0.d0
		ndens_eff(:,:)= 0.d0
		ndens(:,:)    = 0.d0
		ns_tmp(:,:)   = 0.d0
		totvol = 0.d0
		dz_extra_extra = 0.d0
		soilage_avg(:) = 0.d0
		soilage_max = 100.d0
		rho_dpm_rpm = 19.6
		rho_bio_hum = 118.d0
		cs_min = 1.0e-6
		compress_eff_dr(:) = 0.d0
		compress_eff_bh(:) = 0.d0

		do i = 1,n_sp
			dzsoil_ext(:) = dzsoil

		! Calculate effective z's and carbon densities
			do j = 1,dim_cslyr
				if (j==1) then
					compress_eff_dr(j) = 1.0 - 0.03 * 9.8 * 0.0 * dzsoil_ext(j)
					compress_eff_bh(j) = 1.0 - 0.01 * 9.8 * 0.0 * dzsoil_ext(j)
				else
					compress_eff_dr(j) = 1.0 - 0.03 * 9.8 * ( SUM( dz_eff(1:(j-1)) ) * 0.0 * &
																				dzsoil_ext(j) )
					compress_eff_bh(j) = 1.0 - 0.01 * 9.8 * ( SUM( dz_eff(1:(j-1)) ) * 0.0 * &
																				dzsoil_ext(j) )
				end if
				compress_eff_dr(j) = MIN( MAX(compress_eff_dr(j), 0.62) , 1.0)
				compress_eff_bh(j) = MIN( MAX(compress_eff_bh(j), 0.87) , 1.0)

				dz_eff(j) = dzsoil(j) + SUM( pc(i,j,1:2)  ) / ( rho_dpm_rpm       & 
																	/ compress_eff_dr(j) )&
									+ SUM( pc(i,j,3:4) ) / ( rho_bio_hum       &
																	/ compress_eff_bh(j) )&
								+ (compress_eff_dr(j) - compress_dr(j)) * SUM(cs(i,j,1:2)) & 
								/ rho_dpm_rpm                                            &
								+ (compress_eff_bh(j) - compress_bh(j)) * SUM(cs(i,j,3:4)) &
								/ rho_bio_hum
			end do
			

			dz_eff(dim_cslyr)     = dz_eff(dim_cslyr)     + dz_extra(i)
			dzsoil_ext(dim_cslyr) = dzsoil_ext(dim_cslyr) + dz_extra(i)
		  
			dz_extra_extra = SUM(dz_eff(:)) - SUM(dzsoil_ext(:))
			dzsoil_ext(dim_cslyr) = dzsoil_ext(dim_cslyr) + dz_extra_extra
			dz_extra(i) = dz_extra(i) + dz_extra_extra
		  
			do j = 1,dim_cslyr
			  z_eff(j)     = SUM(dz_eff(1:j)) - 0.5 * dz_eff(j)
			  zsoil_ext(j) = SUM(dzsoil_ext(1:j)) - 0.5 * dzsoil_ext(j)
		  
			  do z = 1,4
				cdens_eff(j,z) = ( cs(i,j,z) + pc(i,j,z)  ) / dz_eff(j)
				if(z<3)ndens_eff(j,z) = ns_dpm_rpm(i,j,z) / dz_eff(j)
			  end do
			end do


		! Interpolate the carbon densities to the correct layer points
		! Note that at the moment we use cdens with C pools 
		! of 1:4 (DPM, RPM, BIO,HUM)
		do z = 1,4
			cdens(1,z) = cdens_eff(1,z) + (z_eff(1) - zsoil_ext(1)) /                 &
						(z_eff(1) - z_eff(2)) * (cdens_eff(2,z) - cdens_eff(1,z))
			cdens(dim_cslyr,z) = cdens_eff(dim_cslyr,z) +                         &
						(z_eff(dim_cslyr) - zsoil_ext(dim_cslyr)) /              &
						(z_eff(dim_cslyr) - z_eff(dim_cslyr-1)) *                &
						(cdens_eff(dim_cslyr-1,z) - cdens_eff(dim_cslyr,z))
			if (z<3) then
				ndens(1,z) = ndens_eff(1,z) + (z_eff(1) - zsoil_ext(1)) /               &
							(z_eff(1) - z_eff(2)) * (ndens_eff(2,i) - ndens_eff(1,i))
				ndens(dim_cslyr,z) = ndens_eff(dim_cslyr,z) +                       &
							(z_eff(dim_cslyr) - zsoil_ext(dim_cslyr)) /              &
							(z_eff(dim_cslyr) - z_eff(dim_cslyr-1)) *                &
							(ndens_eff(dim_cslyr-1,z) - ndens_eff(dim_cslyr,z))
			end if 
			soilage_tmp(1,z) = soilcarb_age_pool(i,1,z) + (z_eff(1) - zsoil_ext(1))&
					/ (z_eff(1) - z_eff(2)) * (soilcarb_age_pool(i,2,z) -       &
						soilcarb_age_pool(i,1,z))
			soilage_tmp(dim_cslyr,z) = soilcarb_age_pool(i,dim_cslyr,z) +      &
						(z_eff(dim_cslyr) - zsoil_ext(dim_cslyr)) /              &
						(z_eff(dim_cslyr) - z_eff(dim_cslyr-1)) *                &
						(soilcarb_age_pool(i,dim_cslyr-1,z) -                   &
						soilcarb_age_pool(i,dim_cslyr,z))
		end do

		do j = 2,dim_cslyr-1
			if ( z_eff(j) < zsoil_ext(j) ) then
				do z = 1,4
					cdens(j,z) = cdens_eff(j,z) + 0.75 * (z_eff(j) - zsoil_ext(j)) /      &
							(z_eff(j) - z_eff(j+1)) * (cdens_eff(j+1,z) - cdens_eff(j,z))&
							+ 0.25 * (z_eff(j) - zsoil_ext(j)) /                         &
							(z_eff(j) - z_eff(j-1)) * (cdens_eff(j-1,i) - cdens_eff(j,z))
					if (z<3) then
					ndens(j,z) = ndens_eff(j,z) + 0.75 * (z_eff(j) - zsoil_ext(j)) /    &
							(z_eff(j) - z_eff(j+1)) * (ndens_eff(j+1,i) - ndens_eff(j,z))&
							+ 0.25 * (z_eff(j) - zsoil_ext(j)) /                         &
							(z_eff(j) - z_eff(j-1)) * (ndens_eff(j-1,z) - ndens_eff(j,z))
					end if
					soilage_tmp(j,z) = soilcarb_age_pool(i,j,z) + 0.75 * (z_eff(j) -   &
							zsoil_ext(j)) / (z_eff(j) - z_eff(j+1)) *                    &
							(soilcarb_age_pool(i,j+1,z) - soilcarb_age_pool(i,j,z))&
							+ 0.25 * (z_eff(j) - zsoil_ext(j)) /                         &
							(z_eff(j) - z_eff(j-1)) * (soilcarb_age_pool(i,j-1,z) -   &
							soilcarb_age_pool(i,j,z))
				end do
			else ! z_eff >= zsoil_ext
				do z = 1,4
					cdens(j,z) = cdens_eff(j,z) + 0.75 * (z_eff(j) - zsoil_ext(j)) /      &
							(z_eff(j) - z_eff(j-1)) * (cdens_eff(j-1,z) - cdens_eff(j,z))&
							+ 0.25 * (z_eff(j) - zsoil_ext(j)) /                          &
							(z_eff(j) - z_eff(j+1)) * (cdens_eff(j+1,z) - cdens_eff(j,z))
					if (z<3) then
						ndens(j,z) = ndens_eff(j,z) + 0.75 * (z_eff(j) - zsoil_ext(j)) /    &
							(z_eff(j) - z_eff(j-1)) * (ndens_eff(j-1,z) - ndens_eff(j,z))&
							+ 0.25 * (z_eff(j) - zsoil_ext(j)) /                          &
							(z_eff(j) - z_eff(j+1)) * (ndens_eff(j+1,z) - ndens_eff(j,z))
					end if
					soilage_tmp(j,z) = soilcarb_age_pool(i,j,z) + 0.75 * (z_eff(j) -   &
							zsoil_ext(j)) / (z_eff(j) - z_eff(j-1)) *                    &
							(soilcarb_age_pool(i,j-1,z) - soilcarb_age_pool(i,j,z))&
							+ 0.25 * (z_eff(j) - zsoil_ext(j)) /                         &
							(z_eff(j) - z_eff(j+1)) * (soilcarb_age_pool(i,j+1,z) -   &
							soilcarb_age_pool(i,j,z))
				end do
			end if
		end do

		! Make sure volume in layer does not exceed allowed volume in layer.
		do j = 1,dim_cslyr
			do z = 1,4
				cdens(j,z) = MAX( cdens(j,z), cs_min / dzsoil_ext(j) ) 
			if(z<3) then
				ndens(j,z) = MAX( ndens(j,z), cs_min /         &
														(300.0 * dzsoil_ext(j)) )
			end if
				soilage_tmp(j,z) = MIN( MAX( soilage_tmp(j,z), 0.001 ), soilage_max )
			end do

			totvol = SUM(cdens(j,1:2)) / ( rho_dpm_rpm / compress_dr(j) )              &
				+ SUM(cdens(j,3:4)) / ( rho_bio_hum / compress_bh(j) )
			if ( totvol > 1.0 ) then
				cdens(j,:) = cdens(j,:) / totvol
				ndens(j,:) = ndens(j,:) / totvol
			end if

			! Soil C and N in layer
			cs_tmp(j,:) = cdens(j,:) * dzsoil_ext(j)
			ns_tmp(j,:) = ndens(j,:) * dzsoil_ext(j)

		end do !dim_cslyr

		! Normalise to correct amount of soil C and N, and average age
		do z=1,4
			cs(i,:,z) = cs_tmp(:,z) * SUM( cs(i,:,z) + ( pc(i,:,z)  ) ) /    &
						SUM( cs_tmp(:,z) )
			if (z<3) then
				ns_dpm_rpm(i,:,z) = ns_tmp(:,z) * SUM( ns_dpm_rpm(i,:,z) ) /            &
								SUM (ns_tmp(:,z))
			end if
			if( SUM(cdens_eff(:,z)*dz_eff(:)) > cs_min ) then
				soilage_avg(z) = SUM(cdens_eff(:,z)*dz_eff(:)*soilcarb_age_pool(i,:,z)) &
							/ SUM(cdens_eff(:,z)*dz_eff(:))
			else
				soilage_avg(z) = 0.0
			end if
		end do

		do j = 1,dim_cslyr
			do z = 1,4
				cs(i,j,z) = MAX(cs_min,cs(i,j,z))
			end do
		end do

		!! Carbon quantity should be conserved and therefore so should soil age
		do z=1,4
			soilage_tmp(:,z) = soilage_tmp(:,z) * SUM( cs(i,:,z) ) /                  &
							SUM( soilage_tmp(:,z) * cs(i,:,z) ) * soilage_avg(z)
		end do
			soilcarb_age_pool(i,:,:) = soilage_tmp(:,:)

		end do 
		
	end subroutine decay_withacc

	subroutine calc_organic_c(n_sp, dim_cslyr, dzsoil, cs_dpm_dz,cs_rpm_dz,cs_bio_dz,cs_hum_dz, &
		soil_minvol_soilt_1pt,                   &
		soil_minmass_soilt_1pt,                  &
		compress_dr, compress_bh, f_org)
		
		INTEGER, INTENT(IN) :: n_sp, dim_cslyr
		REAL(8), INTENT(IN), DIMENSION(n_sp,dim_cslyr) :: cs_dpm_dz,cs_rpm_dz,cs_bio_dz,cs_hum_dz
		REAL(8), INTENT(IN), DIMENSION(dim_cslyr) :: compress_dr, & 
											compress_bh, dzsoil
		REAL(8), INTENT(INOUT), DIMENSION(dim_cslyr) :: soil_minvol_soilt_1pt, &
												soil_minmass_soilt_1pt, f_org
		INTEGER :: i, j, z 
		REAL(8), DIMENSION(dim_cslyr) :: bd_org, var_org
		REAL(8), DIMENSION(dim_cslyr,4) :: cs_s
		REAL(8) :: cs_min, rho_dpm_rpm, rho_bio_hum, dzsoil_extra_1pt	
		
		rho_dpm_rpm = 19.6
		rho_bio_hum = 118.d0
		cs_min = 1.0e-6
		dzsoil_extra_1pt = 0.0

		do j=1, dim_cslyr
			cs_s(j,1) = SUM(cs_dpm_dz(:,j))
			cs_s(j,2) = SUM(cs_rpm_dz(:,j))
			cs_s(j,3) = SUM(cs_bio_dz(:,j))
			cs_s(j,4) = SUM(cs_hum_dz(:,j))
		end do

		do j = 1,dim_cslyr
			f_org(j) = ( (cs_s(j,1) + cs_s(j,2)) / (rho_dpm_rpm / compress_dr(j)) +     &
				(cs_s(j,3) + cs_s(j,4)) / (rho_bio_hum / compress_bh(j)) ) / dzsoil(j)
		if (f_org(j) > 1.0) then
			f_org(j) = 1.0
		end if

		if (f_org(j) >= cs_min) then
			bd_org(j) =  (1 / 0.56) * SUM( cs_s(j,:) ) / ( f_org(j) * dzsoil(j) )
		else
			bd_org(j) =  (1 / 0.56)  / ( dzsoil(j) )
		end if

		if (bd_org(j) > 210.0) then
			bd_org(j) = 210.0
		end if

			dzsoil_extra_1pt = dzsoil_extra_1pt + (cs_s(j,1) + cs_s(j,2)) / &
			(rho_dpm_rpm / compress_dr(j)) + &
			(cs_s(j,3) + cs_s(j,4)) / (rho_bio_hum / compress_bh(j))
		end  do

		! correct f_org(dim_cslyr)
		f_org(dim_cslyr) = ( (cs_s(dim_cslyr,1) + cs_s(dim_cslyr,2)) / &
							(rho_dpm_rpm / compress_dr(dim_cslyr)) + &
							(cs_s(dim_cslyr,3) + cs_s(dim_cslyr,4)) / &
							(rho_bio_hum / compress_bh(dim_cslyr)) ) / &
							( dzsoil(dim_cslyr) + dzsoil_extra_1pt )
		if (f_org(dim_cslyr) > 1.0) then
			f_org(dim_cslyr) = 1.0
		end if


		if (f_org(dim_cslyr) >= cs_min) then
			bd_org(dim_cslyr) =  (1 / 0.56) * SUM( cs_s(dim_cslyr,:) ) /               &
					( f_org(dim_cslyr) * ( dzsoil(dim_cslyr) + dzsoil_extra_1pt) )
		else
			bd_org(dim_cslyr) =  (1 / 0.56)  / ( dzsoil(dim_cslyr) + dzsoil_extra_1pt )
		end if

		if (bd_org(dim_cslyr) > 210.0) then
			bd_org(dim_cslyr) = 210.0
		end if
		


		! Calculate mineral volume fraction and mineral mass fraction
		soil_minvol_soilt_1pt = 1.0 - f_org
		! remove magic numbers here
		soil_minmass_soilt_1pt = soil_minvol_soilt_1pt * 1500.0 /                      &
						(soil_minvol_soilt_1pt * 1500.0 + bd_org * f_org)

		!var_org = 0.0

		!CASE ( 'bexp' )
		!  do n = 1,dim_cslyr
		!    var_org = bd_org(j) * 0.0304 + 1.53
		!    var_out(j) = var_org * f_org(j) + var_in(j) * (1 - f_org(j))
		!  end  do

		!CASE ( 'sathh' )
		!  do n = 1,dim_cslyr
		!    var_org = EXP( bd_org(j) * 0.023 - 5.08 )
		!    var_out(j) = ( var_org ** f_org(j) ) * ( var_in(j)** (1 - f_org(j)) )
		!  end  do

		!CASE ( 'satcon' )
		!  do n = 1,dim_cslyr
		!    var_org = EXP( bd_org(j) * (-0.0532) -6.63 ) * 1000.0
		!    var_out(j) = ( var_org ** f_org(j) ) * ( var_in(j)** (1 - f_org(j)) )
		!  end  do

		!CASE ( 'sm_sat' )
		!  do n = 1,dim_cslyr
		!    var_org = 1.0 - bd_org(j) / 1260.0
		!    var_out(j) = var_org * f_org(j) + var_in(j) * (1 - f_org(j))
		!  end  do

		!CASE ( 'hcon' )
		!  do n = 1,dim_cslyr
		!    var_out(j) = ( 0.06 ** f_org(j) ) * ( var_in(j)** (1 - f_org(j)) )
		!  end  do

		!CASE ( 'hcap' )
		!  do n = 1,dim_cslyr
		!    var_org = bd_org(j) * 2500.0 / 1.26
		!    var_out(j) = var_org * f_org(j) + var_in(j) * (1 - f_org(j))
		!  end  do
		
	end subroutine calc_organic_c

!-------------------------------------------------------------------------------------

end module soil_cnp_subroutines