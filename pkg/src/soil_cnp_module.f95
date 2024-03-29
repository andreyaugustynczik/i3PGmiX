module soil_cnp

use soil_cnp_subroutines

IMPLICIT NONE
public :: run_icbm, mod5c, mod5c_year, compAWENH, yasso_outflux, yasso_outflux_year, &
            yasso_nav, run_rothC, yasso_wetland 

contains

    !Subroutines for ICBM/2N--------------------------------------------------------------
    subroutine run_icbm(krmax, klmax, komax, Yl_C, Yr_C, O_C, Yl_N, Yr_N, O_N, Litter_C, &
                        hc, el, er, qbc, qh, qil, qir, Yl_input, Yr_input, O_Coutflux, &
                        Yl_Coutflux, Yr_Coutflux, Nav, TotalCarbo, TotalNitro, fol_litter, &
                        temp, pr, excessSW, asw, aswmax, soil_class, n_depo, dmC, n_sp)
        implicit none		 
        integer, intent(in):: n_sp, soil_class
        integer :: i
        real(kind=8), dimension(n_sp), intent(in):: krmax, klmax, qil, qir, fol_litter, &
                                                    hc, el, er, qbc, qh, Yl_input, Yr_input
        real(kind=8), intent(out):: O_C, O_N, Nav, TotalCarbo, TotalNitro, O_Coutflux	
        real(kind=8), dimension(n_sp), intent(out):: Yl_C, Yr_C, Yl_N, Yr_N, Litter_C, &
                                                    Yl_Coutflux, Yr_Coutflux											
        real(kind=8), intent(in):: komax, dmC, temp, pr, n_depo, excessSW, asw, aswmax						
        real(kind=8), dimension(n_sp) :: Litter_Coutflux, humification_Litter, &
                                        humification_l, humification_r, Yl_Noutflux, humification_N_l, &
                                        humification_N_r, Yr_Noutflux
        real(kind=8) :: leaching_org, leaching_min, water_ratio, kl, kr, ko, O_Noutflux								
        
        water_ratio = asw / aswmax
        
        do i=1, n_sp
            !Adjust decomposition rates based on iLand implementation (Seidl et al. 2012)
            kr = krmax(i) * (1.d0/(1.d0+30.d0*exp(-8.5d0*(water_ratio)))) * exp(5.508032845d0 - 308.56d0/(273.15d0+temp-227.13d0)) 
            kl = klmax(i) * (1.d0/(1.d0+30.d0*exp(-8.5d0*(water_ratio)))) * exp(5.508032845d0 - 308.56d0/(273.15d0+temp-227.13d0)) 

            !Compute monthly out fluxes
            Yl_Coutflux(i) = kl * (1.d0 - hc(i)) * Yl_C(i)  
            Yr_Coutflux(i) = kr * (1.d0 - hc(i)) * Yr_C(i) 
            Litter_Coutflux(i) = kl * (1.d0 - hc(i)) * Litter_C(i) 
            humification_Litter(i) = kl * hc(i) * Litter_C(i)
            humification_l(i) = kl * hc(i) * Yl_C(i)
            humification_r(i) = kr * hc(i) * Yr_C(i)
            
            !Compute N out fluxes
            Yl_Noutflux(i) = kl * ((1.d0 - hc(i)) / (1.d0 - el(i))) * (Yl_N(i) - el(i) * (Yl_C(i) / qbc(i)))
            humification_N_l(i) = kl * hc(i) * (Yl_N(i) / qh(i))
            humification_N_r(i) = kr * hc(i) * (Yr_N(i) / qh(i))
            Yr_Noutflux(i) = kr * ((1.d0 - hc(i)) / (1.d0 - er(i))) * (Yr_N(i) - er(i) * (Yr_C(i) / qbc(i)))

            !Now calculate the end-of-month carbon and nitrogen pools				
            Yr_C(i) = Yr_C(i) + Yr_input(i) - Yr_Coutflux(i) - humification_r(i)  
            Yl_C(i) = Yl_C(i) + Yl_input(i) - Yl_Coutflux(i) - humification_l(i) 
            Litter_C(i) = Litter_C(i) + fol_litter(i) - Litter_Coutflux(i) - humification_Litter(i)

            Yr_N(i) = Yr_N(i) + ((Yr_input(i)) / (2.d0 * qir(i))) - Yr_Noutflux(i) - humification_N_r(i)
            Yl_N(i) = Yl_N(i) + ((Yl_input(i)) / (2.d0 * qil(i))) - Yl_Noutflux(i) - humification_N_l(i)
        end do

        !Repeat for the old carbon and nitrogen pools
        ko = komax * (1.d0/(1.d0+30.d0*exp(-8.5d0*(water_ratio)))) * exp(5.508032845d0 - 308.56d0/(273.15d0+temp-227.13d0)) 

        O_Coutflux = ko * O_C  
        O_Noutflux = ko * O_N

        O_C = O_C + sum(humification_l(:)) + sum(humification_r(:)) - O_Coutflux
        O_N = O_N + sum(humification_N_r(:)) + sum(humification_N_l(:)) - O_Noutflux

        TotalCarbo = sum(Yr_C(:) + Yl_C(:)) + O_C
        TotalNitro = sum(Yr_N(:) + Yl_N(:)) + O_N

        !Calculate N leaching
        call compute_leaching( soil_class, leaching_org, leaching_min, excessSW, asw )
        
        !Get available soil N
        Nav = sum(Yr_Noutflux(:)) * (1.d0-leaching_org)  + sum(+ Yl_Noutflux(:)) * (1.d0-leaching_org) + &
                        O_Noutflux * (1.d0 - leaching_min)  + n_depo


    end subroutine run_icbm
    ! End subroutines for ICBM/2N---------------------------------------------------------

    subroutine run_rothC(temp, temp_year, ASW, asw_max, snow_water, lai, month, soil_class, &
        Q10_soil, Yr_C, Yl_C, dpm_ratio,  &
        resp_soc, soc, soc_int, n_sp, biom_loss_foliage, biom_loss_root, mF, &
        mort_stress, mort_thinn, biom_foliage, &
        stems_n, mR, biom_root, dmC, mS, biom_stem, n_depo, cn_bio, cn_hum, cp_bio, cp_hum, &
        k_bio, k_hum, fn_l, lit_cn, lit_cp, cn_veg, cp_veg, ns, ns_int, &
        ps, ps_int, kaps_roth_bio, kaps_roth_hum, &
        kaps_roth_dpm, kaps_roth_rpm, npr_soil, lch_dis, excessSW, Nav, Pav, TotalCarbo,   &
        TotalCarbo_dz_10cm, TotalCarbo_dz_35cm, TotalCarbo_dz_100cm, TotalCarbo_dz_300cm, & 
        Total_resp, TotalNitro,&
        TotalP, Yl_input, Yr_input, cs_dpm, cs_rpm, cp_dpm, cp_rpm, cs_bio, &
        cs_hum, cs_dpm_dz, cs_rpm_dz, &
        cs_bio_dz, cs_hum_dz, n_inorg_T, minl_n_tot, &
        immob_n_tot, n_gas, n_lch, p_inorg_T, minl_p_tot, immob_p_tot, p_lch, &
        peat_io, anoxic_io, dim_cslyr, &
        pc, cs, dz_extra, ns_dpm_rpm, soilcarb_age_pool, resp_s, resp_frac, dcs, dzsoil, lit_frac, &
        mix_s, mix_term, zz_l, soilage_max, compress_dr, compress_bh, soil_minvol_soilt_1pt, &
        soil_minmass_soilt_1pt, f_org)

        integer, intent(in):: n_sp, soil_class, month        
        integer :: i,j,z
        real(kind=8), intent(in) :: Q10_soil
        integer, intent(inout) :: peat_io, anoxic_io, dim_cslyr
        real(8), dimension(n_sp,dim_cslyr,2), intent(inout) :: ns_dpm_rpm
        real(kind=8), dimension(n_sp), intent(in) :: dpm_ratio, kaps_roth_dpm,kaps_roth_rpm,kaps_roth_bio,kaps_roth_hum
        real(kind=8), dimension(n_sp), intent(in) :: cn_bio, cn_hum , cp_bio, cp_hum, lit_cn, &
                                            lit_cp, cn_veg, cp_veg, lai, biom_loss_foliage, biom_loss_root, mF, &
                                            mort_stress, mort_thinn, biom_foliage, stems_n, mR, biom_root, mS, &
                                            biom_stem, Yl_input, Yr_input

        real(kind=8), dimension(n_sp), intent(out) ::  Yr_C, Yl_C
        real(kind=8), dimension(dim_cslyr), intent(out) ::  lit_frac, dzsoil
        real(kind=8), dimension(n_sp,dim_cslyr,4), intent(inout) ::  pc, cs, soilcarb_age_pool, resp_s, &
                mix_term, mix_s, resp_frac, dcs

        real (kind=8), intent(in) :: n_depo, ASW, asw_max, snow_water, dmC, temp, temp_year, excessSW
        real (kind=8), intent(in) :: soc_int, ns_int, ps_int
        real (kind=8), intent(inout), dimension(dim_cslyr)  :: compress_dr, compress_bh
        real (kind=8), intent(inout):: soilage_max
        real (kind=8), dimension(n_sp,dim_cslyr), intent(inout) ::cs_dpm_dz, cs_rpm_dz, cs_bio_dz, cs_hum_dz
        real (kind=8), dimension(dim_cslyr), intent(inout) ::zz_l, soil_minvol_soilt_1pt, soil_minmass_soilt_1pt, f_org
        real (kind=8), intent(in) :: lch_dis, npr_soil
        real (kind=8), dimension(n_sp) ::  resp_dpm, resp_rpm, resp_bio, resp_hum, dec_dpm_rpm, cn_dpm, cn_rpm
        real (kind=8), intent(inout), dimension(n_sp) :: cs_dpm, cs_rpm, cp_dpm, cp_rpm, cs_bio, cs_hum, dz_extra
        real (kind=8), intent(out), dimension(n_sp) :: n_inorg_T, &
                minl_n_tot, immob_n_tot, n_lch, n_gas, p_inorg_T, &
                minl_p_tot, immob_p_tot, p_lch
        real (kind=8), dimension(n_sp) :: ns_dpm, ns_rpm, ns_bio, ns_hum 
        real (kind=8), dimension(4) :: kaps_roth
        real (kind=8) :: soc, resp_soc, alpha_moist, alpha_temp, nminl_gas, k_dpm, k_rpm, k_bio, k_hum, &
                leaching_min, leaching_org, t_soil, t_soil_surf, isunfrozen, litc_norm, tau_lit, tau_resp
        real (kind=8), dimension(n_sp) :: minl_n_dpm_pot, minl_n_rpm_pot, minl_n_bio_pot, minl_n_hum_pot, minl_n_pot_tot, fn_l
        real (kind=8), dimension(n_sp) :: minl_n_dpm, minl_n_rpm, minl_n_bio, minl_n_hum
        real (kind=8), dimension(n_sp) :: immob_n_dpm_pot, immob_n_rpm_pot, immob_n_bio_pot, immob_n_hum_pot, immob_n_pot_tot
        real (kind=8), dimension(n_sp) :: immob_n_dpm, immob_n_rpm, immob_n_bio, immob_n_hum
        real (kind=8), dimension(n_sp) ::  ps_dpm, ps_rpm, ps_bio, ps_hum, fp
        real (kind=8), dimension(n_sp) :: minl_p_dpm_pot, minl_p_rpm_pot, minl_p_bio_pot, minl_p_hum_pot, minl_p_pot_tot
        real (kind=8), dimension(n_sp) :: minl_p_dpm, minl_p_rpm, minl_p_bio, minl_p_hum
        real (kind=8), dimension(n_sp) :: immob_p_dpm_pot, immob_p_rpm_pot, immob_p_bio_pot, immob_p_hum_pot, immob_p_pot_tot
        real (kind=8), dimension(n_sp) :: immob_p_dpm, immob_p_rpm, immob_p_bio, immob_p_hum
        real (kind=8), intent(out) :: Nav, Pav, TotalCarbo, Total_resp, TotalNitro, TotalP, ns, ps 
        real (kind=8), intent(out) :: TotalCarbo_dz_10cm, TotalCarbo_dz_35cm, TotalCarbo_dz_100cm, TotalCarbo_dz_300cm

        isunfrozen = 1.d0
        litc_norm  = 1.d0
        tau_lit    = 5.d0
        tau_resp   = 1.d0


        ! Define 4 soil layers (this should be later coded to be dynamic)
        ! dim_cslyr : is the number of soil layers that can differ based on the depth 
        dzsoil(1) = 0.1
        dzsoil(2) = 0.25
        dzsoil(3) = 0.65
        dzsoil(4) = 2.d0

        ! Call Moisture Rate Modifier Subroutine
        CALL moisture_rate_modifier(ASW, asw_max, alpha_moist, peat_io)

        ! Call Temperature Rate Modifier Subroutine
        CALL temperature_rate_modifier(t_soil, temp, temp_year, ASW, snow_water, lai, month, soil_class, alpha_temp, Q10_soil)

        IF (t_soil < 0.d0) isunfrozen = 0.d0
        ! Carbon Pools Calculation
        do i=1, n_sp
            k_dpm = dpm_ratio(i)                    !Distribution coef. of litter to DPM pool
            k_rpm = 1.d0-dpm_ratio(i)       !Distribution coef. of litter to RPM pool
            k_bio =  0.46d0                 !Distribution coef. to Biomass pool
            k_hum =  0.54d0                 !Distribution coef. to Humus pool
            kaps_roth(1)=k_dpm
            kaps_roth(2)=k_rpm
            kaps_roth(3)=k_bio
            kaps_roth(4)=k_hum

            nminl_gas = 0.01d0              !Fraction of net mineralization of N that is lost as gas

            ! C pools - based on a fixed ratio from litter

            cs_dpm(i) =  cs_dpm(i) + (Yl_input(i))
            resp_dpm(i) =  kaps_roth_dpm(i) * cs_dpm(i) * alpha_temp * alpha_moist  
            cs_dpm(i) = cs_dpm(i) - resp_dpm(i)

            cs_rpm(i) =  cs_rpm(i) + (Yr_input(i))
            resp_rpm(i) = kaps_roth_rpm(i) * cs_rpm(i) * alpha_temp * alpha_moist
            cs_rpm(i) = cs_rpm(i) - resp_rpm(i)

            dec_dpm_rpm(i) =  resp_dpm(i) + resp_rpm(i)

            cs_bio(i) =  cs_bio(i) + (k_bio * dec_dpm_rpm(i))
            resp_bio(i) = kaps_roth_bio(i) * cs_bio(i) * alpha_temp * alpha_moist
            cs_bio(i) = cs_bio(i) - resp_bio(i)

            cs_hum(i) =  cs_hum(i) + (k_hum * dec_dpm_rpm(i))
            resp_hum(i) = (kaps_roth_hum(i)) * cs_hum(i) * alpha_temp * alpha_moist
            cs_hum(i) = cs_hum(i) - resp_hum(i)

            !-----------------------------------------------------------------------------
            ! Calculate vertical profile of litter inputs.
            !-----------------------------------------------------------------------------
                    
            lit_frac(1) = dzsoil(1) *  EXP( -tau_lit * 0.5 * dzsoil(1) ) / litc_norm
            do j = 2,dim_cslyr
                lit_frac(j) = dzsoil(j) * EXP( -tau_lit *                                   &
                        ( SUM(dzsoil(1:j-1)) + 0.5 * dzsoil(j) )  ) / litc_norm
            end do

            !-----------------------------------------------------------------------------
            ! Calculate vertical profile of Soil C and resp
            !-----------------------------------------------------------------------------

            do j=1, dim_cslyr
                cs(i,j,1)= cs_dpm(i) * zz_l(j)
                cs(i,j,2)= cs_rpm(i) * zz_l(j)
                cs(i,j,3)= cs_bio(i) * zz_l(j)
                cs(i,j,4)= cs_hum(i) * zz_l(j)
            end do

            do z=1, 4
                resp_s(i,1,z) = kaps_roth(z) * cs(i,1,z) * alpha_temp * alpha_moist  &
                            * (EXP(-0.5 * dzsoil(1) * tau_resp))
                do j = 2,dim_cslyr
                    resp_s(i,j,z) = kaps_roth(z) * cs(i,j,z) * alpha_temp * alpha_moist  &
                                * (EXP( - (SUM(dzsoil(1:(j-1)))              &
                                + 0.5 * dzsoil(j)) * tau_resp))
                end do
            end do

            !-------------------------------------------------------------------------
            ! Diagnose the net local carbon flux into the soil
            !-------------------------------------------------------------------------
            do j = 1,dim_cslyr
                pc(i,j,1) = (Yl_input(i) * lit_frac(j)) - resp_s(i,j,1)
                pc(i,j,2) = (Yr_input(i)* lit_frac(j)) - resp_s(i,j,2)
                pc(i,j,3) = (k_bio * dec_dpm_rpm(i)) - resp_s(i,j,3)
                pc(i,j,4) = (k_hum * dec_dpm_rpm(i)) - resp_s(i,j,4)
            end do

            !Now calculate the end-of-month carbon pools
            Yl_C(i) =  cs_dpm(i)
            Yr_C(i) =  cs_rpm(i)


            soc =    soc_int + (sum(cs_dpm(:) + cs_rpm(:) + cs_bio(:) + cs_hum(:)))/n_sp

            !---------------------------------------------------------------------------
            ! Calculate a mixing term
            !---------------------------------------------------------------------------
            CALL dz_soil_cs(n_sp, dim_cslyr, zz_l, dzsoil, cs_dpm, cs_rpm, cs_bio, cs_hum, &
                cs_dpm_dz, cs_rpm_dz, cs_bio_dz, cs_hum_dz, t_soil, mix_term, mix_s, month, dz_extra)

            !-------------------------------------------------------------------------
            ! Save current value of soil carbon
            !-------------------------------------------------------------------------
            dcs(:,:,1) = cs_dpm_dz(:,:)
            dcs(:,:,2) = cs_rpm_dz(:,:)
            dcs(:,:,3) = cs_bio_dz(:,:)
            dcs(:,:,4) = cs_hum_dz(:,:)


            ! Nitrogen Pools Calculation

            !Initial N pool (BIO + HUM) are diagnostic from equivalent C soil pool
            cn_dpm(i) = (dpm_ratio(i) * cn_veg(i))
            cn_rpm(i) = ((1-dpm_ratio(i)) * cn_veg(i))
            ns_bio(i) = cs_bio(i) / cn_bio(i)
            ns_hum(i) = cs_hum(i) / cn_hum(i)

            !First calculate the potential mineralized and immobilized N
            minl_n_dpm_pot(i)  = resp_dpm(i) / cn_dpm(i)
            minl_n_rpm_pot(i)  = resp_rpm(i) / cn_rpm(i)
            minl_n_bio_pot(i)  = resp_bio(i) / cn_bio(i)
            minl_n_hum_pot(i)  = resp_hum(i) / cn_hum(i)
            minl_n_pot_tot(i)  = minl_n_dpm_pot(i) + minl_n_rpm_pot(i) + minl_n_bio_pot(i) + minl_n_hum_pot(i)

            immob_n_dpm_pot(i) = (k_bio * resp_dpm(i) / cn_bio(i)) + (k_hum * resp_dpm(i) / cn_hum(i))
            immob_n_rpm_pot(i) = (k_bio * resp_rpm(i) / cn_bio(i)) + (k_hum * resp_rpm(i) / cn_hum(i))
            immob_n_bio_pot(i) = (k_bio * resp_bio(i) / cn_bio(i)) + (k_hum * resp_bio(i) / cn_hum(i))
            immob_n_hum_pot(i) = (k_bio * resp_hum(i) / cn_bio(i)) + (k_hum * resp_hum(i) / cn_hum(i))
            immob_n_pot_tot(i) = immob_n_dpm_pot(i) + immob_n_rpm_pot(i) + immob_n_bio_pot(i) + immob_n_hum_pot(i)

            ! Calculate the demand vs supply for organic and inorganic N
            ! include unfrozen frac 
            if ((immob_n_pot_tot(i) - minl_n_pot_tot(i))> (n_inorg_T(i)) * isunfrozen) then
                fn_l(i) = (((minl_n_bio_pot(i) + minl_n_hum_pot(i) - immob_n_bio_pot(i) - immob_n_hum_pot(i)) + &
                        n_inorg_T(i)) * isunfrozen) / ((immob_n_dpm_pot(i) + immob_n_rpm_pot(i) - &
                        minl_n_dpm_pot(i) - minl_n_rpm_pot(i)))
                fn_l(i) = MIN(MAX(fn_l(i),0.d0),1.d0)
            else
                fn_l(i) = 1.d0
            end if

            !Final minl. and immob. based on limiting N
            minl_n_dpm(i)  = minl_n_dpm_pot(i) * fn_l(i)
            minl_n_rpm(i)  = minl_n_rpm_pot(i) * fn_l(i)
            minl_n_bio(i)  = minl_n_bio_pot(i)
            minl_n_hum(i)  = minl_n_hum_pot(i)

            immob_n_dpm(i)  = immob_n_dpm_pot(i) * fn_l(i)
            immob_n_rpm(i)  = immob_n_rpm_pot(i) * fn_l(i)
            immob_n_bio(i)  = immob_n_bio_pot(i)
            immob_n_hum(i)  = immob_n_hum_pot(i)

            immob_n_tot(i) = immob_n_dpm(i) + immob_n_rpm(i) + immob_n_bio(i) + immob_n_hum(i)
            minl_n_tot(i) = minl_n_dpm(i) + minl_n_rpm(i) + minl_n_bio(i) + minl_n_hum(i)

            ! C decomposition correction based on N content for microbial activities (Craine et al., 2017)
            resp_dpm(i) = resp_dpm(i) * fn_l(i)
            resp_rpm(i) = resp_rpm(i) * fn_l(i)
            resp_bio(i) = resp_bio(i) * fn_l(i)
            resp_hum(i) = resp_hum(i) * fn_l(i)

            resp_soc = (sum(resp_dpm(:) + resp_rpm(:) + resp_bio(:) + resp_hum(:)))/n_sp


            ns_dpm(i) = dpm_ratio(i) * ((Yl_C(i)+Yr_C(i)) / lit_cn(i)) - minl_n_dpm(i)
            ns_rpm(i) = (1-dpm_ratio(i)) * ((Yl_C(i)+Yr_C(i)) / lit_cn(i)) - minl_n_rpm(i)
            ns_bio(i) = k_bio * immob_n_tot(i) - minl_n_bio(i)
            ns_hum(i) = k_hum * immob_n_tot(i) - minl_n_hum(i)


            n_inorg_T(i) = n_inorg_T(i) + (minl_n_dpm(i) + minl_n_rpm(i) + minl_n_bio(i) + minl_n_hum(i)) - &
                    (immob_n_dpm(i) + immob_n_rpm(i) + immob_n_bio(i) + immob_n_hum(i)) + &
                    n_depo


            ! Call Nitrogen Gas and Leaching Subroutine
            call nitrogen_gas(n_sp, n_inorg_T, n_gas, nminl_gas, minl_n_tot, immob_n_tot)


            call compute_leaching( soil_class, leaching_org, leaching_min, excessSW, asw )

            n_lch(i)   = n_inorg_T(i)  * (1.d0-leaching_org) / npr_soil

            if ((n_inorg_T(i)-n_lch(i)) .GT. 0.d0) then
                n_inorg_T(i) = n_inorg_T(i) - n_lch(i)
            else
                n_lch(i) = 0.d0
            end if

            ns = ns_int + sum(ns_dpm(:) + ns_rpm(:) + ns_bio(:) + ns_hum(:)+ n_inorg_T(:))/n_sp


            ns_dpm_rpm(:,:,1) = ns_dpm(i)
            ns_dpm_rpm(:,:,2) = ns_rpm(i)

            IF (peat_io .eq. int(1)) THEN

                CALL soilcarb_age(n_sp, dim_cslyr, dz_extra, mix_s, resp_s, resp_frac, &
                    soilcarb_age_pool, dcs, Yr_C, Yl_C, lit_frac, dpm_ratio, dzsoil, soilage_max, month)
                CALL decay_withacc(n_sp, dim_cslyr, dz_extra, pc, cs, soilcarb_age_pool, &
                    ns_dpm_rpm, soilage_max, dzsoil, compress_dr, compress_bh)
            END IF 

            !Phosphorus Pools Calculation

            ! Phosphorus Pools Calculation

            !Initial P pool (BIO + HUM) are diagnostic from equivalent C soil pool
            cp_dpm(i) = (dpm_ratio(i) * cp_veg(i))
            cp_rpm(i) = ((1-dpm_ratio(i)) * cp_veg(i))
            ps_bio(i) = cs_bio(i) / cp_bio(i)
            ps_hum(i) = cs_hum(i) / cp_hum(i)

            !First calculate the potential mineralized and immobilized N
            minl_p_dpm_pot(i)  = resp_dpm(i) / cp_dpm(i)
            minl_p_rpm_pot(i)  = resp_rpm(i) / cp_rpm(i)
            minl_p_bio_pot(i)  = resp_bio(i) / cp_bio(i)
            minl_p_hum_pot(i)  = resp_hum(i) / cp_hum(i)
            minl_p_pot_tot(i)  = minl_p_dpm_pot(i) + minl_p_rpm_pot(i) + minl_p_bio_pot(i) + minl_p_hum_pot(i)

            immob_p_dpm_pot(i) = (k_bio * resp_dpm(i) / cp_bio(i)) + (k_hum * resp_dpm(i) / cp_hum(i))
            immob_p_rpm_pot(i) = (k_bio * resp_rpm(i) / cp_bio(i)) + (k_hum * resp_rpm(i) / cp_hum(i))
            immob_p_bio_pot(i) = (k_bio * resp_bio(i) / cp_bio(i)) + (k_hum * resp_bio(i) / cp_hum(i))
            immob_p_hum_pot(i) = (k_bio * resp_hum(i) / cp_bio(i)) + (k_hum * resp_hum(i) / cp_hum(i))
            immob_p_pot_tot(i) = immob_p_dpm_pot(i) + immob_p_rpm_pot(i) + immob_p_bio_pot(i) + immob_p_hum_pot(i)

            ! Calculate the demand vs supply for organic and inorganic N
            if ((immob_p_pot_tot(i) - minl_p_pot_tot(i))> (p_inorg_T(i))) then
                fp(i) = (((minl_p_bio_pot(i) + minl_p_hum_pot(i) - immob_p_bio_pot(i) - immob_p_hum_pot(i)) + &
                    p_inorg_T(i))) / ((immob_p_dpm_pot(i) + immob_p_rpm_pot(i) - &
                    minl_p_dpm_pot(i) - minl_p_rpm_pot(i)))
                fp(i) = MIN(MAX(fp(i),0.d0),1.d0)
            else
                fp(i) = 1.d0
            end if

            !Final minl. and immob. based on limiting N
            minl_p_dpm(i)  = minl_p_dpm_pot(i) * fp(i)
            minl_p_rpm(i)  = minl_p_rpm_pot(i) * fp(i)
            minl_p_bio(i)  = minl_p_bio_pot(i)
            minl_p_hum(i)  = minl_p_hum_pot(i)

            immob_p_dpm(i)  = immob_p_dpm_pot(i) * fp(i)
            immob_p_rpm(i)  = immob_p_rpm_pot(i) * fp(i)
            immob_p_bio(i)  = immob_p_bio_pot(i)
            immob_p_hum(i)  = immob_p_hum_pot(i)

            immob_p_tot(i) = immob_p_dpm(i) + immob_p_rpm(i) + immob_p_bio(i) + immob_p_hum(i)
            minl_p_tot(i) = minl_p_dpm(i) + minl_p_rpm(i) + minl_p_bio(i) + minl_p_hum(i)

            ps_dpm(i) = dpm_ratio(i) * ((Yl_C(i)+Yr_C(i)) / lit_cp(i)) - minl_p_dpm(i)
            ps_rpm(i) = (1-dpm_ratio(i)) * ((Yl_C(i)+Yr_C(i)) / lit_cp(i)) - minl_p_rpm(i)
            ps_bio(i) = k_bio * immob_p_tot(i) - minl_p_bio(i)
            ps_hum(i) = k_hum * immob_p_tot(i) - minl_p_hum(i)


            p_inorg_T(i) = p_inorg_T(i) + (minl_p_dpm(i) + minl_p_rpm(i) + minl_p_bio(i) + minl_p_hum(i)) - &
                    (immob_p_dpm(i) + immob_p_rpm(i) + immob_p_bio(i) + immob_p_hum(i)) 


            ! Call Phosphorus  Leaching Subroutine

            call compute_leaching( soil_class, leaching_org, leaching_min, excessSW, asw )

            p_lch(i)   = p_inorg_T(i)  * (1.d0-leaching_org) / npr_soil

            if ((p_inorg_T(i)-p_lch(i)) .GT. 0.d0) then
                p_inorg_T(i) = p_inorg_T(i) - p_lch(i)
            else
                p_lch(i) = 0.d0
            end if

            ps = ps_int + sum(ps_dpm(:) + ps_rpm(:) + ps_bio(:) + ps_hum(:)+ p_inorg_T(:))/n_sp

        end do 

        Pav = sum(ps_dpm(:) + ps_rpm(:) + ps_bio(:) + ps_hum(:)) - sum(p_lch(:))

        Nav = sum(ns_dpm(:) + ns_rpm(:) + ns_bio(:) + ns_hum(:)) - sum(n_lch(:)) - sum(n_gas(:)) +n_depo

        IF (peat_io .eq. int(1)) THEN

            call calc_organic_c(n_sp, dim_cslyr, dzsoil, cs_dpm_dz,cs_rpm_dz,cs_bio_dz,cs_hum_dz, &
                soil_minvol_soilt_1pt, soil_minmass_soilt_1pt, compress_dr, compress_bh, f_org)

            TotalCarbo_dz_10cm = sum(cs_dpm_dz(:,1) + cs_rpm_dz(:,1) + cs_bio_dz(:,1) + cs_hum_dz(:,1))
            TotalCarbo_dz_35cm = sum(cs_dpm_dz(:,2) + cs_rpm_dz(:,2) + cs_bio_dz(:,2) + cs_hum_dz(:,2))
            TotalCarbo_dz_100cm = sum(cs_dpm_dz(:,3) + cs_rpm_dz(:,3) + cs_bio_dz(:,3) + cs_hum_dz(:,3))
            TotalCarbo_dz_300cm = sum(cs_dpm_dz(:,4) + cs_rpm_dz(:,4) + cs_bio_dz(:,4) + cs_hum_dz(:,4))

        END IF


        TotalCarbo = soc



        !TotalCarbo_dz_10cm =  f_org(1)
        !TotalCarbo_dz_35cm =  f_org(2)
        !TotalCarbo_dz_100cm = f_org(3)
        !TotalCarbo_dz_300cm = f_org(4)

        !TotalCarbo_dz_10cm = sum(cs(:,1,:) )
        !TotalCarbo_dz_35cm = sum(cs(:,2,:) )
        !TotalCarbo_dz_100cm = sum(cs(:,3,:) )
        !TotalCarbo_dz_300cm = sum(cs(:,4,:) )


        !TotalCarbo_dz_10cm = sum(soilcarb_age_pool(:,1,1) + soilcarb_age_pool(:,1,2) + &
        !                        soilcarb_age_pool(:,1,3) + soilcarb_age_pool(:,1,4))
        !TotalCarbo_dz_35cm = sum(soilcarb_age_pool(:,2,1) + soilcarb_age_pool(:,2,2) + &
        !soilcarb_age_pool(:,2,3) + soilcarb_age_pool(:,2,4))
        !TotalCarbo_dz_100cm = sum(soilcarb_age_pool(:,3,1) + soilcarb_age_pool(:,3,2) + &
        !soilcarb_age_pool(:,3,3) + soilcarb_age_pool(:,3,4))
        !TotalCarbo_dz_300cm = sum(soilcarb_age_pool(:,4,1) + soilcarb_age_pool(:,4,2) + &
        !soilcarb_age_pool(:,4,3) + soilcarb_age_pool(:,4,4))

        Total_resp = resp_soc
        TotalNitro = ns
        TotalP     = ps


    end subroutine run_rothC
! End subroutines for rothC---------------------------------------------------------


    !Subroutines for Yasso ---------------------------------------------------------------

    SUBROUTINE mod5c(theta,time,temp,prec,init,b,d,leac,xt,steadystate_pred)
    IMPLICIT NONE
        !********************************************* &
        ! GENERAL DESCRIPTION FOR ALL THE MEASUREMENTS
        !********************************************* &
        ! returns the model prediction xt for the given parameters
        ! 1-16 matrix A entries: 4*alpha, 12*p

        ! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION

        ! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2

        ! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2

        ! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2

        ! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H

        ! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)

        ! 33-35 Woody parameters: theta_1, theta_2, r

        REAL (kind=8),DIMENSION(35),INTENT(IN) :: theta ! parameters
        REAL (kind=8),INTENT(IN) :: time,d,leac ! time,size,leaching
        !REAL (kind=8),DIMENSION(12),INTENT(IN) :: temp ! monthly mean temperatures
        REAL (kind=8),INTENT(IN) :: temp ! monthly mean temperatures
        REAL (kind=8),INTENT(IN) :: prec ! annual/monthly precipitation
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: init ! initial state
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: b ! infall
        REAL (kind=8),DIMENSION(5),INTENT(OUT) :: xt ! the result i.e. x(t)
        INTEGER, INTENT(IN) :: steadystate_pred
        ! LOGICAL,OPTIONAL,INTENT(IN) :: steadystate_pred ! set to true if ignore 'time' and compute solution
        ! in steady-state conditions (which sould give equal solution as if time is set large enough)
        REAL (kind=8),DIMENSION(5,5) :: A,At,mexpAt
        INTEGER :: i
        REAL (kind=8),PARAMETER :: pi = 3.141592653589793d0
        REAL (kind=8) :: tem,temN,temH,size_dep
        REAL (kind=8),DIMENSION(5) :: z1,z2
        REAL (kind=8),PARAMETER :: tol = 10E-12
        LOGICAL :: ss_pred

        ! IF(PRESENT(steadystate_pred)) THEN
            ! ss_pred = steadystate_pred
        ! ENDIF
        IF(steadystate_pred == 1) THEN
            ss_pred = .TRUE.
        ELSE
            ss_pred = .FALSE.
        ENDIF

        !#########################################################################
        ! Compute the coefficient matrix A for the differential equation

        ! temperature annual cycle approximation
        ! te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
        ! te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
        ! te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
        ! te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

        ! DO i = 1,4 ! Average temperature dependence
        !     tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
        !     temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
        !     temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
        ! END DO

        ! Set up climate dependence factors
        tem = 0.d0
        temN = 0.d0
        temH = 0.d0

        ! Monthly temperature dependence
        tem = tem+EXP(theta(22)*temp+theta(23)*temp**2.d0)
        temN = temN+EXP(theta(24)*temp+theta(25)*temp**2.d0)
        temH = temH+EXP(theta(26)*temp+theta(27)*temp**2.d0)
        
        ! Precipitation dependence
        tem = tem*(1.d0-EXP(theta(28)*prec/83.3333d0))
        temN = temN*(1.d0-EXP(theta(29)*prec/83.3333d0))
        temH = temH*(1.d0-EXP(theta(30)*prec/83.3333d0))
        
        ! Monthly temperature dependence
        !DO i = 1,12
        !    tem = tem+EXP(theta(22)*temp(i)+theta(23)*temp(i)**2.d0)
        !    temN = temN+EXP(theta(24)*temp(i)+theta(25)*temp(i)**2.d0)
        !    temH = temH+EXP(theta(26)*temp(i)+theta(27)*temp(i)**2.d0)
        !END DO
        
        !tem = tem*(1.d0-EXP(theta(28)*prec/1000.d0))/12
        !temN = temN*(1.d0-EXP(theta(29)*prec/1000.d0))/12
        !temH = temH*(1.d0-EXP(theta(30)*prec/1000.d0))/12

        ! Size class dependence -- no effect if d == 0.0
        size_dep = MIN(1.d0,(1.d0+theta(33)*d+theta(34)*d**2.d0)**(-ABS(theta(35))))

        ! check rare case where no decomposition happens for some compartments
        ! (basically, if no rain)
        IF (tem <= tol) THEN
            xt = init + b*time
            return
        END IF

        ! Calculating matrix A (will work ok despite the sign of alphas)
        DO i = 1,3
            A(i,i) = -ABS(theta(i))*tem*size_dep
        END DO
        A(4,4) = -ABS(theta(4))*temN*size_dep

        A(1,2) = theta(5)*ABS(A(2,2))
        A(1,3) = theta(6)*ABS(A(3,3))
        A(1,4) = theta(7)*ABS(A(4,4))
        A(1,5) = 0.0 ! no mass flows from H -> AWEN
        A(2,1) = theta(8)*ABS(A(1,1))
        A(2,3) = theta(9)*ABS(A(3,3))
        A(2,4) = theta(10)*ABS(A(4,4))
        A(2,5) = 0.0
        A(3,1) = theta(11)*ABS(A(1,1))
        A(3,2) = theta(12)*ABS(A(2,2))
        A(3,4) = theta(13)*ABS(A(4,4))
        A(3,5) = 0.0
        A(4,1) = theta(14)*ABS(A(1,1))
        A(4,2) = theta(15)*ABS(A(2,2))
        A(4,3) = theta(16)*ABS(A(3,3))
        A(4,5) = 0.0
        A(5,5) = -ABS(theta(32))*temH ! no size effect in humus
        DO i = 1,4
            A(5,i) = theta(31)*ABS(A(i,i)) ! mass flows AWEN -> H (size effect is present here)
        END DO

        ! Leaching (no leaching for humus)
        DO i = 1,4
            ! A(i,i) = A(i,i)+leac*climate(2)/1000.0
            A(i,i) = A(i,i)+leac*prec/1000.d0
        END DO

        !#########################################################################
        ! Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init

        IF(ss_pred) THEN
            ! Solve DE directly in steady state conditions (time = infinity)
            ! using the formula 0 = x'(t) = A*x + b => x = -A**-1*b
            CALL solve(-A, b, xt)
        ELSE
            ! Solve DE in given time
            z1 = MATMUL(A,init) + b * 12.d0 ! correct to the full monthly litter input
            At = A/12.d0 !A*time !At = A*t
            CALL matrixexp(At,mexpAt)
            z2 = MATMUL(mexpAt,z1) - b * 12.d0 ! correct to the full monthly litter input
            CALL solve(A,z2,xt) ! now it can be assumed A is non-singular
        ENDIF

    END SUBROUTINE mod5c


    SUBROUTINE mod5c_year(theta,time,temp,prec,init,b,d,leac,xt,steadystate_pred)
        IMPLICIT NONE
        !********************************************* &
        ! GENERAL DESCRIPTION FOR ALL THE MEASUREMENTS
        !********************************************* &
        ! returns the model prediction xt for the given parameters
        ! 1-16 matrix A entries: 4*alpha, 12*p

        ! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION

        ! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2

        ! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2

        ! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2

        ! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H

        ! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)

        ! 33-35 Woody parameters: theta_1, theta_2, r

        REAL (kind=8),DIMENSION(35),INTENT(IN) :: theta ! parameters
        REAL (kind=8),INTENT(IN) :: time,d,leac ! time,size,leaching
        REAL (kind=8),DIMENSION(12),INTENT(IN) :: temp ! monthly mean temperatures
        REAL (kind=8),INTENT(IN) :: prec ! annual precipitation
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: init ! initial state
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: b ! infall
        REAL (kind=8),DIMENSION(5),INTENT(OUT) :: xt ! the result i.e. x(t)
        INTEGER, INTENT(IN) :: steadystate_pred
        ! LOGICAL,OPTIONAL,INTENT(IN) :: steadystate_pred ! set to true if ignore 'time' and compute solution
        ! in steady-state conditions (which sould give equal solution as if time is set large enough)
        REAL (kind=8),DIMENSION(5,5) :: A,At,mexpAt
        INTEGER :: i
        REAL (kind=8),PARAMETER :: pi = 3.141592653589793
        REAL (kind=8) :: tem,temN,temH,size_dep
        REAL (kind=8),DIMENSION(5) :: z1,z2
        REAL (kind=8),PARAMETER :: tol = 1E-12
        LOGICAL :: ss_pred

        ! IF(PRESENT(steadystate_pred)) THEN
            ! ss_pred = steadystate_pred
        ! ENDIF
        IF(steadystate_pred == 1) THEN
            ss_pred = .TRUE.
        ELSE
            ss_pred = .FALSE.
        ENDIF

        !#########################################################################
        ! Compute the coefficient matrix A for the differential equation

        ! temperature annual cycle approximation
        ! te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
        ! te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
        ! te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
        ! te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

        ! DO i = 1,4 ! Average temperature dependence
        !     tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
        !     temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
        !     temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
        ! END DO

        ! Set up climate dependence factors
        tem = 0.0
        temN = 0.0
        temH = 0.0

        ! Monthly temperature dependence
        DO i = 1,12
            tem = tem+EXP(theta(22)*temp(i)+theta(23)*temp(i)**2.0)
            temN = temN+EXP(theta(24)*temp(i)+theta(25)*temp(i)**2.0)
            temH = temH+EXP(theta(26)*temp(i)+theta(27)*temp(i)**2.0)
         END DO

        ! Precipitation dependence
        tem = tem*(1.0-EXP(theta(28)*prec/1000.0))/12
        temN = temN*(1.0-EXP(theta(29)*prec/1000.0))/12
        temH = temH*(1.0-EXP(theta(30)*prec/1000.0))/12

        ! Size class dependence -- no effect if d == 0.0
        size_dep = MIN(1.0,(1.0+theta(33)*d+theta(34)*d**2.0)**(-ABS(theta(35))))

        ! check rare case where no decomposition happens for some compartments
        ! (basically, if no rain)
        IF (tem <= tol) THEN
            xt = init + b*time
            return
        END IF

        ! Calculating matrix A (will work ok despite the sign of alphas)
        DO i = 1,3
            A(i,i) = -ABS(theta(i))*tem*size_dep
        END DO
        A(4,4) = -ABS(theta(4))*temN*size_dep

        A(1,2) = theta(5)*ABS(A(2,2))
        A(1,3) = theta(6)*ABS(A(3,3))
        A(1,4) = theta(7)*ABS(A(4,4))
        A(1,5) = 0.0 ! no mass flows from H -> AWEN
        A(2,1) = theta(8)*ABS(A(1,1))
        A(2,3) = theta(9)*ABS(A(3,3))
        A(2,4) = theta(10)*ABS(A(4,4))
        A(2,5) = 0.0
        A(3,1) = theta(11)*ABS(A(1,1))
        A(3,2) = theta(12)*ABS(A(2,2))
        A(3,4) = theta(13)*ABS(A(4,4))
        A(3,5) = 0.0
        A(4,1) = theta(14)*ABS(A(1,1))
        A(4,2) = theta(15)*ABS(A(2,2))
        A(4,3) = theta(16)*ABS(A(3,3))
        A(4,5) = 0.0
        A(5,5) = -ABS(theta(32))*temH ! no size effect in humus
        DO i = 1,4
            A(5,i) = theta(31)*ABS(A(i,i)) ! mass flows AWEN -> H (size effect is present here)
        END DO

        ! Leaching (no leaching for humus)
        DO i = 1,4
            ! A(i,i) = A(i,i)+leac*climate(2)/1000.0
            A(i,i) = A(i,i)+leac*prec/1000.0
        END DO

        !#########################################################################
        ! Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init

        IF(ss_pred) THEN
            ! Solve DE directly in steady state conditions (time = infinity)
            ! using the formula 0 = x'(t) = A*x + b => x = -A**-1*b
            CALL solve(-A, b, xt)
        ELSE
            ! Solve DE in given time
            z1 = MATMUL(A,init) + b
            At = A*time !At = A*t
            CALL matrixexp(At,mexpAt)
            z2 = MATMUL(mexpAt,z1) - b
            CALL solve(A,z2,xt) ! now it can be assumed A is non-singular
        ENDIF                                   
    
    END SUBROUTINE mod5c_year



        !>
    !! @par Copyright
    !! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
    !! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
    !! file COPYING in the root of the source tree for this code.
    !! Where software is supplied by third parties, it is indicated in the headers of the routines.
    !!


    SUBROUTINE yasso_wetland (n_sp, time, init, tmp_ave,prec, Lit, parsAWEN, WoodLitterSize, Yasso_out, &
			frac_aboveground)
    !     !! DSG: 11.03.2013
    !     !! modified Yasso soil C model using an Euler scheme and readable coding style
    !     !! JSBACH specific modification of model structure: AWEN pools are separated into aboveground and belowground pools
    !     !! 
    !     !! Herbivory flux is treated like normal litter, because the herbivory fluxes in JSBACH are not evaluated, yet.
    !     !! If the represenation of herbivory in JSBACH is evaluated I recommend to treat herbivory flux in YASSO different from general
    !     !! litter.
    !     !!
    !     !! The routine needs as input 1. daily litter flux
    !     !!                            2. Yasso pools
    !     !!                            3. PFT specific parameters
    !     !!
    !     !! output is on a daily timestep: respiration
    !     !! output pools: AWEN + humus 
    !     !! The routine must be called twice, first for non-woody litter
    !     !! then for woody litter
    !     !!-------------------------------------------------------------------------------------------------------------------------
    !     !! Yasso - Soil carbon model (Jari Liski)
    !     !!
    !     !! The code for the yasso model has been developed and is provided by the Finnish Environment Institute SYKE. It is 
    !     !! distributed as part of the MPI Earth System Model under the MPI-M license conditions as stated in the header of this module. 
    !     !!
    !     !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !     !! heterotrophic soil respiration.
    !     !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !     !!
    !     !! implementation: Tea Thum (FMI), Petri Räisänen (FMI), Daniel Goll (MPI-M)
    !     !!-------------------------------------------------------------------------------------------------------------------------

        integer, intent(in):: n_sp
        integer :: i
        real (kind=8), DIMENSION(5),  intent(in)  :: init, frac_aboveground                    
        real (kind=8), INTENT(IN)  :: WoodLitterSize,time
        REAL (kind=8),INTENT(IN) :: tmp_ave, Lit 
        REAL (kind=8),INTENT(IN) :: prec 
        real (kind=8), DIMENSION(5),  intent(in)  :: parsAWEN        
        real (kind=8), DIMENSION(5), intent(out) :: Yasso_out         

        ! anoxic & peat decomposition modifier
        LOGICAL      :: anoxic_flag, peat_flag             
        
        ! Yasso pools  
        real (kind=8), dimension(n_sp)  :: YCpool_acid_ag, YCpool_water_ag, YCpool_ethanol_ag, &
                        YCpool_nonsoluble_ag, YCpool_acid_bg, YCpool_water_bg,              &
                        YCpool_ethanol_bg, YCpool_nonsoluble_bg, YCpool_humus 

        ! Yasso fluxes (internally used, therefore annual fluxes)  
        real (kind=8), dimension(n_sp) :: Cflx_litter_2_acid, Cflx_litter_2_water, Cflx_litter_2_ethanol, &
                        Cflx_litter_2_nonsoluble, Cflx_litter_2_humus

        ! Aboveground
        real (kind=8), dimension(n_sp) :: Cflx_from_acid_ag, Cflx_from_water_ag, Cflx_from_ethanol_ag, &
                        Cflx_from_nonsoluble_ag, Cflx_2_acid_ag, Cflx_2_water_ag, &
                        Cflx_2_ethanol_ag, Cflx_2_nonsoluble_ag
        ! Belowground
        real (kind=8), dimension(n_sp) :: Cflx_from_acid_bg, Cflx_from_water_bg, Cflx_from_ethanol_bg, &
                        Cflx_from_nonsoluble_bg, Cflx_from_humus, Cflx_2_acid_bg, &
                        Cflx_2_water_bg, Cflx_2_ethanol_bg, Cflx_2_nonsoluble_bg, &
                        Cflx_2_humus, Cflx_2_humusAG, Cflx_2_humusBG, YCpool_humus_in 
        
        ! Respiration is output of Yasso therefore a daily flux 
        real (kind=8), dimension(n_sp) :: soilResp_rateYasso, soilResp_rateLitterAG, soilResp_rateLitterBG
        
        ! Decomposition rates [1/a]
        real (kind=8), dimension(n_sp) :: d_acid, d_water, d_ethanol, d_nonsoluble, d_humus
        
        ! Dcomposition rate for N litter green (needed for N cycle)
        real (kind=8), dimension(n_sp) :: d_litter_green, Pseudo_litter_green

        ! Species-specific parameters
        real (kind=8), dimension(n_sp) :: frac_litter_2_acid, frac_litter_2_water, frac_litter_2_ethanol, &
                        frac_litter_2_nonsoluble, frac_litter_2_humus
        
        !Decomposition parameters
        real (kind=8), dimension(n_sp) :: d_anoxic_ag, d_anoxic_bg, &
                        d_peat_acrotelm, d_peat_catotelm
        real (kind=8) :: N_2_A, N_2_W, N_2_E, AWEN_2_H, N_2_H

        ! climatic drivers
        real (kind=8)     :: precip, temp, d_temp, d_precip, d_size, dt
            
        ! Yasso parameters
        ! decomposition rates in [1/yr] 
        real (kind=8), PARAMETER   :: ref_decomp_rate_acid       = -0.72d0   
        real (kind=8), PARAMETER   :: ref_decomp_rate_water      = -5.9d0     
        real (kind=8), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28d0    
        real (kind=8), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031d0   
        real (kind=8), PARAMETER   :: ref_decomp_rate_humus      = -0.0016d0  
        ! parameters for temperature dependence
        real (kind=8), PARAMETER   :: temp_p1                    = 0.095d0    
        real (kind=8), PARAMETER   :: temp_p2                    = -0.0014d0  
        ! parameter for precipitation dependence
        real (kind=8), PARAMETER   :: precip_p1                  = -1.21d0    
        ! parameters for size dependence 
        real (kind=8), PARAMETER   :: size_p1                    = -1.71d0    
        real (kind=8), PARAMETER   :: size_p2                    = 0.86d0     
        real (kind=8), PARAMETER   :: size_p3                    = -0.306d0   
        !Decomposition rates 
        real (kind=8), PARAMETER   :: A_2_W                      =  0.99d0
        real (kind=8), PARAMETER   :: A_2_E                      =  0.0d0
        real (kind=8), PARAMETER   :: A_2_N                      =  0.0d0
        real (kind=8), PARAMETER   :: W_2_A                      = 0.48d0
        real (kind=8), PARAMETER   :: W_2_E                      = 0.0d0
        real (kind=8), PARAMETER   :: W_2_N                      = 0.015d0
        real (kind=8), PARAMETER   :: E_2_A                      = 0.01d0
        real (kind=8), PARAMETER   :: E_2_W                      = 0.0d0
        real (kind=8), PARAMETER   :: E_2_N                      = 0.95d0
        

        dt = time/12.0d0


        N_2_A                      = 0.83d0
        N_2_W                      = 0.01d0
        N_2_E                      = 0.02d0
        AWEN_2_H                   = 0.0045d0
        N_2_H                      = 0.0045d0

        ! set up anoxic and peatland flags
        anoxic_flag = .TRUE.
        peat_flag = .TRUE.



        do i=1, n_sp

            ! set up anoxic decomposition rate modifiers
            d_anoxic_bg(i) = 1.d0
            d_anoxic_ag(i) = 1.d0
            if (anoxic_flag) then
                d_anoxic_bg = 0.35d0		
                d_anoxic_ag = 0.35d0		
            end if

            ! set up peat decomposition rate modifiers
            d_peat_acrotelm(i) = 1.d0
            d_peat_catotelm(i) = 1.d0
            if (peat_flag) then
                d_peat_acrotelm(i) = 1.d0         
                d_peat_catotelm(i) = 0.125d0		
                N_2_A                      = 0.66d0
                N_2_W                      = 0.01d0
                N_2_E                      = 0.02d0
                N_2_H                      = 0.17d0

            end if 

            frac_litter_2_acid(i)       = parsAWEN(1)
            frac_litter_2_water(i)      = parsAWEN(2)
            frac_litter_2_ethanol(i)    = parsAWEN(3)
            frac_litter_2_nonsoluble(i) = parsAWEN(4)
            frac_litter_2_humus(i)      = 0.0d0

            
            YCpool_acid_ag(i)        = init(1)
            YCpool_water_ag(i)       = init(2)
            YCpool_ethanol_ag(i)     = init(3)
            YCpool_nonsoluble_ag(i)  = init(4)
            YCpool_acid_bg(i)        = init(1)
            YCpool_water_bg(i)       = init(2)
            YCpool_ethanol_bg(i)     = init(3)
            YCpool_nonsoluble_bg(i)  = init(4)

            YCpool_humus(i)          = init(5)

            ! get the sum of the AWEN pools before decomposition
            ! This for N pools 
            Pseudo_litter_green(i) = YCpool_acid_ag(i)  &
                        + YCpool_water_ag(i)            &
                        + YCpool_ethanol_ag(i)          &
                        + YCpool_nonsoluble_ag(i)       &
                        + YCpool_acid_bg(i)             &
                        + YCpool_water_bg(i)            &
                        + YCpool_ethanol_bg(i)          &
                        + YCpool_nonsoluble_bg(i)

            ! Change units of climatic forcing variables
            precip = (prec  / 1000.d0)*dt  ! mm/month -> m/a
            temp = tmp_ave                  ! in C

            ! Calculate annual litter influxes [mol(c)/m2/a]
            Cflx_litter_2_acid(i)       = (frac_litter_2_acid(i)        *Lit)/dt
            Cflx_litter_2_water(i)      = (frac_litter_2_water(i)       *Lit)/dt
            Cflx_litter_2_ethanol(i)    = (frac_litter_2_ethanol(i)     *Lit)/dt
            Cflx_litter_2_nonsoluble(i) = (frac_litter_2_nonsoluble(i)  *Lit)/dt
            Cflx_litter_2_humus(i)      = (frac_litter_2_humus(i)       *Lit)/dt

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Calculate decomposition rates 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.d0)  
            d_size   = MIN(1.d0,(1.d0 + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.d0)**size_p3)

            ! in peatland precip modifier is inactive
            d_precip = 1.d0


            ! decomposition rates accounting for temperature, precipitation, litter size, peatland, and nutrient limitation 
            d_acid(i)       =   ref_decomp_rate_acid       *d_temp *d_precip *d_size *d_peat_acrotelm(i)
            d_water(i)      =   ref_decomp_rate_water      *d_temp *d_precip *d_size *d_peat_acrotelm(i)
            d_ethanol(i)    =   ref_decomp_rate_ethanol    *d_temp *d_precip *d_size *d_peat_acrotelm(i)
            d_nonsoluble(i) =   ref_decomp_rate_nonsoluble *d_temp *d_precip *d_size *d_peat_acrotelm(i)
            d_humus(i)      =   ref_decomp_rate_humus      *d_temp *d_precip *d_peat_catotelm(i)        

            

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculate fluxes
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! loss fluxes (negative values): 
            Cflx_from_acid_ag(i)         = YCpool_acid_ag(i)       * d_acid(i)        * d_anoxic_ag(i)
            Cflx_from_water_ag(i)        = YCpool_water_ag(i)      * d_water(i)       * d_anoxic_ag(i)           
            Cflx_from_ethanol_ag(i)      = YCpool_ethanol_ag(i)    * d_ethanol(i)     * d_anoxic_ag(i)        
            Cflx_from_nonsoluble_ag(i)   = YCpool_nonsoluble_ag(i) * d_nonsoluble(i)  * d_anoxic_ag(i)     

            Cflx_from_acid_bg(i)         = YCpool_acid_bg(i)       * d_acid(i)        * d_anoxic_bg(i)
            Cflx_from_water_bg(i)        = YCpool_water_bg(i)      * d_water(i)       * d_anoxic_bg(i)           
            Cflx_from_ethanol_bg(i)      = YCpool_ethanol_bg(i)    * d_ethanol(i)     * d_anoxic_bg(i)        
            Cflx_from_nonsoluble_bg(i)   = YCpool_nonsoluble_bg(i) * d_nonsoluble(i)  * d_anoxic_bg(i)     

            Cflx_from_humus(i)           = YCpool_humus(i)         * d_humus(i)       * d_anoxic_bg(i)  

            ! gain fluxes (positive values): 
            Cflx_2_acid_ag(i)           = ABS(Cflx_from_water_ag(i)    * W_2_A  &    
                                + Cflx_from_ethanol_ag(i)    * E_2_A  &
                                + Cflx_from_nonsoluble_ag(i) * N_2_A) 

            Cflx_2_water_ag(i)          = ABS(Cflx_from_acid_ag(i)       * A_2_W  &
                                + Cflx_from_ethanol_ag(i)    * E_2_W  &
                                + Cflx_from_nonsoluble_ag(i) * N_2_W) 

            Cflx_2_ethanol_ag(i)        = ABS(Cflx_from_acid_ag(i)       * A_2_E  &
                                + Cflx_from_water_ag(i)      * W_2_E  &
                                + Cflx_from_nonsoluble_ag(i) * N_2_E) 

            Cflx_2_nonsoluble_ag(i)     = ABS(Cflx_from_acid_ag(i)       * A_2_N  &
                                + Cflx_from_water_ag(i)      * W_2_N  &
                                + Cflx_from_ethanol_ag(i)    * E_2_N) 

            Cflx_2_acid_bg(i)           = ABS(Cflx_from_water_bg(i)      * W_2_A  &    
                                + Cflx_from_ethanol_bg(i)    * E_2_A  &
                                + Cflx_from_nonsoluble_bg(i) * N_2_A) 

            Cflx_2_water_bg(i)          = ABS(Cflx_from_acid_bg(i)       * A_2_W  &
                                + Cflx_from_ethanol_bg(i)    * E_2_W  &
                                + Cflx_from_nonsoluble_bg(i) * N_2_W) 

            Cflx_2_ethanol_bg(i)        = ABS(Cflx_from_acid_bg(i)       * A_2_E  &
                                + Cflx_from_water_bg(i)      * W_2_E  &
                                + Cflx_from_nonsoluble_bg(i) * N_2_E) 

            Cflx_2_nonsoluble_bg(i)     = ABS(Cflx_from_acid_bg(i)       * A_2_N  &
                                + Cflx_from_water_bg(i)      * W_2_N  &
                                + Cflx_from_ethanol_bg(i)    * E_2_N) 

                            
            Cflx_2_humus(i)          = ABS(Cflx_from_acid_ag(i)       &
                                + Cflx_from_water_ag(i)               &
                                + Cflx_from_ethanol_ag(i)             &
                                + Cflx_from_acid_bg(i)                &
                                + Cflx_from_water_bg(i)               &
                                + Cflx_from_ethanol_bg(i)             &
                                ) * AWEN_2_H                          &
                                + ABS(Cflx_from_nonsoluble_ag(i)      &
                                + Cflx_from_nonsoluble_bg(i)          &
                                ) * N_2_H  

            Cflx_2_humusag(i)        = ABS(Cflx_from_acid_ag(i)       &
                                + Cflx_from_water_ag(i)               &
                                + Cflx_from_ethanol_ag(i)             &
                                + Cflx_from_nonsoluble_ag(i)          &
                                ) * AWEN_2_H

            Cflx_2_humusbg(i)        = ABS(Cflx_from_acid_bg(i)       &
                                + Cflx_from_water_bg(i)               &
                                + Cflx_from_ethanol_bg(i)             &
                                + Cflx_from_nonsoluble_bg(i)          &
                                ) * AWEN_2_H

            ! the remaining fractions of the loss fluxes enter the atmosphere as respiration
            soilResp_rateYasso(i)    =  (Cflx_from_acid_ag(i) + Cflx_from_acid_bg(i))   &
                            * (1.d0 - A_2_W - A_2_E - A_2_N - AWEN_2_H)                 &
                            + (Cflx_from_water_ag(i) + Cflx_from_water_bg(i))           &      
                            * (1.d0 - W_2_A - W_2_E - W_2_N - AWEN_2_H)                 &
                            + (Cflx_from_ethanol_ag(i) + Cflx_from_ethanol_bg(i))       &     
                            * (1.d0 - E_2_A - E_2_W - E_2_N - AWEN_2_H)                 &
                            + (Cflx_from_nonsoluble_ag(i) + Cflx_from_nonsoluble_bg(i)) &  
                            * (1.d0 - N_2_A - N_2_W - N_2_E - N_2_H)                    &
                            + Cflx_from_humus(i) 

            ! litter ag & bg respiration (needed for N cycle)                       
            soilResp_rateLitterAG(i)   =  (Cflx_from_acid_ag(i) )                       &
                            * (1.d0 - A_2_W - A_2_E - A_2_N - AWEN_2_H)                 &
                            + (Cflx_from_water_ag(i) )                                  &
                            * (1.d0 - W_2_A - W_2_E - W_2_N - AWEN_2_H)                 &
                            + (Cflx_from_ethanol_ag(i) )                                &
                            * (1.d0 - E_2_A - E_2_W - E_2_N - AWEN_2_H)                 &
                            + (Cflx_from_nonsoluble_ag(i) )                             &
                            * (1.d0 - N_2_A - N_2_W - N_2_E - AWEN_2_H)         

            soilResp_rateLitterBG(i)   =  (Cflx_from_acid_bg(i) )                       &
                            * (1.d0 - A_2_W - A_2_E - A_2_N - AWEN_2_H)                 &
                            + (Cflx_from_water_bg(i) )                                  &
                            * (1.d0 - W_2_A - W_2_E - W_2_N - AWEN_2_H)                 &
                            + (Cflx_from_ethanol_bg(i) )                                &
                            * (1.d0 - E_2_A - E_2_W - E_2_N - AWEN_2_H)                 &
                            + (Cflx_from_nonsoluble_bg(i) )                             &
                            * (1.d0 - N_2_A - N_2_W - N_2_E - AWEN_2_H)
            
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Update Yasso pools
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            YCpool_acid_ag(i)           = MAX(0.d0,                                     &
                            YCpool_acid_ag(i)                                           & 
                            + (Cflx_from_acid_ag(i)                                     & 
                                + Cflx_2_acid_ag(i)                                     & 
                                + frac_aboveground(i) * Cflx_litter_2_acid(i)) *dt)         
            YCpool_water_ag(i)          = MAX(0.d0,                                     &
                            YCpool_water_ag(i)                                          &
                            + (Cflx_from_water_ag(i)                                    &
                                + Cflx_2_water_ag(i)                                    &
                                + frac_aboveground(i) * Cflx_litter_2_water(i)) *dt)
            YCpool_ethanol_ag(i)        = MAX(0.d0,                                     &
                            YCpool_ethanol_ag(i)                                        &
                            + (Cflx_from_ethanol_ag(i)                                  &
                                + Cflx_2_ethanol_ag(i)                                  &
                                + frac_aboveground(i) * Cflx_litter_2_ethanol(i)) *dt)
            YCpool_nonsoluble_ag(i)     = MAX(0.d0,                                     &
                            YCpool_nonsoluble_ag(i)                                     &
                            + (Cflx_from_nonsoluble_ag(i)                               &
                                + Cflx_2_nonsoluble_ag(i)                               &
                                + frac_aboveground(i) * Cflx_litter_2_nonsoluble(i)) *dt)


            YCpool_acid_bg(i)           = MAX(0.d0,                                     &
                            YCpool_acid_bg(i)                                           & 
                            + (Cflx_from_acid_bg(i)                                     & 
                                + Cflx_2_acid_bg(i)                                     & 
                                + (1.d0 - frac_aboveground(i)) *                           & 
                                Cflx_litter_2_acid(i)) *dt)                            
            YCpool_water_bg(i)          = MAX(0.d0,                                     & 
                            YCpool_water_bg(i)                                          &
                            + (Cflx_from_water_bg(i)                                    &
                                + Cflx_2_water_bg(i)                                    &
                                + (1.d0 - frac_aboveground(i)) *                           &
                                Cflx_litter_2_water(i)) *dt )                                
            YCpool_ethanol_bg(i)        = MAX(0.d0,                                     &
                            YCpool_ethanol_bg(i)                                        &
                            + (Cflx_from_ethanol_bg(i)                                  &
                                + Cflx_2_ethanol_bg(i)                                  &
                                + (1.d0 - frac_aboveground(i)) *                           &
                                    Cflx_litter_2_ethanol(i)) *dt)
            YCpool_nonsoluble_bg(i)     = MAX(0.d0,                                     &
                            YCpool_nonsoluble_bg(i)                                     &
                            + (Cflx_from_nonsoluble_bg(i)                               &
                                + Cflx_2_nonsoluble_bg(i)                               &
                                + (1.d0 - frac_aboveground(i)) *                           &
                                Cflx_litter_2_nonsoluble(i)) *dt)
            YCpool_humus_in(i) = YCpool_humus(i)
            YCpool_humus(i)          = MAX(0.d0,                                        &
                            YCpool_humus(i)                                             &
                            + (Cflx_from_humus(i)                                       &  
                                + Cflx_2_humus(i)                                       &
                                + Cflx_litter_2_humus(i)) *dt)

            ! compute d_litter_green; this is given by the change in AWEN pools
            if (Pseudo_litter_green(i) .GT. 0.d0) then
                d_litter_green(i) = (YCpool_acid_ag(i)   &
                            + YCpool_water_ag(i)         &
                            + YCpool_ethanol_ag(i)       &
                            + YCpool_nonsoluble_ag(i)    &
                            + YCpool_acid_bg(i)          &
                            + YCpool_water_bg(i)         &
                            + YCpool_ethanol_bg(i)       &
                            + YCpool_nonsoluble_bg(i))   &
                            / Pseudo_litter_green(i) - 1.d0
            else
                d_litter_green(i) = 0.d0
            end if 
            

            Yasso_out(1)           = YCpool_acid_bg(i) 
            Yasso_out(2)           = YCpool_water_bg(i) 
            Yasso_out(3)           = YCpool_ethanol_bg(i)
            Yasso_out(4)           = YCpool_nonsoluble_bg(i) 

        end do 

    end subroutine yasso_wetland                                                            



    !#########################################################################
    ! Functions for solving the diff. equation, adapted for the Yasso case
    SUBROUTINE matrixexp(A,B)
        IMPLICIT NONE
        ! Approximated matrix exponential using Taylor series with scaling & squaring
        ! Accurate enough for the Yasso case
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n,n),INTENT(OUT) :: B
        REAL (kind=8),DIMENSION(n,n) :: C,D
        REAL (kind=8) :: p,normiter
        INTEGER :: i,q,j
        q = 10 ! #terms in Taylor
        B = 0.d0
        DO i = 1,n
            B(i,i) = 1.d0
        END DO
        normiter = 2.d0 ! Amount of scaling & squaring
        j = 1
        CALL matrixnorm(A, p)
        DO
            IF (p<normiter) THEN
                EXIT
            END IF
            normiter = normiter*2.d0
            j = j+1
        END DO
        !write(*,*) normiter
        C = A/normiter ! scale
        B = B+C
        D = C
        DO i = 2,q ! compute Taylor expansion
            D = MATMUL(C,D)/REAL(i)
            B = B+D
        END DO
        DO i = 1,j ! square
            B = MATMUL(B,B)
        END DO
    END SUBROUTINE matrixexp

    SUBROUTINE matrixnorm(A,B)
        !returns elementwise (i.e. Frobenius) norm of a square matrix
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),INTENT(OUT) :: b
        INTEGER :: i
        b = 0.d0
        DO i = 1,n
            b = b+SUM(A(:,i)**2.d0)
        END DO
        b = SQRT(b)
    END SUBROUTINE matrixnorm


    SUBROUTINE solve(A, b, x)
        ! Solve linear system A*x = b
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n),INTENT(IN) :: b
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: x
        REAL (kind=8),DIMENSION(n,n) :: U
        REAL (kind=8),DIMENSION(n) :: c
        INTEGER :: i

        ! transform the problem to upper diagonal form
        CALL pgauss(A, b, U, c)

        ! solve U*x = c via back substitution
        x(n) = c(n)/U(n,n)
        DO i = n-1,1,-1
            x(i) = (c(i) - DOT_PRODUCT(U(i,i+1:n),x(i+1:n)))/U(i,i)
        END DO
    END SUBROUTINE solve

    SUBROUTINE pgauss(A, b, U, c)
        ! Transform the lin. system to upper diagonal form using gaussian elimination
        ! with pivoting
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n),INTENT(IN) :: b
        REAL (kind=8),DIMENSION(n,n),INTENT(OUT) :: U
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: c
        INTEGER :: k, j
        REAL,PARAMETER :: tol = 1E-12

        U = A
        c = b
        DO k = 1,n-1
            CALL pivot(U,c,k) ! do pivoting (though may not be necessary in our case)
            IF (ABS(U(k,k)) <= tol) THEN
                write(*,*) 'Warning!!! Matrix is singular to working precision!'
            END IF
            U(k+1:n,k) = U(k+1:n,k)/U(k,k)
            DO j = k+1,n
                U(j,k+1:n) = U(j,k+1:n) - U(j,k)*U(k,k+1:n)
            END DO
            c(k+1:n) = c(k+1:n) - c(k)*U(k+1:n,k)
        END DO
    END SUBROUTINE pgauss

    SUBROUTINE pivot(A, b, k)
        ! perform pivoting to matrix A and vector b at row k
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(INOUT) :: A
        REAL (kind=8),DIMENSION(n),INTENT(INOUT) :: b
        INTEGER,INTENT(IN) :: k
        INTEGER :: q, pk

        !write(*,*) 'Pivot elements are: ', A(k:n,k)
        q = MAXLOC(ABS(A(k:n,k)),1)
        !write(*,*) q
        IF (q > 1) THEN
            pk = k-1+q
            A(k:pk:pk-k,:) = A(pk:k:k-pk,:)
            b(k:pk:pk-k) = b(pk:k:k-pk)
        END IF
        !write(*,*) 'Pivot elements are: ', A(k:n,k)
    END SUBROUTINE pivot

    SUBROUTINE compAWENH(Lit,AWENH,parsAWEN)
        IMPLICIT NONE
        REAL (kind=8),DIMENSION(5),INTENT(OUT) :: AWENH
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: parsAWEN
        REAL (kind=8),INTENT(IN) :: Lit
        AWENH(1) = parsAWEN(1)*Lit
        AWENH(2) = parsAWEN(2)*Lit
        AWENH(3) = parsAWEN(3)*Lit
        AWENH(4) = parsAWEN(4)*Lit
        AWENH(5) = 0.d0
    END SUBROUTINE compAWENH

    ! Get stock changes from Yasso
    subroutine yasso_outflux (soilC, Yr_C, Yl_C, O_C, Yr_input, Yb_input, Yl_input, Yl_Coutflux, &
                            Yr_Coutflux, O_Coutflux, n_sp, nm, ii)

        integer, intent(in):: nm, n_sp, ii
        integer :: i
        real(kind=8), dimension(nm,n_sp,4,5), intent(inout) :: soilC
        real(kind=8), dimension(n_sp), intent(in) :: Yr_input, Yb_input, Yl_input
        real(kind=8), dimension(n_sp), intent(out) :: Yl_Coutflux, Yr_Coutflux, Yr_C, Yl_C
        real(kind=8), intent(out) :: O_C, O_Coutflux

        do i =1, n_sp

            Yl_Coutflux(i) = sum(soilC(ii-1,i,1,1:4)) - sum(soilC(ii,i,1,1:4)) + Yl_input(i)  
            Yr_Coutflux(i) = sum(soilC(ii-1,i,2:3,1:4)) - sum(soilC(ii,i,2:3,1:4)) + Yr_input(i) + Yb_input(i)
            
            Yr_C(i) = sum(soilC(ii,i,2:3,1:4))
            Yl_C(i) = sum(soilC(ii,i,1,1:4))

        end do

        O_Coutflux = soilC(ii-1,1,4,5) - soilC(ii,1,4,5)

        soilC(ii,1,4,5) = soilC(ii,1,4,5) + sum(soilC(ii,:,1:3,5))
        soilC(ii,:,1:3,5) = 0.d0

        O_C = soilC(ii,1,4,5)

    end subroutine  yasso_outflux 

        ! Get stock changes from Yasso at yearly time step
    subroutine yasso_outflux_year (soilC, Yr_C, Yl_C, O_C, Yr_input, Yb_input, Yl_input, Yl_Coutflux, &
                                    Yr_Coutflux, O_Coutflux, n_sp, nm, ii, step_aux)

        integer, intent(in):: nm, n_sp, ii, step_aux
        integer :: i
        real(kind=8), dimension(nm,n_sp,4,5), intent(inout) :: soilC
        real(kind=8), dimension(n_sp), intent(in) :: Yr_input, Yb_input, Yl_input
        real(kind=8), dimension(n_sp), intent(out) :: Yl_Coutflux, Yr_Coutflux, Yr_C, Yl_C
        real(kind=8), intent(out) :: O_C, O_Coutflux

        do i =1, n_sp

            Yl_Coutflux(i) = sum(soilC(ii-step_aux,i,1,1:4)) - sum(soilC(ii,i,1,1:4)) + Yl_input(i)  
            Yr_Coutflux(i) = sum(soilC(ii-step_aux,i,2:3,1:4)) - sum(soilC(ii,i,2:3,1:4)) + Yr_input(i) + Yb_input(i)

            Yr_C(i) = sum(soilC(ii,i,2:3,1:4))
            Yl_C(i) = sum(soilC(ii,i,1,1:4))

        end do

        O_Coutflux = soilC(ii-step_aux,1,4,5) - soilC(ii,1,4,5)

        soilC(ii,1,4,5) = soilC(ii,1,4,5) + sum(soilC(ii,:,1:3,5))
        soilC(ii,:,1:3,5) = 0.d0

        O_C = soilC(ii,1,4,5)

    end subroutine  yasso_outflux_year 

    ! Get available nitrogen based on stoichiometry pof carbon decomposition
    subroutine yasso_nav (Yr_C, Yl_C, O_C, Yl_Coutflux, Yr_Coutflux, O_Coutflux, humification_N_l, humification_N_r, & 
                        Yr_Noutflux, Yl_Noutflux, O_Noutflux, Yr_N, Yl_N, O_N, hc, qh, el, qbc, qir, qil, er, &
                        Yr_input, Yl_input, Yb_input, soil_class, excessSW, asw, n_depo, Nav, n_sp)
    
        integer, intent(in):: n_sp, soil_class
        integer :: i
        real(kind=8), dimension(n_sp), intent(in) :: hc, qh, el, qbc, er, Yr_input, Yl_input, Yb_input, qir, qil, &
                                                    Yl_Coutflux, Yr_Coutflux,  Yr_C, Yl_C
        real(kind=8), intent(in) :: excessSW, asw, n_depo, O_C, O_Coutflux
        real(kind=8), dimension(n_sp), intent(out) :: Yl_Noutflux, Yr_Noutflux, humification_N_l, &
                                                        humification_N_r, Yr_N, Yl_N
        real(kind=8), intent(out) :: Nav,  O_N,  O_Noutflux
        real(kind=8) :: leaching_min, leaching_org

        do i =1, n_sp
            
            humification_N_l(i) = Yl_Coutflux(i)/(Yl_C(i) + 0.01d0) * hc(i) * (Yl_N(i) / qh(i))
            Yl_Noutflux(i) = max((Yl_Coutflux(i)/(Yl_C(i) + 0.01d0)) * (1-hc(i)) / (1 - el(i)) * &
                    (Yl_N(i) - el(i) * (Yl_C(i) / qbc(i))),0.d0)

            humification_N_r(i) = Yr_Coutflux(i)/(Yr_C(i) + 0.01d0) * hc(i) * (Yr_N(i) / qh(i))
            Yr_Noutflux(i) = max((Yr_Coutflux(i)/(Yr_C(i) + 0.01d0)) * (1-hc(i)) / (1 - er(i)) * &
                        (Yr_N(i) - er(i) * (Yr_C(i) / qbc(i))),0.d0)

            !Now calculate the end-of-month carbon and nitrogen pools					
            Yr_N(i) = Yr_N(i) + ((Yr_input(i) + Yb_input(i)) / (2.d0 * qir(i))) - Yr_Noutflux(i) - humification_N_r(i)
            Yl_N(i) = Yl_N(i) + ((Yl_input(i)) / (2.d0 * qil(i))) - Yl_Noutflux(i) - humification_N_l(i)
    
        end do

        O_Noutflux = (O_Coutflux/O_C) * O_N

        O_N = O_N - O_Noutflux + sum(humification_N_r(:)) + sum(humification_N_l(:))
                
        call compute_leaching( soil_class, leaching_org, leaching_min, excessSW, asw )
    
        Nav = sum(Yr_Noutflux(:)) * (1.d0-leaching_org)  + sum(+ Yl_Noutflux(:)) * (1.d0-leaching_org) + &
                    O_Noutflux * (1.d0 - leaching_min)  + n_depo		
        
    end subroutine  yasso_nav 
    ! End subroutines for Yasso -----------------------------------------------------------


end module soil_cnp