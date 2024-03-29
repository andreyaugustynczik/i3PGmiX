library(r3PG)
library(testthat)

context("Basic Model runs work")

d_site$to <- "2050-12"
d_thinning$BA[1] <- 1
d_thinning$stems_n[1] <- 0
d_thinning$Vol[1] <- -999

test_that("basic model run", {
    out <- run_3PG(
      site = d_site,
      species = d_species[1,],
      climate = d_climate[1:12,],
      thinning = d_thinning,
      parameters = d_parameters,
      size_dist = d_sizeDist,
      settings = list(light_model = 2, transp_model = 1, phys_model = 1,soil_model=2, method_jmaxlim=3,nut_lim=0,
                      wind_dist=1, beetle_dist=1,dist_start=5*12,
        correct_bias = 1, calculate_d13c = 0, canopy_cond=1, gpp_model=2, maint_resp=2),
      check_input = TRUE, df_out = FALSE)

    testthat::expect_true(class(out) == "array")
})


