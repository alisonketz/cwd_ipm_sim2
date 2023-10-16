############################################################################################
############################################################################################
############################################################################################
###
### Model Statement
###
############################################################################################
############################################################################################
############################################################################################

modelcode <- nimbleCode({

  ##############################
  ### Priors
  ##############################

  ##############################
  ### Force of infection model
  ##############################

  tau_age_foi  ~ dgamma(1, 1)
  tau1_age_foi <- .0000001 * tau_age_foi
  foi_age_effect[1] ~ dnorm(0, tau1_age_foi)
  foi_age_effect[2] ~ dnorm(0, tau1_age_foi)
  for (i in 3:n_ageclass) {
    foi_age_effect[i]~dnorm(2 * foi_age_effect[i-1] - foi_age_effect[i-2], tau_age_foi)
  }
  foi_age_mu <- mean(foi_age_effect[1:n_ageclass])

  ############################################################
  ############################################################
  ### Age/Period Survival Model
  ############################################################
  ############################################################

  ####################################
  ### Susceptibles survival intercept
  ####################################

  beta0_sus_temp ~ dnorm(0, .1)
  sus_mix ~ dunif(-1, 1)
  beta0_survival_sus <- beta0_sus_temp * sus_mix

  ##################################
  ### Infected survival intercept
  ##################################

  beta0_inf_temp ~ dnorm(0, .1)
  inf_mix ~ dunif(-1, 1)
  beta0_survival_inf <- beta0_inf_temp * inf_mix

  ########################################
  ### Priors for Age Effects Survival
  ########################################

  ### Age effects
  for (k in 1:nknots_age) {
    ln_b_age_survival[k] ~ dnorm(0, tau_age_survival)
    b_age_survival[k] <- exp(ln_b_age_survival[k])
  }
  tau_age_survival ~ dgamma(1, 1)

  for (t in 1:nT_age) {
    age_effect_survival_temp[t] <- inprod(b_age_survival[1:nknots_age],
                                     Z_age[t, 1:nknots_age])
  }
  mu_age_effect_survival_temp <- mean(age_effect_survival_temp[1:nT_age])
  
  for (t in 1:nT_age) {
    age_effect_survival[t] <-  age_effect_survival_temp[t] -
                               mu_age_effect_survival_temp
  }

  ########################################
  ### Priors for Period Effects Survival
  ########################################

  mix_survival ~ dunif(-1, 1)
  ln_sk_period ~ dnorm(0, sd = 1)
  sdk_period <- exp(mix_survival * ln_sk_period)
  tauk_period <- 1 / sdk_period^2
  stauk_period <- sqrt(tauk_period)
  sda_period ~ T(dnorm(0, sd = 1), 0, Inf)#<- 1/sqrt(taua_period)
  taua_period <- 1 / sda_period^2
  for (i in 1:(nknots_period)) {
    alpha_period[i] ~ dnorm(0, 1)
    alphau_period[i] <- sda_period * alpha_period[i]
  }
  ratioinf_period <- sdk_period / sda_period #ratio of variability

  period_effect_surv[1:nT_period] <- kernel_conv(
    nT = nT_period,
    Z = Z_period[1:nT_period, 1:nknots_period],
    stauk = stauk_period,
    nconst = nconst,
    tauk = tauk_period,
    nknots = nknots_period,
    alphau = alphau_period[1:nknots_period]
  )


  #######################################################################
  #######################################################################
  ## Likelihoods of Joint Model
  #######################################################################
  #######################################################################


  #######################################################################
  ###
  ###   User defined distribution for likelihood for
  ###   infected deer mortalities for radio marked deer that
  ###   enter the study as test positive at capture
  ###
  ###   d_fit_icap_mort
  ###
  ###   Overleaf Equation (12)
  ###
  #######################################################################

    y_icap_mort ~ dIcapMort(
      n_samples = n_fit_icap_mort,
      e = icap_mort_e_age[1:n_fit_icap_mort],
      r = icap_mort_r_age[1:n_fit_icap_mort],
      s = icap_mort_s_age[1:n_fit_icap_mort],
      age2date = icap_mort_age2date[1:n_fit_icap_mort],
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age],
      period_effect_surv = period_effect_surv[1:nT_period]
      )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   infected deer mortalities for radio marked deer that
    ###   enter the study as test positive at capture
    ###
    ###   d_fit_icap_cens
    ###
    ###   Overleaf Equation (10)
    ###
    #######################################################################

    y_icap_cens ~ dIcapCens(
      n_samples = n_fit_icap_cens,
      e = icap_cens_e_age[1:n_fit_icap_cens],
      r = icap_cens_r_age[1:n_fit_icap_cens],
      age2date = icap_cens_age2date[1:n_fit_icap_cens],
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age],
      period_effect_surv = period_effect_surv[1:nT_period]
      )

    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   infected deer mortalities for radio marked deer that
    ###   enter the study as test negative at capture
    ###
    ###   d_fit_idead
    ###
    ###   Overleaf Equation (24)
    ###
    #######################################################################

    y_idead ~ dNegCapPosMort(
        n_samples = n_fit_idead,
        e = idead_e_age[1:n_fit_idead],
        r = idead_r_age[1:n_fit_idead],
        s = idead_s_age[1:n_fit_idead],
        age2date = idead_age2date[1:n_fit_idead],
        beta0_sus = beta0_survival_sus,
        beta0_inf = beta0_survival_inf,
        age_effect_surv = age_effect_survival[1:nT_age],
        period_effect_surv = period_effect_surv[1:nT_period],
        foi_age_effect = foi_age_effect[1:n_ageclass],
        age_lookup = age_lookup[1:nT_age]
        )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   uninfected radio-marked deer mortalities:
    ###   test neg at cap and tested mort
    ###
    ###   d_fit_sus_mort_posttest
    ###
    ###   Overleaf Equation (11)
    ###
    #######################################################################

    y_sus_mort_posttest ~ dSusMortTest(
        n_samples = n_fit_sus_mort_posttest,
        e = sus_mort_posttest_e_age[1:n_fit_sus_mort_posttest],
        r = sus_mort_posttest_r_age[1:n_fit_sus_mort_posttest],
        s = sus_mort_posttest_s_age[1:n_fit_sus_mort_posttest],
        age2date = sus_mort_posttest_age2date[1:n_fit_sus_mort_posttest],
        beta0_sus = beta0_survival_sus,
        age_effect_surv = age_effect_survival[1:nT_age],
        period_effect_surv = period_effect_surv[1:nT_period],
        foi_age_effect =  foi_age_effect[1:n_ageclass],
        age_lookup = age_lookup[1:nT_age]
        )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   Uninfected radio-marked deer right censored:
    ###   Test neg at cap and censoring
    ###
    ###   d_fit_sus_cens_postno
    ###
    ###   Overleaf Equation (9)
    ###
    #######################################################################

    y_sus_cens_postno ~ dSusCensNo(
          n_samples = n_fit_sus_cens_postno,
          e = sus_cens_postno_e_age[1:n_fit_sus_cens_postno],
          r = sus_cens_postno_r_age[1:n_fit_sus_cens_postno],
          age2date = sus_cens_postno_age2date[1:n_fit_sus_cens_postno],
          beta0_sus = beta0_survival_sus,
          beta0_inf = beta0_survival_inf,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )

    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   Uninfected radio-marked deer right censored:
    ###   Test neg at cap and censoring
    ###
    ###   d_fit_sus_cens_posttest
    ###   Overleaf Equation 2
    ###
    ###
    #######################################################################

    y_sus_cens_posttest ~ dSusCensTest(
          n_samples = n_fit_sus_cens_posttest,
          e = sus_cens_posttest_e_age[1:n_fit_sus_cens_posttest],
          r = sus_cens_posttest_r_age[1:n_fit_sus_cens_posttest],
          age2date = sus_cens_posttest_age2date[1:n_fit_sus_cens_posttest],
          beta0_sus = beta0_survival_sus,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   uninfected radio-marked deer mortalities:
    ###   test neg at cap and no test at mortality, no recap
    ###
    ###   d_fit_sus_mort_postno
    ###
    ###   Overleaf Equation (8)
    ###
    #######################################################################

    y_sus_mort_postno ~ dSusMortNoTest(
          n_samples = n_fit_sus_mort_postno,
          e = sus_mort_postno_e_age[1:n_fit_sus_mort_postno],
          r = sus_mort_postno_r_age[1:n_fit_sus_mort_postno],
          s = sus_mort_postno_s_age[1:n_fit_sus_mort_postno],
          age2date = sus_mort_postno_age2date[1:n_fit_sus_mort_postno],
          beta0_sus = beta0_survival_sus,
          beta0_inf = beta0_survival_inf,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )

    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   uninfected deer that were test neg at capture,
    ###   then test negative at recap, that are right censored,
    ###   and have been tested post censoring
    ###
    ###   d_fit_rec_neg_cens_posttest
    ###
    ###   Overleaf Equation (14)
    ###
    #######################################################################

    y_rec_neg_cens_posttest ~ dRecNegCensTest(
          n_samples = n_fit_rec_neg_cens_posttest,
          e = rec_neg_cens_posttest_e_age[1:n_fit_rec_neg_cens_posttest],
          r = rec_neg_cens_posttest_r_age[1:n_fit_rec_neg_cens_posttest],
          age2date = rec_neg_cens_posttest_age2date[1:n_fit_rec_neg_cens_posttest],
          beta0_sus = beta0_survival_sus,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   uninfected deer that were test neg at capture,
    ###   then test negative at recap, that are right censored,
    ###   and have not been tested post censoring
    ###
    ###   d_fit_rec_neg_cens_postno
    ###
    ###   Overleaf Equation (16)
    ###
    #######################################################################

    y_rec_neg_cens_postno ~ dRecNegCensPostNo(
          n_samples = n_fit_rec_neg_cens_postno,
          e = rec_neg_cens_postno_e_age[1:n_fit_rec_neg_cens_postno],
          r = rec_neg_cens_postno_r_age[1:n_fit_rec_neg_cens_postno],
          dn1 = rec_neg_cens_postno_rec_age[1:n_fit_rec_neg_cens_postno],
          age2date = rec_neg_cens_postno_age2date[1:n_fit_rec_neg_cens_postno],
          beta0_sus = beta0_survival_sus,
          beta0_inf = beta0_survival_inf,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   uninfected deer that were test neg at capture,
    ###   then test negative at recap,
    ###   that die
    ###
    ###   d_fit_rec_neg_mort
    ###
    ###   Overleaf Equation (18)
    ###
    #######################################################################

    y_rec_neg_mort ~ dRecNegMort(
          n_samples = n_fit_rec_neg_mort,
          e = rec_neg_mort_e_age[1:n_fit_rec_neg_mort],
          r = rec_neg_mort_r_age[1:n_fit_rec_neg_mort],
          s = rec_neg_mort_s_age[1:n_fit_rec_neg_mort],
          age2date = rec_neg_mort_age2date[1:n_fit_rec_neg_mort],
          beta0_sus = beta0_survival_sus,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )

    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   infected deer that were test neg at capture,
    ###   then test positive at recap,
    ###   than these were right censored
    ###
    ###   d_fit_rec_pos_cens
    ###
    ###   Overleaf Equation (20)
    ###
    #######################################################################

    y_rec_pos_cens ~ dRecPosCens(n_samples = n_fit_rec_pos_cens,
          e = rec_pos_cens_e_age[1:n_fit_rec_pos_cens],
          r = rec_pos_cens_r_age[1:n_fit_rec_pos_cens],
          dn = rec_pos_cens_rec_age[1:n_fit_rec_pos_cens],
          age2date = rec_pos_cens_age2date[1:n_fit_rec_pos_cens],
          beta0_sus = beta0_survival_sus,
          beta0_inf = beta0_survival_inf,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )

    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   deer that were test neg at capture,
    ###   then test positive at recap,
    ###   than die
    ###
    ###   d_fit_rec_pos_mort
    ###
    ###   Overleaf Equation (22)
    ###
    #######################################################################

    y_rec_pos_mort ~ dRecPosMort(
          n_samples = n_fit_rec_pos_mort,
          e = rec_pos_mort_e_age[1:n_fit_rec_pos_mort],
          r = rec_pos_mort_r_age[1:n_fit_rec_pos_mort],
          s = rec_pos_mort_s_age[1:n_fit_rec_pos_mort],
          dn = rec_pos_mort_rec_age[1:n_fit_rec_pos_mort],
          age2date = rec_pos_mort_age2date[1:n_fit_rec_pos_mort],
          beta0_sus = beta0_survival_sus,
          beta0_inf = beta0_survival_inf,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )

    #######################################################################
    #######################################################################
    #######################################################################
    ###
    ###   Derived parameters: 
    ###   Annual survival estimates "May 15 - May 14"
    ###   (not using hazards from Jan - May of the first year)
    ###
    ###
    #######################################################################
    #######################################################################
    #######################################################################

    sn_sus[1:n_age, 1:n_year]  <- calc_surv_aah(nT_age = nT_age,
        nT_period = nT_period,
        nT_age_short = nT_age_short,
        nT_age_surv_aah = nT_age_surv_aah,
        beta0 = beta0_survival_sus,
        age_effect = age_effect_survival[1:nT_age],
        period_effect = period_effect_surv[1:nT_period],
        yr_start_age = yr_start_age[1:n_yr_start_age],
        yr_start_pop = yr_start_pop[1:n_year],
        n_year = n_year,
        n_age = n_age)

    sn_inf[1:n_age, 1:n_year]  <- calc_surv_aah(nT_age = nT_age,
        nT_period = nT_period,
        nT_age_short = nT_age_short,
        nT_age_surv_aah = nT_age_surv_aah,
        beta0 = beta0_survival_inf,
        age_effect = age_effect_survival[1:nT_age],
        period_effect = period_effect_surv[1:nT_period], 
        yr_start_age = yr_start_age[1:n_yr_start_age],
        yr_start_pop = yr_start_pop[1:n_year],
        n_year = n_year,
        n_age = n_age)

    psi[1:n_age, 1:n_year] <- calc_infect_prob(age_lookup = age_lookup[1:nT_age],
                        n_age = n_age,
                        yr_start = yr_start_age[1:n_yr_start_age],
                        age_foi = foi_age_effect[1:n_ageclass],
                        nT_period = nT_period,
                        n_year = n_year,
                        nT_age_surv_aah = nT_age_surv_aah
                        )

})#end model statement
