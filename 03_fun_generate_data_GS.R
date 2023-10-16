 ageperiod_surv_foi_sim_data <- function(
    beta0_survival_sus,
    beta0_survival_inf,
    foi_age_effect,
    age_effect,
    period_effect,
    nT_age,
    nT_period,
    prop_recap = 0.03
    ) {
  
  set.seed(12345)
  
  # sample size of individuals from each year
  n_years <- 5
  n_per_yr <- 250
  n_yr <- rep(n_per_yr, n_years)
  n_ind <- sum(n_yr)

  ### Intercepts for mortality
  beta0_sus <- beta0_survival_sus
  beta0_inf <- beta0_survival_inf

  ### maximum age interval
  nT_age <- nT_age

  ### maximum period interval
  nT_period <- nT_period

  ########################################
  # Mortality hazards
  ########################################

  age_effect_surv <- age_effect
  period_effect_surv <- period_effect

  ########################################
  # Age effects for the mortality hazard
  ########################################
  # age <- seq(1, nT_age - 1, by = 1)
  # age_effect <- age_effect
  # age_effect_surv <- age_effect
  
  ########################################
  # Age effects for the foi log hazard
  ########################################

  # logistic
  # age_effect_foi <- 1/(1+(1/.01-1)*exp(-0.7/52*age))
  # age_effect_foi <- age_effect_foi - mean(age_effect_foi)

  # png("figures/age_effects_foi_generating.png")
  # plot(age_effect_foi)
  # dev.off()


  #####################################################
  ### Age at entry with staggered entry,
  ### drawn from a stable age distribution
  ### based on the age effect hazard (but considering only susceptibles)
  ### for n_years 'years', and for the age that
  ### individuals that go from birth to censoring
  #####################################################

  hazard_se <- -exp((beta0_sus + age_effect_surv)) # not quite right because it only accounts for sus survival
  stat_se <- rep(NA, nT_age - 1)
  stat_se[1] <- (1 - exp(hazard_se[1]))
  for (j in 2:(nT_age - 1)) {
    stat_se[j] <- (1 - exp(hazard_se[j])) * exp(sum(hazard_se[1:(j - 1)]))
  }
  left_age <- nimble::rcat(n_ind, stat_se)
  maxages <- nT_age - left_age

  ###########################################
  # Staggered entry period effects
  ###########################################

  weeks_entry <- 12
  left_yr <- matrix(NA,n_years,n_per_yr)
  for(i in 1:n_years){
	left_yr[i,] <- nimble::rcat(n_per_yr,
          prob = rep(1 / weeks_entry, weeks_entry)) + (i - 1) * 52
  }
  left_period <- c(t(left_yr))

  ### dealing with maximum possible period individuals can
  ### be alive and in the study
  maxtimes <- nT_period - left_period
  maxtimes <- ifelse(maxtimes <= maxages, maxtimes, maxages)
 
  ########################################
  ### Period effects for the survival log hazard 
  ########################################

  # period <- seq(1, nT_period - 1, by = 1)
  # period_effect_surv <- 1 * sin(2/52 * pi * (period) + 1)
  # period_effect_surv <- period_effect_surv - mean(period_effect_surv)

  # png("figures/period_effects_survival_generating.png")
  # plot(period_effect_surv)
  # dev.off()

  ########################################
  ### Force of infection, log hazards 
  ########################################


  ########################################
  ### Period effects for the foi log hazard 
  ########################################

  # logistic
  # period_effect_foi <- 1/(1+(1/.001-1)*exp(-0.5/52*period))
  # period_effect_foi <- period_effect_foi - mean(period_effect_foi)

  # png("figures/period_effects_foi_generating.png")
  # plot(period_effect_foi)
  # dev.off()
  ########################################
  ### Log hazards 
  ########################################
  # account for lesser foi age effects for individuals older at capture 
  foi_age_discount <- seq(2.5,1.25,length.out=nT_age)

  hazard_sus <- hazard_inf <- hazard_foi <- matrix(0, n_ind, nT_age)
  for (i in 1:n_ind) {
    hazard_sus[i,left_age[i]:(left_age[i] + maxtimes[i] - 1)] <-  exp(beta0_sus+
                      age_effect_surv[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_surv[left_period[i]:(left_period[i] + maxtimes[i] - 1)])

    hazard_inf[i,left_age[i]:(left_age[i] + maxtimes[i] - 1)] <- exp(beta0_inf+
                      age_effect_surv[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_surv[left_period[i]:(left_period[i] + maxtimes[i] - 1)])
	
	# account for lower foi hazard prior to study period
	foi_age_effect_discount <- foi_age_effect
	foi_age_effect_discount[1:left_age[i]] <- foi_age_effect[1:left_age[i]]*foi_age_discount[(nT_age-left_age[i]+1):nT_age]
	
	hazard_foi[i,1:(left_age[i] + maxtimes[i] - 1)] <- exp(
                      foi_age_effect_discount[1:(left_age[i] + maxtimes[i] - 1)])
  }


  ########################################
  ### Probability of infection
  ########################################

 # incidence <- matrix(0, n_ind, nT_age)
 # for (i in 1:n_ind) {
 #   for (j in left_age[i]:(left_age[i] + maxages[i] - 1)) {
 #     if (j==left_age[i]) {
 #       incidence[i, j] <- (1 - exp(-exp(hazard_foi[i, j])))
 #     } else {
 #       incidence[i, j] <-
 #         (1 - exp(-exp(hazard_foi[i, j])))*exp(-sum(exp(hazard_foi[i, left_age[i]:(j - 1)])))
 #     }
 #   }
 #   incidence[i, left_age[i] + maxages[i]] <-
 #         exp(-sum(exp(hazard_foi[i, left_age[i]:(left_age[i] + maxages[i] - 1)])))
 # }

  incidence <- matrix(0, n_ind, nT_age)
  for (i in 1:n_ind) {
    for (j in 1:(left_age[i] + maxtimes[i] - 1)) {
        incidence[i, j] <- (1 - exp(-hazard_foi[i, j]))
      } 
  }


  ########################################
  ### Calculating infection status
  ########################################
  inf_age <- rep(NA, n_ind)
  inf_status <- matrix(0, n_ind, nT_age)
  for(i in 1:n_ind) {
    inf_status[i,] <- rbinom(nT_age,1,incidence[i,1:nT_age])
    if(max(inf_status[i,]) == 0) next
      inf_age[i] <- min(which(inf_status[i,] == 1))
      inf_status[i,inf_age[i]:nT_age] <- 1
  }

  ########################################
  ### Probability of mortality
  ########################################

  survival_prob <- matrix(0, n_ind, nT_age)
  for (i in 1:n_ind) {
    for (j in left_age[i]:(left_age[i] + maxtimes[i] - 1)) {
      if (j==left_age[i]) {
        survival_prob[i, j] <- 
          1 - ((1 - inf_status[i, j]) * exp(-hazard_sus[i, j]) + 
				  inf_status[i, j] * exp(-hazard_inf[i, j]))
      } else {
	      survival_prob[i, j] <- 
          (1 - ((1 - inf_status[i, j])*exp(-hazard_sus[i, j]) + 
				  inf_status[i, j] * exp(-hazard_inf[i, j]))) *
			    exp(-(sum((1 - inf_status[i, left_age[i]:(j - 1)]) * hazard_sus[i, left_age[i]:(j - 1)]) +
				      sum(inf_status[i, left_age[i]:(j - 1)] * hazard_inf[i, left_age[i]:(j - 1)])))
	    }
    }
    survival_prob[i, left_age[i] + maxtimes[i]] <-
		 (1 - inf_status[i, left_age[i] + maxtimes[i]])*exp(-sum(hazard_sus[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)])) + 
			  inf_status[i, left_age[i] + maxtimes[i]]*exp(-sum(hazard_inf[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)]))
  }


  ########################################
  ### Calculating right censoring
  ########################################
  p_censor <- 0.0002
  fail_int <- rep(0, n_ind)
  right_age <- rep(0, n_ind)
  rt_censor <- rep(0, n_ind)
  collar_loss_int <- rep(NA, n_ind)
  test1 <- left_age
  pos1 <- inf_during_study <- rep(0, n_ind)
  pos2 <- pos3 <- test2 <- test3 <- rec_age <- data_type <- rep(NA, n_ind)
 
  for(i in 1:n_ind) {
    fail_int[i] <- which(rmultinom(1, 1, survival_prob[i,1:nT_age]) == 1)
    right_age[i] <- ifelse(fail_int[i] >= left_age[i] + maxtimes[i],
                           left_age[i] + maxtimes[i],
                           fail_int[i] + 1)
	#loss <- rbinom(length(left_age[i]:right_age[i]),1,p_censor)
	#if(max(loss) == 1){
	#	collar_loss_int[i] <- min(which(loss == 1)) + left_age[i]
	#	right_age[i] <- min(collar_loss_int[i],right_age[i])
	#}
    #rt_censor[i] <- ifelse(fail_int[i] >= (left_age[i] + maxtimes[i]) | is.na(collar_loss_int[i]), 1, 0)
    rt_censor[i] <- ifelse(fail_int[i] < (left_age[i] + maxtimes[i]), 0, 1)
    if(!is.na(inf_age[i])){
      if(inf_age[i] >= right_age[i]){
        inf_age[i] <- NA
      } else {
        if(inf_age[i] < left_age[i]) {
          pos1[i] <- 1
          }
      }
    }
  }
  inf_during_study[which(inf_age>left_age & inf_age<right_age)] <- 1
  #(quick <- which(right_age-left_age < 2))
  #(right_age-left_age)[quick]
  

	########################################
	# Deal with those that got infected during study period
	# could be idead, rec_pos_mort, rec_pos_cens, sus_cens_postno, 
	# endlive, sus_mort_post_no, rec_neg_cens_postno 
	# infected dead seem over-represented, and infected censored under-represented
	########################################
	cand_dead <- which(pos1 == 0 & inf_during_study == 1 & rt_censor == 0)
	cand_cens <- which(pos1 == 0 & inf_during_study == 1 & rt_censor == 1)
	# allocate dead ones (most go in the "known infected" bins)
	prob <- c(7,8,108)/sum(c(7,8,108))
	adj_prob <- cumsum(rev(prob))
	inf_dead_class <- sample(rep(3:1, diff(floor(length(cand_dead) * c(0,adj_prob)))))
	for(i in 1:length(cand_dead)){
		j <- cand_dead[i]
		if(inf_dead_class[i] == 3){ 
			data_type[j] <- "idead"
			test2[j] <- right_age[j]
			pos2[j] = 1
		}
 		if(inf_dead_class[i] == 2){ 
			data_type[j] <- "rec_pos_mort"
			rec_age[j] <- sample(inf_age[j]:(right_age[j]-1),1)
			test2[j] <- rec_age[j]
			pos2[j] = 1
		}
 		if(inf_dead_class[i] == 1){ 
			data_type[j] <- "sus_mort_postno"
		}
   }
	# allocate censored ones (sus_cens_postno and endlive are the same here)
	# Make sure to put one in the rec_pos_cens group
	prob <- c(353,12)/sum(c(353,12))
	adj_prob <- cumsum(prob)
	inf_cens_class <- c(3,3,sample(rep(1:2, diff(floor(length(cand_cens[-c(1:2)]) * c(0,adj_prob))))))
	for(i in 1:length(cand_cens)){
		j <- cand_cens[i]
 		if(inf_cens_class[i] == 1){ 
			data_type[j] <- "sus_cens_postno"
		}
 		if(inf_cens_class[i] == 2){ 
			data_type[j] <- "rec_neg_cens_postno"
			rec_age[j] <- sample(left_age[j]:(inf_age[j]-1),1)
			test2[j] <- rec_age[cand_cens[j]]
			pos2[j] = 0
		}
 		if(inf_cens_class[i] == 3){ 
			data_type[j] <- "rec_pos_cens"
			rec_age[j] <- sample(inf_age[j]:(right_age[j]-1),1)
			test2[j] <- rec_age[j]
			pos2[j] = 1
		}
   }

	########################################
	# Deal with those that not uninfected
	# some died, some censored, some posttested, some not posttest
	# more dead than expected, and fewer censored
	########################################
	cand_dead <- which(pos1 == 0 & inf_during_study == 0 & rt_censor == 0)
	cand_cens <- which(pos1 == 0 & inf_during_study == 0 & rt_censor == 1)
	# allocate dead ones 
	prob <- c(428,7,5)/sum(c(428,7,5))
	adj_prob <- cumsum(prob)
	sus_dead_class <- sample(rep(1:3, diff(floor(length(cand_dead) * c(0,adj_prob)))))
	for(i in 1:length(cand_dead)){
		j <- cand_dead[i]
		if(sus_dead_class[i] == 1){ 
			data_type[j] <- "sus_mort_posttest"
			test2[j] <- right_age[j]
			pos2[j] = 0
		}
 		if(sus_dead_class[i] == 2){ 
			data_type[j] <- "sus_mort_postno"
		}
 		if(sus_dead_class[i] == 3){ 
			data_type[j] <- "rec_neg_mort"
			rec_age[j] <- sample(left_age[j]:(right_age[j]-1),1)
			test2[j] <- rec_age[j]
			pos2[j] = 0
			test3[j] <- right_age[j]
			pos3[j] = 0
		}
   }
	# allocate censored ones (sus_cens_postno and endlive are the same here)
	prob <- c(353,66,12)/sum(c(353,66,12))
	adj_prob <- cumsum(prob)
	sus_cens_class <- c(3,3,sample(rep(c(2,1,4), diff(floor(length(cand_cens[-c(1:2)]) * c(0,adj_prob))))))
	for(i in 1:length(cand_cens)){
		j <- cand_cens[i]
		if(sus_cens_class[i] == 1){ 
			data_type[j] <- "sus_cens_posttest"
			test2[j] <- right_age[j]
			pos2[j] = 0
		}
 		if(sus_cens_class[i] == 2){ 
			data_type[j] <- "sus_cens_postno"
		}
 		if(sus_cens_class[i] == 3){ 
			data_type[j] <- "rec_neg_cens_posttest"
			rec_age[j] <- sample(left_age[j]:(right_age[j]-1),1)
			test2[j] <- rec_age[j]
			pos2[j] = 0
			test3[j] <- right_age[j]
			pos3[j] = 0
		}
 		if(sus_cens_class[i] == 4){ 
			data_type[j] <- "rec_neg_cens_postno"
			rec_age[j] <- sample(left_age[j]:(right_age[j]-1),1)
			test2[j] <- rec_age[j]
			pos2[j] = 0
		}
	}

	########################################
	# Deal with those infected at capture 
	# Based on real data, these are way over-represented
	########################################
	cand_dead <- which(pos1 == 1 & rt_censor == 0)
	cand_cens <- which(pos1 == 1 & rt_censor == 1)
	# dead ones 
	for(i in cand_dead){
		data_type[i] <- "icap_mort"
	}
	for(i in cand_cens){
		data_type[i] <- "icap_cens"
	}
  
table(data_type)
  ########################################
  ### setup data into overall data frame
  ########################################

  right_period <- right_age - left_age + left_period
  cwd_mort <- pos1*(1-rt_censor) + pos2*(1-rt_censor)
  cwd_mort[is.na(cwd_mort)] <- 0
  
  df_fit <- data.frame(id = 1:n_ind,
                      left_age = left_age,
                      right_age  = right_age,
                      left_period = left_period,
                      right_period = right_period,
                      rt_censor = rt_censor,
                      cwd_cap = pos1,
					            cwd_mort = cwd_mort,
                      inf_age = inf_age,
					  inf_during = inf_during_study,
					  test2_age = test2,
					  pos2 = pos2,
					  test3_age = test3,
					  pos3 = pos3,
					  rec_age = rec_age,
					  data_type	= data_type				  
                      )

  df_fit$inf_period = df_fit$right_age - df_fit$inf_age + df_fit$left_period

  ### setting up in terms of e/r/s
  df_fit$e_age <- df_fit$left_age
  df_fit$r_age <- df_fit$right_age
  df_fit$s_age <- df_fit$right_age
  df_fit$r_age[which(df_fit$rt_censor == 0)] <- df_fit$r_age[which(df_fit$rt_censor == 0)] - 1
  df_fit$s_age[which(df_fit$rt_censor == 1)] <- NA


  df_fit$e_period <- df_fit$left_period
  df_fit$r_period <- df_fit$right_period
  df_fit$s_period <- df_fit$right_period
  df_fit$r_period[which(df_fit$rt_censor == 0)] <- df_fit$r_period[which(df_fit$rt_censor == 0)] - 1
  df_fit$s_period[which(df_fit$rt_censor == 1)] <- NA

  df_fit$age2date <- df_fit$e_period - df_fit$e_age
  ########################################
  ### setup all data types
  ########################################
  cols_to_keep <- which(colnames(df_fit) %in% c("e_age","r_age","s_age","e_period","r_period","s_period",
                      "rt_censor","cwd_cap","cwd_mort","test2_age","pos2","test3_age","pos3","rec_age","age2date"))
  
  d_fit_sus_cens_posttest <- df_fit[which(data_type == "sus_cens_posttest"),cols_to_keep]
  d_fit_sus_cens_postno <- df_fit[which(data_type == "sus_cens_postno"),cols_to_keep]
  d_fit_sus_mort_posttest <- df_fit[which(data_type == "sus_mort_posttest"),cols_to_keep]
  d_fit_sus_mort_postno <- df_fit[which(data_type == "sus_mort_postno"),cols_to_keep]
  d_fit_icap_cens <- df_fit[which(data_type == "icap_cens"),cols_to_keep]
  d_fit_icap_mort <- df_fit[which(data_type == "icap_mort"),cols_to_keep]
  d_fit_rec_neg_cens_posttest <- df_fit[which(data_type == "rec_neg_cens_posttest"),cols_to_keep]
  d_fit_rec_neg_cens_postno <- df_fit[which(data_type == "rec_neg_cens_postno"),cols_to_keep]
  d_fit_rec_neg_mort <- df_fit[which(data_type == "rec_neg_mort"),cols_to_keep]
  d_fit_rec_pos_cens <- df_fit[which(data_type == "rec_pos_cens"),cols_to_keep]
  d_fit_rec_pos_mort <- df_fit[which(data_type == "rec_pos_mort"),cols_to_keep]
  d_fit_idead <- df_fit[which(data_type == "idead"),cols_to_keep]



  ################################################
  ###
  ### remove the over-abundance infected
  ### at capture, build in 4 options
  ###   no removal, 
  ###   random sample,
  ###   weighted sample by left_age/sum(left_age)
  ###   weighted sample by (1-1/left_age)/sum(1-1/left_age))
  ###
  #################################################


  ########################################
  ### Return values
  ########################################

  return(list(n_ind = n_ind,
              nT_age = nT_age,
              nT_period = nT_period,
              beta0_sus = beta0_sus,
              beta0_inf = beta0_inf,
              period_effect_surv = period_effect_surv,
              age_effect_surv = age_effect_surv,
              foi_age_effect = foi_age_effect,
              hazard_sus = hazard_sus,
              hazard_inf = hazard_inf,
              hazard_foi = hazard_foi,
              survival_prob = survival_prob,
              incidence = incidence,
              right_age = right_age,
              left_age = left_age,
              right_period = right_period,
              left_period = left_period,
              rt_censor = rt_censor,
              inf_age = inf_age,
              inf_status = inf_status,
              test1 = test1,
              test2 = test2,
			  test3 = test3,
              pos1 = pos1,
              pos2 = pos2,
			  pos3 = pos3,
              )
         )
}
