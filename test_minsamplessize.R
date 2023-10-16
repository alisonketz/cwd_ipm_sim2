minsamplesize <- c()

source("03_fun_generate_data.R")


for(i in 0:99){
dat <- ageperiod_surv_foi_sim_data(
            beta0_survival_sus = beta0_survival_sus,
            beta0_survival_inf = beta0_survival_inf,
            foi_age_effect = foi_age_effect,
            age_effect = age_effect_true,
            period_effect = period_effect_true,
            nT_age = nT_age,
            nT_period = nT_period,
            processnum = i
            )

n_fit <- nrow(dat$df_fit)
n_fit_sus_cens_posttest <- nrow(dat$d_fit_sus_cens_posttest)
n_fit_sus_cens_postno <- nrow(dat$d_fit_sus_cens_postno)
n_fit_sus_mort_posttest <- nrow(dat$d_fit_sus_mort_posttest)
n_fit_sus_mort_postno <- nrow(dat$d_fit_sus_mort_postno)
n_fit_icap_cens <- nrow(dat$d_fit_icap_cens)
n_fit_icap_mort <- nrow(dat$d_fit_icap_mort)
n_fit_rec_neg_cens_posttest <- nrow(dat$d_fit_rec_neg_cens_posttest)
n_fit_rec_neg_cens_postno <- nrow(dat$d_fit_rec_neg_cens_postno)
n_fit_rec_neg_mort <- nrow(dat$d_fit_rec_neg_mort)
n_fit_rec_pos_cens <- nrow(dat$d_fit_rec_pos_cens)
n_fit_rec_pos_mort <- nrow(dat$d_fit_rec_pos_mort)
n_fit_idead <- nrow(dat$d_fit_idead)

### test that these data combined are the same
### dimension as the overall generated data
# n_fit_sus_cens_posttest+
# n_fit_sus_cens_postno+
# n_fit_sus_mort_posttest+
# n_fit_sus_mort_postno+
# n_fit_icap_cens+
# n_fit_icap_mort+
# n_fit_rec_neg_cens_posttest+
# n_fit_rec_neg_cens_postno+
# n_fit_rec_neg_mort+
# n_fit_rec_pos_cens+
# n_fit_rec_pos_mort+
# n_fit_idead

datatypes_out <- data.frame(datatype = c("overall","sus_cens_posttest",
"sus_cens_postno",
"sus_mort_posttest",
"sus_mort_postno",
"icap_cens",
"icap_mort",
"rec_neg_cens_posttest",
"rec_neg_cens_postno",
"rec_neg_mort",
"rec_pos_cens",
"rec_pos_mort",
"idead"),
number_samples = c(n_fit,
n_fit_sus_cens_posttest,
n_fit_sus_cens_postno,
n_fit_sus_mort_posttest,
n_fit_sus_mort_postno,
n_fit_icap_cens,
n_fit_icap_mort,
n_fit_rec_neg_cens_posttest,
n_fit_rec_neg_cens_postno,
n_fit_rec_neg_mort,
n_fit_rec_pos_cens,
n_fit_rec_pos_mort,
n_fit_idead))

minsamplesize <- c(minsamplesize,min(datatypes_out$number_samples))
}
