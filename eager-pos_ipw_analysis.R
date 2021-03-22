
##############################################################################################################  
#
# Project: Compute the per-protocol effect for the EAGeR trial (outcome: pregnancy)
#
# Purpose: Run IPW analyses
#
# Author: Jacqueline Rudolph
#
# Last Update: 11 Mar 2021
#
############################################################################################################## 


packages <- c("survival", "tidyverse", "tidyselect", "data.table", "parallel", "geepack", "splines")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

`%nin%` = Negate(`%in%`)

# When do we admin censor?
maxt <- 26

# Number of bootstrap resamples
nboot <- 200


# Read in and prepare data ------------------------------------------------

eager <- read.table(file="../data/eager_long_hcg.txt", sep="\t", header=TRUE) %>% 
  mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
  select(-outcome) %>% 
  group_by(id) %>% 
  mutate(compliance_1=lag(compliance, default=0),
         compliance_2=lag(compliance, n=2, default=0),
         nausea_1=lag(nausea, default=0),
         nausea_2=lag(nausea, n=2, default=0),
         bleed_1=lag(bleed, default=0),
         bleed_2=lag(bleed, n=2, default=0)) %>% 
  ungroup()


# Bootstrap ---------------------------------------------------------------

boot.res <- data.frame(
  method=c("Flexible", "Smooth"),
  boot_num=rep(NA, 2),
  r1=rep(NA, 2),
  r0=rep(NA, 2),
  rd=rep(NA, 2),
  stringsAsFactors=FALSE
)

bootrep <- function (r) {
  boot.res$boot_num <- r
  set.seed(r)
  
  firstobs <- eager[eager$week == 1, ]
  samp <- table(firstobs[sample(1:nrow(firstobs), nrow(firstobs), replace=T), (names(eager) == "id")])
  
  # The step below pulls in the original data for boot=0; otherwise grabs all records for the resampled observations
  boot <- NULL
  if(r==0){
    boot <- eager %>% 
      rename(bid = id)
  } else{
    for(zzz in 1:max(samp)){ 
      cc <- eager[eager$id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid <- paste0(cc$id, zzz)
      boot <- rbind(boot, cc)
    }
    boot <- select(boot, -id)
  }

  
# TMLE-IPW ----------------------------------------------------------------
  # "Flexible" IPW
  
  # Compliance
  # Denominator models (to match TMLE, run one per time point)
  boot$ps_trt <- rep(NA, dim(boot)[1])
  for (wk in 1:maxt) {
    dat <- filter(boot, week==wk)
    boot$ps_trt[boot$week==wk] <- glm(compliance ~ treatment + compliance_1 + compliance_2 +
                                        bs(age, df=3) + bs(BMI, df=3) + smoke + eligibility + 
                                        bleed + bleed_1 + bleed_2 + nausea + nausea_1 + nausea_2,
                                        family=binomial(link="logit"), data=dat)$fitted.values
  }
  
  boot <- boot %>% 
    group_by(bid) %>% 
    mutate(cum_ps_trt = cumprod(ps_trt),
           cum_comp = cumprod(compliance),
           unstab_trt_wt1 = (treatment*cum_comp)/(0.5*cum_ps_trt),
           unstab_trt_wt0 = ((1-treatment)*cum_comp)/(0.5*cum_ps_trt))
  
  mean1 <- mean(boot$unstab_trt_wt1)
  mean0 <- mean(boot$unstab_trt_wt0)
  
  boot <- boot %>% 
    mutate(stab_trt_wt1 = unstab_trt_wt1/mean1,
           stab_trt_wt0 = unstab_trt_wt0/mean0)
  
  # Drop out
  boot$ps_drop <- glm((1-drop) ~ treatment + compliance + compliance_1 + compliance_2 + 
                        bs(age, df=3) + bs(BMI, df=3) + smoke + eligibility + 
                        bleed + bleed_1 + bleed_2 + 
                        nausea + nausea_1 + nausea_2 + 
                        as.factor(week) + week:compliance + week:compliance_1 + week:compliance_2 +
                        week:bleed + week:bleed_1 + week:bleed_2 +
                        week:nausea + week:nausea_1 + week:nausea_2,
                        family=binomial(link="logit"), data=boot)$fitted.values
  
  boot <- boot %>% 
    group_by(bid) %>% 
    mutate(cum_ps_drop = cumprod(ps_drop),
           cum_nodrop = cumprod((1-drop)),
           unstab_drop_wt1 = cum_nodrop/(cum_ps_drop))
  
  mean1 <- mean(boot$unstab_drop_wt1)

  boot <- boot %>% 
    mutate(stab_drop_wt1 = unstab_drop_wt1/mean1,
           tot_wt1 = stab_trt_wt1*stab_drop_wt1,
           tot_wt0 = stab_trt_wt0*stab_drop_wt1)
  
  # Weighted survival
  sfit <- summary(survfit(Surv((week-1), week, delta) ~ 1, data=boot, weights=tot_wt1))
  boot.res$r1[1] <- 1 - min(sfit$surv)
  
  sfit <- summary(survfit(Surv((week-1), week, delta) ~ 1, data=boot, weights=tot_wt0))
  boot.res$r0[1] <- 1 - min(sfit$surv)


# Smoothed IPW ------------------------------------------------------------
  # Using a cumulative exposure model
  
  # Run models by randomized treatment
  trt1 <- filter(boot, treatment==1)
  trt0 <- filter(boot, treatment==0)
  
  # Model the probability of complying by treatment arm
  ps.mod <- function(dat) {
    dat2 <- dat
    dat2$ps_den <- glm(compliance ~ compliance_1 + compliance_2 + 
                         bs(age, df=3) + bs(BMI, df=3) + smoke + eligibility + 
                         bleed + bleed_1 + bleed_2 + 
                         nausea + nausea_1 + nausea_2 + 
                         as.factor(week),
                         family=binomial(link="logit"), data=dat2)$fitted.values
    dat2$ps_num <- glm(compliance ~ compliance_1 + compliance_2 + as.factor(week),
                         family=binomial(link="logit"), data=dat2)$fitted.values
    dat2 <- dat2 %>% 
      group_by(bid) %>% 
      mutate(num = compliance*ps_num + (1-compliance)*(1-ps_num),
             den = compliance*ps_den + (1-compliance)*(1-ps_den),
             wt_t = cumprod(num/den)) %>% 
      ungroup()
    
    return(dat2)
  }
  
  wt1 <- ps.mod(trt1)
  wt0 <- ps.mod(trt0)
  
  # Combine and create the cumulative exposure variables
    # X1: took aspirin vs. everything else
    # X2: took nothing vs. took something
  boot2 <- bind_rows(wt1, wt0) %>% 
    group_by(bid) %>% 
    mutate(took_asp = as.numeric(treatment==1 & compliance==1),
           took_nothing = 1 - compliance,
           cum_asp = cumsum(took_asp),
           cum_nothing = cumsum(took_nothing),
           cum_avg_asp = cum_asp/week,
           cum_avg_nothing = cum_nothing/week) %>% 
    ungroup()
  
  # Model the probability of not dropping out
  boot2$drop_den <- glm((1-drop) ~ treatment + compliance + compliance_1 + compliance_2 + 
                        bs(age, df=3) + bs(BMI, df=3) + smoke + eligibility + 
                        bleed + bleed_1 + bleed_2 + 
                        nausea + nausea_1 + nausea_2 + 
                        as.factor(week),
                        family=binomial(link="logit"), data=boot2)$fitted.values
  boot2$drop_num <- glm((1-drop) ~ treatment + compliance + as.factor(week),
                        family=binomial(link="logit"), data=boot2)$fitted.values
  
  boot2 <- boot2 %>% 
    group_by(bid) %>% 
    mutate(cum_drop = cumprod(drop_num/drop_den)) %>% 
    ungroup()
  boot2$wt_c <- ifelse(boot2$drop == 1, 0, boot2$cum_drop)
  boot2 <- mutate(boot2, wt = wt_t*wt_c)
  
  # Model the outcome
  fit <- glm(delta==0 ~ cum_avg_asp + week:cum_avg_asp +
                        cum_avg_nothing + week:cum_avg_nothing +
                        as.factor(week), family=binomial(link="logit"), weights = wt, data=boot2)
  
  # Risk if always took aspirin
  dat1 <- data.frame(week = c(1:26), cum_avg_asp = rep(1, 26), cum_avg_nothing = rep(0, 26))
    dat1$pred <- predict(fit, newdata=dat1, type="response")
    risk1 <- dat1 %>% 
      mutate(s = cumprod(pred),
             r = 1 - s)
  boot.res$r1[2] <- risk1$r[26]

  # Risk if always took placebo
  dat0 <- data.frame(week = c(1:26), cum_avg_asp = rep(0, 26), cum_avg_nothing = rep(0, 26))
  dat0$pred <- predict(fit, newdata=dat0, type="response")
    risk0 <- dat0 %>% 
      mutate(s = cumprod(pred),
             r = 1 - s)
  boot.res$r0[2] <- risk0$r[26]
  
  # Estimate risk difference
  boot.res$rd <- boot.res$r1 - boot.res$r0  
  
  return(boot.res)
}

cores <- detectCores() - 1
all.boot <- mclapply(0:nboot, function(tt) {bootrep(tt)}, mc.cores=cores, mc.set.seed=FALSE)
all.boot <- do.call(rbind, all.boot)

res <- filter(all.boot, boot_num==0)

boot.summ <- all.boot %>% 
  group_by(method) %>% 
  summarize(se = sd(rd))

res <- left_join(res, boot.summ, by="method") %>% 
  mutate(lower = rd - 1.96*se,
         upper = rd + 1.96*se)

write.table(res, file="../results/results_ipw.txt", sep="\t", row.names=FALSE)


