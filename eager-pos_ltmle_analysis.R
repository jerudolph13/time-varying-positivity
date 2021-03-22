
##############################################################################################################  
#
# Project: Compute the per-protocol effect for the EAGeR trial (outcome: pregnancy)
#
# Purpose: Run LTMLE analysis
#
# Author: Jacqueline Rudolph
#
# Last Update: 18 Feb 2020
#
############################################################################################################## 


packages <- c("survival", "ltmle", "tidyverse", "tidyselect", "data.table", "parallel", "splines",
              "SuperLearner", "nnet", "gam", "earth", "glmnet", "e1071", "arm")

for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

`%nin%` = Negate(`%in%`)

#when do we admin censor?
maxt <- 26

#do we want to run LTMLE using GLM? (1=yes, 0=no)
param <- 0

#set seed for SuperLearner
set.seed(123)


# Read in and prepare data ------------------------------------------------

eager <- read.table(file="../data/eager_wide_hcg.txt", sep="\t", header=TRUE)

#Baseline variables
  #Put splines on continuous vars (for parametric implementation)
  age_spline <- bs(eager$age, df=3)
    age1 <- age_spline[ , 1]
    age2 <- age_spline[ , 2]
    age3 <- age_spline[ , 3]
  bmi_spline <- bs(eager$BMI, df=3)
    bmi1 <- bmi_spline[ , 1]
    bmi2 <- bmi_spline[ , 2]
    bmi3 <- bmi_spline[ , 3]
  splines <- data.frame(id=eager$id, age1, age2, age3, bmi1, bmi2, bmi3)
  
  eager_param <- merge(splines, eager, by="id") %>% 
    dplyr::select(-c(age, BMI))
  
  base_param <- paste(names(eager_param)[2:9], collapse="+")
  base_nonparam <- paste(names(eager)[2:5], collapse="+")

#Time-varying confounders
  Lnodes <- as.character()
  bleed <- vars_select(names(eager), starts_with("bleed"))
  nausea <- vars_select(names(eager), starts_with("nausea"))
  for (i in 1:maxt){
    Lnodes[2*i-1] <- nausea[i]
    Lnodes[2*i]<- bleed[i]
  }
  
#Exposure    
  Anodes <- vars_select(names(eager), starts_with(c("treatment", "compliance")))
  # idx <- vars_select(names(eager), starts_with("compliance"))
  # eager[ , idx] <- ifelse(eager[ , idx] > 5/7, 1, 0)
  # eager_param[ , idx] <- ifelse(eager_param[ , idx] > 5/7, 1, 0)
  
#Censoring
  Cnodes <- vars_select(names(eager), starts_with("drop"))
  for (i in 1:maxt) {
    eager[ , Cnodes[i]] <- factor(eager[ , Cnodes[i]], levels=c(0, 1), labels=c("uncensored", "censored"))
    eager_param[ , Cnodes[i]] <- factor(eager_param[ , Cnodes[i]], levels=c(0, 1), labels=c("uncensored", "censored"))
  }

#Outcome
  Ynodes <- vars_select(names(eager), starts_with("delta"))


# Specify models ----------------------------------------------------------
  
  base <- ifelse(param==1, base_param, base_nonparam)
  
  LYformula <- rep(NA, 3*maxt)
    for (i in 1:maxt){
      if (i==1) {
        #Formulae for L(1)
        names(LYformula)[i] <- nausea[i]
          LYformula[i] <- paste("Q.kplus1 ~ treatment", base, sep="+")  
        names(LYformula)[i+1] <- bleed[i]
          LYformula[i+1] <- paste("Q.kplus1 ~ treatment", base, sep="+")
        #Formula for Y(1)
        names(LYformula)[i+2] <- Ynodes[i]
          LYformula[i+2] <- paste("Q.kplus1 ~ treatment", Anodes[i+1], bleed[i], nausea[i], base, sep="+")
      } else if (i==2) {
        #Formula for L(t), t>1
        names(LYformula)[3*i - 2] <- nausea[i]
          LYformula[3*i - 2] <- paste("Q.kplus1 ~ treatment", Anodes[i], nausea[i-1], base, sep="+")
        names(LYformula)[3*i - 1] <- bleed[i]
          LYformula[3*i - 1] <- paste("Q.kplus1 ~ treatment", Anodes[i], bleed[i-1], base, sep="+")
        #Formula for Y(t), t>1
        names(LYformula)[3*i] <- Ynodes[i]
          LYformula[3*i] <- paste("Q.kplus1 ~ treatment", Anodes[i+1], Anodes[i], 
                                  bleed[i], bleed[i-1],
                                  nausea[i], nausea[i-1], base, sep="+")
      } else {
        names(LYformula)[3*i - 2] <- nausea[i]
        LYformula[3*i - 2] <- paste("Q.kplus1 ~ treatment", Anodes[i], nausea[i-1], base, sep="+")
        names(LYformula)[3*i - 1] <- bleed[i]
        LYformula[3*i - 1] <- paste("Q.kplus1 ~ treatment", Anodes[i], bleed[i-1], base, sep="+")
        #Formula for Y(t), t>1
        names(LYformula)[3*i] <- Ynodes[i]
        LYformula[3*i] <- paste("Q.kplus1 ~ treatment", Anodes[i+1], Anodes[i], Anodes[i-1], 
                                bleed[i], bleed[i-1], bleed[i-2],
                                nausea[i], nausea[i-1], nausea[i-2], base, sep="+")
      }
    }

  ACformula <- rep(NA, 2*maxt+1)
    for (i in 1:maxt){
      if (i==1) {
        ACformula[1] <- paste(Anodes[i], "~ eligibility")
        ACformula[2] <- paste(Anodes[i+1], "~", Anodes[1], "+", bleed[i], "+", nausea[i], "+", base, sep="") 
        ACformula[3] <- paste(Cnodes[i], "~", Anodes[1], "+", Anodes[i+1], "+", bleed[i], "+", nausea[i], "+", 
                              base, sep="")
      } else if (i==2) {
        ACformula[2*i] <- paste(Anodes[i+1], "~", Anodes[1], "+", Anodes[i], "+", 
                                bleed[i], "+", bleed[i-1], "+", 
                                nausea[i], "+", nausea[i-1], "+",
                                base, sep="")
        ACformula[2*i+1] <- paste(Cnodes[i], "~",  Anodes[1], "+", Anodes[i+1], "+", Anodes[i], "+",
                                  bleed[i], "+", bleed[i-1], "+",
                                  nausea[i], "+", nausea[i-1], "+",
                                  base, sep="")
      } else {
        ACformula[2*i] <- paste(Anodes[i+1], "~", Anodes[1], "+", Anodes[i], "+", Anodes[i-1], "+",
                                bleed[i], "+", bleed[i-1], "+", bleed[i-2], "+",
                                nausea[i], "+", nausea[i-1], "+", nausea[i-2], "+",
                                base, sep="")
        ACformula[2*i+1] <- paste(Cnodes[i], "~",  Anodes[1], "+", Anodes[i+1], "+", Anodes[i], "+", Anodes[i-1], "+",
                                  bleed[i], "+", bleed[i-1], "+", bleed[i-2], "+",
                                  nausea[i], "+", nausea[i-1], "+", nausea[i-2], "+",
                                  base, sep="") 
      }
    }


# Run LTMLE ---------------------------------------------------------------
    
#Specify our interventions
  #Set to treatment and then continuously compliant
  abar1 <- rep(1, maxt+1)
  
  #Set to placebo and then continuously compliant
  abar0 <- c(0, rep(1, maxt))

#Data frame to hold results
res <- data.frame(estimand=c("R1","R0","RD"), estimate=rep(NA,3), se=rep(NA,3))
  
#Parametric LTMLE
  if (param==1){
    param_ltmle <- ltmle(dplyr::select(eager_param, -id), Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                         #Qform=LYformula,
                         #gform=ACformula,
                         abar=list(treatment=abar1, control=abar0), 
                         survivalOutcome=TRUE,
                         variance.method = "tmle",
                         SL.library=NULL)
    summ <- summary(param_ltmle)  
    res$estimate[1] <- summ$effect.measures$treatment$estimate
     res$se[1] <- summ$effect.measures$treatment$std.dev
    res$estimate[2] <- summ$effect.measures$control$estimate
     res$se[2] <- summ$effect.measures$control$std.dev
    res$estimate[3] <- summ$effect.measures$ATE$estimate
     res$se[3] <- summ$effect.measures$ATE$std.dev
    res$lower <- res$estimate - 1.96*res$se
    res$upper <- res$estimate + 1.96*res$se

    #write.table(res, file="../results/results_param_hcg.txt", sep="\t", row.names=FALSE)
  }

#Non-parametric LTMLE  
  if (param==0) {
    SL.lib <- c("SL.glm", "SL.gam", "SL.earth", "SL.nnet", "SL.bayesglm") #, "SL.glmnet")#, "SL.svm")
      #svm and glmnet would not run due to sparse data

    nonparam_ltmle <- ltmle(dplyr::select(eager, -id), Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                            #Qform=LYformula,
                            #gform=ACformula,
                            abar=list(treatment=abar1, control=abar0), 
                            survivalOutcome=TRUE, 
                            variance.method = "tmle",
                            SL.library=SL.lib)
    summ <- summary(nonparam_ltmle)  
    res$estimate[1] <- summ$effect.measures$treatment$estimate
      res$se[1] <- summ$effect.measures$treatment$std.dev
    res$estimate[2] <- summ$effect.measures$control$estimate
      res$se[2] <- summ$effect.measures$control$std.dev
    res$estimate[3] <- summ$effect.measures$ATE$estimate
      res$se[3] <- summ$effect.measures$ATE$std.dev
    res$lower <- res$estimate - 1.96*res$se
    res$upper <- res$estimate + 1.96*res$se
      
    write.table(res, file="../results/results_sl.txt", sep="\t", row.names=FALSE)
  }    


# Positivity --------------------------------------------------------------

#Let's look at positivity using propensity score overlap
  #Cumulative probabilities of exposure and censoring from first intervention scenario
  g1 <- param_ltmle$cum.g[ , , 1]
  g_val <- as.data.frame(g1[ , seq(2,52, 2)]) #I believe this is pulling out the exposure probabilities
  g_long <- g_val %>% 
    pivot_longer(cols=starts_with("V"),
                 names_to = "week",
                 names_prefix = "V",
                 values_to = "prob_g") %>% 
    select(-week)
  
  eager_small <- eager_param[ , vars_select(names(eager_param), starts_with(c("id", "treatment", "compliance")))]
  eager_long <- pivot_longer(eager_small, 
                             cols=starts_with("compliance"),
                             names_to="week",
                             names_prefix="compliance",
                             values_to = "compliance")
  eager_g <- bind_cols(eager_long, g_long)
  #write_csv(eager_g, file="../results/eager_ltmle_ps.csv")
  

  