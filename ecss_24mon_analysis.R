## -------------------------------------------------------------------------- ##
## ecss_24mon_analysis.R ---------------------------------------------------- ##
## Author: Peter Norwood, UCLA ---------------------------------------------- ##
## Purpose: clean the datasets to be able to run each analysis -------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## load packages
library(tidyverse)
library(lme4)
library(readxl)
library(zscorer)
library(geepack)
library(lmerTest)
library(blme)
library(knitr)
library(nlme)
library(shiny)
library(optimx)
library(lubridate)
library(psych)

## set working directory
#setwd("~/ECSS_24Mon_Paper")

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 1: DATA LOADING AND MANIPULATION

## This first section loads in the full six month dataset, conducts
## some data manipulation (adding new variables, etc), and also
## loads in "identifier" datasets from the 15 and 24 month surveys.
## This identifier information is PID, time, intervention, etc that we
## later use to join with the previous information.


## six month data
dat_6mon <- read_excel("ecss_complete_dataset_aug2020.xlsx")

## some basic data manipulation
dat_6mon <- dat_6mon %>%
  mutate(month=case_when(time==1~0,
                         time==2~3,
                         time==3~6)) %>%
  group_by(pid) %>%
  ## set clinic from baseline, sometimes clinic later is missing so we 
  ## always use clinic baseline
  mutate(clinic_baseline=clinic[month==0]) %>%
  ungroup() %>%
  ## calculating if more than 4 antenatal visists were taken, one of the outcomes
  mutate(howmanytimeshospitalpast3m=as.numeric(howmanytimeshospitalpast3m),
         howmanytimeshospitalpast3m=ifelse(is.na(howmanytimeshospitalpast3m),
                                           0,howmanytimeshospitalpast3m),
         ## counting up antenatal visits
         antenatalvisitsrecall=as.numeric(antenatalvisitsrecall),
         antenatalvisitsrecall=ifelse(antenatalvisitsrecall==10,0,antenatalvisitsrecall),
         antenatalvisitsrecall=ifelse(antenatalvisitsrecall==98,NA,antenatalvisitsrecall),
         antenatal_gt4=ifelse(antenatalvisitsrecall>3,1,0),
         ## calculating breastfeeding at six months (0 if no breastfeeding at 3 months)
         bf_upto6mons=ifelse(bf_upto3mons==0,0,bf_upto6mons),
         ## calculating exclusive breastfeeding
         exclusive_bf_3mon=ifelse(bf_upto3mons==1 & formula_3mons!=1 & nestum_3mons!=1,1,0),
         exclusive_bf_6mon=ifelse(bf_upto6mons==1 & formula_6mons!=1 & nestum_6mons!=1,1,0)) %>%
  group_by(pid) %>%
  ungroup() %>%
  ## creating an intervention char variable to be clear with tables
  mutate(intv_char=ifelse(intv==1,"intv","control"))


## grab pid, intervention, and clinic for an "essential" dataset that we can
## use to grab this information later
essential <- dat_6mon %>% select(pid,intv,clinic_baseline,time) %>% 
  group_by(pid) %>%
  select(pid,intv,clinic_baseline) %>%
  unique()


## 24 "identifier" datasets with PID information
dat_24mon_identifier <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                                   sheet=4)

dat_24mon_identifier <- dat_24mon_identifier %>% select("[Submission Id]",
                                                        "Participant ID",
                                                        "PID check",
                                                        "Current clinic")

## 24 month phone interviews "identifiers"
dat_24mon_phone_identifier <- read_excel("24_months_phone_survey_completed_20210304-091337.xlsx",
                                         sheet=4)

dat_24mon_phone_identifier <- dat_24mon_phone_identifier %>% select("[Submission Id]",
                                                                    "Participant ID",
                                                                    "PID check")

## 15 month "identifier" datasets with PID information
dat_15mon_identifier <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                                   sheet=4)

dat_15mon_identifier <- dat_15mon_identifier %>% select("[Submission Id]",
                                                        "Participant ID",
                                                        "PID check",
                                                        "Current clinic")

## 15 month phone interviews "identifiers
dat_15mon_phone_identifier <- read_excel("15_months_phone_survey_completed_20210304-122838.xlsx",
                                         sheet=4)

dat_15mon_phone_identifier <- dat_15mon_phone_identifier %>% select("[Submission Id]",
                                                                    "Participant ID",
                                                                    "PID check")


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 2: BREASTFEDING AT 6 MONTHS OUTCOME

## This section tests whether the intervention increases breastfeeding
## rates at six months.

## create a breastfeeding dataset specifically
bf_dat <- dat_6mon %>% 
  filter(month==6) %>% 
  select(pid,intv,intv_char,clinic_baseline,month,bf_upto6mons)

## check missing values for breastfeeding variable
## there are three in the intervention gruop
table(bf_dat$intv_char,bf_dat$bf_upto6mons,useNA="always")

## raw #s for the table
bf_dat %>%
  filter(!is.na(bf_upto6mons)) %>%
  group_by(intv,month) %>%
  summarise(n_bf=sum(bf_upto6mons),
            n=n(),
            perc=round(100*mean(bf_upto6mons),1))

## model for breastfeeding
bf_model <- bglmer(bf_upto6mons~(1|clinic_baseline)+intv,
                   family=binomial(),
                   data=bf_dat)

summary(bf_model)

## asy covariance matrix
bf_cov <- vcov(bf_model)

## fixed effect vector
bf_fixed <- fixed.effects(bf_model)

## test
bf_contrast_6mon <- c(1,1)-c(1,0)
bf_est_6mon <- bf_contrast_6mon %*% bf_fixed
bf_se_6mon <- sqrt(t(bf_contrast_6mon) %*% (bf_cov %*% bf_contrast_6mon))
bf_Z_6mon <- bf_est_6mon/bf_se_6mon
round(1-pnorm(bf_Z_6mon[1,1]),3) ## 0.032

## save smaller dataset for correlation calculation
bf_corr <- bf_dat %>% select(pid,bf_upto6mons) %>% na.omit()

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 3: LOW BIRTH WEIGHT OUTCOME

## This section tests whether the intervention decreases the probability
## of low-birth weight

## create low birth weight dataset
lbw_dat <- dat_6mon %>% 
  mutate(lbw_ind=ifelse(birth_weight<=2.5,1,0)) %>%
  group_by(pid) %>%
  fill(lbw_ind,.direction="updown") %>%
  filter(month==3) %>%
  select(pid,month,clinic_baseline,intv,intv_char,birth_weight,lbw_ind)

## check raw counts of lbw indicators
table(lbw_dat$intv_char,lbw_dat$lbw_ind,useNA="always")

## raw numbers for table
lbw_dat %>% group_by(intv) %>%
  filter(!is.na(lbw_ind)) %>%
  summarise(n_lbw=sum(lbw_ind),
            n=n(),
            perc=round(100*n_lbw/n,1))

## lbw model
lbw_model <- bglmer(lbw_ind~(1|clinic_baseline)+intv,
                    family=binomial,
                    data=lbw_dat)

## lbw model estimates and covariance matrix
lbw_cov <- vcov(lbw_model)
lbw_fixed <- fixed.effects(lbw_model)

## test
lbw_contrast <- c(1,0)-c(1,1)
lbw_est <- lbw_contrast %*% lbw_fixed
lbw_se <- sqrt(t(lbw_contrast) %*% (lbw_cov %*% lbw_contrast))
lbw_Z <- lbw_est/lbw_se
round(1-pnorm(lbw_Z[1,1]),3) ## 0.397

## save smaller dataset for correlation calculation
lbw_corr <- lbw_dat %>% select(pid,lbw_ind) %>% na.omit()

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 4: NO ALCOHOL DURING PREGNANCY

## This section tests whether the intervention decreases the probability
## of drinking during pregnancy

## create alcohol dataset
alc_dat <- dat_6mon %>%
  mutate(alcoholduringpregnancyafterl=ifelse(alcoholduringpregnancyafterl=="NA",
                                             NA,alcoholduringpregnancyafterl),
         noalc_preg=ifelse(alcoholduringpregnancyafterl=="1",1,0)) %>%
  filter(month==0) %>%
  select(pid,intv,intv_char,clinic_baseline,month,noalc_preg)

table(alc_dat$intv_char,alc_dat$noalc_preg,useNA="always")

## raw numbers for table
alc_dat %>%
  group_by(intv) %>%
  summarise(n_noalc=sum(noalc_preg),
            n=n(),
            perc=round(100*n_noalc/n,1))

## no alcohol during pregnancy model
alc_model <- bglmer(noalc_preg~(1|clinic_baseline)+intv,
                    family=binomial,
                    data=alc_dat)

## lbw model estimates and covariance matrix
alc_cov <- vcov(alc_model)
alc_fixed <- fixed.effects(alc_model)

## test
alc_contrast <- c(1,1)-c(1,0)
alc_est <- alc_contrast %*% alc_fixed
alc_se <- sqrt(t(alc_contrast) %*% (alc_cov %*% alc_contrast))
alc_Z <- alc_est/alc_se
round(1-pnorm(alc_Z[1,1]),3) ## 0.247

## save smaller dataset for correlation calculation
alc_corr <- alc_dat %>% select(pid,noalc_preg) %>% na.omit()


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 5: EVER DEPRESSED

## This section tests whether the intervention decreases the probability
## of being depressed across the time period


## read in the 24 month data
epds_24mon <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                             sheet=32)

## calculate epds variable
epds_24mon <- epds_24mon %>%
  mutate(epds_laugh=as.numeric(`Laugh`),
         epds_enjoy=as.numeric(`Enjoyment`),
         epds_self_blame=as.numeric(`Self blame`),
         epds_anxious=as.numeric(`Anxious or worried`),
         epds_panic=as.numeric(`Panicky`),
         epds_piled=as.numeric(`Things piled up`),
         epds_diff_sleep=as.numeric(`Difficulty sleeping`),
         epds_sad=as.numeric(`Sad or miserable`),
         epds_cry=as.numeric(`Crying`),
         epds_harm=as.numeric(`Self harm`),
         epds=epds_laugh+epds_enjoy+epds_self_blame+epds_anxious+
           epds_panic+epds_piled+epds_diff_sleep+epds_sad+
           epds_cry+epds_harm) %>%
  select("[Submission Id]",starts_with("epds"))

## join with the identifier variable to get intv, clinic, etc
epds_24mon <- left_join(dat_24mon_identifier,epds_24mon) %>%
  mutate(pid=`Participant ID`,
         time=5) %>%
  inner_join(essential) %>%
  select(pid,intv,time,clinic_baseline,epds)


## read in 15 month data
epds_15mon <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                             sheet=31)

## calculate epds
epds_15mon <- epds_15mon %>%
  mutate(epds_laugh=as.numeric(`Laugh`),
         epds_enjoy=as.numeric(`Enjoyment`),
         epds_self_blame=as.numeric(`Self blame`),
         epds_anxious=as.numeric(`Anxious or worried`),
         epds_panic=as.numeric(`Panicky`),
         epds_piled=as.numeric(`Things piled up`),
         epds_diff_sleep=as.numeric(`Difficulty sleeping`),
         epds_sad=as.numeric(`Sad or miserable`),
         epds_cry=as.numeric(`Crying`),
         epds_harm=as.numeric(`Self harm`),
         epds=epds_laugh+epds_enjoy+epds_self_blame+epds_anxious+
           epds_panic+epds_piled+epds_diff_sleep+epds_sad+
           epds_cry+epds_harm) %>%
  select("[Submission Id]",starts_with("epds"))

## join with the identifier dataset
epds_15mon <- left_join(dat_15mon_identifier,epds_15mon) %>%
  mutate(pid=`Participant ID`,
         time=4) %>%
  inner_join(essential) %>%
  select(pid,intv,time,clinic_baseline,epds)

## previous data (baseline - 6 month)
epds_prev <- dat_6mon %>% select(pid,intv,time,clinic_baseline,epds)

## full epds dataset
epds_dat <- rbind(epds_prev,epds_15mon,epds_24mon)

## create likely depressed and not depressed variables
epds_dat <- epds_dat %>%
  mutate(epds=as.numeric(epds),
         likely_dep=ifelse(epds>=13,1,0),
         not_dep=1-likely_dep,
         month=case_when(time==1~0,
                         time==2~3,
                         time==3~6,
                         time==4~15,
                         time==5~24))

## repot number of missing values at each time point
epds_dat %>% group_by(intv,month) %>%
  summarise(num_missing=sum(is.na(epds))) %>% kable()

## calculate info dataset to use to later join
epds_info <- epds_dat %>% 
  select(pid,intv,clinic_baseline) %>%
  unique()

## filter to obs where they have all 5 observations
everdep_dat <- epds_dat %>% 
  filter(!is.na(likely_dep)) %>%
  group_by(pid) %>%
  mutate(n=n()) %>%
  filter(n==5) %>%
  summarise(ever_dep=max(likely_dep))

## join with intv/clinic info
everdep_dat <- left_join(everdep_dat,epds_info)

## ever depressed count
everdep_dat %>%
  group_by(intv) %>%
  summarise(n_dep=sum(ever_dep),
            n=n(),
            perc=round(100*n_dep/n,2))

## model
everdep_model <- bglmer(ever_dep~(1|clinic_baseline)+intv,
                         family=binomial(),
                         data=everdep_dat)

summary(everdep_model)

## fixed effects and covariance matrix
everdep_cov <- vcov(everdep_model)
everdep_fixed <- fixed.effects(everdep_model)

## test
everdep_contrast <- c(1,0) - c(1,1)
everdep_est <- everdep_contrast %*% everdep_fixed
everdep_se <- sqrt(t(everdep_contrast) %*% (everdep_cov %*% everdep_contrast))
everdep_Z <- everdep_est/everdep_se
1-pnorm(everdep_Z[1,1]) ## 0.5823702

## save smaller dataset for correlation calculation
everdep_corr <- everdep_dat %>% select(pid,ever_dep) %>% na.omit()

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 6: EVER HAZ or WAZ <-2

## This section tests whether the intervention decreases the probability
## of being stunted or severely underweight

## 24 month growth dataset
growth_24mon <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                           sheet=7)

## calculate weights and heights at 24 months
growth_24mon <- growth_24mon %>%
  mutate(`Baby DoB`=as.numeric(`Baby DoB`),
         dob=as.Date(`Baby DoB`,origin = "1899-12-30"),
         date_rec=as.Date(`Received Date (GMT/UTC)`),
         days_old=date_rec-dob,
         months_old=round(days_old/30),
         weight1=as.numeric(`Infant weight DC 1`),
         weight2=as.numeric(`Infant weight DC 2`),
         height1=as.numeric(`Infant length DC 1`),
         height2=as.numeric(`Infant length DC 2`),
         weight_kg=0.5*(weight1+weight2),
         height_cm=0.5*(height1+height2),
         sex_obs=ifelse(`Sex updated`=="1","2",
                        ifelse(`Sex updated`=="2","1",NA)))


## calculate wazs at 24 months
wazs <- c()
for(i in 1:nrow(growth_24mon)){
  
  waz <- try(zscorer::getWGS(sexObserved = growth_24mon$sex_obs[i],
                             firstPart=growth_24mon$weight_kg[i],
                             secondPart=growth_24mon$months_old[i],
                             index="wfa"))
  
  if(is.numeric(waz)){
    wazs[i]=waz
  }else{
    wazs[i]=NA
  }
  
}
## save to 24 month datasets
growth_24mon$waz <- wazs

## calculate haz at 24 months
hazs <- c()
for(i in 1:nrow(growth_24mon)){
  
  haz <- try(zscorer::getWGS(sexObserved = growth_24mon$sex_obs[i],
                             firstPart=growth_24mon$height_cm[i],
                             secondPart=growth_24mon$months_old[i],
                             index="hfa"))
  
  if(is.numeric(haz)){
    hazs[i]=haz
  }else{
    hazs[i]=NA
  }
  
}
## save to the dataset
growth_24mon$haz <- hazs


## add "identifier" information to the 24 month dataset
growth_24mon <- left_join(dat_24mon_identifier,growth_24mon) %>%
  mutate(pid=`Participant ID`,
         time=5) %>%
  inner_join(essential) %>%
  select(pid,intv,time,clinic_baseline,waz,haz)

## growth from the phone interviews at 24 months
growth_24mon_phone <- read_excel("24_months_phone_survey_completed_20210304-091337.xlsx",
                                 sheet=6)

## calculate heights and weights
growth_24mon_phone <- growth_24mon_phone %>%
  mutate(`Baby DoB`=as.numeric(`Baby DoB`),
         dob=as.Date(`Baby DoB`,origin = "1899-12-30"),
         date_rec=as.Date(`Received Date (GMT/UTC)`),
         days_old=date_rec-dob,
         months_old=round(days_old/30),
         weight1=as.numeric(`Infant weight DC 1`),
         weight2=as.numeric(`Infant weight DC 2`),
         height1=as.numeric(`Infant length DC 1`),
         height2=as.numeric(`Infant length DC 2`),
         weight_kg=0.5*(weight1+weight2),
         height_cm=0.5*(height1+height2),
         sex_obs=ifelse(`Sex updated`=="1","2",
                        ifelse(`Sex updated`=="2","1",NA)))

## calculate waz values from phone survey at 24 months
phone_wazs <- c()
for(i in 1:nrow(growth_24mon_phone)){
  
  waz <- try(zscorer::getWGS(sexObserved = growth_24mon_phone$sex_obs[i],
                             firstPart=growth_24mon_phone$weight_kg[i],
                             secondPart=growth_24mon_phone$months_old[i],
                             index="wfa"))
  
  if(is.numeric(waz)){
    phone_wazs[i]=waz
  }else{
    phone_wazs[i]=NA
  }
  
}
## save
growth_24mon_phone$waz <- phone_wazs

## calculate haz scores from the phone survey at 24 months
phone_hazs <- c()
for(i in 1:nrow(growth_24mon_phone)){
  
  haz <- try(zscorer::getWGS(sexObserved = growth_24mon_phone$`Sex updated`[i],
                             firstPart=growth_24mon_phone$height_cm[i],
                             secondPart=growth_24mon_phone$months_old[i],
                             index="hfa"))
  
  if(is.numeric(haz)){
    phone_hazs[i]=haz
  }else{
    phone_hazs[i]=NA
  }
  
}
## save to dataset
growth_24mon_phone$haz <- phone_hazs

## add identifier information to 24 month phone survey
growth_24mon_phone <- left_join(dat_24mon_phone_identifier,
                                growth_24mon_phone) %>%
  mutate(pid=`Participant ID`,
         time=5) %>%
  inner_join(essential) %>%
  select(pid,intv,time,clinic_baseline,waz,haz)


## fill in phone values for those missing
for(i in 1:nrow(growth_24mon)){
  
  if(is.na(growth_24mon$waz[i])){
    temp_waz <- growth_24mon_phone %>% filter(pid==growth_24mon$pid[i])
    if(nrow(temp_waz)!=0){
      growth_24mon$waz[i] <- temp_waz$waz
    }
  }
  
  if(is.na(growth_24mon$haz[i])){
    temp_haz <- growth_24mon_phone %>% filter(pid==growth_24mon$pid[i])
    if(nrow(temp_haz)!=0){
      growth_24mon$haz[i] <- temp_waz$haz
    }
  }
  
}


## 15 month haz and waz scores -- same steps for new dataset
growth_15mon <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                           sheet=6)

growth_15mon <- growth_15mon %>%
  mutate(`Baby DoB`=as.numeric(`Baby DoB`),
         dob=as.Date(`Baby DoB`,origin = "1899-12-30"),
         date_rec=as.Date(`Received Date (GMT/UTC)`),
         days_old=date_rec-dob,
         months_old=round(days_old/30),
         weight1=as.numeric(`Infant weight DC 1`),
         weight2=as.numeric(`Infant weight DC 2`),
         height1=as.numeric(`Infant length DC 1`),
         height2=as.numeric(`Infant length DC 2`),
         weight_kg=0.5*(weight1+weight2),
         height_cm=0.5*(height1+height2),
         sex_obs=ifelse(`Sex updated`=="1","2",
                        ifelse(`Sex updated`=="2","1",NA)))


## calculate wazs at 15 months
wazs <- c()
for(i in 1:nrow(growth_15mon)){
  
  waz <- try(zscorer::getWGS(sexObserved = growth_15mon$sex_obs[i],
                             firstPart=growth_15mon$weight_kg[i],
                             secondPart=growth_15mon$months_old[i],
                             index="wfa"))
  
  if(is.numeric(waz)){
    wazs[i]=waz
  }else{
    wazs[i]=NA
  }
  
}
## save
growth_15mon$waz <- wazs

## calculate hazs at 15 months
hazs <- c()
for(i in 1:nrow(growth_15mon)){
  
  haz <- try(zscorer::getWGS(sexObserved = growth_15mon$sex_obs[i],
                             firstPart=growth_15mon$height_cm[i],
                             secondPart=growth_15mon$months_old[i],
                             index="hfa"))
  
  if(is.numeric(haz)){
    hazs[i]=haz
  }else{
    hazs[i]=NA
  }
  
}
## save
growth_15mon$haz <- hazs

## save identifier information for main survey at 15 months
growth_15mon <- left_join(dat_15mon_identifier,growth_15mon) %>%
  mutate(pid=`Participant ID`,
         time=4) %>%
  inner_join(essential) %>%
  select(pid,intv,time,clinic_baseline,waz,haz)

## growth from the phone interviews at 15 months
growth_15mon_phone <- read_excel("15_months_phone_survey_completed_20210304-122838.xlsx",
                                 sheet=6)

## we had issues loading in these date variables, because it was only a few,
## i just manually inputted them in
growth_15mon_phone$bb_dob <- c("1/6/2019",
                               "1/6/2019",
                               "2/26/2019",
                               "2/22/2019",
                               "1/4/2019",
                               "2/21/2019",
                               "1/9/2019",
                               "2/1/2019",
                               "2/16/2019",
                               "2/25/2019",
                               "12/27/2018",
                               "3/7/2019",
                               "1/29/2019",
                               "2/19/2019",
                               "1/28/2019",
                               "1/1/2019",
                               "2/22/2019",
                               "1/5/2019",
                               "2/1/2019",
                               "1/18/2019",
                               "4/3/2019",
                               "1/15/2019",
                               "4/2/2019")


## calculate heights and weights
growth_15mon_phone <- growth_15mon_phone %>%
  mutate(`Baby DoB`=as.numeric(`Baby DoB`),
         #dob_pre=substr(`Baby DoB`,1,7),
         dob=mdy(bb_dob),
         date_rec=as.Date(`Received Date (GMT/UTC)`),
         days_old=date_rec-dob,
         months_old=round(days_old/30),
         weight1=as.numeric(`Infant weight DC 1`),
         weight2=as.numeric(`Infant weight DC 2`),
         height1=as.numeric(`Infant length DC 1`),
         height2=as.numeric(`Infant length DC 2`),
         weight_kg=0.5*(weight1+weight2),
         height_cm=0.5*(height1+height2),
         sex_obs=ifelse(`Sex updated`=="1","2",
                        ifelse(`Sex updated`=="2","1",NA)))

## calculate waz scores from 15 month phone survey
phone_wazs <- c()
for(i in 1:nrow(growth_15mon_phone)){
  
  waz <- try(zscorer::getWGS(sexObserved = growth_15mon_phone$sex_obs[i],
                             firstPart=growth_15mon_phone$weight_kg[i],
                             secondPart=growth_15mon_phone$months_old[i],
                             index="wfa"))
  
  if(is.numeric(waz)){
    phone_wazs[i]=waz
  }else{
    phone_wazs[i]=NA
  }
  
}
## save
growth_15mon_phone$waz <- phone_wazs

## calcualte haz scores from 15 month phone survey
phone_hazs <- c()
for(i in 1:nrow(growth_15mon_phone)){
  
  haz <- try(zscorer::getWGS(sexObserved = growth_15mon_phone$`Sex updated`[i],
                             firstPart=growth_15mon_phone$height_cm[i],
                             secondPart=growth_15mon_phone$months_old[i],
                             index="hfa"))
  
  if(is.numeric(haz)){
    phone_hazs[i]=haz
  }else{
    phone_hazs[i]=NA
  }
  
}
## save
growth_15mon_phone$haz <- phone_hazs

## add identifier info
growth_15mon_phone <- left_join(dat_15mon_phone_identifier,growth_15mon_phone) %>%
  mutate(pid=`Participant ID`,
         time=4) %>%
  inner_join(essential) %>%
  select(pid,intv,time,clinic_baseline,waz,haz)


## fill in phone values for those missing at 15 months
for(i in 1:nrow(growth_15mon)){
  
  if(is.na(growth_15mon$waz[i])){
    temp_waz <- growth_15mon_phone %>% filter(pid==growth_15mon$pid[i])
    if(nrow(temp_waz)!=0){
      growth_15mon$waz[i] <- temp_waz$waz
    }
  }
  
  if(is.na(growth_15mon$haz[i])){
    temp_haz <- growth_15mon_phone %>% filter(pid==growth_15mon$pid[i])
    if(nrow(temp_haz)!=0){
      growth_15mon$haz[i] <- temp_waz$haz
    }
  }
  
}


## 3 and 6 month dataset
growth_pre <- dat_6mon %>% select(pid,intv,time,clinic_baseline,waz,haz)

## combine datasets to a "full" dataset with 3,6,15,24 month growth data
## create less than -2 sd variables
growth_full <- rbind(growth_pre,growth_15mon,growth_24mon) %>%
  mutate(waz_lt2=ifelse(waz<=-2,1,0),
         waz_gt2=1-waz_lt2,
         haz_lt2=ifelse(haz<=-2,1,0),
         haz_gt2=1-haz_lt2) %>%
  mutate(month=case_when(time==1~0,
                         time==2~3,
                         time==3~6,
                         time==4~15,
                         time==5~24))

## growth essential info
growth_info <- growth_full %>%
  select(pid,intv,clinic_baseline) %>%
  unique()

## check missing values for waz
growth_full %>% 
  filter(month>0) %>%
  group_by(intv,month) %>%
  summarise(waz_missing=sum(is.na(waz_lt2)),
            haz_missing=sum(is.na(haz_lt2))) %>%
  kable()

## -------------------------------------------------------------------------- ##
## every have waz <-2sd calculations
complete_waz <- growth_full %>% 
  filter(month>0 & !is.na(waz_lt2)) %>%
  group_by(pid) %>%
  mutate(n=n()) %>%
  filter(n==4) %>%
  summarise(ever_waz_lt2=max(waz_lt2))

## join with intv/clinic info
waz_dat <- left_join(complete_waz,growth_info)

## every waz <-2 count
waz_dat %>%
  group_by(intv) %>%
  summarise(n_dep=sum(ever_waz_lt2),
            n=n(),
            perc=round(100*n_dep/n,2))

waz_model <- bglmer(ever_waz_lt2~(1|clinic_baseline)+intv,
                    family=binomial(),
                    data=waz_dat)

summary(waz_model)

## fixed effects and covariance matrix
waz_cov <- vcov(waz_model)
waz_fixed <- fixed.effects(waz_model)

## test
waz_contrast <- c(1,0) - c(1,1)
waz_est <- waz_contrast %*% waz_fixed
waz_se <- sqrt(t(waz_contrast) %*% (waz_cov %*% waz_contrast))
waz_Z <- waz_est/waz_se
1-pnorm(waz_Z[1,1]) ## 0.07805849

## -------------------------------------------------------------------------- ##


## HAZ final dataset and modeling
complete_haz <- growth_full %>% 
  filter(month>0 & !is.na(haz_lt2)) %>%
  group_by(pid) %>%
  mutate(n=n()) %>%
  filter(n==4) %>%
  summarise(ever_haz_lt2=max(haz_lt2))

## join with intv/clinic info
haz_dat <- left_join(complete_haz,growth_info)

## ever haz <-2 count
haz_dat %>%
  group_by(intv) %>%
  summarise(n_dep=sum(ever_haz_lt2),
            n=n(),
            perc=round(100*n_dep/n,2))

## haz model
haz_model <- bglmer(ever_haz_lt2~(1|clinic_baseline)+intv,
                    family=binomial(),
                    data=haz_dat)

summary(haz_model)

## fixed effects and covariance matrix
haz_cov <- vcov(haz_model)
haz_fixed <- fixed.effects(haz_model)

## test
haz_contrast <- c(1,0) - c(1,1)
haz_est <- haz_contrast %*% haz_fixed
haz_se <- sqrt(t(haz_contrast) %*% (haz_cov %*% haz_contrast))
haz_Z <- haz_est/haz_se
1-pnorm(haz_Z[1,1]) ## 0.3932384

## -------------------------------------------------------------------------- ##

## remove the previous datasets that keep giving not-imporant warnings
growth_15mon <- NULL
growth_24mon <- NULL
growth_24mon_phone <- NULL
growth_15mon_phone <- NULL
temp_haz <- NULL
temp_waz <- NULL

## save waz and haz datasets for correlation
waz_corr <- waz_dat %>% select(pid,ever_waz_lt2) %>% na.omit()
haz_corr <- haz_dat %>% select(pid,ever_haz_lt2) %>% na.omit()

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 7: Child support grant by month 6

## This section tests whether the intervention increases the probability
## that a mother get a CSG for her child

## create child support grant variable and dataset
csg_dat <- dat_6mon %>%
  filter(time==3) %>%
  mutate(csg=ifelse(csg_0=="1",1,0))

table(csg_dat$csg,useNA="always")

csg_dat %>%
  group_by(intv) %>%
  summarise(n_csg=sum(csg),
            n=n(),
            perc=100*mean(csg))


## child support grant model
csg_model <- bglmer(csg~(1|clinic_baseline)+intv,
                    family=binomial(),
                    data=csg_dat)

summary(csg_model)

## save model info
csg_cov <- vcov(csg_model)
csg_fixed <- fixed.effects(csg_model)

## run test
csg_contrast <- c(1,1)-c(1,0)
csg_est <- csg_contrast %*% csg_fixed
csg_se <- sqrt(t(csg_contrast) %*% (csg_cov %*% csg_contrast))
csg_Z <- csg_est/csg_se
pnorm(csg_Z[1,1]) ## 0.5312229

## save smaller dataset for correlation
csg_corr <- csg_dat %>% select(pid,csg) %>% na.omit()

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 8: Immunizations up to date

## This section tests whether the intervention increases the probability
## that a child has its immunizations up to date at 6, 15, and 24 months

## load 15 month survey
immune_15mon <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                           sheet=16)

## join with identifier information
immune_15mon <- inner_join(immune_15mon,dat_15mon_identifier) %>%
  select(c("PID check","Current clinic","Immunizations up to date")) %>%
  mutate(time=4,
         intv=NA,
         pid=`PID check`)


## 15 month phone surveys
immune_15mon_phone <- read_excel("15_months_phone_survey_completed_20210304-122838.xlsx",
                                 sheet=9)

## join with identifier information
immune_15mon_phone <- inner_join(immune_15mon_phone,dat_15mon_phone_identifier) %>%
  select(c("PID check","Immunizations up to date")) %>%
  mutate(time=4,
         intv=NA,
         pid=`PID check`)

## fill in phone values for those missing
for(i in 1:nrow(immune_15mon)){
  
  if(is.na(immune_15mon$`Immunizations up to date`[i])){
    temp <- immune_15mon_phone %>% filter(pid==immune_15mon$pid[i])
    if(nrow(temp)!=0){
      immune_15mon$`Immunizations up to date`[i] <- temp$`Immunizations up to date`
    }
  }
  
}

## 24 month survey
immune_24mon <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                           sheet=18)

## join with identifier information
immune_24mon <- inner_join(immune_24mon,dat_24mon_identifier) %>%
  select(c("PID check","Current clinic","Immunizations up to date")) %>%
  mutate(time=5,
         intv=NA,
         pid=`PID check`)

## phone surveys
immune_24mon_phone <- read_excel("24_months_phone_survey_completed_20210304-091337.xlsx",
                                 sheet=11)

## join with identifier information
immune_24mon_phone <- inner_join(immune_24mon_phone,dat_24mon_phone_identifier) %>%
  select(c("PID check","Immunizations up to date")) %>%
  mutate(time=5,
         intv=NA,
         pid=`PID check`)

## fill in phone values for those missing
for(i in 1:nrow(immune_24mon)){
  
  if(is.na(immune_24mon$`Immunizations up to date`[i])){
    temp <- immune_24mon_phone %>% filter(pid==immune_24mon$pid[i])
    if(nrow(temp)!=0){
      immune_24mon$`Immunizations up to date`[i] <- temp$`Immunizations up to date`
    }
  }
  
}

## combine all of them 
immune_pre <- dat_6mon %>% 
  select(c(pid,intv,time,clinic_baseline,immunizationsuptodate_0))

colnames(immune_pre) <- c("pid","intv","time",
                          "clinic_baseline","immune_upto_date")

immune_15mon <- immune_15mon %>% mutate(clinic_baseline=NA) %>%
  select(pid,intv,time,clinic_baseline,"Immunizations up to date")

colnames(immune_15mon) <- c("pid","intv","time",
                            "clinic_baseline","immune_upto_date")

immune_24mon <- immune_24mon %>% mutate(clinic_baseline=NA) %>%
  select(pid,intv,time,clinic_baseline,"Immunizations up to date")

colnames(immune_24mon) <- c("pid","intv","time",
                            "clinic_baseline","immune_upto_date")


## full dataset
immune_full <- rbind(immune_pre,immune_15mon,immune_24mon) %>%
  group_by(pid) %>%
  fill(clinic_baseline,.direction="updown") %>%
  fill(intv,.direction="updown") %>%
  filter(time>2) %>%
  mutate(immune_upto_date=as.numeric(immune_upto_date)) %>%
  mutate(month=case_when(time==1~0,
                         time==2~3,
                         time==3~6,
                         time==4~15,
                         time==5~24))

## check missing values
immune_full %>%
  group_by(intv,month) %>%
  summarise(n_na=sum(is.na(immune_upto_date)))

## immune ever calculations and model
immune_info <- immune_full %>%
  select(pid,intv,clinic_baseline) %>%
  unique()

## immunization at all time points dataset
immune_dat <- immune_full %>% 
  filter(month>3 & !is.na(immune_upto_date)) %>%
  group_by(pid) %>%
  mutate(n=n()) %>%
  filter(n==3) %>%
  summarise(ever_not_immune=min(immune_upto_date)) %>%
  mutate(always_upto=1-ever_not_immune)

## join with intv/clinic info
immune_dat <- left_join(immune_dat,immune_info)

## always immune count
immune_dat %>%
  group_by(intv) %>%
  summarise(n_upto=sum(always_upto),
            n=n(),
            perc=round(100*n_upto/n,2))

## immunization model
immune_model <- bglmer(always_upto~(1|clinic_baseline)+intv,
                       family=binomial(),
                       data=immune_dat)

summary(immune_model)

## fixed effects and covariance matrix
immune_cov <- vcov(immune_model)
immune_fixed <- fixed.effects(immune_model)

## test
immune_contrast <- c(1,1) - c(1,0)
immune_est <- immune_contrast %*% immune_fixed
immune_se <- sqrt(t(immune_contrast) %*% (immune_cov %*% immune_contrast))
immune_Z <- immune_est/immune_se
1-pnorm(immune_Z[1,1]) ## 0.6766569

## save smaller dataset for correlation
immune_corr <- immune_dat %>% select(pid,always_upto)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 9: No hospitalizations at 6, 15, and 24 Months

## This section tests whether the intervention decreases the probability
## that a child goes to the hospital in the first 24 months

## 15 month data
hos_15mon <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                        sheet=23)

## add identifier information
hos_15mon <- hos_15mon %>%
  mutate(hos=`How many times hospital past 3 months`,
         hos=as.numeric(hos),
         hos_no=ifelse(is.na(hos)|hos<1,1,0)) %>%
  inner_join(dat_15mon_identifier) %>%
  mutate(pid=`PID check`,
         intv=NA,
         clinic_baseline=NA,
         time=4) %>%
  select(pid,intv,time,clinic_baseline,hos,hos_no)

## 24 month data
hos_24mon <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                        sheet=25)

## add identifier information
hos_24mon <- hos_24mon %>%
  mutate(hos=`How many times hospital past 3 months`,
         hos=as.numeric(hos),
         hos_no=ifelse(is.na(hos)|hos<1,1,0)) %>%
  inner_join(dat_24mon_identifier) %>%
  mutate(pid=`PID check`,
         intv=NA,
         clinic_baseline=NA,
         time=5) %>%
  select(pid,intv,time,clinic_baseline,hos,hos_no)

hos_pre <- dat_6mon %>%
  mutate(hos=howmanytimesherbalistpast3,
         hos=as.numeric(hos),
         hos_no=ifelse(is.na(hos)|hos<1,1,0)) %>%
  select(pid,intv,time,clinic_baseline,hos,hos_no)


hos_full <- rbind(hos_pre,
                  hos_15mon,hos_24mon) %>%
  group_by(pid) %>%
  fill(clinic_baseline,.direction="updown") %>%
  fill(intv,.direction="updown") %>%
  ungroup() %>%
  mutate(month=case_when(time==1~0,
                         time==2~3,
                         time==3~6,
                         time==4~15,
                         time==5~24)) %>%
  filter(month>3)

hos_full %>% 
  group_by(month,intv) %>% 
  summarise(sum_na=sum(is.na(hos_no)))


## ever hospitalized data
hos_dat <- hos_full %>%
  group_by(pid,clinic_baseline,intv) %>%
  mutate(visit=1-hos_no) %>%
  summarise(hos_ever=max(visit))

## get raw #s for table
hos_dat %>% group_by(intv) %>%
  summarise(hos_ever=sum(hos_ever),
            n=n(),
            prop=round(100*hos_ever/n,1))

## modeling
hos_model <- bglmer(hos_ever~(1|clinic_baseline)+intv,
                    family=binomial(),
                    data=hos_dat)
summary(hos_model)

## test 
hos_cov <- vcov(hos_model)
hos_fixed <- fixed.effects(hos_model)

hos_contrast <- c(1,0)-c(1,1)
hos_est <- hos_contrast %*% hos_fixed
hos_se <- sqrt(t(hos_contrast) %*% (hos_cov %*% hos_contrast))
hos_Z <- hos_est/hos_se
1-pnorm(hos_Z[1,1]) ## 0.2988315

## save smaller data
hos_corr <- hos_dat %>% ungroup() %>% select(pid,hos_ever)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 10: Adherence to Antenatal visits (>=4) for all mothers

## This section tests whether the intervention increases the probability
## that a mother goes to 4 or more antenatal doctors appointments

## dataset
antenatal_dat <- dat_6mon %>%
  mutate(antenatalvisitsrecall=ifelse(antenatalvisitsrecall=="NA",NA,antenatalvisitsrecall),
         antenatalvisitsrecall=as.numeric(antenatalvisitsrecall),
         antenatalvisitsrecall=ifelse(antenatalvisitsrecall==10,0,antenatalvisitsrecall),
         antenatalvisitsrecall=ifelse(antenatalvisitsrecall==98,0,antenatalvisitsrecall),
         antenatal_gt4=ifelse(antenatalvisitsrecall>3,1,0)) %>%
  filter(time==2)

## check missing values
antenatal_dat %>% 
  group_by(intv) %>%
  summarise(sum_na=sum(is.na(antenatal_gt4)))

## raw counts
antenatal_dat %>% 
  filter(!is.na(antenatal_gt4)) %>%
  group_by(intv) %>%
  summarise(gt4=sum(antenatal_gt4),
            n=n(),
            perc=round(100*gt4/n,1)) 


## modeling antenatal visits
antenatal_model <- bglmer(antenatal_gt4~(1|clinic)+intv,
                          family=binomial(),
                          data=antenatal_dat)

summary(antenatal_model) 

## save model info
antenatal_cov <- vcov(antenatal_model)
antenatal_fixed <- fixed.effects(antenatal_model)

## test
antenatal_contrast <- c(1,1)-c(1,0)
antenatal_est <- antenatal_contrast %*% antenatal_fixed
antenatal_se <- sqrt(t(antenatal_contrast) %*% (antenatal_cov %*% antenatal_contrast))
antenatal_Z <- antenatal_est/antenatal_se
1-pnorm(antenatal_Z[1,1]) ## 0.5388112

## save smaller dataset for correlation
antenatal_corr <- antenatal_dat %>% select(pid,antenatal_gt4)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 11: Adherence to child HIV preventive measures (HIV+ mothers)

## This section tests whether the intervention increases the probability
## that a mother (a) exclusively breastfeeds for six months, (b) gives bactrim
## at birth, (c) gives nvp at birth, (d) tests a child at birth, and 
## (e) goes to the clinic to get the test

## grab baseline/3/6month data for identifier
adher_info <- dat_6mon %>% 
  filter(time==1) %>%
  select(pid,intv,clinic_baseline,hivpos)

## exclusive breastfeeding
excbf_dat <- dat_6mon %>%
  filter(hivpos==1) %>%
  group_by(pid) %>%
  fill(exclusive_bf_3mon,.direction="updown") %>%
  fill(exclusive_bf_6mon,.direction="updown") %>%
  ungroup() %>%
  filter(time==1) %>%
  select(pid,exclusive_bf_6mon)

## bactrim
bactrim_dat <- dat_6mon %>%
  filter(hivpos==1 & time==2) %>%
  select(pid,bactrimsyrup)

## nvp
nvp_dat <- dat_6mon %>%
  filter(hivpos==1 & time==2) %>%
  select(pid,nvpfrombirthbaby)

## test for hiv
test_dat <- dat_6mon %>%
  filter(hivpos==1 & time==2) %>%
  select(pid,childtestedatbirth)

## go back to get test
pcr_dat <- dat_6mon %>%
  filter(hivpos==1 & time==2) %>%
  select(pid,pcrresultfrombirth)

## join these together
child_adher_dat <- left_join(adher_info,excbf_dat)
child_adher_dat <- left_join(child_adher_dat,bactrim_dat)
child_adher_dat <- left_join(child_adher_dat,nvp_dat)
child_adher_dat <- left_join(child_adher_dat,test_dat)
child_adher_dat <- left_join(child_adher_dat,pcr_dat)

## create summed variable
child_adher_dat <- child_adher_dat %>%
  filter(hivpos==1) %>%
  mutate(bactrim_num=as.numeric(bactrimsyrup),
         bactrim_binary = ifelse(bactrim_num==1 | bactrim_num==2,1,0),
         nvp_num=as.numeric(nvpfrombirthbaby),
         nvp_binary = ifelse(nvp_num==1 | nvp_num==2,1,0),
         test_num=as.numeric(childtestedatbirth),
         test_binary=ifelse(test_num==1,1,0),
         pcr_num=as.numeric(pcrresultfrombirth),
         pcr_binary=ifelse(pcr_num==1 | pcr_num==2 | pcr_num==3,1,0)) %>%
  mutate(child_sum=exclusive_bf_6mon+bactrim_binary+nvp_binary+test_binary+pcr_binary,
         child_all=ifelse(child_sum==5,1,0))

## run summaries for table
child_adher_dat %>%
  filter(!is.na(child_sum)) %>%
  group_by(intv) %>%
  summarise(mean=mean(child_sum),
            sd=sd(child_sum)) %>%
  kable(digits=1)

## child adherence model
child_adher_model <- blmer(child_sum~(1|clinic_baseline)+intv,
                     data=child_adher_dat)

summary(child_adher_model)
qqnorm(residuals(child_adher_model)); qqline(residuals(child_adher_model))

## save model info
child_adher_cov <- vcov(child_adher_model)
child_adher_fixed <- fixed.effects(child_adher_model)

## test for difference in #
child_adher_contrast <- c(1,1)-c(1,0)
child_adher_est <- child_adher_contrast %*% child_adher_fixed
child_adher_se <- sqrt(t(child_adher_contrast) %*% (child_adher_cov %*% child_adher_contrast))
child_adher_Z <- child_adher_est/child_adher_se
round(1-pnorm(child_adher_Z[1,1]),3) ## 0.241

## save for correlation
child_adher_corr <- child_adher_dat %>% select(pid,child_sum)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 12: Adherence to ARV medication (HIV+ mothers)

## This section tests whether the intervention increases the probability
## that a mother is adherent to her arv medication at all timepoints

## we need to pull directly from the original surveys rather than combined
## dataset due to data cleaning issues

## baseline arv answers
arv_baseline <- read.csv("ECSS_baseline.csv")

arv_baseline <- arv_baseline %>%
  mutate(pid=`PID.check`,
         intv=NA,
         clinic_baseline=NA,
         arvadher=`Adherence.last.30.days`,
         time=1,
         month=0) %>%
  select(pid,clinic_baseline,arvadher,time,month)

table(arv_baseline$arvadher,useNA="always")

## 3 month arv answers
arv_3mon <- read.csv("ECSS_3mo.csv")

arv_3mon <- arv_3mon %>%
  mutate(pid=`PID.check`,
         intv=NA,
         clinic_baseline=NA,
         arvadher=`Adherence.last.30.days`,
         time=2,
         month=3) %>%
  select(pid,clinic_baseline,arvadher,time,month)

table(arv_3mon$arvadher,useNA="always")

## 6 month arv answers
arv_6mon <- read.csv("ECSS_6mo.csv")

arv_6mon <- arv_6mon %>%
  mutate(pid=`PID.check`,
         intv=NA,
         clinic_baseline=NA,
         arvadher=`Adherence.last.30.days`,
         time=3,
         month=6) %>%
  select(pid,clinic_baseline,arvadher,time,month)

table(arv_6mon$arvadher,useNA="always")

## 15 month arv answers
arv_15mon <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                             sheet=28)

arv_15mon <- arv_15mon %>%
  inner_join(dat_15mon_identifier) %>%
  mutate(pid=`PID check`,
         intv=NA,
         clinic_baseline=NA,
         time=4,
         arvadher=`Adherence last 30 days`,
         month=15) %>%
  select(pid,clinic_baseline,arvadher,time,month)

table(arv_15mon$arvadher,useNA="always")

## 24 month arv answers
arv_24mon <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                             sheet=29)

arv_24mon <- arv_24mon %>%
  inner_join(dat_24mon_identifier) %>%
  mutate(pid=`PID check`,
         intv=NA,
         clinic_baseline=NA,
         time=5,
         arvadher=`Adherence last 30 days`,
         month=24) %>%
  select(pid,clinic_baseline,arvadher,time,month)

table(arv_24mon$arvadher,useNA="always")

arv_dat <- rbind(arv_baseline,
                  arv_3mon,
                  arv_6mon,
                  arv_15mon,
                  arv_24mon)

## classify a binary arv adherence variable to be 1 if 
## the response is "very good" or "excellent" and 0 otherwise
arv_dat <- arv_dat %>%
  mutate(arvadher_binary=case_when(arvadher %in% c("5","6")~1,
                                   arvadher=="-"~NA_real_,
                                   is.na(arvadher)~NA_real_,
                                   arvadher %in% c("1","2","3","4","95")~0)) %>%
  mutate(arvadher_baseline=ifelse(time==1,arvadher,NA),
         arvadher_3mon=ifelse(time==2,arvadher,NA),
         arvadher_6mon=ifelse(time==3,arvadher,NA),
         arvadher_15mon=ifelse(time==4,arvadher,NA),
         arvadher_24mon=ifelse(time==5,arvadher,NA),
         arvadher_binary_baseline=ifelse(time==1,arvadher_binary,NA),
         arvadher_binary_3mon=ifelse(time==2,arvadher_binary,NA),
         arvadher_binary_6mon=ifelse(time==3,arvadher_binary,NA),
         arvadher_binary_15mon=ifelse(time==4,arvadher_binary,NA),
         arvadher_binary_24mon=ifelse(time==5,arvadher_binary,NA)) %>%
  group_by(pid) %>%
  fill(arvadher_baseline,.direction="updown") %>%
  fill(arvadher_3mon,.direction="updown") %>%
  fill(arvadher_6mon,.direction="updown") %>%
  fill(arvadher_15mon,.direction="updown") %>%
  fill(arvadher_24mon,.direction="updown") %>%
  fill(arvadher_binary_baseline,.direction="updown") %>%
  fill(arvadher_binary_3mon,.direction="updown") %>%
  fill(arvadher_binary_6mon,.direction="updown") %>%
  fill(arvadher_binary_15mon,.direction="updown") %>%
  fill(arvadher_binary_24mon,.direction="updown") %>%
  ungroup() %>%
  filter(time==1) %>%
  select(pid,
         arvadher_baseline,arvadher_3mon,arvadher_6mon,
         arvadher_15mon,arvadher_24mon,arvadher_binary_baseline,arvadher_binary_3mon,arvadher_binary_6mon,
         arvadher_binary_15mon,arvadher_binary_24mon)

## join with adherence baseline info
arv_dat <- left_join(arv_dat,adher_info)

## filter to hiv positive mothers and add the sum
arv_dat <- arv_dat %>%
  filter(hivpos==1) %>%
  mutate(arv_sum=arvadher_binary_baseline+arvadher_binary_3mon+arvadher_binary_6mon+
           arvadher_binary_15mon+arvadher_binary_24mon)

## run summaries
arv_dat %>% 
  filter(!is.na(arv_sum)) %>%
  group_by(intv) %>%
  summarise(round(mean(arv_sum),1),
            round(sd(arv_sum),1))

## arv adherence model
arv_model <- blmer(arv_sum~(1|clinic_baseline)+intv,
                   data=arv_dat)
summary(arv_model)
qqnorm(residuals(arv_model)); qqline(residuals(arv_model))

## asy covariance matrix
arv_cov <- vcov(arv_model)

## fixed effect vector
arv_fixed <- fixed.effects(arv_model)

## test for difference in #
arv_contrast <- c(1,1)-c(1,0)
arv_est <- arv_contrast %*% arv_fixed
arv_se <- sqrt(t(arv_contrast) %*% (arv_cov %*% arv_contrast))
arv_Z <- arv_est/arv_se
round(1-pnorm(arv_Z[1,1]),3) ## 0.008

## save smaller dataset for correlation
arv_corr <- arv_dat %>% select(pid,arv_sum)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## BLOCK 13: Developmental milestones

## This section tests whether the intervention increases the mean number
## of developmental milestones hit at 6, 15, and 24 months

## save identifier information
dev_info <- dat_6mon %>% 
  select(pid,intv,clinic_baseline)

## six month data
dev_6mon <- dat_6mon %>%
  filter(time==3) %>%
  filter(who1swscriteria!="NA") %>%
  mutate(dev1=ifelse(who1swscriteria=="3",1,0),
         dev2=ifelse(who2hkccriteria=="3",1,0),
         dev3=ifelse(who3swacriteria=="3",1,0),
         dev4=ifelse(who4wwacriteria=="3",1,0),
         dev5=ifelse(who5sacriteria=="3",1,0),
         dev6=ifelse(finemotorpaper_griffiths=="3",1,0),
         dev7=ifelse(personalsocial_griffiths_mirror=="3",1,0),
         dev_sum=dev1+dev2+dev3+dev4+dev5+dev6+dev7,
         dev_sum5=dev1+dev2+dev3+dev4+dev5)

## check missing values
table(dev_6mon$dev_sum,useNA="always")

## 15 month developmental milestones
dev_15mon <- read_excel("EC_RCT_15_months__20210304-105302.xlsx",
                        sheet=7)
dev_15mon <- dev_15mon %>% inner_join(dat_15mon_identifier)

dev_15mon_phone <- read_excel("15_months_phone_survey_completed_20210304-122838.xlsx",
                              sheet=7)
dev_15mon_phone <- inner_join(dev_15mon_phone,dat_15mon_phone_identifier)
dev_15mon_phone$`Current clinic` <- NA

## join together
dev_15mon <- rbind(dev_15mon,dev_15mon_phone)
dev_15mon <- dev_15mon %>% mutate(pid=`PID check`)

## add other info
dev_15mon <- left_join(dev_15mon,dev_info)

## remove missing values
dev_15mon <- dev_15mon %>%
  filter(`WHO 1 SWS criteria`!="-")

## remove duplicates
dev_15mon <- unique(dev_15mon)

## calculate the sum
dev_15mon <- dev_15mon %>%
  mutate(dev1=ifelse(`WHO 1 SWS criteria`=="3",1,0),
         dev2=ifelse(`WHO 2 HKC criteria`=="3",1,0),
         dev3=ifelse(`WHO3 SWA criteria`=="3",1,0),
         dev4=ifelse(`WHO 4 WWA criteria`=="3",1,0),
         dev5=ifelse(`WHO5 SA criteria`=="3",1,0),
         dev6=ifelse(`WHO6 WA criteria`=="3",1,0),
         dev7=ifelse(Speech=="3",1,0),
         dev8=ifelse(`RTHC Social_Hearing`=="1",1,0),
         dev_sum=dev1+dev2+dev3+dev4+dev5+dev6+dev7+dev8,
         all8=ifelse(dev_sum==8,1,0))

## check missing values
table(dev_15mon$dev_sum,useNA="always")


## developmental milestones 24 months
dev_24mon <- read_excel("EC_RCT_24_months_20210304-080857.xlsx",
                        sheet=9)

## phone surveys
dev_24mon_phone <- read_excel("24_months_phone_survey_completed_20210304-091337.xlsx",
                              sheet=8)

## join with identifier dataset
dev_24mon <- dev_24mon %>% inner_join(dat_24mon_identifier)
dev_24mon_phone <- inner_join(dev_24mon_phone,dat_24mon_phone_identifier)
dev_24mon_phone$`Current clinic` <- NA

## join together
dev_24mon <- rbind(dev_24mon,dev_24mon_phone)
dev_24mon <- dev_24mon %>% mutate(pid=`PID check`)

## add baseline info
dev_24mon <- left_join(dev_24mon,dev_info)

## remove duplicates
dev_24mon <- unique(dev_24mon)

## remove missing values
dev_24mon <- dev_24mon %>%
  filter(`WHO 1 SWS criteria`!="-")

## calculate the sum
dev_24mon <- dev_24mon %>%
  mutate(dev1=ifelse(`WHO 1 SWS criteria`=="3",1,0),
         dev2=ifelse(`WHO 2 HKC criteria`=="3",1,0),
         dev3=ifelse(`WHO3 SWA criteria`=="3",1,0),
         dev4=ifelse(`WHO 4 WWA criteria`=="3",1,0),
         dev5=ifelse(`WHO5 SA criteria`=="3",1,0),
         dev6=ifelse(`WHO6 WA criteria`=="3",1,0),
         dev7=ifelse(Speech=="3",1,0),
         dev8=ifelse(`RTHC Social_Hearing`=="1",1,0),
         dev_sum=dev1+dev2+dev3+dev4+dev5+dev6+dev7+dev8)

## check missing values/check missing values
table(dev_24mon$dev_sum,useNA="always")

## combining both of these
dev_6mon_small <- dev_6mon %>%
  mutate(dev_sum_6mon=dev_sum) %>%
  select(pid,intv,clinic_baseline,dev_sum_6mon)

dev_15mon_small <- dev_15mon %>% 
  mutate(dev_sum_15mon=dev_sum) %>%
  select(pid,intv,clinic_baseline,dev_sum_15mon)

dev_24mon_small <- dev_24mon %>% 
  mutate(dev_sum_24mon=dev_sum) %>%
  select(pid,intv,clinic_baseline,dev_sum_24mon)

## join them
dev_dat <- inner_join(dev_6mon_small,
                      dev_15mon_small)

dev_dat <- inner_join(dev_dat,
                      dev_24mon_small)

## create sum variable and then truncate for values under 10
dev_dat<- dev_dat %>%
  mutate(sum_all=dev_sum_6mon+dev_sum_15mon+dev_sum_24mon) %>%
  mutate(sum_all=ifelse(sum_all<10,10,sum_all))

## means for table 
dev_dat %>% 
  group_by(intv) %>% 
  summarise(n=n(),
            mean=mean(sum_all),
            sd=sd(sum_all)) %>%
  kable(digits=1)

## developmental milestones model
dev_model <- blmer(sum_all~(1|clinic_baseline)+intv,
                   data=dev_dat)
summary(dev_model)

hist(residuals(dev_model))
qqnorm(residuals(dev_model));qqline(residuals(dev_model))

## asy covariance matrix
dev_cov <- vcov(dev_model)

## fixed effect vector
dev_fixed <- fixed.effects(dev_model)

## all test
dev_contrast <- c(1,1)-c(1,0)
dev_est <- dev_contrast %*% dev_fixed
dev_se <- sqrt(t(dev_contrast) %*% (dev_cov %*% dev_contrast))
dev_Z <- dev_est/dev_se
round(1-pnorm(dev_Z[1,1]),3) ## 0.068

## save smaller dataset for correlation
dev_corr <- dev_dat %>% select(pid,sum_all) %>% na.omit()

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## Calculate correlation between outcomes

## join all the datasets together
corr_dat <- full_join(bf_corr,lbw_corr)
corr_dat <- full_join(corr_dat,alc_corr)
corr_dat <- full_join(corr_dat,everdep_corr)
corr_dat <- full_join(corr_dat,haz_corr)
corr_dat <- full_join(corr_dat,waz_corr)
corr_dat <- full_join(corr_dat,csg_corr)
corr_dat <- full_join(corr_dat,immune_corr)
corr_dat <- full_join(corr_dat,hos_corr)
corr_dat <- full_join(corr_dat,antenatal_corr)
corr_dat <- full_join(corr_dat,child_adher_corr)
corr_dat <- full_join(corr_dat,arv_corr)
corr_dat <- full_join(corr_dat,dev_corr)

## correct so all variables are coded so higher values are better
corr_dat <- corr_dat %>%
  mutate(not_lbw=1-lbw_ind,
         not_dep=1-ever_dep,
         no_haz=1-ever_haz_lt2,
         no_waz=1-ever_waz_lt2,
         no_hos=1-hos_ever) %>%
  select(pid,bf_upto6mons,not_lbw,noalc_preg,not_dep,no_haz,no_waz,csg,always_upto,
         no_hos,antenatal_gt4,child_sum,arv_sum,sum_all)

## loop through all variables
corr_mat <- matrix(NA,nrow=13,ncol=13)
tet_mat <- matrix(NA,nrow=13,ncol=13)

## is_binary function
is_binary <- function(x){
  unique_x <- unique(x)
  return(length(unique_x)==2)
}

for(i in 1:13){
  for(j in 1:13){
    
    temp <- corr_dat[,c(1,i+1,j+1)]
    temp <- na.omit(temp)
    temp <- as.matrix(temp)
    
    pearson_corr <- cor(temp[,2],temp[,3])
    
    a_binary <- is_binary(temp[,2])
    b_binary <- is_binary(temp[,3])
    
    if(a_binary & b_binary){
      tet_corr <- tetrachoric(table(temp[,2],temp[,3]))
      tet_corr <- tet_corr[1]$rho
      corr_mat[i,j] <- max(abs(tet_corr),abs(pearson_corr))
    }else{
      corr_mat[i,j] <- pearson_corr[1]
    }
    
  }
}

corrs <- as.vector(corr_mat[upper.tri(corr_mat)])
mean(abs(corrs)) ## 0.108 --> use 0.2 as correlation per the text

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

## run z-scores together
Z_scores <- c(bf_Z_6mon[1,1],lbw_Z[1,1], alc_Z[1,1], everdep_Z[1,1], haz_Z[1,1], waz_Z[1,1], csg_Z[1,1], 
              immune_Z[1,1], hos_Z[1,1], antenatal_Z[1,1], child_adher_Z[1,1], arv_Z[1,1], dev_Z[1,1])

contrast <- rep(1,13)
num <- contrast %*% Z_scores

## total variable
total_cov <- matrix(0.1,13,13)
diag(total_cov) <- rep(1,13)

## different standard errors based on covariance matrix
se1 <- sqrt(t(contrast) %*% total_cov %*% contrast)
se2 <- sqrt(t(contrast) %*% corr_mat %*% contrast)

## test statistics
Z_stat1 <- num/se1
Z_stat2 <- num/se2

## p-values
1-pnorm(Z_stat1[1,1])
1-pnorm(Z_stat2[1,1])



Z_scores>0
sum(Z_scores>0)

## using a sign test
## 0.011
dbinom(11,13,0.5)+dbinom(12,13,0.5)+dbinom(13,13,0.5)

