## -------------------------------------------------------------------------- ##
## table1.R ----------------------------------------------------------------- ##
## Author: Peter Norwood, UCLA ---------------------------------------------- ##
## Purpose: Table 1 for the 24 month ECSS paper ----------------------------- ##
## -------------------------------------------------------------------------- ##

## load packages
library(tidyverse)
library(lme4)
library(readxl)
library(patchwork)
library(blme)
library(knitr)

## set working directory
setwd("~/ECSS_24Mon_Paper")

## expit function
expit <- function(x){
  exp(x)/(1+exp(x))
}
expit <- Vectorize(expit)

## load data and basic data cleaning
dat <- read_excel("ecss_complete_dataset_aug2020.xlsx")

## remove missing clinics
## create month variable
dat <- dat %>%
  filter(clinic!=99) %>%
  mutate(month=case_when(time==1~0,
                         time==2~3,
                         time==3~6),
         clinic=as.factor(clinic),
         intv_char=ifelse(intv==1,"intv","control"))

## baseline table
baseline <- dat %>% filter(time==1)

## totals for each group
table(baseline$intv,useNA="always")

## LIST OF VARIABLES

## -------------------------------------------------------------------------- ##
# Age

## check missing values
table(baseline$age,useNA="always")

## by group
baseline %>% 
  mutate(age=as.numeric(age)) %>%
  group_by(intv) %>%
  summarise(mean=mean(age),
            sd=sd(age)) %>%
  kable(digits=2)

## overall
baseline %>% 
  mutate(age=as.numeric(age)) %>%
  summarise(mean=mean(age),
            sd=sd(age)) %>%
  kable(digits=2)

## model testing difference
age_fit <- blmer(as.numeric(age)~(1|clinic)+intv,data=baseline)
## act as if it's a z-test
sum_age_fit <- summary(age_fit)
2*(1-pnorm(abs(sum_age_fit$coefficients[2,3]))) ## 0.5255542

## -------------------------------------------------------------------------- ##
# Pregnancy recruit or post-birth

## -------------------------------------------------------------------------- ##
# Education

## check missing values
table(baseline$education,useNA="always")

## by group
baseline %>% 
  mutate(education=as.numeric(education)) %>%
  group_by(intv) %>%
  summarise(mean=mean(education),
            sd=sd(education)) %>%
  kable(digits=2)

## overall
baseline %>% 
  mutate(education=as.numeric(education)) %>%
  summarise(mean=mean(education),
            sd=sd(education)) %>%
  kable(digits=2)

## model testing difference
edu_fit <- blmer(as.numeric(education)~(1|clinic)+intv,data=baseline)
## act as if it's a z-test
sum_edu_fit <- summary(edu_fit)
2*(1-pnorm(abs(sum_edu_fit$coefficients[2,3]))) ## 0.012

## -------------------------------------------------------------------------- ##
# Employment

## employment using karl's question
## by group
baseline <- baseline %>% mutate(employ_ind=ifelse(participantanyemployment=="5",0,1))
table(baseline$employ_ind,useNA="always")

baseline %>% 
  group_by(intv) %>%
  summarise(n_emp=sum(employ_ind),
            n_total=n(),
            prop=100*mean(employ_ind)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n_emp=sum(employ_ind),
            n_total=n(),
            prop=100*mean(employ_ind)) %>%
  kable(digits=2)

## model testing difference
emp_fit <- bglmer(employ_ind~(1|clinic)+intv,family=binomial(),
                  data=baseline)

## z test
summary(emp_fit) ## 0.0236


## -------------------------------------------------------------------------- ##
# Household income

## -------------------------------------------------------------------------- ##
# Child grants

## check missing values
table(baseline$intv,baseline$csg,useNA="always")

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(csg),
            prop=100*mean(csg)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(csg),
            prop=100*mean(csg)) %>%
  kable(digits=2)

## model testing difference
csg_fit <- bglmer(as.numeric(csg)~(1|clinic)+intv,family=binomial(),
                  data=baseline)

## z test
summary(csg_fit) ## 0.111

## -------------------------------------------------------------------------- ##
# Previous pregnancy

## check missing values
table(baseline$intv,baseline$numberofpregnancies,useNA="always")

## create yes/no variable
baseline <- baseline %>%
  mutate(prev_preg=ifelse(numberofpregnancies==0,0,1))

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(prev_preg),
            prop=100*mean(prev_preg)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(prev_preg),
            prop=100*mean(prev_preg)) %>%
  kable(digits=2)

## model testing difference
preg_fit <- bglmer(prev_preg~(1|clinic)+intv,family=binomial(),
                  data=baseline)

## z test
summary(preg_fit) ## 0.111

## -------------------------------------------------------------------------- ##

# Water on premises

## check missing values
table(baseline$intv,baseline$sourceofdrinkingwater,useNA="always")

## create yes/no variable
baseline <- baseline %>%
  mutate(wonsite=ifelse(sourceofdrinkingwater %in% c("1","2"),1,0))

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(wonsite),
            prop=100*mean(wonsite)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(wonsite),
            prop=100*mean(wonsite)) %>%
  kable(digits=2)

## model testing difference
wonsite_fit <- bglmer(wonsite~(1|clinic)+intv,family=binomial(),
                   data=baseline)

## z test
summary(wonsite_fit) ## 0.111

## -------------------------------------------------------------------------- ##

# flush Toilet

## check missing values
table(baseline$intv,baseline$maintypeoftoiletfacility,useNA="always")

## create yes/no variable
baseline <- baseline %>%
  mutate(flush_toilet=ifelse(maintypeoftoiletfacility %in% c("1","2"),1,0))

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(flush_toilet),
            prop=100*mean(flush_toilet)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(flush_toilet),
            prop=100*mean(flush_toilet)) %>%
  kable(digits=2)

## model testing difference
wonsite_fit <- bglmer(wonsite~(1|clinic)+intv,family=binomial(),
                      data=baseline)

## z test
summary(wonsite_fit) ## 0.111

## -------------------------------------------------------------------------- ##
# Electricity

## -------------------------------------------------------------------------- ##
# Distance to health care if it has it. We have distances charted other places

## -------------------------------------------------------------------------- ##
# Depression scale score

## check missing values
table(baseline$intv,baseline$epds,useNA="always")

## by group
baseline %>% 
  mutate(epds=as.numeric(epds)) %>%
  group_by(intv) %>%
  summarise(mean=mean(epds),
            sd=sd(epds)) %>%
  kable(digits=2)

## overall
baseline %>% 
  mutate(epds=as.numeric(epds)) %>%
  summarise(mean=mean(epds),
            sd=sd(epds)) %>%
  kable(digits=2)

## model testing difference
epds_fit <- blmer(as.numeric(epds)~(1|clinic)+intv,data=baseline)
## act as if it's a z-test
sum_epds_fit <- summary(epds_fit)
2*(1-pnorm(abs(sum_epds_fit$coefficients[2,3]))) ## 0.03695627

## -------------------------------------------------------------------------- ##
# Depression > 13

## create new variable
baseline <- baseline %>%
  mutate(epds=as.numeric(epds),
         epds_gt13=ifelse(epds>13,1,0))

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(epds_gt13),
            prop=100*mean(epds_gt13)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(epds_gt13),
            prop=100*mean(epds_gt13)) %>%
  kable(digits=2)

## model testing difference
epds13_fit <- bglmer(epds_gt13~(1|clinic)+intv,family=binomial(),
                      data=baseline)

## z test
summary(epds13_fit) ## 0.955
## -------------------------------------------------------------------------- ##
# Alcohol use/not

## -------------------------------------------------------------------------- ##
# Interpersonal partner violence (ever)

## recode variable
baseline <- baseline %>%
  mutate(ipv_ever=ifelse(violenceever=="1",1,0),
         ipv_current=ifelse(violencecurrentpartner=="1",1,0))

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(ipv_ever),
            prop=100*mean(ipv_ever)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(ipv_ever),
            prop=100*mean(ipv_ever)) %>%
  kable(digits=2)

## model testing difference
ipv_ever_fit <- bglmer(ipv_ever~(1|clinic)+intv,family=binomial(),
                     data=baseline)

## z test
summary(ipv_ever_fit) ## 0.830

## -------------------------------------------------------------------------- ##

## violence current partner

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(ipv_current),
            prop=100*mean(ipv_current)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(ipv_current),
            prop=100*mean(ipv_current)) %>%
  kable(digits=2)

## model testing difference
ipv_current_fit <- bglmer(ipv_current~(1|clinic)+intv,family=binomial(),
                       data=baseline)

## z test
summary(ipv_current_fit) ## 0.862

## -------------------------------------------------------------------------- ##
# other household members

## create new variable
baseline <- baseline %>%
  mutate(numberofadultmen=as.numeric(numberofadultmen),
         numberofadultwomen=as.numeric(numberofadultwomen),
         numberofboys=as.numeric(numberofboys),
         numberofgirls=as.numeric(numberofgirls),
         num_household=numberofadultmen+numberofadultwomen+numberofboys+numberofgirls,
         other_household=ifelse(num_household>0,1,0))


table(baseline$intv,baseline$num_household,useNA="always")

## by group
baseline %>% 
  filter(!is.na(num_household)) %>%
  group_by(intv) %>%
  summarise(mean=mean(num_household),
            sd=sd(num_household)) %>%
  kable(digits=2)

## overall
baseline %>% 
  filter(!is.na(num_household)) %>%
  summarise(mean=mean(num_household),
            sd=sd(num_household)) %>%
  kable(digits=2)

## model testing difference
num_household_fit <- blmer(num_household~(1|clinic)+intv,data=baseline)
## act as if it's a z-test
sum_num_fit <- summary(num_household_fit)
2*(1-pnorm(abs(sum_num_fit$coefficients[2,3]))) ## 0.8897593

## -------------------------------------------------------------------------- ##
# Mother in charge of childcare or mother left for Western Cape

## -------------------------------------------------------------------------- ##
# Father in household

## recode variable
baseline <- baseline %>%
  mutate(father_household=ifelse(fatherstayingwithyou=="1",1,0))

## by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(father_household),
            prop=100*mean(father_household)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(father_household),
            prop=100*mean(father_household)) %>%
  kable(digits=2)

## model testing difference
father_fit <- bglmer(father_household~(1|clinic)+intv,family=binomial(),
                          data=baseline)

## z test
summary(father_fit) ## 0.862

## -------------------------------------------------------------------------- ##
# alcohol before/after and after pregnancy
  
## recode variable
baseline <- baseline %>%
  mutate(alc_before=ifelse(alcoholinpregnancybeforeknow=="1",0,1),
         alc_after=ifelse(alcoholduringpregnancyafterl=="1",0,1))

## before - by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(alc_before),
            prop=100*mean(alc_before)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(alc_before),
            prop=100*mean(alc_before)) %>%
  kable(digits=2)

## model testing difference
alc_before_fit <- bglmer(alc_before~(1|clinic)+intv,family=binomial(),
                     data=baseline)

## z test
summary(alc_before_fit) ## 0.706

## after - by group
baseline %>% 
  group_by(intv) %>%
  summarise(n=sum(alc_after),
            prop=100*mean(alc_after)) %>%
  kable(digits=2)

## overall
baseline %>% 
  summarise(n=sum(alc_after),
            prop=100*mean(alc_after)) %>%
  kable(digits=2)

## model testing difference
alc_after_fit <- bglmer(alc_after~(1|clinic)+intv,family=binomial(),
                         data=baseline)

## z test
summary(alc_after_fit) ## 0.494 
