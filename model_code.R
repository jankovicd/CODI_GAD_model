
library(BCEA)
library(SHELF)
library(reshape2)

rm(list=ls(all=TRUE))                                    

nums<-function (x) (as.numeric(as.character((x))))

t_names <- c("NI", "UDC", "SDC", "UDI", "SDI", "SNoDI", "M")
n_treat <- length(t_names)

setwd("//.../GitHub model")

coda_folder <- data_folder <- "inputs"
tabl_folder <- "Results - tables/"
figs_folder <- "Results - figures/"

########### import GAD-7 scores ###########

# post-treatment GAD-7 scores
gad7_scr <- read.csv(file=paste0(coda_folder,"/coda.csv",sep=""), header=TRUE)[-c(1:4),2:8]

colnames(gad7_scr)<-t_names

n_samps<-nrow(gad7_scr) # number of random samples (number of iterations in meta-analysis)

######## model + pop specification #########

age <- 47.6                 # age point estimate. In Revicki mean age was 47.6, sd 13.7
sex <- 2                    # single sex (female). In Revicki 72.4% female.
n_years <- 101-floor(age)   # model time horizon
n_cycle <- n_years*4        # number of cycles
dr <- 1.035                 # discount rate

# gad severity states by gad-7 score
state1 <- 0: 4 
state2 <- 5: 9
state3 <-10:14
state4 <-15:21

#each scenario is ran separately
b_assum<-2 #1 assumes patients in baseline are constant 
#2 assumes baseline goes down (base case)
t_assum<-1 #1 assumes treatment effect is constant indefinitely (base case)
#2 assumes treatment effect drops off after 1 year
#3 assumes treatment effect drops gradually between 1 and 10 years

########### set costs and qalys ###########

######## utilities ######

uts1 <- matrix(c(0.72,0.1,47,0.64,0.1,101,0.60,0.1,102,0.53,0.1,46),3,4) #utility, sd and sample size reported in Revicki (by severity level)
uts2 <- t(matrix(c(uts1[1,],c((uts1[2,]/sqrt(uts1[3,]))^2)),4,2)) #derive beta parameters

falfa <- function (m,v) ((((1-m)/v)-(1/m))*m^2) #derive beta parameter 1 from mean and sd
fbeta <- function (a,m) (a*(1/m-1)) #derive beta parameter 2 from mean and sd

qaly_mat <- array(NA,c(n_samps,4,n_treat),list(1:n_samps, 1:4, colnames(gad7_scr)))

# random utility score for each state
for (i in 1:4){
  alfa<-falfa(uts2[1,i],uts2[2,i])
  beta<-fbeta(alfa,uts2[1,i])
  
  for (j in 1:n_treat){
    qaly_mat[,i,j]<-rbeta(n_samps,alfa, beta)/4
  }
}

#age-related qaly decrements
qaly_dec<-read.csv(file=paste(data_folder,"/utility_decrements.csv",sep=""), header=TRUE)

qaly_age<-qaly_rda<-array(NA,c(n_samps,101,n_treat)) # matrix with QALYs and QALY decrements for different ages
for(j in 1:n_treat){ qaly_age[,1,j] <- rnorm(n_samps, qaly_dec[1,3], qaly_dec[1,4]) }
qaly_rda[,1,] <- 0

for (a in 2:101){
  for(b in 1:n_treat){
    qaly_age[,a,b] <- rnorm(n_samps, qaly_dec[a,3], qaly_dec[a,4])
    qaly_rda[,a,b] <- qaly_age[,a,b] - qaly_age[,a-1,b]
  }
}

rm(alfa, beta, falfa, fbeta, uts1, uts2, qaly_dec, qaly_age)

######## costs #######

#import inflation rates
cpi05 <- read.csv(file=paste(data_folder,"/CPI05.csv",sep=""), header=FALSE)[-(1:7),]   
cpi09 <- read.csv(file=paste(data_folder,"/CPI.csv",sep=""), header=FALSE)[-(1:7),c(1,3)]
cpi05 <- cpi05[grep("May",cpi05[,1]),]
cpi09 <- cpi09[grep("May",cpi09[,1]),]
colnames(cpi05) <- colnames(cpi09) <- c("year","cpi")
cpi <- rbind(cpi05,cpi09[grep("16", cpi09[,1]):nrow(cpi09),]); rm(cpi05, cpi09)

sel_cost <- 4 #1 for costs from VL, 2 for costs from Kumar, 3 for no state cost, 4 for McCrone (base case)

cost_psychol <- 53 #cost of psychologist's time

# unit cost for resources reported in Vera-Llonch
uc_vl <- c(31, # per 9.22 minute appointment (£204 per hour incl. direct care). gp. £39 per hour for nurse. £59 per hour for nurse f2f. 
         109, # cost of psychiatrist per hour,
         cost_psychol, # cost of clinical psychologist per hour (band 7),
         222, # Emergency room average treatment cost, assuming patient arrives alive.
         #VB01Z	Emergency Medicine, Any Investigation with Category 5 Treatment	338 
         #VB02Z	Emergency Medicine, Category 3 Investigation with Category 4 Treatment	338 
         #VB03Z	Emergency Medicine, Category 3 Investigation with Category 1-3 Treatment	252 
         #VB04Z	Emergency Medicine, Category 2 Investigation with Category 4 Treatment	227 
         #VB05Z	Emergency Medicine, Category 2 Investigation with Category 3 Treatment	184 
         #VB06Z	Emergency Medicine, Category 1 Investigation with Category 3-4 Treatment	130 
         #VB07Z	Emergency Medicine, Category 2 Investigation with Category 2 Treatment	163 
         #VB08Z	Emergency Medicine, Category 2 Investigation with Category 1 Treatment	155 
         #VB09Z	Emergency Medicine, Category 1 Investigation with Category 1-2 Treatment	106 
         0,#other
         6*cumprod(1+nums(cpi[grep("13",cpi[,1]):nrow(cpi),2])/100)[7], 
         #6.909099, #6(2013/14 price)*1.016*1.012*1.025*1.041*1.027*1.022 - blood count (2012 prices, add inflation)
         58, #code RD51A
         1.837497,#thyroid function test in 2008/09 adjusted 1.34*1.024*1.048*1.048*1.031*1.027*1.016*1.012*1.025*1.041*1.027*1.022
         1603)#inpatient days 2017/18 tariff

#unit costs for resources reported in Kumar
uc_ku <- c(31, # per 9.22 minute appointment (£204 per hour incl. direct care). gp. £39 per hour for nurse. £59 per hour for nurse f2f.
         222, # Emergency room average treatment cost, assuming patient arrives alive.
         1603)#inpatient days

# unit costs for resources reported in McChrone
uc_mc <- t(matrix(c(1.96,0.016,
                  0.848,0.003,
                  0.295,0.001,
                  0.32,0.0007),2,4))

#Resource use in Kumar. Kumar reports confidence intervals, so fit gamma for use in sensitivity analysis
ru_gamma <- rbind(fitdist(c(1.1,1.2,1.3), c(0.025,0.5,0.975), lower=0)$Gamma,
                fitdist(c(1.5,1.7,1.9), c(0.025,0.5,0.975), lower=0)$Gamma,
                fitdist(c(1.9,2.2,2.5), c(0.025,0.5,0.975), lower=0)$Gamma,
                fitdist(c(2.0,2.4,2.8), c(0.025,0.5,0.975), lower=0)$Gamma)

#Primary care visits (mean) (per month) in Vera Llonch
ru_vl <- t(matrix(c(0.44, 1.03, 1.26, 1.80,
                  ##Specialist visits (mean) (per month)
                  #Psychiatrist
                  0.42, 0.48, 0.48, 0.49,
                  #Psychologist
                  0.48, 0.52, 1.03, 1.37,
                  #Emergency room
                  0.14, 0.26, 0.37, 0.56,
                  #Other
                  0.33, 0.37, 0.58, 0.52,
                  ##Other outpatient services (mean) (per month)
                  #Blood counts
                  0.35, 0.38, 0.50, 0.43,
                  #Electrocardiogram
                  0.33, 0.35, 0.33, 0.18,
                  #Thyroid function
                  0.33, 0.33, 0.36, 0.35,
                  #Inpatient days (mean)
                  0.12, 0.18, 0.37, 0.49),4,9))

#resource use in Kumar
ru_ku <- t(matrix(c(1.2, 1.7, 2.2, 2.4,
                  0.014,0.019,0.025,0.027,
                  0.014,0.019,0.025,0.027),4,3))


uc <- list(vl=uc_vl,ku=uc_ku,no=0,mc=uc_mc)[[sel_cost]]
ru <- list(vl=ru_vl,ku=ru_ku,no=0,mc=0)[[sel_cost]]

#mean(apply(ru*uc,2,sum))

cost_mat <- array(NA, c(n_samps, 4, n_treat), list(1:n_samps, 1:4, colnames(gad7_scr)))

for (t in 1:n_treat){#done separately for each treatment to ensure random samples are slightly different
  
  if(sel_cost <3 ){
    ruse_mat <- array(NA,c(n_samps, nrow(ru), ncol(ru)), list(1:n_samps, 1:nrow(ru), 1:ncol(ru)))#sample resource use
  }
  
  if(sel_cost == 1){  
    
    for (i in 1:ncol(ru)){
      
      for(j in 1:nrow(ru)){ ruse_mat[,j,i] <- runif(n_samps, ru[j,i] * 0.75, ru[j,i] * 1.25) * 3 } #uniform +/- 25% for now
      
      cost_mat[,i,t] <- ruse_mat[,1,i]*uc[1] + ruse_mat[,2,i]*uc[2] + ruse_mat[,3,i]*uc[3] +
        ruse_mat[,4,i]*uc[4] + ruse_mat[,5,i]*uc[5] + ruse_mat[,6,i]*uc[6] +
        ruse_mat[,7,i]*uc[7] + ruse_mat[,8,i]*uc[8] + ruse_mat[,9,i]*uc[9]     
    }
    
    
  } else if(sel_cost==2) {
    
    for (i in 1:ncol(ru)){
      
      for(j in 1:nrow(ru)){ruse_mat[,j,i]<-runif(n_samps,ru[j,i]*0.75, ru[j,i]*1.25)*3} #uniform +/- 25% for now
      
      cost_mat[,i,t]<-ruse_mat[,1,i]*uc[1] + ruse_mat[,2,i]*uc[2] + ruse_mat[,3,i]*uc[3]
      
    }
    
  } else if(sel_cost==3){
    
    cost_mat<-array(0,c(n_samps,4,n_treat),list(1:n_samps, 1:4, colnames(gad7_scr)))
    
  } else if(sel_cost==4){
    
    for (i in 1:4){
      cost_mat[,i,t]<-(rgamma(n_samps,uc[i,1],uc[i,2])/2)*cumprod(1+nums(cpi[grep("05",cpi[,1]):nrow(cpi),2])/100)[15]
    }
  }
  
  rm(ruse_mat)
  
}

cost_int <- matrix(NA,n_samps,n_treat)
cost_dig <- 0 # cost of the digital psychological intervention. Explore other costs.
cost_con <- 1.5*cost_psychol # cost of contact in digital interventions (90min of contact with psychologist)
n_drugs <- c(1547338,1225544,601118,114328,113647)#number of prescriptions issued in january 2020 per type of SSRI - from opensource prescribing
#sertraline, citalopram, fluoxetine, escitalopram, paroxetine
c_drugs <- c(18.40,12.35,13.13,(43.94+56.29)/2,16.9)#cost for 13 x 28 tablets as per drug tariff
#sertraline (mean for 50mg and 100mg), citalopram (20mg), fluoxetine (20mg), escitalopram (mean for 10mg and 20mg), paroxetine 20mg
cost_dr <- sum(c_drugs*n_drugs/sum(n_drugs))
cost_di <- 12*1.26 #dispensed 13 times per year
cost_gp <- 42.6

cost_int[,1] <- 0
cost_int[,2] <- 0
cost_int[,3] <- 0 + cost_psychol/3 #assume 20min to deliver
cost_int[,4] <- cost_dig
cost_int[,5] <- cost_dig + 1.5*cost_psychol #assume 1.5 hours to deliver
cost_int[,6] <- 1.5*cost_psychol # 1.5 hours with clinical psychologist (band 7) per person (group therapy)
cost_int[,7] <- cost_dr + cost_di + cost_gp*7.5 #meds in first year of treatment based on number of GP visits, number of Rxs and cost of meds.

cost_med <- c(0, ((cost_dr + cost_di+cost_gp*4) * (1/(1.035^c(0:3)))))#cost in first 5 years, discounted. Assumes treatment will last 5 years. cost in years 2:5 lower as less monitoring required.

rm(ru, uc, cost_dig, cost_con, cost_dr, cost_di, cost_gp)

######## mortality rates ########

mort_ons <- read.csv(file=paste(data_folder,"/mort_data_ons.csv",sep=""), header=FALSE)   #### isnt there better data than this??? THE THING OF HAVING A PROBABILITY OF 1 OF DYING FROM 85 ABOVE IS A BIT WEIRD
mort_ons <- mort_ons[8:nrow(mort_ons), c(1, 3, 9)]
rownames(mort_ons) <- 0:(nrow(mort_ons) - 1)
colnames(mort_ons) <- c("lower_age", "male", "female")
mort_all <- matrix(c(1:nrow(mort_ons), rep(NA,2*nrow(mort_ons))), nrow(mort_ons), 3)
mort_all[,2] <- 1 - exp(log(1 - nums(mort_ons[,2]))*0.25)  ### 1YEAR PROBABILITY TO RATE AND BACK TO 3 MONTH PROBABILITY
mort_all[,3] <- 1 - exp(log(1 - nums(mort_ons[,3]))*0.25)
rm(mort_ons)

mort_cyc <- matrix(NA, n_samps, n_cycle)#mortality per cycle
mort_cyc[,1] <- 0
for (i in 2:n_cycle){
  if(sex==1){
    mort_cyc[,i] <- mort_all[floor(age + (i-1)/4),2]
  }else{
    mort_cyc[,i] <- mort_all[floor(age + (i-1)/4),3]
  }
}

# add excess mortality
mortex_mn <- c(0.11, 0.13, 0.172, 0.236) #mean from study    #### WHICH STUDY??
mortex_ns <- c(688, 315, 203, 144) #n from study
mortex_as <- mortex_mn*mortex_ns #alpha parameter for the beta distribution
mortex_bs <- mortex_ns-mortex_as #beta parameter for the beta distribution

mort_exs <- array(NA, c(n_samps, 4, n_treat))
for (j in 1:n_treat){
  temp<-rbeta(n_samps, mortex_as[1], mortex_bs[1])
  mort_exs[,1,j] <- temp/temp
  mort_exs[,2,j] <- rbeta(n_samps,mortex_as[2], mortex_bs[2])/temp
  mort_exs[,3,j] <- rbeta(n_samps,mortex_as[3], mortex_bs[3])/temp
  mort_exs[,4,j] <- rbeta(n_samps,mortex_as[4], mortex_bs[4])/temp
}

rm(mortex_mn, mortex_ns, mortex_as, mortex_bs, temp)

########### graph parameters ##########

lnwds <- c(1,1,2,1,2,2,3)
lntys <- c(1,3,3,5,5,1,1)
# lnwds<-c(2,1,1,1,1,2,2)
# lntys<-c(1,1,2,3,4,2,3)
x_yrs <- (0:(n_cycle-1))/4

f_gph_name<-function(a)
  (paste0(figs_folder,a,"_b",b_assum,"t",t_assum,".jpeg"))
f_tab_name<-function(a)
  (paste0(tabl_folder,a,"_b",b_assum,"t",t_assum,".csv"))

########### run the model (predict GAD-7 scores) ############

#if running all scenarios at once, remove # in the next 4 lines
#rm(b_assum, t_assum)
#
#for(b_assum in 1:2){
#  for(t_assum in 1:3){

counter<-array(NA,c(n_samps,n_cycle,n_treat),list(1:n_samps, 1:n_cycle, colnames(gad7_scr)))
counter[,1,]<-as.matrix(ceiling(nums(gad7_scr[,1])))
for(t in 1:n_treat){
  counter[,2,t]<-as.matrix(ceiling(nums(gad7_scr[,t])))
}

for (t in 1:n_treat){counter[,3:n_cycle,t]<-counter[,2,t] }

if (b_assum==2){
  
    samplen<-sample(which(nums(counter[,2,1])>4),round(n_samps*0.3,digits=0)) # sample random 30% of simulations
    #select 50% of the random sample above because 15% improve in the first year
    sample1<-samplen[1:((length(samplen)*0.5))] # 15% of simulations (50% of the random sample above) improve after one year
    #select further 1/3 of the random sample because 10% improve in the second year
    sample2<-samplen[(length(sample1)+1):round(length(samplen)*(0.5+1/3),digits=0)] # 10% of simulations (33% of the random sample above) improve after two years
    #select the rest of the random sample samplen because 5% improve in the third year
    sample3<-samplen[(length(sample1)+length(sample2)+1):(length(samplen))] # 5% of simulations (16.7% of the random sample above) improve after THREE years
    counter[sample1, 5:n_cycle,1]<-nums(counter[sample1,4,1])-5  # move patients into lower state by reducing GAD7 score by 5. 
    counter[sample2, 9:n_cycle,1]<-nums(counter[sample2,8,1])-5  # move patients into lower state by reducing GAD7 score by 5. 
    counter[sample3,13:n_cycle,1]<-nums(counter[sample3,12,1])-5 # move patients into lower state by reducing GAD7 score by 5.
    
    for (t in 2:n_treat){

      counter[sample1, 5:n_cycle,t]<-ifelse(counter[sample1, 5:n_cycle,t]>4,nums(counter[sample1, 4,t])-5,counter[sample1, 5:n_cycle,t])  # move patients into lower state by reducing GAD7 score by 5.
      counter[sample2, 9:n_cycle,t]<-ifelse(counter[sample2, 9:n_cycle,t]>4,nums(counter[sample2, 8,t])-5,counter[sample2, 9:n_cycle,t])  # move patients into lower state by reducing GAD7 score by 5.
      counter[sample3,13:n_cycle,t]<-ifelse(counter[sample3,13:n_cycle,t]>4,nums(counter[sample3,12,t])-5,counter[sample3,13:n_cycle,t]) # move patients into lower state by reducing GAD7 score by 5.

    }

    rm(samplen,sample1,sample2,sample3)
        
  }

if (t_assum==2){
  for (t in 2:n_treat){
    counter[,5:n_cycle,t]<-counter[,5:n_cycle,1]
    }
  } else if (t_assum==3){
    for (t in 2:n_treat){
      for (i in 5:40){
        diffs<-(nums(counter[,i,1])-nums(counter[,i,t]))/36
        counter[,i,t]<-nums(counter[,i,t])+(diffs*(i-4))
      }
      counter[,41:n_cycle,t]<-counter[,41:n_cycle,1]
      rm(diffs)
  }
}

#sum(nums(unlist(counter[,2,]))<0)#count gad scores below 0 after 1st cycle
#sum(nums(unlist(counter[,n_cycle,]))<0)#count gad scores below 0 after last cycle

########### plot outcomes over time ##########

if(t_assum==1){

#count patients in each state, after each cycle
prob_sta<-array(NA,c(4,n_cycle,n_treat),
                list(c("No anxiety", "Mild anxiety", "Moderate anxiety", "Severe anxiety"),1:n_cycle,t_names))
prob_sta[1,,]<-apply(counter,c(2,3),function (x) mean(x %in% state1))  # NOT SURE IF WE SHOULD USE THE MID POINT OR THE END OF DISTRIBUTION
prob_sta[2,,]<-apply(counter,c(2,3),function (x) mean(x %in% state2))
prob_sta[3,,]<-apply(counter,c(2,3),function (x) mean(x %in% state3))
prob_sta[4,,]<-apply(counter,c(2,3),function (x) mean(x %in% state4))

write.csv(prob_sta[,2,], file =paste0(tabl_folder,"cycle2_proportion_in_each_state_b",b_assum,"t",t_assum,".csv"), row.names=FALSE)

#graph of proportion in each state with no intervention
jpeg(f_gph_name("basline_states_over_time"), width = 600, height = 370)  # bland Altman plot
par(mar=c(4,4,0.2,0.2))
plot(x_yrs,prob_sta[1,,1],
     #ylim=c(0,0.9),xlim=c(0,10),#full sample
     ylim=c(0,1),xlim=c(0,10),#subgroup analysis
     ylab="Proportion in state",xlab="Time (years)",
     type="l",lty=2)
lines(x_yrs,prob_sta[2,,1],lty=2,lwd=2)
lines(x_yrs,prob_sta[3,,1],lwd=1)
lines(x_yrs,prob_sta[4,,1],lwd=2)
legend("topleft", c("No Anxiety","Mild", "Moderate", "Severe"),lty=c(2,2,1,1), lwd=c(1,2,1,2), horiz=TRUE, bty = "n")#subgroup analysis
dev.off()

}
 
gad_mean<-apply(counter, c(2,3), mean)
 
jpeg(f_gph_name("mean_gad7_over_time"), width = 600, height = 370)  
par(mar=c(0.2,4,0.2,0.2))
plot(x_yrs,nums(gad_mean[,1]),ylim=c(3,12),xlim=c(0,10),
    ylab="Mean GAD-7 score",
    xlab="Years after starting treatment",
    #ylab="",
    #xlab="",
     type="l",lwd=lnwds[1],lty=lntys[1])
for(i in 1:n_treat){
  lines(x_yrs,nums(gad_mean[,i]),lwd=lnwds[i],lty=lntys[i])
}
legend("topleft",t_names,lwd=lnwds,lty=lntys,bty = "n", horiz=TRUE, cex=0.9)
dev.off()
 
m_gad_chg<-gad_mean[,1]-gad_mean[,]
 
if(b_assum==1){
  jpeg(f_gph_name("mean_te_over_time_small"), width = 280, height = 262)
  par(mar=c(4,0.2,0.2,0.2))
} else {
  jpeg(f_gph_name("mean_te_over_time_small"), width = 280, height = 218)
  par(mar=c(0.2,0.2,0.2,0.2))
}
  
plot(x_yrs,nums(m_gad_chg[,1]),
     ylim=c(0,7),xlim=c(0,10),#full sample
     ylab="Mean reduction in GAD-7 score",
     xlab="Years after starting treatment",
     type="l",lwd=lnwds[1],lty=lntys[1])
for(i in 2:n_treat){
  lines(x_yrs,nums(m_gad_chg[,i]),lwd=lnwds[i],lty=lntys[i])
}
 if(t_assum==1){
   legend(0,8,t_names,lwd=c(1,1,2,1,2,2,3),lty=c(1,3,3,5,5,1,1),bty = "n",horiz=TRUE, cex=0.9)
 }
dev.off()

  #} #remove hash if running all scenarios at once


######## apply QALYs and costs ########

qaly_count<-qaly_count_a<-cost_count_a<-array(NA,c(n_samps,n_cycle,n_treat),list(1:n_samps, 1:n_cycle, t_names))

for (i in 1:n_cycle){
  for (j in 1:n_treat){
    
    qaly_count[,i,j]<-ifelse(as.numeric(counter[,i,j])<=4,qaly_mat[,1,j],
                             ifelse(as.numeric(counter[,i,j])<=9,qaly_mat[,2,j],
                                    ifelse(as.numeric(counter[,i,j])<=14,qaly_mat[,3,j],qaly_mat[,4,j])))
    cost_count_a[,i,j]<-ifelse(as.numeric(counter[,i,j])<=4,cost_mat[,1,j],
                             ifelse(as.numeric(counter[,i,j])<=9,cost_mat[,2,j],
                                    ifelse(as.numeric(counter[,i,j])<=14,cost_mat[,3,j],cost_mat[,4,j])))
  }}

for (i in 1:n_cycle){
  
  qaly_count[,i,][counter[,i,]%in%state1]<-qaly_mat[,1,][counter[,i,]%in%state1]
  qaly_count[,i,][counter[,i,]%in%state2]<-qaly_mat[,2,][counter[,i,]%in%state2]
  qaly_count[,i,][counter[,i,]%in%state3]<-qaly_mat[,3,][counter[,i,]%in%state3]
  qaly_count[,i,][counter[,i,]%in%state4]<-qaly_mat[,4,][counter[,i,]%in%state4]
  
  cost_count_a[,i,][counter[,i,]%in%state1]<-cost_mat[,1,][counter[,i,]%in%state1]
  cost_count_a[,i,][counter[,i,]%in%state2]<-cost_mat[,2,][counter[,i,]%in%state2]
  cost_count_a[,i,][counter[,i,]%in%state3]<-cost_mat[,3,][counter[,i,]%in%state3]
  cost_count_a[,i,][counter[,i,]%in%state4]<-cost_mat[,4,][counter[,i,]%in%state4]  
}

#adjust for age related qaly decrement and cost change
qaly_count_a[,1:4,] <- qaly_count[,1:4,]

for (i in 5:n_cycle){
  
  cycl_age<-floor(age+(i-1)/4)
  
  if(cycl_age>floor(age+(i-2)/4)){
    
    for (k in 1:n_samps){qaly_count_a[k,i,]<-qaly_count_a[k,i-1,]+qaly_rda[k,cycl_age,]/4}
    
  } else {
    
    for (k in 1:n_samps){qaly_count_a[k,i,]<-qaly_count_a[k,i-1,]}
    
  }
  
}

rm(qaly_count,cycl_age)

######## discount #######

qaly_count_d<-cost_count_d<-array(NA,c(n_samps,n_cycle,n_treat),list(1:n_samps, 1:n_cycle, t_names))

 for (i in 1:n_cycle){
   qaly_count_d[,i,]<-qaly_count_a[,i,]/(dr^((i-1)/4))
   cost_count_d[,i,]<-cost_count_a[,i,]/(dr^((i-2)/4))
 }

rm(qaly_count_a,cost_count_a)

######## apply mortality ########

mort_sta<-surv_sta<-array(NA,c(n_samps,n_cycle,n_treat),list(1:n_samps, 1:n_cycle, colnames(gad7_scr)))
surv_sta[,1,]<-1

for (b in 1:n_cycle){
  for (c in 1:n_treat){
    mort_sta[,b,c]<-ifelse(counter[,b,c] <=4, mort_cyc[,b],
                           ifelse(counter[,b,c] <=9, mort_exs[,2,c]*mort_cyc[,b],
                                  ifelse(counter[,b,c] <=14, mort_exs[,3,c]*mort_cyc[,b],mort_exs[,4,c]*mort_cyc[,b])))
  }
  if(b>1){
    surv_sta[,b,]<-surv_sta[,b-1,]*(1-mort_sta[,b,])
  }
}

rm(mort_sta)

qaly_tot <- qaly_count_d * surv_sta
cost_tot <- cost_count_d * surv_sta

cost_int2 <- array(0, c(n_samps, n_cycle, n_treat))
cost_int2[,1,] <- cost_int
for (i in 1:n_samps){
  cost_int2[i,1:20,7]<-cost_int2[i,1:20,7]+(cost_med[ceiling((1:20)/4)]/4)*surv_sta[i,1:20,7]
}
cost_int3<-cost_int2
rm(cost_int2)


####### cummulative outcomes over time #######

 cum_lyrs_over_time <- apply(surv_sta, c(1,3), cumsum) / 4
 cum_qaly_over_time <- apply(qaly_tot, c(1,3), cumsum)
 cum_cost_over_time <- apply(cost_tot, c(1,3), cumsum)
 tot_cost_over_time <- cost_tot + cost_int3
 cum_ctot_over_time <- apply(tot_cost_over_time, c(1,3), cumsum)
 
jpeg(f_gph_name("cum_lys_over_time"), width = 450, height = 260)
par(mar=c(2,4,0.2,0.2))
plot(x_yrs,apply(cum_lyrs_over_time[,,1],1,mean),xlim=c(0,50),ylim=c(0,40),type="l",lwd=2,
     xlab="Years after treatment",ylab="Cumulative mean LY gain")
for(k in 2:n_treat){
  lines(x_yrs,apply(cum_lyrs_over_time[,,k],1,mean),lty=lntys[k],lwd=lnwds[k])
}
legend(-0.2,42,t_names,lwd=lnwds,lty=lntys,cex=1.1,bty='n')
dev.off()

jpeg(f_gph_name("cum_qalys_over_time"), width = 450, height = 260)  # bland Altman plot
par(mar=c(2,4,0.2,0.2))
plot(x_yrs,apply(cum_qaly_over_time[,,1],1,mean),xlim=c(0,50),ylim=c(0,15),type="l",lwd=2,
     xlab="Years after starting treatment",ylab="Cumulative mean QALY gain")
for(k in 2:n_treat){
  lines(x_yrs,apply(cum_qaly_over_time[,,k],1,mean),lty=lntys[k],lwd=lnwds[k])
}
legend(-0.2,15.1,t_names,lwd=lnwds,lty=lntys,cex=0.75,bty='n')
dev.off()

jpeg(f_gph_name("cum_cost_over_time"), width = 450, height = 280)  # bland Altman plot
par(mar=c(4,4,0.2,0.2))
plot(x_yrs,apply(cum_cost_over_time[,,1],1,mean),xlim=c(0,50),ylim=c(0,17000),type="l",lwd=2,
     xlab="Years after starting treatment", ylab="Cumulative mean healthcare cost")
for(k in 2:n_treat){
  lines(x_yrs,apply(cum_cost_over_time[,,k],1,mean),lty=lntys[k],lwd=lnwds[k])
}
legend(-0.2,26000,t_names,lwd=lnwds,lty=lntys,cex=0.75,bty='n')
dev.off()

jpeg(f_gph_name("cum_ctot_over_time"), width = 450, height = 280)  # bland Altman plot
par(mar=c(4,4,0.2,0.2))
plot(x_yrs,apply(cum_ctot_over_time[,,1],1,mean),xlim=c(0,50),ylim=c(0,17000),type="l",lwd=2,
     xlab="Years after starting treatment", ylab="Cumulative mean total cost")
for(k in 2:n_treat){
  lines(x_yrs,apply(cum_ctot_over_time[,,k],1,mean),lty=lntys[k],lwd=lnwds[k])
}
legend(-0.2,26000,t_names,lwd=lnwds,lty=lntys,cex=0.75,bty='n')
dev.off()

 nmb_time<-matrix(NA,n_cycle,n_treat)
 
 for (i in 1:n_cycle){
   for (t in 1:n_treat){
       
      nmb_time[i,t]<-mean(15000*cum_qaly_over_time[i,,t]-tot_cost_over_time[,i,t])
       
   }
 }
 
 rm(cum_lyrs_over_time,cum_qaly_over_time,cum_cost_over_time,cum_ctot_over_time)

 jpeg(f_gph_name("nmb_5years"), width = 600, height = 370)
 par(mar=c(4,4,0.2,0.2))
 plot(x_yrs[1:20],nmb_time[1:20,1]/1000, type="l",
      lwd=lnwds[1], lty=lntys[1],
      ylim=c(0,50),
      xlab="Years after treatment", ylab="Net benefit (£1000s)")
 for(t in 2:n_treat){
   lines(x_yrs,nmb_time[,t]/1000, lwd=lnwds[t], lty=lntys[t])
 }
 dev.off()
 
 jpeg(f_gph_name("nmb_50years"), width = 600, height = 370)
 par(mar=c(4,4,0.2,0.2))
 plot(x_yrs,nmb_time[,1]/1000, type="l",lwd=lnwds[1], lty=lntys[1],
      ylim=c(0,200),xlab="Years after treatment", ylab="Net benefit (£1000s)")
 for(t in 2:n_treat){
   lines(x_yrs,nmb_time[,t]/1000, lwd=lnwds[t], lty=lntys[t])
 }
 dev.off()
 
 
####### add costs and qalys #######

qaly_sum_all<-apply(qaly_tot,c(1,3),sum)
lyrs_sum_all<-apply(surv_sta,c(1,3),sum)/4
cost_sum_hcr<-apply(cost_tot,c(1,3),sum)
cost_sum_int<-apply(cost_int3,c(1,3),sum)
cost_sum_tot<-cost_sum_hcr+cost_sum_int
lambda<-seq(0,30000,1000)

write.csv(qaly_sum_all, file =f_tab_name("qalys"), row.names=FALSE)
write.csv(cost_sum_tot, file =f_tab_name("costs"), row.names=FALSE)

####### import results (use if analysing pre-saved results) ########

qaly_sum_all <- as.matrix(read.csv(file=paste0(tabl_folder,"costs_b2t1.csv"), header=TRUE))
cost_sum_tot <- as.matrix(read.csv(file=paste0(tabl_folder,"qalys_b2t1.csv"), header=TRUE))

n_samps<-nrow(cost_sum_tot)
lambda<-seq(0,30000,1000)

######## tabulate results #########

#if excluding digital controls
# n_treat<-5
# t_names<-t_names[-c(2:3)]
# cost_sum_tot<-cost_sum_tot[,-c(2:3)]
# qaly_sum_all<-qaly_sum_all[,-c(2:3)]

nmb_mod1<-array(NA,c(n_samps,length(lambda),n_treat))
pce_mod1<-matrix(NA,length(lambda), n_treat)

for (l in 1:length(lambda)){
  nmb_mod1[,l,]<-as.matrix(qaly_sum_all*lambda[l]-cost_sum_tot)
  }

for (l in 1:length(lambda)){
  maxnmb<-apply(nmb_mod1[,l,],1,max)
  for (t in 1:n_treat){
    pce_mod1[l,t]<-mean(nmb_mod1[,l,t]==maxnmb)
  }
}

temp<-matrix(NA,nrow(pce_mod1),2)
temp[,1]<-lambda

for (i in 1:nrow(pce_mod1)){
  temp[i,2]<-which.max(pce_mod1[i,])
}

jpeg(f_gph_name("ceac_colours"), width = 600, height = 370)
rnd_col<-c("grey","red","black","blue","red","black","blue")
#rnd_col<-c("grey","blue","red","black","blue")
par(mar=c(4,4,0.2,0.2))
plot(lambda, pce_mod1[,1], ylim=c(0,0.5), xlim=c(lambda[1], lambda[length(lambda)]), type="l",
     lwd=lnwds[1], lty=lntys[1], col=rnd_col[1], xlab="Opportunity cost (£/QALY)", ylab="Probability cost-effective")
for (t in 2:n_treat){lines(lambda, pce_mod1[,t], lwd=lnwds[t], lty=lntys[t], col=rnd_col[t])}
legend(-2,0.61,t_names,lwd=lnwds,lty=lntys,col=rnd_col,bty = "n",horiz=TRUE, cex=0.9)
dev.off()

results_mod1<-matrix(NA,11,n_treat)

results_mod1[ 1,]<-paste(round(apply(qaly_sum_all,2,mean),digits=2)," (",
                         round(apply(qaly_sum_all,2,quantile, probs=0.025),digits=2)," to ",
                         round(apply(qaly_sum_all,2,quantile, probs=0.975),digits=2),")",sep="")
results_mod1[ 2,]<-paste(round(apply(lyrs_sum_all,2,mean),digits=2)," (",
                         round(apply(lyrs_sum_all,2,quantile, probs=0.025),digits=2)," to ",
                         round(apply(lyrs_sum_all,2,quantile, probs=0.975),digits=2),")",sep="")
results_mod1[ 3,]<-paste(format(round(apply(cost_sum_hcr,2,mean),digits=0),big.mark=",",scientific=FALSE)," (",
                         format(round(apply(cost_sum_hcr,2,quantile, probs=0.025),digits=0),big.mark=",",scientific=FALSE)," to ",
                         format(round(apply(cost_sum_hcr,2,quantile, probs=0.975),digits=0),big.mark=",",scientific=FALSE),")",sep="")
results_mod1[ 4,]<-paste(format(round(apply(cost_sum_int,2,mean),digits=0),big.mark=",",scientific=FALSE))
results_mod1[4,7]<-paste(format(round(apply(cost_sum_int,2,mean)[7],digits=0),big.mark=",",scientific=FALSE)," (",
                         format(round(apply(cost_sum_int,2,quantile, probs=0.025)[7],digits=0),big.mark=",",scientific=FALSE)," to ",
                         format(round(apply(cost_sum_int,2,quantile, probs=0.975)[7],digits=0),big.mark=",",scientific=FALSE),")",sep="")
results_mod1[ 5,]<-paste(format(round(apply(cost_sum_tot,2,mean),digits=0),big.mark=",",scientific=FALSE)," (",
                         format(round(apply(cost_sum_tot,2,quantile, probs=0.025),digits=0),big.mark=",",scientific=FALSE)," to ",
                         format(round(apply(cost_sum_tot,2,quantile, probs=0.975),digits=0),big.mark=",",scientific=FALSE),")",sep="")
results_mod1[ 6,]<-paste(format(round(apply(nmb_mod1[, 1,],2,mean),digits=0),big.mark=",",scientific=FALSE)," (",
                         format(round(apply(nmb_mod1[, 1,],2,quantile, probs=0.025),digits=0),big.mark=",",scientific=FALSE)," to ",
                         format(round(apply(nmb_mod1[, 1,],2,quantile, probs=0.975),digits=0),big.mark=",",scientific=FALSE),")",sep="")
results_mod1[ 8,]<-paste(format(round(apply(nmb_mod1[,16,],2,mean),digits=0),big.mark=",",scientific=FALSE)," (",
                         format(round(apply(nmb_mod1[,16,],2,quantile, probs=0.025),digits=0),big.mark=",",scientific=FALSE)," to ",
                         format(round(apply(nmb_mod1[,16,],2,quantile, probs=0.975),digits=0),big.mark=",",scientific=FALSE),")",sep="")
results_mod1[10,]<-paste(format(round(apply(nmb_mod1[,31,],2,mean),digits=0),big.mark=",",scientific=FALSE)," (",
                         format(round(apply(nmb_mod1[,31,],2,quantile, probs=0.025),digits=0),big.mark=",",scientific=FALSE)," to ",
                         format(round(apply(nmb_mod1[,31,],2,quantile, probs=0.975),digits=0),big.mark=",",scientific=FALSE),")",sep="")
results_mod1[ 7,]<-round(pce_mod1[ 1,],digits=3)
results_mod1[ 9,]<-round(pce_mod1[16,],digits=3)
results_mod1[11,]<-round(pce_mod1[31,],digits=3)

rownames(results_mod1)<-c("QALYs", "LYs", "Cost of HC", "Intervention cost", "Total cost",
                          "NB, threshold=0","PCE, threshold=0", "NB, threshold=15K","PCE, threshold=15K",
                          "NB, threshold=30K","PCE, threshold=30K")
colnames(results_mod1)<-t_names

results<-t(results_mod1)

write.csv(results, file =f_tab_name("t_results"), row.names=FALSE)

#}
#}

############ Present results ############

####### NMB at £15000 over time ########

#CEAC
par(mar=c(4,4,0.2,0.2))
plot (lambda,pce_mod1[,1],ylim=c(0,0.6),lwd=2,"l",
      ylab="Probability cost-effective",
      xlab="Opportunity cost")
lines(lambda,pce_mod1[,2],lwd=1,lty=1)
lines(lambda,pce_mod1[,3],lwd=1,lty=2)
lines(lambda,pce_mod1[,4],lwd=1,lty=3)
lines(lambda,pce_mod1[,5],lwd=1,lty=4)
lines(lambda,pce_mod1[,6],lwd=2,lty=2)
lines(lambda,pce_mod1[,7],lwd=2,lty=3)
#legend(0,0.65,t_names,lwd=lnwds,lty=lntys)

temp<-rep(NA,n_samps)
for (i in 1:n_samps){
  temp[i]<-length(unique(nmb_mod1[i,1,]))
}

######### EVPI #########

max_nmb<-nmb_max<-matrix(NA,n_samps, length(lambda))
evi_ind<-rep(NA,length(lambda))

pop<-c(250000,500000)
tim<-c(1,5,10)
evi_mult<-matrix(NA,length(pop),length(tim))
for(i in 1:length(pop)){
  for (j in 1:length(tim)){
    evi_mult[i,j]<-cumsum(pop[i]/(dr^(0:(tim[j]-1))))[tim[j]]
  }
}
evi_mult<-unlist(melt(evi_mult)[3])

rm(pop, tim)

evi_pop<-matrix(NA,length(lambda), length(evi_mult))

for (l in 1:length(lambda)){
  max_nmb [,l] <-apply(nmb_mod1[,l,],1,max)
  nmb_max [,l] <-nmb_mod1[,l,which.max(apply(nmb_mod1[,l,],2,mean))]
  evi_ind [l]  <- mean(max_nmb[,l]-nmb_max[,l])
  evi_pop [l,] <- evi_ind[l]*evi_mult
}

jpeg(f_gph_name("popEVPI"), width = 600, height = 370)  # bland Altman plot
par(mar=c(4,4,0.2,0.2))
plot(lambda,evi_pop[,3]/1000000000, 
     ylim=c(5,25),
     xlab="Opportunity cost (£/QALY)",
     ylab="Population EVPI (£ billion)",
     type="l")
dev.off()

round(min(evi_pop[,3])/1000000000,digits=1)
lambda[which.min(evi_pop[,3])]

######### 1 time horizon and population ###########

max_nmb<-nmb_max<-matrix(NA,n_samps, length(lambda))
evi_ind<-rep(NA,length(lambda))

pop<-250000
tim<-5
evi_mult<-cumsum(pop/(dr^(0:(tim-1))))[tim]

evi_pop<-rep(NA,length(lambda))

for (l in 1:length(lambda)){
  max_nmb [,l] <-apply(nmb_mod1[,l,],1,max)
  nmb_max [,l] <-nmb_mod1[,l,which.max(apply(nmb_mod1[,l,],2,mean))]
  evi_ind [l]  <- mean(max_nmb[,l]-nmb_max[,l])
  evi_pop [l] <- evi_ind[l]*evi_mult
}

#jpeg(f_gph_name("popEVPI_no_control"), width = 600, height = 370)  # bland Altman plot
par(mar=c(4,4,0.2,0.2))
plot(lambda,evi_pop/1000000000, 
     ylim=c(5,25),
     xlab="Opportunity cost (£/QALY)",
     ylab="Population EVPI (£ billion)",
     type="l")
#dev.off()

round(min(evi_pop)/1000000000,digits=1)
lambda[which.min(evi_pop)]


