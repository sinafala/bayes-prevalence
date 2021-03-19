######################################
#
# OCAP - Ohio COVID-19 Prevalence
# Annotated Code
# Antibody Tests - Past Infection
# 
# https://arxiv.org/abs/2011.09033
#
# 11/23/20
#
# Notes:
# - Key parts of code to illustrate
#   data processing, model, & summary
# - Omitted some trivial loading of
#   data since it cannot be shared
# - May not run due to omitted steps
#
######################################

#####################################
#Preamble
#####################################
#load packages
library(nimble)
library(coda)
library(mcmcplots)
library(plyr)
library(reshape)

#####################################
#Read in data
#####################################
#Omitted for brevity since data are restricted

#####################################
#Process data
#####################################
#calculate strata
#age groups
ocap.data$age.group=cut(ocap.data$age,c(18,44,64,max(ocap.data$age)),include.lowest=TRUE)

#gender 1=M 2=F
ocap.data$sex.group=NA
ocap.data$sex.group[ocap.data$s15==1]='Male'
ocap.data$sex.group[ocap.data$s15==2]='Female'
ocap.data$sex.group=factor(ocap.data$sex.group,levels=c('Male','Female'))

#LC testing data 1=postiive 2=negative 3=equivocal
ocap.data$LC_IgM=NA
ocap.data$LC_IgM[ocap.data$igm_qual==1]='Positive'
ocap.data$LC_IgM[ocap.data$igm_qual==2]='Negative'

ocap.data$LC_IgG=NA
ocap.data$LC_IgG[ocap.data$igg_qual==1]='Positive'
ocap.data$LC_IgG[ocap.data$igg_qual==2]='Negative'

#OSU testing data 1=detected 2=not detected 3=other
ocap.data$OSU_IgG=NA
ocap.data$OSU_IgG[ocap.data$igg_osu==1]='Positive'
ocap.data$OSU_IgG[ocap.data$igg_osu==2]='Negative'

#####################################
#Data for analysis
#####################################
#select variables from survey needed for analysis
ocap.vars=c('PIN','GEOID','age.group','sex.group','LC_IgM','LC_IgG','OSU_IgG','test_pattern')
#select variables required to be complete
ocap.cc=c('PIN','GEOID','age.group','sex.group','test_pattern')

#create indicator of which tests are observed for each subject
#indicator of missing test
ocap.data$test_pattern=NA
#all 3 observed
ocap.data$test_pattern[!(is.na(ocap.data$LC_IgM) & is.na(ocap.data$LC_IgG) & is.na(ocap.data$OSU_IgG))]=1

#patterns with a subset of tests observed
#LC_IgM LC_IgG
ocap.data$test_pattern[(!is.na(ocap.data$LC_IgM) & !is.na(ocap.data$LC_IgG) & is.na(ocap.data$OSU_IgG))]=2
#LC_IgM OSU_IgG
ocap.data$test_pattern[(!is.na(ocap.data$LC_IgM) & is.na(ocap.data$LC_IgG) & !is.na(ocap.data$OSU_IgG))]=3
#LC_IgG OSU_IgG
ocap.data$test_pattern[(is.na(ocap.data$LC_IgM) & !is.na(ocap.data$LC_IgG) & !is.na(ocap.data$OSU_IgG))]=4
#LC_IgM
ocap.data$test_pattern[(!is.na(ocap.data$LC_IgM) & is.na(ocap.data$LC_IgG) & is.na(ocap.data$OSU_IgG))]=5
#LC_IgG
ocap.data$test_pattern[(is.na(ocap.data$LC_IgM) & !is.na(ocap.data$LC_IgG) & is.na(ocap.data$OSU_IgG))]=6
#OSU_IgG
ocap.data$test_pattern[(is.na(ocap.data$LC_IgM) & is.na(ocap.data$LC_IgG) & !is.na(ocap.data$OSU_IgG))]=7


#complete survey data limited to variables needed
ocap.analysis=ocap.data[complete.cases(ocap.data[,ocap.cc]),ocap.vars]

#merge in region ID
ocap.analysis=join(ocap.analysis,region.data[,c('GEOID','REGION_ML')],by='GEOID',type='left')

#add population - standardized log population
ocap.analysis=join(ocap.analysis,pop.data[,c('GEOID','st.log.pop')],by='GEOID',type='left')

#sort by pattern for convenient indexing in MCMC
ocap.analysis=ocap.analysis[order(ocap.analysis$test_pattern),]

#constants for MCMC
#number of tests
num.tests=3
#number of subjects
subjects=length(unique(ocap.analysis$PIN))
#number of regions
regions=length(unique(ocap.analysis$REGION_ML))
#number of tracts
tracts=length(unique(ocap.analysis$GEOID))

#create indices defining groups by which tests observed
nt=rep(0,7)
nt[1]=sum(ocap.analysis$test_pattern==1)
nt[2]=sum(ocap.analysis$test_pattern==2)+nt[1]
nt[3]=sum(ocap.analysis$test_pattern==3)+nt[2]
nt[4]=sum(ocap.analysis$test_pattern==4)+nt[3]
nt[5]=sum(ocap.analysis$test_pattern==5)+nt[4]
nt[6]=sum(ocap.analysis$test_pattern==6)+nt[5]
nt[7]=sum(ocap.analysis$test_pattern==7)+nt[6]

#map test results to result patterns
#columns 8,12,16,20 all negative results
test.results=array(0,dim=c(subjects,23))

#1,1,1
test.results[,1]=(ocap.analysis$LC_IgG=='Positive' & ocap.analysis$OSU_IgG=='Positive' & ocap.analysis$LC_IgM=='Positive')

#1,1,0
test.results[,2]=(ocap.analysis$LC_IgG=='Positive' & ocap.analysis$OSU_IgG=='Positive' & ocap.analysis$LC_IgM=='Negative')

#1,0,1
test.results[,3]=(ocap.analysis$LC_IgG=='Positive' & ocap.analysis$OSU_IgG=='Negative' & ocap.analysis$LC_IgM=='Positive')

#0,1,1
test.results[,4]=(ocap.analysis$LC_IgG=='Negative' & ocap.analysis$OSU_IgG=='Positive' & ocap.analysis$LC_IgM=='Positive')

#1,0,0
test.results[,5]=(ocap.analysis$LC_IgG=='Positive' & ocap.analysis$OSU_IgG=='Negative' & ocap.analysis$LC_IgM=='Negative')

#0,1,0
test.results[,6]=(ocap.analysis$LC_IgG=='Negative' & ocap.analysis$OSU_IgG=='Positive' & ocap.analysis$LC_IgM=='Negative')

#0,0,1
test.results[,7]=(ocap.analysis$LC_IgG=='Negative' & ocap.analysis$OSU_IgG=='Negative' & ocap.analysis$LC_IgM=='Positive')

#0,0,0
test.results[,8]=(ocap.analysis$LC_IgG=='Negative' & ocap.analysis$OSU_IgG=='Negative' & ocap.analysis$LC_IgM=='Negative')

#patterns with at least 1 test result not observed
#1,NA,1
test.results[,9]=(ocap.analysis$LC_IgG=='Positive' & is.na(ocap.analysis$OSU_IgG) & ocap.analysis$LC_IgM=='Positive')

#1,NA,0
test.results[,10]=(ocap.analysis$LC_IgG=='Positive' & is.na(ocap.analysis$OSU_IgG) & ocap.analysis$LC_IgM=='Negative')

#0,NA,1
test.results[,11]=(ocap.analysis$LC_IgG=='Negative' & is.na(ocap.analysis$OSU_IgG) & ocap.analysis$LC_IgM=='Positive')

#0,NA,0
test.results[,12]=(ocap.analysis$LC_IgG=='Negative' & is.na(ocap.analysis$OSU_IgG) & ocap.analysis$LC_IgM=='Negative')

#NA,1,1
test.results[,13]=(is.na(ocap.analysis$LC_IgG) & ocap.analysis$OSU_IgG=='Positive' & ocap.analysis$LC_IgM=='Positive')

#NA,1,0
test.results[,14]=(is.na(ocap.analysis$LC_IgG) & ocap.analysis$OSU_IgG=='Positive' & ocap.analysis$LC_IgM=='Negative')

#NA,0,1
test.results[,15]=(is.na(ocap.analysis$LC_IgG) & ocap.analysis$OSU_IgG=='Negative' & ocap.analysis$LC_IgM=='Positive')

#NA,0,0
test.results[,16]=(is.na(ocap.analysis$LC_IgG) & ocap.analysis$OSU_IgG=='Negative' & ocap.analysis$LC_IgM=='Negative')

#1,1,NA
test.results[,17]=(ocap.analysis$LC_IgG=='Positive' & ocap.analysis$OSU_IgG=='Positive' & is.na(ocap.analysis$LC_IgM))

#1,0,NA
test.results[,18]=(ocap.analysis$LC_IgG=='Positive' & ocap.analysis$OSU_IgG=='Negative' & is.na(ocap.analysis$LC_IgM))

#0,1,NA
test.results[,19]=(ocap.analysis$LC_IgG=='Negative' & ocap.analysis$OSU_IgG=='Positive' & is.na(ocap.analysis$LC_IgM))

#0,0,NA
test.results[,20]=(ocap.analysis$LC_IgG=='Negative' & ocap.analysis$OSU_IgG=='Negative' & is.na(ocap.analysis$LC_IgM))

#NA,NA,1
test.results[,21]=(is.na(ocap.analysis$LC_IgG) & is.na(ocap.analysis$OSU_IgG) & ocap.analysis$LC_IgM=='Positive')

#1,NA,NA
test.results[,22]=(ocap.analysis$LC_IgG=='Positive' & is.na(ocap.analysis$OSU_IgG) & is.na(ocap.analysis$LC_IgM))

#NA,1,NA
test.results[,23]=(is.na(ocap.analysis$LC_IgG) & ocap.analysis$OSU_IgG=='Positive' & is.na(ocap.analysis$LC_IgM))


#fixed effects matrix
X=model.matrix(~1+age.group+sex.group+st.log.pop,data=ocap.analysis,contrasts=list(age.group='contr.sum',sex.group='contr.sum'))


#####################################
#Model
#####################################
#nimble model code
model_code=nimbleCode({
  
  #data model
  #all 3 tests observed
  for(i in 1:nt[1]){
    
    #test result pattern
    t[i,1:8] ~ dmulti(p[i,1:8],1)
    
    #IgG results
    p12[i,1]<- d[i]*(s1*s2+r12.pos)+(1-d[i])*((1-c1)*(1-c2)+r12.neg) #++
    p12[i,2]<- d[i]*(s1*(1-s2)-r12.pos)+(1-d[i])*((1-c1)*c2-r12.neg) #+-
    p12[i,3]<- d[i]*((1-s1)*s2-r12.pos)+(1-d[i])*(c1*(1-c2)-r12.neg) #-+
    p12[i,4]<- d[i]*((1-s1)*(1-s2)+r12.pos)+(1-d[i])*(c1*c2+r12.neg) #--
    
    #IgM results
    p3[i,1]<- d[i]*(s3)+(1-d[i])*((1-c3)) #+
    p3[i,2]<- d[i]*((1-s3))+(1-d[i])*(c3) #-
    
    #combined probabilities
    p[i,1]<- p12[i,1]*p3[i,1] #+++
    p[i,2]<- p12[i,1]*p3[i,2] #++-
    p[i,3]<- p12[i,2]*p3[i,1] #+-+
    p[i,4]<- p12[i,3]*p3[i,1] #-++
    p[i,5]<- p12[i,2]*p3[i,2] #+--
    p[i,6]<- p12[i,3]*p3[i,2] #-+-
    p[i,7]<- p12[i,4]*p3[i,1] #--+
    p[i,8]<- p12[i,4]*p3[i,2] #---
    
  }
  
  #models for patterns with at least 1 test unobserved
  
  for(i in (nt[1]+1):nt[2]){
    
    #test result pattern
    t[i,9:12] ~ dmulti(p[i,9:12],1)
    
    #IgG results
    p1[i,1]<- d[i]*(s1)+(1-d[i])*((1-c1)) #+
    p1[i,2]<- d[i]*((1-s1))+(1-d[i])*(c1) #-
    
    #IgM results
    p3[i,1]<- d[i]*(s3)+(1-d[i])*((1-c3)) #+
    p3[i,2]<- d[i]*((1-s3))+(1-d[i])*(c3) #-
    
    #combined probabilities
    p[i,9]<- p1[i,1]*p3[i,1] #++
    p[i,10]<- p1[i,1]*p3[i,2] #+-
    p[i,11]<- p1[i,2]*p3[i,1] #-+
    p[i,12]<- p1[i,2]*p3[i,2] #--
    
  }
  
  #Note: no subjects had this pattern
  #for(i in (nt[2]+1):nt[3]){
  
    #test result pattern
    #t[i,13:16] ~ dmulti(p[i,13:16],1)
  
    #IgG results
    #p2[i,1]<- d[i]*(s2)+(1-d[i])*((1-c2)) #+
    #p2[i,2]<- d[i]*((1-s2))+(1-d[i])*(c2) #-
  
    #IgM results
    #p3[i,1]<- d[i]*(s3)+(1-d[i])*((1-c3)) #+
    #p3[i,2]<- d[i]*((1-s3))+(1-d[i])*(c3) #-
  
    #combined probabilities
    #p[i,13]<- p2[i,1]*p3[i,1] #++
    #p[i,14]<- p2[i,1]*p3[i,2] #+-
    #p[i,15]<- p2[i,2]*p3[i,1] #-+
    #p[i,16]<- p2[i,2]*p3[i,2] #--
  
  #}
  
  for(i in (nt[3]+1):nt[4]){
    
    #test result pattern
    t[i,17:20] ~ dmulti(p[i,17:20],1)
    
    #IgG results
    p12[i,1]<- d[i]*(s1*s2+r12.pos)+(1-d[i])*((1-c1)*(1-c2)+r12.neg) #++
    p12[i,2]<- d[i]*(s1*(1-s2)-r12.pos)+(1-d[i])*((1-c1)*c2-r12.neg) #+-
    p12[i,3]<- d[i]*((1-s1)*s2-r12.pos)+(1-d[i])*(c1*(1-c2)-r12.neg) #-+
    p12[i,4]<- d[i]*((1-s1)*(1-s2)+r12.pos)+(1-d[i])*(c1*c2+r12.neg) #--
    
    #combined probabilities
    p[i,17]<- p12[i,1] #++
    p[i,18]<- p12[i,2] #+-
    p[i,19]<- p12[i,3] #-+
    p[i,20]<- p12[i,4] #--
    
  }
  
  for(i in (nt[4]+1):nt[5]){
    
    #test result pattern
    t[i,21] ~ dbern(p[i,21])
    
    #IgM results
    p3[i,1]<- d[i]*(s3)+(1-d[i])*((1-c3)) #+
    
    #combined probabilities
    p[i,21]<- p3[i,1] #+
    
  }
  
  
  for(i in (nt[5]+1):nt[6]){
    
    #test result pattern
    t[i,22] ~ dbern(p[i,22])
    
    #IgG results
    p1[i,1]<- d[i]*(s1)+(1-d[i])*((1-c1)) #+
    
    #combined probabilities
    p[i,22]<- p1[i,1] #+
    
  }
  
  #Note: no subjects had this pattern
  #for(i in (nt[6]+1):nt[7]){
  
    #test result pattern
    #t[i,23] ~ dbern(p[i,23])
  
    #IgG results
    #p2[i,1]<- d[i]*(s2)+(1-d[i])*((1-c2)) #+
  
    #combined probabilities
    #p[i,23]<- p2[i,1] #+
  
  #}
  
  #process model
  for(i in 1:n){   
    
    #prevalence model
    d[i] ~ dbern(pi[i])
    logit(pi[i])<- inprod(X[i,1:xp],beta[1:xp])+alpha[region[i]]+b[tract[i]]
    
  }
  
  
  #random region effects
  for(i in 1:r){
    alpha[i] ~ dnorm(mu,sd=sigma)
  }
  
  #tract random effects
  for(i in 1:ct){
    b[i] ~ dnorm(0,sd=tau)
  }
  
  #prior model
  #sensitivity and specificity of LC IgG
  s1 ~ dbeta(alpha.s1,beta.s1)
  c1 ~ dbeta(alpha.c1,beta.c1)
  #sensitivity and specificity of OSU IgG
  s2 ~ dbeta(alpha.s2,beta.s2)
  c2 ~ dbeta(alpha.c2,beta.c2)
  #sensitivity and specificity of LC IgM
  s3 ~ dbeta(alpha.s3,beta.s3)
  c3 ~ dbeta(alpha.c3,beta.c3)
  
  #covariance between results of test 1 and 2 given the disease status is positive
  r12.pos ~ dunif(0, (min(s1,s2)-(s1*s2)) )
  #covariance between results of test 1 and 2 given the disease status is negative
  r12.neg ~ dunif(0, (min(c1,c2)-(c1*c2)) )
  
  #fixed effects
  for(i in 1:xp){
    beta[i] ~ dnorm(0,sd=3)
  }
  
  #overall intercept
  mu ~ dnorm(log(.03/.97),sd=1)
  
  #region effect SD
  sigma ~ dunif(0,5)
  
  #tract effect SD
  tau ~ dunif(0,5)
  
})

#####################################
#NIMBLE - MCMC
#####################################
#constants
#number of fixed effect parameters - removing intercept
xp=dim(X[,-1])[2]
#initial values for probabilities
p.init=array(NA,dim=c(subjects,23))
p.init[,1:8]=1/8
p.init[,9:20]=1/4
p.init[,21:23]=1/2

#model inputs
constants=list(n=subjects,xp=xp,r=regions,ct=tracts,nt=nt,
               region=as.numeric(ocap.analysis$REGION_ML),tract=as.numeric(as.factor(ocap.analysis$GEOID)),
               alpha.s1=109,beta.s1=13,alpha.s2=96,beta.s2=39,alpha.s3=9,beta.s3=11,
               alpha.c1=1066,beta.c1=4,alpha.c2=1074,beta.c2=16,alpha.c3=54,beta.c3=0.1)

#test results and design matrix without intercept
model.data=list(t=test.results,X=X[,-1])

#initial values
inits=list(p=p.init,d=rep(0,subjects),alpha=rep(0,regions),sigma=1,beta=rep(0,xp),
           b=rep(0,tracts),tau=1,mu=-3)

#load limited data to run code
load('ocap_antibody_limited_data.Rda')

#build the model
model <- nimbleModel(model_code, constants,model.data,inits)
compiled_model <- compileNimble(model,resetFunctions = TRUE)

#set up samplers
mcmc_conf <- configureMCMC(model,monitors=c('alpha','sigma','beta','s1','c1','s2','c2','s3','c3','r12.pos','r12.neg','tau','mu','d'),useConjugacy = TRUE)

#build MCMC
mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(mcmc, project = model,resetFunctions = TRUE)

#run the model 
st<-Sys.time()
samples=runMCMC(compiled_mcmc,inits=inits,
                nchains = 1, nburnin=250000,niter = 500000,samplesAsCodaMCMC = TRUE,thin=20,
                summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st

#trace plots
mcmcplot(samples)


#####################################
#Statewide strata data
#####################################
#statewide data for poststratification
strata.data=join(t.pop.data[,c('GEOID','age.group','sex.group','pop')],region.data[,c('GEOID','region','st.log.pop')],by='GEOID',type='left')
#statewide design matrix
X=model.matrix(~1+age.group+sex.group+st.log.pop,data=strata.data,contrasts=list(age.group='contr.sum',sex.group='contr.sum'))

#####################################
#Calculate posterior
#####################################
#mcmc object to dataframe
data.samples=as.data.frame(samples)

#calculate response rate by region
response.data$response.rate=response.data$RESPONDERS/response.data$TOTAL

#constants
#number of MCMC samples
reps=dim(data.samples)[1]
#number of regions in state
regions=length(unique(strata.data$region))
#number of census tracts in state
n.tracts=length(unique(strata.data$GEOID))

#tract index for state
tracts=as.numeric(as.factor(strata.data$GEOID))

#non-response RR
RR=c(0.9,1,1.1,1.25,1.50,2,2.5,3)
num=length(RR)

#initialize storage
cases=array(NA,dim=c(reps,regions,num))
total.cases=array(NA,dim=c(reps,num))
rate=array(NA,dim=c(reps,num))

#calculate posterior by region and strata for each RR
for(k in 1:num){
  for(i in 1:reps){
    
    #tract SD
    ct.sd=data.samples[i,'tau']
    
    #tract RE
    b.ct=rnorm(n.tracts,0,ct.sd)
    
    #beta
    beta=data.samples[i,c('beta[1]','beta[2]','beta[3]','beta[4]')]
    
    for(j in 1:regions){
      
      #indicator of strata and region
      ind=(strata.data$region==j)
      
      #calculate probability estimate
      mu=X[ind,-1]%*%t(beta)+data.samples[i,paste0('alpha[',j,']')]+b.ct[tracts[ind]]
      
      #back transform
      est.p=exp(mu)/(1+exp(mu))
      
      #calculate adjusted probability
      adj.p=est.p*response.data$response.rate[response.data$region==j]+RR[k]*est.p*(1-response.data$response.rate[response.data$region==j])
      
      #calculate expected number of cases
      cases[i,j,k]=sum(adj.p*strata.data$pop[ind])
      
    }
  }
  #sum cases over region and strata
  total.cases[,k]=rowSums(cases[,,k])
  
  #compute rate for state
  rate[,k]=total.cases[,k]/sum(strata.data$pop)
}

#posterior mean and HPD CI for each RR
post.mean=colMeans(rate)
cred.int=HPDinterval(as.mcmc(rate))

