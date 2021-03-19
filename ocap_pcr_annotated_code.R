######################################
#
# OCAP - Ohio COVID-19 Prevalence
# Annotated Code
# PCR Test - Current Infection
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

#PCR test 3=negative 4=positive
ocap.data$PCR=NA
ocap.data$PCR[ocap.data$pcr_taqpath_covid_19==4]='Positive'
ocap.data$PCR[ocap.data$pcr_taqpath_covid_19==3]='Negative'

#####################################
#Data for analysis
#####################################
#select variables from survey
ocap.vars=c('PIN','GEOID','age.group','sex.group','PCR')

#complete survey data limited to variables needed
ocap.analysis=ocap.data[complete.cases(ocap.data[,ocap.vars]),ocap.vars]

#merge in region ID
ocap.analysis=join(ocap.analysis,region.data[,c('GEOID','REGION_ML')],by='GEOID',type='left')

#add population - standardized log population
ocap.analysis=join(ocap.analysis,pop.data[,c('GEOID','st.log.pop')],by='GEOID',type='left')

#constants for MCMC
#number of subjects
subjects=length(unique(ocap.analysis$PIN))
#number of regions
regions=length(unique(ocap.analysis$REGION_ML))
#number of tracts
tracts=length(unique(ocap.analysis$GEOID))

#test result
test.result=NA
test.result[ocap.analysis$PCR=='Positive']=1
test.result[ocap.analysis$PCR=='Negative']=0

#fixed effects matrix
X=model.matrix(~1+age.group+sex.group+st.log.pop,data=ocap.analysis,contrasts=list(age.group='contr.sum',sex.group='contr.sum'))

#####################################
#Model
#####################################
#nimble model code
model_code=nimbleCode({
  
  #test result model
  for(i in 1:n){
    
    #test result
    t[i] ~ dbern(p[i])
    p[i]<- d[i]*s+(1-d[i])*(1-c)
    
    #prevalence model
    d[i] ~ dbern(pi[i])
    logit(pi[i])<- inprod(X[i,1:xp],beta[1:xp])+alpha[region[i]]+b[tract[i]]
    
  }
  
  #process model
  #random region effects
  for(i in 1:r){
    alpha[i] ~ dnorm(mu,sd=sigma)
  }
  
  #tract random effects
  for(i in 1:ct){
    b[i] ~ dnorm(0,sd=tau)
  }
  
  
  #fixed effects
  for(i in 1:xp){
    beta[i] ~ dnorm(0,sd=3)
  }
  
  mu ~ dnorm(log(0.01/0.99),sd=0.5)
  
  sigma ~ dunif(0,5)
  
  #tract effects
  tau ~ dunif(0,5)
  
  #test characteristics
  s ~ dbeta(890,110)
  c ~ dbeta(995,5)
  
})

#####################################
#NIMBLE
#####################################
#constants
#number of fixed effect parameters - removing intercept
xp=dim(X[,-1])[2]

#model inputs
constants=list(n=subjects,xp=xp,r=regions,ct=tracts,
               region=as.numeric(ocap.analysis$REGION_ML),tract=as.numeric(as.factor(ocap.analysis$GEOID)))

#test results and design matrix without intercept
model.data=list(t=test.result,X=X[,-1])

#initial values
inits=list(alpha=rep(0,regions),sigma=1,beta=rep(0,xp),d=rep(0,subjects),p=rep(0.1,subjects),
           b=rep(0,tracts),tau=1,mu=-3)

#load limited data to run code
load('ocap_pcr_limited_data.Rda')

# Build the model.
model <- nimbleModel(model_code, constants,model.data,inits)
compiled_model <- compileNimble(model,resetFunctions = TRUE)

# Set up samplers.
mcmc_conf <- configureMCMC(model,monitors=c('alpha','sigma','beta','tau','mu','c','s'),useConjugacy = TRUE)

#build MCMC
mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(mcmc, project = model,resetFunctions = TRUE)

# Run the model 
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


