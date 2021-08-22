#Lynn Waterhouse
#Based on J. Hoenig structure written 26 Aug 2019
#created Aug 30, 2019
#modified January 3, 2020 to match what John wrote.
#Updated Dec 7, 2020 to not use fish lengths, just spine data
#Updated Jan 5, 2020 for simulations

#set working directory
#setwd("C:/Users/Lynn Waterhouse/Desktop/BayesSpine_originalmodel_notgoodBayes")
setwd("C:/Users/water/Desktop/Jan2021_BRIAN_datascience_spineonly_sims/300fish_unifages_0.2sdgrowth_100sims")

#set seed so reproducible
set.seed=12345

#Sample will have fish ages 1 to 10
min.age<-1   #min age 1
max.age<-10  #max.age 10
all.ages<-seq(min.age,max.age,by=1)  #vector of all ages
no.perage<-30  #how many fish per age do we want (roughly)

#growth von-Bert K=0.4, Linf=4, T0=-0.1
#Length(age_t)=Linf*(1-exp(K_growth*(age_t-t0)))
#Length(age_t)=Linf*(1-exp(K_growth*(Obs.ring+Miss.ring-t0)))
K_growth<- 0.4
Dinf<- 4
t0<- -0.1
#when generating lengths add random error growth.err~rnorm(0,3)
#mean 0, sd 3 (based on Then et al. paper)
gr.mean=0
gr.sd=0.2  #Med sd growth

#assume vascularization is % of overall spine. we'll assume 75%
#rad_vasc=0.6*rad_spine
p.vasc<-0.6  #based on mean from Guelson's YFT data which was 0.5456
#add random error ~rnorm(0,.3)
va.mean<-0
va.sd<-sqrt(0.003)  #variance by estimated age ranged from 0.001 to 0.0036 from Guelson's YFT data

set.seed(12345)
nsims<-100
for(my.sim in 1:nsims){

  
  ##############################
  ##############################
  ####### Generate data ########
  ##############################
  remove(data.age, data.max, rings_obs, rings_miss, true_age,
         rad_spine, rad_spine.age, rad_spine_exp, growth.err, fish.cens,
         fish.ring.min, fish.ring.miss,calc.age, calc.age.med,
         calc.age.mode,true.mat,true.mat2,toy_basic,ring.rad, ring.rad.obs,
         repp, rdatatemp, obs.mat, obs.mat2, naive_basic, jags.file, 
         gelman, dat, a)
  #create lookup matrix of age, exp(length), exp(spine dia.)
  data.age<-NULL
  data.age<-as.data.frame(matrix(data=rep(NA,length(all.ages)*2),
                                 ncol=2,byrow=TRUE))
  colnames(data.age)<-c("age","rad_spine.age")
  #fill ages
  data.age$age<-all.ages
  for (i in 1:length(all.ages)){
    data.age$rad_spine.age[i]=Dinf*(1-exp(-K_growth*(data.age$age[i]-t0)))
  }
  
  
  #create matrix for simulated data
  data.max<-NULL
  data.max<-as.data.frame(matrix(data=rep(NA,no.perage*max.age*7),
                                 nrow=no.perage*max.age, ncol=7,byrow=TRUE))
  colnames(data.max)<-c("true_age","rad_spine","rad_spine_exp","growth.err","rad_vasc","rings_obs",
                        "rings_miss")
  #fill ages
  #randomly sample 200 fish ages from possible ages
  data.max$true_age<-sample(all.ages,size=length(all.ages)*no.perage, replace=T) 
  ring.rad<-NULL
  ring.rad<-as.data.frame(matrix(data=rep(NA,no.perage*max.age*max.age),
                                 nrow=no.perage*max.age, ncol=max.age,byrow=TRUE))
  rings_miss2<-NULL
  
  #fill lengths
  attach(data.max)
  attach(data.age)
  
  for (i in (1:dim(data.max)[1])){
    #expected spine
    rad_spine_exp[i]=Dinf*(1-exp(-K_growth*(true_age[i]-t0))) 
    #individual error on growth (spine)
    growth.err[i]=rnorm(1,mean=gr.mean,sd=gr.sd)
    #individual spine
    rad_spine[i]=rad_spine_exp[i]+growth.err[i]  
    #dummy term in logit space to get variability added to ratio of vascularization to spine 
    rad_vasc_p_log=log(p.vasc/(1-p.vasc))+rnorm(1,mean=va.mean,va.sd)
    #idividual spine radius
    rad_vasc[i]= (exp(rad_vasc_p_log)/(1+exp(rad_vasc_p_log)))*rad_spine[i]
    #p.vasc*rad_spine[i]+rnorm(1,mean=va.mean,va.sd) former way of doing this
    #add for loop for creating a matrix that has the ring radius for each individual at each age
    for(j in 1:true_age[i]){   
      ring.rad[i,j]<-rad_spine.age[j]*rad_spine[i]/rad_spine.age[true_age[i]]
    }
    #how many rings are going to be missing (missing is <= vascularized region)
    rings_miss[i]<-max(which(ring.rad[i,]<=rad_vasc[i])) #using differing growth rates (using this)
    #rings_miss2[i]<-max(which(data.age$rad_spine.age<=rad_vasc[i])) #no growth considered
    #make zeros when no rings missing
    if(is.infinite(rings_miss[i])) rings_miss[i]<-0 #USING THIS
    #how many rings do we see for each individual
    rings_obs[i]=true_age[i]-rings_miss[i]
    
    #if(is.infinite(rings_miss2[i])) rings_miss2[i]<-0
  }

  
  data.max$growth.err<-growth.err
  data.max$rad_spine<-rad_spine
  data.max$rad_spine_exp<-rad_spine_exp
  data.max$rad_vasc<-rad_vasc
  data.max$rings_obs<-rings_obs
  data.max$rings_miss<-rings_miss
  head(data.max)  
  
  #generate observed ring radius - nullify any ages for fish with vascularization
  #rings disappear if smaller than or equal to vascularized region
  #keep rings that are bigger than (outside of) vascularized region
  ring.rad.obs<-NULL
  ring.rad.obs<-as.data.frame(matrix(data=rep(NA,no.perage*max.age*max.age),
                                     nrow=no.perage*max.age, ncol=max.age,byrow=TRUE)) #empty matrix
  for (i in 1:dim(ring.rad)[1]){
    for (j in 1:true_age[i]){
      if(ring.rad[i,j]<=rad_vasc[i]){ring.rad.obs[i,j]<-NA    } #NA if ring is vascularized that far
      if(ring.rad[i,j]>rad_vasc[i]){ring.rad.obs[i,j]<- ring.rad[i,j]   } #actual value if no vascularization
    }
  }
  
  #create vector of 1st ring observed per fish
  fish.ring.min<-NULL #vector of minimum ring size observed
  for(i in 1:dim(ring.rad.obs)[1]){
    if(sum(!is.na(ring.rad.obs[i,1:true_age[i]]))>0){
      temp<-min(ring.rad.obs[i,1:true_age[i]],na.rm=T)}
    if(sum(!is.na(ring.rad.obs[i,1:true_age[i]]))==0){
      temp<-NA
    }
    fish.ring.min<-c(fish.ring.min,temp)
  }
  length(fish.ring.min) #should be nfish long
  
  ##Histogram of 1st ring radius seen on each spine
  #hist(fish.ring.min,breaks=100,xlab="diameter 1st ring")
  #summary(fish.ring.min) 
  
  #create a variable that filters between fish with rings obscured and those without
  #use the minimum ring diamater observed as a 'tau' variable
  m.tau<-min(fish.ring.min)  #smallest ring diameter observed
  m.tau
  fish.ring.miss<-NULL #this will be 0 if not censored, 1 if censored
  for(i in 1:length(fish.ring.min)){
    if(rad_vasc[i]<m.tau){fish.ring.miss[i]<-0} #not censored
    if(rad_vasc[i]>=m.tau){fish.ring.miss[i]<-1} #censored
  }
  
  #look at how many we think are censored versus not
  table(fish.ring.miss) 
  
  #make a new variable that is NA is censored, 0 if not
  fish.cens<-fish.ring.miss
  fish.cens[which(fish.cens==1)]<-NA
  table(fish.cens)
  summary(fish.cens)
  
  
  data.max$rings_miss<-rings_miss
  data.max$fish.ring.min<-fish.ring.min
  data.max$fish.ring.miss<-fish.ring.miss
  data.max$fish.cens<-fish.cens
  head(data.max) 
  
  write.csv(data.max,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                                "sim.data.nsim",my.sim,".csv",sep=""))
  
  #library(ggplot2)
  ##Look at some of the plots
  #plot(true_age,rad_spine)  #age vs. spine plot
  #plot(true_age,rad_vasc)  #age vs. vasc plot
  #plot(rad_spine,rad_vasc)  #spine vs. vasc plot
  #ggplot(data.max, aes(x=true_age, y=rad_spine)) +
  #  geom_point(aes(color=fish.cens))
  #plot(data.max$rings_obs,data.max$rad_spine) #observed rngs vs. spine plot
  #plot(data.max$true_age,data.max$rings_miss) #age vs. missing rings
  
  #ggplot(data.max, aes(x=true_age, y=rings_miss)) +
  #  geom_count()+ scale_size_area()
  
  #plot(data.max$true_age,data.max$rings_obs)
  #ggplot(data.max, aes(x=true_age, y=rings_obs)) +
  # geom_count()+ scale_size_area()
  
  #hist(true_age,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11))
  #hist(data.max$rings_obs,breaks=c(0,1,2,3,4,5,6,7,8,9))
  #hist(rad_vasc)
  #hist(rad_spine)
  #hist(rad_vasc/rad_spine)
  
  #plot(rad_vasc,rad_spine)
  
  #check summary of vascularization:spine radius
  summary(rad_vasc/rad_spine)
  

  ########################################################################
  ########################################################################
  ########################################################################
  ################################
  #Bayes model to fit missing data
  ########################################################################
  #now fit using the posteriors
  #attach(data.maxt)
  attach(data.age)
  nfish<-dim(data.max)[1]
  
  summary(rad_vasc/rad_spine)
  #hist((rad_vasc/rad_spine))
  #hist(rlnorm(1000,mean=log(mean(rad_vasc/rad_spine)),sd=sd(rad_vasc/rad_spine)))
  
  save.image(,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                         "Simulation.data.toy_bayes.spine.nsim",my.sim,".RDATA",
                         sep=""))
  

  ########################################################################
  ########################################################################
  ########################################################################
  ################################
  #Bayes model to fit missing data
  ########################################################################
  #now fit using the posteriors
  #attach(data.maxt)
  attach(data.age)
  nfish<-dim(data.max)[1]
  
  #summary(rad_vasc/rad_spine)
  #hist((rad_vasc/rad_spine))
  #hist(rlnorm(1000,mean=log(mean(rad_vasc/rad_spine)),sd=sd(rad_vasc/rad_spine)))
  
  save.image(,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                         "Simulation.data.toy_bayes.spine.nsim",my.sim,".RDATA",
                         sep=""))
  
  # Age probability Matrix
  # makes it only possible for the fish to be an age equal to the number of observed rings or greater
  max.age2<-max(rings_obs)*2
  repp<-matrix(rep(NA,(max.age2+1)*nfish),nrow=nfish)
  for (i in 1:nfish){
    try<-c( c(rep(0, (rings_obs[i]) )),
            rep( (1/(max.age2-rings_obs[i]+1)),(max.age2-rings_obs[i]+1) ))
    repp[i,]<-try
    # CHANGE TO MAKE --- add an if statement - so that if we know no rings are missing it accounts for this
  }
  
  ###################################################################
  # below builds the model file:
  ###################################################################
  
  cat("
      model{
      m.Dinf~dunif(0,30)                 #prior for Dinf -- VB 
      m.K~dunif(0,1)                       #prior for growth K -- VB 
      m.t0~dunif(-1,1)                     #prior for length age zero -- VB  
      m.blank~dunif(0,1)                  #prior for proportion vascularized
      
      spine.mu<-0                #prior for mean of spine error
      sigma.spine ~ dunif(0, 1)           # Prior for sd of spine error
      sigma2.spine <- pow(sigma.spine, 2)  #variance of spine error
      tau.spine <- pow(sigma.spine, -2)    #one over variance of spine error
      
      vasc.mu<-0                #prior for mean of vasc error
      sigma.vasc ~ dunif(0, 1)           # Prior for sd of vasc error
      sigma2.vasc <- pow(sigma.vasc, 2)  #variance of vasc error
      tau.vasc <- pow(sigma.vasc, -2)    #one over variance of vasc error
      
      for(i in 1:nfish){
      #expected spine of actual fish
      mu.spine[i]<-m.Dinf*(1-exp(-m.K*(true.age[i]-m.t0))) #VB equation
      exp.spine[i]~dnorm(mu.spine[i]+spine.mu, tau.spine) #random error
      true.age[i]~dcat(c(repp[i, 2:(max.age2+1)])) #discrete uniform prior on ages of fish
      #true.age[i]<-exp.rings_miss[i]+rings_obs[i]
      }
      
      for(iii in 1:nfish){
      #vascularization
      #add in logit error
      #err.rad_vasc[iii]~dnorm(vasc.mu, tau.vasc) #error for vascularization
      #m.blank.lg[iii]<-log(m.blank/(1-m.blank))+err.rad_vasc[iii] #add error in logit space
      #m.blank.sp[iii]<-exp(m.blank.lg[iii])/(1+exp(m.blank.lg[iii])) #back transform out of logit space
      #exp.rad_vasc[iii]~dnorm(m.blank.sp[iii]*exp.spine_rad[iii],0) #expected vasularized radius
      mean.rad_vasc[iii]<-m.blank*exp.spine[iii]
      exp.rad_vasc[iii]~dnorm(mean.rad_vasc[iii],tau.vasc) #if we ignore the random error
      }
      
      for(iV in 1:nfish){
      #age #CHANGED 01.03.20- added temp, changes poiss statment
         #temp[iV]<-ifelse((1-(m.a+m.b*exp.spine_rad[iV])/(m.Linf))>=1e-100, 
         #   (1-(m.a+m.b*exp.spine_rad[iV])/(m.Linf)), 
         #    exp(-m.K*(rings_obs[iV]-m.t0))) 
         #swap to rad_vasc for obliterated
      temp[iV]<-ifelse((1-(exp.rad_vasc[iV])/(m.Dinf))>=1e-100, 
      (1-(exp.rad_vasc[iV])/(m.Dinf)), 
      exp(m.K*m.t0) ) #for obliterated    
      
      #temp2[iV]<-log(temp[iV])/(-m.K) + m.t0 - rings_obs[iV]
      temp2[iV]<-log(temp[iV])/(-m.K)+m.t0 #if using obliterated region
      
      exp.rings_miss[iV]~dpois(temp2[iV]) 
      #expected missing rings, come from poisson with mean of missing rings for fish that age
      }
      
      #end model    
      }",file="toy_bayes.spine.txt")

  dat=list("rings_obs"=rings_obs, "exp.rad_vasc"=rad_vasc, "exp.spine"=rad_spine, 
           "max.age2"=max.age2,"nfish"=nfish,"repp"=repp,"exp.rings_miss"=fish.cens)

  parameters<-c("m.Dinf","m.K","m.t0","sigma.spine",
                "exp.rings_miss","m.blank","sigma.vasc") 
  model="toy_bayes.spine.txt"
  
  library(runjags)
  library(R2jags)
  library(gtools)
  library(gdata)
  ###############################  
  ###############################
  ###############################
  #Jags parameters for later on
  my.chains = 8		  #number of chains (4) 8
  my.iter = 200000 	#number of iterations in each chain (25,000) 200,000
  my.burnin = 4000 	#burnin per chain (5000) 5000
  my.thin=40    #thining (40) 400
  #At the end how many MCMC results will we have
  my.mcmc.count<-my.chains*(my.iter-my.burnin)/my.thin 
  my.mcmc.count
  ################################
  # run the model in JAGS, using default settings
  
  toy_basic <- jags(data=dat, inits=NULL, parameters.to.save=parameters, 
                    model.file=model, n.chains = my.chains, n.iter = my.iter, 
                    n.burnin = my.burnin, n.thin=my.thin, DIC = TRUE, 
                    progress.bar="text")
  
  attach.jags(toy_basic)
  
  ##future modification
  ##need to look into how to add constraint to a parameter so that true age is sum of the other 2 things
  
  save.image(,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                         "Simulation.toy_bayes.spine.nsim",my.sim,".RDATA", sep=""))
  
  
  
  
  
  ###################################################################
  ###################################################################
  # Looking at Bayes Results
  
  #Check out VB growth params
  #K_growth<- 0.4
  quantile(m.K,c(0.5, 0.025, 0.975))
  #Dinf<- 4
  quantile(m.Dinf,c(0.5, 0.025, 0.975))
  #t0<- -0.1
  quantile(m.t0,c(0.5, 0.025, 0.975))
  
  #check out a,b, and m.blank
  quantile(m.blank,c(0.5, 0.025, 0.975)) #.6

  write.csv(rbind(
    c("K",quantile(m.K,c(0.5, 0.025, 0.975))),
    c("Dinf",quantile(m.Dinf,c(0.5, 0.025, 0.975))),
    c("t0",quantile(m.t0,c(0.5, 0.025, 0.975))),
    c("beta",quantile(m.blank,c(0.5, 0.025, 0.975)))
  ), file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                "VBparams_toys_basic.nsim",my.sim,".csv",sep=""))  
  
  #function to calculate mode
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  #summary stats on posteriors
  rings_post.med<-apply(exp.rings_miss,2,median)
  rings_post.mode<-apply(exp.rings_miss,2,getmode)
  rings_post.mean<-apply(exp.rings_miss,2,mean)
  rings_post.UCI<-apply(exp.rings_miss,2,quantile,probs=.975)
  rings_post.LCI<-apply(exp.rings_miss,2,quantile,probs=.025)
  rings_post.UCI2<-apply(exp.rings_miss,2,quantile,probs=.9)
  rings_post.LCI2<-apply(exp.rings_miss,2,quantile,probs=.1)
  
  #how many times did it get it right - over all fish
  table(rings_post.med==rings_miss)
  table(rings_post.mode==rings_miss)
  table(round(rings_post.mean,digits=0)==rings_miss)
  table(floor(rings_post.mean)==rings_miss)
  table(rings_miss>=rings_post.LCI&rings_miss<=rings_post.UCI)
  table(rings_miss2>=rings_post.LCI&rings_miss<=rings_post.UCI2)
  
  write.csv(rbind(
    c("median",table(rings_post.med==rings_miss)),
    c("mode",table(rings_post.mode==rings_miss)),
    c("round_mean",table(round(rings_post.mean,digits=0)==rings_miss)),
    c("floor_mean",table(floor(rings_post.mean)==rings_miss)),
    c("95CI",table(rings_miss>=rings_post.LCI&rings_miss<=rings_post.UCI)),
    c("80CI",table(rings_miss>=rings_post.LCI2&rings_miss<=rings_post.UCI2)),
    c("naive_model",table(rings_miss==0))
  ), file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                "Overall_AgeRightWrong_BRIAN.nsim",my.sim,".csv",sep=""))  
  

  
  
  #how many times did it get it right - over censored fish
  table(rings_post.med[is.na(fish.cens)]==rings_miss[is.na(fish.cens)])
  table(rings_post.mode[is.na(fish.cens)]==rings_miss[is.na(fish.cens)])
  
  #how many times did it get it right - over fish we "knew" age
  table(rings_post.med[which(fish.cens==0)]==rings_miss[which(fish.cens==0)])
  table(rings_post.mode[which(fish.cens==0)]==rings_miss[which(fish.cens==0)])
  
  #How many times was the truth in Bayes CI
  table(rings_miss<=rings_post.UCI&rings_miss>=rings_post.LCI)
  
  #Look at some histograms
  true_age<-data.max$true_age
  #hist(true_age, main="true age")
  #hist(rings_post.mode+rings_obs,main="Posterior Mode Estimated Age")
  #hist((rings_post.mode+rings_obs)-true_age,main="Estimated Age-True Age",
  #     xlab="Estimated Age - True Age")
  
  #hist((rings_obs)-true_age,main="Observed Ring Age-True Age",
  #     xlab="Observed Ring  Age - True Age")
  
  #Calculate the final age (missing rings plus observed rings)
  calc.age<-matrix(rep(NA,length(exp.rings_miss)),ncol=ncol(exp.rings_miss),
                   nrow=nrow(exp.rings_miss))
  for (nn in 1:length(rings_miss)){
    calc.age[,nn]<-exp.rings_miss[,nn]+rings_obs[nn]
  }
  calc.age.med<-apply(calc.age,2,median) #median posterior estimate age
  calc.age.mode<-apply(calc.age,2,getmode) #mode of posterior estimated age
  
  
  
  #look at how it did by age- using median posterior
  calc.age.med<-round(rings_post.med)+rings_obs #need to round the results if using median
  #We will use the mode since it is a Poisson distribution
  true.mat<-matrix(rep(NA,max(calc.age.mode)*max(true_age)),
                   ncol=max(calc.age.mode), nrow=max(true_age))
  true.mat2<-true.mat
  for(i in 1:max(true_age)){
    for(j in 1:max(calc.age.mode)){
      true.mat[i,j]<-sum(calc.age.mode[which(true_age==i)]==j)
      true.mat2[i,j]<-true.mat[i,j]/sum(true_age==i)
    }
  }
  true.mat
  round(true.mat2, digits=3)
  round(true.mat2, digits=2)
  
  write.csv(true.mat, file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                            "BayesPostAgecol_TrueAgerow.nsim",my.sim,".csv"), row.names = TRUE)
  write.csv(true.mat2, file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                             "propsBayesPostAgecol_TrueAgerow.nsim",my.sim,".csv"), row.names = TRUE)
  
  ##look at how it did by age-- for just rings observations
  ##This would be what happens if you naively assume the ring count is the true age
  obs.mat<-matrix(rep(NA,max(rings_obs)*max(true_age)),
                  ncol=max(rings_obs), nrow=max(true_age))
  obs.mat2<-obs.mat
  for(i in 1:max(true_age)){
    for(j in 1:max(rings_obs)){
      obs.mat[i,j]<-sum(rings_obs[which(true_age==i)]==j)
      obs.mat2[i,j]<-obs.mat[i,j]/sum(true_age==i)
    }
  }
  
  obs.mat
  round(obs.mat2,digits=2)
  
  
  write.csv(obs.mat,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                               "NaiveAgecol_TrueAgerow.nsim",my.sim,".csv"), row.names = TRUE)
  write.csv(obs.mat2,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                                "propsNaiveAgecol_TrueAgerow.nsim",my.sim,".csv"), row.names = TRUE)
  

  ######################################################################
  library(mcmcplots)
  #mcmc plots
  #mcmcplot(toy_basic) #if you run this it takes forever
  #Gelman-rubin statistics
  #want this to be below 1.1
  # citation: https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
  #setwd(mywd)
  jags.file<-toy_basic
  mcmcjags.file<-as.mcmc(jags.file)
  n.var <- coda::nvar(mcmcjags.file)
  gelman <- matrix(NA, nrow=n.var, ncol=2)
  for (v in 1:n.var) {gelman[v,] <- coda::gelman.diag(mcmcjags.file[,v],autoburnin = FALSE )$psrf }
  which(gelman[,1]>1.05)
  
  func.gelman<-function(jags.file){
    mcmcjags.file<-as.mcmc(jags.file)
    n.var <- coda::nvar(mcmcjags.file)
    gelman <- matrix(NA, nrow=n.var, ncol=2)
    for (v in 1:n.var) {gelman[v,] <- coda::gelman.diag(mcmcjags.file[,v],autoburnin = FALSE )$psrf }
    gelman.all <- gelman[which(!is.nan(gelman[,1])),] # Remove dummy variables (show up as NA) 
    gelman_short <- gelman[order(gelman[,1],decreasing=T),] 
    if(n.var>10) gelman_short <- gelman_short[1:10,] 
    gelman_fail <- c(length(which(gelman[,1]>1.01)), length(which(gelman[,1]>1.05)), 
                     length(which(gelman[,1]>1.1))) 
    gelman_fail_CI <- c(length(which(gelman[,2]>1.01)), length(which(gelman[,2]>1.05)), 
                        length(which(gelman[,2]>1.1)))
    #failing is >1.05, so middle number
    use.DIC=jags.file$BUGSoutput$DIC
    return(list("gelman_fail"=gelman_fail,"gelman_fail_CI"=gelman_fail_CI,
                "gelman_short"=gelman_short,
                "gelman.all"=gelman.all,"numbervar"=n.var,"DIC"=use.DIC))
  }
  
  a<-func.gelman(toy_basic)
  a
  
  sink(paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
             "gelmanrubin_toy_basic.nsim",my.sim,".csv",sep=""))
  print(a)
  sink()
  
  
  
  ############################################################################
  table(rings_post.med)
  table((rings_post.med+rings_obs)-true_age) 
  table(fish.cens)  #reminder of how many fish are censored, 0=not censored, rest are NA=yes censor
  
  #where is it getting it wrong
  table(rings_miss[which((rings_post.med != rings_miss))])
  
  #where is it getting it right
  table(rings_miss[which((rings_post.med == rings_miss))])
  
  #reminder of what were original missing values
  table(rings_miss)
  
  
  
  ############################################################################
  ############################################################################
  #### Let's do the naive thing and fit VB to the observed ring counts #######
  ############################################################################
  table(rings_obs) #observed 'ages' if we use observed rings
  table(true_age) #compare with true ages
  
  table(rings_obs,true_age) #rows observed age, columns true age
  obs.mat #rows true age, columns observed age
  
  write.csv(table(rings_obs,true_age),file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                                                 "ObsAgeRingObsrow_TrueAgecol.nsim",my.sim,".csv"), 
            row.names = TRUE)
  write.csv(table(true_age,rings_obs),file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
                                                 "ObsAgeRingObscol_TrueAgerow.nsim",my.sim,".csv"), 
            row.names = TRUE)
  ###################################################################
  # below builds the model file:
  ###################################################################
  
  cat("
      model{
      m.Dinf~dunif(0,30)                   #prior for Dinf -- VB 
      m.K~dunif(0,1)                         #prior for growth K -- VB 
      m.t0~dunif(-1,1)                       #prior for length age zero -- VB  
      
      spine.mu<-0                           #prior for mean of spine error
      sigma.spine~ dunif(0, 10)            #Prior for sd of spine error
      sigma2.spine <- pow(sigma.spine, 2)  #variance of spine error
      tau.spine <- pow(sigma.spine, -2)    #one over variance of spine error
      
      
      for(i in 1:nfish){
      #expected spine of actual fish
      mu.spine[i]<-m.Dinf*(1-exp(-m.K*(rings_obs[i]-m.t0))) #VB equation
      exp.spine[i]~dnorm(mu.spine[i]+spine.mu, tau.spine) #random error
      }
      #end model    
      }",file="toy_bayes.spine.naive.txt")

  dat=list("rings_obs"=rings_obs,"nfish"=nfish,"exp.spine"=rad_spine)
  
  parameters<-c("m.Dinf","m.K","m.t0","sigma.spine")
  model="toy_bayes.spine.naive.txt"
  
  library(runjags)
  library(R2jags)
  library(gtools)
  library(gdata)
  ###############################
  #Jags parameters for later on
  my.chains = 8		  #number of chains (4) 8
  my.iter = 200000 	#number of iterations in each chain (25,000) 200,000
  my.burnin = 4000 	#burnin per chain (5000) 5000
  my.thin=40    #thining (40) 400
  #At the end how many MCMC results will we have
  my.mcmc.count<-my.chains*(my.iter-my.burnin)/my.thin 
  my.mcmc.count
  ################################
  # run the model in JAGS, using default settings
  
  naive_basic <- jags(data=dat, inits=NULL, parameters.to.save=parameters, 
                      model.file=model, n.chains = my.chains, n.iter = my.iter, 
                      n.burnin = my.burnin, n.thin=my.thin, DIC = TRUE, 
                      progress.bar="text")
  
  attach.jags(naive_basic)
  
  ###################################################################
  ###################################################################
  # Looking at Bayes Results
  
  #Check out VB growth params
  #K_growth<- 0.4
  quantile(m.K,c(0.5, 0.025, 0.975))
  #Dinf<- 100
  quantile(m.Dinf,c(0.5, 0.025, 0.975))
  #t0<- -0.1
  quantile(m.t0,c(0.5, 0.025, 0.975))
  
  write.csv(rbind(
    c("K",quantile(m.K,c(0.5, 0.025, 0.975))),
    c("Dinf",quantile(m.Dinf,c(0.5, 0.025, 0.975))),
    c("t0",quantile(m.t0,c(0.5, 0.025, 0.975)))
      ), file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
         "VBparams_naive_basic.nsim",my.sim,".csv",sep=""))
  
  ######################################################################
  library(mcmcplots)
  #mcmc plots
  #mcmcplot(naive_basic) #if you run this it takes forever
  #Gelman-rubin statistics
  #want this to be below 1.1
  # citation: https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
  #setwd(mywd)
  jags.file<-naive_basic
  mcmcjags.file<-as.mcmc(jags.file)
  n.var <- coda::nvar(mcmcjags.file)
  gelman <- matrix(NA, nrow=n.var, ncol=2)
  for (v in 1:n.var) {gelman[v,] <- coda::gelman.diag(mcmcjags.file[,v],autoburnin = FALSE )$psrf }
  which(gelman[,1]>1.05)
  
  #func.gelman from before
  
  a<-func.gelman(naive_basic)
  a
  
  sink(paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),
             "gelmanrubin_naive_basic.nsim",my.sim,".csv",sep=""))
  print(a)
  sink()
  
  
  save.image(,file=paste(format(Sys.time(), "%Y.%m.%d.%I.%p."),"Simulation.toy_bayes.spine.nsim",
                         my.sim,".RDATA",sep=""))

}

