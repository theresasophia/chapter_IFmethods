######################################################################
## R code generating graphs and numbers for the Analysis of 
## Outbreak Data section of the book chapter "Iterated filtering 
## methods for Markov process epidemic models"  by T. Stocks in the
## Handbook of Infectious Disease Data Analysis.
##
## Author: Theresa Stocks <http://www.su.se/english/profiles/tstoc-1.219526>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 7/12/2017
######################################################################


rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows

###Install packages
library(pomp)
library(dplyr)
library(readr)
library(dplyr)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
require(doMPI)
require(doParallel)
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")
require(doParallel)
cores <- detectCores()
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)
 

#parallelization of the code
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

#read in the data which was simulated from the model
dat <- read.table("sim_data.txt")

#plot the data
ggplot(dat,mapping=aes(x=time,y=cases))+
  geom_bar(stat="identity")+guides(color=FALSE)+ ylab("Weekly new reported cases")+ 
  xlab("Time (weeks)")+theme(text = element_text(size=45))  

#determinsitc SIR transmission model
sir.skel <-  Csnippet("DS = - Beta*I/N*S;
                       DI = Beta*I/N*S - gamma*I;
                       DR = gamma*I;
                       DH = Beta*I/N*S;")

#stochastic SIR transmission model
sir.step <- Csnippet("double rate[2];
                      double dN[2];
                     rate[0] = Beta*I/N;
                     rate[1] = gamma;
                     reulermultinom(1, S, &rate[0], dt, &dN[0]);
                     reulermultinom(1, I, &rate[1], dt, &dN[1]);
                     S += -dN[0];
                     I += dN[0] - dN[1];
                     R += dN[1];
                     H += dN[0];
                     ")

#intial values 
init <- Csnippet("S = N-1;
                   I = 1;
                   R = 0;
                   H = 0;")

#observation model
dmeas <-Csnippet("lik =  dpois(cases, H, 1);
                 lik = (give_log) ? lik : exp(lik);")

rmeas <- Csnippet("cases = rpois(H);")


#transmformation of all bounded parameters to unbounded
fromEst <- Csnippet("TBeta = exp(Beta);
                     Tgamma = exp(gamma);")
toEst <- Csnippet(" TBeta = log(Beta);
                    Tgamma = log(gamma);")

#parameter vector
params <- c(Beta=1, gamma=1/2, N=10000)

#pomp object
pomp(data=dat,     
     times = 'time',
     t0 = 1,
     rprocess = euler.sim(step.fun = sir.step, delta.t = 1/10),
     skeleton = vectorfield(sir.skel),
     initializer = init,
     statenames = c("S","I","R","H"),
     zeronames = c("H"),
     paramnames=names(params),
     toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     dmeasure = dmeas,
     rmeasure = rmeas,
     params = params
) -> sir

##  generating the data from above
# #simulate data 
# plot(simulate(sir,obs=TRUE,seed=123),type='l')
# x <-simulate(sir,obs=TRUE,seed=123)
# dat <- data.frame(time=seq(1:length(x[,,])),cases=x[,,])
# write.table(dat,file="~/Dropbox/Bookchapter/Rcode/sir_simple/sim_data.txt" )


#fixing the parameters we do not want to estimate (here population size N)
sir_fixed_params<- params["N"]

#global search: inital guesses are drawn from hypercube
sir_box <- rbind(
  Beta=c(0,2),
  gamma= c(0,1)
)


########################################################################################## 
# Fitting the model with stochastic underlying transmission model to the simulated data  #   
# with iterated filtering                                                                #
##########################################################################################

#number of particles is N=1000 and number of mif2 iterations is M=100
stew(file="sim_data_mif.rda",{
  t_global <- system.time({
    mifs_global <- foreach(i=1:10,.packages='pomp', .combine=function(...){x<-c(...);saveRDS(x, file = "tmp.rds");x}, .options.multicore=mcopts) %dopar% {
      mif2(
        sir,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=100,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02,gamma=0.02)
      )
    }
  })
},seed=1270401374,kind="L'Ecuyer")

# diagnostic plots of the iterated filtering search
mifs_global %>%
  conv.rec(c("loglik", "nfail","Beta","gamma"))%>%
  melt() %>%
  mutate(variable = factor(variable)) %>%
  mutate(variable = recode(variable, Beta = "beta")) %>%
  mutate(variable = recode(variable, gamma = "gamma"))-> df
df %>%
  subset(iteration>0)%>%
  ggplot(aes(x=iteration,y=value,color=variable, group=L1))+
  geom_line(size=1)+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+scale_colour_grey( start = 0.3, end = 0.3)+
  facet_wrap(~variable,scales="free_y",ncol=2,labeller = label_parsed)+
theme(text = element_text(size=45))

#evaluation of the logliklihood of the final mif outputs with particle filter
stew(file="mif_eval.rda",{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:10,.packages='pomp',.errorhandling="remove",.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=2000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

#choose the vector with highest loglik as MLE
best <- which.max(liks_global[,1])


coef(sir)<-coef(mifs_global[[best]])

#############################################################################
# calculations of PROFILE LIKELIHOODS                                        #
#############################################################################

##PROFILE DESIGNS for beta
profileDesign(
  Beta=seq(from=0.9,to=1.15,length=30),
  lower=c(gamma=0,sir_fixed_params ),
  upper=c(gamma=1, sir_fixed_params ),
  nprof=2
) -> pd_Beta
pd <- pd_Beta


#independent mif2 searches in order to construct a profile likelihood
stew(file="profile_beta_mif.rda",{
  t_global <- system.time({
    mifs_global <- foreach(p=iter(pd,"row"),.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir,
        start=unlist(p),
        Np=1000,
        Nmif=100,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(gamma=0.02)
      )
    }
  })
},seed=1270401374,kind="L'Ecuyer")


#calculation of pfilter for ecah mif2 output,.errorhandling=pass to find the filtering failures
stew(file="profile_beta_lik1.rda",{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:dim(pd)[1],.errorhandling="pass",.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

#extract filtering failures
rownames(liks_global)<- seq(1:nrow(liks_global))
ind<-as.numeric(rownames(liks_global[which(liks_global[,"call"] == "NULL"), ]))

# repeat exact same calculations but now .errorhandling=remove to find the filtering failures
stew(file="profile_beta_lik2.rda",{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:dim(pd)[1],.errorhandling="remove",.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


results_beta<- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global[-ind],coef)))


##PROFILE DESIGNS for gamma
profileDesign(
  gamma=seq(from=0.4,to=0.6,length=30),
  lower=c(Beta=0,sir_fixed_params ),
  upper=c(Beta=2, sir_fixed_params ),
  nprof=2
) -> pd_gamma

pd <- pd_gamma



#mif2 search
stew(file="profile_gamma_mif.rda",{
  
  t_global <- system.time({
    mifs_global <- foreach(p=iter(pd,"row"),.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir,
        start=unlist(p),
        Np=1000,
        Nmif=100,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02)
      )
    }
  })
},seed=1270401374,kind="L'Ecuyer")

#evaluate particle filter and have errorhandling=pass to identify location of filtering failure
stew(file="profile_gamma_lik1.rda",{
  
  t_global_eval <- system.time({
    liks_global_gamma <- foreach(i=1:dim(pd)[1],.errorhandling="pass",.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

#extract rownumber of filtering failure
rownames(liks_global_gamma)<- seq(1:nrow(liks_global_gamma))
ind<-as.numeric(rownames(liks_global_gamma[which(liks_global_gamma[,"se"] == "NULL"), ]))

#reeat exact same calculations but errorhandling=remove
stew(file="profile_gamma_lik2.rda",{
  
  t_global_eval <- system.time({
    liks_global_gamma2 <- foreach(i=1:dim(pd)[1],.errorhandling="remove",.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


results_gamma<- data.frame(logLik=liks_global_gamma2[,1],logLik_se=liks_global_gamma2[,2],t(sapply(mifs_global[-ind],coef)))

####### now with the results we do the monte carlo adjusted profile likelihood cf Ionides 2017

mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000){
  smooth_fit <- loess(lp ~ parameter,span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0,length(parameter))
  weight[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"] ) 
  a <- unname(coef(quadratic_fit)["a"] ) 
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a) 
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp,parameter=parameter,confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(
         parameter=parameter_grid,
         smoothed=smoothed_loglik,
         quadratic=predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
       ),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), se=sqrt(se_total_squared)
  )
}


lambda_spacetime <- 1 
spacetime_mcap_beta <- mcap(lp=results_beta$logLik,parameter=results_beta$Beta,lambda=lambda_spacetime)
spacetime_mcap_gamma <- mcap(lp=results_gamma$logLik,parameter=results_gamma$gamma,lambda=lambda_spacetime)

plot.profile <- function(spacetime_mcap,true,xlab){
  data_res<- data.frame(param=spacetime_mcap$parameter, lp=spacetime_mcap$lp)
  data_lp<- data.frame(x=spacetime_mcap$fit$parameter, b=spacetime_mcap$fit$smoothed,d=spacetime_mcap$fit$quadratic)
  
  ggplot(data=data_res,aes(x=param,y=lp))+geom_point(size=2)+geom_line(data=data_lp,aes(x=x,y=b),color="grey25",size=2.5) +
    geom_vline(xintercept = spacetime_mcap$ci[1],linetype=2,size=2) +
    geom_vline(xintercept = spacetime_mcap$ci[2],linetype=2,size=2)+geom_line(data=data_lp,aes(x=x,y=d),color="grey50",size=1.5,linetype=2.5)+
    geom_line(aes(x=param,y= max(spacetime_mcap$fit$smoothed,na.rm=T)-spacetime_mcap$delta),size=1.5)+ geom_vline(xintercept = true,size=2) +
    xlab(xlab)+ ylab("Loglikelihood")+theme(text = element_text(size=75)) +
    scale_y_continuous(limits = c(-160, -154))
}

plot.profile(spacetime_mcap_beta,1,expression(beta))
plot.profile(spacetime_mcap_gamma,0.5,expression(gamma))


