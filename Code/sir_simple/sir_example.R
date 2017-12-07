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
 

#read in the data which was simulated from the model
dat <- read.table("sim_data.txt")

#plot the data
cairo_ps(file="~/Dropbox/Bookchapter/Latex/data.eps", width=17, height=9)
ggplot(dat,mapping=aes(x=time,y=cases))+
  geom_bar(stat="identity")+guides(color=FALSE)+ ylab("Weekly new reported cases")+ 
  xlab("Time (weeks)")+theme(text = element_text(size=45))  
dev.off()

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

#check if the stochastic model is working by simulating from it
plot(simulate(sir))
#check if the right data is loaded into the pomp object by plotting it
plot(sir)


##  generating the data from above
# #simulate data 
# plot(simulate(sir,obs=TRUE,seed=123),type='l')
# x <-simulate(sir,obs=TRUE,seed=123)
# dat <- data.frame(time=seq(1:length(x[,,])),cases=x[,,])
# write.table(dat,file="~/Dropbox/Bookchapter/Rcode/sir_simple/sim_data.txt" )

#check if the determistic skeleton is working as well 
trajectory(sir,params=params,as=TRUE) -> x
plot(x$H,type='l')

#parallelization of the code
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

#fixing the parameters we do not want to estimate (here population size N)
sir_fixed_params<- params["N"]

#global search: inital guesses are drawn from hypercube
sir_box <- rbind(
  Beta=c(0,2),
  gamma= c(0,1)
)



########################################################################################## 
#Fitting the model with deterministic underlying transmission model to the simulated data#
##########################################################################################
#number of searches
n<- 1000

#trajectory matching where errors are removed
stew(file="sim_dat_traj.rda",{
  
  t_global <- system.time({
    guesses <- as.data.frame(apply(sir_box,1,function(x)runif(n,x[1],x[2])))
    mle_traj <- foreach(guess=iter(guesses,"row"),.errorhandling="remove",.packages='pomp', .combine=c, 
                        .options.multicore=mcopts) %dopar% {
                          traj.match(sir,
                                     start=c(unlist(guess),sir_fixed_params),
                                     est=c( "Beta", "gamma"),
                                     method="Nelder-Mead",reltol=1e-8,maxit=4000, transform=TRUE)
                          
                        }
  })
},seed=1270401374,kind="L'Ecuyer")


#extract the failures and set n to the successful number of searches
n=as.numeric(length(mle_traj))
#extract the logliklihood and choose the MLE as the parameter vector with highest loglik
f <- function(x) {summary(mle_traj[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)

#MLE and loglik of MLE
round(coef(mle_traj[[best]]),2)
summary(mle_traj[[best]])$loglik

# check if all searches converged (convergence is indicated with 0)
h<- function(x) {summary(mle_traj[[x]])$convergence}
conv <- apply(x,2,FUN=h)
length(which(conv==0))

#extract the coeficients
g <- function(x) {coef(mle_traj[[x]])}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
coef <- apply(x,2,FUN=g)
g<-t(coef[,which(conv==0)])
results_global_all <- data.frame(logLik=liks,gamma=coef["gamma",],Beta=coef["Beta",])
results_global <- results_global_all[which(conv==0),]
guesses_con <- guesses[which(conv==0),] 


# starting values and corresponding estimates which are in a 90 % quantile of the loglik
quan<-0.90
all <- ldply( list(guess=guesses_con[which(results_global$logLik > quantile(results_global$logLik,probs = quan)),], 
                   results_global=subset(results_global, logLik > quantile(results_global$logLik,probs = quan)) ), .id="type")
plot(liks[which(conv==0)])

#pair plot- in grey the starting value of the algorithm and in red the estimated value 
pairs(~logLik+gamma+Beta,data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
sum(liks < quantile(results_global$logLik,probs = quan))


#plot 95% in-sample prediction interval of the model evaluated at the MLE 
# (data black solid line, mean white solid line, prediction interval red shading)
trajectory(sir,params=coef(mle_traj[[best]]),as=TRUE) -> x

nfu<-c(qpois(0.975,x$H))
nfl<-c(qpois(0.025,x$H))
df2 <- data.frame(lower= nfl, upper=nfu)
df2$time <- dat$time
df2$type <- "95% prediction interval"
df2$type <- factor(df2$type, levels = c("95% prediction interval"))
df2$labelmedian <- "Mean"
df2$labelmedian <- factor(df2$labelmedian, levels = c("Mean"))
df2$mean <- c(x$H)
dfcases <-melt(dat, id.vars = c("time"), variable.name = "cases")
df2$cases <- dat$cases
df2$labelcases <- "Data"
df2$labelcases <- factor(df2$labelcases, levels = c("Data"))

#pdf(file="~/Dropbox/Cuba/pred_without.pdf", width=17, height=9)
ggplot(df2) + 
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill=type), alpha = 0.3) + 
  geom_line(aes(x = time, y = mean, color = labelmedian),color='white',size=1)+
  geom_line(aes(x = time, y = cases, color = labelcases),color='black') +
  ylab("Weekly new cases")+ xlab("Time")+  labs(color="")+ ylab("Number of new cases")+ 
  xlab("Time (weeks)")+  labs(color="")+theme_bw()+
  theme(legend.position="none",text = element_text(size=15))
#dev.off()


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
cairo_ps(file="~/Dropbox/Bookchapter/Latex/mif.eps", width=17, height=9)
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
dev.off()


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

#simulate 10 realizations form the model evaluated at the MLE
sim = simulate(sir, nsim=1000, states=TRUE,obs=TRUE,seed=1234)

periods <- (length(dat$time))/52
axis.spots <- (0:(periods))*52+2;
axis.labels <- as.character((2001):(2001+periods));

## x is time, y is variable
quantile_obs_fun <- function(x){
  quantile(sim$obs[,,x],probs=c(0.025,0.5,0.975))
}

quantile_states_fun <- function(x,y){
  quantile(sim$states[,,x][y,],probs=c(.025,0.5,.975))
}

# quantiles,time, variable
quantile_obs <- array(0,c(3,length(dat$time),1),dimnames=list(c("2.5%","50%","97.5%"),NULL, c("Children")))
quantile_obs[,,"Children"] <- mapply(quantile_obs_fun,seq(1:length(dat$time)))


quantile_states <- array(0,c(3,length(dat$time),1),dimnames=list(c("2.5%","50%","97.5%"),NULL, c("Children")))
quantile_states[,,"Children"]<-mapply(quantile_states_fun,seq(1:length(dat$time)),"H")


#preparing the data fram
library(tidyr)
df2 <- as.data.frame(t(cbind(quantile_obs[ , , 1])))
df4 <- as.data.frame(t(cbind(quantile_states[ , , 1])))
names(df4)[1:3] <- c("lowernnb", "midnnb", "uppernnb")
df2$age <- factor(rep(c("children"), each = dim(quantile_obs)[2] ),
                  levels = c("children"))

df2$time <- rep(seq.int(dim(quantile_obs)[2] ), 1)
names(df2)[1:3] <- c("lower", "mid", "upper")
head(df2)
df2$lowernnb <-df4$lowernnb
df2$midnnb <-df4$midnnb
df2$uppernnb <-df4$uppernnb
df2$type <- "95% PI process model"
df2$type <- factor(df2$type, levels = c("95% PI process model"))
df2$type1 <- "95% PI  model"
df2$type1 <- factor(df2$type1, levels = c("95% PI  model"))
df2$labelmedian <- "Median total model"
df2$labelmedian <- factor(df2$labelmedian, levels = c("Median total model"))
df2$labelmediannnb <- "Median process model"
df2$labelmediannnb <- factor(df2$labelmediannnb, levels = c("Median process model"))
dfcases <-melt(dat[,], id.vars = c("time"), variable.name = "cases")
df2$cases <- dfcases$value
df2$labelcases <- "Data"
df2$labelcases <- factor(df2$labelcases, levels = c("Data"))

ggplot(df2) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill=type1), alpha = 0.2) +
  geom_ribbon(aes(x = time, ymin = lowernnb, ymax = uppernnb, fill=type), alpha = 0.5) +
  geom_line(aes(x = time, y = mid, color = labelmedian),color='white',size=1) +
  geom_line(aes(x = time, y = midnnb, color = labelmediannnb),color='red',linetype=3,size=1) +
  geom_line(aes(x = time, y = cases, color = labelcases),color='black') +
  facet_wrap( ~age, ncol=1, scales =  "free_y") +
  scale_x_continuous( breaks = axis.spots,labels = axis.labels)+ ylab("Weekly new cases")+
  xlab("Time (weeks)")+  labs(color="")+
  theme_bw()+theme(text = element_text(size=30))+ theme(legend.position="none")

sims <- simulate(sir,params=coef(mifs_global[[best]]),
                 nsim=1000,as.data.frame=TRUE,include.data=TRUE)

ggplot(sims,mapping=aes(x=time,y=H,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)


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
dim(pd)

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

#check for convergence
mifs_global[2] %>%
  conv.rec(c("loglik","nfail","Beta", "gamma")) %>%
  melt() %>%subset(iteration>0)%>%
  ggplot(aes(x=iteration,y=value,color=variable,group=L1))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()

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
ind
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
dim(pd)


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

#diagnostic plot
mifs_global[1] %>%
  conv.rec(c("loglik","nfail","Beta", "gamma")) %>%
  melt() %>%subset(iteration>0)%>%
  ggplot(aes(x=iteration,y=value,color=variable,group=L1))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()

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
ind

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

cairo_ps(file="~/Dropbox/Bookchapter/Latex/profile_beta.eps", width=17, height=14)
plot.profile(spacetime_mcap_beta,1,expression(beta))
dev.off()
cairo_ps(file="~/Dropbox/Bookchapter/Latex/profile_gamma.eps", width=17, height=14)
plot.profile(spacetime_mcap_gamma,0.5,expression(gamma))
dev.off()

