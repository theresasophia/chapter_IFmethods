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
# close graphics windows
###Install tgern doMPI in R 
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
library(dplyr)
library(doParallel)
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

# measurement model 
dmeas <- Csnippet("if (ISNA(cases)) {
                  lik = (give_log) ? 0 : 1;
                  } else {
                  lik =  dnbinom_mu(cases, 1/od, H, 1) ;
                  lik = (give_log) ? lik : exp(lik);}")

rmeas <-  Csnippet("cases = rnbinom_mu(1/od,H);")



# transmission model is Markovian SIRS with 3 age classes, seasonal forcing and overdispersion
sir.step <- Csnippet("double rate[7];
                     double dN[7];
                     double Beta1;
                     
                     Beta1 = beta1*(1 + beta11 * cos(M_2PI/52*t + phi));
                     rate[0] = mu*N;
                     rate[1] = Beta1*I/N;
                     rate[2] = mu;
                     rate[3] = gamma;
                     rate[4] = mu;
                     rate[5] = omega;
                     rate[6] = mu;
                     
                     dN[0] = rpois(rate[0]*dt); // births are Poisson
                     reulermultinom(2, S, &rate[1], dt, &dN[1]);
                     reulermultinom(2, I, &rate[3], dt, &dN[3]);
                     reulermultinom(2, R, &rate[5], dt, &dN[5]);
                     S += dN[0] - dN[1] - dN[2] + dN[5];
                     I += dN[1] - dN[3] - dN[4];
                     R += dN[3] - dN[6] - dN[5];
                     H += dN[1];
")


# ------------ deterministic skeleton-----------------------------
sir.skel <- "
                     double rate[7];
                     double term[7];
                     double Beta1;
                     
                     
                     Beta1 = beta1*(1 + beta11 * cos(M_2PI/52*t + phi)); //seasonal forcing
                     
                     rate[0] = mu*N;
                     rate[1] = Beta1*(I)/N;;
                     rate[2] = mu;
                     rate[3] = gamma;
                     rate[4] = mu;
                     rate[5] = omega;  
                     rate[6] = mu;  
                     
                     
                     
                     // compute the several terms
                     term[0] = rate[0];
                     
                     term[1] = rate[1] * S;
                     term[2] = rate[2] * S;
                     
                     term[3] = rate[3] * I;
                     term[4] = rate[4] * I;
                     
                     term[5] = rate[5] * R;
                     term[6] = rate[6] * R;
                     
                     DS = term[0] - term[1] - term[2] + term[5];
                     DI = term[1]          - term[3] - term[4];
                     DR = term[3]          - term[6] - term[5];
                     DH = term[1];
                     
                     " 

# read in the data
# add at t=0 a row of NAs to not have problems with the accumulator variable since
# t0 is much less than t[1]
read.table("total_rota.txt")->dat_temp
dat_temp[,-(3:5)]%>%
  rbind(data.frame(time=0,cases=NA)) %>%
  arrange(time) -> dat


# define parameters (without betas)
params_fixed <- c(gamma=1, mu=1/(78.86912*52), N=82372825, omega=1/(1*52))
first_data <- c(y=dat$cases[2])


# initializer
init <- function(params, t0, ...) {
  x0 <- c(S=0,I=0,R=0,H=0)
  y <- params[c("y")]
  x0["I"] <- y[1]/((params["gamma"]+params["mu"]))
  x0["S"] <- (params["mu"]*params["N"]-(params["gamma"]+params["mu"])*x0["I"]+
                 params["omega"]*(params["N"]*params["mu"]/params["mu"]-x0["I"]))/(params["mu"]+params["omega"])

  x0["R"] <- (params["N"]*params["mu"]/params["mu"]-x0["S"]-x0["I"])

  round(x0) 
}

#help parameters with different data ie the mean data
mean_data <- c(y=mean(dat$cases[-1])) 
help_param <- c(params_fixed,mean_data)
# analytic guess for the betas
beta_ana <-  function(params){
  beta_ana <- c(beta1=0)
  I_bar <- init(params)["I"]
  beta_ana["beta1"] <- ((params["gamma"]+params["mu"])*init(params)["I"]*params["N"])/(init(params)["S"]*init(params)["I"])

  return(beta_ana)
}

# paramtervector with betas and inital data 
params <- c(beta_ana(help_param), beta11=0.15, phi=0.1, params_fixed,first_data,od=0.3, sigma=0.05)

#transformation of parameter space for unbounded search
toEst<- Csnippet("
                 Tbeta1  = log(beta1);
                 Tbeta11 = logit(beta11);
                 Tsigma = log(sigma);
                 Tphi    = logit(phi/(M_2PI));
                 Tod = log(od);")

fromEst <-Csnippet("
                   Tbeta1  = exp(beta1);
                   Tsigma = exp(sigma);
                   Tbeta11 = expit(beta11);
                   Tphi    = M_2PI*expit(phi);
                   Tod = exp(od);")

#pomp object
pomp(data = dat,
     times="time",
     t0=1-6*52,
     dmeasure = dmeas,
     rmeasure = rmeas,
     rprocess = euler.sim(step.fun = sir.step, delta.t = 1/10),
     statenames = c("S", "I", "R", "H"),
     paramnames = names(params),
     zeronames=c("H"),
     skeleton=vectorfield(Csnippet(sir.skel)),
     initializer=init,
     toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     params = params
) -> sir


################### GLOBAL SEARCH ############################
require(doParallel)
cores <- 20
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)


#hypercube with sensible parameters
sir_box <- rbind(
  beta1=c(0,2),
  beta11=c(0.11,0.16),
  phi=c(0.01,0.3),
  od=c(0.001,0.3),
  sigma=c(0.001,0.2)
)


# parameters which are assued to be fixed and known and dont have to estimated
sir_fixed_params <- c(params_fixed, y=dat$cases[2])

stew(file="mle_gamma_data_agg.rda",{
  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=function(...){x<-c(...);
    saveRDS(x, file = "tmp.rds");x}, .options.multicore=mcopts) %dopar% {
      mif2(
        sir,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=5000,
        Nmif=300,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(beta1=0.002,beta11=0.001,phi=0.01,od=0.01, sigma=0.01)
      )
    }
  })
},seed=1270401374,kind="L'Ecuyer")


#particle filter evaulations of the output vetor of the mif2 searches, remove filtering failures
stew(file="mle_gamma_data_lik_agg.rda",{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.errorhandling= "remove",.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=5000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

# choose the set of parameters with hightest loglik as MLE (since errorhabdling=remove, we have to reset the number)
best <- which.max(liks_global[,1])
best <- 22



#plot the results
coef(sir) <- coef(mifs_global[[best]])
sim = simulate(sir, nsim=1000, states=TRUE,obs=TRUE,seed=1234)

periods <- (length(dat$time[-1]))/52
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
df2 <- as.data.frame(t(cbind(quantile_obs[ , -1, 1])))
df4 <- as.data.frame(t(cbind(quantile_states[ , -1, 1])))
names(df4)[1:3] <- c("lowernnb", "midnnb", "uppernnb")
df2$age <- factor(rep(c("children"), each = dim(quantile_obs)[2] - 1),
                  levels = c("children"))

df2$time <- rep(seq.int(dim(quantile_obs)[2] - 1), 1)
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
dfcases <-melt(dat[-1,], id.vars = c("time"), variable.name = "cases")
df2$cases <- dfcases$value
df2$labelcases <- "Data"
df2$labelcases <- factor(df2$labelcases, levels = c("Data"))

ggplot(df2) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill=type1), alpha = 0.15) +
  geom_ribbon(aes(x = time, ymin = lowernnb, ymax = uppernnb, fill=type), alpha = 0.35) +
  geom_line(aes(x = time, y = mid, color = labelmedian),color='white',size=1) +
  #geom_line(aes(x = time, y = midnnb, color = labelmediannnb),color='red',linetype=3,size=1) +
  geom_line(aes(x = time, y = cases, color = labelcases),color='black',size=1) +
#  facet_wrap( ~age, ncol=1, scales =  "free_y") +
  scale_x_continuous( breaks = axis.spots,labels = axis.labels)+ ylab("Weekly new cases")+ xlab("Time (weeks)")+
  labs(color="")+scale_fill_grey( start = 0.1, end = 0.2)+
theme(text = element_text(size=45))+ theme(legend.position="none")


