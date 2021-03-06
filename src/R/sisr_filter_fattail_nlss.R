#################################################################################################
#
#  FAT-TAILED NON-LINEAR STATE SPACE MODEL 
#
#  SEQUENTIAL IMPORTANCE SAMPLING WITH RESAMPLING
#
#  MODEL:
#  y[t]      =  theta[t]  +  sqrt(lambda[t])*v[t]
#  theta[t]  =  beta*(theta[t-1]/(1+theta[t-1]^2))  +  w[t]
#
#  lambda[t]    ~ IG(nu/2,nu/2)
#  v[t]         ~ N(0,v)
#  y | theta[t] ~ T_nu(theta[t],v)    student T distribution with nu degrees of freedom
#  theta[t] | theta[t-1]  ~  N(beta*(theta[t-1]/(1+theta[t-1]^2)),w)
#
#################################################################################################
# Author :   Rounak Meyur
#            Power and Energy Center Virginia Tech
#            Blacksburg, Virginia 24060
#            Email: rounakm8@vt.edu
#################################################################################################


# True parameters of the model
# -----------------------------
nu = 2
beta = 0.9
v = 0.04
w = 0.01
Tbig = 200

# User defined functions for fast computation
# --------------------------------------------
G           <- function(x){x/(1+x^2)}
quant025    <- function(x){quantile(x,0.025)}
quant975    <- function(x){quantile(x,0.975)}

# Simulating the observation and state data
# ------------------------------------------
set.seed(12345)
sim.y         <- rep(0,Tbig)
sim.theta     <- rep(0,Tbig)
sim.lambda    <- 1/rgamma(Tbig,nu/2,nu/2)
sim.theta0    <- 0.1

g             <- G(sim.theta0)
sim.theta[1]  <- rnorm(1,g*beta,sqrt(w))
sim.y[1]      <- rnorm(1,sim.theta[t],sqrt(sim.lambda[1]*v))
for(t in 2:Tbig)
{
  g            <- G(sim.theta[t-1])
  sim.theta[t] <- rnorm(1,g*beta,sqrt(w))
  sim.y[t]     <- rnorm(1,sim.theta[t],sqrt(sim.lambda[t]*v))
}

# Plot the time series
# ---------------------
par(mfrow=c(2,2))
plot(sim.y,type="l",xlab="time,t",ylab="Observations,y")
plot(sim.theta,type="l",xlab="time,t",ylab=TeX("State,$\\theta$"))

plot(sim.theta,sim.y,type="p",pch=19,xlab=TeX("State,$\\theta$"),ylab="Observations,y")
plot(sim.theta[1:199],sim.theta[2:200],type="p",pch=19,xlab=TeX("State,$\\theta_{t-1}$"),
     ylab=TeX("State,$\\theta_{t}$"))


# Sequential importance sampling with resampling (SISR)
# -----------------------------------------------------
smc.sisr <- function(Tbig,M,y,theta0,wt0)
{
  est.theta <- NULL
  ess <- NULL
  
  theta <- theta0
  wt <- wt0
  
  for(t in 1:Tbig)
  {
    m.theta <- beta*t(apply(cbind(theta),1,G))
    theta <- rnorm(M,m.theta,sqrt(w))
    
    wt <- dt(y[t],df=nu,ncp=theta)
    wt <- wt/sum(wt)
    
    # Estimated states
    est.theta <- rbind(est.theta,sample(theta,size=M,
                                        replace=T,prob=wt))
    # Estimated sample size
    ess <- rbind(ess,1/sum(wt^2))
  }
  return(list("theta"=est.theta,"ess"=ess))
}


# Estimate the states
# --------------------
set.seed(12545)
M <- 1000
m0 <- 0
C0 <- 100

theta0 <- rnorm(M,m0,sqrt(C0))
wt0 <- rep(1/M,M)

est <- smc.sisr(Tbig,M,sim.y,theta0,wt0)
# smt <- smoothing(Tbig,M,est$theta)
# sisr.stheta = apply(smt,1,mean)
sisr.mtheta = apply(est$theta,1,mean)
sisr.ltheta = apply(est$theta,1,quant025)
sisr.utheta = apply(est$theta,1,quant975)
sisr.es <- est$ess


# Plot the estimates
# ------------------
par(mfrow=c(2,1))
plot(sisr.mtheta,xlab="Time",ylab="State estimates",
     main=expression(theta_t),type = 'l',col="blue")
lines(sim.theta,col="red",lty=2)
# lines(sisr.stheta,col=1,lty=1)
plot(sisr.es,type="l",ylab="Effective sample size")
