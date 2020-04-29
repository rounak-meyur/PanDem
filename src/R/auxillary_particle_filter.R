#################################################################################################
#
#  NON-LINEAR NORMAL DYNAMIC MODEL 
#
#  AUXILIARY PARTICLE FILTER 
#
#################################################################################################
# Author :   Rounak Meyur
#            Power and Energy Center Virginia Tech
#            Blacksburg, Virginia 24060
#            Email: rounakm8@vt.edu
#################################################################################################

# True parameters of the model
a = 1/20
b = 1/2
c = 25
d = 8
beta <- c(b,c,d)
omega = 1.2
v = 10
w = 1
Tbig = 200

# User defined functions for fast computation
G <- function(arg){c(arg[1],arg[1]/(1+arg[1]^2),cos(omega*arg[2]))}
quant025    <- function(x){quantile(x,0.025)}
quant975    <- function(x){quantile(x,0.975)}

# Simulating the observation and state data
set.seed(12345)
sim.y         <- rep(0,Tbig)
sim.theta     <- rep(0,Tbig)
sim.theta0    <- 0.1
sim.t0        <- 0.0

g             <- G(c(sim.theta0,sim.t0))
sim.theta[1]  <- rnorm(1,sum(g*beta),sqrt(w))
sim.y[1]      <- rnorm(1,a*sim.theta[1]^2,sqrt(v))
for(t in 2:Tbig)
{
  g            <- G(c(sim.theta[t-1],t-1))
  sim.theta[t] <- rnorm(1,sum(g*beta),sqrt(w))
  sim.y[t]     <- rnorm(1,a*sim.theta[t]^2,sqrt(v))
}

# Plot the time series
par(mfrow=c(2,2))
plot(sim.y,type="l",xlab="time,t",ylab="Observations,y")
plot(sim.theta,type="l",xlab="time,t",ylab=TeX("State,$\\theta$"))

plot(sim.theta,sim.y,type="p",pch=19,xlab=TeX("State,$\\theta$"),ylab="Observations,y")
plot(sim.theta[1:199],sim.theta[2:200],type="p",pch=19,xlab=TeX("State,$\\theta_{t-1}$"),
     ylab=TeX("State,$\\theta_{t}$"))


# Auxillary particle filtering
# -----------------------------
smc.apf <- function(Tbig,M,y,theta0,wt0)
{
  # Initialize the arrays of estimates
  est.theta <- NULL
  ess <- NULL
  
  theta <- theta0
  wt <- wt0
  
  # Iterate over the time instants
  for (t in 1:Tbig)
  {
    X <- t(apply(cbind(theta,t-1),1,G))
    m.theta <- apply(X*beta,1,sum)
    theta.prior <- rnorm(M,m.theta,sqrt(w))
    wt <- dnorm(y[t],a*theta.prior^2,sqrt(v))*wt
    k <- sample(M,size=M,replace=T,prob=wt/sum(wt))
    
    theta <- rnorm(M,m.theta[k],sqrt(w))
    wt <- dnorm(y[t],a*theta^2,sqrt(v))/dnorm(y[t],
                                              a*m.theta[k]^2,sqrt(v))  
    wt <- wt/sum(wt)
    
    # Find estimated state through resampling
    est.theta <- rbind(est.theta,sample(theta, M, 
                                        replace = T, prob = wt))
    # Estimated sample size
    ess <- rbind(ess,1/sum(wt^2))
  }
  return(list("theta"=est.theta,"ess"=ess))
}


# Smoothing function
# -------------------
smoothing <- function(Tbig,M,theta)
{
  theta.smooth <- array(0,dim=c(Tbig,M))
  theta.smooth[Tbig,] <- theta[Tbig,]
  for(t in (Tbig-1):1)
  {
    X <- t(apply(cbind(theta[t,],t-1),1,G))
    m.theta <- apply(X*beta,1,sum)
    wt <- dnorm(theta.smooth[t+1,],m.theta,sqrt(v))
    wt <- wt/sum(wt)
    theta.smooth[t,] <- sample(theta[t,],size=M,
                               replace=T,prob=wt)
  }
  return(theta.smooth)
}


# Estimate the states
# --------------------
set.seed(12545)
M <- 1000
m0 <- 0
C0 <- 100

theta0 <- rnorm(M,m0,sqrt(C0))
wt0 <- rep(1/M,M)

est <- smc.apf(Tbig,M,sim.y,theta0,wt0)
smt <- smoothing(Tbig,M,est$theta)
apf.stheta = apply(smt,1,mean)
apf.mtheta = apply(est$theta,1,mean)
apf.ltheta = apply(est$theta,1,quant025)
apf.utheta = apply(est$theta,1,quant975)
apf.es <- est$ess


# Plot the estimates
# ------------------
par(mfrow=c(2,1))
plot(apf.mtheta,xlab="Time",ylab="State estimates",
     main=expression(theta_t),type = 'l', col="blue")
lines(sim.theta,col="red",lty=2)
lines(apf.stheta,col=1,lty=1)
plot(apf.es,type="l",ylab="Efective sample size")
