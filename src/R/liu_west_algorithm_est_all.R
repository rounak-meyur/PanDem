#################################################################################################
#
#  NON-LINEAR NORMAL DYNAMIC MODEL 
#
#  AUXILIARY PARTICLE FILTER + PARAMETER LEARNING 
#
#################################################################################################
#
# Author : Hedibert Freitas Lopes
#          The University of Chicago Booth School of Business
#          5807 South Woodlawn Avenue
#          Chicago, Illinois, 60637
#          Email : hlopes@Chicagobooth.edu
#
#################################################################################################
# Modified : Rounak Meyur
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





# Liu and West filter
# -------------------

# Define discount factor
# -------------------
delta   <- 0.99
h2      <- 1-((3*delta-1)/(2*delta))^2
A       <- sqrt(1-h2)

# Sample theta and phi from particle approximations
# -------------------------------------------------
sample.theta.phi <- function(M,y,t,theta,phi)
{
  # Compute phi mean and phi variance
  m.phi <- apply(phi,2,mean)
  v.phi <- var(phi)
  
  # Compute auxillary indices
  X <- t(apply(cbind(theta,t-1),1,G))
  theta.prior <- apply(X*phi[,1:3],1,sum)
  phi.prior <- A*phi+(1-A)*matrix(m.phi,M,5,byrow=T)
  
  # Compute auxillary indices
  p <- dnorm(y[t],a*theta.prior^2,exp(phi.prior[,4]/2))
  k <- sample(1:M,size=M,replace=T,prob=p)
  
  # Sample theta and phi
  phi.post <- phi.prior[k,] + matrix(rnorm(5*M),M,5)%*%chol(h2*v.phi)
  X <- t(apply(cbind(theta[k],t-1),1,G))
  theta.post <- rnorm(M,apply(X*phi.post[,1:3],1,sum),exp(phi.post[,5]/2))
  weight <- dnorm(y[t],a*theta.post^2,exp(phi.post[,4]/2))/
    dnorm(y[t],a*theta.prior[k]^2,exp(phi.prior[k,4]/2))
  weight <- weight/sum(weight)
  ind    = sample(1:M,size=M,replace=T,prob=weight)
  
  return(list("theta"=theta.post[ind],"phi"=phi.post[ind,]))
}

# Sequential Monte Carlo
# -----------------------
smc.liu <- function(Tbig,M,y,theta0,phi0)
{
  # Initialize the arrays of estimates
  est.theta <- NULL
  est.phi <- array(0,c(M,5,Tbig))
  
  theta <- theta0
  phi <- phi0
  
  # Iterate over the time instants
  for (t in 1:Tbig)
  {
    samples <- sample.theta.phi(M,y,t,theta,phi)
    
    theta <- samples$theta
    phi <- samples$phi
    
    # Update
    est.theta <- rbind(est.theta,theta)
    est.phi[,,t] <- phi 
  }
  mtheta = apply(est.theta,1,mean)
  ltheta = apply(est.theta,1,quant025)
  utheta = apply(est.theta,1,quant975)
  mphi  = matrix(0,Tbig,5)
  lphi  = matrix(0,Tbig,5)
  uphi  = matrix(0,Tbig,5)
  for (i in 1:3)
  {
    mphi[,i] = apply(est.phi[,i,],2,mean)
    lphi[,i] = apply(est.phi[,i,],2,quant025)
    uphi[,i] = apply(est.phi[,i,],2,quant975)
  }
  for (i in 4:5)
  {
    mphi[,i] = apply(exp(est.phi[,i,]),2,mean)
    lphi[,i] = apply(exp(est.phi[,i,]),2,quant025)
    uphi[,i] = apply(exp(est.phi[,i,]),2,quant975)
  }
  
  return(list("theta_mu"=mtheta,"theta_l"=ltheta,"theta_u"=utheta,
              "phi_mu"=mphi,"phi_l"=lphi,"phi_u"=uphi))
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

# Liu and West filter
# -------------------
set.seed(8642)
M <- 10000

# Prior hyperparameters
# ---------------------
a0  = 3
A0  = 20
b0  = 3
B0  = 2
c0  = c(0.5,25,8)
C0  = c(0.1,16,2)
m0  = 0
V0  = 5

phi0 <- cbind(rnorm(M,c0[1],sqrt(C0[1])),rnorm(M,c0[2],sqrt(C0[2])),rnorm(M,c0[3],sqrt(C0[3])),
              log(1/rgamma(M,a0,A0)),log(1/rgamma(M,b0,B0)))
theta0 <- rnorm(M,m0,sqrt(V0))


# Estimate the states and parameters
# -----------------------------------
est <- smc.liu(200,M,sim.y,theta0,phi0)


mx <- est$theta_mu
lx <- est$theta_l
ux <- est$theta_u

mpars <- est$phi_m
lpars <- est$phi_l
upars <- est$phi_u

fil.error <- abs(sim.theta-mx)


# Plot the estimates and quantiles
# --------------------------------
pars.true <- c(beta,v,w)

layout(matrix(c(1,2), 2, 1, byrow = TRUE))
plot(mx,xlab="Time",ylab="State estimates",
     main=TeX("State $\\theta_t$ estimates after filtering"),
     type = 'b',col="blue",pch=20,cex.main=1.5)
lines(sim.theta,col=1,lty="dotted")
plot(fil.error,type="h",ylab="Error in filtering estimation",
     xlab="Time",main="Filtering error",cex.main=1.5)


names = c("Parameter: b","Parameter: c","Parameter: d","Parameter: v","Parameter: w")
par(mfrow=c(1,3))
for (i in 1:3){
  ts.plot(mpars[,i],ylim=range(lpars[,i],upars[,i]),ylab="Estimated parameter",
          main=names[i],col="black")
  lines(lpars[,i],col="blue")
  lines(upars[,i],col="blue")
  abline(h=pars.true[i],col="red",lty="dotted")
}


par(mfrow=c(1,2))
for (i in 4:5){
  ts.plot(mpars[,i],ylim=range(lpars[,i],upars[,i]),ylab="Estimated parameter",
          main=names[i],col="black")
  lines(lpars[,i],col="blue")
  lines(upars[,i],col="blue")
  abline(h=pars.true[i],col="red",lty="dotted")
}
