Posterior Simulation: Sequential Monte Carlo
================

## Nonlinear time series model

The following nonlinear time series model is analyzed.

![](README_files/figure-gfm/equations-1.png)<!-- -->![](README_files/figure-gfm/equations-2.png)<!-- -->

where and \(w_t\thicksim N(0,w)\) are mutually independent Gaussian
random variables.

``` r
# Function to simulate observations for model in example 1
simulate.example1 = function(Tbig,a,b,c,d,omega,v,w)
{
  y = rep(NA,Tbig)
  theta = rep(NA,Tbig)
  theta[1] = rnorm(1,mean=0,sd=sqrt(w))
  y[1] = a*theta[1]^2 + rnorm(1,mean=0,sd=sqrt(v))
  for(t in 2:Tbig)
  {
    theta[t] = b*theta[t-1] + c*theta[t-1]/(1+theta[t-1]^2) + d*cos(omega*t) + rnorm(1,mean=0,sd=sqrt(w))
    y[t] = a*theta[t]^2 + rnorm(1,mean=0,sd=sqrt(v))
  }
  list("y"=y,"theta"=theta)
}
```

The function is used to generate \(y_t\) and \(\theta_t\) for
\(t=1:200\) and with particular values of the parameters
\((a=\frac{1}{20},b=\frac{1}{2},c=25,d=8,\omega=1.2,v=10\textrm{ and }w=1)\).

``` r
a = 1/20
b = 1/2
c = 25
d = 8
omega = 1.2
v = 10
w = 1
Tbig = 200

set.seed(12345) # Set seed for reproducability
sim.data = simulate.example1(Tbig,a,b,c,d,omega,v,w)
sim.y = sim.data$y
sim.theta = sim.data$theta

# Plot the time series
par(mfrow=c(2,2))
plot(sim.y,type="l",xlab="time,t",ylab="Observations,y")
plot(sim.theta,type="l",xlab="time,t",ylab=TeX("State,$\\theta$"))

plot(sim.theta,sim.y,type="p",pch=19,xlab=TeX("State,$\\theta$"),ylab="Observations,y")
plot(sim.theta[1:199],sim.theta[2:200],type="p",pch=19,xlab=TeX("State,$\\theta_{t-1}$"),
     ylab=TeX("State,$\\theta_{t}$"))
```

![](README_files/figure-gfm/simulate-1.png)<!-- -->

The first plot shows the variation of observations \(y_t\) with state
\(\theta_t\) and the second one plots \(`\theta_t`\) versus
\(`\theta_{t-1}`\).

The first goal is to find the probability density functions
\(`p(\theta_t^{(m)}|\theta_{(t-1)}^{(m)})`\) and
\(`p(y_t|\theta_{(t)}^{(m)})`\). From the state evolution equation, with
evolution error \(`w_t`\) following a Gaussian distribution
\(`N(0,w)`\), we can evaluate the density functions to be a Gaussian
distribution \(`N\big(a\theta_{t}^2,v\big)`\) and
\(`N\bigg(b\theta^{(m)}_{t-1}+c\dfrac{\theta^{(m)}_{t-1}}{1+\theta_{t-1}^{(m)2}}+d\cos(\omega t),w\bigg)`\).
The importance density
\(`g_t(\theta_t^{(m)}|\theta_{0:(t-1)}^{(m)},y_{1:t})`\) is considered
to be Gaussian distribution \(`N(\theta_{t-1},w)`\).

``` r
#####################################################################
# Computes the expection of current state given previous state:     #
# E[theta_t|theta_{t-1}]=b*theta+c*theta/(1+theta^2)+d*cos(omega*t) #
mu.theta = function(theta,t) b*theta+c*theta/(1+theta^2)+d*cos(omega*t)

#####################################################################
# Computes the expection of current observation given present state:#
# E[y_t|theta_t]=a*theta^2                                          #
mu.y = function(theta) a*theta^2                                    

#####################################################################
# Returns a sample from the importance distribution and the weight  #
# to scale the target probability distribution to the importance    #
imp.sampling <- function(M,y,t,theta.previous,weight.prev)
{
  theta.current <- rep(NA,M)
  weight.current <- rep(NA,M)
  for(m in 1:M)
  {
    theta.current[m] <- rnorm(1,mean=mu.theta(theta.previous[m],t),sd=sqrt(w))
    num1 <- dnorm(y[t],mu.y(theta.current[m]),sqrt(v))
    num2 <- dnorm(theta.current[m],mu.theta(theta.previous[m],t),sqrt(w))
    den <- dnorm(theta.current[m],mean=theta.previous[m],sd=sqrt(w))
    weight.current[m] <- weight.prev[m]*num1*num2/den
  }
  weight.current <- weight.current/sum(weight.current)
  return(list("theta"=theta.current,"weight"=weight.current))
}

#####################################################################
# Performs filtering operation with importance sampling to estimate
# the states for times t=1 to t=T.
smc.sis <- function(Tbig,M,y)
{
  # Initialize the arrays of estimates
  theta <- array(NA,dim=c(M,Tbig))
  est.theta <- array(NA,dim=c(M,Tbig))
  weight <- array(NA,dim=c(M,Tbig))
  
  # Initial state and weight
  theta[,1] <- rnorm(M,mean=0,sd=sqrt(w))
  weight[,1] <- rep(1/M,M)
  # Find estimated state
  est.theta[,1] <- sample(x = theta[,1], M, 
                          replace = T, prob = weight[,1])
  
  for(t in 2:Tbig)
  {
    # Importance sampling step
    samples <- imp.sampling(M,y,t,est.theta[,t-1],weight[,t-1])
    theta[,t] <- samples$theta
    weight[,t] <- samples$weight
    # Find estimated state through resampling
    est.theta[,t] <- sample(x = theta[,t], M, 
                          replace = T, prob = weight[,t])
  }
  theta.mean <- apply(est.theta,2,mean)
  theta.var <- apply(est.theta,2,var)
  return(list("mu"=theta.mean,"sigma2"=theta.var))
}
```

``` r
sis.est.theta <- smc.sis(200,500,sim.y)
theta.mean <- sis.est.theta$mu
theta.var <- sis.est.theta$sigma2

ll = min(c(min(theta.mean-2*sqrt(theta.var)),min(sim.theta),sim.y))
ul = max(c(max(theta.mean+2*sqrt(theta.var)),max(sim.theta),sim.y))

plot(theta.mean,type="l",xlab="time",ylab="State",ylim=c(ll,ul),
     lwd=2)
lines(sim.theta,col="red")
lines(theta.mean+1.96*sqrt(theta.var),lty=4,lwd=0.5,col="blue")
lines(theta.mean-1.96*sqrt(theta.var),lty=4,lwd=0.5,col="blue")
```

![](README_files/figure-gfm/estimate-1.png)<!-- -->

``` r
#points(sim.y,lwd=2)
```

We now analyze the auxillary particle filter approach for posterior
simulation using sequential Monte Carlo. The first

``` r
##########################################################################
# Samples M auxillary variables from the set {1,2,...,M}. The probablity #
# of each index is proportional to w_k*p(y|mu_k) where w_k is the weight #
# corresponding to the kth sample of previous iteration and mu_k is the  #
# mean of the distribution of theta_t|theta_{t-1}^{(k)}.

aux.var <- function(M,y,t,theta.previous,weight.previous)
{
  # Get probablity of auxillary particles
  prob <- rep(NA,M)
  for(k in 1:M)
  {
    mu.k <- mu.theta(theta.previous[k],t)
    prob[k] <- dnorm(y[t],mean=mu.y(mu.k),sd=sqrt(v))*weight.previous[k]
  }
  prob <- prob/sum(prob)
  
  # Sample auxillary variable
  return(sample(M,size=M,replace=T,prob=prob))
}

##########################################################################
# Performs importance sampling with auxillary particle sampling

imp.sampling.aux <- function(M,y,t,theta.previous,weight.prev,aux.var)
{
  theta.current <- rep(NA,M)
  weight.current <- rep(NA,M)
  for(m in 1:M)
  {
    k <- aux.var[m]
    theta.current[m] <- rnorm(1,mean=mu.theta(theta.previous[k],t),
                              sd=sqrt(w))
    mu.aux <- mu.theta(theta.previous[k],t)
    
    num <- dnorm(y[t],mean=mu.y(theta.current[m]),sd=sqrt(v))
    den <- dnorm(y[t],mean=mu.y(mu.aux),sd=sqrt(v))
    weight.current[m] <- num/den
  }
  weight.current <- weight.current/sum(weight.current)
  return(list("theta"=theta.current,"weight"=weight.current))
}

smc.apf <- function(Tbig,M,y)
{
  # Initialize the arrays of estimates
  theta <- array(NA,dim=c(M,Tbig))
  weight <- array(NA,dim=c(M,Tbig))
  est.theta <- array(NA,dim=c(M,Tbig))
  
  # Initial time instant
  theta[,1] <- rnorm(M,mean=0,sd=sqrt(w))
  weight[,1] <- rep(1/M,M)
  est.theta[,1] <- sample(x = theta[,1], M, 
                          replace = T, prob = weight[,1])
  
  # Iterate over the time instants
  for (t in 2:Tbig)
  {
    # Find probabilities for auxilliary variables
    aux.k <- aux.var(M,y,t,est.theta[,t-1],weight)
    # Importance sampling step
    samples <- imp.sampling.aux(M,y,t,est.theta[,t-1],weight[,t-1],aux.k)
    theta[,t] <- samples$theta
    weight[,t] <- samples$weight
    
    # Find estimated state through resampling
    est.theta[,t] <- sample(x = theta[,t], M, 
                          replace = T, prob = weight[,t])
  }
  theta.mean <- apply(est.theta,2,mean)
  theta.var <- apply(est.theta,2,var)
  return(list("mu"=theta.mean,"sigma2"=theta.var))
}
```

``` r
apf.est.theta <- smc.apf(200,500,sim.y)
apf.theta.mean <- apf.est.theta$mu
apf.theta.var <- apf.est.theta$sigma2

ll = min(c(min(apf.theta.mean-2*sqrt(apf.theta.var)),min(sim.theta),sim.y))
ul = max(c(max(apf.theta.mean+2*sqrt(apf.theta.var)),max(sim.theta),sim.y))

plot(apf.theta.mean,type="l",xlab="time",ylab="State",ylim=c(ll,ul),
     lwd=2)
lines(sim.theta,col="red")
lines(apf.theta.mean+1.96*sqrt(apf.theta.var),lty=4,lwd=0.5,col="blue")
lines(apf.theta.mean-1.96*sqrt(apf.theta.var),lty=4,lwd=0.5,col="blue")
```

![](README_files/figure-gfm/estimate_aux-1.png)<!-- -->

``` r
#points(sim.y,lwd=2)
```
