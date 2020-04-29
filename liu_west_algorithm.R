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
G          = function(arg){c(arg[1],arg[1]/(1+arg[1]^2),cos(1.2*arg[2]))}
pri        = function(x,x0,sV0){dnorm(x,x0,sV0)}
likelihood = function(y,x,sig){dnorm(y,x^2/20,sig)}
rlike      = function(x,sig){rnorm(1,x^2/20,sig)}
post       = function(y,x,x0,V0){pri(x,x0,sV0)*likelihood(y,x)}
quant025   = function(x){quantile(x,0.025)}
quant975   = function(x){quantile(x,0.975)}

# Simulating the data
set.seed(12345)
n         = 200
sig2      = 10
tau2      = 1
sig       = sqrt(sig2)
tau       = sqrt(tau2)
beta      = c(0.5,25,8)
pars.true = c(beta,sig2,tau2)
y         = rep(0,n)
x         = rep(0,n)
x0        = 0.1
zero      = 0.0
g         = G(c(x0,zero))
x[1]      = rnorm(1,sum(g*beta),tau)
y[1]      = rlike(x[1],sig)
for (t in 2:n){
  g        = G(c(x[t-1],t-1))
  x[t] = rnorm(1,sum(g*beta),tau)
  y[t]     = rlike(x[t],sig)
}

pdf(file="nonlinear-data.pdf",width=10,height=7)
par(mfrow=c(2,2))
ts.plot(y,xlab="time",ylab="",main=expression(y[t]))
ts.plot(x,xlab="time",ylab="",main=expression(x[t]))
plot(x,y,xlab=expression(x[t]),ylab=expression(y[t]))
plot(x[1:(n-1)],x[2:n],xlab=expression(x[t-1]),ylab=expression(x[t]))
dev.off()

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

# Liu and West filter
# -------------------
set.seed(8642)
N       = 10000
delta   = 0.99
h2      = 1-((3*delta-1)/(2*delta))^2
a       = sqrt(1-h2)
pars    = cbind(rnorm(N,c0[1],sqrt(C0[1])),rnorm(N,c0[2],sqrt(C0[2])),rnorm(N,c0[3],sqrt(C0[3])),
                log(1/rgamma(N,a0,A0)),log(1/rgamma(N,b0,B0)))
xs  = rnorm(N,m0,sqrt(V0))
ESS     = NULL
xss = NULL
parss   = array(0,c(N,5,n))
ws      = NULL
par(mfrow=c(1,1))
for (t in 1:n){
  m.phi         = apply(phi,2,mean)
  v.par         = var(phi)
  phi.prior     = a*phi+(1-a)*matrix(m.phi,N,5,byrow=T)
  X             = t(apply(cbind(theta,t-1),1,G))
  theta.prior   = apply(X*phi[,1:3],1,sum)
  p             = dnorm(y[t],0.05*theta.prior^2,exp(phi.prior[,4]/2))
  k             = sample(1:N,size=N,replace=T,prob=p)
  phi.post      = phi.prior[k,] + matrix(rnorm(5*N),N,5)%*%chol(h2*vpar)
  X             = t(apply(cbind(theta[k],t-1),1,G))
  theta.post    = rnorm(N,apply(X*phi.post[,1:3],1,sum),exp(phi.post[,5]/2))
  w             = dnorm(y[t],a*theta.post^2,exp(ms1[,4]/2))/likelihood(y[t],mus[k],exp(ms[k,4]/2))
  w      = w/sum(w)
  ind    = sample(1:N,size=N,replace=T,prob=w)
  xs = xt[ind]
  pars   = ms1[ind,]
  xss     = rbind(xss,xs)
  parss[,,t]  = pars
}

mx = apply(xss,1,mean)
lx = apply(xss,1,quant025)
ux = apply(xss,1,quant975)
mpars  = matrix(0,n,5)
lpars  = matrix(0,n,5)
upars  = matrix(0,n,5)
for (i in 1:3){
  mpars[,i] = apply(parss[,i,],2,mean)
  lpars[,i] = apply(parss[,i,],2,quant025)
  upars[,i] = apply(parss[,i,],2,quant975)
}
for (i in 4:5){
  mpars[,i] = apply(exp(parss[,i,]),2,mean)
  lpars[,i] = apply(exp(parss[,i,]),2,quant025)
  upars[,i] = apply(exp(parss[,i,]),2,quant975)
}

par(mfrow=c(2,3))
plot(x,mx,xlab="TRUE",ylab="Posterior mean",main=expression(x[t]),pch=16)
abline(0,1,col=2,lwd=2)
names = c("alpha","beta","gamma","sigma2","tau2")
for (i in 1:5){
  ts.plot(lpars[,i],ylim=range(lpars[,i],upars[,i]),ylab="",main=names[i])
  lines(lpars[,i],col=1)
  lines(mpars[,i],col=1)
  lines(upars[,i],col=1)
  abline(h=pars.true[i],col=2)
}


