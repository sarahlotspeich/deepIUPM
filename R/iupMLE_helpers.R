##########################################
# Simulation study                12/29/15
#
# Author: Pedro Baldoni
#
# Functions:
#   -loglik: loglikelihood function
#   -gloglik: gradient of loglik
#   -sim.data: simulates an experiment
#   -mle: computes the MLE's
#   -hec,hec.jac: equality constraint and gradient of the equality constraint for the profile likelihood confidence interval (not used anymore)
#   -hin,hin.jac: inequality constraint and gradient of the inequality constraint for the profile likelihood confidence interval (not used anymore)
#   -mle.profile: computes the restricted MLE's
#   -CI.profile: computes the profile Confidence Interval
##########################################
#rm(list=ls())
#install.packages("alabama")
#library(alabama)
#setwd("~sarahlotspeich/Dropbox/UNC/Collab-CFAR/Hudgens-IUPM-Paper/")
################
###Simulation###
################
loglik<-function(l,n,M,xM,m1,count,...)
{
  if(sum(l<0)>0){ll = 0}
  else{
    #l: vector of parameters
    #n: # of lineages (rows)
    #M: # of wells originally sequenced (columns)
    #xM: # of positive wells originally sequenced (xM<M)
    #m: # of wells sequenced (m<xM)
    #count: vector of counts (length=n)

    ll = 0
    for(i in 1:n)
    {
      ll = ll + (count[i]*log(1-exp(-l[i])) - l[i]*(M-xM+m1-count[i])) #Observed loglikelihood
    }
    ll = ll + (xM-m1)*log(1-exp(-sum(l))) #Missing loglikelihood
  }
  return(-ll)
}
gloglik<-function(l,n,M,xM,m1,count,...)
{
  gradient = NULL
  for(i in 1:n)
  {
    gradient[i] = count[i]*(exp(-l[i])/(1-exp(-l[i])))-(M-xM+m1-count[i])+(xM-m1)*(exp(-sum(l))/(1-exp(-sum(l))))
  }
  return(-gradient)
}
heq <- function(l,That1,...) {#Equality constraint
  h <- rep(NA, 1)
  h[1] <- sum(l) - That1
  h
}
heq.jac <- function(l,...) {#Equality gradient
  j <- matrix(NA, 1, length(l))
  j[1, ] <- rep(1,length(l))
  j
}
hin <- function(l,...) {#Inequality constraint
  h <- rep(NA, 1)
  for(i in 1:length(l)){h[i] = l[i]}
  h
}
hin.jac <- function(l,...) {#Inequality gradient
  j = diag(length(l))
  j
}
mle.profile<-function(data,That,mle.old,epsilon=10^-6,it=10000)
{
  grid = That - sum(mle.old)
  M = ncol(data)
  n = nrow(data)
  xM = M-sum(colSums(data)==0,na.rm=T)
  m = xM-sum(is.na(colSums(data)))
  count = rowSums(data,na.rm=T)

  U = rbind(diag(n),rep(-1,n),rep(1,n))
  C = c(rep(epsilon,n),-(That+epsilon),+(That-epsilon))

  if(That==sum(mle.old+grid/n)){theta = (mle.old+grid/n)}
  else
  {
    theta = (mle.old+grid/n)+c((That-sum(mle.old+grid/n)),rep(0,(length(mle.old)-1)))
    #if(That==sum((mle.old+grid/n)+(That-sum(mle.old+grid/n))/n)){theta = (mle.old+grid/n)+(That-sum(mle.old+grid/n))/n}
    #else{theta = rep(That/n,n)}
  }
  #theta = mle.old+grid/n

  lambda = constrOptim(theta=theta,f=loglik,grad=gloglik,ui=U,ci=C,control=list(maxit=it),
                       n=n,M=M,xM=xM,m1=m,count=count)$par
  #lambda =  constrOptim.nl(par=theta,fn=loglik,gr=gloglik,heq=heq,heq.jac=heq.jac,hin=hin,hin.jac=hin.jac,
  #                    n=n,M=M,xM=xM,m1=m,count=count,That1=That)$par

  return(lambda)
}
CI.profile<-function(data,alpha,grid)
{
  M = ncol(data)
  n = nrow(data)
  xM = M-sum(colSums(data)==0,na.rm=T)
  m = xM-sum(is.na(colSums(data)))
  count = rowSums(data,na.rm=T)

  cutChisq = qchisq(1-alpha,1)

  lhat = mle(data)$mle
  That = sum(lhat)

  l0 = - loglik(l=lhat,n=n,M=M,xM=xM,m1=m,count=count) #loglik(MLE)

  diff.lower<-0
  diff.upper<-0
  lambda.lower<-lhat
  lambda.upper<-lhat

  That.lower = That
  That.upper = That

  i<-1
  while(diff.lower<cutChisq | diff.upper<cutChisq)
  {
    # cat(rep('# ',28),'\n')
    # cat('T_hat (LB, lower bound): ',That.lower,'.\n')
    # cat('T_hat (UB, upper bound): ',That.upper,'.\n')
    # cat('\n')
    # cat('2*(loglik0 - loglik_LB): ',diff.lower,' (it should be grater than ',cutChisq,').\n')
    # cat('2*(loglik0 - loglik_UB): ',diff.upper,' (it should be grater than ',cutChisq,').\n')
    # cat(rep('# ',28),'\n')

    That.lower.new<-ifelse((That.lower-grid)<0,0,That.lower-grid)
    That.upper.new<-That.upper+grid

    if(diff.lower>cutChisq){
      diff.lower<-diff.lower
    } else{
      lambda.lower.new = mle.profile(data=data,That=That.lower.new,mle.old=lambda.lower)
      diff.lower<-2*(l0+loglik(l=lambda.lower.new,n=n,M=M,xM=xM,m1=m,count=count))
      lambda.lower<-lambda.lower.new
    }

    if(diff.upper>cutChisq){
      diff.upper<-diff.upper
    } else{
      lambda.upper.new<-mle.profile(data=data,That=That.upper.new,mle.old=lambda.upper)
      diff.upper<-2*(l0+loglik(l=lambda.upper.new,n=n,M=M,xM=xM,m1=m,count=count))
      lambda.upper<-lambda.upper.new
    }
    i<-i+1
    That.lower = That.lower.new
    That.upper = That.upper.new
  }
  CI<-c(sum(lambda.lower),sum(lambda.upper))
  return(CI)
}



###Example
#n=6
#M=32
#l=rep(0.5,n)
#p=0.1
#
#data = sim.data(n,M,l,p)
#
#mle(data)$mle
#That = sum(mle(data)$mle)
#mle(data)$cov #Inverse of -d2l/dtheta2
#mle.profile(data,That+0.01,mle(data)$mle) #MLE subject restricted to sum(lambda)=1.5
#CI.profile(data=data,alpha=0.05,grid=0.01)
