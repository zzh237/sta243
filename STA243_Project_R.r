
## function to run 10 times for a given theta.opt and sigma and given algorithm
r1 = function(c,sigma,datafunc,algofunc){
  data = datafunc(c,sigma)
  theta.opt = data$theta
  theta0=c()
  ## based on Fig3, the theta.opt_i - sqrt(c/10)<=theta0_i<=theta.opt_i + sqrt(c/10)
  dif = runif(1,c*sqrt(1/10)-0.5,c*sqrt(1/10)-0.5)
  for(i in 1:length(theta.opt)){
    p = sample(c(1,-1),1)
    value = theta.opt[i] + p*dif
    theta0 = append(theta0, value)
  } 
  norm(theta.opt,'2')
  norm(theta0,'2')
  norm((theta.opt-theta0),'2')
  # call the function
  mod=algofunc(data,theta0,sigma)
  theta.t = mod$theta.t
  theta.hat = mod$theta.hat
  end = length(theta.t)
  
  print(end)
  err.opt = c()
  err.stat = c()
  for(i in 1:end){
    theta = theta.t[[i]]
    e1=log(norm(theta.hat-theta,'2'))
    e2=log(norm(theta.opt-theta,'2'))
    #print(e1)
    #print(e2)
    err.opt=c(err.opt, e1)
    err.stat=c(err.stat, e2)
  }
  return(list(err.opt=err.opt,err.stat=err.stat))
}



## Linear Regression with Missing Covariates

LRMCdata<-function(c,sigma){
  n = 1000;k=2;d=10;sigma2 =sigma^2;p = 0.2
  y = matrix(0,nrow=n, ncol=1)
  gauss = rnorm(d)
  length = norm(gauss,'2')
  c=c
  theta.opt=c*gauss/length
  #print(norm(theta.opt,'2'))
  theta = matrix(1:(k*d),nrow=k,byrow = T)
  theta[1,] = theta.opt
  theta[2,] = -theta.opt
  pi = c(1,1)/k
  Z = sample(1:k,n,prob=pi,replace = T)
  X = matrix(0,nrow=n,ncol=d)
  y = matrix(0,nrow=n, ncol=1)
  for(i in 1:n){
    xi = rnorm(d,mean=0,sd=1)
    mask = sample(c(NA,1), prob=c(p, 1-p), replace=TRUE, size=d)
    X[i,] = xi*mask
    vi = rnorm(1, mean=0,sd=1)
    y[i,] = sum(X[i,]*theta[Z[i],], na.rm = TRUE) +vi
  }
  return(list(X=X,y=y,theta=theta.opt))
}

## EM for LRMC
EMLRMC<-function(data,theta0,sigma0,alpha=0.01,eta=0.0001){
  X=data$X
  y=data$y
  sigma20=sigma0^2
  n = nrow(X); d = ncol(X)
  theta.t = list()
  repeat({
    first.term=0
    second.term=0
    for(i in 1:n){
      xi=X[i,]
      yi=y[i,]
      NAindex <-which(is.na(xi))
      numNA = length(NAindex)
      
      if(numNA>0)
      {
        NonNAindex <- which(!is.na(xi))
        xi.obs=xi[NonNAindex]
        xi.mis=xi[NAindex]
        theta0.obs=theta0[NonNAindex]
        theta0.mis=theta0[NAindex]
        z.obs=c(xi.obs,yi)
        U.theta=cbind(-theta0.mis%*%t(theta0.obs),theta0.mis)*1/(norm(theta0.mis,type="2")^2+sigma20)
        u=c(U.theta%*%z.obs,xi.obs)
        p1=matrix(1, numNA, numNA)
        p2=U.theta%*%z.obs%*%t(xi.obs)
        p3=xi.obs%*%t(z.obs)%*%t(U.theta)
        p4=xi.obs%*%t(xi.obs)
        #conditional second moment matrix
        csmm=rbind(cbind(p1,p2),cbind(p3,p4))
        second.term=second.term+yi*u
        first.term=first.term+csmm
      }
      else
      { 
        second.term=second.term+yi*xi #d*1 vector
        first.term=first.term+xi%*%t(xi) #d*d matrix
      }
      
    }
    theta=solve(first.term)%*%second.term
    err = norm(theta - theta0,"2")
    theta.t[[length(theta.t)+1]] = theta
    theta0=theta
    #print(err)
    if(err<eta)
      break
  })
  return(list(theta.hat=theta,theta.t=theta.t))
}



## First-Order-EM for LRMC

FEMLRMC<-function(data,theta0,sigma0,alpha=0.01,eta=0.0001){
  X=data$X
  y=data$y
  sigma20=sigma0^2
  n = nrow(X); d = ncol(X)
  theta.t = list()
  repeat({
    gradient=0
    for(i in 1:n){
      xi=X[i,]
      yi=y[i,]
      NAindex <-which(is.na(xi))
      numNA = length(NAindex)
      if(numNA>0)
      {
        NonNAindex <- which(!is.na(xi))
        xi.obs=xi[NonNAindex]
        xi.mis=xi[NAindex]
        theta0.obs=theta0[NonNAindex]
        theta0.mis=theta0[NAindex]
        z.obs=c(xi.obs,yi)
        U.theta=cbind(-theta0.mis%*%t(theta0.obs),theta0.mis)*1/(norm(theta0.mis,type="2")^2+sigma20)
        u=c(U.theta%*%z.obs,xi.obs)
        p1=matrix(1, numNA, numNA)
        p2=U.theta%*%z.obs%*%t(xi.obs)
        p3=xi.obs%*%t(z.obs)%*%t(U.theta)
        p4=xi.obs%*%t(xi.obs)
        #conditional second moment matrix
        csmm=rbind(cbind(p1,p2),cbind(p3,p4))
        gradient=gradient+yi*u-csmm%*%theta0
      }
      else
      {
        gradient=gradient+yi*xi-xi%*%t(xi)%*%theta0
      }
      
    }
    update=1/n*gradient
    theta=theta0+alpha*update
    err = norm(theta - theta0,"2")
    theta.t[[length(theta.t)+1]] = theta
    theta0=theta
    #print(err)
    if(err<eta)
      break
  })
  return(list(theta.hat=theta,theta.t=theta.t))
}


## Fig8(a)
rs = lapply(rep(2,10), function(c){
  sigma = 1
  out=r1(c,sigma,LRMCdata,EMLRMC)
})

end = 0
for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  end_ = length(err.opt)
  end = max(end,end_)
}

end=10
x = seq(1:end-1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-7,1),
     main='EM, Missing Data Regression')

for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  err.stat=rs[[i]]$err.stat
  end = length(err.opt)-1
  end=10
  x = seq(1:end)
  points(x, err.opt[1:end], col="blue", pch="o")
  lines(x, err.opt[1:end], col="blue", lty=1)
  points(x, err.stat[1:end], col="red", pch="*")
  lines(x, err.stat[1:end], col="red", lty=2)
}
legend(1,-5,legend=c("Opt. error","Stat. error"), col=c("blue","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)


## Fig8(b)
rs = lapply(rep(2,10), function(c){
  sigma = 1
  out=r1(c,sigma,LRMCdata,FEMLRMC)
})

end = 0
for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  end_ = length(err.opt)
  end = max(end,end_)
}

end=500
x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-7.5,1),
     main='First-order EM, Missing Data Regression',cex=0.3,lwd=0.3)

for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  err.stat=rs[[i]]$err.stat
  end = length(err.opt)-1
  end =500
  x = seq(1:end)
  points(x, err.opt[1:end], col="blue", pch="o",cex=0.3,lwd=0.3)
  lines(x, err.opt[1:end], col="blue", lty=1, cex=0.3,lwd=0.3)
  points(x, err.stat[1:end], col="red", pch="*",cex=0.3,lwd=0.3)
  lines(x, err.stat[1:end], col="red", lty=2,cex=0.3,lwd=0.3)
}
legend(1,-5,legend=c("Opt. error","Stat. error"), col=c("blue","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)










## fig6 suplement for missing data covariates
end = 0
rss = list()
snrs = c(0.5,0.75,1,1.8,2.5)
repnum=10
for(i in snrs){
  rs = lapply(rep(i,repnum), function(i){
    sigma = 1
    out=r1(i,sigma,LRMCdata,EMLRMC)
  })
  for(j in 1:repnum){
    err.opt=rs[[j]]$err.opt
    end_ = length(err.opt)
    end = max(end,end_)
  } 
  rss[[length(rss)+1]]=rs
}


#end=100
x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-10,1),
     main='EM, Mixture of Missing Data Regression',cex=0.3,lwd=0.3)

colors=c('red','blue','green','black','cyan')
for(j in 1:length(snrs)){
  rs=rss[[j]]
  color=colors[j]
  for(i in 1:repnum){
    err.opt=rs[[i]]$err.opt
 #   end=100
    end = length(err.opt)-1
    x = seq(1:end)
    points(x, err.opt[1:end], col=color, pch="o",cex=1,lwd=1)
    lines(x, err.opt[1:end], col=color, lty=1, cex=1,lwd=1)
  }
}
legend(20,-2,legend=snrs, col=colors,
       pch=rep("o",5),lty=rep(1,5), ncol=1)



## fig6 suplement for missing data covariates
end = 0
rss = list()
snrs = c(0.5,0.75,1,1.8,2.5)
repnum=10
for(i in snrs){
  rs = lapply(rep(i,repnum), function(i){
    sigma = 1
    out=r1(i,sigma,LRMCdata,FEMLRMC)
  })
  for(j in 1:repnum){
    err.opt=rs[[j]]$err.opt
    end_ = length(err.opt)
    end = max(end,end_)
  } 
  rss[[length(rss)+1]]=rs
}

#end=100
x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-10,1),
     main='First-Order-EM, Mixture of Missing Data Regression',cex=0.3,lwd=0.3)

colors=c('red','blue','green','black','cyan')
for(j in 1:length(snrs)){
  rs=rss[[j]]
  color=colors[j]
  for(i in 1:repnum){
    err.opt=rs[[i]]$err.opt
    #   end=100
    end = length(err.opt)-1
    x = seq(1:end)
    points(x, err.opt[1:end], col=color, pch="o",cex=1,lwd=1)
    lines(x, err.opt[1:end], col=color, lty=1, cex=1,lwd=1)
  }
}
legend(90,-2,legend=snrs, col=colors,
       pch=rep("o",5),lty=rep(1,5), ncol=1)

## Figure3 
sigma=1
d=10
theta0 = matrix(runif(d),nrow=d,ncol=1)

rs3 = sapply(c(1,1.5,2,2.5,3,3.5,4), function(a){
  sigma=1
  data=LRMCdata(a,sigma)
  theta.opt=data$theta
  #print(theta.opt)
  mod=EMLRMC(data,theta0,sigma)
  theta=mod$theta
  r=norm(theta.opt-theta,'2')
  #print(norm(theta0-theta,'2'))
  print(r)
  r
})




## GMM Data Generation
GMMdata<-function(c,sigma){
  # c decide the norm of the theta.opt
  n = 1000;k=2;d=10;sigma2 = sigma^2;
  gauss = rnorm(d)
  length = norm(gauss,'2')
  c=c
  ## construct the theta star that length is c
  theta.opt=c*gauss/length
  
  theta = matrix(1:(k*d),nrow=k,byrow = T)
  theta[1,] = theta.opt
  theta[2,] = -theta.opt
  pi = c(1,1)/k
  Z = sample(1:k,n,prob=pi,replace = T)
  X = matrix(0,nrow=n,ncol=d)
  for(i in 1:n){
    X[i,] = rnorm(d,mean=theta[Z[i],],sd=sqrt(sigma2))
  }
  return(list(X=X,theta=theta.opt))
}


## EM-Algorithm for GMM
#X: the data matrix, n X d. n is the number of observations, d is the dimension of features.
#eta: the tolerance of difference of the theta values.  
EMGMM<-function(data,theta0,sigma0,alpha=0.01,eta=0.00001){
  X=data$X
  sigma20=sigma0^2
  n = nrow(X); d = ncol(X)
  theta.t = list()
  repeat({
    first.term=0
    second.term=0
    for(i in 1:n){
      xi=X[i,]
      w.n=exp(-norm(theta0-xi,"2")^2/(2*sigma20))
      w.p=exp(-norm(theta0+xi,"2")^2/(2*sigma20))
      w=w.n*(w.n+w.p)^(-1)
      first.term=first.term+w*xi
      second.term=second.term+xi
    }
    theta=2/n*first.term-1/n*second.term
    err = norm(theta - theta0,"2")
    theta.t[[length(theta.t)+1]] = theta
    theta0=theta
    #print(err)
    if(err<eta)
      break
  })
  return(list(theta.hat=theta,theta.t=theta.t))
}


## First-Order EM-Algorithm for GMM
FEMGMM<-function(data,theta0,sigma0,alpha=0.01,eta=0.0001){
  X=data$X
  sigma20=sigma0^2
  n = nrow(X); d = ncol(X)
  theta.t = list()
  repeat({
    gradient=0
    for(i in 1:n){
      xi=X[i,]
      w.n=exp(-norm(theta0-xi,"2")^2/(2*sigma20))
      w.p=exp(-norm(theta0+xi,"2")^2/(2*sigma20))
      w=w.n*(w.n+w.p)^(-1)
      gradient=gradient+(2*w-1)*xi
    }
    update=1/n*gradient-theta0
    theta=theta0+alpha*update
    err = norm(theta - theta0,"2")
    theta.t[[length(theta.t)+1]] = theta
    theta0=theta
    #print(err)
    if(err<eta)
      break
  })
  return(list(theta.hat=theta,theta.t=theta.t))
}


## Fig5(a)
rs = lapply(rep(2,10), function(c){
  sigma = 1
  out=r1(c,sigma,GMMdata,EMGMM)
})

end = 0
for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  end_ = length(err.opt)
  end = max(end,end_)
}

x = seq(1:end-1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-8,1),
     main='EM, Mixture of Gaussian',cex=1,lwd=1)

for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  err.stat=rs[[i]]$err.stat
  end = length(err.opt)-1
  x = seq(1:end)
  points(x, err.opt[1:end], col="blue", pch="o",cex=1,lwd=1)
  lines(x, err.opt[1:end], col="blue", lty=1, cex=1,lwd=1)
  points(x, err.stat[1:end], col="red", pch="*",cex=1,lwd=1)
  lines(x, err.stat[1:end], col="red", lty=2,cex=1,lwd=1)
}
legend(1,-6,legend=c("Opt. error","Stat. error"), col=c("blue","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)


## Fig 5(b)
rs = lapply(rep(2,10), function(c){
  sigma = 1
  out=r1(c,sigma,GMMdata,FEMGMM)
})

end = 0
for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  end_ = length(err.opt)
  end = max(end,end_)
}

end=500
x = seq(1:end-1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-6,1),
     main='First-order EM, Mixture of Gaussian',cex=0.3,lwd=0.3)

for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  err.stat=rs[[i]]$err.stat
  end = length(err.opt)-1
  end=500
  x = seq(1:end)
  points(x, err.opt[1:end], col="blue", pch="o",cex=0.3,lwd=0.3)
  lines(x, err.opt[1:end], col="blue", lty=1, cex=0.3,lwd=0.3)
  points(x, err.stat[1:end], col="red", pch="*",cex=0.3,lwd=0.3)
  lines(x, err.stat[1:end], col="red", lty=2,cex=0.3,lwd=0.3)
}
legend(1,-4,legend=c("Opt. error","Stat. error"), col=c("blue","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)


## Fig6
end = 0
rss = list()
snrs = c(0.5,0.75,1,1.8,2.5)
repnum=10
for(i in snrs){
  rs = lapply(rep(i,repnum), function(i){
    sigma = 1
    out=r1(i,sigma,GMMdata,EMGMM)
  })
  for(j in 1:repnum){
    err.opt=rs[[j]]$err.opt
    end_ = length(err.opt)
    end = max(end,end_)
  } 
  rss[[length(rss)+1]]=rs
}


x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-10,1),
     main='EM, Mixture of Gaussian',cex=0.3,lwd=0.3)

colors=c('red','blue','green','black','cyan')
for(j in 1:length(snrs)){
  rs=rss[[j]]
  color=colors[j]
  for(i in 1:repnum){
    err.opt=rs[[i]]$err.opt
    end = length(err.opt)-1
    x = seq(1:end)
    points(x, err.opt[1:end], col=color, pch="o",cex=1,lwd=1)
    lines(x, err.opt[1:end], col=color, lty=1, cex=1,lwd=1)
  }
}
legend(70,-2,legend=snrs, col=colors,
       pch=rep("o",5),lty=rep(1,5), ncol=1)





end = 0
rss = list()
snrs = c(0.5,0.75,1,1.8,2.5)
repnum=10
for(i in snrs){
  rs = lapply(rep(i,repnum), function(i){
    sigma = 1
    out=r1(i,sigma,GMMdata,FEMGMM)
  })
  for(j in 1:repnum){
    err.opt=rs[[j]]$err.opt
    end_ = length(err.opt)
    end = max(end,end_)
  } 
  rss[[length(rss)+1]]=rs
}


x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-10,1),
     main='First-Order-EM, Mixture of Gaussian',cex=0.3,lwd=0.3)

colors=c('red','blue','green','black','cyan')
for(j in 1:length(snrs)){
  rs=rss[[j]]
  color=colors[j]
  for(i in 1:repnum){
    err.opt=rs[[i]]$err.opt
    end = length(err.opt)-1
    x = seq(1:end)
    points(x, err.opt[1:end], col=color, pch="o",cex=1,lwd=1)
    lines(x, err.opt[1:end], col=color, lty=1, cex=1,lwd=1)
  }
}
legend(70,-2,legend=snrs, col=colors,
       pch=rep("o",5),lty=rep(1,5), ncol=1)






#Figure3
sigma2=1
d=10
theta0=runif(d)
rs = sapply(c(1,1.5,2,2.5,3,3.5,4), function(a){
  X=GMMdata(a)
  theta.opt=X$theta
  #print(theta.opt)
  mod=EMGMM(X,theta0,sigma2)
  theta.hat=mod$theta.hat
  r=norm(theta.opt-theta.hat,'2')
  r
  print(r)
})

print(rs)





## Mixture of regression data

MRdata<-function(c,sigma){
  sigma2 = sigma^2
  n = 1000;k=2;d=10;sigma2 = 1;
  gauss = rnorm(d)
  length = norm(gauss,'2')
  c=c
  theta.opt=c*gauss/length
  theta = matrix(1:(k*d),nrow=k,byrow = T)
  theta[1,] = theta.opt
  theta[2,] = -theta.opt
  pi = c(1,1)/k
  Z = sample(1:k,n,prob=pi,replace = T)
  X = matrix(0,nrow=n,ncol=d)
  y = matrix(0,nrow=n, ncol=1)
  for(i in 1:n){
    xi = rnorm(d,mean=0,sd=1)
    vi = rnorm(1, mean=0,sd=1)
    X[i,] = xi
    y[i,] = xi%*%theta[Z[i],]+vi
  }
  return(list(X=X,y=y,theta=theta.opt))
}


##  EM for Mixture of Regressions
EMMR<-function(data,theta0,sigma0,alpha=0.01,eta=0.0001){
  X=data$X
  y=data$y
  sigma20=sigma0^2
  n = nrow(X); d = ncol(X)
  theta.t=list()
  repeat({
    first.term=0
    second.term=0
    for(i in 1:n){
      xi=X[i,]
      yi=y[i,]
      w.n=exp(-(yi-xi%*%theta0)^2/(2*sigma20))
      w.p=exp(-(yi+xi%*%theta0)^2/(2*sigma20))
      w=w.n*(w.n+w.p)^(-1)
      first.term=first.term+xi%*%t(xi)
      second.term=second.term+(2*w-1)*yi*xi
    }
    theta=solve(first.term)%*%second.term
    err = norm(theta - theta0,"2")
    theta.t[[length(theta.t)+1]] = theta
    theta0=theta
    #print(err)
    if(err<eta)
      break
  })
  return(list(theta.hat=theta,theta.t=theta.t))
}


## first-order EM for Mixture of Regressions
FEMMR<-function(data,theta0,sigma0,alpha=0.01,eta=0.0001){
  X=data$X
  y=data$y
  sigma20=sigma0^2
  n = nrow(X); d = ncol(X)
  theta.t=list()
  repeat({
    gradient=0
    for(i in 1:n){
      xi=X[i,]
      yi=y[i,]
      w.n=exp(-(yi-xi%*%theta0)^2/(2*sigma20))
      w.p=exp(-(yi+xi%*%theta0)^2/(2*sigma20))
      w=w.n/(w.n+w.p)
      gradient=gradient+(2*w-1)*yi*xi-xi%*%t(xi)%*%theta0
    }
    update=1/n*gradient
    theta=theta0+alpha*update
    err = norm(theta - theta0,"2")
    theta.t[[length(theta.t)+1]] = theta
    theta0=theta
    #print(err)
    if(err<eta)
      break
  })
  return(list(theta.hat=theta,theta.t=theta.t))
}



## Fig7(a)
## based on the theorem, the theta.opt_i - sqrt(c/10)<=theta0_i<=theta.opt_i + sqrt(c/10)
rs = lapply(rep(2,10), function(c){
  sigma = 1
  out=r1(c,sigma,MRdata,EMMR)
})

end = 0
for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  end_ = length(err.opt)
  end = max(end,end_)
}


x = seq(1:end-1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-8,1),
     main='EM, Mixture of Regression',cex=1,lwd=1)

for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  err.stat=rs[[i]]$err.stat
  end = length(err.opt)-1
  x = seq(1:end)
  points(x, err.opt[1:end], col="blue", pch="o",cex=1,lwd=1)
  lines(x, err.opt[1:end], col="blue", lty=1, cex=1,lwd=1)
  points(x, err.stat[1:end], col="red", pch="*",cex=1,lwd=1)
  lines(x, err.stat[1:end], col="red", lty=2,cex=1,lwd=1)
}
legend(1,-6,legend=c("Opt. error","Stat. error"), col=c("blue","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)


## Fig 7(b)
rs = lapply(rep(2,10), function(c){
  sigma = 1
  out=r1(c,sigma,MRdata,FEMMR)
})

end = 0
for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  end_ = length(err.opt)
  end = max(end,end_)
}

end=400
x = seq(1:end-1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-8,1),
     main='First-order EM, Mixture of Regression',cex=0.3,lwd=0.3)

for(i in 1:10){
  err.opt=rs[[i]]$err.opt
  err.stat=rs[[i]]$err.stat
  end = length(err.opt)-1
  end=400
  x = seq(1:end)
  points(x, err.opt[1:end], col="blue", pch="o",cex=0.3,lwd=0.3)
  lines(x, err.opt[1:end], col="blue", lty=1, cex=0.3,lwd=0.3)
  points(x, err.stat[1:end], col="red", pch="*",cex=0.3,lwd=0.3)
  lines(x, err.stat[1:end], col="red", lty=2,cex=0.3,lwd=0.3)
}
legend(1,-6,legend=c("Opt. error","Stat. error"), col=c("blue","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)




## fig6 suplement for regression
end = 0
rss = list()
snrs = c(0.5,0.75,1,1.8,2.5)
repnum=10
for(i in snrs){
  rs = lapply(rep(i,repnum), function(i){
    sigma = 1
    out=r1(i,sigma,MRdata,EMMR)
  })
  for(j in 1:repnum){
    err.opt=rs[[j]]$err.opt
    end_ = length(err.opt)
    end = max(end,end_)
  } 
  rss[[length(rss)+1]]=rs
}

end=2118
x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-10,1),
     main='First-Order-EM, Mixture of Regression',cex=0.3,lwd=0.3)

colors=c('red','blue','green','black','cyan')
for(j in 1:length(snrs)){
  rs=rss[[j]]
  color=colors[j]
  for(i in 1:repnum){
    err.opt=rs[[i]]$err.opt
    end=100
    end = length(err.opt)-1
    x = seq(1:end)
    points(x, err.opt[1:end], col=color, pch="o",cex=1,lwd=1)
    lines(x, err.opt[1:end], col=color, lty=1, cex=1,lwd=1)
  }
}


legend(0,-6,legend=snrs, col=colors,
       pch=rep("o",5),lty=rep(1,5), ncol=1)


### fig6 first order em
end = 0
rss = list()
snrs = c(0.5,0.75,1,1.8,2.5)
repnum=10
for(i in snrs){
  rs = lapply(rep(i,repnum), function(i){
    sigma = 1
    out=r1(i,sigma,MRdata,FEMMR)
  })
  for(j in 1:repnum){
    err.opt=rs[[j]]$err.opt
    end_ = length(err.opt)
    end = max(end,end_)
  } 
  rss[[length(rss)+1]]=rs
}

end=100
x = seq(1:end+1)
plot(x,xlab="Iteration #", ylab="log error",ylim=c(-10,1),
     main='First-Order-EM, Mixture of Regression',cex=0.3,lwd=0.3)

colors=c('red','blue','green','black','cyan')
for(j in 1:length(snrs)){
  rs=rss[[j]]
  color=colors[j]
  for(i in 1:repnum){
    err.opt=rs[[i]]$err.opt
    end=100
    end = length(err.opt)-1
    x = seq(1:end)
    points(x, err.opt[1:end], col=color, pch="o",cex=1,lwd=1)
    lines(x, err.opt[1:end], col=color, lty=1, cex=1,lwd=1)
  }
}
legend(90,-2,legend=snrs, col=colors,
       pch=rep("o",5),lty=rep(1,5), ncol=1)




#Fig3

res=list()
delta=list()
delta[[1]]=c(1,1.5,2,2.5,3,3.5,4)-1
delta[[2]]=c(1,1.5,2,2.5,3,3.5,4)+0.5
theta_norm = c(1,1.5,2,2.5,3,3.5,4)
for(j in 1:length(delta)){
  delta_j = delta[[j]]
  rs = 
    mapply(function(c,e)
    {
    
    sigma=1
    
    data=MRdata(c,sigma)
    theta.opt=data$theta
    print('theta norm')
    print(norm(theta.opt,'2'))

    dif = e*sqrt(1/10)
    for(i in 1:length(theta.opt)){
      p = sample(c(1,-1),1)
      value = theta.opt[i] + p*dif
      theta0 = c(theta0, value)
    } 
    
    norm(theta0,'2')
    print('difference')
    print(norm((theta.opt-theta0),'2'))
    Sys.sleep(1)
    #print(theta.opt)
    mod=FEMMR(data,theta0,sigma)
    theta.hat=mod$theta.hat
    r=norm(theta.opt-theta.hat,'2')
    
    r
    print('error')
    print(r)
  }, theta_norm, delta_j)
  
  
  res[[length(res)+1]]=rs
}

res1=res[[1]]
res2=res[[2]]

value=c()## is 3*7 vector
for(i in 1:length(theta_norm)){
  value=c(value,res1[[i]],res2[[i]])
}
theta_norm_label=c(rep("1", 2), rep("1.5", 2),rep("2", 2),rep("2.5", 2),rep("3", 2),rep("3.5", 2),rep("4", 2))
delta_label=rep(c("equal" , "over") , 7)

data <- data.frame(theta_norm_label,delta_label,value)
library("ggplot2")
# Grouped
ggplot(data, aes(fill=delta_label, y=value, x=theta_norm_label)) + 
  geom_bar(position="dodge", stat="identity")










