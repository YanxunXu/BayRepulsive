##########################
const <- function(theta, tau){
  a = 1/(4*tau^2)
  b = 1/(theta^2)
  c = sqrt(a^2 + 2*a*b)
  sum = a + b + c
  lambda = rep(0, 20)
  for (i in 1:20){
    lambda[i] = sqrt(2*a/sum)*(b/sum)^(i-1)
  }
  ind<-combn(lambda,K)
  return(sum(apply(ind,2,prod)))
}

update.C<-function(Z,K,theta,tau){
  #given features, which is in a matrix Z; number of cell type K;
  #two hyperparameters theta and tau(theta is the theta in your paper, tau is sigma_{q} in your paper)
  #compute the matrix C, while det(C) is the likelihood of DPP
  C<-matrix(1,nrow=K,ncol=K)
  q<-rep(0,K)
  for(i in 1:K){
    #compute q(X)
    q[i]<-exp(-sum((Z[i,]-mu0)^2)/(2*tau*tau))
  }
  for(i in 1:K){
    #compute C(x,x')
    for(j in i:K){
      C[i,j]<-exp(-sum((Z[i,]-Z[j,])^2)/(theta^2))
    }
    for(j in 1:(i-1)){
      C[i,j]<-C[j,i]
    }
  }
  C<-t(t(C)*q)*q
  C+0.00001*diag(K)
}

likeli.theta<-function(theta){
  #function to minimize in order to update theta
  C<-update.C(Z.now,K,theta,tau)
  Likeli<-log(det(C))-log(const(theta,tau))
  -Likeli-dgamma(theta,20,2,log=TRUE)
}

likeli.tau<-function(tau){
  #function to minimize in order to update tau
  C<-update.C(Z.now,K,theta,tau)
  Likeli<-log(det(C))-log(const(theta,tau))
  -Likeli-dgamma(tau,2,2,log=TRUE)
}
#updata W
fnw<-function(x){
  mu<-t(as.matrix(x))%*%Z.now
  s<- sum((mu-data.now[i,])^2)/(2*sigma.square)
  s
}


fnw_new<- function(x){
  w_temp = exp(x-max(x))/sum(exp(x-max(x)))
  mu<-t(as.matrix(w_temp))%*%Z.now
  s<- sum((mu-data.now[i,])^2)/(2*sigma.square)
  s
}

grw<-function(x){
  mu<-t(as.matrix(x))%*%Z.now
  ans<-rep(0,K)
  for(j in 1:K){
    s<-0
    for(k in 1:Nfeature){
      s<-s+(mu[k]-data.now[i,k])*Z.now[j,k]/sigma.square
    }
    ans[j]<-s
  }
  ans
}

heqw<-function(x){
  h<-rep(NA,1)
  h[1]<-sum(x)-1
  h
}

heq.jacw<-function(x){
  m<-matrix(NA,1,length(x))
  m[1,]<-rep(1, length(x))
  m
}

hin<-function(x){
  h<-rep(NA,length(x))
  for(i in 1:length(x)){
    h[i]<-x[i]-0.00001
  }
  h
}

hin.jac<-function(x){
  diag(length(x))
}


#update Z
fnz1<-function(x){
  Z0<-Z.now
  Z0[,i]<-x
  mu<-W%*%Z0
  s<-sum((mu[,i]-data.now[,i])^2)/(2*sigma.square)
  C<-update.C(Z0,K,theta,tau)
  s-log(det(C))*temperature
}

fnz<-function(x){
  Z0<-Z.now
  Z0[ind[i],i]<-x
  mu<-W%*%Z0
  s<-sum((mu[,i]-data.now[,i])^2)/(2*sigma.square)
  C<-update.C(Z0,K,theta,tau)
  s-log(det(C))*temperature
}

grz1<-function(x){
  Z0<-Z.now
  Z0[,i]<-x
  mu<-W%*%Z0
  C<-update.C(Z0,K,theta,tau)
  C.inv<-solve(C)
  ans<-rep(0,K)
  for(k in 1:K){
    s<-0
    for(j in 1:Nobs){
      s<-s+(mu[j,i]-data.now[j,i])*W[j,k]/sigma.square
    }
    s1<-0
    for(l in 1:K){
      s1<-s1+C.inv[k,l]*C[k,l]*((Z0[k,i]-mu0[i])/(tau^2)+(2*(Z0[k,i]-Z0[l,i]))/(theta^2))
    }
    ans[k]<-s+2*s1*temperature
  }
  ans
}

grz<-function(x){
  Z0<-Z.now
  Z0[ind[i],i]<-x
  mu<-W%*%Z0
  C<-update.C(Z0,K,theta,tau)
  C.inv<-solve(C)
  s<-0
  for(j in 1:Nobs){
    s<-s+(mu[j,i]-data.now[j,i])*W[j,ind[i]]/sigma.square
  }
  s1<-0
  for(l in 1:K){
    s1<-s1+C.inv[ind[i],l]*C[ind[i],l]*((x-mu0[i])/(tau^2) + (2*(x-Z0[l,i])/(theta^2)))
    #print(s1)
  }
  ans<-s+2*s1*temperature
  ans
}

#update sigma.square
update.sigma<-function(W,Z,DATA){
  mu<-DATA-(W%*%Z)
  return(sum(mu^2)/(Nobs*Nfeature))
}
#############


#' BayRepulsive_unknown is the NMF function used when the number of subclones is unknown
#'
#' Takes in the observed data matrix, the range of number of subclones, the number of features and samples
#' Gives the estimation of the NMF, including the estimated number of subclones.
#' @usage BayRepulsive_unknown(Datause, K_min, K_max, Nobs, Nfeature, Niter = 100, epsilon = 0.0001,tau = 100, seed = 1 )
#' @param Datause  The observed data matrix. Each row is a sample.
#' @param K_min The minimum number of subclones
#' @param K_max The Maximum number of subclones
#' @param Nobs The number of samples, i.e., the number of rows of the Datause
#' @param Nfeature The number of features, i.e., the number of columns of the Datause
#' @param Niter The number of maximum iterations
#' @param epsilon break if the L2 distance of the two estiamted proportion matrix in row is less than epsilon
#' @param tau The hyperparameter for DPP, a large number is prefered, default value is 100
#' @param seed The random seed, default as 1
#' @return W, the estiamted signature matrix.
#' @return Z, and the estiamted number of subclones.
#' @return K, the estimated number of subclones.
#' @details Given an observed matrix, whose rows are mixed samples of unknown number of subclones,
#' we give an estimation of number of subclones, and do the NMF.
#'
#' This function will create a bunch of globel variables, named "Datause", "Nobs", "Nfeature", "sigma0",
#' "mu0", "K", "Theta", "W.star", "Z.star",  "W_temp", "sigma.square", "data.now", "Z.now", "i".
#' Thus, users should avoid these variable names when using BayRepulsive_unknown, if they don't want the variables to be overwritten.
#' Especially "i", which is commonly used in loops.
#' @source BayRepulsive: A Bayesian Repulsive Deconvolution Model for Inferring Latent Biologic Structure
#' @export
BayRepulsive_unknown<-function(Datause, K_min, K_max, Nobs, Nfeature, Niter = 100, epsilon = 0.0001, tau = 100, seed = 1 ){
  library(mvtnorm)
  library(alabama)
  library(psych)
  library(optimx)

  N_K<- K_max - K_min + 1
  C<-rep(0,N_K)
  result<-NULL
  result$W<-array(NA,dim = c(Nobs,K_max,N_K))
  result$Z<-array(NA,dim = c(K_max,Nfeature,N_K))
  assign('Datause', Datause, envir = .GlobalEnv)
  assign('Nobs', Nobs, envir = .GlobalEnv)
  assign('Nfeature', Nfeature, envir = .GlobalEnv)
  assign('sigma0', apply(Datause,2,var), envir = .GlobalEnv)
  assign('mu0', apply(Datause,2,mean), envir = .GlobalEnv)
  assign('temperature', Nobs * Nfeature, envir = .GlobalEnv)
  assign('tau', tau, envir = .GlobalEnv)

  for(K in K_min:K_max){
    set.seed(seed)
    assign('K', K, envir = .GlobalEnv)
    #initial data
    cat(paste0("Number of subclones ", as.character(K), "...", sep=" "));

    assign('Theta', matrix(0,nrow = Nobs, ncol = K), envir = .GlobalEnv)
    assign('W.star', matrix(0,nrow = Nobs, ncol = K), envir = .GlobalEnv)

    for(i in 1:Nobs){
      .GlobalEnv$Theta[i,]<-rgamma(K,1/K,1)
      .GlobalEnv$W.star[i,]<-.GlobalEnv$Theta[i,]/sum(.GlobalEnv$Theta[i,])
    }

    assign('Z.star', pmax(rmvnorm(K, mu0, diag(sigma0)),0) , envir = .GlobalEnv)
    assign('W', W.star, envir = .GlobalEnv)
    assign('W_temp', W, envir = .GlobalEnv)
    assign('W_former', W, envir = .GlobalEnv)
    assign('sigma.square', 1, envir = .GlobalEnv)



    assign('data.now', .GlobalEnv$Datause, envir = .GlobalEnv)
    assign('Z.now', Z.star, envir = .GlobalEnv)


    for(iter in 1:Niter){
      #theta<-auglag(par = theta, fn = likeli.theta, hin = hin, hin.jac = hin.jac,control.outer = list(trace=F))$par
      assign('theta', optimize(f = likeli.theta, interval = c(0, 100))$minimum, envir = .GlobalEnv)
      for(i in 1:Nobs){
        assign('i',i,envir = .GlobalEnv)
        .GlobalEnv$W_temp[i,]<-nlminb(.GlobalEnv$W_temp[i,], fnw_new, lower = -10, upper = 10)$par
        .GlobalEnv$W[i,] <- exp(.GlobalEnv$W_temp[i,]-max(.GlobalEnv$W_temp[i,]))/sum(exp(.GlobalEnv$W_temp[i,]-max(.GlobalEnv$W_temp[i,])))
      }
      for(i in 1:Nfeature){
        assign('i',i,envir = .GlobalEnv)
        z<-.GlobalEnv$Z.now[,i]
        z<-pmax(z,0)
        #z0<-optim(par = z,fn=fnz1,gr = grz1,lower = rep(0,K),method = "L-BFGS-B",control=list(trace=0))$par
        z0<-nlminb(z, fnz1, lower = 0)$par
        .GlobalEnv$Z.now[,i]<-z0
      }
      assign('sigma.square', update.sigma(.GlobalEnv$W, .GlobalEnv$Z.now, .GlobalEnv$data.now), envir = .GlobalEnv)

      if(iter>1){
        if(sum((.GlobalEnv$W-.GlobalEnv$W_former)^2)<epsilon){
          break
        }
      }
      .GlobalEnv$W_former<-.GlobalEnv$W
    }
    temp<-(Datause) - .GlobalEnv$W%*%.GlobalEnv$Z.now
    C[K-1]<-sum(temp^2)
    result$W[,1:K,K-1]<-.GlobalEnv$W
    result$Z[1:K,,K-1]<-.GlobalEnv$Z.now
    cat("is done!"); cat("\n")
  }
  k_star<-which.max(diff(diff(C))) + 2
  W<-result$W[,1:k_star,k_star-1]
  Z<-result$Z[1:k_star,,k_star-1]

  return(list(W=W, Z=Z, k =k_star))
}






#' BayRepulsive_known is the NMF function used when the number of subclones is known
#'
#' Takes in the observed data matrix, the number of subclones, the number of features and samples.
#' Gives the estiamted NMF result.
#' @usage BayRepulsive_known(Datause, K, Nobs, Nfeature,Niter = 100, epsilon = 0.0001, tau = 100, seed = 1 )
#' @param Datause  The observed data matrix. Each row is a sample.
#' @param K The number of subclones
#' @param Nobs The number of samples, i.e., the number of rows of the Datause
#' @param Nfeature The number of features, i.e., the number of columns of the Datause
#' @param Niter The number of maximum iterations
#' @param epsilon break if the L2 distance of the two estiamted proportion matrix in row is less than epsilon
#' @param tau The hyperparameter for DPP, a large number is prefered, default value is 100
#' @param seed The random seed, default as 1
#' @return W, the estiamted signature matrix.
#' @return Z, and the estiamted number of subclones.
#' @return C, an measure of performance for NMF.
#' @details Given an observed matrix, whose rows are mixed samples of unknown number of subclones,
#' we give an estimation of number of subclones, and do the NMF.
#'
#' This function will create a bunch of globel variables, named "Datause", "Nobs", "Nfeature", "sigma0",
#' "mu0", "K", "Theta", "W.star", "Z.star",  "W_temp", "sigma.square", "data.now", "Z.now", "i".
#' Thus, users should avoid these variable names when using BayRepulsive_unknown, if they don't want the variables to be overwritten.
#' Especially "i", which is commonly used in loops.
#' @source BayRepulsive: A Bayesian Repulsive Deconvolution Model for Inferring Latent Biologic Structure
#' @export
BayRepulsive_known<-function(Datause, K, Nobs, Nfeature, Niter = 100, epsilon = 0.0001, tau = 100, seed = 1 ){
  library(mvtnorm)
  library(alabama)
  library(psych)
  library(optimx)
  assign('Datause', Datause, envir = .GlobalEnv)
  assign('K', K, envir = .GlobalEnv)
  assign('Nobs', Nobs, envir = .GlobalEnv)
  assign('Nfeature', Nfeature, envir = .GlobalEnv)
  assign('sigma0', apply(Datause,2,var), envir = .GlobalEnv)
  assign('mu0', apply(Datause,2,mean), envir = .GlobalEnv)
  assign('temperature', Nobs * Nfeature, envir = .GlobalEnv)
  assign('tau', tau, envir = .GlobalEnv)
  set.seed(seed)
  #initial data

  assign('Theta', matrix(0,nrow = Nobs, ncol = K), envir = .GlobalEnv)
  assign('W.star', matrix(0,nrow = Nobs, ncol = K), envir = .GlobalEnv)

  for(i in 1:Nobs){
    .GlobalEnv$Theta[i,]<-rgamma(K,1/K,1)
    .GlobalEnv$W.star[i,]<-.GlobalEnv$Theta[i,]/sum(.GlobalEnv$Theta[i,])
  }

  assign('Z.star', pmax(rmvnorm(K, mu0, diag(sigma0)),0) , envir = .GlobalEnv)
  assign('W', W.star, envir = .GlobalEnv)
  assign('W_temp', W, envir = .GlobalEnv)
  assign('W_former', W, envir = .GlobalEnv)
  assign('sigma.square', 1, envir = .GlobalEnv)



  assign('data.now', .GlobalEnv$Datause, envir = .GlobalEnv)
  assign('Z.now', Z.star, envir = .GlobalEnv)


  for(iter in 1:Niter){
    #theta<-auglag(par = theta, fn = likeli.theta, hin = hin, hin.jac = hin.jac,control.outer = list(trace=F))$par
    assign('theta', optimize(f = likeli.theta, interval = c(0, 100))$minimum, envir = .GlobalEnv)
    for(i in 1:Nobs){
      assign('i',i,envir = .GlobalEnv)
      .GlobalEnv$W_temp[i,]<-nlminb(.GlobalEnv$W_temp[i,], fnw_new, lower = -10, upper = 10)$par
      .GlobalEnv$W[i,] <- exp(.GlobalEnv$W_temp[i,]-max(.GlobalEnv$W_temp[i,]))/sum(exp(.GlobalEnv$W_temp[i,]-max(.GlobalEnv$W_temp[i,])))
    }
    for(i in 1:Nfeature){
      assign('i',i,envir = .GlobalEnv)
      z<-.GlobalEnv$Z.now[,i]
      z<-pmax(z,0)
      #z0<-optim(par = z,fn=fnz1,gr = grz1,lower = rep(0,K),method = "L-BFGS-B",control=list(trace=0))$par
      z0<-nlminb(z, fnz1, lower = 0)$par
      .GlobalEnv$Z.now[,i]<-z0
    }
    assign('sigma.square', update.sigma(.GlobalEnv$W, .GlobalEnv$Z.now, .GlobalEnv$data.now), envir = .GlobalEnv)

    if(iter>1){
      if(sum((.GlobalEnv$W-.GlobalEnv$W_former)^2)<epsilon){
        break
      }
    }
    .GlobalEnv$W_former<-.GlobalEnv$W
  }
  temp<-(Datause) - W%*%Z.now
  C<-sum(temp^2)

  return(list(W=W, Z=Z.now, C = C))
}






