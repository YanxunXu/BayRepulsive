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
  test.env$ind<-combn(lambda,test.env$K)
  return(sum(apply(test.env$ind,2,prod)))
}

update.C<-function(Z,K,theta,tau){
  #given features, which is in a matrix Z; number of cell type K;
  #two hyperparameters theta and tau(theta is the theta in your paper, tau is sigma_{q} in your paper)
  #compute the matrix C, while det(C) is the likelihood of DPP
  C<-matrix(1,nrow=K,ncol=K)
  q<-rep(0,K)
  for(i in 1:K){
    #compute q(X)
    q[i]<-exp(-sum((Z[i,]-test.env$mu0)^2)/(2*tau*tau))
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
  C<-update.C(test.env$Z.now,test.env$K,theta,test.env$tau)
  Likeli<-log(det(C))-log(const(theta,test.env$tau))
  -Likeli-dgamma(theta,test.env$a.theta,test.env$b.theta,log=TRUE)
}

#updata W

fnw<-function(x){
  mu<-t(as.matrix(x))%*%test.env$Z.now
  s<- sum((mu-test.env$data.now[test.env$i,])^2)/(2*test.env$sigma.square)
  s
}

fnw_new<- function(x){
  w_temp = exp(x-max(x))/sum(exp(x-max(x)))
  mu<-t(as.matrix(w_temp))%*%test.env$Z.now
  s<- sum((mu-test.env$data.now[test.env$i,])^2)/(2*test.env$sigma.square)
  s
}

grw<-function(x){
  mu<-t(as.matrix(x))%*%test.env$Z.now
  ans<-rep(0,test.env$K)
  for(j in 1:test.env$K){
    s<-0
    for(k in 1:test.env$Nfeature){
      s<-s+(mu[k]-test.env$data.now[test.env$i,k])*test.env$Z.now[j,k]/test.env$sigma.square
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
  Z0<-test.env$Z.now
  Z0[,test.env$i]<-x
  mu<-test.env$W%*%Z0
  s<-sum((mu[,test.env$i]-test.env$data.now[,test.env$i])^2)/(2*test.env$sigma.square)
  C<-update.C(Z0,test.env$K,test.env$theta,test.env$tau)
  s-log(det(C))*test.env$temperature
}

fnz<-function(x){
  Z0<-test.env$Z.now
  Z0[ind[test.env$i],test.env$i]<-x
  mu<-test.env$W%*%Z0
  s<-sum((mu[,i]-test.env$data.now[,i])^2)/(2*test.env$sigma.square)
  C<-update.C(Z0,test.env$K,test.env$theta,test.env$tau)
  s-log(det(C))*test.env$temperature
}

grz1<-function(x){
  Z0<-test.env$Z.now
  Z0[,test.env$i]<-x
  mu<-test.env$W%*%Z0
  C<-update.C(Z0,test.env$K,test.env$theta,test.env$tau)
  C.inv<-solve(C)
  ans<-rep(0,test.env$K)
  for(k in 1:test.env$K){
    s<-0
    for(j in 1:test.env$Nobs){
      s<-s+(mu[j,test.env$i]-data.now[j,test.env$i])*test.env$W[j,k]/test.env$sigma.square
    }
    s1<-0
    for(l in 1:test.env$K){
      s1<-s1+C.inv[k,l]*C[k,l]*((Z0[k,i]-test.env$mu0[i])/(test.env$tau^2)+(2*(Z0[k,test.env$i]-Z0[l,test.env$i]))/(test.env$theta^2))
    }
    ans[k]<-s+2*s1*test.env$temperature
  }
  ans
}

grz<-function(x){
  Z0<-test.env$Z.now
  Z0[ind[test.env$i],test.env$i]<-x
  mu<-test.env$W%*%Z0
  C<-update.C(Z0,test.env$K,test.env$theta,test.env$tau)
  C.inv<-solve(C)
  s<-0
  for(j in 1:test.env$Nobs){
    s<-s+(mu[j,test.env$i]-data.now[j,test.env$i])*test.env$W[j,ind[test.env$i]]/test.env$sigma.square
  }
  s1<-0
  for(l in 1:K){
    s1<-s1+C.inv[ind[test.env$i],l]*C[ind[test.env$i],l]*((x-test.env$mu0[test.env$i])/(test.env$tau^2) + (2*(x-Z0[l,test.env$i])/(test.env$theta^2)))
    #print(s1)
  }
  ans<-s+2*s1*test.env$temperature
  ans
}

#update sigma.square
update.sigma<-function(W,Z,DATA){
  mu<-DATA-(W%*%Z)
  return(sum(mu^2)/(Nobs*Nfeature))
}

#-------------
#' BayRepulsive_known is a deconvolution function designed for inferring tumor heterogeneity, when the number of subclones is known.
#'
#' This function gives deconvolution results of the observed matrix, along with the square sum of the residuals.
#' @usage BayRepulsive_known(DATA, K, Nobs, Nfeature,
#'                   Niter = 100, epsilon = 0.0001, tau = 100,
#'                   a.theta = 0.01, b.theta = 0.01, seed = 1 )
#' @param DATA  The observed data matrix. Each row represents a feature (gene); each column represents a sample.
#' @param K The number of subclones.
#' @param Nobs The number of samples, i.e., the number of columns of the \code{DATA}.
#' @param Nfeature The number of features, i.e., the number of rows of the \code{DATA}.
#' @param Niter The number of maximum iterations.
#' @param epsilon Tolerance for convergence. We determine whether to break based on the estimated weight matrix.
#' We decide to break if the distance induced by L2 norm between two successive estimated weight matrices is less than epsilon.
#' @param tau The hyperparameter for DPP. A large number is preferred. See more in \strong{Details}.
#' @param a.theta The hyperparameter for DPP. See more in \strong{Details}.
#' @param b.theta The hyperparameter for DPP. See more in \strong{Details}.
#' @param seed The random seed.
#' @return A list of following components:
#' \tabular{lllllllllll}{
#' \tab \code{Z} \tab  \tab  \tab    The estimated subclone-specific expression matrix.\tab  \tab  \tab \tab  \tab  \tab    \cr
#' \tab \code{W} \tab  \tab  \tab    The estimated weight matrix. \tab  \tab  \tab \tab  \tab  \tab   \cr
#' \tab \code{C} \tab  \tab  \tab    Square sum of the residuals used as a measure of performance.\tab  \tab  \tab  \tab  \tab  \tab
#' }
#' @details Given an observed matrix, whose columns are mixed samples of known number of subclones, the function returns the
#' deconvolution results.
#'
#' The deconvolution model is based on the assumption that
#' \deqn{Y = ZW+E.}
#' Here \eqn{Y} is the observed matrix \code{DATA};
#' \eqn{Z} is the subclone-specific expression matrix;
#' \eqn{W} is the weight matrix; \eqn{E} is the matrix whose entries are independent white noises, with unknown variance \eqn{\sigma^{2}}{\sigma^2}.
#' We assume each column of \eqn{W}, \eqn{W_j}{Wj} has a prior \eqn{W_j\sim\mathrm{Dir}(\alpha)}{Wj~Dir(\alpha)}, where \eqn{\alpha} is a vector with elements 1.
#' We also assume an improper uniform prior for \eqn{\sigma^2}: \eqn{\sigma^2\sim\mathrm{Uniform}(0, 10^6)}{\sigma^2~U(0,10^6)}.
#' We use a fixed-size determinant point process \insertCite{kulesza2011k}{BayRepulsive}
#' as a prior for the subclone-specific expression matrix \eqn{Z}.
#' Suppose there are \eqn{K} subclones and let \eqn{Z_k}{Zk} be the expression profile of subclone \eqn{k}.
#' Mean zero multivariate normal density functions are commonly used as quality functions in DPP.
#' Since the subclone-specific expression matrix is nonnegative, we cosider a transformation,
#' \eqn{\tilde{Z}_k = Z_k-\bar{Y}}{\tilde{Z}k = Zk - \bar{Y}}, where the vector \eqn{\bar{Y}} is the mean of average expression level in the \code{DATA} of each gene.
#' The prior of (\eqn{\tilde{Z}_1,\dots, \tilde{Z}_K}{\tilde{Z}1,...,\tilde{Z}K}) is proportional to the determinant of a \eqn{K\times K}{K*K} matrix \eqn{L},
#' defined by \eqn{L_{ij} = q(\tilde{Z}_i)\phi(\tilde{Z}_i, \tilde{Z}_j)q(\tilde{Z}_j)}{Lij = q(\tilde{Z}i)\phi(\tilde{Z}i,\tilde{Z}j)q(\tilde{Z}j)},
#' where \eqn{\phi(\tilde{Z}_i, \tilde{Z}_j) = \exp\{ - \frac{\|\tilde{Z}_i-\tilde{Z}_j\|^2}{\theta^2}\}}{\phi(\tilde{Z}i,\tilde{Z}j)=exp\{||\tilde{Z}i-\tilde{Z}j||^2 / \theta^2\}}
#' and \eqn{q(\tilde{Z}_j)}{q(\tilde{Z}j)} is the density function of a multivariate normal distribution with mean being the zero vector and variance being \eqn{\tau^2 I}.
#' Here \eqn{\tau} is the parameter \code{tau} in the function,
#' and \eqn{\theta} is an unknown parameter with a prior \eqn{\theta\sim\mathrm{Gamma}(a_{\theta},b_{\theta})}{theta~Gamma(a.theta,b.theta)}.
#' Here \eqn{a_\theta}{a.theta} and \eqn{b_\theta}{b.theta} are the parameters \code{a.theta} and \code{b.theta} in the function.
#'
#' @source BayRepulsive: A Bayesian Repulsive Deconvolution Model for Inferring Tumor Heterogeneity
#' @examples
#' rm(list=ls())
#' library(BayRepulsive)
#' data(CCLE)
#' set.seed(1)
#' Nobs     <- 24
#' Nfeature <- 100
#' K0       <- 3
#' ### randomly generate weight matrix W for 24 mixing samples
#' W        <- matrix(0,nrow = K0, ncol = Nobs)
#' for(i in 1:Nobs){
#'   Theta <- rgamma(K0,1/K0,1)
#'   W[,i] <- Theta/sum(Theta)
#' }
#' ### add some noise
#' error    <- t(matrix(rnorm(Nfeature * Nobs, mean = 0, sd = 0.5), nrow = Nobs))
#' DATA     <- CCLE$Z%*%W + error
#' ### Note: please make sure that there are no negative values after adding the noise
#' result1  <- BayRepulsive_known(DATA = DATA, K = K0, Nobs = Nobs,
#'                                Nfeature = Nfeature)
#' cor(as.vector(result1$W), as.vector(W))
#'


#' @keywords functions
#' @references
#' \insertRef{kulesza2011k}{BayRepulsive}
#' @export
BayRepulsive_known<-function(DATA, K, Nobs, Nfeature, Niter = 100, epsilon = 0.0001, tau = 100,
                             a.theta = 0.01, b.theta = 0.01, seed = 1 ){
  require(mvtnorm)
  require(alabama)
  require(psych)
  require(optimx)
  Datause <- t(DATA)
  assign('test.env', new.env(), envir = .GlobalEnv)
  assign('Datause', Datause, envir = test.env)
  assign('K', K, envir = test.env)
  assign('Nobs', Nobs, envir = test.env)
  assign('a.theta', a.theta, envir = test.env)
  assign('b.theta', b.theta, envir = test.env)
  assign('Nfeature', Nfeature, envir = test.env)
  assign('sigma0', apply(Datause,2,var), envir = test.env)
  assign('mu0', apply(Datause,2,mean), envir = test.env)
  assign('temperature', Nobs * Nfeature, envir = test.env)
  assign('tau', tau, envir = test.env)
  set.seed(seed)
  #initial data

  assign('Theta', matrix(0,nrow = Nobs, ncol = K), envir = test.env)
  assign('W.star', matrix(0,nrow = Nobs, ncol = K), envir = test.env)

  for(i in 1:test.env$Nobs){
    test.env$Theta[i,]<-rgamma(K,1/K,1)
    test.env$W.star[i,]<-test.env$Theta[i,]/sum(test.env$Theta[i,])
  }

  assign('Z.star', pmax(rmvnorm(K, test.env$mu0, diag(test.env$sigma0)),0) , envir = test.env)
  assign('W', test.env$W.star, envir = test.env)
  assign('W_temp', test.env$W, envir = test.env)
  assign('W_former', test.env$W, envir = test.env)
  assign('sigma.square', 1, envir = test.env)



  assign('data.now', test.env$Datause, envir = test.env)
  assign('Z.now', test.env$Z.star, envir = test.env)


  for(iter in 1:Niter){
    #theta<-auglag(par = theta, fn = likeli.theta, hin = hin, hin.jac = hin.jac,control.outer = list(trace=F))$par
    assign('theta', optimize(f = likeli.theta, interval = c(0, 100))$minimum, envir = test.env)
    for(i in 1:Nobs){
      assign('i',i,envir = test.env)
      test.env$W_temp[i,]<-nlminb(test.env$W_temp[i,], fnw_new, lower = -10, upper = 10)$par
      test.env$W[i,] <- exp(test.env$W_temp[i,]-max(test.env$W_temp[i,]))/sum(exp(test.env$W_temp[i,]-max(test.env$W_temp[i,])))
    }
    for(i in 1:Nfeature){
      assign('i',i,envir = test.env)
      z<-test.env$Z.now[,i]
      z<-pmax(z,0)
      #z0<-optim(par = z,fn=fnz1,gr = grz1,lower = rep(0,K),method = "L-BFGS-B",control=list(trace=0))$par
      z0<-nlminb(z, fnz1, lower = 0)$par
      test.env$Z.now[,i]<-z0
    }
    assign('sigma.square', update.sigma(test.env$W, test.env$Z.now, test.env$data.now), envir = test.env)

    if(iter>1){
      if(sum((test.env$W-test.env$W_former)^2)<epsilon){
        break
      }
    }
    test.env$W_former<-test.env$W
  }
  temp<-(Datause) - test.env$W%*%test.env$Z.now
  C<-sum(temp^2)
  W<- t(test.env$W)
  Z<- t(test.env$Z.now)
  rm(list=ls(test.env), envir = test.env)
  return(list( W= W, Z=Z, C = C))
}


#--------------------
#' BayRepulsive_unknown is a deconvolution function designed for inferring tumor heterogeneity,
#' when the number of subclones is unknown.
#'
#' This function gives the estimated number of subclones, along with deconvolution results of the observed matrix.
#' @usage BayRepulsive_unknown(DATA, K_min, K_max, Nobs, Nfeature,
#'                       Niter = 100, epsilon = 0.0001, tau = 100,
#'                       a.theta = 0.01, b.theta = 0.01, seed = 1 )
#' @param DATA  The observed data matrix. Each row represents a feature (gene); each column represents a sample.
#' @param K_min The minimum number of subclones.
#' @param K_max The Maximum number of subclones.
#' @param Nobs The number of samples, i.e., the number of columns of the \code{DATA}.
#' @param Nfeature The number of features, i.e., the number of rows of the \code{DATA}.
#' @param Niter The number of maximum iterations.
#' @param epsilon Tolerance for convergence. We determine whether to break based on the estimated weight matrix.
#' We decide to break if the distance induced by L2 norm between two successive estimated weight matrices is less than epsilon.
#' @param tau The hyperparameter for DPP. A large number is preferred. See \code{\link{BayRepulsive_known}} for more detials.
#' @param a.theta The hyperparameter for DPP. See \code{\link{BayRepulsive_known}} for more detials.
#' @param b.theta The hyperparameter for DPP. See \code{\link{BayRepulsive_known}} for more detials.
#' @param seed The random seed.
#'
#' @return A list of following components:
#' \tabular{llllllllllllll}{
#' \tab \code{Z} \tab  \tab  \tab    The estimated subclone-specific expression matrix.\tab  \tab  \tab \tab  \tab  \tab \tab  \tab  \tab    \cr
#' \tab \code{W} \tab  \tab  \tab    The estimated weight matrix. \tab  \tab  \tab \tab  \tab  \tab  \tab  \tab  \tab  \cr
#' \tab \code{K_hat} \tab  \tab  \tab    The number of estimated subclones.\tab  \tab  \tab  \tab  \tab  \tab \tab  \tab  \tab
#' }
#' @details Given an observed matrix, whose columns are mixed samples of unknown number of subclones,
#' this function gives an estimation of number of subclones along with deconvolution results.
#'
#' We first use the algorithm in \code{\link{BayRepulsive_known}} to fit the data for every possible number of subclones.
#' Let \eqn{S(k)} denote the square sum of the residuals when the number of subclones is fixed at \eqn{k}.
#' We define the second-order finite difference \eqn{\Delta^2 S(k)}{Delta^2 S(k)} of the residual by
#' \eqn{\Delta^2S(k) = [S(k+1)-S(k)]-[S(k)-S(k-1)]}{Delta^2(k) = [S(k+1)-S(k)]-[S(k)-S(k-1)]}, for \eqn{k=K_{min}+1,\dots,K_{max}-1}{k=Kmin+1,...,Kmax-1}.
#' Then the optimal \eqn{\hat{K}} estimated by BayRepulsive is
#' \deqn{\hat{K}=\mathrm{arg max}_k \Delta^2 S(k).}{\hat{K} = argmin Delta^2 S(k).}
#' And the deconvolution results are the corresponding results when the number of subclones is fixed at \eqn{\hat{K}}.
#' @source BayRepulsive: A Bayesian Repulsive Deconvolution Model for Inferring Tumor Heterogeneity
#' @examples
#' rm(list=ls())
#' library(BayRepulsive)
#' data(CCLE)
#' set.seed(1)
#' Nobs     <- 24
#' Nfeature <- 100
#' K0       <- 3
#' ### randomly generate weight matrix W for 24 mixing samples
#' W        <- matrix(0,nrow = K0, ncol = Nobs)
#' for(i in 1:Nobs){
#'   Theta <- rgamma(K0,1/K0,1)
#'   W[,i] <- Theta/sum(Theta)
#' }
#' ### add some noise
#' error    <- t(matrix(rnorm(Nfeature * Nobs, mean = 0, sd = 0.5), nrow = Nobs))
#' DATA     <- CCLE$Z%*%W + error
#' ### Note: please make sure that there are no negative values after adding the noise
#' result1  <- BayRepulsive_unknown(DATA = DATA, K_min = 2, K_max = 6, Nobs = Nobs,
#'                                  Nfeature = Nfeature)
#' cor(as.vector(result1$W), as.vector(W))
#' @keywords functions
#' @export
#' @seealso \code{\link{BayRepulsive_known}}
BayRepulsive_unknown<-function(DATA, K_min, K_max, Nobs, Nfeature, Niter = 100, epsilon = 0.0001, tau = 100,
                               a.theta = 0.01, b.theta = 0.01, seed = 1 ){
  require(mvtnorm)
  require(alabama)
  require(psych)
  require(optimx)
  Datause <- t(DATA)

  N_K<- K_max - K_min + 1
  C<-rep(0,N_K)
  result<-NULL
  assign('test.env', new.env(), envir = .GlobalEnv)
  result$W<-array(NA,dim = c(Nobs,K_max,N_K))
  result$Z<-array(NA,dim = c(K_max,Nfeature,N_K))
  assign('Datause', Datause, envir = test.env)
  assign('Nobs', Nobs, envir = test.env)
  assign('Nfeature', Nfeature, envir = test.env)
  assign('sigma0', apply(Datause,2,var), envir = test.env)
  assign('mu0', apply(Datause,2,mean), envir = test.env)
  assign('temperature', Nobs * Nfeature, envir = test.env)
  assign('tau', tau, envir = test.env)
  assign('a.theta', a.theta, envir = test.env)
  assign('b.theta', b.theta, envir = test.env)

  for(K in K_min:K_max){
    set.seed(seed)
    assign('K', K, envir = test.env)
    #initial data
    cat(paste0("Number of subclones ", as.character(K), "...", sep=" "));

    assign('Theta', matrix(0,nrow = Nobs, ncol = K), envir = test.env)
    assign('W.star', matrix(0,nrow = Nobs, ncol = K), envir = test.env)

    for(i in 1:Nobs){
      assign('i',i,envir = test.env)
      test.env$Theta[test.env$i,]<-rgamma(K,1/K,1)
      test.env$W.star[test.env$i,]<-test.env$Theta[test.env$i,]/sum(test.env$Theta[test.env$i,])
    }

    assign('Z.star', pmax(rmvnorm(K, test.env$mu0, diag(test.env$sigma0)),0) , envir = test.env)
    assign('W', test.env$W.star, envir = test.env)
    assign('W_temp', test.env$W, envir = test.env)
    assign('W_former', test.env$W, envir = test.env)
    assign('sigma.square', 1, envir = test.env)



    assign('data.now', test.env$Datause, envir = test.env)
    assign('Z.now', test.env$Z.star, envir = test.env)


    for(iter in 1:Niter){
      #theta<-auglag(par = theta, fn = likeli.theta, hin = hin, hin.jac = hin.jac,control.outer = list(trace=F))$par
      assign('theta', optimize(f = likeli.theta, interval = c(0, 100))$minimum, envir = test.env)
      for(i in 1:Nobs){
        assign('i',i,envir = test.env)
        test.env$W_temp[i,]<-nlminb(test.env$W_temp[i,], fnw_new, lower = -10, upper = 10)$par
        test.env$W[i,] <- exp(test.env$W_temp[i,]-max(test.env$W_temp[i,]))/sum(exp(test.env$W_temp[i,]-max(test.env$W_temp[i,])))
      }
      for(i in 1:Nfeature){
        assign('i',i,envir = test.env)
        z<-test.env$Z.now[,i]
        z<-pmax(z,0)
        #z0<-optim(par = z,fn=fnz1,gr = grz1,lower = rep(0,K),method = "L-BFGS-B",control=list(trace=0))$par
        z0<-nlminb(z, fnz1, lower = 0)$par
        test.env$Z.now[,i]<-z0
      }
      assign('sigma.square', update.sigma(test.env$W, test.env$Z.now, test.env$data.now), envir = test.env)

      if(iter>1){
        if(sum((test.env$W-test.env$W_former)^2)<epsilon){
          break
        }
      }
      test.env$W_former<-test.env$W
    }
    temp<-(Datause) - test.env$W%*%test.env$Z.now
    C[K-1]<-sum(temp^2)
    result$W[,1:K,K-1]<-test.env$W
    result$Z[1:K,,K-1]<-test.env$Z.now
    cat("is done!"); cat("\n")
  }
  k_star<-which.max(diff(diff(C))) + 2
  W<-result$W[,1:k_star,k_star-1]
  Z<-result$Z[1:k_star,,k_star-1]
  rm(list=ls(test.env), envir = test.env)
  return(list(W=t(W), Z=t(Z), K_hat =k_star))
}








