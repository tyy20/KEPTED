#' Kernel embedding of probability test for elliptical distribution
#'
#' This function gives a test on whether the data is elliptically distributed
#' based on kernel embedding of probability. See Tang and Li (2023+) for
#' details. Gaussian kernels and product-type inverse quadratic kernels are
#' considered.
#'
#' @param X A matrix with n rows and d columns.
#' @param eps The regularization constant added to the diagonal to avoid
#' singularity. Default value is \code{1e-6}.
#' @param kerU The type of kernel function on \code{U}. Currently supported
#' options are \code{"Gaussian"} and \code{"PIQ"}.
#' @param kerTheta The type of kernel function on \code{Theta}. Currently
#' supported options are \code{"Gaussian"} and \code{"PIQ"}.
#' @param gamma.U The tuning parameter \code{gamma} in the kernel function
#' k_U(u1,u2). If \code{gamma.U=0}, the recommended procedure of selecting
#' tuning parameter will be applied. Otherwise, the value given in
#' \code{gamma.U} will be directly used as the tuning parameter. Default value
#' is \code{gamma.U=0}. See "Details" for more information.
#' @param gamma.Theta The tuning parameter \code{gamma} in the kernel function
#' k_Theta(theta1,theta2). If \code{gamma.Theta=0}, the recommended procedure
#' of selecting tuning parameter will be applied. Otherwise, the value given in
#' \code{gamma.Theta} will be directly used as the tuning parameter. Default
#' value is \code{gamma.Theta=0}. See "Details" for more information.
#'
#' @details
#' The Gaussian kernel is defined as k(z1,z2)=exp(-gamma*||z1-z2||^2), and the
#' Product-type Inverse-Quadratic (\code{PIQ}) kernel is defines as
#' k(z1,z2)=Prod_j(1/(1+gamma*(z1_j-z2_j)^2)). The recommended procedure of
#' selecting tuning parameter is given as in the simulation section of Tang and
#' Li (2023+), where we set
#' 1/sqrt(gamma)=(n(n-1)/2)^(-1)*sum_\{1<=i<j<=n\}||Z_i-Z_j||.
#'
#' @return A list of the following:
#' \item{stat}{The value of the test statistic.}
#' \item{pval}{The p-value of the test.}
#' \item{lambda}{The \code{n} eigenvalues in the approximated asymptotic
#' distribution.}
#' \item{gamma.U}{The tuning parameter \code{gamma.U} used in the test. Same as
#' the input if its input is nonzero.}
#' \item{gamma.Theta}{The tuning parameter \code{gamma.Theta} used in the test.
#' Same as the input if its input is nonzero.}
#'
#' @references
#' \cite{Tang, Y. and Li, B. (2023+), “A nonparametric test for elliptical
#' distribution based on kernel embedding of probabilities,”
#' \url{https://arxiv.org/abs/2306.10594}}
#'
#' @note
#' In the arguments, \code{eps} refers to a regularization constant added to
#' the diagonal. When the dimension is high, we recommend increasing \code{eps}
#' to avoid singularity.
#'
#' @examples
#' set.seed(313)
#' n=200
#' d=3
#'
#' ## Null Hypothesis
#' X=matrix(rnorm(n*d),nrow=n,ncol=d)
#' EllKEPT(X)
#'
#' ## Alternative Hypothesis
#' X=matrix(rchisq(n*d,2)-2,nrow=n,ncol=d)
#' EllKEPT(X)
#'
#' @export

EllKEPT=function(X,eps=1e-6,kerU="Gaussian",kerTheta="Gaussian",gamma.U=0,gamma.Theta=0){
  n=nrow(X)
  d=ncol(X)

  mu.hat=colMeans(X)
  Sigma.hat=cov(X)*(n-1)/n


  Sigma.hat.sqrt=sqrtm(Sigma.hat+diag(eps,d))
  Sigma.hat.inv=solve(Sigma.hat+diag(eps,d))
  Sigma.hat.sqrt.inv=sqrtm(Sigma.hat.inv+diag(eps,d))


  X.hat=X
  for(i in 1:n){
    X.hat[i,]=Sigma.hat.sqrt.inv%*%(X[i,]-mu.hat)
  }

  U=sqrt(rowSums(X.hat^2))
  V=X.hat
  for(i in 1:n){
    V[i,]=X.hat[i,]/U[i]
  }


  Theta=matrix(NA,n,d-1)
  for(i in 1:n){
    Theta[i,]=Rec2Polar(X.hat[i,])$Theta
  }



  if(gamma.U==0){
    sum.gamma.U=0
    for(k in 1:(n-1)){
      for(l in (k+1):n){
        sum.gamma.U=sum.gamma.U+abs(U[k]-U[l])
      }
    }
    gamma.U=(sum.gamma.U/(n*(n-1)/2))^(-2)
  }


  if(kerU=="Gaussian"){
    kernel.U=function(u1,u2){
      exp(-gamma.U*(u1-u2)^2)
    }
  }else if(kerU=="PIQ"){
    kernel.U=function(u1,u2){
      1/(1+gamma.U*(u1-u2)^2)
    }
  }

  if(gamma.Theta==0){
    sum.gamma.Theta=0
    for(k in 1:(n-1)){
      for(l in (k+1):n){
        sum.gamma.Theta=sum.gamma.Theta+sqrt(sum((Theta[k,]-Theta[l,])^2))
      }
    }
    gamma.Theta=(sum.gamma.Theta/(n*(n-1)/2))^(-2)
  }


  if(kerTheta=="Gaussian"){
    kernel.Theta=function(theta1,theta2){
      exp(-gamma.Theta*sum((theta1-theta2)^2))
    }
  }else if(kerTheta=="PIQ"){
    kernel.Theta=function(theta1,theta2){
      prod(1/(1+gamma.Theta*(theta1-theta2)^2))
    }
  }



  ### Compute kernel matrices

  K.U=matrix(NA,n,n)
  K.Theta=matrix(NA,n,n)
  for(i in 1:n){
    for(j in i:n){
      K.U[i,j]=kernel.U(U[i],U[j])
      K.U[j,i]=K.U[i,j]
      K.Theta[i,j]=kernel.Theta(Theta[i,],Theta[j,])
      K.Theta[j,i]=K.Theta[i,j]
    }
  }

  ### Compute centralized kernel matrix for \Theta


  norm.const=2*pi
  for(j in 1:(d-2)){
    norm.const=norm.const*integrate(
      function(x){(cos(x))^(d-1-j)},
      lower=-pi/2,
      upper=pi/2,
      abs.tol=1e-12)$value
  }

  if(kerTheta=="Gaussian"){
    E.kernel.Theta=rep(1,n)
    for(i in 1:n){
      for(j in 1:(d-2)){
        E.kernel.Theta[i]=E.kernel.Theta[i]*
          integrate(
            function(x){
              exp(-gamma.Theta*(Theta[i,j]-x)^2)*(cos(x))^(d-1-j)
            },
            lower=-pi/2,
            upper=pi/2,
            abs.tol=1e-12
          )$value
      }
      E.kernel.Theta[i]=E.kernel.Theta[i]*
        integrate(
          function(x){
            exp(-gamma.Theta*(Theta[i,d-1]-x)^2)
          },
          lower=-pi,
          upper=pi,
          abs.tol=1e-12
        )$value
      E.kernel.Theta[i]=E.kernel.Theta[i]/norm.const
    }


    EE.kernel.Theta=1
    if(d>2){
      for(j in 1:(d-2)){
        EE.kernel.Theta=EE.kernel.Theta*
          cubintegrate(
            function(x){
              exp(-gamma.Theta*(x[1]-x[2])^2)*(cos(x[1]))^(d-1-j)*(cos(x[2]))^(d-1-j)
            },
            lower=c(-pi/2,-pi/2),
            upper=c(pi/2,pi/2),
            absTol=1e-12
          )$integral
      }
    }
    EE.kernel.Theta=EE.kernel.Theta*
      cubintegrate(
        function(x){
          exp(-gamma.Theta*(x[1]-x[2])^2)
        },
        lower=c(-pi,-pi),
        upper=c(pi,pi),
        absTol=1e-12
      )$integral
    EE.kernel.Theta=EE.kernel.Theta/norm.const^2
  }else if(kerTheta=="PIQ"){
    E.kernel.Theta=rep(1,n)
    for(i in 1:n){
      for(j in 1:(d-2)){
        E.kernel.Theta[i]=E.kernel.Theta[i]*
          integrate(
            function(x){
              1/(1+gamma.Theta*(Theta[i,j]-x)^2)*(cos(x))^(d-1-j)
            },
            lower=-pi/2,
            upper=pi/2,
            abs.tol=1e-12
          )$value
      }
      E.kernel.Theta[i]=E.kernel.Theta[i]*
        integrate(
          function(x){
            1/(1+gamma.Theta*(Theta[i,d-1]-x)^2)
          },
          lower=-pi,
          upper=pi,
          abs.tol=1e-12
        )$value
      E.kernel.Theta[i]=E.kernel.Theta[i]/norm.const
    }


    EE.kernel.Theta=1
    if(d>2){
      for(j in 1:(d-2)){
        EE.kernel.Theta=EE.kernel.Theta*
          cubintegrate(
            function(x){
              1/(1+gamma.Theta*(x[1]-x[2])^2)*(cos(x[1]))^(d-1-j)*(cos(x[2]))^(d-1-j)
            },
            lower=c(-pi/2,-pi/2),
            upper=c(pi/2,pi/2),
            absTol=1e-12
          )$integral
      }
    }
    EE.kernel.Theta=EE.kernel.Theta*
      cubintegrate(
        function(x){
          1/(1+gamma.Theta*(x[1]-x[2])^2)
        },
        lower=c(-pi,-pi),
        upper=c(pi,pi),
        absTol=1e-12
      )$integral
    EE.kernel.Theta=EE.kernel.Theta/norm.const^2
  }



  K.tilde.Theta=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      K.tilde.Theta[i,j]=K.Theta[i,j]-E.kernel.Theta[i]-E.kernel.Theta[j]+EE.kernel.Theta
    }
  }


  one.n=matrix(1,n,1)
  test.stat.hs=as.numeric(1/n*t(one.n)%*%(K.U*K.tilde.Theta)%*%one.n)



  if(kerU=="Gaussian"){
    kernel.U.deriv=function(u1,u2){
      exp(-gamma.U*(u1-u2)^2)*2*gamma.U*(u1-u2)
    }
  }else if(kerU=="PIQ"){
    kernel.U.deriv=function(u1,u2){
      2*gamma.U*(u1-u2)/(1+gamma.U*(u1-u2)^2)^2
    }
  }

  if(kerTheta=="Gaussian"){
    kernel.Theta.deriv=function(k,theta1,theta2){
      exp(-gamma.Theta*sum((theta1-theta2)^2))*2*gamma.Theta*(theta1[k]-theta2[k])
    }
  }else if(kerTheta=="PIQ"){
    kernel.Theta.deriv=function(k,theta1,theta2){
      kernel.Theta(theta1,theta2)*
        2*gamma.Theta*(theta1[k]-theta2[k])/
        (1+gamma.Theta*(theta1[k]-theta2[k])^2)
    }
  }







  Psi.1=matrix(0,n*n,d)
  Psi.2=matrix(0,n*n,d*d)

  K.U.inv=solve(K.U+diag(eps,n))
  K.Theta.inv=solve(K.Theta+diag(eps,n))
  Id=diag(1,n)


  Sigma.inf.mat=Sigma.hat.sqrt%x%Sigma.hat+Sigma.hat%x%Sigma.hat.sqrt+diag(eps,d^2)
  Sigma.inf.mat.eig=eigen(Sigma.inf.mat)
  D.Sigma.inf.inv.mat=diag(1/Sigma.inf.mat.eig$values)
  P.Sigma.inf.mat=Sigma.inf.mat.eig$vectors

  Sigma.inf.inv=P.Sigma.inf.mat%*%D.Sigma.inf.inv.mat%*%t(P.Sigma.inf.mat)
  Sigma.inf.inv=Re(Sigma.inf.inv)


  for(j in 1:n){

    normlize.Xj=Sigma.hat.sqrt.inv%*%(X[j,]-mu.hat)
    normlize.Xj.Sigmainv=Sigma.hat.inv%*%(X[j,]-mu.hat)
    norm.sq.Xj=as.numeric(t(normlize.Xj)%*%normlize.Xj)
    norm.Xj=sqrt(norm.sq.Xj)
    U1.star.Xj=-t(normlize.Xj.Sigmainv)/norm.Xj
    U2.star.Xj=-(t(normlize.Xj.Sigmainv))%x%(t(normlize.Xj.Sigmainv))/(2*norm.Xj)
    V1.star.Xj=-(normlize.Xj%*%U1.star.Xj/norm.sq.Xj+Sigma.hat.sqrt.inv/norm.Xj)
    V2.star.Xj=-(normlize.Xj%*%U2.star.Xj/norm.sq.Xj+
                   ((t(X[j,]-mu.hat))%x%diag(1,d))%*%Sigma.inf.inv/norm.Xj)

    k.deriv.U.Uj=matrix(NA,n,1)
    for(i in 1:n){
      k.deriv.U.Uj[i,1]=kernel.U.deriv(U[i],U[j])
    }
    k.deriv.Theta.Thetaj=matrix(NA,n,d-1)
    for(i in 1:n){
      for(k in 1:(d-1)){
        k.deriv.Theta.Thetaj[i,k]=kernel.Theta.deriv(k,Theta[i,],Theta[j,])
      }
    }

    J.Xj=PolarDerivative(V[j,])


    Psi.1=Psi.1+
      (K.U.inv%*%k.deriv.U.Uj%*%U1.star.Xj)%x%(Id[,j]-one.n/n)+
      Id[,j]%x%(K.Theta.inv%*%k.deriv.Theta.Thetaj%*%J.Xj%*%V1.star.Xj)
    Psi.2=Psi.2+
      (K.U.inv%*%k.deriv.U.Uj%*%U2.star.Xj)%x%(Id[,j]-one.n/n)+
      Id[,j]%x%(K.Theta.inv%*%k.deriv.Theta.Thetaj%*%J.Xj%*%V2.star.Xj)
  }
  Psi.1=Psi.1/n
  Psi.2=Psi.2/n


  Gamma.hat=matrix(NA,n^2,n)
  for(i in 1:n){
    Gamma.hat[,i]=Id[,i]%x%(Id[,i]-one.n/n)+Psi.1%*%(X[i,]-mu.hat)+
      Psi.2%*%as.vector((X[i,]-mu.hat)%*%t(X[i,]-mu.hat)-Sigma.hat)
  }
  Gamma.hat=Gamma.hat/sqrt(n)

  Gamma.tilde.base=Gamma.hat
  K.U.sqrt=sqrtm(K.U+diag(eps,n))
  K.Theta.sqrt=sqrtm(K.Theta+diag(eps,n))
  for(i in 1:n){
    Gamma.tilde.base[,i]=matrix(K.Theta.sqrt%*%matrix(Gamma.hat[,i],n,n)%*%K.U.sqrt,n^2,1)
  }

  Gamma.tilde=t(Gamma.tilde.base)%*%Gamma.tilde.base

  lambda.hs=eigen(Gamma.tilde+diag(eps,n),only.values=TRUE)$values


  p.val.hs=max(imhof(test.stat.hs,lambda.hs)$Qq,0)

  list(stat=test.stat.hs,
       pval=p.val.hs,
       lambda=lambda.hs,
       gamma.U=gamma.U,
       gamma.Theta=gamma.Theta)
}










