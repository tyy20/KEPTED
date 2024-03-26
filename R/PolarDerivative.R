#' Derivative of the polar coordinate transformation
#'
#' This function compute the Jacobian matrix of the polar transformation
#' theta=g(v), i.e., the transformation from the the rectangular coordinate
#' representation of the directional vector to its angular representation.
#'
#' @param v A \code{d}-dimensional directional vector of length 1.
#'
#' @details
#' See Lemma 3 of Tang and Li (2024).
#'
#' @return The Jacobian matrix of the polar transformation theta=g(v), with
#' \code{d-1} rows and \code{d} columns.
#'
#' @references
#' \cite{Tang, Y. and Li, B. (2024), “A nonparametric test for elliptical
#' distribution based on kernel embedding of probabilities,”
#' \url{https://arxiv.org/abs/2306.10594}}
#'
#' @examples
#' X=c(3,1,3)
#' V=X/sqrt(sum(X^2))
#' PolarDerivative(V)
#'
#' @export

PolarDerivative=function(v){
  d=length(v)
  S.sq=rep(NA,d)
  S.sq[d]=v[d]^2
  for(i in (d-1):1){
    S.sq[i]=S.sq[i+1]+v[i]^2
  }
  S=sqrt(S.sq)
  g.deriv=matrix(0,d-1,d)
  for(i in 1:(d-2)){
    g.deriv[i,i]=S[i+1]/S.sq[i]
    for(k in (i+1):d){
      g.deriv[i,k]=-v[i]*v[k]/(S.sq[i]*S[i+1])
    }
  }
  g.deriv[d-1,d-1]=v[d]/S.sq[d-1]
  g.deriv[d-1,d]=-v[d-1]/S.sq[d-1]
  g.deriv
}
