#' Rectangular to polar coordinates
#'
#' Given a \code{d}-dimensional vector \code{X} in rectangular coordinate, this
#' function compute its polar coordinate \code{(R,Theta)}, where \code{R} is the
#' length of \code{X} and the \code{(d-1)}-dimensional vector \code{Theta}
#' contains the \code{d-1} angles of \code{X}.
#'
#' @param X A vector in rectangular coordinate. Suppose the dimension of
#' \code{X} is \code{d}.
#'
#' @details
#' The formula corresponds to theta=g(v) as in Lemma 1 of Tang and Li (2024).
#'
#' @return A list of the following:
#' \item{R}{The length of \code{X}.}
#' \item{Theta}{A vector of length \code{d-1}, containing the angles of \code{X}.}
#'
#' @references
#' \cite{Tang, Y. and Li, B. (2024), “A nonparametric test for elliptical
#' distribution based on kernel embedding of probabilities,”
#' \url{https://arxiv.org/abs/2306.10594}}
#'
#' @examples
#' X=c(3,1,3)
#' Rec2Polar(X)
#'
#' @export

Rec2Polar=function(X){
  d=length(X)
  r=sqrt(sum(X^2))
  theta=rep(NA,d-1)
  for(i in 1:(d-2)){
    theta[i]=atan(X[i]/sqrt(sum(X[(i+1):d]^2)))
  }
  theta[d-1]=atan2(X[d-1],X[d])
  list(R=r,Theta=theta)
}
