#' Polar to rectangular coordinates
#'
#' Given a polar coordinate representation \code{(R,Theta)} of a
#' \code{d}-dimensional vector \code{X}, where \code{R} is the length of
#' \code{X} and the \code{(d-1)}-dimensional vector \code{Theta} contains the
#' \code{d-1} angles of \code{X}, this function compute \code{X} in its
#' rectangular coordinate representation.
#'
#' @param R The length of \code{X}.
#' @param Theta A vector of length \code{d-1}, containing the angles of \code{X}.
#'
#' @details
#' The formula corresponds to v=rho(theta) as in Lemma 1 of Tang and Li (2024).
#' See also Anderson (2003).
#' Note that when \code{d=2}, \code{V} will be \code{(sin(Theta),cos(Theta))}.
#'
#' @return A list of the following:
#' \item{X}{A vector in rectangular coordinate.}
#' \item{V}{The directional vector of \code{X}. Note that \code{V} is always on
#' the unit sphere.}
#'
#' @references
#' \cite{Tang, Y. and Li, B. (2024), “A nonparametric test for elliptical
#' distribution based on kernel embedding of probabilities,”
#' \url{https://arxiv.org/abs/2306.10594}}
#' \cite{Anderson, T. W. (2003). An Introduction to Multivariate Statistical
#' Analysis. John Wiley & Suns, Inc. Huboken, New Jersey.}
#'
#' @examples
#' R=2
#' Theta=c(pi/6,pi/3)
#' Polar2Rec(R,Theta)
#'
#' @export

Polar2Rec=function(R,Theta){
  d=length(Theta)+1
  if(d==2){
    v=c(sin(Theta),cos(Theta))
  }else{
    v=rep(0,d)
    v[1]=sin(Theta[1])
    for(i in 2:(d-1)){
      v[i]=prod(cos(Theta[1:(i-1)]))*sin(Theta[i])
    }
    v[d]=prod(cos(Theta[1:(d-1)]))
  }
  x=R*v
  list(X=x,V=v)
}
