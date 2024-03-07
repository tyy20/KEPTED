#' Gaussian kernel computation
#'
#' Computing the values of Gaussian kernel functions.
#'
#' @param gamma A number, the bandwidth parameter in the Gaussian kernel.
#' @param z1 A vector, the first input of the Gaussian kernel.
#' @param z2 A vector, the second input of the Gaussian kernel.
#'
#' @details
#' The Gaussian kernel is defined as k(z1,z2)=exp(-gamma*||z1-z2||^2).
#'
#' @return A number, the value of the Gaussian kernel function.
#'
#'
#' @examples
#' gamma=0.02
#' z1=c(3,1,3)
#' z2=c(8,1,9)
#' kerGauss(gamma,z1,z2)
#'
#' @export


kerGauss=function(gamma,z1,z2){
  exp(-gamma*sum((z1-z2)^2))
}
