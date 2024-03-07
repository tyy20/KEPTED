#' Product-type Inverse-Quadratic (PIQ) kernel computation
#'
#' Computing the values of Product-type Inverse-Quadratic (PIQ) kernel functions.
#'
#' @param gamma A number, the bandwidth parameter in the PIQ kernel.
#' @param z1 A vector, the first input of the PIQ kernel.
#' @param z2 A vector, the second input of the PIQ kernel.
#'
#' @details
#' The Product-type Inverse-Quadratic (PIQ) kernel is defined as
#' k(z1,z2)=Prod_j(1/(1+gamma*(z1_j-z2_j)^2)).
#'
#' @return A number, the value of the PIQ kernel function.
#'
#'
#' @examples
#' gamma=0.02
#' z1=c(3,1,3)
#' z2=c(8,1,9)
#' kerPIQ(gamma,z1,z2)
#'
#' @export


kerPIQ=function(gamma,z1,z2){
  prod(1/(1+gamma*(z1-z2)^2))
}
