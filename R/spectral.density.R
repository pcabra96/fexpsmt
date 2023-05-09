spectral.density <- function(eta = c(), p = 0, q = 0, d_long=0, x){
  
  if(p > 0) {
    phi <- eta[1:(p)]
    px <- outer(x, 1:p)
    Rar <- cos(px) %*% phi
    Iar <- sin(px) %*% phi
    
    far <- (1-Rar)^2 + Iar^2
  } else {
    phi <- numeric(0)
    far <- 1
  }
  
  if(q > 0) {
    psi <- eta[(p+1):(p+q)]
    px <- outer(x, 1:q)
    Rma <- cos(px) %*% psi
    Ima <- sin(px) %*% psi
    
    fma <- (1+Rma)^2 + Ima^2
  } else {
    psi <- numeric(0)
    fma <- 1
  }
  if(d_long>0){
    d <- eta[p+q+1]
    lm = sapply(x, function(y) abs(2*sin(y*(1/2)))^(-2*d))
  }
  else{
    psi = 0
    lm <- 1
  }
  
  spec <- lm*(fma/far)
  
  spectrum = (1/(2*pi))*spec
  
  return(spectrum)
}