# generate Halton sequence

generate_Halton_pts <- function(n = 100, seeds = c(0,0), bases = c(2,3)) {

  d <- length(bases)
  pts <- c()
  for (i in 1:d) {
    b <- bases[i]
    u <- seeds[i]
    u_plus_k <- rep(u, n) + seq(0, n - 1)
    xk <- (u_plus_k %% b)/b;
    for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
      xk <- xk + (floor(u_plus_k/(b^j)) %% b)/(b^(j+1));
    }
    pts <- cbind(pts,xk)
  }
  return(pts)
}
