# This script is for checking equivalency of theoretical and empirical power of Tq and DOT.
# Empirical power is obtained via simulating regressions of Y on each Xi where X is discrete 
# (with values 1, 2, 3) with a given correlation matrix and Y ~ Beta*X + error.
library(MASS)
library(Matrix)
library(CompQuadForm)
library(GenOrd)
Fu = function(x,l)  davies(x, lambda = l, lim=1e7, acc=1e-8)$Qq
Finv = function(p, l, br=c(1e-3, 10000)) {
  G = function(x) Fu(x,l) - p
  return( uniroot(G,br)$root )
}
options(width=180)

ptol <- 1e-8
alpha <- 0.05
N <- 1000
p <- 5
s <- 1
R <- matrix(1, p, p)
repeat {
    R[lower.tri(R)] <- sort(2*rbeta(p*(p-1)/2, 0.1, 0.6) - 1) 
    R <- (R * lower.tri(R)) + t(R * lower.tri(R))
    diag(R) <- 1
    R <- as.matrix(nearPD(R, corr=TRUE, posd.tol = ptol)$mat)
    Beta <- runif(p, -0.1, 0.05) # random coefficients
    Cors <- c((Beta %*% R) / sqrt(as.numeric(c(Beta) %*% R %*% Beta + s^2)))
    x <- Cors; s.R <- t.R <- t.V <- R
    for(i in 1:p) {
        for(v in 1:p) {
            xi <- x[i]; xj = x[v]; Sij = R[i,v]
            t.V[i,v] <- 0.5*(Sij^2*xi*xj +  (1-xi^2-xj^2)*(2*Sij - xi*xj)) / N
            t.R[i,v] <- N * t.V[i,v] / ((1-xi^2)*(1-xj^2))
        }
    }
    t.R <- as.matrix(nearPD(t.R, corr=TRUE, posd.tol = ptol)$mat)
    f <- runif(p, 0.1, 0.9) # if z is SNP vector with values {0,1,2}, f = (1 - mean(z-1)) / 2 (Zaykin 2004)
    # Assume Hardy-Weinberg equilibrium
    ff <- list(); for(i in 1:p) ff[[i]] <- c(f[i]^2, 2*f[i]*(1-f[i]) + f[i]^2)
    Bah <- tryCatch(corrcheck(ff), error = function(e) NA)
    Rd <- tryCatch(contord(ff, t.R), error = function(e) NA)
    if(identical(NA, Bah) == FALSE && identical(NA, Rd) == FALSE ) {
        xg.sim <- tryCatch(ordsample(n=N, ff, Rd), error = function(e) NA)
        if(identical(NA, xg.sim) == FALSE) break
    }
}
Cors <- c((Beta %*% Rd) / sqrt(as.numeric(c(Beta) %*% Rd %*% Beta + s^2)))
ee <- eigen(Rd)
eivec <- ee$vectors
eigva <- ee$values
sD <- diag(sqrt(1/eigva)); pc <- eivec %*% sD %*% t(eivec)
Nonc <- ( (sqrt(N)*Cors) %*% eivec %*% sD )^2
crit.val.m <- Finv(alpha, ee$values)
(pwt.YY <- imhof(crit.val.m, lambda = eigva, delta = Nonc, epsrel = 1e-11, limit = 2e5)$Qq)
(pwt.XX <- 1 - pchisq(qchisq(1-alpha, df=p), df = p, ncp = sum(Nonc)))

oo <- 1000
St <- matrix(nrow=oo, ncol=p)
X2 <- YY <- pwe.YY <- pwe.XX <- rep(0, oo)
Y <- vec.1 <- rep(1, p)
for(j in 1:oo) {
    xg.sim <- scale(ordsample(n=N, ff, Rd))
    Ys <- (xg.sim %*% Beta) + s*rnorm(N)
    ee <- eigen(Rd); eivec <- ee$vectors; eigva <- ee$values
    sD <- diag(sqrt(1/eigva)); pc <- eivec %*% sD %*% t(eivec)
    Cors <- c((Beta %*% Rd) / sqrt(as.numeric(c(Beta) %*% Rd %*% Beta + s^2)))
    Nonc <- ( (sqrt(N)*Cors) %*% eivec %*% sD )^2
    crit.val.m <- Finv(alpha, ee$values)
    for(v in 1 : p) { Y[v] <- sqrt(N)*cor(Ys, xg.sim[,v]) }
    Sgm.est <- cor(xg.sim)
    ee <- eigen(Sgm.est)
    eivec <- ee$vectors
    eigva <- ee$values
    sD <- diag(sqrt(1/eigva))
    pc <- eivec %*% sD %*% t(eivec)
    St[j,] <- Y
    Y2 <- Y^2
    YY[j] <- Y2 %*% vec.1
    X2[j] <- (Y %*% pc)^2 %*% vec.1
    pwe.XX[j] <- 1 - pchisq(X2[j], df = p)
    pwe.YY[j] <- davies(YY[j], lambda = ee$values, lim=1e7, acc=1e-8)$Qq
    if(j %% 10 == 0) cat(".", j, ".", sep=""); if(j %% 100 == 0) cat("\r")
}
cat("\n")

( e.yy <- ecdf(pwe.YY)(alpha) ) ;
pwt.YY
( e.xx <- ecdf(pwe.XX)(alpha) ) ;
pwt.XX

# p = 10
# > ( e.yy <- ecdf(pwe.YY)(alpha) ) ; pwt.YY
# [1] 0.281
# [1] 0.2767971
# > ( e.xx <- ecdf(pwe.XX)(alpha) ) ; pwt.XX
# [1] 0.337
# [1] 0.3429447

# p = 10
# > ( e.yy <- ecdf(pwe.YY)(alpha) ) ; pwt.YY
# [1] 0.653
# [1] 0.6573614
# > ( e.xx <- ecdf(pwe.XX)(alpha) ) ; pwt.XX
# [1] 0.822
# [1] 0.805308


# p = 5
# > ( e.yy <- ecdf(pwe.YY)(alpha) ) ; pwt.YY
# [1] 0.775
# [1] 0.7729011
# > ( e.xx <- ecdf(pwe.XX)(alpha) ) ; pwt.XX
# [1] 0.831
# [1] 0.8333806

# p = 5
# > ( e.yy <- ecdf(pwe.YY)(alpha) ) ; pwt.YY
# [1] 0.932
# [1] 0.9376426
# > ( e.xx <- ecdf(pwe.XX)(alpha) ) ; pwt.XX
# [1] 0.886
# [1] 0.8964152

# p = 15
# > ( e.yy <- ecdf(pwe.YY)(alpha) ) ; pwt.YY
# [1] 0.78
# [1] 0.7673209
# > ( e.xx <- ecdf(pwe.XX)(alpha) ) ; pwt.XX
# [1] 0.93
# [1] 0.9212168


