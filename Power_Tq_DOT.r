# This script is to check equivalency of theoretical and empirical power of Tq and DOT.
# Empirical power is obtained via simulating regressions of Y on each Xi where X is MVN 
# and Y ~ Beta*X + error. Errors have Laplace distribution with variance 1. 

library(MASS)
library(Matrix)
library(CompQuadForm)

Fu = function(x,l)  davies(x, lambda = l, lim=1e7, acc=1e-8)$Qq
Finv = function(p, l, br=c(1e-3, 10000)) {
  G = function(x) Fu(x,l) - p
  return( uniroot(G,br)$root )
}

Laplace <- function(n, sdev=1, b=sqrt(0.5)) { # R.V. from Laplace distribution
    u = runif(n, -0.5, 0.5)
    x = -b*sign(u)*log(1-2*abs(u))
    sdev*(x/sqrt(2 * b^2))
}

alpha <- 1e-7 # alpha level for power
Beta <- c(0, 0, 0, 0.3) # true regression coefficients
p <- length(Beta) # number of SNPs
N <- 500
sigma <- 1
R <- matrix(1, p, p)
#
if(p == 2) {
    p.sub <- p # number of proxi SNPs to be tested, p.sub <= p
    R[1,2] <- R[2,1] <- 0.5
} else {
    p.sub <- p - 1 # number of proxi SNPs to be tested, p.sub <= p
    R[lower.tri(R)] <- 2*rbeta(p*(p-1)/2, 0.9, 0.9) - 1
    R <- (R * lower.tri(R)) + t(R * lower.tri(R))
    diag(R) <- 1
    R <- nearPD(R, corr=TRUE, posd.tol = 1e-5)$mat
}
#
(R <- as.matrix(R))
solve(R)
#
Cors <- c((Beta %*% R) / sqrt(as.numeric(c(Beta) %*% R %*% Beta + sigma^2)))
Cors
x <- Cors
t.M <- R + R^2*(x %*% t(x))/2 - R*(x^2 %*% t(x^2)) - (x %*% t(x))/2; diag(t.M) <- 1
t.M
(t.M <- as.matrix(nearPD(t.M, corr=TRUE)$mat))
t.R <- R
for(i in 1:p) {
    for(j in 1:p) {
        xi <- x[i]
        xj <- x[j]
        Sij <- R[i,j]
        Vij <- 0.5 * (Sij^2*xi*xj +  (1-xi^2-xj^2)*(2*Sij - xi*xj)) / N
        t.R[i,j] <- N * Vij / ((1-xi^2)*(1-xj^2))
    }
}
t.R
(t.R <- as.matrix(nearPD(t.R, corr=TRUE)$mat))
Sgm <- t.R[1:p.sub, 1:p.sub]
Mf <- sqrt(N)
ee <- eigen(Sgm); eivec <- ee$vectors; eigva <- ee$values
sD <- diag(sqrt(1/eigva)); pc <- eivec %*% sD %*% t(eivec)
Nonc <- ( (Mf*Cors[1:p.sub]) %*% eivec %*% sD )^2
crit.val.m <- Finv(alpha, ee$values)
(pow.yy <- imhof(crit.val.m, lambda = eigva, delta = Nonc, epsrel = 1e-11, limit = 2e5)$Qq)
(pow.xx <- 1 - pchisq(qchisq(1-alpha, df=p.sub), df = p.sub, ncp = sum(Nonc)))

# Now do regression simulations and compare theoretical and empirical results
oo <- 5000
pows.YY <- pows.XX <- rep(0, oo)
Y <- vec.1 <- rep(1, p.sub)
j <- 1
for(j in 1:oo) {
    xg.sim <- mvrnorm(n = N, mu = rep(0, p), Sigma = R, emp=FALSE)
    Cv <- cov(xg.sim)
    Sgm.est <- cov2cor(Cv[1:p.sub, 1:p.sub])
    Ys <- (xg.sim %*% Beta) + Laplace(N) # rlogis(N, scale = sqrt(3)/pi) # rnorm(N)
    for(v in 1 : p.sub) {
        #Y[v] <- summary(lm(Ys ~ xg.sim[,v]))$coefficients[2,3] # <-- Slower way
        Y[v] <- sqrt(N)*cor(Ys, xg.sim[,v]) # <-- Faster way
    }
    Y2 <- Y^2
    YY <- Y2 %*% vec.1
    ee <- eigen(Sgm.est)
    eivec <- ee$vectors
    eigva <- ee$values
    sD <- diag(sqrt(1/eigva))
    pc <- eivec %*% sD %*% t(eivec)
    X2 <- (Y %*% pc)^2 %*% vec.1
    pc.cor <- pc
    pows.XX[j] <- 1 - pchisq(X2, df = p.sub)
    pows.YY[j] <- davies(YY, lambda = ee$values, lim=1e7, acc=1e-8)$Qq
    if(j %% 100 == 0) cat(".", j, ".", sep=""); if(j %% 1000 == 0) cat("\r")
}
cat("\n")
cat("Empirical powers, Tq =", ecdf(pows.YY)(alpha), "DOT =", ecdf(pows.XX)(alpha), "\n")
cat("Theoretical powers, Tq =", pow.yy, "DOT =", pow.xx, "\n")

### High correlation results in the induced effect sizes that are similar and Tq wins
# alpha <- 1e-7
# p.sub <- p <- 2
# N <- 500
# Beta <- c(0.3, 0)
# R[1,2] <- R[2,1] <- 0.99
# > Cors
# [1] 0.2873479 0.2844744
# Empirical powers, Tq = 0.8823 DOT = 0.8134 
# Theoretical powers, Tq = 0.8608337 DOT = 0.7969669 

### Low correlation results in the induced effect sizes that differ and DOT wins
# Note how at R=0.5 the two powers switch almost exactly
# alpha <- 1e-7
# p.sub <- p <- 2
# N <- 500
# Beta <- c(0.3, 0)
# R[1,2] <- R[2,1] <- 0.5
# > Cors
# [1] 0.2873479 0.1436739
# Empirical powers, Tq = 0.7236 DOT = 0.8088 
# Theoretical powers, Tq = 0.72062 DOT = 0.7972284 
