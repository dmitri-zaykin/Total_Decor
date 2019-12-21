# This script is for power study of DOT and TQ

library(MASS)
library(Matrix)
library(CompQuadForm)


Fu = function(x,l)  davies(x, lambda = l, lim=1e6, acc=1e-7)$Qq
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
ha <- c(0, 0.5)
h0 <- rep( 0, 100-length(ha) )
Beta <- c(h0, ha) # true regression coefficients
p <- length(Beta) # number of SNPs
N <- 500
sigma <- 1
R <- matrix(1, p, p)
#
cnt <- 0
oo <- 100
sumCr <- matrix(0, nrow=oo, ncol=p)
tqpow <- dotpow <- dif <- minR <- maxR <- rep(0,oo)
for(k in 1:oo) {
    Beta <- sample(Beta)
    #p.sub <- p-1
    p.sub <- round(0.2*p)
    R[lower.tri(R)] <- 2*rbeta(p*(p-1)/2, 0.25, 0.25) - 1
    R <- (R * lower.tri(R)) + t(R * lower.tri(R))
    diag(R) <- 1
    R <- nearPD(R, corr=TRUE, posd.tol = 1e-8)$mat
    R <- as.matrix(R)
    #rho <- 0.5
    #S <- (1-rho)*diag(1,p) + rho * matrix(1,p,p)
    #dlt <- 5
    #v <- u <- runif(p, -dlt, dlt)
    #B <- S + u %*% t(v)
    #B <- cov2cor(B)
    #B[lower.tri(B)] <- B[lower.tri(B)] + (2*rbeta(p*(p-1)/2, 0.9, 0.9) - 1)
    #B <- (B * lower.tri(B)) + t(B * lower.tri(B))
    #diag(B) <- 1
    #R <- cov2cor(B)
    Cors <- c((Beta %*% R) / sqrt(as.numeric(c(Beta) %*% R %*% Beta + sigma^2)))
    sumCr[k,] <- x <- Cors
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
    t.R <- as.matrix(nearPD(t.R, corr=TRUE)$mat)
    Sgm <- t.R[1:p.sub, 1:p.sub]
    Mf <- sqrt(N)
    ee <- eigen(Sgm); eivec <- ee$vectors; eigva <- ee$values
    sD <- diag(sqrt(1/eigva)); pc <- eivec %*% sD %*% t(eivec)
    Nonc <- ( (Mf*Cors[1:p.sub]) %*% eivec %*% sD )^2
    crit.val.m <- Finv(alpha, ee$values)
    pow.yy <- imhof(crit.val.m, lambda = eigva, delta = Nonc, epsrel = 1e-11, limit = 2e5)$Qq
    pow.xx <- 1 - pchisq(qchisq(1-alpha, df=p.sub), df = p.sub, ncp = sum(Nonc))
    if(pow.xx > pow.yy) cnt <- cnt+1
    #cat(pow.yy, pow.xx, "\n")
    dotpow[k] <- pow.xx
    tqpow[k] <- pow.yy
    dif[k] <- abs(pow.yy-pow.xx)
    cr <- R[lower.tri(R)]
    minR[k] <- min(cr)
    maxR[k] <- max(cr)
}
cnt/oo
summary(R[lower.tri(R)])
mean(tqpow)
mean(dotpow)
mean(dif)
mean(minR); mean(maxR)
Cr <- apply(sumCr, 2, sum)
sd(Cr); sd(Cr)/mean(Cr)
#summary(Cr)
