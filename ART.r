# Sample code for ART, ART-A, and RTP
# https://www.biorxiv.org/content/10.1101/667238v1
#
library(MASS)
library(mvtnorm)

# RTP via Equ. 2
RTP <- function(Z, k, L) { # Z <- sum(-log(P[1:k]))
   integrate(function(x,y,m,n) 1-pgamma(log(qbeta(x,m+1,n-m))*m+y,m),0,1,Z,k,L)$va
}

ART.A <- function(P, k, L) { # "ART-A"
    wgt <- rep(1,k)
    z <- P
    z[1] <- ( 1 - P[1] )^L
    for(j in 2:L) z[j] <- ((1-P[j]) / (1-P[j-1]))^((L-(j-1)))
    p <- (1-z)[1:k]
    k = length(p)
    sumZ <- rep(0, k)
    y <- qnorm(p)
    z <- y
    gz <- z[1] * wgt[1]
    sumZ[1] <- gz
    for(i in 2:k) {
        gz <- p[ 1 : i ]
        for(j in 1:i) gz[j] <- z[j] * wgt[j]
        sumZ[i] <- sum(gz)
    }
    Lo = diag(k); Lo[lower.tri(Lo)] <- 1
    pSg <- Lo %*% diag(wgt[1:k]^2) %*% t(Lo)
    pCr <- cov2cor(pSg)
    sZ <- sumZ
    for(i in 1:k) {
        sZ[i] <- sumZ[i] / sqrt(diag(pSg))[i]
    }
    ppZ <- pmvnorm(lower = rep(-Inf,k), upper = rep(max(sZ), k), sigma = pCr)[1]
    c(ppZ, which.min(sZ))
}

ART <- function(lW, Pk, k, L) {
   d = (k-1)*(digamma(L+1) - digamma(k))
   ak = (k-1)*log(Pk) - lW + qgamma(1-pbeta(Pk, k, L-k+1), shape=d)
   1 - pgamma(ak, shape=k+d-1)
}

## read in summary statistics and LD matrix
d <- read.csv("Z_score.csv", header=TRUE, as.is=TRUE)
S <- read.csv("LD_matrix.csv", header=FALSE, as.is=TRUE)
## the number of statistics
Z <- d$Z
L <- length(Z)
## transform Z-scores to chi-squares
y <- Z^2
## calculate TQ-statistic
yy <- y %*% rep(1, L)
## eigen decomposition of the LD matrix
ee <- eigen(S); eivec <- ee$vectors; eigva <- ee$values
pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
## calculate decorreated statistics
x <- (Z %*% pc)^2
k <- round(0.5*L)
#k <- L
px <- 1-pchisq(x, df=1)
P <- sort(px)
P.rtp <- RTP(sum(-log(P[1:k])), k, L)
P.art <- ART(sum(log(P[1:(k-1)])), P[k], k, L)
P.arta <- ART.A(P, k, L)
cat(P.art, P.rtp, P.arta[1], P.arta[2], "\n")
