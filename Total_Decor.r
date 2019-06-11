# Sample code for the DOT:
# https://doi.org/10.1101/665133
#
library(MASS)
library(Matrix)
library(CompQuadForm)
library(mvtnorm)

## read in summary statistics and LD matrix
d <- read.csv("Z_score.csv", header=TRUE, as.is=TRUE)
S <- read.csv("LD_matrix.csv", header=FALSE, as.is=TRUE)
## the number of statistics
Z <- d$Z
n <- length(Z)
## transform Z-scores to chi-squares
y <- Z^2
## calculate TQ-statistic
yy <- y %*% rep(1, n)
## eigen decomposition of the LD matrix
ee <- eigen(S); eivec <- ee$vectors; eigva <- ee$values
pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
## calculate decorreated statistics
x <- (Z %*% pc)^2
## calculate DOT-statistic
xx <- x %*% rep(1, n)
## calculate DOT P-value
DOT <- 1 - pchisq(as.numeric(xx), df = n, ncp=0)
## calculate TQ O-value
TQ <- imhof(yy, lambda = eigva, delta = rep(0, n), epsrel = 1e-11, limit = 2e5)$Qq
## calculate average pair-wise LD
rho.av2 <- sqrt(mean(S[lower.tri(S)]^2))
## two way to approaximate TQ P-value
P.approx.1 <- 1-pchisq( (yy - (1 - rho.av2)*(n - 1) ) / (rho.av2*(n - 1) + 1), df=1)
Sgm.e <- (1-rho.av2)*diag(1,n,n) + rho.av2 * matrix(1,n,n)
ev <- eigen(Sgm.e)$values
P.approx.2 <- imhof(yy, lambda = eigen(Sgm.e)$values, delta = rep(0, n), epsrel = 1e-11, limit = 2e5)$Qq
cat("P-values by DOT and TQ", DOT, TQ, "\n")
cat("TQ Approximation by 1 degree of freedom chi-square:", P.approx.1, "\n")
cat("TQ Approximation by reducing LD to equicorrelation matrix:", P.approx.2, "\n")


#######################################################################
## Explore top results for TOX3
#######################################################################
## correlation between rs999737 (anchor SNP) and the rest of the SNPs in the block
c = S[,9]

## keep only relavent information
tmp.dat <- data.frame(rs = d$rs,
                      pos = d$pos,
                      cor = c,
                      Z = Z,
                      YY = y,
                      XX = as.vector(x))

## identify the highest values of the decorrelated statistics
tmp.dat[order(tmp.dat$XX),]
## 1  rs3095611 52483447  0.3356984 -0.8180895  0.66927051 3.4435309
## 8  rs8051542 52534167  0.6181170  3.4219796 11.70994441 4.3985414
## 4    rs43143 52522449 -0.4738636 -3.4036086 11.58455152 4.9062960 <-- 3rd
## 12 rs8046979 52555742 -0.3467780 -2.1149178  4.47287742 7.9250069 <-- 2nd
## 9  rs4784220 52535810  1.0000000  3.8715836 14.98915963 8.6790381 <-- 1st

## -----------------------------
## the highest (1sr) value of XX
## -----------------------------
tmp <- data.frame(wgts = pc[9,],
                  Z = Z,
                  cor = c,
                  pos = d$pos,
                  rs = d$rs)
tmp$comps <- pc[9,] * Z
## ## check whether we grabed the write weights
## (sum(tmp$wgts * tmp$Z))^2
## ## [1] 8.679038
## identify the sign of XX value prior to squaring
sum(tmp$wgts * tmp$Z)
## [1] 2.946021 <-- positive

tmp[order(tmp$comps),]                  #the highest positive values in the linear combination
## 6  -0.030608001 -2.7462054 -0.4010791 52528727 rs9926539  0.084055858
## 9   1.289591855  3.8715836  1.0000000 52535810 rs4784220  4.992762686 <--- ancor SNP, associated with breast cancer 


## -----------------------------
## the second highest value of XX
## -----------------------------
tmp <- data.frame(wgts = pc[12,],
                  Z = Z,
                  cor = c,
                  pos = d$pos,
                  rs = d$rs)
tmp$comps <- pc[12,] * Z
## ## check whether we grabed the write weights
## (sum(tmp$wgts * tmp$Z))^2
## ## [1] 7.925007
sum(tmp$wgts * tmp$Z)
## [1] -2.815139 <-- negative

tmp[order(tmp$comps),]                  #the highest positive values
##            wgts          Z        cor      pos        rs         comps
## 12  8.043394133 -2.1149178 -0.3467780 52555742 rs8046979 -17.011117652 <--- SNP with the highest contribution
## 4   0.149748810 -3.4036086 -0.4738636 52522449   rs43143  -0.509686339
## 7   0.044117485 -3.3588794 -0.5227874 52532355 rs1111481  -0.148185310
## 1   0.002893467 -0.8180895  0.3356984 52483447 rs3095611  -0.002367115


## -----------------------------
## the third highest value of XX
## -----------------------------
tmp <- data.frame(wgts = pc[4,],
                  Z = Z,
                  cor = c,
                  pos = d$pos,
                  rs = d$rs)
tmp$comps <- pc[4,] * Z
## ## check whether I grabbed the write weights
## (sum(tmp$wgts * tmp$Z))^2
## ## [1] 4.906296
sum(tmp$wgts * tmp$Z)
## [1] -2.215016  <-- negative

tmp[order(tmp$comps),]                  #the highest positive values
##           wgts           Z        cor      pos         rs         comps
## 4   4.12287930 -3.4036086 -0.4738636 52522449   rs43143 -14.03266745 <--- SNP with the highest contribution
## 11 -0.32945077  1.9551947  0.3571002 52543645 rs9302556  -0.64414041
