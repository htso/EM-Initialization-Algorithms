# loop-runs-script.R
# Aug 26, 2012

library(Hext)
library(InitExp)

data(simdat2single)
dat.nm ="Single2"
simdat = simdat2.1
x = list(x=as.matrix(simdat), TT=TT)
class(x) = "Hext.data"
mtype = rep(5, NCOL(simdat))
p = c(0, 0, 0, 0, length(mtype))
D = 0

maxit = 500
tol = 1e-02
var.reduct = 0.5
spherical = TRUE
wdf=2

Nruns = 3
M = M.truth
K = K.truth
MM = c(M, M+1, M+2)
KK = c(K, K+1, K+2)

system.time(ans1 <- loop.over.hext.plot(fname="-test", init.method="random", x=x, mtype=mtype, 
                                        K=KK, K.truth=K.truth, M=MM, D=D, 
                                        maxit=maxit, verbose=FALSE))

system.time(ans2 <- loop.over.hext(init.method="unifseg", x=x, mtype=mtype, K=KK, M=MM, D=D, 
                                   maxit=maxit, verbose=TRUE))

ans3 = loop.over.hext(init.method="maxsep", x=x, mtype=mtype, K=KK, M=MM, D=D, maxit=maxit, verbose=TRUE)
ans4 = loop.over.hext(init.method="LL",     x=x, mtype=mtype, K=KK, M=MM, D=D, maxit=maxit, verbose=TRUE)
ans5 = loop.over.hext(init.method="segkm",  x=x, mtype=mtype, K=KK, M=MM, D=D, maxit=maxit, verbose=TRUE)
ans6 = loop.over.hext(init.method="PCA",    x=x, mtype=mtype, K=KK, M=MM, D=D, maxit=maxit, verbose=TRUE)







