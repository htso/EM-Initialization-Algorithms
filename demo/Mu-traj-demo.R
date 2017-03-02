library(Hext)
library(InitExp)

data(simdat2fullcov)
dat.nm ="Fullcov"
simdat = simdat2.3
x = list(x=as.matrix(simdat), TT=TT)
class(x) = "Hext.data"
mtype = rep(5, NCOL(simdat))
p = c(0, 0, 0, 0, length(mtype))

Nruns = 5
K = 3
M = 5
maxit = 20
tol = 1e-01
var.reduct=0.5
spherical=TRUE
wdf=4
D = 0

system.time(mu.paths <- mu.trajectory(Nruns=Nruns, init.method = "random", x=x, mtype=mtype, 
                         K=K, M=M, D=D, maxit=maxit, tol=tol, 
                         wdf=wdf, var.reduct=var.reduct, spherical=spherical,
                         verbose=FALSE))

X11();par(mar=c(1,1,4,1))
plot(simdat, type="p", col="grey", xlab="", ylab="",
     main=paste("EM estimates of Gaussian Means vs True Means"))
points(mu.paths$mu[[1]], col="blue", pch=18, cex=2)
points(mu.paths$mu[[2]], col="green", pch=17, cex=2)
points(mu.paths$mu[[3]], col="brown", pch=20, cex=2)
points(mu.paths$mu[[4]], col="grey", pch=16, cex=2)
points(mu.paths$mu[[5]], col="pink", pch=15, cex=2)
points(mu.paths$mu[[6]], col="black", pch=14, cex=2)
ctrs = M5.truth$mu.ll
for ( k in 1:K ) {
  x.m1 = ctrs[[k]][[1]][1] 
  y.m1 = ctrs[[k]][[1]][2]
  x.m2 = ctrs[[k]][[2]][1] 
  y.m2 = ctrs[[k]][[2]][2]
  points(x=x.m1, y=y.m1, pch=19, cex=3, col="red")
  points(x=x.m2, y=y.m2, pch=19, cex=3, col="red")
}

X11();hist(mu.paths$ans[,3], breaks=30, xlim=c(-9600, -8200), xlab="",
               main=paste("Loglikelihood Distribution for", dat.nm))




