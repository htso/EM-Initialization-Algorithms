#' Records the position of the means of each Gaussian mixture state for clustering display
#'
#' @param Nruns number of runs
#' @param init.method initialization method
#' @param x data structure for Hext models
#' @param mtype vector to indicate the type of variables in obs 
#' @param maxit maximum number of iterations, default = 1000
#' @param tol tolerance, default = 1e-14
#' @param K number of hidden states
#' @param M number of mixtures for each state
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param wdf degree of freedom for Wishart
#' @param var.reduct var.reduct
#' @param fudge fudge factor
#' @param spherical boolean
#' @param verbose boolean 
#' @param TT length of observation sequence, if vector, it has the length of each sequence in x
#' @description see example 
#' @return list of three elements : mu.init, the initial position of the gaussian means, mu.final, the final position,
#'            ans, the result of one call to single.run(..)
#' @export
mu.trajectory = function(Nruns, init.method, x, mtype, K, M, D, maxit=1000, tol=1e-04, 
                         wdf=1, var.reduct=0.5, fudge=0.1, verbose=TRUE, spherical=TRUE) {
  p = NCOL(x$x)
  ans = matrix(NA, nrow=Nruns, ncol=4)
  colnames(ans) = c("K", "M", paste(init.method, c("LL", "Niter"), sep="."))
  mu.final = rep(list(list()), K*M)
  mu.init = rep(list(list()), K*M)
  for ( i in 1:(K*M)) {
    mu.final[[i]] = matrix(NA, ncol=p, nrow=Nruns )
    mu.init[[i]] = matrix(NA, ncol=p, nrow=Nruns )
  }
  ans[,1] = K
  ans[,2] = M
  for ( i in 1:Nruns ) {
    cat("run #:", i, " by init method:", init.method, " about to start.\n")
    res = single.run(x=x, mtype=mtype, K=K, M=M, D=D, 
                     init.method=init.method, wdf=wdf, tol=tol, maxit=maxit, 
                     var.reduct=var.reduct, fudge=fudge, verbose=verbose)
    ans[i,3] = res$LL
    ans[i,4] = res$last.it
    for ( k in 1:K ) {
      for ( m in 1:M ) {
        mu.final[[m + (k-1)*M]][i,] = res$final$model$parms.emission$M5$mu.ll[[k]][[m]]
        mu.init[[m + (k-1)*M]][i,] = res$init$b$M5$mu.ll[[k]][[m]]
      }
    }
    cat("done.\n")
  }
  list(ans=ans, mu=mu.final, mu.init=mu.init)
}

#' Plot the trajectories of the gaussian means over each EM iterations until convergence
#'
#' @param dat data set
#' @param n.it number of iterations, to be use to loop over mu
#' @param mu list with the same structure as model$parms.emission$M5$mu.ll, see details
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param LLtrace LLtrace
#' @description paramter mu is a list where each element holds the paramter of one iteration, and each element 
#'              has the gaussian means of each mixture component in each state.
#' @return None
#' @export
mu.trace = function(dat, n.it, K, M, LLtrace, mu){
  mycol = c("red", "blue", "brown","green",  "pink", "black")
  windows()
  plot(dat, type="p", col="grey", main=paste("LL:", round(LLtrace[n.it],2)))
  for ( k in 1:K ) {
    for ( i in 2:n.it ) {
      x0.m1 = ctrs[[i-1]][[k]][[1]][1] 
      y0.m1 = ctrs[[i-1]][[k]][[1]][2]
      x1.m1 = ctrs[[i]][[k]][[1]][1]
      y1.m1 = ctrs[[i]][[k]][[1]][2]
      x0.m2 = ctrs[[i-1]][[k]][[2]][1] 
      y0.m2 = ctrs[[i-1]][[k]][[2]][2]
      x1.m2 = ctrs[[i]][[k]][[2]][1]
      y1.m2 = ctrs[[i]][[k]][[2]][2]
      lines(x=c(x0.m1, x1.m1), y=c(y0.m1, y1.m1), col=mycol[k], lwd=2)
      lines(x=c(x0.m2, x1.m2), y=c(y0.m2, y1.m2), col=mycol[k], lwd=2)
      if ( i == 2 ) {
        points(x=x0.m1, y=y0.m1, pch=1, col=mycol[6])
        points(x=x0.m2, y=y0.m2, pch=1, col=mycol[6])
      }
      if ( i == Nit ) {
        points(x=x1.m1, y=y1.m1, pch=19, col=mycol[k])
        points(x=x1.m2, y=y1.m2, pch=19, col=mycol[k])
      }
    }
  }
}

#' Measure the decoding performance of a hidden Markov model
#'
#' @param ss.fit predicted state sequence from a decoding algorithm, e.g. Viterbi
#' @param ss true sequence of hidden states
#' @param digit number of digit, default 3
#' @description If the decoding is perfect, then the decoded state sequence should match with the true state
#'          sequence exactly. If some of the observations should be in state i, but is mapped to 
#'          state j, then that is a decoding error. 
#' @return a matrix of size K x K.out, where K is the number of true states in the obs sequence, K.out is the number of states
#'         decoded by the decoding algorithm. It could be interpreted as follow. If the decoding is perfect,
#'          then every row and every column would have *exactly* one non-zero element. This means every actual
#'          state is mapped to *exactly* one decoded state. However, if a row has multiple entries, 
#'          that means the HMM thinks some of state i look different and assign them to more than
#'          one state. If a column has multiple entries, that means the variability of 
#'          certain state is so great that some of its observations are classified into separate
#'          decoded states.
#' @export
match.confusion = function(ss, ss.fit, digit=3) {
  n = length(ss)
  n.fit = length(ss.fit)
  if ( n != n.fit ) {
    error("sequence length mismatch.")
  }
  K = length(unique(ss))
  KK = max(ss.fit)
  confuse = matrix(0, nrow=K, ncol=KK)
  # Here I assume ss has states numbered from 1 to K. No gap in between.
  for ( i in 1:K ) {
    ix = which(ss==i)
    tbl = table(ss.fit[ix])
    ss.no = as.numeric(names(tbl))
    for ( j in 1:length(ss.no)) {
      confuse[i, ss.no[j]] = tbl[j] / sum(tbl)
    }
  }
  round(confuse, digit)
}
#' Plot the Viterbi states sequence by sequence
#'
#' @param SS state sequences
#' @param TT vector of length of each sequence
#' @description plot the decoded state sequence
#' @return None
#' @export
Vit.SS.plot = function(SS, TT) {
  windows()
  n = length(TT)
  Tcum = c(0, cumsum(TT))
  if ( n < 9 ) {
    nr = 2
    nc = ceiling(n/2)
  } else {
    nr = 3
    nc = ceiling(n/3)
  }
  par(mfrow=c(nr, nc))
  for ( i in 1:n ) {
    plot(SS[(Tcum[i]+1):Tcum[i+1]], type="s", xlab="Time", ylab="State", main=paste("Seq:", i))
  }
}

#' Plot the loglikelihood trace of the fitted HMM
#'
#' @param mod HMM model
#' @param n.trunc n.trunc
#' @description plot the LL 
#' @return None
#' @export
LL.plot = function(mod, n.trunc=5) {
  X11()
  n = mod$last.it
  plot(mod$loglik1[n.trunc:n], type='l', col="red", ylab="Log Likelihood", xlab="iteration")
  lines(mod$loglik2[n.trunc:n], type='l', col="blue")
  title(paste("Best LL:", round(mod$loglik1[n],2), "  N.iter:", n))
}

#' One pass through the HMM fit-and-decode cycle
#'
#' @param mod HMM model
#' @param x data structure for Hext models
#' @param mtype vector to indicate the type of variables in obs 
#' @param TT vector of the lengths of observation sequences
#' @param K number of states for the run
#' @param M number of mixtures for the run
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param init.method initialization method
#' @param wdf degree of freedom for Wishart
#' @param tol tolerance
#' @param var.reduct variance reduction factor
#' @param fudge fudge factor
#' @param verbose boolean
#' @param spherical spherical
#' @param maxit maximum number of iterations
#' @description This function takes the observation sequences, initialize the HMM parameters, call hextfit to 
#'              fit a model, then decode the hidden state sequence with the Viterbi algorithm.
#' @return list of the following items :
#'          sseq : the decoded state sequences,
#'          LL : loglikelihood of the HMM model,
#'          AIC, AICc, BIC : the AIC, AICC, and BIC of the model,
#'          init : the initial value of the paramters generated by the initialization method,
#'          res : the fitted model,
#'          last.it : the state of last iteration.
#' @export
single.run = function(x, mtype, K, M, D, init.method="random", wdf=1, tol=1e-01, maxit=1000, 
                      var.reduct=1, fudge=0.1, spherical=TRUE, verbose=TRUE) {
  p = c(0, 0, 0, 0, length(mtype))
  D = 0
  obs = x$x
  if ( init.method == "km") {
    init =  init.GM.by.kmeans(obs, K=K, M=M, var.reduct=var.reduct, verbose=verbose)
  } else if (init.method == "maxsep") {
    init = init.GM.max.sep(obs, K=K, M=M, var.reduct=var.reduct, verbose=verbose)
  } else if (init.method == "unifseg") {
    init = init.GM.unif.seg(x, K=K, M=M, var.reduct=var.reduct, verbose=verbose)
  } else if ( init.method == "random") {
    init = init.GM.random(obs, K=K, M=M, var.reduct=var.reduct, verbose=verbose)
  } else if (init.method == "LL")  {
    init = init.GM.by.LL(x=x, K=K, M=M, D=D, mtype=mtype, Nset=30, 
                         var.reduct=var.reduct, verbose=verbose)
  } else if (init.method == "PCA")  {
    init = init.GM.by.PCA(x, K=K, M=M, A=NULL, Pi=NULL, 
                          spherical=spherical, 
                          var.reduct=var.reduct, 
                          fudge=fudge,
                          verbose=verbose)
  } else if (init.method == "segkm")  {    
    # reason for this try() wrapper is Seg,KMeans may fail and we need to 
    # keep calling it until function returns valid result. 
    res = try(init.GM.Seg.KMeans(x=x, K=K, M=M, D=D, mtype=mtype, Pi=NULL, A=NULL, var.reduct=var.reduct, verbose=verbose), TRUE)
    while(inherits(res, "try-error")) {
      res = try(init.GM.Seg.KMeans(x=x, K=K, M=M, D=D, mtype=mtype, Pi=NULL, A=NULL, var.reduct=var.reduct, verbose=verbose), TRUE)
    }
    init = res
  } else {
    stop("initialization method not recognized.")
  }
  mod = hmmspec.gen(K=K, D=D, M=M, mtype=mtype,
                    linit=init$Pi, ltrans=init$A, parms.emission=init$b, 
                    mstep=mstep.generic, 
                    dens.emission=generic.density, 
                    Log=TRUE, Intv=FALSE, wdf=wdf)
  res = hextfit(x, mod, tol=tol, maxit=maxit, verbose=verbose)
  pred = Viterbi.hext(res$model, x)
  sseq = pred$s
  LL = res$loglik1[res$last.it]
  AIC = res$Aic
  AICc = res$Aicc
  BIC = res$Bic
  return(list(sseq=sseq, LL=LL, AIC=AIC, AICc=AICc, BIC=BIC, init=init, final=res, 
              last.it=res$last.it))
}

#' Plotting function for 2-dim models
#'
#' @param mod HMM model
#' @param x data matrix for Hext models, where each row is an observation at time t, columns are features
#' @param mtype vector to indicate the type of variables in obs 
#' @param TT vector of the lengths of observation sequences
#' @param K vector of length Nruns, each component is the number of states for the ith run
#' @param M vector of length Nruns, each component is the number of mixtures for the ith run
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param init.method initialization method
#' @param wdf degree of freedom of the Wishart
#' @param tol tolerance
#' @param var.reduct variance reduction factor
#' @param fudge fudge factor
#' @param verbose boolean
#' @param spherical boolean 
#' @param maxit maximum number of iterations for each run
#' @description Plot function to write the plot to a PDF file. See loop.over.hext(..) for details.
#' @return None
#' @export
loop.over.hext.plot = function(init.method, x, mtype, K, K.truth, M, D, 
                               wdf=1, tol=1e-01, maxit=1000, var.reduct=0.5, fudge=0.1,
                               spherical=TRUE, verbose=TRUE, fname=NULL) {
  require(ellipse)
  if ( !is.null(fname) ) {
    pdf(paste("model-init-by-", init.method, "-", "K", K.truth, fname, ".pdf", sep=""))  
  }
  Nruns = length(K)
  IC = matrix(NA, ncol=5, nrow=Nruns)
  for ( i in 1:Nruns ) {
    cat("K", K[i], " -- M", M[i], "\n")
    res = single.run(x=x, mtype=mtype, K=K[i], M=M[i], D=D, var.reduct=var.reduct, fudge=fudge,
                     init.method=init.method, wdf=wdf, tol=tol, maxit=maxit, verbose=verbose )
    IC[i,1] = res$AIC
    IC[i,2] = res$AICc
    IC[i,3] = res$BIC
    IC[i,4] = -res$LL
    IC[i,5] = res$last.it-1
    # the following ploting codes work best if the data is 2-dimensional, otherwise
    # it only plots the first and second dimension
    if ( is.null(fname)) X11()
    par(mfrow=c(2, ceiling(K[i]/2)), mar=c(1,1,2,1))
    for ( j in 1:K[i] ) {
      plot(x$x[,1:2], type="p", col="grey", main=paste("Iter:", i, " State ", j, " init means & cov"))
      for ( m in 1:M[i] ) {
        lines(ellipse(x=res$init$b$M5$sigma.ll[[j]][[m]], 
                      centre=res$init$b$M5$mu.ll[[j]][[m]], which=c(1,2)), col=j, lwd=2)
      }
    }
    if ( is.null(fname)) X11()
    par(mfrow=c(2, ceiling(K[i]/2)), mar=c(1,1,2,1))
    for ( j in 1:K[i] ) {
      plot(x$x[,1:2], type="p", col="grey", main=paste("Iter:", i, " State ", j, " final means & cov"))
      for ( m in 1:M[i] ) {
        lines(ellipse(x=res$final$model$parms.emission$M5$sigma.ll[[j]][[m]], 
                      centre=res$final$model$parms.emission$M5$mu.ll[[j]][[m]], which=c(1,2)), col=j, lwd=2)
      }
    }
    if ( is.null(fname)) X11()
    par(mfrow=c(1,1), mar=c(1,1,2,1))
    plot(res$sseq, type="s", main=paste("Iter:", i, " State Sequence (AIC:", round(res$AIC,0), " BIC:", round(res$BIC,0), ")"))
  }
  if ( is.null(fname)) X11()
  par(mfrow=c(2,2), mar=c(2,2,2,2))
  plot(K, IC[,1], type="b", main="AIC vs K", col="red")
  plot(K, IC[,2], type="b", main="AICc vs K", col="pink")
  plot(K, IC[,3], type="b", main="BIC vs K", col="blue")
  plot(K, IC[,4], type="b", main="negative LL vs K", col="black")
  if (!is.null(fname)) dev.off()
  return(IC)
}
#' Loop over the single.run function
#'
#' @param init.method initialization algorithm to use
#' @param x data matrix for Hext models, where each row is an observation at time t, columns are features
#' @param mtype vector to indicate the type of variables in obs 
#' @param K vector of length Nruns, each component is the number of states for the ith run
#' @param M vector of length Nruns, each component is the number of mixtures for the ith run
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param wdf degree of freedom for Wishart
#' @param tol tolerance
#' @param var.reduct variance reduction factor
#' @param fudge fudge factor
#' @param verbose print out status of the run, boolean
#' @param spherical spherical
#' @param maxit maximum number of iterations for each run
#' @description See single.run(...) for details
#' @return list of Lik, K, M, sigmal.ll, mu.ll
#' @export
loop.over.hext = function(init.method, x, mtype, K, M, D, maxit=1000, tol=1e-04, 
                          wdf=1, var.reduct=0.5, fudge=0.1, verbose=TRUE, spherical=TRUE) {
  Nruns = length(K)
  Lik = matrix(NA, nrow=Nruns, ncol=4)
  colnames(Lik) = paste(init.method, c("AIC", "AICc", "BIC", "LL"), sep=".")
  sigma.ll = list()
  mu.ll = list()
  for ( i in 1:Nruns ) {
    cat("run #:", i, " init method:", init.method, "\n")
    res = single.run(x=x, mtype=mtype, K=K[i], M=M[i], D=D, 
                     init.method=init.method, wdf=wdf, tol=tol, maxit=maxit, 
                     fudge=fudge, var.reduct=var.reduct, spherical=spherical, verbose=verbose)
    sigma.ll[[i]] = res$final$model$parms.emission$M5$sigma.ll
    mu.ll[[i]] = res$final$model$parms.emission$M5$mu.ll
    Lik[i,1] = res$AIC
    Lik[i,2] = res$AICc
    Lik[i,3] = res$BIC
    Lik[i,4] = res$LL
    cat("\n")
  }
  return(list(Lik=Lik, K=K, M=M, sigma.ll=sigma.ll, mu.ll=mu.ll))
}




