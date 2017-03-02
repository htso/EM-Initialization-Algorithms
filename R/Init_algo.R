#' Non-informative common start method
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @description see my paper for detailed description
#' @details  Reference : HTK Toolkit, Ch. 8, HINT module
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.Noninform = function(x, K, M, Pi=NULL, A=NULL, var.reduct=1, dir.a=1, spherical=TRUE, verbose=TRUE) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("Non-informative method.\n")    
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  obs = x$x
  p = NCOL(obs)
  mu.global = apply(obs, 2, mean)
  var.global = var(obs)
  cov.m = diag(p)
  if (spherical==TRUE) {
    diag(cov.m) = diag(var.global) * var.reduct
  } else {
    cov.m = var.global
    diag(cov.m) = diag(cov.m) * var.reduct
  }
  
  for ( i in 1:K ) {
    for ( j in 1:M ) {
      cov1.m = cov.m
      diag(cov1.m) = diag(cov1.m)*0.05
      b$M5$mu.ll[[i]][[j]] = mu.global + rmvnorm(1, mean=rep(0, p), sigma=cov1.m)
      b$M5$sigma.ll[[i]][[j]] = cov.m  
    }
    b$M5$c.l[[i]] = rep(1/M, M)
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  if ( is.null(Pi) ) {
    dir.parm = rep(dir.a, K)
    Pi = log(as.vector(rdirichlet(1, dir.parm)))
  }
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  if ( is.null(A) ) {
    dir.parm = rep(dir.a, K)
    A = log(rdirichlet(K, dir.parm))
  }
  if ( verbose == TRUE)
    cat("DONE.\n")
  list(Pi=Pi, A=A, b=b)
}

#' Uniform common start method
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @description see my paper for detailed description
#' @details  Reference : HTK Toolkit, Ch. 8, HINT module
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.unifcommon = function(x, K, M, Pi=NULL, A=NULL, var.reduct=1, dir.a=1, spherical=TRUE, verbose=TRUE) {
  if ( verbose==TRUE)  
    cat("Uniform common start method.\n")
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  dirich = rep(1, M)
  
  obs = x$x
  p = NCOL(obs)
  TT = x$TT
  # Since TT is a vector of varying lengths for each sequence, 
  # i need to divide each length into K equal sized portion.
  len = ceiling( TT / K )  # len is a vector
  seg.ix = NULL
  # create a linear index into obs, which consists of multiple sequences
  for ( i in 1:length(len) ) {
    tmp = rep(1:K, each=len[i])
    tmp2 = tmp[1:TT[i]]
    seg.ix = c(seg.ix, tmp2)
  }
  
  for ( i in 1:K ) {
    seg = obs[which(seg.ix==i),]
    cov.m = diag(p)
    mu.global = apply(seg, 2, mean)
    var.global = var(seg)
    if (spherical==TRUE) {
      diag(cov.m) = diag(var.global) * var.reduct
    } else {
      cov.m = var.global
      diag(cov.m) = diag(cov.m) * var.reduct
    }
    
    for ( j in 1:M ) {
      b$M5$mu.ll[[i]][[j]] = mu.global 
      b$M5$sigma.ll[[i]][[j]] = cov.m  
    }
    b$M5$c.l[[i]] = rep(1/M, M)
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  if ( is.null(Pi) ) {
    dir.parm = rep(dir.a, K)
    Pi = log(as.vector(rdirichlet(1, dir.parm)))
  }
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  if ( is.null(A) ) {
    dir.parm = rep(dir.a, K)
    A = log(rdirichlet(K, dir.parm))
  }
  if ( verobse == TRUE)
    cat("DONE.\n")
  list(Pi=Pi, A=A, b=b)
}


#' Segmental k-means algorithm
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param mtype vector to indicate the type of variables in obs 
#' @description see my paper for detailed description
#' @details  Reference : Juang & Rabiner, Mixture Autoregressive HMM for Speech Signals, IEEE
#'           Transactions on Acoustics, Speech, and Signal Processing, Vol. ASSP-33, No 6, Dec 1985.
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.Seg.KMeans = function(x, K, M, D, mtype, Pi=NULL, A=NULL, 
                              tol=1e-2, maxit=20, dir.a=3, var.reduct=1, verbose=FALSE) {
  if ( verbose == TRUE)
    cat("Segmental k-means method.\n")
  obs = x$x
  p = NCOL(obs)
  TT = x$TT
  Ttot = sum(TT)
  if ( verbose == TRUE)
    cat("....enter infinite loop.\n")
  while (TRUE) {
    init = init.GM.random(obs=obs, K=K, M=M, dir.a=dir.a, var.reduct=var.reduct, verbose=verbose)
    mod = hmmspec.gen(K=K, D=D, M=M, mtype=mtype,
                      linit=init$Pi, ltrans=init$A, parms.emission=init$b,
                      Log=TRUE, Intv=FALSE)
    pred = Viterbi.hext(mod, x)
    ss = pred$s
    if ( all(table(ss) >= 3) & length(unique(ss)) == K )
      break
  }
  if ( verbose == TRUE)
    cat("....exit infinite loop.\n")
  LL = pred$loglik
  
  LLlast = 0
  LLbest = -1e300
  init.best = NULL
  it = 1
  while( abs(LL - LLlast) > tol & it <= maxit ) {
    for ( i in 1:K ) {
      ix = which(ss == i)
      if ( length(ix) > 3 ) {
        mx = kmeans(obs[ix,], M)
        ctr = mx$centers
        mxi = mx$cluster
        for ( j in 1:M ) {
          init$b$M5$mu.ll[[i]][[j]] = ctr[j,] 
          mat.tmp = var(obs[mxi==j,])
          init$b$M5$sigma.ll[[i]][[j]] = mat.tmp
          frac = sum(mxi==j)/length(mxi)
          init$b$M5$c.l[[i]][[j]] = frac
        }
      } else {
        stop("Segmental k-means failed. Some states are sparse.")
      }
    }
    mod = hmmspec.gen(K=K, D=D, M=M, mtype=mtype,
                      linit=init$Pi, ltrans=init$A, parms.emission=init$b, 
                      Log=TRUE, Intv=FALSE)
    pred = Viterbi.hext(mod, x)
    ss = pred$s
    LLlast = LL
    LL = pred$loglik
    if ( LL > LLbest ) {
      LLbest = LL
      init.best = init
    }
    it = it + 1
    if (verbose == TRUE) {
      cat("it:", it, " - ", round(LL,2), "\n")
    }
  }
  if ( verbose == TRUE)
    cat("DONE.\n")
  return(init.best)
}

#' Random initialization for Pi, A, and b
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param mtype vector to indicate the type of variables in obs 
#' @description The means of the gaussian mixture are generated by uniform sampling 
#'               from the observed range in each dimension. The variances for each
#'               spherical covariance matrix is picked uniformly from half the 
#'               min variance to the maximum variance. 
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.random = function(obs, K, M, var.reduct=1, dir.a=1, verbose=FALSE) {
  require(MCMCpack)
  if ( verbose == TRUE)
    cat("Random initialization method.\n")
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  
  p = NCOL(obs)
  
  # rough estimate of the min and max of variance-covariances from observations  
  VCOV = var(obs)
  Var.act = diag(VCOV)
  Var.max = max(Var.act)
  Var.min = min(Var.act)
  mean.vec = apply(obs, 2, mean)
  max.rng = apply(obs, 2, max)
  min.rng = apply(obs, 2, min)
  # generate M*K uniformly distributed random points inside the data range
  n = M * K
  pts = matrix(NA, ncol=p, nrow=n)
  # OLD METHOD : simulate points in a p-dimensional box bounded by empirical ranges,
  # this does not take into account shape of the data.
  # NOTE : need a function to simulate from multivariate empirical distribution.
  #for ( i in 1:p ) {
  #  pts[,i] = runif(n, max=max.rng[i], min=min.rng[i])
  #}
  # generate a variance for each of the M*K random points
  # where each variance is uniformly sampled from half the
  # minimum observed variance to the maximum obs variance.
  # var.random = matrix(runif(n, max=Var.max, min=Var.min*0.5), ncol=K)
  #alternatively, just set to the average variance
  # var.ave = mean(Var.act)
  
  # NEW METHOD : simulate from a multivariate guassian with means equal the empirical means
  # and sigma equal the empirical covariance matrix.
  # Note that the underlying distribution is assumed to be a mixture of gaussian,
  # so this simulation is only approximately correct
  pts = rmvnorm(n, mean=mean.vec, sigma=VCOV)
  
  dir.parm = rep(3, M)
  for ( i in 1:K ) {
    for ( j in 1:M ) {
      # This assignment is OK since these M*K points are random
      b$M5$mu.ll[[i]][[j]] = pts[i+(j-1)*K,]   
      # assume spherical covariance for each point
      # >>>>> CHOICE 1 : 
      # b$M5$sigma.ll[[i]][[j]] = var.reduct*var.ave*diag(p)
      # >>>>> CHOICE 2 : 
      # b$M5$sigma.ll[[i]][[j]] = var.random[j,i]*diag(p)  
      # >>>>> CHOICE 3 :
      # RHS below : scalar * vector * vector = vector
      cov.tmp = var.reduct*Var.act*diag(p)
      if ( is.positive.definite(cov.tmp) == TRUE ) {
        b$M5$sigma.ll[[i]][[j]] = cov.tmp  
      } else {
        mat.tmp = make.positive.definite(cov.tmp)
        b$M5$sigma.ll[[i]][[j]] = mat.tmp
      }
    }
    # use dirichlet sampling for the mixture proportions
    b$M5$c.l[[i]] = as.vector(rdirichlet(1, dir.parm))
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  dir.parm = rep(dir.a, K)
  Pi = log(as.vector(rdirichlet(1, dir.parm)))
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  dir.parm = rep(dir.a, K)
  A = log(rdirichlet(K, dir.parm))
  if ( verbose == TRUE)
    cat("DONE.\n")
  list(Pi=Pi, A=A, b=b)
}

#' Best Loglikelihood algorithm
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param mtype vector to indicate the type of variables in obs 
#' @param Nset number of random sets from which to choose the best one
#' @param dir.a dir.a
#' @param var.reduct variance reduction factor
#' @description Choose initial values by picking the best set from Nset sets of random 
#'         values. This function calls init.GM.random to generate Nset random sets.
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.by.LL = function(x, K, M, D, mtype, Nset=10, dir.a=1, var.reduct=1, verbose=TRUE) {
  if (verbose == TRUE) 
    cat("Best Likelihood method.\n")    
  obs = x$x
  TT = x$TT
  startval = rep(list(list()), Nset)
  LL = rep(NA, Nset)
  for ( i in 1:Nset ) {
    # Step 1 : generate N sets of starting values for {pi, A, b}
    tmp0 = init.GM.random(obs, K=K, M=M, dir.a=dir.a, var.reduct=var.reduct, verbose=verbose)
    mod = hmmspec.gen(K=K, D=D, M=M, 
                      linit=tmp0$Pi, ltrans=tmp0$A, parms.emission=tmp0$b, 
                      mtype=mtype, dens.emission=generic.density, mstep=mstep.generic, 
                      Log=TRUE, Intv=FALSE, wdf=1)
    startval[[i]] = tmp0
    # Step 2 : run viterbi on each
    viterb = Viterbi.hext(mod, x)
    LL[i] = viterb$loglik
  }
  # Step 3 : choose the set with the highest likelihood
  i.best = which.max(LL)
  if (verbose == TRUE) 
    cat("DONE.\n")    
  return(startval[[i.best]])
}

#' Initialization by PCA
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param mtype vector to indicate the type of variables in obs 
#' @param Nset number of random sets from which to choose the best one
#' @param dir.a dir.a
#' @param var.reduct variance reduction factor
#' @description This algorithm uses principal component analysis to discern local maxima 
#' in the data. It divides the data into K equal sized segments, then send
#' each segment into Cluster.by.PCA to find M clustering points in the data
#' space, which becomes the mu's of the gaussian mixtures in M5. The 
#' covariance matrices of the components are calculated from the PCA-segmented
#' data whenever possible, and are made spherical if the segmented data is
#' insufficient to come up with a positive definite covariance matrix. If the 
#' number of points returned by Cluster.by.PCA is less than M, then the deficit
#' is supplied with random vectors centered around the PCA-clustered points. 
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.by.PCA = function(x, K, M, Pi=NULL, A=NULL, spherical=TRUE, var.reduct=1, fudge=0.1, dir.a=1, verbose=FALSE) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("PCA Initialization method.\n")    
  require(corpcor)
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  
  obs = x$x
  TT = x$TT
  p = NCOL(obs)
  n = NROW(obs)
  
  if ( p <= 3 ) {
    m.star = p
  } else {
    m.star = ceiling(log(p)) # m.star should scale like log(.)
  }
  
  dirich = rep(1, M)
  
  # Since TT is a vector of varying lengths for each sequence, 
  # i need to divide each length into K equal sized portion.
  len = ceiling( TT / K )  # len is a vector
  seg.ix = NULL
  # create a linear index into obs, which consists of multiple sequences
  for ( i in 1:length(len) ) {
    tmp = rep(1:K, each=len[i])
    tmp2 = tmp[1:TT[i]]
    seg.ix = c(seg.ix, tmp2)
  }
  
  for ( i in 1:K ) {
    seg = obs[which(seg.ix==i),]
    seg.var = diag(var(seg))
    ave.var = mean(seg.var)
    # Use the version that doesn't require ks package for Ubuntu
    Cpts = Cluster.by.PCA1(seg, m.star=m.star, K=M, 0.1, verbose=verbose)
    n.ctrs = NROW(Cpts) # n.ctrs should be exactly M or less.
    if ( verbose == TRUE ) {
      cat("seg:",i," size:", NROW(seg), "\t")
      cat("PCA found:", n.ctrs, " centers.\t")
      cat("seg var:", seg.var, "\n")
    }
    if ( n.ctrs == 1) {
      copies = make.noisy.copies(Cpts, n=(M-1), var.s=ave.var*fudge)
      Cpts = rbind(Cpts, copies)
    } else if ( (n.ctrs > 1) & (n.ctrs < M) ) {
      ix1 = sample(n.ctrs, size = M - n.ctrs, replace=TRUE) # replace needs to be set to TRUE here
      cc = table(ix1)
      ix.set = as.numeric(names(cc))
      copies = NULL
      for ( ii in 1:length(ix.set) ) {
        tmp = make.noisy.copies(Cpts[ix.set[ii],], n=cc[ii], var.s=ave.var*fudge )
        copies = rbind(copies, tmp)
      }
      Cpts = rbind(Cpts, copies)
    }
    if (verbose==TRUE) {
      cat("no of vectors in Cpts:", NROW(Cpts), "\n")
    }
    # just pick the first M points from Cpts
    icls = data2cluster(Cpts, seg)
    cnts = table(icls)
    if ( verbose == TRUE ) {
      cat("members in each cluster:", cnts, "\n")
    }
    ix2 = which(cnts > 5) # must have at least 5 pts for cov matrix calc
    M.adj = length(ix2)
   
    if ( spherical == FALSE ) {  
      for ( j in 1:M.adj ) {
        b$M5$mu.ll[[i]][[j]] = Cpts[j,] 
        ii = as.numeric(names(ix2[j]))
        tmp = seg[icls==ii,]
        cov.m = var(tmp)
        if (verbose == TRUE) cat("cov.m:", cov.m, "\n")
        diag(cov.m) = diag(cov.m)*var.reduct
        if ( is.positive.definite(cov.m) == TRUE ) {
          b$M5$sigma.ll[[i]][[j]] = cov.m  
        } else {
          mat.tmp = make.positive.definite(cov.m)
          b$M5$sigma.ll[[i]][[j]] = mat.tmp
        }
      }
      if ( M.adj < M ) {
        cov.m = matrix(0, ncol=p, nrow=p)
        diag(cov.m) = seg.var * var.reduct
        if (verbose == TRUE) cat("diag(cov.m):", diag(cov.m), "\n")
        for ( j in (M.adj+1):M ) {
          b$M5$mu.ll[[i]][[j]] = Cpts[j,]
          if ( is.positive.definite(cov.m) == TRUE ) {
            b$M5$sigma.ll[[i]][[j]] = cov.m
          } else {
            mat.tmp = make.positive.definite(cov.m)
            b$M5$sigma.ll[[i]][[j]] = mat.tmp
          }
        }
      }
    } else {
      cov.m = matrix(0, ncol=p, nrow=p)
      # use a spherical covariance matrix for each mixture component,
      for ( j in 1:M ) {
        b$M5$mu.ll[[i]][[j]] = Cpts[j,]
        diag(cov.m) = seg.var * var.reduct
        if ( is.positive.definite(cov.m) == TRUE ) {
          b$M5$sigma.ll[[i]][[j]] = cov.m
        } else {
          mat.tmp = make.positive.definite(cov.m)
          b$M5$sigma.ll[[i]][[j]] = mat.tmp
        }
      }
    }
    frac = rep(0, M)
    frac[] = 1/M
    b$M5$c.l[[i]] = frac
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  if ( is.null(Pi) ) {
    dir.parm = rep(dir.a, K)
    Pi = log(as.vector(rdirichlet(1, dir.parm)))
  }
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  if ( is.null(A) ) {
    dir.parm = rep(dir.a, K)
    A = log(rdirichlet(K, dir.parm))
  }
  if (verbose == TRUE) 
    cat("DONE.\n")    
  list(Pi=Pi, A=A, b=b)
}

#' Uniform segmentation method
#'
#' @param x data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param D vector of the size of the alphabets in each column of the multivariate discrete obs matrix
#' @param mtype vector to indicate the type of variables in obs 
#' @param Nset number of random sets from which to choose the best one
#' @param dir.a dir.a
#' @param var.reduct variance reduction factor
#' @description This uniform segmentation method is described in the HTK Toolkit (Section 8.2, p.133)
#' and makes sense only for left-right models. The codes below check whether a 
#' covariance matrix is positive definite, if not it computes the nearest p.d. matrix
#' using the algorithm in NJ Higham (1988, Linear Algebra Appl. 103:103-118).
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.unif.seg = function(x, K, M, Pi=NULL, A=NULL, var.reduct=1, dir.a=1, verbose=TRUE) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("Uniform segmentation method.\n")    
  require(corpcor)
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  
  dirich = rep(1, M)
  
  obs = x$x
  TT = x$TT
  # Since TT is a vector of varying lengths for each sequence, 
  # i need to divide each length into K equal sized portion.
  len = ceiling( TT / K )  # len is a vector
  seg.ix = NULL
  # create a linear index into obs, which consists of multiple sequences
  for ( i in 1:length(len) ) {
     tmp = rep(1:K, each=len[i])
     tmp2 = tmp[1:TT[i]]
     seg.ix = c(seg.ix, tmp2)
  }

  for ( i in 1:K ) {
    seg = obs[which(seg.ix==i),]
    km = kmeans(seg, M)
    mu = km$centers
    vq = km$cluster
    frac = rep(0, M)
    for ( j in 1:M ) {
      b$M5$mu.ll[[i]][[j]] = mu[j,] 
      cov.m = var(seg[vq==j,])
      diag(cov.m) = diag(cov.m) * var.reduct
      if ( is.positive.definite(cov.m) == TRUE ) {
        b$M5$sigma.ll[[i]][[j]] = cov.m  
      } else {
        mat.tmp = make.positive.definite(cov.m)
        b$M5$sigma.ll[[i]][[j]] = mat.tmp
      }
      frac[j] = sum(vq==j) / NROW(seg)
    }
    b$M5$c.l[[i]] = frac
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  if ( is.null(Pi) ) {
    dir.parm = rep(dir.a, K)
    Pi = log(as.vector(rdirichlet(1, dir.parm)))
  }
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  if ( is.null(A) ) {
    dir.parm = rep(dir.a, K)
    A = log(rdirichlet(K, dir.parm))
  }
  if (verbose == TRUE) 
    cat("DONE.\n")    
  list(Pi=Pi, A=A, b=b)
}

#' Initialization by k-means method
#'
#' @param obs data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param dir.a dir.a
#' @param var.reduct variance reduction factor
#' @description This methods runs k-means twice, first to segment the entire data set into K groups,
#'       then segment each group into M subgroups corresponding to the mixture components.
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.by.kmeans = function(obs, K, M, Pi=NULL, A=NULL, dir.a=1, var.reduct=1, verbose=TRUE) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("Initialization by k-means.\n")    
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  
  VQ = Kmeans.Vector.Quantization(obs, K)
  p = NCOL(obs)
  dirich = rep(1, M)
  for ( i in 1:K ) {
    ix = which(VQ$codebk == i)
    # if there are more than 5 data points in this cluster, then we're OK
    if ( length(ix) > 5 ) {
      tmp = obs[ix,]
      mx = kmeans(tmp, M)
      ctr = mx$centers
      mxi = mx$cluster
      for ( j in 1:M ) {
        b$M5$mu.ll[[i]][[j]] = ctr[j,] 
        b$M5$sigma.ll[[i]][[j]] = var(tmp[mxi==j,])
      }
    } else if ( length(ix) > 2 ) {
      # In case we have fewer than 5 data points in this cluster,
      # use ad hoc method to estimate the mean and 
      # covariances; requires further thought.....
      for ( j in 1:M ) {
        mu0 = mean(obs[ix,])
        # perturb the cluster mean by a little noise
        mu1 = mu0 * runif(1, min=0.9, max=1.1)
        var0 = var(obs[ix,])
        # perturb the cluster variance by a little noise
        var1 = var0 * var.reduct
        b$M5$mu.ll[[i]][[j]] = rep(mu1, p)
        b$M5$sigma.ll[[i]][[j]] = var1*diag(p)
      }
    } else {
      stop("k-means fails.")
    }
    trial = as.vector(rdirichlet(1, dirich))
    while( any(trial < 0.1) ) {
      trial = as.vector(rdirichlet(1, dirich))
    }
    b$M5$c.l[[i]] = trial
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  if ( is.null(Pi) ) {
    dir.parm = rep(dir.a, K)
    Pi = log(as.vector(rdirichlet(1, dir.parm)))
  }
  
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  if ( is.null(A) ) {
    dir.parm = rep(dir.a, K)
    A = log(rdirichlet(K, dir.parm))
  }
  if (verbose == TRUE) 
    cat("DONE.\n")    
  return(list(Pi=Pi, A=A, b=b))
}

#' Initialization by Max Sep Algorithm
#'
#' @param obs data used in hextfit(..), must have components : x$x and x$TT
#' @param K number of hidden state
#' @param M number of gaussian mixtures
#' @param dir.a dir.a
#' @param var.reduct variance reduction factor
#' @description This is the maximal separation algorithm described in my paper (cite)
#'       which finds the points farthest away from each other. It calls max.separate.set
#' @details Note that M5 has the same structure as in gen.init.       
#' @return list of Pi, initial probability of the hidden states, A, the transition matrix, and b, the emission matrix
#' @export
init.GM.max.sep = function( obs, K, M, Pi=NULL, A=NULL, dir.a=1, var.reduct=1, verbose=TRUE ) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("Max Sep method.\n")    
  require(fields)
  require(MCMCpack)
  b = list(M1=list(prob.l=list(), labels=NULL),
           M2=list(mu.l=list(), sigma.l=list()),
           M3=list(shape.l=list(), scale.l=list()),
           M4=list(shape1.l=list(), shape2.l=list()),
           M5=list(c.l=list(), mu.ll=rep(list(list()),K), sigma.ll=rep(list(list()),K)))
  
  p = NCOL(obs)
  n = M * p * 100
  
  # rough estimate of the min and max of variance-covariances from observation data matrix
  Vmat = var(obs)
  Var.act = diag(Vmat)
  Var.max = max(Var.act)
  Var.min = min(Var.act)
  mean.vec = apply(obs, 2, mean)
  max.vec = apply(obs, 2, max)
  min.vec = apply(obs, 2, min)
  
  len = ceiling(NROW(obs) / K)
  ix = rep(1:K, each=len) 
  ix = ix[1:NROW(obs)]
  
  for ( i in 1:K ) {
    # c_jm :
    # make sure no component weight is too small
    dir.parm = rep(1, M)
    cjm = as.vector(rdirichlet(1, dir.parm))
    while( any(cjm < 0.1) ) {
      cjm = as.vector(rdirichlet(1, dir.parm))
    }
    b$M5$c.l[[i]] = cjm
    # mu_jm :
    # before i can think of better way to cover the 0-dim space, just 
    # generate random vector according to mean.vec
    seg = obs[which(ix==i),]
    Smax = max.separate.N(seg, M)
    for ( m in 1:M ) {
      b$M5$mu.ll[[i]][[m]] = obs[Smax[m],] #test this code !!!!!!!
    }
    # sigma_jm :
    # simple stuff : spherical with same radius for all mixture states
    cov.m = matrix(0, ncol=p, nrow=p)
    diag(cov.m) = Var.act * var.reduct
    for ( m in 1:M ) {
      b$M5$sigma.ll[[i]][[m]] = cov.m
    }
  }
  
  # initial prob generated from Dirichlet distribution with equal alpha
  if ( is.null(Pi) ) {
    dir.parm = rep(dir.a, K)
    Pi = log(as.vector(rdirichlet(1, dir.parm)))
  }
  # initialize each row of the transition prob with Dirichlet distribution 
  # where alpha are set 1.0 
  if ( is.null(A) ) {
    dir.parm = rep(dir.a, K)
    A = log(rdirichlet(K, dir.parm))
  }
  if (verbose == TRUE) 
    cat("DONE.\n")    
  return(list(Pi=Pi, A=A, b=b))
}




#' Clustering by PCA, version 2
#'
#' @param dat n x p matrix
#' @param K the number of cluster points to return
#' @param m.star the number of principal components to consider, should be no bigger than 5
#' @param cutoff.frac the level on the kernel density below which peaks are not considered
#' @param verbose boolean
#' @description similar to Cluster.by.PCA(..)
#' @return matrix
#' @export
Cluster.by.PCA1 = function (dat, m.star, K, cutoff.frac=0.1, verbose=FALSE) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("PCA Clustering method (version 2).\n")
  n = NROW(dat)
  p = NCOL(dat)
  # number of principal components should not exceed dimension of the data
  if ( m.star > p ) { m.star = p  }
  # prcomp is preferred over princomp because it uses SVD, 
  # which has better numeric accuracy ?
  # cat("dat size:", dim(dat), "\n")
  pca = prcomp(dat, retx=TRUE, center=TRUE, scale=FALSE)
  # pick the first m.star principal vectors
  PC = as.matrix(pca$rotation[,1:m.star]) # as.matrix is needed in case m.star=1
  proj = matrix(NA, nrow=n, ncol=m.star)
  fhat = rep(list(list()), m.star)
  pk.loc = rep(list(list()), m.star)
  if ( verbose == TRUE ) {
    nn = ceiling(sqrt(m.star))
    X11(); par(mfrow=c(nn,nn))  
  }
  
  # Go each principal component direction searching for local maxima
  for ( i in 1:m.star) {
    proj[,i] = as.vector(dat %*% PC[,i])
    # use the build-in kernel density estimation function
    fhat[[i]] = density(x=proj[,i], bw="sj")
    
    if( verbose == TRUE  ) {
      plot(fhat[[i]]$x, fhat[[i]]$y, type="l", main=paste("m.star:", i))
    }
    # find all local maxima over a neighborhood 
    # of 5 points (5 is deemed sufficient)
    pk = PK.Identify(fhat[[i]]$y, 5) 
    ix = which(pk == 1)
    # unlikely that ix is length zero?? but...
    if (length(ix) == 0 ) { stop("Couldn't find any peak in kernel density profile.") }
    pk.levels = fhat[[i]]$y[ix]
    pk.highest = max(pk.levels)
    if (verbose == TRUE ) {
      cat("pk.levels:", pk.levels, "\n")  
    }
    
    # filter out all the small peaks, which is defined as %cutoff.frac of the
    # highest peak in the kernel density
    ix.sub = which( pk.levels > cutoff.frac * pk.highest)
    if (verbose == TRUE ) {
      cat("significant peaks (ix.sub):", ix.sub, "\n")  
    }
    # note that ix.sub is at least length 1.

    # assign a vector to the pk.loc list
    pk.loc[[i]] = fhat[[i]]$x[ix[ix.sub]]
  }
  #pt0 = pk.loc[1] * PC[,i] + pk.loc[1] * PC2
  n.vec = sapply(pk.loc, function(ll)length(ll))
  #cat("n.vec:", n.vec, "\n")
  # create the product set :
  n.ll = rep(list(list()), m.star)
  for ( j in 1:m.star ) {
    n.ll[[j]] = 1:n.vec[j]
  }
  product.set = expand.grid(n.ll) # no. rows = # permutation, no. cols = m.star 
  
  N.pts = NROW(product.set)
  if ( verbose == TRUE ) {
    cat("no in product set:", N.pts, "\n")
  }
  pk.mat = matrix(NA, ncol=m.star, nrow=N.pts)
  for ( i in 1:N.pts ) {
    for ( j in 1:m.star ) {
      pk.mat[i,j] = pk.loc[[j]][product.set[i,j]]
    }
  }
  # pk.mat : N.pts x m.star
  # PC : p x m.star
  PP = pk.mat %*% t(PC)
  # result : (N.pts x m.star) * t(p x m.star) = N.pts x p
  #cat("PP:", NROW(PP), "\n")
  
  # Calculate the distance from each point in PP to all the data points in
  # the original coordinate. 
  # CAUTION : This matrix could be huge!!
  D = matrix(NA, ncol=N.pts, nrow=n)
  # D : each column corresponds to a point in PP and each row is the distance of
  # that cluster point to all the actual data points.
  for ( i in 1:N.pts ) {
    D[,i] = dist.one(PP[i,], dat) 
  }
  # For each cluster point, find the 10th percentile, which is the distance
  # needed to go out to encounter 1/10 of the data points. This is a better
  # measure of the proximity to data than the above, which could be distorted
  # by a few outliers.
  qtl.vec = apply(D, 2, quantile, probs=0.1) # pick 10%, with curse of dimensionality consideration 
  if ( verbose == TRUE ) {
    cat("qtl.vec:", qtl.vec, "\n")
  }
  if ( N.pts > K ) {
    # pick the cluster points which are closest to some data
    tmp = sort(qtl.vec, decreasing=FALSE)[1:K]
    ix = which(qtl.vec %in% tmp )
    candidates = PP[ix,]
  } else {
    # if there are less than K points, just return all.
    candidates = PP
  }
  if (verbose == TRUE) 
    cat("DONE.\n")    
  return(candidates)
}

#' Clustering by PCA, version 1
#'
#' @param dat n x p matrix
#' @param K the number of cluster points to return
#' @param m.star the number of principal components to consider, should be no bigger than 5
#' @param cutoff.frac the level on the kernel density below which peaks are not considered
#' @param verbose boolean
#' @description similar to Cluster.by.PCA(..)
#' @return candidate matrix
#' @export
Cluster.by.PCA = function (dat, m.star=4, K, cutoff.frac=0.1, verbose=FALSE) {
  require(MCMCpack)
  if (verbose == TRUE) 
    cat("PCA Clustering method (version1 ).\n")    
  require(ks)
  require(Fundifun)
  n = NROW(dat)
  p = NCOL(dat)
  # number of principal components should not exceed dimension of the data
  if ( m.star > p ) { m.star = p  }
  # prcomp is preferred over princomp because it uses SVD, 
  # which has better numeric accuracy ?
  # cat("dat size:", dim(dat), "\n")
  pca = prcomp(dat, retx=TRUE, center=TRUE, scale=FALSE)
  # pick the first m.star principal vectors
  
  PC = as.matrix(pca$rotation[,1:m.star]) # as.matrix is needed in case m.star=1
  proj = matrix(NA, nrow=n, ncol=m.star)
  fhat = rep(list(list()), m.star)
  pk.loc = rep(list(list()), m.star)
  nn = ceiling(sqrt(m.star))
  #X11(); par(mfrow=c(nn,nn))
  # Go each principal component direction searching for local maxima
  for ( i in 1:m.star) {
    proj[,i] = as.vector(dat %*% PC[,i])
    fhat[[i]] = kde(x=proj[,i], h=hpi(proj[,i]))
    # find all local maxima over a neighborhood 
    # of 5 points (5 is deemed sufficient)
    pk = PK.Identify(fhat[[i]]$estimate, 5) 
    ix = which(pk == 1)
    # unlikely that ix is length zero?? but...
    if (length(ix) == 0 ) { stop("Couldn't find any peak in kernel density profile.") }
    pk.levels = fhat[[i]]$estimate[ix]
    pk.highest = max(pk.levels)
    # filter out all the small peaks, which is defined as %cutoff.frac of the
    # highest peak in the kernel density
    ix.sub = which( pk.levels > cutoff.frac * pk.highest)
    # note that ix.sub is at least length 1.
    if( length(ix.sub) <= 5 ) {
      ix = ix[ix.sub]
    } else {
      stop("Too many local maxima.")
      #plot(fhat[[i]]$eval.points, fhat[[i]]$estimate, type="l", main=paste("PC",i))
    }
    # assign a vector to the pk.loc list
    pk.loc[[i]] = fhat[[i]]$eval.points[ix]
  }
  #pt0 = pk.loc[1] * PC[,i] + pk.loc[1] * PC2
  n.vec = sapply(pk.loc, function(ll)length(ll))
  #cat("n.vec:", n.vec, "\n")
  # create the product set :
  n.ll = rep(list(list()), m.star)
  for ( j in 1:m.star ) {
    n.ll[[j]] = 1:n.vec[j]
  }
  product.set = expand.grid(n.ll) # no. rows = # permutation, no. cols = m.star 
  
  N.pts = NROW(product.set)
  pk.mat = matrix(NA, ncol=m.star, nrow=N.pts)
  for ( i in 1:N.pts ) {
    for ( j in 1:m.star ) {
      pk.mat[i,j] = pk.loc[[j]][product.set[i,j]]
    }
  }
  # pk.mat : N.pts x m.star
  # PC : p x m.star
  PP = pk.mat %*% t(PC)
  # result : (N.pts x m.star) * t(p x m.star) = N.pts x p
  #cat("PP:", NROW(PP), "\n")
  
  # Calculate the distance from each point in PP to all the data points in
  # the original coordinate. 
  # CAUTION : This matrix could be huge!!
  D = matrix(NA, ncol=N.pts, nrow=n)
  # D : each column corresponds to a point in PP and each row is the distance of
  # that cluster point to all the actual data points.
  for ( i in 1:N.pts ) {
    D[,i] = dist.one(PP[i,], dat) 
  }
  # For each cluster point, find the 10th percentile, which is the distance
  # needed to go out to encounter 1/10 of the data points. This is a better
  # measure of the proximity to data than the above, which could be distorted
  # by a few outliers.
  qtl.vec = apply(D, 2, quantile, probs=0.1) # pick 10%, with curse of dimensionality consideration  
  if ( N.pts > K ) {
    # pick the cluster points which are closest to some data
    tmp = sort(qtl.vec, decreasing=FALSE)[1:K]
    ix = which(qtl.vec %in% tmp )
    candidates = PP[ix,]
  } else {
    # if there are less than K points, just return all.
    candidates = PP
  }
  if (verbose == TRUE) 
    cat("DONE.\n")    
  return(candidates)
}

#' Assigns a cluster number (index) to each data point
#'
#' @param Cpts n.ctr x p matrix of cluster points returned by Clusters.by.PCA, 
#'        where each row is a point; columns are dimension of the data  
#' @param dat n x p matrix of the data, rows are observation and columns dimension
#' @description Assigns a cluster number(index) to a data point, where the center
#'     of the clusters is given by Cpts. 
#' @return vector of the same length as the number of observations, which tells 
#'          you which cluster center the data point is closest to.
#' @export
data2cluster = function( Cpts, dat ) {
  n = NROW(dat)
  n.ctr = NROW(Cpts)
  Dij = matrix(NA, ncol=n.ctr, nrow=n)
  for ( i in 1:n.ctr ) {
    Dij[,i] = dist.one(as.vector(Cpts[i,]), dat)
  }
  iclass = apply(Dij, 1, which.min)
  return(iclass)
}

#' Vector quantization by k-means
#'
#' @param obs observation matrix
#' @param K number of codewords
#' @description Specifying a K, apply k-means to data. The center of the clusters is mu, while the cluster 
#'        index is returned as the codebook. 
#' @return list of codebk, the codebook, mu.m, the cluster center, and sigma.l, the variance of the clusters.
#' @export
Kmeans.Vector.Quantization = function(obs, K) {
  km = stats::kmeans(obs, K)
  mu.m = km$centers
  codebk = km$cluster
  sigma = rep(list(list()), K)
  for ( i in 1:K ) {
    ix = which(codebk == i)
    sigma[[i]] = var(obs[ix,])
  }
  return(list(codebk=codebk, mu.m=mu.m, sigma.l=sigma))
}


#' Test the PCA algorithm 
#' 
#' @param obs observation matrix
#' @param K K
#' @param M M
#' @description Call Cluster.by.PCA multiple (K) times
#' @return list of Centers, Memb, where Centers is a list of K lists, each being the center returned by Cluster.by.PCA 
#' @export
PCA.algo.test = function(obs, K, M) {
  require(corpcor)
  p = NCOL(obs)
  n = NROW(obs)
  if ( p <= 2 ) {
    m.star = p
  } else {
    m.star = ceiling((p+1)/2) # better estimate ?  
  }
  
  len = ceiling(n / K)
  ix = rep(1:K, each=len) 
  ix = ix[1:n]
  Centers = rep(list(list()), K)
  Memb = rep(list(list()), K)
  for ( i in 1:K ) {
    seg = obs[which(ix==i),]
    seg.var = diag(var(seg))
    Cpts = Cluster.by.PCA(seg, m.star=m.star, K=M, 0.1)
    Centers[[i]] = Cpts
    n.ctrs = NROW(Cpts) # n.ctrs should be exactly M or less.
    icls = data2cluster(Cpts, seg)
    cnts = table(icls)
    Memb[[i]] = cnts
  }
  return(list(Centers=Centers, Memb=Memb))
}

#' Make copies of x by adding noise
#'
#' @param x either a vector of length p or a matrix with p columns
#' @param n number of copies to make
#' @param var.s scalar variance for spherical multivariate noise
#' @description See code
#' @return n x p matrix where each column is a noisy copy of x
#' @export
make.noisy.copies = function(x, n, var.s) {
  require(mvtnorm)
  if (is.matrix(x) == TRUE ) {
    x.v = x[1,]
  } else if ( is.vector(x) == TRUE ) {
    x.v = x
  } else {
    stop("x is neither a vector nor a matrix.")
  }
  p = length(x.v)
  x.copies = matrix(rep(x.v, each=n), ncol=p)
  epsi = rmvnorm(n, mean=rep(0,p), sigma=var.s*diag(p))
  noisy = x.copies + epsi
  return(noisy)
}

#' find peak
#'
#' @param x plain vector of values.
#' @param win.width window size
#' @description identify the location of peaks
#' @return integer vector the length of x, where locations of 1s mark the occurence of a peak.
#' @export
PK.Identify = function(x, win.width=21) {
  require(zoo)
  z = zoo(x)
  pk.mark = pk.scoring(z, win.width)
  y = ifelse( pk.mark == win.width, 1, 0)
  return(y)
}

#' Internal function to show the score from the first pass
#'
#' @param z a zoo object
#' @param win.width window size
#' @description Note requires rollapply from zoo package
#' @return vector of votes, same length as z
#' @export
pk.scoring = function(z, win.width=21) {
  require(zoo)
  inx = rollapply(z, width=win.width, pk.find, by=1, align="left")
  # initialize every time pt starting with zero vote.
  pk.mark = rep(0, length(z))
  # following gathers the votes from every time point so that the highest votes
  # go to the extrema.
  for( i in 1:length(inx) ) {
    pk.mark[inx[i] + (i-1)] =  pk.mark[inx[i] + (i-1)] + 1
  }
  return(pk.mark)
}

#' Locates the maximum value in a vector x and returns the location.
#'
#' @param x numeric vector
#' @description Low level function being called in a rollapply loop by pk.scoring
#' @return vector of peak locations
#' @export
pk.find = function(x) {
  loc=which(x==max(x)) # alternative is to use which.max
  # The following condition check is absolutely needed. It's possible to have
  # **identical** values in x, despite the floating point representation of real numbers. 
  # Just look at the SP500 Close on 1997-1-27 and 1-28!!
  # Pick the first one. Any better choice here ?
  if ( length(loc) > 1 )
    return(loc[1])
  else
    return(loc)
}
