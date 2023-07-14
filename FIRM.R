# production function and marginal products
fY     <- function(K,L,TFP) TFP*K^alpha*L^(1-alpha);
MPK    <- function(K,L,TFP) TFP*alpha*(K/L)^(alpha-1);
MPKinv <- function(MPKin,L,TFP) L*(MPKin/(TFP*alpha))^(1/(alpha-1));
MPL    <- function(K,L,TFP) TFP*(1-alpha)*(K/L)^alpha;

# invert the firm problem: finds r for given V (and LS, TFP and taxes)
rdemand <- function(assetsupply, maxiter = 20, tol = 1e-6, verbose = F) {

  K2      = K;
  uck2    = uck;
  r2      = r;
  qTob2   = qTob;
  
  error   = Inf;
  iter    = 0;
  
  repeat {
    
    iter           = iter + 1;
    error_old      = error;
    qTob_old       = qTob2;
    
    K2             = assetsupply/qTob2;
    #K2[1]          = K[1]; #predetermined
    uck2           = MPK(K2,LD,TFP);
    qTob2          = (1-tauprof)*uck2 + tauprof*delta + (1-delta);
    
    error = sum(abs(qTob2-qTob_old));
    
    if (verbose == T) cat(paste("Iteration:\t",iter,"\t\tError:\t",error,"\n"));
    
    if (iter > maxiter)    {if (verbose) cat("No convergence!!\n"); break;}
    if (error < tol)       {if (verbose) cat("Convergence!!\n"); break;}
    if (error > error_old) {if (verbose) cat("Increasing error: stop at previous step!\n"); qTob2 = qTob_old; break;}
    
  }
  
  r2[1:(tend-1)] = qTob2[2:tend]-1;
  r2[tend] = r2[tend-1];
  
  return(r2)
}

