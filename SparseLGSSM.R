#########################################################################################
######### Writing class SparseLGSSM (Sparse Linear Gaussian State Space Model) ##########
#########################################################################################

# Dependencies --------------------------------------------------------------------------
library(mvtnorm)
library(matrixcalc)
library(glmnet)
library(glasso)
library(R.utils)


# Generic methods -----------------------------------------------------------------------
sim_par <- function(x, ...) {
  UseMethod('sim_par')
}



sim_data <- function(x, ...) {
  UseMethod('sim_data')
}



kalman_filtering <- function(x, ...) {
  UseMethod('kalman_filtering')
}



kalman_smoothing <- function(x, ...) {
  UseMethod('kalman_smoothing')
}



E_step <- function(x, ...) {
  UseMethod('E_step')
}



M_step <- function(x, ...) {
  UseMethod('M_step')
}



EM <- function(x, ...) {
  UseMethod('EM')
}



tune_lambda <- function(x, ...) {
  UseMethod('tune_lambda')
}



# Methods for the SparseLGSSM class -----------------------------------------------------
create_LGSSM <- function(data, p, N, type='LWF') {
  # If data is missing, i.e. we run simulation
  if(missing(data)) {
    model <- structure(
      list(
        p = p,
        N = N,
        type = type
      ),
      class='SparseLGSSM'
    )
  }
  else { # If data is provided by the user
    X0 <- data$X0
    Y <- scale(data$Y)
    N <- nrow(Y)
    p <- ncol(Y)
    
    model <- structure(
      list(
        data = list(X=X0, Y=Y),
        N = N,
        p = p,
        type = type
      ),
      class='SparseLGSSM'
    )
  }
  
  return(model)
}



alt_search_B <- function(theta, M1, M2, M3, lambda_1, N) {
  tmp1 <- kronecker(M3, theta)
  tmp1 <- chol(tmp1)
  A <- sqrt(2 * N) * tmp1
  y <- sqrt(N / 2) * tmp1 %*% as.vector(t(M2) %*% solve(M3))
  beta <- glmnet(x=A, y=y, lambda=lambda_1, family='gaussian',
                 alpha=1, standardize=F, intercept = F,
                 standardize.response = F)$beta
  
  matrix(as.numeric(beta), ncol=ncol(M1))
}



alt_search_theta <- function(B, theta_prev, M1, M2, M3, lambda_2) {
  S <- M1 - B %*% M2 + B %*% M3 %*% t(B)
  S <- (S + t(S)) / 2
  glasso(s=S, rho=lambda_2)$wi
} 



make_W <- function(p, lambda_1, lambda_2) {
  W <- matrix(0, ncol = 2*p, nrow = 2*p)
  W[1:p, 1:p] <- lambda_2
  #diag(W[1:p, 1:p]) <- 0
  W[(p+1):(2*p), 1:p] <- W[1:p, (p+1):(2*p)] <- lambda_1/2
  W
}



glasso_suf_stat <- function(M1, M2, M3) {
  tmp1 <- cbind(M1, t(M2)/2)
  tmp2 <- cbind(M2/2, M3)
  M <- rbind(tmp1, tmp2)
  
  if (!is.positive.semi.definite(M))
    stop('M is not positive semi-definite')
  
  M
}



compute_improvement <- function(x, model, B, Lambda, Theta, Sigma) {
    x <- rbind(x, c(
      norm(B - model$par_est$B, type='F') / norm(model$par_est$B, type='F'),
      norm(Lambda - model$par_est$Lambda, type='F') / norm(model$par_est$Lambda, type='F'),
      norm(Theta - model$par_est$Theta, type='F') / norm(model$par_est$Theta, type='F'),
      norm(Sigma - model$par_est$Sigma, type='F') / norm(model$par_est$Sigma, type='F')))
    
    x[is.nan(x)] <- 0
    x
}



log_lik_sample <- function(X, B, Theta) {
  N <- nrow(X)
  tmp <- 0
  for (t in 2:nrow(X)) {
    tmp <- tmp + (- 0.5 * as.numeric(determinant(Theta, logarithm=T)$modulus) +
      0.5 * t(X[t, ] - B %*% X[t-1, ]) %*% Theta %*% (X[t, ] - B %*% X[t-1, ])) / N
  }
  as.numeric(tmp)
}



sim_par.SparseLGSSM <- function(model, rho2=0.1, s1=0.75, s2=0.75) {
  type <- model$type
  p <- model$p
  if (!(type=='LWF' | type=='AMP'))
    stop('Model penalization scheme must be LWF or AMP')
  
  # Simulating concentration matrix with sparsity s1 until we get a valid one
  tmp <- T
  while (tmp) {
    A <- matrix(rnorm(p^2, mean=0, sd=0.2), ncol=p, nrow=p)
    Theta <- t(A) %*% A
    diag(Theta) <- diag(Theta) + 1
    
    # Simulate entries to set to zero
    ind <- expand.grid(1:p, 1:p)
    ind <- ind[sample(1:p^2, size = floor(s2*p^2/1.3)), ]
    exc <- (ind[, 1] == ind[, 2])
    ind <- as.matrix(ind[!exc, ])
    
    for (i in 1:nrow(ind)) {
      foo <- ind[i, ]
      Theta[foo[1], foo[2]] <- Theta[foo[2], foo[1]] <- 0
    }
    
    Theta <- (Theta + t(Theta)) / 2
    tmp <- !(is.positive.definite(Theta))
  }
  
  Sigma <- solve(Theta)
  Sigma <- (Sigma + t(Sigma)) / 2
  
  model$par_true$Theta <- Theta
  model$par_true$Sigma <- Sigma
  
  if (type=='LWF') {
    # Simulating sparse Lambda matrix ---------------------
    Lambda <- matrix(rnorm(p^2, mean=0, sd=2), ncol=p, nrow=p)
    ind <- expand.grid(1:p, 1:p)
    ind <- ind[sample(1:p^2, size = floor(s1*p^2)), ]
    
    for (i in 1:nrow(ind)) {
      foo <- as.numeric(ind[i, ])
      Lambda[foo[1], foo[2]] <- 0
    }
    
    # Making sure that B has spectral radius less than 1
    B <- Sigma %*% Lambda
    tmp <- eigen(B)$values
    max(sqrt(Re(tmp)^2 + Im(tmp)^2))
    B <- B / (1.001 * max(sqrt(Re(tmp)^2 + Im(tmp)^2)))
    Lambda <- Theta %*% B
    
    Lambda[Lambda < 1e-11] <- 0
    
    model$par_true$Lambda <- Lambda
    model$par_true$B <- B
  }
  else {
    # Simulating sparse B matrix
    B <- matrix(rnorm(p^2, mean=0, sd=2), ncol=p, nrow=p)
    ind <- expand.grid(1:p, 1:p)
    ind <- ind[sample(1:p^2, size = floor(s1*p^2)), ]

    for (i in 1:nrow(ind)) {
      foo <- as.numeric(ind[i, ])
      B[foo[1], foo[2]] <- 0
    }
    
    tmp <- eigen(B)$values
    max(sqrt(Re(tmp)^2 + Im(tmp)^2))
    B <- B / (1.001 * max(sqrt(Re(tmp)^2 + Im(tmp)^2)))
    
    model$par_true$B <- B
    model$par_true$Lambda <- model$par_true$Theta %*% B
  }
  
  model$par_true$Omega <- rho2 * mean(Sigma) * diag(p)
  model$par_true$C <- diag(p)
  
  return(model)
}



sim_data.SparseLGSSM <- function(model) {
  if(is.null(model$par_true))
    stop('True parameters must be provided when simulating data')
  
  # Extracting model parameters
  p <- model$p
  N <- model$N
  
  Omega <- model$par_true$Omega
  Sigma <- model$par_true$Sigma
  C <- model$par_true$C
  B <- model$par_true$B
  
  # Simulating data
  Y <- X <- matrix(NA, nrow=N, ncol=p)
  
  X[1, ] <- rmvnorm(n=1, mean=rep(0, p), sigma=Sigma)
  Y[1, ] <- X[1, ] + rmvnorm(n=1, mean= rep(0, p), sigma=Omega)
  
  eps <- rmvnorm(n=N-1, mean=rep(0, p), sigma=Sigma)
  nu <- rmvnorm(n=N-1, mean=rep(0, p), sigma=Omega)
  
  for (t in 2:N) {
    X[t, ] <- t(B %*% X[t-1, ]) + eps[t-1, ]
    Y[t, ] <- t(C %*% X[t, ]) + nu[t-1, ]
  }
  
  X <- scale(X, center = F, scale = F)
  Y <- scale(Y, center = F, scale = F)
  
  model$data <- list(X=X, Y=Y)
  return(model)
}



kalman_filtering.SparseLGSSM <- function(model) {
  # Extracting working parameter estimates and data
  p <- model$p
  N <- model$N
  
  Omega <- model$par_est$Omega
  Sigma <- model$par_est$Sigma
  B <- model$par_est$B
  C <- model$par_est$C
  Y <- model$data$Y
  
  # Variables for storing prediction and filtering estimates
  X_pred <- matrix(NA, ncol=p, nrow=N) # Note: first row always NAs
  X_filt <- matrix(NA, ncol=p, nrow=N)
  P_pred <- array(NA, dim = c(p, p, N)) # Note: first matrix always NAs
  P_filt <- array(NA, dim=c(p, p, N))

  # Computing initial filtering estimates
  X_filt[1, ] <- model$data$X[1, ]
  P_filt[, , 1] <- 0
  
  # Alternating between prediction and filtering
  for (t in 2:N) {
    # Prediction
    X_pred[t, ] <- B %*% X_filt[t-1, ]
    P_pred[, , t] <- B %*% tcrossprod(P_filt[, , t-1], B) + Sigma
    P_pred[, , t] <- (P_pred[, , t] + t(P_pred[, , t])) / 2
    # Filtering
    tmp <- solve(C %*% tcrossprod(P_pred[, , t], C) + Sigma)
    X_filt[t, ] <- X_pred[t, ] + P_pred[, , t] %*% tmp %*% (Y[t, ] - C %*% X_pred[t, ])
    P_filt[, , t] <- P_pred[, , t] - tcrossprod(P_pred[, , t], C) %*% tmp %*% C %*% P_pred[, , t]
    P_filt[, , t] <- (P_filt[, , t] + t(P_filt[, , t])) / 2
  }
  
  # Saving results to model object and returning
  model$X_pred <- X_pred
  model$X_filt <- X_filt
  model$P_pred <- P_pred
  model$P_filt <- P_filt
  
  return(model)
}



kalman_smoothing.SparseLGSSM <- function(model) {
  # Extracting model parameters and data
  p <- model$p
  N <- model$N
  
  Omega <- model$par_est$Omega
  Sigma <- model$par_est$Sigma
  B <- model$par_est$B
  C <- model$par_est$C
  
  # Doing Kalman filtering if not already done
  if (is.null(model$X_pred)) {
    model <- kalman_filtering(model)
  }
  
  X_pred <- model$X_pred
  X_filt <- model$X_filt
  P_pred <- model$P_pred
  P_filt <- model$P_filt
  
  # Variables for storing smoothed estimates
  X_smooth <- matrix(NA, ncol=p, nrow=N)
  P_smooth <- array(NA, dim=c(p, p, N))
  # Cov(X_t-1, X_t | Y_0:N) makes sense (here) for t = 2,...,N so first slash in array empty
  Cov_cross_smooth <- array(NA, dim=c(p, p, N))
  
  # First smoothed variable is simply the filtered variable
  X_smooth[N, ] <- X_filt[N, ]
  P_smooth[, , N] <- P_filt[, , N]
  
  # Running the smoother using predicted and filtered estimates
  for (t in (N-1):1) {
    G <- P_filt[, , t] %*% B %*% solve(P_pred[, , t+1])
    X_smooth[t, ] <- X_filt[t, ] + G %*% (X_smooth[t+1, ] - X_pred[t+1, ])
    P_smooth[, , t] <- P_filt[, , t] + G %*% tcrossprod(P_smooth[, , t+1] - P_pred[, , t+1], G)
    P_smooth[, , t] <- (P_smooth[, , t] + t(P_smooth[, , t])) / 2
    Cov_cross_smooth[, , t] <- G %*% P_smooth[, , t+1]
  }
  
  # Saving results to model object and returning
  model$X_smooth <- X_smooth
  model$P_smooth <- P_smooth
  model$Cov_cross_smooth <- Cov_cross_smooth
  
  return(model)
}



E_step.SparseLGSSM <- function(model) {
  # Perform forwards Kalman prediction and filtering pass
  model <- kalman_filtering(model)
  
  # Perform backwards Kalman smoothing pass
  model <- kalman_smoothing(model)
  
  # Extract the information that we need
  N <- model$N
  p <- model$p
  X_smooth <- model$X_smooth
  P_smooth <- model$P_smooth
  Cov_cross_smooth <- model$Cov_cross_smooth

  #Compute the expected sufficient statistics M1, M2 and M3
  M1 <- M2 <- matrix(0, ncol=p, nrow=p)

  for (t in 2:N) {
    tmp <- outer(X_smooth[t, ], X_smooth[t, ])
    M1 <- M1 + (tmp + t(tmp)) / 2 + P_smooth[, , t]
  }

  for (t in 1:(N-1)) {
    M2 <- 2 * outer(X_smooth[t, ], X_smooth[t+1, ]) + 2 * Cov_cross_smooth[, , t]
  }
  
  M1 <- M1 / N
  M2 <- M2 / N
  
  M3 <- (t(X_smooth[1, , drop=F]) %*% X_smooth[1, , drop=F] - 
           t(X_smooth[N, , drop=F]) %*% X_smooth[N, , drop=F] -
           P_smooth[, , N]) / N + M1
  
  # Save the resulting matrices
  model$M1 <- M1
  model$M2 <- M2
  model$M3 <- M3
  
  return(model)
}



M_step.SparseLGSSM <- function(model) {
  # Checking model penalization scheme
  type <- model$type
  
  if(!(type=='LWF' | type=='AMP'))
    stop('Must choose either LWF or AMP')
  
  # Extracting sufficient statistics from E-step
  M1 <- model$M1
  M2 <- model$M2
  M3 <- model$M3
  lambda_1 <- model$lambda_1
  lambda_2 <- model$lambda_2
  p <- model$p
  N <- model$N
  
  # Doing M-step with glasso if LWF is chosen
  if(type=='LWF') {
    # Prepare glasso
    M <- glasso_suf_stat(M1, M2, M3)
    W <- make_W(p=p, lambda_1=lambda_1, lambda_2=lambda_2)
    
    # Running glasso
    Th <- glasso(s=M, rho=W)$wi

    # Extracting parameter estimates
    Theta <- Th[1:p, 1:p]
    Lambda <- -Th[1:p, (p+1):(2*p)]
    Sigma <- solve(Theta)
    Sigma <- (Sigma + t(Sigma)) / 2
    B <- Sigma %*% Lambda
    
    # Saving relative improvement from last estimate
    model$convergence$impr_rel <- compute_improvement(model=model, x=model$convergence$impr_rel,
                                                      B=B, Lambda=Lambda, Theta=Theta, Sigma=Sigma)

    # Saving current parameter updates
    model$par_est$Sigma <- Sigma
    model$par_est$B <- B
    model$par_est$Theta <- Theta
    model$par_est$Lambda <- Lambda
  }
  
  if (type=='AMP') {
    B <- model$par_est$B
    theta <- model$par_est$Theta
    tmp <- Inf
    
    while (tmp > 1e-6) {
      theta_new <- alt_search_theta(B=B, theta_prev=theta, M1=M1, M2=M2, M3=M3, lambda_2=lambda_2)
      B_new <- alt_search_B(theta=theta, M1=M1, M2=M2, M3=M3, lambda_1=lambda_1, N=N)
      
      if (norm(B, type='F')==0) {
        tmp <- norm(theta - theta_new, type='F') / norm(theta, type='F')
      } else {
        tmp <- max(norm(B - B_new, type='F') / norm(B, type='F'),
                   norm(theta - theta_new, type='F') / norm(theta, type='F'))
      }
      
      theta <- theta_new
      B <- B_new
    }
    
    Sigma <- solve(theta)
    Lambda <- theta %*% B
    
    # Saving relative improvement from last estimate
    model$convergence$impr_rel <- compute_improvement(model=model, x=model$convergence$impr_rel,
                                                      B=B, Lambda=Lambda, Theta=theta, Sigma=Sigma)
    
    model$par_est$B <- B
    model$par_est$Theta <- theta
    model$par_est$Sigma <- Sigma
    model$par_est$Lambda <- Lambda
  }
  
  return(model)
}



EM.SparseLGSSM <- function(model, lambda_1, lambda_2, max_iter=50,
                           eps=1e-6, type, par_init, print_progress=F) {
  if (missing(type))
    type <- model$type
  p <- model$p
  model$lambda_1 <- lambda_1
  model$lambda_2 <- lambda_2
  
  if (!(type %in% c('LWF', 'AMP')))
    stop('Either LWF or AMP parametrization must be chosen')
  
  # Creating initial values for the optimization
  if(missing(par_init)) {
    # We simply use the true Omega and C
    model$par_init$Omega <- model$par_true$Omega
    model$par_init$C <- model$par_true$C
    
    # We use a simulated B matrix with small entries
    model$par_init$B <- matrix(rnorm(p^2, sd=0.1), ncol=p)
    
    # We use some diagonal for Theta 
    model$par_init$Theta <- diag(rnorm(p, mean=10, sd=1))
    model$par_init$Sigma <- solve(model$par_init$Theta)
    
    model$par_init$Lambda <- model$par_init$Theta %*% model$par_init$B
  }
  
  # Using initial parameters as starting values
  model$par_est <- model$par_init
  
  # Saving relative distance to last estimate
  model$convergence$impr_rel <- matrix(Inf, ncol=4, dimnames = list(NULL, c('B', 'Lambda', 'Theta', 'Sigma')))
  model$convergence$it <- 1
  
  while (model$convergence$it < max_iter &
         max(model$convergence$impr_rel[model$convergence$it, ], na.rm = T) > eps) {
    model <- E_step(model)
    model <- M_step(model)
    
    if (print_progress==T) {
      message(cat(paste('EM iteration', model$convergence$it, 'done')))
    }
    
    model$convergence$it <- model$convergence$it + 1
    if (model$convergence$it > max_iter) {
      stop('Maximum number of iterations reached and EM did not converge')
    }
  }
  
  model$convergence$it <- model$convergence$it - 1
  model$convergence$impr_rel <- model$convergence$impr_rel[-1, ]
  
  # Deleting irrelevant quantities
  model$M1 <- model$M2 <- model$M3 <- model$X_pred <-
    model$X_filt <- model$X_smooth <- model$P_filt <-
    model$P_pred <- model$P_smooth <- model$Cov_cross_smooth <- NULL
  
  return(model)
}


tune_lambda.SparseLGSSM <- function(model, fit_type, B1=100, B2=100,
                                    lambda_1_min, lambda_1_max,
                                    lambda_2_min, lambda_2_max,
                                    grid_size) {
  # Grid of lambdas
  lambda <- as.matrix(expand.grid(exp(seq(log(lambda_1_min), log(lambda_1_max), length=grid_size)),
                                  exp(seq(log(lambda_2_min), log(lambda_2_max), length=grid_size))))
  
  # Computing risk for each lambda
  risk <- rep(NA, nrow(lambda))
  for (i in 1:length(risk)) {
    message(cat(paste('Running EM for lambda number', i)))
    lambda_1 <- lambda[i, ][1]
    lambda_2 <- lambda[i, ][2]
    
    set.seed(123)
    
    # Averaging over EM estimates
    tmp <- numeric(B1)
    for (j in 1:B1) {
      model <- sim_data(model)
      
      par_est <- NULL
      par_est$B <- NA
      par_est$Theta <- NA
      
      try(expr = {
        withTimeout(
          expr = {
            par_est <- EM(model, lambda_1=lambda_1, lambda_2=lambda_2, type=fit_type)$par_est
            message(cat(paste('EM fitting nr', j, 'done')))
          },
          onTimeout = 'error', elapsed = 3, timeout = 3
        )},
        silent = T)

      if (!is.na(par_est$B)) {
        foo <- numeric(B2)
        for (k in 1:B2) {
          X <- sim_data(model)$data$X
          foo[k] <- log_lik_sample(X=X, B=par_est$B, Theta=par_est$Theta)
        }
        tmp[j] <- mean(foo, na.rm=T)
      } else {
        tmp[j] <- NA
      }
    }
    risk[i] <- mean(tmp, na.rm=T)
  }
  
  return(lambda[which.min(risk), ])
}