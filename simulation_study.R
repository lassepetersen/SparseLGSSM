# Loading the SparseLGSSM class
source('SparseLGSSM.R')
library(reshape2)
library(tidyverse)
library(xtable)
theme_set(theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=20)))



# Convergence properties --------------------------------------------------------------------------
if (F) {
  set.seed(1)
  model_conv1 <- create_LGSSM(p=40, N=200, type='LWF')
  model_conv1 <- sim_par(model_conv1, rho2 = 0.1, s1=0.75, s2=0.75)
  model_conv2 <- create_LGSSM(p=40, N=200, type='AMP')
  model_conv2 <- sim_par(model_conv2, rho2 = 0.1, s1=0.75, s2=0.75)
  
  model_conv1 <- sim_data(model_conv1)
  model_conv1 <- EM(model_conv1, lambda_1=0.001, lambda_2=0.1,
                    print_progress = T, eps = 1e-9)
  
  model_conv2 <- sim_data(model_conv2)
  model_conv2 <- EM(model_conv2, lambda_1=0.001, lambda_2=0.1,
                    print_progress = T, eps = 1e-9)
  
  tmp1 <- as.data.frame(model_conv1$convergence$impr_rel)
  tmp2 <- as.data.frame(model_conv2$convergence$impr_rel)
  
  tmp1 <- cbind(tmp1, Iteration=1:model_conv1$convergence$it, type='LWF')
  tmp2 <- cbind(tmp2, Iteration=1:model_conv2$convergence$it, type='AMP')
  
  tmp <- rbind(tmp1, tmp2)
  tmp <- melt(tmp, measure.vars = 1:4)
  colnames(tmp)[3] <- 'Variable'
  
  ggplot(tmp, aes(x=Iteration, y=log(value), color=Variable)) +
    geom_point(size=1) +
    facet_wrap(~ type) +
    ylab(expression(log(RD))) + 
    ggtitle('EM Convergence') +
    scale_color_manual(labels = c("B", expression(Lambda), expression(Theta), expression(Sigma)),
                       values=c('forestgreen', 'firebrick1', 'blue2', 'darksalmon'))
  
    rm(model_conv1, model_conv2, tmp1, tmp2, tmp)
}



# Choosing optimal lambdas for the case p = 3 and N = 100 -----------------------------------------
set.seed(12345670)
p <- 3
N <- 100

model_LWF <- create_LGSSM(p = p, N = N, type = 'LWF')
model_LWF <- sim_par(model_LWF, rho2 = 0.1, s1 = 0.6, s2 = 0.2)

set.seed(1234)
model_AMP <- create_LGSSM(p = p, N = N, type = 'AMP')
model_AMP <- sim_par(model_AMP, rho2 = 0.1, s1 = 0.70, s2 = 0.3)

#### WARNING - The following grid search can take several hours - WARNING ###
if (F) {
  lambda1 <- tune_lambda(model_LWF,
                         fit_type = 'LWF',
                         B1 = 20,
                         B2 = 100,
                         lambda_1_min = 0.00001,
                         lambda_1_max = 1,
                         lambda_2_min = 0.001,
                         lambda_2_max = 5,
                         grid_size = 20)
  save(lambda1, file='lambda1')
}

#### WARNING - The following grid search can take several hours - WARNING ###
if (F) {
  lambda2 <- tune_lambda(model_AMP,
                         fit_type = 'AMP',
                         B1 = 20,
                         B2 = 100,
                         lambda_1_min = 0.00001,
                         lambda_1_max = 1,
                         lambda_2_min = 0.001,
                         lambda_2_max = 5,
                         grid_size = 20)
  save(lambda2, file='lambda2')
}

#### WARNING - The following grid search can take several hours - WARNING ###
if (F) {
  lambda3 <- tune_lambda(model_AMP,
                         fit_type = 'LWF',
                         B1 = 20,
                         B2 = 100,
                         lambda_1_min = 0.00001,
                         lambda_1_max = 1,
                         lambda_2_min = 0.001,
                         lambda_2_max = 5,
                         grid_size = 20)
  save(lambda3, file='lambda3')
}

#### WARNING - The following grid search can take several hours - WARNING ###
if (F) {
  lambda4 <- tune_lambda(model_LWF,
                         fit_type = 'AMP',
                         B1 = 20,
                         B2 = 100,
                         lambda_1_min = 0.00001,
                         lambda_1_max = 1,
                         lambda_2_min = 0.001,
                         lambda_2_max = 5,
                         grid_size = 20)
  save(lambda4, file='lambda4')
}



# Comparing the estimators ------------------------------------------------------------------------
if (F) {
  load('lambda1')
  load('lambda2')
  load('lambda3')
  load('lambda4')
}

# Bootstrapping
if (F) {
  n_sim <- 1000
  Lambda1 <- Lambda4 <- matrix(NA, ncol = p*p, nrow = n_sim)
  B2 <- B3 <- matrix(NA, ncol = p*p, nrow = n_sim)
  Theta1 <- Theta2 <- Theta3 <- Theta4 <- matrix(NA, ncol = p*p, nrow = n_sim)
  
  
  for (i in 1:n_sim) {
    print(i)
    
    Lambda <- rep(NA, p^2)
    Theta <- rep(NA, p^2)
    
    try(expr = {
    model_LWF <- sim_data(model_LWF)
    
    withCallingHandlers(expr = {
      setTimeLimit(elapsed = 3)
      par_est <- EM(model_LWF, type='LWF', lambda_1 = lambda1[1], lambda_2 = lambda1[2])$par_est
    }, error = function(e) stop('EM did not converge'))
    
    Lambda <- as.numeric(par_est$Lambda);
    Theta <- as.numeric(par_est$Theta);
    }, 
    silent = T)
    
    Lambda1[i, ] <- Lambda
    Theta1[i, ] <- Theta
  }
  
  for (i in 1:n_sim) {
    print(i)
    
    B <- rep(NA, p^2)
    Theta <- rep(NA, p^2)
    
    try(expr = {
      model_AMP <- sim_data(model_AMP)
      
      withCallingHandlers(expr = {
        setTimeLimit(elapsed = 3)
        par_est <- EM(model_AMP, type='AMP', lambda_1 = lambda2[1], lambda_2 = lambda2[2])$par_est
      }, error = function(e) stop('EM did not converge'))
      
      B <- as.numeric(par_est$B);
      Theta <- as.numeric(par_est$Theta);
    }, 
    silent = T)
    
    B2[i, ] <- B
    Theta2[i, ] <- Theta
}
  
  for (i in 1:n_sim) {
    print(i)
    
    B <- rep(NA, p^2)
    Theta <- rep(NA, p^2)
    
    try(expr = {
      model_AMP <- sim_data(model_AMP)
      
      withCallingHandlers(expr = {
        setTimeLimit(elapsed = 3)
        par_est <- EM(model_AMP, type='LWF', lambda_1 = lambda3[1], lambda_2 = lambda3[2])$par_est
      }, error = function(e) stop('EM did not converge'))
      
      B <- as.numeric(par_est$B);
      Theta <- as.numeric(par_est$Theta);
    }, 
    silent = T)
    
    B3[i, ] <- B
    Theta3[i, ] <- Theta
  }
  
  for (i in 1:n_sim) {
    print(i)
    
    Lambda <- rep(NA, p^2)
    Theta <- rep(NA, p^2)
    
    try(expr = {
      model_LWF <- sim_data(model_LWF)
      
      withCallingHandlers(expr = {
        setTimeLimit(elapsed = 3)
        par_est <- EM(model_LWF, type='AMP', lambda_1 = lambda4[1], lambda_2 = lambda4[2])$par_est
      }, error = function(e) stop('EM did not converge'))
      
      Lambda <- as.numeric(par_est$Lambda);
      Theta <- as.numeric(par_est$Theta);
    }, 
    silent = T)
    
    Lambda4[i, ] <- Lambda
    Theta4[i, ] <- Theta
  }
  
  bootstraps <- list(Theta1=Theta1, Theta2=Theta2, Theta3=Theta3, Theta4=Theta4,
                     B2=B2, B3 = B3, Lambda1=Lambda1, Lambda4=Lambda4)
  save(bootstraps, file='bootstraps')
}

# Creating summaries for the thesis

create_latex_matrix <- function(x) {
  x=xtable(x,align=rep("",ncol(x)+1))
  print(x, floating=FALSE, tabular.environment="pmatrix", 
        hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
}

# 1
create_latex_matrix(matrix(colMeans(bootstraps$Theta1 == 0, na.rm = T), ncol = 3))
create_latex_matrix(matrix(colMeans(bootstraps$Lambda1 == 0, na.rm = T), ncol = 3))

# 2
create_latex_matrix(matrix(colMeans(bootstraps$Theta2 == 0, na.rm = T), ncol = 3))
create_latex_matrix(matrix(colMeans(bootstraps$B2 == 0, na.rm = T), ncol = 3))

# 3
create_latex_matrix(matrix(colMeans(bootstraps$Theta3 == 0, na.rm = T), ncol = 3))
create_latex_matrix(matrix(colMeans(bootstraps$B3 == 0, na.rm = T), ncol = 3))

# 4
create_latex_matrix(matrix(colMeans(bootstraps$Theta4 == 0, na.rm = T), ncol = 3))
create_latex_matrix(matrix(colMeans(bootstraps$Lambda4 == 0, na.rm = T), ncol = 3))


# The true graph structure
find_common_structure <- function(theta, x) {
  mat <- cbind(theta, x)
  mat <- na.omit(mat) != 0
  ind <- table(apply(mat, 1, paste, collapse = "/"))
  print(max(ind))
  ind <- which.max(ind)
  tmp <- as.logical(strsplit(names(ind), "/")[[1]])
  
  tmp1 <- matrix(tmp[1:9], ncol = 3)
  tmp2 <- matrix(tmp[10:18], ncol = 3)
  
  list(tmp1, tmp2)
}

# 1
find_common_structure(bootstraps$Theta1, bootstraps$Lambda1)

# 2
find_common_structure(bootstraps$Theta2, bootstraps$B2)

# 3
find_common_structure(bootstraps$Theta3, bootstraps$B3)

# 4
find_common_structure(bootstraps$Theta4, bootstraps$Lambda4)

# Comparing noise sizes rho2=0.1, rho2=0.05, rho2=0.01 --------------------------------------------
set.seed(12345670)
p <- 3
N <- 100
model_LWF_1 <- create_LGSSM(p = p, N = N, type = 'LWF')
model_LWF_1 <- sim_par(model_LWF_1, rho2 = 0.05, s1 = 0.6, s2 = 0.2)

set.seed(12345670)
model_LWF_2 <- create_LGSSM(p = p, N = N, type = 'LWF')
model_LWF_2 <- sim_par(model_LWF_2, rho2 = 0.01, s1 = 0.6, s2 = 0.2)

if (F) {
  lambda5 <- tune_lambda(model_LWF_1,
                         fit_type = 'LWF',
                         B1 = 20,
                         B2 = 100,
                         lambda_1_min = 0.00001,
                         lambda_1_max = 1,
                         lambda_2_min = 0.001,
                         lambda_2_max = 5,
                         grid_size = 20)
  save(lambda5, file='lambda5')
}

if (F) {
  lambda6 <- tune_lambda(model_LWF_2,
                         fit_type = 'LWF',
                         B1 = 20,
                         B2 = 100,
                         lambda_1_min = 0.00001,
                         lambda_1_max = 1,
                         lambda_2_min = 0.001,
                         lambda_2_max = 5,
                         grid_size = 20)
  save(lambda6, file='lambda6')
}

# Boostrapping
Lambda5 <- Lambda6 <- matrix(NA, ncol = p*p, nrow = n_sim)
Theta5 <- Theta6 <- matrix(NA, ncol = p*p, nrow = n_sim)

set.seed(123)
for (i in 1:n_sim) {
  print(i)
  
  Lambda <- rep(NA, p^2)
  Theta <- rep(NA, p^2)
  
  try(expr = {
    model_LWF_1 <- sim_data(model_LWF_1)
    
    withCallingHandlers(expr = {
      setTimeLimit(elapsed = 3)
      par_est <- EM(model_LWF_1, type='LWF', lambda_1 = lambda5[1], lambda_2 = lambda5[2])$par_est
    }, error = function(e) stop('EM did not converge'))
    
    Lambda <- as.numeric(par_est$Lambda);
    Theta <- as.numeric(par_est$Theta);
  }, 
  silent = T)
  
  Lambda5[i, ] <- Lambda
  Theta5[i, ] <- Theta
}

set.seed(12345)
for (i in 1:n_sim) {
  print(i)
  
  Lambda <- rep(NA, p^2)
  Theta <- rep(NA, p^2)
  
  try(expr = {
    model_LWF_2 <- sim_data(model_LWF_2)
    
    withCallingHandlers(expr = {
      setTimeLimit(elapsed = 3)
      par_est <- EM(model_LWF_2, type='LWF', lambda_1 = lambda6[1], lambda_2 = lambda6[2])$par_est
    }, error = function(e) stop('EM did not converge'))
    
    Lambda <- as.numeric(par_est$Lambda);
    Theta <- as.numeric(par_est$Theta);
  }, 
  silent = T)
  
  Lambda6[i, ] <- Lambda
  Theta6[i, ] <- Theta
}

# 5
create_latex_matrix(matrix(colMeans(bootstraps$Theta5 == 0, na.rm = T), ncol = 3))
create_latex_matrix(matrix(colMeans(bootstraps$Lambda5 == 0, na.rm = T), ncol = 3))

find_common_structure(bootstraps$Theta5, bootstraps$Lambda5)

# 6
create_latex_matrix(matrix(colMeans(bootstraps$Theta6 == 0, na.rm = T), ncol = 3))
create_latex_matrix(matrix(colMeans(bootstraps$Lambda6 == 0, na.rm = T), ncol = 3))

find_common_structure(bootstraps$Theta6, bootstraps$Lambda6)

