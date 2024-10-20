# Description ------------------------------
# Date: 2024-10-01
# Purpose  
# All functions used for full likelihood model
# Finite Mixture project


# Functions ---------------------------------------------------------------

# setup the pi_matrix with all possible terms

init.pi.matrix <- function(conc_vars = NULL, n_grps = NULL) {
  main_terms <- c("(Intercept)",conc_vars)
  pi_matrix <- matrix(NA, nrow = n_grps - 1, ncol = length(main_terms),
                      dimnames = list(paste0("g_",seq(2,n_grps)),
                                      main_terms))
  return(pi_matrix)
}

init.theta.matrix <- function(n_grps = NULL, n_thetas = 3, 
                              n_times) {
  grp_vector        <- paste0("g_",2:n_grps) # group terms (g_2, g_3, ...)
  rep_vector        <- paste0("t_", seq(2, n_times)) # time terms (t_2, t_3, ...)
  interaction_terms <- c(outer(rep_vector, grp_vector, paste, sep = ":"))
  main_terms        <- c(rep_vector, grp_vector)
  all_terms         <- c("(Intercept)",main_terms, interaction_terms)
  
  all_thetas        <- paste0("theta_",seq(0, n_thetas - 1))
  
  theta_matrix      <- matrix(NA, nrow = length(all_thetas), ncol = length(all_terms),
                              dimnames = list(all_thetas, all_terms))
  return(theta_matrix)
}

# get.probs()
# Probabilities of group membership is calculated using 
# multinomial logistic function  
# to enable regressing >2 latent groups
get.pi <- function(dat = NULL, pi_matrix = NULL) {
  # retrieve common vars
  conc_vars           <- c_env$conc_vars
  id_col              <- c_env$id_col
  
  # get the data frame with only concomitant variables
  # then setup the model matrix and calculate the odds 
  df_only_conc_vars   <- unique(dat[,c(id_col,conc_vars), drop = FALSE]) # if repeated measurements, this gets rid of duplicates
  model_formula       <- as.formula(paste("~ 1", paste(conc_vars, collapse = " + "))) 
  M                   <- model.matrix(model_formula, data = df_only_conc_vars) 
  logits              <- exp(M %*% t(pi_matrix))          # relative probability P(g = k) / P(g = K) 
  sum_logits          <- rowSums(logits)                  # sum of exp(x * beta_j)
  p_base              <- 1 / (1 + sum_logits)             # P(g = K)
  p_matrix            <- cbind(p_base, p_base * logits)   # P(g = k)  
  c_names             <- paste0("g_",seq(1,nrow(pi_matrix)+1))
  colnames(p_matrix)  <- c_names
  return(p_matrix)                                        # matrix of probabilites per cluster and group
}

# This is the "main" function.
# @in theta index and group
# @out theta values (theta_0, theta_1 and so forth)
get.thetas <- function(g_vector = NULL, dat = NULL, 
                       theta_matrix = NULL) {
  # access common variables
  k                 <- get("k", c_env)
  time_col          <- get("time_col", c_env)
  
  theta_str         <- paste0("theta_",seq_len(nrow(theta_matrix))) # theta string to index theta_matrix
  
  # if only group number is supplied, replicate to nrow data frame
  if (length(g_vector) == 1) g_vector <- rep(g_vector, nrow(dat)) else {
    if (length(g_vector) != nrow(dat)) stop("g_vector needs to be same length as data frame")
  }
  
  # data frame to be used in model matrix
  # adds a factor variable with k levels that is expanded to dummy variables
  t_df <- cbind(dat,
                g = factor(g_vector, levels = seq_len(k)))
  t_df[,time_col] <- as.factor(t_df[,time_col])
  form <- as.formula(paste("~ ", time_col, "* g"))
  M <- model.matrix(form, data = t_df) # model matrix with dummies
  if(any(colnames(M) != gsub("_","",colnames(theta_matrix)))) stop("Colnames of design matrix and theta matrix do not match")
  
  # Each row of design matrix is multiplied by each theta vector
  theta <- M %*% t(theta_matrix)
  return(theta)
}



# linear transformation of y, t and group
# # h(.) approximates g(.) through Taylor expansion
h.fun <- function(g = NULL, dat = NULL, theta_matrix = NULL,
                  epsilon = 0) {
  thetas  <- get.thetas(g, dat, theta_matrix)
  loc     <- dat$y - thetas[,1] + epsilon
  h       <- loc + thetas[,2] * loc + thetas[,3] * loc ^ 2
  return(h)
}


# given thetas, group and h(y), solves for y
# where h(y) = residual = y - theta_0 + theta_1 (y - theta_0) + theta_2 (y - theta_0) ^ 2
# let z = y - theta_0
# then h(z) = z + theta_1 * z + theta_2 * z^2 + , ... , + theta_m * z^m
#           = (1 + theta_1) * z + theta_2 * z^2 ... 
# increasing order of terms:
# 0 = -e_term , (1 + theta_0) , theta_1 , theta_2, ... , theta_j
# which is passed to polyroot()
solve.y <- function(res = NULL, g = NULL, dat = NULL, 
                    theta_matrix = NULL) {
  
  thetas <- get.thetas(g, dat, theta_matrix) 
  
  M <- cbind(res, thetas)
  
  # each row is passed to polyroot where the real numbered root is returned
  z <- apply(M, 1, \(i) {
    
    # if residual is 0 return 0 instead of zero length vector
    if (i[1] == 0) {
      root <- polyroot(c(1 + i[3], i[4])) 
      if (length(root) == 0) root <- 0 # if zero length, set to 0
      return(Re(root[1]))
    } else {
      root <- polyroot(c(-i[1], 1 + i[3], i[4])) 
      return(Re(root[1]))
    }})
  
  # z = y - theta_0
  y <- thetas[,1] + z 
  return(y)
}

# Likelihood function

l.fun <- function(pars = NULL, pars_list = NULL, theta_matrix = NULL,
                  pi_matrix = NULL, dat = NULL) {
  
  # get common vars
  k       <- get("k", c_env)
  id_col  <- get("id_col", c_env)
  
  # assign theta_matrix and pi_matrix
  # (Intercept)"   t_2"         "t_3"         "g_2"       "t_2:g_2"     "t_3:g_2" 
  new_pars_list   <- relist(pars, pars_list)
  theta_matrix    <- new_pars_list$thetas
  pi_matrix       <- new_pars_list$pis
  
  # H.fun approximates g(y).
  # An offset/epsilon of 0.01 is passed to the H.fun to enable use of the CDF
  # to estimate the probability.
  h_matrix <- sapply(seq_len(c_env$k), \(g) h.fun(g = g, dat = dat,
                                                  theta_matrix = theta_matrix))
  
  fij_matrix <- dnorm(h_matrix)
  
  # derivative of h.fun  
  d_matrix_u <- sapply(seq_len(c_env$k), \(g) h.fun(g = g, dat = dat,
                                                  theta_matrix = theta_matrix,
                                                  epsilon = 0.001))
  d_matrix_l <- sapply(seq_len(c_env$k), \(g) h.fun(g = g, dat = dat,
                                                  theta_matrix = theta_matrix,
                                                  epsilon = -0.001))
                    
                   
  pen_matrix <- (dnorm(d_matrix_u) - dnorm(d_matrix_l)) / (0.001 * 2)
  # if positive, no penalty
  pen_matrix[pen_matrix > 0] <- 0
  # if negative incur penalty
  pen_matrix[pen_matrix < 0] <- 1e-6 * -pen_matrix[pen_matrix < 0]
  
  # add penalty to density matrix
  fij_matrix <- fij_matrix - pen_matrix
  
  # repeated measurements are multiplied per individual
  fi_matrix <- aggregate(fij_matrix,
                         by = list(id = dat[,c_env$id_col]),
                         prod)
  fi_matrix <- fi_matrix[, colnames(fi_matrix) != c_env$id_col] # remove id column
  
  # matrix of group probabilites
  pi <- get.pi(dat = dat, pi_matrix = pi_matrix)
  
  # element multiplication of pi and fi_matrix
  pi_fi <- pi * fi_matrix
  
  # each row is summed and logged
  li <- log(rowSums(pi_fi))
  l <- -sum(li)
  
  print(l)
  return(l)
}




