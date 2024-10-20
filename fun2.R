# Description ------------------------------
# Date: 2024-10-01
# Purpose  
# All functions used for full likelihood model
# Finite Mixture project


# Functions ---------------------------------------------------------------

# each row is replicated k times and assigned a group
repl.dat <- function(k, dat) {
  row_seq         <- seq_len(nrow(dat)) 
  rep_row_seq     <- rep(row_seq, each = k)
  t_dat           <- dat[rep_row_seq,]
  t_dat[,"g"]     <- factor(rep(seq_len(k), times = nrow(dat)))
  rownames(t_dat) <- seq_len(nrow(t_dat))
  return(t_dat)
}

# setup the pi_matrix with all possible terms

init.pi.pars <- function(pi_model_matrix, k) {
  n_pars <- ncol(model_matrix_pi)
  pi_par_matrix <- matrix(nrow = n_pars, ncol = k,
                          dimnames = list(colnames(pi_model_matrix),
                                          paste0("g_",seq(1, k))))
  return(pi_par_matrix)
}

init.theta.pars <- function(theta_model_matrix, n_thetas) {
    n_pars <- ncol(theta_model_matrix)
    theta_par_matrix <- matrix(ncol = n_thetas, nrow = n_pars,
                            dimnames = list(colnames(theta_model_matrix),
                                            paste0("theta_",seq(0, n_thetas-1))))
  return(theta_par_matrix)
}

# get.probs()
# Probabilities of group membership is calculated using 
# multinomial logistic function  
# to enable regressing >2 latent groups
get.pi <- function(conc_vars = NULL, dat = NULL, rep_var = "t", pi_pars = NULL) {
  form_pi <- formula(paste(c("~ 1","g",conc_vars), collapse = "+"))
  pi_vars <- c("g","id",conc_vars)
  # model matrix for group probabilities. 
  # collapses nrow by removing repeated measurements
  pi_dat <- dat
  pi_dat$g            <- 1 # the model matrix needs g = 1
  pi_dat              <- unique(pi_dat[, colnames(pi_dat) %in% pi_vars])
  model_matrix_pi    <- model.matrix(form_pi, pi_dat)
  pi_pars            <- cbind(c(0,0),
                              c(1,0.5))
  
  pi_matrix          <- exp(model_matrix_pi %*% pi_pars)
  
  pi_matrix          <- (1 / (1 + rowSums(pi_matrix[,seq(2,k), drop = FALSE]))) * pi_matrix
  colnames(pi_matrix) <- paste0("g_",seq_len(k))
  return(pi_matrix)
}

# This is the "main" function.
# @in theta index and group
# @out theta values (theta_0, theta_1 and so forth)
get.thetas <- function(theta_matrix, theta_model_matrix) {
  
  thetas <- theta_model_matrix %*% theta_pars


  
  
  
  }



# linear transformation of y, t and group
# # h(.) approximates g(.) through Taylor expansion
H.fun <- function(g = NULL, dat = NULL, theta_matrix = NULL,
                  offset = 0) {
  thetas  <- get.thetas(g, dat, theta_matrix)
  loc     <- dat$y - thetas[,1] + offset
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
solve.y <- function(res = NULL, thetas = NULL,
                    dat) {
  
  M     <- cbind(res, thetas)
  
  # each row is passed to polyroot where the real numbered root is returned
  z <- apply(M, 1, \(i) {
      root <- polyroot(c(i[1], 1 + i[3], i[4])) # (1 + theta_1)z, theta_2 z^2
  })
  z <- unlist(z)

  y <- Re(thetas[,1] + z )
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
  # to estimate the density.
  browser()
  plot.h.y(dat, theta_matrix)
  H_matrix_u <<- sapply(seq_len(c_env$k), \(g) H.fun(g = g, dat = dat,
                                                  theta_matrix = theta_matrix,
                                                  offset = 0.01))
  H_matrix_l <<- sapply(seq_len(c_env$k), \(g) H.fun(g = g, dat = dat,
                                                    theta_matrix = theta_matrix,
                                                    offset = -0.01))
  Fij_matrix_h <- pnorm(H_matrix_u)
  Fij_matrix_l <- pnorm(H_matrix_l)
  
  fij_matrix <- Fij_matrix_h - Fij_matrix_l
  
  # there are 0 probabilities, these are set to minimum value
  # There are also iterations where the parameters results in cells 
  # with higher Fij_matrix_l than Fij_matrix_h --> negative values
  # 
  
  # NaNs are generated, this can be avoided by setting all negative values
  # to xmin. In this version, all 0 values are set to minimum value.
  fij_matrix[fij_matrix <= 0] <- .Machine$double.xmin
  
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
  if (!is.finite(l)) {
    l <- Inf }  # Return Inf to prevent invalid results
  print(l)
  return(l)
}


plot.h.y <- function(dat, theta_matrix) {
  
  H_matrix_u <- sapply(seq_len(c_env$k), \(g) H.fun(g = g, dat = dat,
                                                    theta_matrix = theta_matrix,
                                                    offset = 0.01))
  H_matrix_l <- sapply(seq_len(c_env$k), \(g) H.fun(g = g, dat = dat,
                                                    theta_matrix = theta_matrix,
                                                    offset = -0.01))
  H_matrix <- sapply(seq_len(c_env$k), \(g) H.fun(g = g, dat = dat,
                                                  theta_matrix = theta_matrix,
                                                  offset = 0))
  
  # data frame with boolean values if pnorm difference is positive
  p_b_matrix <- pnorm(H_matrix_u) - pnorm(H_matrix_l) >= 0
  colnames(p_b_matrix) <- c("g_1","g_2")
  p_b_df <- data.frame(p_b_matrix) |> 
    pivot_longer(everything(), names_to = "g",
                 values_to = "pos_p")
  
  
  t_df <- data.frame(H_matrix) |> 
    rename(any_of(c("g_1" = "X1", "g_2" = "X2"))) |> 
    bind_cols(true_g = dat$true_g,
              y = dat$y,
              t = dat$t) |> 
    pivot_longer(c("g_1","g_2"),
                 names_to = "g",
                 values_to = "h.y") |> 
    mutate(g = str_extract(g, "\\d")) |> 
    mutate(true_g_b = if_else(g == true_g, "true","false")) |> 
    bind_cols(pos_p = p_b_df$pos_p)
  
  plot.obj <- ggplot(t_df, aes(x = y, y = h.y, color = true_g_b)) +
    geom_point() +
    facet_wrap(~ pos_p)
  
  print(plot.obj)
}

