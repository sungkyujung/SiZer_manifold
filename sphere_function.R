library(scales)
library(expm)
library(pracma)
library(vMF)
library(plotrix)
library(reticulate)

np <- import("numpy")
sci <- import("scipy")

generate_grid_sphere = function(grid_num){
  lon = pi*(-1+sqrt(5))*(-grid_num:grid_num)
  lat = 2*(-grid_num:grid_num)/(2*grid_num+1)
  
  x = cos(lon)*sqrt(1-lat^2)
  y = sin(lon)*sqrt(1-lat^2)
  z = lat
  
  grid_sph <- cbind(x, y, z)
  return(grid_sph)
}

t_stat_cal = function(x){
  wing <- x[1:2]
  est_var_inv <- matrix(c(x[4], -x[5], -x[5], x[3]), nrow=2)/(x[3]*x[4]-x[5]^2)
  t_stat <- t(wing)%*%est_var_inv%*%wing
  return(t_stat)
}

get_quantile_for_hess = function(sigma, n, alpha, m){
  n <- as.integer(n)
  dim <- ncol(sigma)[1]
  grid <- dim/3
  alpha_star <- 1-(1-alpha)^(1/m)

  data_c <- np$random$multivariate_normal(mean=rep(0, dim), cov=sigma, size=n, tol=1e-20)
  eigen_bdd <- array(NA, dim = c(grid, 4))
  
  idx <- 1
  for(i in 1:grid){
    data <- t(data_c[, idx:(idx+2)])
    com <- data[1, ] + data[2, ]
    pm <- sqrt((data[1, ] - data[2, ])^2 + 4*data[3, ]^2)
    
    eigen_max_boot <- Re((com + pm)/2)
    eigen_min_boot <- Re((com - pm)/2)
    
    q_max <- quantile(eigen_max_boot, probs=c(alpha_star/4, 1-alpha_star/4))
    q_min <- quantile(eigen_min_boot, probs=c(alpha_star/4, 1-alpha_star/4))
    
    eigen_bdd[i, ] <- as.vector(c(q_max, q_min))
    idx <- idx+3
  }
  return(eigen_bdd)
}

get_boots_for_grad = function(rot_data, k, idx, grid, n, grad_1, grad_2){
  rot_data_b <- rot_data[, idx, ]
  k_rot_data_b <- rot_data_b*k
  c_exp_ky1_b <- CpvMF(2, k)*exp(k_rot_data_b[1, , ])
  
  grad_mat_1_b <- -k_rot_data_b[3, , ]*c_exp_ky1_b
  grad_mat_2_b <- k_rot_data_b[2, , ]*c_exp_ky1_b
  grad_1_b <- colMeans(grad_mat_1_b)
  grad_2_b <- colMeans(grad_mat_2_b)
  
  t_stat_b <- NULL
  
  for (j in 1:grid){
    est_var_b <- cov(cbind(grad_mat_1_b[, j], grad_mat_2_b[, j]))/n
    t_stat_b <- c(t_stat_b, t_stat_cal(c(grad_1_b[j]-grad_1[j], grad_2_b[j]-grad_2[j], est_var_b[1, 1], est_var_b[2, 2], est_var_b[1, 2])))
  }
  
  return(max(abs(t_stat_b)))
}

get_arrow = function(grid_pt, grad, scope){
  pt <- grid_pt
  
  endpt <- c(1, grad[2]/scope, -grad[1]/scope)
  endpt <- endpt/norm(endpt, "2")
  
  if(pt[1]!=-1){
    endpt <- endpt + endpt[1]*c(pt[1]-1, pt[2], pt[3]) - endpt[2]*pt[2]/(1+pt[1])*c(1+pt[1], pt[2], pt[3]) - endpt[3]*pt[3]/(1+pt[1])*c(1+pt[1], pt[2], pt[3])
  }else{
    endpt <- c(-endpt[1:2], endpt[3])
  }
  
  start <- pt[1:2]/sqrt(1+abs(pt[3]))
  end <- endpt[1:2]/sqrt(1+abs(endpt[3]))
  
  return(end-start)
}

get_grad = function(grid_pt, grad){
  pt <- grid_pt
  grad_res <- c(0, -grad[2], grad[1])
  
  if(pt[1]!=-1){
    grad_res <- grad_res - grad_res[2]*pt[2]/(1+pt[1])*c(1+pt[1], pt[2], pt[3]) - grad_res[3]*pt[3]/(1+pt[1])*c(1+pt[1], pt[2], pt[3])
  }else{
    grad_res <- c(-grad_res[1:2], grad_res[3])
  }
  return(grad_res)
}

sphere_sizer = function(data, grid_num, seq_kappa, B, alpha, method){
  n <- size(data)[2]
  grid_pt <- generate_grid_sphere(grid_num)
  grid <- size(grid_pt)[1]
  
  rot_data <- array(NA, dim = c(3, n, grid))
  
  for(idx in 1:grid){
    pt <- grid_pt[idx, ]
    if(pt[1]!=-1){
      rot_data[, , idx] <- data + matrix(c(pt[1]-1, -pt[2], -pt[3]))%*%data[1, ] + pt[2]/(1+pt[1])*matrix(c(1+pt[1], -pt[2], -pt[3]))%*%data[2, ] + pt[3]/(1+pt[1])*matrix(c(1+pt[1], -pt[2], -pt[3]))%*%data[3, ]
    }else{
      rot_data[, , idx] <- rbind(-data[1:2, ], data[3, ])
    }
  }
  
  par_length <- length(seq_kappa)
  t_stat <- array(NA, dim = c(grid, par_length))
  arrow_1 <- array(NA, dim = c(grid,  par_length))
  arrow_2 <- array(NA, dim = c(grid,  par_length))
  ESS <- array(NA, dim = c(grid,  par_length))
  m <- rep(NA, par_length)
  
  sample_idx <- matrix(sample(1:n, n*B, replace=T), nrow=B)
  quant_boot <- matrix(nrow=par_length, ncol=B)
  grad_1 <- array(0, dim = c(grid, par_length))
  grad_2 <- array(0, dim = c(grid, par_length))
  
  result <- as.list(NA)
  
  for (i in 1:par_length){
    kappa <- seq_kappa[i]
    kappa_rot_data <- rot_data*kappa
    c_exp_ky1 <- CpvMF(2, kappa)*exp(kappa_rot_data[1, , ])

    grad_mat_1 <- -kappa_rot_data[3, , ]*c_exp_ky1
    grad_mat_2 <- kappa_rot_data[2, , ]*c_exp_ky1
    grad_1[, i] <- colMeans(grad_mat_1)
    grad_2[, i] <- colMeans(grad_mat_2)
    
    hess_mat_11 <- (-kappa_rot_data[1, , ]+kappa_rot_data[3, , ]^2)*c_exp_ky1
    hess_mat_22 <- (-kappa_rot_data[1, , ]+kappa_rot_data[2, , ]^2)*c_exp_ky1
    hess_mat_12 <- -kappa_rot_data[2, , ]*kappa_rot_data[3, , ]*c_exp_ky1
    hess_mat_cat <- array(NA, dim = c(3*grid, n))
    hess_11 <- colMeans(hess_mat_11)
    hess_22 <- colMeans(hess_mat_22)
    hess_12 <- colMeans(hess_mat_12)
    
    color <- rep('green', grid)
    
    ESS[, i] <- colSums(exp(kappa_rot_data[1, , ]))/exp(kappa)
    ESS_pos <- ESS[, i]
    m[i] <- n/mean(ESS_pos[ESS_pos>=5])
    
    idx <- 1
    for(j in 1:grid){
      hess_mat_cat[idx, ] <- hess_mat_11[, j]
      idx <- idx+1
      hess_mat_cat[idx, ] <- hess_mat_22[, j]
      idx <- idx+1
      hess_mat_cat[idx, ] <- hess_mat_12[, j]
      idx <- idx+1
    }
    
    hess_var_est <- np$cov(hess_mat_cat)
    eigen_bdd <- get_quantile_for_hess(hess_var_est/n, B, alpha, m[i])

    for (j in 1:grid){
      est_var <- cov(cbind(grad_mat_1[, j], grad_mat_2[, j]))/n
      t_stat[j, i] <- t_stat_cal(c(grad_1[j, i], grad_2[j, i], est_var[1, 1], est_var[2, 2], est_var[1, 2]))
    }
    
    com <- hess_11 + hess_22
    pm <- sqrt((hess_11 - hess_22)^2 + 4*hess_12^2)
    
    eigen_max <- (com + pm)/2
    eigen_min <- (com - pm)/2
    
    color[eigen_max <= eigen_bdd[, 1]] <- 'blue'
    color[eigen_min >= eigen_bdd[, 4]] <- 'yellow'
    color[color == 'green' & eigen_max >= eigen_bdd[, 2] & eigen_min <= eigen_bdd[, 3]] <- 'red'
    color[color == 'green' & eigen_max >= eigen_bdd[, 2]] <- 'orange'
    color[color == 'green' & eigen_min <= eigen_bdd[, 3]] <- 'purple'
    color[ESS_pos<5] <- 'gray'
    
    if(method==3 || method==4){
      quant_boot[i, ] <- apply(sample_idx, 1, function(idx){get_boots_for_grad(rot_data, kappa, idx, grid, n, grad_1[, i], grad_2[, i])})
    }
    
    result[[i]] <- list(kappa = kappa,
                        grid_pt = t(grid_pt),
                        gradient = rbind(grad_1[, i], grad_2[, i]),
                        hessian = rbind(hess_11, hess_22, hess_12),
                        eigen_max = eigen_max,
                        eigen_min = eigen_min,
                        eigen_max_bdd = rbind(eigen_bdd[, 1], eigen_bdd[, 2]),
                        eigen_min_bdd = rbind(eigen_bdd[, 3], eigen_bdd[, 4]),
                        color = color,
                        kde = c_exp_ky1,
                        ESS = ESS_pos,
                        m = m[i])
  }
  
  if(method==1){
    quant <- qchisq(1-alpha, 2)
  }else if(method==2){
    quant <- qchisq((1-alpha)^(1/m), 2)
  }else if(method==3){
    quant <- apply(quant_boot, 1, function(x){quantile(x, probs=1-alpha, type=1)})
  }else if(method==4){
    quant <- apply(quant_boot, 2, function(x){max(x)})
    quant <- quantile(quant, probs=c(1-alpha), type=1)
  }
  
  for (i in 1:par_length){
    if(method==2 || method==3){
      quant_i <- quant[i]
    }else{
      quant_i <- quant
    }
    
    gridsize <- 2/sqrt(grid_num)
    maxgrad <- max(sqrt(grad_1[, i]^2+grad_2[, i]^2))
    scope <- 1*max(maxgrad/gridsize)
    
    result[[i]]$quantile <- quant_i
    result[[i]]$t_statistics <- as.vector(t_stat[, i])
    result[[i]]$signif <- (t_stat[, i] > quant_i)
    result[[i]]$arrow <- apply(rbind(t(grid_pt), grad_1[, i], grad_2[, i]), 2, function(x){get_arrow(x[1:3], x[4:5], scope)})
    result[[i]]$gradient_res <- apply(rbind(t(grid_pt), grad_1[, i], grad_2[, i]), 2, function(x){get_grad(x[1:3], x[4:5])})
  }
  return(result)
}