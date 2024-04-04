library(scales)
library(expm)
library(plotly)
library(reticulate)

np <- import("numpy")
sci <- import("scipy")

t_stat_cal = function(x){
  wing <- x[1:2]
  est_var_inv <- matrix(c(x[4], -x[5], -x[5], x[3]), nrow=2)/(x[3]*x[4]-x[5]^2)
  t_stat <- t(wing)%*%est_var_inv%*%wing
  return(t_stat)
}

get_quantile_for_hess = function(sigma, n, alpha, m){
  n <- as.integer(n)
  dim <- ncol(sigma)[1]
  grid <- sqrt(dim/3)
  alpha_star <- 1-(1-alpha)^(1/m)
  
  data_c <- np$random$multivariate_normal(mean=rep(0, dim), cov=sigma, size=n, tol=1e-20)
  eigen_bdd <- array(NA, dim = c(grid, grid, 4))
  
  idx <- 1
  for(i in 1:grid){
    for(j in 1:grid){
      data <- t(data_c[, idx:(idx+2)])
      com <- data[1, ] + data[2, ]
      pm <- sqrt((data[1, ] - data[2, ])^2 + 4*data[3, ]^2)
      
      eigen_max_boot <- Re((com + pm)/2)
      eigen_min_boot <- Re((com - pm)/2)
      
      q_max <- quantile(eigen_max_boot, probs=c(alpha_star/4, 1-alpha_star/4))
      q_min <- quantile(eigen_min_boot, probs=c(alpha_star/4, 1-alpha_star/4))
      
      eigen_bdd[i, j, ] <- as.vector(c(q_max, q_min))
      idx <- idx+3
    }
  }
  return(eigen_bdd)
}

get_boots_for_grad = function(x1, x2, t1, t2, kappa, idx, grid, n, grad_1, grad_2){
  x1_b <- x1[idx]
  x2_b <- x2[idx]
  
  grid_x1_b <- outer(t1, x1_b, "-")
  grid_x2_b <- outer(t2, x2_b, "-")
  
  bess_const <- (2*pi*besselI(kappa, 0))^2
  
  exp_x1_b <- exp(kappa*cos(grid_x1_b))
  exp_x2_b <- exp(kappa*cos(grid_x2_b))
  pder_1_b <- -kappa*sin(grid_x1_b)*exp_x1_b
  pder_2_b <- -kappa*sin(grid_x2_b)*exp_x2_b
  
  grad_1_b <- 0
  grad_mat_1_b <- array(NA, dim = c(grid, grid, n))
  grad_2_b <- 0
  grad_mat_2_b <- array(NA, dim = c(grid, grid, n))
  for (j in 1:n){
    grad_mat_1_b[, , j] <- pder_1_b[, j]%*%t(exp_x2_b[, j])
    grad_1_b <- grad_1_b+grad_mat_1_b[, , j]
    grad_mat_2_b[, , j] <- exp_x1_b[, j]%*%t(pder_2_b[, j])
    grad_2_b <- grad_2_b+grad_mat_2_b[, , j]
  }
  grad_mat_1_b <- grad_mat_1_b/bess_const
  grad_1_b <- grad_1_b/(n*bess_const)
  grad_mat_2_b <- grad_mat_2_b/bess_const
  grad_2_b <- grad_2_b/(n*bess_const)
  
  t_stat_b <- matrix(nrow=grid, ncol=grid)
  for (j in 1:grid){
    est_var_11_b <- (rowMeans(grad_mat_1_b[j, , ]^2)-rowMeans(grad_mat_1_b[j, , ])^2)/(n-1)
    est_var_22_b <- (rowMeans(grad_mat_2_b[j, , ]^2)-rowMeans(grad_mat_2_b[j, , ])^2)/(n-1)
    est_var_12_b <- (rowMeans(grad_mat_1_b[j, , ]*grad_mat_2_b[j, , ])-rowMeans(grad_mat_1_b[j, , ])*rowMeans(grad_mat_2_b[j, , ]))/(n-1)
    t_stat_b[j, ] <- apply(rbind(grad_1_b[j, ]-grad_1, grad_2_b[j, ]-grad_2, est_var_11_b, est_var_22_b, est_var_12_b), 2, t_stat_cal)
  }
  
  return(max(abs(t_stat_b)))
}

torus_sizer = function(x1, x2, grid, seq_kappa, B, alpha, method){
  n <- length(x1)
  t1 <- (seq(0, 2*pi, length=grid+1))[1:grid]
  grid_x1 <- outer(t1, x1, "-")
  t2 <- (seq(0, 2*pi, length=grid+1))[1:grid]
  grid_x2 <- outer(t2, x2, "-")
  par_length <- length(seq_kappa)
  t_stat <- array(NA, dim = c(grid, grid, par_length))
  signif_save <- array(NA, dim = c(grid, grid, par_length))
  
  sample_idx <- matrix(sample(1:n, n*B, replace=T), nrow=B)
  quant_boot <- matrix(nrow=par_length, ncol=B)
  grad_1 <- array(0, dim = c(grid, grid, par_length))
  grad_2 <- array(0, dim = c(grid, grid, par_length))
  ESS <- array(0, dim = c(grid, grid, par_length))
  m <- rep(0, length=par_length)
  
  result <- as.list(NA)
  
  for (i in 1:par_length){
    kappa <- seq_kappa[i]
    bess_const <- (2*pi*besselI(kappa, 0))^2
    
    exp_x1 <- exp(kappa*cos(grid_x1))
    exp_x2 <- exp(kappa*cos(grid_x2))
    pder_1 <- -kappa*sin(grid_x1)*exp_x1
    pder_2 <- -kappa*sin(grid_x2)*exp_x2
    phes_1 <- (kappa^2*(sin(grid_x1))^2-kappa*cos(grid_x1))*exp_x1
    phes_2 <- (kappa^2*(sin(grid_x2))^2-kappa*cos(grid_x2))*exp_x2
    
    kde <- array(NA, dim = c(grid, grid, n))
    grad_mat_1 <- array(NA, dim = c(grid, grid, n))
    grad_mat_2 <- array(NA, dim = c(grid, grid, n))
    hess_mat_11 <- array(NA, dim = c(grid, grid, n))
    hess_mat_22 <- array(NA, dim = c(grid, grid, n))
    hess_mat_12 <- array(NA, dim = c(grid, grid, n))
    hess_mat_cat <- array(NA, dim = c(3*grid^2, n))
    for (j in 1:n){
      kde[, , j] <- exp_x1[, j]%*%t(exp_x2[, j])
      grad_mat_1[, , j] <- pder_1[, j]%*%t(exp_x2[, j])
      grad_mat_2[, , j] <- exp_x1[, j]%*%t(pder_2[, j])
      hess_mat_11[, , j] <- phes_1[, j]%*%t(exp_x2[, j])
      hess_mat_22[, , j] <- exp_x1[, j]%*%t(phes_2[, j])
      hess_mat_12[, , j] <- pder_1[, j]%*%t(pder_2[, j])
    }
    kde <- kde/bess_const
    grad_mat_1 <- grad_mat_1/bess_const
    grad_mat_2 <- grad_mat_2/bess_const
    hess_mat_11 <- hess_mat_11/bess_const
    hess_mat_22 <- hess_mat_22/bess_const
    hess_mat_12 <- hess_mat_12/bess_const
    kde <- apply(kde, c(1, 2), mean)
    grad_1[, , i] <- apply(grad_mat_1, c(1, 2), mean)
    grad_2[, , i] <- apply(grad_mat_2, c(1, 2), mean)
    hess_11 <- apply(hess_mat_11, c(1, 2), mean)
    hess_22 <- apply(hess_mat_22, c(1, 2), mean)
    hess_12 <- apply(hess_mat_12, c(1, 2), mean)
    
    eigen_max <- matrix(nrow=grid, ncol=grid)
    eigen_min <- matrix(nrow=grid, ncol=grid)
    color <- matrix('green', nrow=grid, ncol=grid)
    
    ESS[, , i] <- exp_x1%*%t(exp_x2)/exp(2*kappa)
    ESS_pos <- ESS[, , i]
    m[i] <- n/mean(ESS_pos[ESS_pos>=5])
    
    idx <- 1
    for(j in 1:grid){
      for(k in 1:grid){
        hess_mat_cat[idx, ] <- hess_mat_11[j, k, ]
        idx <- idx+1
        hess_mat_cat[idx, ] <- hess_mat_22[j, k, ]
        idx <- idx+1
        hess_mat_cat[idx, ] <- hess_mat_12[j, k, ]
        idx <- idx+1
      }
    }
    
    hess_var_est <- np$cov(hess_mat_cat)
    eigen_bdd <- get_quantile_for_hess(hess_var_est/n, B, alpha, m[i])
    
    for (j in 1:grid){
      est_var_11 <- (rowMeans(grad_mat_1[j, , ]^2)-rowMeans(grad_mat_1[j, , ])^2)/(n-1)
      est_var_22 <- (rowMeans(grad_mat_2[j, , ]^2)-rowMeans(grad_mat_2[j, , ])^2)/(n-1)
      est_var_12 <- (rowMeans(grad_mat_1[j, , ]*grad_mat_2[j, , ])-rowMeans(grad_mat_1[j, , ])*rowMeans(grad_mat_2[j, , ]))/(n-1)
      t_stat[j, , i] <- apply(rbind(grad_1[j, , i], grad_2[j, , i], est_var_11, est_var_22, est_var_12), 2, t_stat_cal)
      
      com <- hess_11[j, ] + hess_22[j, ]
      pm <- sqrt((hess_11[j, ] - hess_22[j, ])^2 + 4*hess_12[j, ]^2)
      if(!all(Im(pm)==0)){
        if(sum(abs(Im(pm)))<1e-10){
          pm <- Re(pm)
        }else{
          stop("Error : the complexity of eigenvalues")
        }
      }
      eigen_max[j, ] <- (com + pm)/2
      eigen_min[j, ] <- (com - pm)/2
    }
    
    color[eigen_max <= eigen_bdd[, , 1]] <- 'blue'
    color[eigen_min >= eigen_bdd[, , 4]] <- 'yellow'
    color[color == 'green' & eigen_max >= eigen_bdd[, , 2] & eigen_min <= eigen_bdd[, , 3]] <- 'red'
    color[color == 'green' & eigen_max >= eigen_bdd[, , 2]] <- 'orange'
    color[color == 'green' & eigen_min <= eigen_bdd[, , 3]] <- 'purple'
    color[ESS_pos<5] <- 'gray'
    
    if(method==3 || method==4){
      quant_boot[i, ] <- apply(sample_idx, 1, function(idx){get_boots_for_grad(x1, x2, t1, t2, h, idx, grid, n, grad_1[, , i], grad_2[, , i])})
    }
    
    result[[i]] <- list(kappa = kappa,
                        grid_pt = rbind(rep(t1, times=grid), rep(t2, each=grid)),
                        gradient = rbind(as.vector(grad_1), as.vector(grad_2)),
                        hessian = rbind(as.vector(hess_11), as.vector(hess_22), as.vector(hess_12)),
                        eigen_max = as.vector(eigen_max),
                        eigen_min = as.vector(eigen_min),
                        eigen_max_bdd = rbind(as.vector(eigen_bdd[, , 1]), as.vector(eigen_bdd[, , 2])),
                        eigen_min_bdd = rbind(as.vector(eigen_bdd[, , 3]), as.vector(eigen_bdd[, , 4])),
                        color = as.vector(color),
                        kde = as.vector(kde),
                        ESS = as.vector(ESS_pos),
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
    
    signif <- (t_stat[, , i] > quant_i)
    
    x_gridsize <- t1[2]-t1[1]
    y_gridsize <- t2[2]-t2[1]
    x_maxgrad <- max(grad_1[, , i])
    y_maxgrad <- max(grad_2[, , i])
    scope <- 1.2*max(x_maxgrad/x_gridsize, y_maxgrad/y_gridsize)
    
    arrow_1 <- grad_1[, , i]/scope
    arrow_2 <- grad_2[, , i]/scope
    
    result[[i]]$quantile <- quant_i
    result[[i]]$t_statistics <- as.vector(t_stat[, , i])
    result[[i]]$signif <- as.vector(signif)
    result[[i]]$arrow <- rbind(as.vector(arrow_1), as.vector(arrow_2))
  }
  return(result)
}