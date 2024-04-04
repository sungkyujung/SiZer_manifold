##############################################
#run file 'sphere_function.R' first###########
##############################################
#data: spherical data (matrix, size of 3 times N)
#grid_num: the value of m (scalar) / generate Fibonacci grid of totally 2m+1 points
#seq_kappa: the sequence of the value of kappa (vector)
#B: bootstrap repeatation
#alpha: significance level
#method: 1(Gaussian), 2(simul. Gaussian), 3(simul. bootstrap for grid), 4(simul. bootstrap for grid, kappa)
##############################################
#Line 12-22 is an example: put your data here#
##############################################

n <- 1000
theta <- rnorm(10*n, mean = pi/3, sd = 0.1)
theta <- theta[theta>=0]
theta_len <- length(theta)
rv_unif <- runif(theta_len, min = 0, max = 2*pi)
theta <- theta[sin(theta)>=rv_unif]
theta <- theta[1:n]
phi <- runif(n, min = 0, max = 2*pi)
data <- rbind(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))

seq_kappa <- 10^seq(0, 2, length=200)
res <- sphere_sizer(data, grid_num=800, seq_kappa, B=10000, alpha=0.05, method=2)

##############################################
#Do not edit below this#######################
##############################################
#figure of Hessian############################
##############################################

n <- 100
col_pos <- NULL

for(k in 1:n){
  col_pos <- rbind(col_pos, cbind(rep(k, 8*k-4), 1:(8*k-4)))
}

k <- col_pos[, 1]
i <- col_pos[, 2]
rad_pt <- (2*k-1)/(2*n)
phi_pt <- 2*pi*(2*i-1)/(16*k-8)

rad_lb <- (k-1)/n
rad_ub <- k/n
phi_lb <- 2*pi*(i-1)/(8*k-4)
phi_ub <- 2*pi*i/(8*k-4)

z <- 1-rad_pt^2
x <- sqrt(1-z^2)*cos(phi_pt)
y <- sqrt(1-z^2)*sin(phi_pt)

grid_pt <- res[[1]]$grid_pt
grid_pos <- grid_pt[, grid_pt[3, ]>=0]
col_idx_pos <- NULL

for(i in 1:(4*n^2)){
  pt <- c(x[i], y[i], z[i])
  col_idx_pos <- c(col_idx_pos, which.max(pt%*%grid_pos))
}

z <- -z
grid_neg <- grid_pt[, grid_pt[3, ]<=0]
col_idx_neg <- NULL

for(i in 1:(4*n^2)){
  pt <- c(x[i], y[i], z[i])
  col_idx_neg <- c(col_idx_neg, which.max(pt%*%grid_neg))
}

for(i in 1:length(res)){
  kappa <- res[[i]]$kappa
  grid_pt <- res[[i]]$grid_pt
  arrow <- res[[i]]$arrow
  color <- res[[i]]$color
  ESS <- res[[i]]$ESS
  signif <- res[[i]]$signif
  
#  png(paste("Rplot_", i, ".png"), width=1200, height=600) #If you want to save as files, run this line along with other line
  par(mfrow=c(1, 2))
  
  label.pos <- seq(0, 2*pi, length=13)
  radial.plot(0, label.pos=label.pos, start=0, clockwise=F, show.grid.labels=F, show.grid=F,
              rp.type="s", radial.lim=c(0, 1), point.symbols=16,
              point.col="white", show.centroid=F, grid.col=NULL,
              labels=c("0", "30", "60", "90", "120", "150", "180", "210", "240", "270", "300", "330"))
  
  for(j in 1:(4*n^2)){
    polygon(c(rad_lb[j]*cos(phi_ub[j]), rad_lb[j]*cos(phi_lb[j]), rad_ub[j]*cos(phi_lb[j]), rad_ub[j]*cos(phi_ub[j])), 
            c(rad_lb[j]*sin(phi_ub[j]), rad_lb[j]*sin(phi_lb[j]), rad_ub[j]*sin(phi_lb[j]), rad_ub[j]*sin(phi_ub[j])), 
            col=(color[grid_pt[3, ]>=0])[col_idx_pos[j]], density=NA)
  }
  
  segments(rep(0, 9), rep(0, 9), cos(label.pos), sin(label.pos))
  
  xcir <- sqrt(1-abs(cos(seq(0, pi/2, length=7))))
  ncir <- length(xcir)
  
  for(idx in 1:ncir){
    radial.plot(rep(xcir[idx], 250), seq(0,2*pi,length=250), rp.type="p", radial.lim=c(0, 1), add=T)
  }
  
  xpos <- xcir
  ypos <- rep(0, ncir)
  boxed.labels(xpos, ypos, c("0", "15", "30", "45", "60", "75", "90"), border = FALSE, cex=1)
  
  for(j in 1:length(ESS)){
    if(ESS[j]>=5 & grid_pt[3, j]>=0){
      pt <- grid_pt[, j]
      start <- pt[1:2]/sqrt(1+abs(pt[3]))
      
      if(signif[j]==TRUE){
        end <- start + arrow[, j]
        arrows(start[1], start[2], end[1], end[2], length = 0.05)
      }else{
        rad <- as.numeric(sqrt(start[1]^2+start[2]^2))
        angle <- as.numeric(sign(start[2])*acos(start[1]/rad))
        radial.plot(rad, angle, rp.type="s", point.symbols=16, add=T)
      }
    }
  }
  
  label.pos <- seq(0, 2*pi, length=13)
  radial.plot(0, label.pos=label.pos, start=0, clockwise=F, show.grid.labels=F, show.grid=F,
              rp.type="s", radial.lim=c(0, 1), point.symbols=16,
              point.col="white", show.centroid=F, grid.col=NULL,
              labels=c("0", "30", "60", "90", "120", "150", "180", "210", "240", "270", "300", "330"))
  
  for(j in 1:(4*n^2)){
    polygon(c(rad_lb[j]*cos(phi_ub[j]), rad_lb[j]*cos(phi_lb[j]), rad_ub[j]*cos(phi_lb[j]), rad_ub[j]*cos(phi_ub[j])), 
            c(rad_lb[j]*sin(phi_ub[j]), rad_lb[j]*sin(phi_lb[j]), rad_ub[j]*sin(phi_lb[j]), rad_ub[j]*sin(phi_ub[j])), 
            col=(color[grid_pt[3, ]<=0])[col_idx_neg[j]], density=NA)
  }
  
  segments(rep(0, 9), rep(0, 9), cos(label.pos), sin(label.pos))
  
  xcir <- sqrt(1-abs(cos(seq(0, pi/2, length=7))))
  ncir <- length(xcir)
  
  for(idx in 1:ncir){
    radial.plot(rep(xcir[idx], 250), seq(0,2*pi,length=250), rp.type="p", radial.lim=c(0, 1), add=T)
  }
  
  xpos <- xcir
  ypos <- rep(0, ncir)
  boxed.labels(xpos, ypos, c("180", "165", "150", "135", "120", "105", "90"), border = FALSE, cex=1)
  
  for(j in 1:length(ESS)){
    if(ESS[j]>=5 & grid_pt[3, j]<=0){
      pt <- grid_pt[, j]
      start <- pt[1:2]/sqrt(1+abs(pt[3]))
      
      if(signif[j]==TRUE){
        end <- start + arrow[, j]
        arrows(start[1], start[2], end[1], end[2], length = 0.05)
      }else{
        rad <- as.numeric(sqrt(start[1]^2+start[2]^2))
        angle <- as.numeric(sign(start[2])*acos(start[1]/rad))
        radial.plot(rad, angle, rp.type="s", point.symbols=16, add=T)
      }
    }
  }
  mtext(paste(expression(kappa), "=", toString(round(seq_kappa[i], digits=3), width=11)), side = 1, line = -5, cex = 1.4, outer = TRUE)
#  dev.off() #If you want to save as files, run this line along with other line
}
