##############################################
#run file 'torus_function.R' first############
##############################################
#x: 1st coordiante of torus-valued data (vector)
#y: 2nd coordiante of torus-valued data (vector / same length with x)
#grid: the number of the grid of axis (scalar) / generate equally spaced grid of 'grid times grid'
#seq_kappa: the sequence of the value of kappa (vector)
#B: bootstrap repeatation
#alpha: significance level
#method: 1(Gaussian), 2(simul. Gaussian), 3(simul. bootstrap for grid), 4(simul. bootstrap for grid, kappa)
##############################################
#Line 15-17 is an example: put your data here#
##############################################

x1 <- runif(1000, 0, 2*pi)
x2 <- runif(1000, 0, 2*pi)
seq_kappa <- 10^seq(0, 2, length=200)
res <- torus_sizer(x1, x2, grid=24, seq_kappa, B=10000, alpha=0.05, method=2)

##############################################
#Do not edit below this#######################
##############################################
#figure of gradient###########################
##############################################
for(i in 1:length(res)){
  signif <- res[[i]]$signif
  signif[res[[i]]$signif==FALSE] <- "purple"
  signif[res[[i]]$signif==TRUE] <- "red"
  
#  png(paste("Rplot_", i, ".png"), width=800, height=800) #If you want to save as files, run this line along with other line
  
  plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))
  axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
  axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
  
  idx <- 1
  for(k in 0:(grid-1)){
    for(j in 0:(grid-1)){
      polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*k-1, 2*k+1, 2*k+1, 2*k-1)*pi/grid, col=signif[idx], density=NA)
      idx <- idx+1
    }
  }
  
  grid_pt <- res[[i]]$grid_pt
  arrow <- res[[i]]$arrow
  ESS <- res[[i]]$ESS
  signif <- res[[i]]$signif

  for (j in 1:length(ESS)){
    if(ESS[j]>=5){
      if(signif[j]==TRUE){
        arrows(grid_pt[1, j], grid_pt[2, j], grid_pt[1, j]+arrow[1, j], grid_pt[2, j]+arrow[2, j], length = 0.04, lwd=2)
      }else{
        points(grid_pt[1, j], grid_pt[2, j], pch=20)
      }
    }
  }
#  dev.off() #If you want to save as files, run this line along with other line
}

##############################################
#figure of Hessian############################
##############################################

for(i in 1:length(res)){
  color <- res[[i]]$color

#  png(paste("Rplot_", i, ".png"), width=800, height=800) #If you want to save as files, run this line along with other line
  
  plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))
  axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
  axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
  
  idx <- 1
  for(k in 0:(grid-1)){
    for(j in 0:(grid-1)){
      polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*k-1, 2*k+1, 2*k+1, 2*k-1)*pi/grid, col=color[idx], density=NA)
      idx <- idx+1
    }
  }
  
  for (j in 1:length(ESS)){
    if(ESS[j]>=5){
      if(signif[j]==TRUE){
        arrows(grid_pt[1, j], grid_pt[2, j], grid_pt[1, j]+arrow[1, j], grid_pt[2, j]+arrow[2, j], length = 0.04, lwd=2)
      }else{
        points(grid_pt[1, j], grid_pt[2, j], pch=20)
      }
    }
  }
  #  dev.off() #If you want to save as files, run this line along with other line
}
