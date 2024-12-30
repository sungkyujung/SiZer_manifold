library(MASS)
library(reshape2)
library(ggplot2)

############################################
##########Figure 1 & Movie 1################
############################################

library(ClusTorus)
library(plot3D)

x1 <- SARS_CoV_2$phi
x2 <- SARS_CoV_2$psi

plot(x1, x2, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi), pch=16)
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

set.seed(230329)

x1 <- SARS_CoV_2$phi
x2 <- SARS_CoV_2$psi
grid <- 24
kappa <- 2
res <- torus_sizer(x1, x2, grid, seq_kappa=kappa, B=10000, alpha=0.05, method=2)

grid_pt <- res[[1]]$grid_pt
arrow <- res[[1]]$arrow
ESS <- res[[1]]$ESS
signif <- res[[1]]$signif
color <- res[[1]]$color

plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))##############################################
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

idx <- 1
for(i in 0:(grid-1)){
  for(j in 0:(grid-1)){
    polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*i-1, 2*i+1, 2*i+1, 2*i-1)*pi/grid, col=color[idx], density=NA)
    idx <- idx+1
  }
}

for (j in 1:length(ESS)){
  if(ESS[j]>=5){
    if(signif[j]==TRUE){
      arrows(grid_pt[1, j], grid_pt[2, j], grid_pt[1, j]+arrow[1, j], grid_pt[2, j]+arrow[2, j], length = 0.04, lwd=2)######################################
    }else{
      points(grid_pt[1, j], grid_pt[2, j], pch=20)
    }
  }
}

i <- 1
pt <- res[[i]]$grid_pt
grid_pt_x <- pt[1, ]
grid_pt_y <- pt[2, ]
kde <- res[[i]]$kde
color <- res[[i]]$color
color_adj <- adjustcolor(color, alpha=0.25)

surf3D(matrix(grid_pt_x, nrow=grid, byrow=T), matrix(grid_pt_y, nrow=grid, byrow=T), matrix(kde, nrow=grid, byrow=T), 
       colvar = matrix(color_adj, nrow=grid, byrow=T), colkey = TRUE, box = TRUE, bty = "b", phi = 40, theta = -40)


############################################
##########Movie 3###########################
############################################

library(MASS)
library(tidyverse) 
library(circular)
library(bio3d)

pdb <- read.pdb("9ENJ") 
a <- torsion.pdb(pdb)

data <- cbind(a$phi/180*pi,a$psi/180*pi)
data <- data[-which(is.na(data[,1])|is.na(data[,2])),]

x1 <- data[, 1]+2*pi*(sign(data[, 1])==-1)
x2 <- data[, 2]+2*pi*(sign(data[, 2])==-1)

plot(x1, x2, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi), pch=16)
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

grid <- 24
seq_kappa <- 10^seq(0, 2, length=200)
res <- torus_sizer(x1, x2, grid, seq_kappa, B=10000, alpha=0.05, method=2)

for(i in 1:length(res)){
  grid_pt <- res[[i]]$grid_pt
  arrow <- res[[i]]$arrow
  ESS <- res[[i]]$ESS
  signif <- res[[i]]$signif
  color <- res[[i]]$color
  
  png(paste("Rplot_", i, ".png"), width=800, height=800) #If you want to save as files, run this line along with other line
  par(mar = c(7, 4, 2, 2))
  
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
  mtext(paste(expression(kappa), "=", toString(round(seq_kappa[i], digits=3), width=11)), side = 1, line = -2.5, cex = 1.4, outer = TRUE)
  dev.off() #If you want to save as files, run this line along with other line
}

############################################
##########Figure 4 & Movie 4-5##############
############################################

set.seed(230328)

grid <- 24

seq_kappa <- 4
#seq_kappa <- 10^seq(0, 2, length=200)

x <- runif(1000, 0, 2*pi)
y <- runif(1000, 0, 2*pi)

res <- torus_sizer(x, y, grid, seq_kappa, B=10000, alpha=0.05, method=2)

signif <- res[[1]]$signif
signif[res[[1]]$signif==FALSE] <- "purple"
signif[res[[1]]$signif==TRUE] <- "red"

plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

idx <- 1
for(i in 0:(grid-1)){
  for(j in 0:(grid-1)){
    polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*i-1, 2*i+1, 2*i+1, 2*i-1)*pi/grid, col=signif[idx], density=NA)
    idx <- idx+1
  }
}

grid_pt <- res[[1]]$grid_pt
arrow <- res[[1]]$arrow
ESS <- res[[1]]$ESS
signif <- res[[1]]$signif

for (j in 1:length(ESS)){
  if(ESS[j]>=5){
    if(signif[j]==TRUE){
      arrows(grid_pt[1, j], grid_pt[2, j], grid_pt[1, j]+arrow[1, j], grid_pt[2, j]+arrow[2, j], length = 0.04, lwd=2)
    }else{
      points(grid_pt[1, j], grid_pt[2, j], pch=20)
    }
  }
}

color <- res[[1]]$color

plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

idx <- 1
for(i in 0:(grid-1)){
  for(j in 0:(grid-1)){
    polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*i-1, 2*i+1, 2*i+1, 2*i-1)*pi/grid, col=color[idx], density=NA)
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

############################################
##########Figure 5 & Movie 6-7##############
############################################

set.seed(230329)

grid <- 24
kappa <- 4

rsamp <- runif(2000, 0, 2*pi)
dens <- runif(2000, 0, 1)

idx <- (rsamp<=3*pi/10 & dens<=0.2) | (rsamp>3*pi/10 & rsamp<=7*pi/10 & dens<2/pi*rsamp-0.4) | (rsamp>7*pi/10 & rsamp<=13*pi/10) | (rsamp>13*pi/10 & rsamp<=17*pi/10 & dens+2/pi*rsamp<3.6) | (rsamp>17*pi/10 & dens<=0.4)
x <- (rsamp[idx])[1:1000]
y <- runif(1000, 0, 2*pi)

res <- torus_sizer(x, y, grid, seq_kappa=kappa, B=10000, alpha=0.05, method=2)

signif <- res[[1]]$signif
signif[res[[1]]$signif==FALSE] <- "purple"
signif[res[[1]]$signif==TRUE] <- "red"

plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

idx <- 1
for(i in 0:(grid-1)){
  for(j in 0:(grid-1)){
    polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*i-1, 2*i+1, 2*i+1, 2*i-1)*pi/grid, col=signif[idx], density=NA)
    idx <- idx+1
  }
}

grid_pt <- res[[1]]$grid_pt
arrow <- res[[1]]$arrow
ESS <- res[[1]]$ESS
signif <- res[[1]]$signif

for (j in 1:length(ESS)){
  if(ESS[j]>=5){
    if(signif[j]==TRUE){
      arrows(grid_pt[1, j], grid_pt[2, j], grid_pt[1, j]+arrow[1, j], grid_pt[2, j]+arrow[2, j], length = 0.04, lwd=2)
    }else{
      points(grid_pt[1, j], grid_pt[2, j], pch=20)
    }
  }
}

color <- res[[1]]$color

plot(NULL, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xaxt='n', yaxt='n', xlab = expression(phi), ylab = expression(psi))
axis(1, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
axis(2, at=seq(0, 2*pi, length.out=5), labels=c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))

idx <- 1
for(i in 0:(grid-1)){
  for(j in 0:(grid-1)){
    polygon(c(2*j-1, 2*j-1, 2*j+1, 2*j+1)*pi/grid, c(2*i-1, 2*i+1, 2*i+1, 2*i-1)*pi/grid, col=color[idx], density=NA)
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

############################################
#####Statistics related to Figure 4, 5######
############################################


set.seed(230328)

grid <- 24
kappa <- 4

color_mat <- NULL
signif_mat <- NULL

for(rep in 1:100){
  phi <- runif(1000, 0, 2*pi)
  psi <- runif(1000, 0, 2*pi)
  
  res <- torus_sizer(phi, psi, grid, seq_kappa=kappa, B=10000, alpha=0.05, method=2)
  
  color_err <- NULL
  signif_err <- NULL
  
  for(i in 1:length(res)){
    color <- res[[i]]$color
    signif <- res[[i]]$signif
    
    color_err <- c(color_err, mean(color!='green'))
    signif_err <- c(signif_err, mean(signif=='TRUE'))
  }
  
  color_mat <- rbind(color_mat, color_err)
  signif_mat <- rbind(signif_mat, signif_err)
}

col_m <- mean(color_mat)
col_sd <- sd(color_mat)

sig_m <- mean(signif_mat)
sig_sd <- sd(signif_mat)

sig_m
sig_sd
col_m
col_sd



set.seed(230329)

grid <- 24
kappa <- 4

color_mat <- NULL
signif_pmat <- NULL
signif_mat <- NULL

for(rep in 1:100){
  rsamp <- runif(2000, 0, 2*pi)
  dens <- runif(2000, 0, 1)
  
  idx <- (rsamp<=pi/3 & dens<=0.4) | (rsamp>pi/3 & rsamp<=2*pi/3 & dens<1.8/pi*rsamp-0.2) | (rsamp>2*pi/3 & rsamp<=4*pi/3) | (rsamp>4*pi/3 & rsamp<=5*pi/3 & dens+1.8/pi*rsamp<3.4) | (rsamp>5*pi/3 & dens<=0.4)
  phi <- (rsamp[idx])[1:1000]
  psi <- runif(1000, 0, 2*pi)
  
  res <- torus_sizer(phi, psi, grid, seq_kappa=kappa, B=10000, alpha=0.05, method=2)
  
  color_err <- NULL
  signif_power <- NULL
  signif_err <- NULL
  
  for(i in 1:length(res)){
    grid_x <- res[[i]]$grid_pt[1,]
    color <- res[[i]]$color
    signif <- res[[i]]$signif
    idx <- (grid_x>0.3*pi & grid_x<0.7*pi) | (grid_x>1.3*pi & grid_x<1.7*pi)
    
    color_err <- c(color_err, mean(color!='green'))
    signif_power <- c(signif_power, mean(signif[idx]=='TRUE'))
    signif_err <- c(signif_err, mean(signif[-which(idx)]=='TRUE'))
  }
  
  color_mat <- rbind(color_mat, color_err)
  signif_pmat <- rbind(signif_pmat, signif_power)
  signif_mat <- rbind(signif_mat, signif_err)
}

col_m <- mean(color_mat)
col_sd <- sd(color_mat)

sig_m <- mean(signif_mat)
sig_sd <- sd(signif_mat)

sig_p_m <- mean(signif_pmat)
sig_p_sd <- sd(signif_pmat)

sig_m
sig_sd
col_m
col_sd
sig_p_m
sig_p_sd

############################################
##########Figure 6 & Movie 8################
############################################

set.seed(230330)
n <- 1000

theta <- rnorm(10*n, mean = pi/3, sd = 0.1)
theta <- theta[theta>=0]

theta_len <- length(theta)
rv_unif <- runif(theta_len, min = 0, max = 2*pi)

theta <- theta[sin(theta)>=rv_unif]
theta <- theta[1:n]

phi <- runif(n, min = 0, max = 2*pi)
data <- rbind(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))

seq_kappa <- 20
#seq_kappa <- 10^seq(0, 2, length=200)
res <- sphere_sizer(data, grid_num=800, seq_kappa, B=10000, alpha=0.05, method=2)



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

for(i in 1:length(res)){   #show on viewer
  kappa <- res[[i]]$kappa
  grid_pt <- res[[i]]$grid_pt
  arrow <- res[[i]]$arrow
  color <- res[[i]]$color
  ESS <- res[[i]]$ESS
  signif <- res[[i]]$signif
  
  png(paste("Rplot_", i, ".png"), width=1200, height=600) #save as file
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
  dev.off()  #save as file
}

############################################
##########Figure 7 & Movie 9################
############################################

set.seed(230330)

data_eq <- read.csv('query.csv')
lat <- data_eq[data_eq$lon>=137 | data_eq$lon<=-105, ]$latitude/180*pi
lon <- data_eq[data_eq$lon>=137 | data_eq$lon<=-105, ]$longitude/180*pi
data <- rbind(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
data <- data[, 1:1000]

seq_kappa <- c(12, 13, 14, 15, 16)

res <- sphere_sizer(data, grid_num=800, seq_kappa, B=1000, alpha=0.05, method=2)

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

for(i in 1:length(res)){   #show on viewer
  kappa <- res[[i]]$kappa
  grid_pt <- res[[i]]$grid_pt
  arrow <- res[[i]]$arrow
  color <- res[[i]]$color
  ESS <- res[[i]]$ESS
  signif <- res[[i]]$signif
  
  png(paste("Rplot_", i, ".png"), width=1200, height=600) #save as file
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
  dev.off()  #save as file
}

############################################
##########Figure 8(b)#######################
############################################

set.seed(230329)
x <- mvrnorm(500, c(1, 1), matrix(c(1, 0.5, 0.5, 1), ncol=2))
exp_x <- exp(-(x[, 1]^2+x[, 2]^2)/2)/(2*pi)
z1 <- exp_x*(x[, 1]^2-1)
z2 <- exp_x*(x[, 2]^2-1)
z3 <- exp_x*(x[, 1]*x[, 2])
d2h_1 <- cbind(z1, z2, z3)
colnames(d2h_1) <- c("x1_11", "x1_22", "x1_12")
lam_mat <- var(d2h_1)

x_star <- x-0.2
exp_x <- exp(-(x_star[, 1]^2+x_star[, 2]^2)/2)/(2*pi)
z1 <- exp_x*(x_star[, 1]^2-1)
z2 <- exp_x*(x_star[, 2]^2-1)
z3 <- exp_x*(x_star[, 1]*x_star[, 2])
d2h_2 <- cbind(z1, z2, z3)
colnames(d2h_2) <- c("x2_11", "x2_22", "x2_12")

x_star <- x-0.4
exp_x <- exp(-(x_star[, 1]^2+x_star[, 2]^2)/2)/(2*pi)
z1 <- exp_x*(x_star[, 1]^2-1)
z2 <- exp_x*(x_star[, 2]^2-1)
z3 <- exp_x*(x_star[, 1]*x_star[, 2])
d2h_3 <- cbind(z1, z2, z3)
colnames(d2h_3) <- c("x3_11", "x3_22", "x3_12")

x_star <- x-0.6
exp_x <- exp(-(x_star[, 1]^2+x_star[, 2]^2)/2)/(2*pi)
z1 <- exp_x*(x_star[, 1]^2-1)
z2 <- exp_x*(x_star[, 2]^2-1)
z3 <- exp_x*(x_star[, 1]*x_star[, 2])
d2h_4 <- cbind(z1, z2, z3)
colnames(d2h_4) <- c("x4_11", "x4_22", "x4_12")

d2h_full <- cbind(d2h_1, d2h_2, d2h_3, d2h_4)
cor_melt <- melt(cor(d2h_full))
ggplot(data = cor_melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low = "red", mid='white', high = "blue") 

x <- mvrnorm(500, c(0, 0, 0), lam_mat)

eigenvalues <- NULL
for(i in 1:500){
  hessian <- c(x[i, 1], x[i, 2], x[i, 2], x[i, 3])
  hessian <- matrix(hessian, nrow=2)
  res <- eigen(hessian, symmetric = T, only.values = T)
  eigenvalues <- rbind(eigenvalues, res$values)
}

plot(eigenvalues[, 1], eigenvalues[, 2], pch=16, xlab=expression(lambda^+(x)), ylab=expression(lambda^-(x)))
abline(0, 1)

