## Simulation test for 4 comparison models
## No Proposed model code included

library(FDAimage)
library(MASS)
library(matrixcalc)
library(ivreg)
library(sisVIVE)
library(fields)

set.seed(12345)

#### Data Generation ####
### Setting 1: 
#' (1) endogeneity: sigma_star = 0.3
#' (2) number of invalid IVs: s=2
#' (3) Z_i, Z_j are independent
#' (4) overall strength: 126
#' (5) relative strength: Equal
#' (6) sample size n = 50
#' (7) lambda_1 = 0.1; lambda2 = 0.02


# settings
n <- 50
sigma_star <- 0.3
lambda1 <- 0.1; lambda2 <- 0.02

L=10

# Location information
n1=95
n2=79
u1=seq(0,1,length.out=n1)
v1=seq(0,1,length.out=n2)
uu=rep(u1,each=n2)
vv=rep(v1,times=n1)
loc=as.matrix(cbind(uu,vv))
xx <- loc[, 1]
yy <- loc[, 2]
n_pos <- n1*n2

ind=inVT(V1,Tr1,loc[,1],loc[,2])
ind.inside=ind$ind.inside
length(ind.inside)

# coefficient: alpha, beta 
beta0 <- -3*xx^3 + 5*yy^3
beta1 <- (15 * ((xx - 0.5)^2 + (yy - 0.5)^2))  # 7505*1
beta.true <- cbind(beta0, beta1)
alpha_mat <- matrix(0, ncol = L, nrow = n_pos)  # 7505*10
alpha_1 <- 5 * ((xx - 0.5)^2 + (yy - 0.5)^2)  
alpha_2 <- 3 * ((xx - 1)^2 + (yy - 1)^2)
alpha_mat[,1] <- alpha_1
alpha_mat[,2] <- alpha_2

# coefficient: gamma.true (11)
# (5) relative strength = equal 0.5
gamma <- rep(0.5,L)
gamma0 <- 2
gamma.true <- c(gamma0, gamma)

# error ~ MVN()
mu_e <- c(0, 0, 0)
error_cov <- matrix(c(1, 0, sigma_star, 
                      0, 1, sigma_star, 
                      sigma_star, sigma_star, 1), nrow = 3)
is.positive.definite(error_cov)
samples <- mvrnorm(n = n, mu =  mu_e, Sigma = error_cov)
zeta_1 <- samples[, 1] # 50*1
zeta_2 <- samples[, 2] 
xi <- samples[, 3]   # 50*1

# IV: Z_mat_const (50x11)
mu_Z <- rep(0, L)
Z_matrix <- mvrnorm(n = n, mu = mu_Z, Sigma = diag(1, L)) #50x10
Z_mat_const <- cbind(rep(1,L),Z_matrix) #50x11

# Exposure: D (50)
D <- as.vector(Z_mat_const%*%gamma.true)+xi
# overall strength
# t(D)%*%D/2/sigma_star^2

# components in eta
lamda_zeta1 <- sqrt(lambda1)*zeta_1  # 50*1
lamda_zeta2 <- sqrt(lambda2)*zeta_2
phi1 <- 0.56 * sin(2 * pi * xx/n1)  # 7505*1
phi2 <- 0.61 * cos(2 * pi * yy/n2)

# individual level data
Y <- c()
c.eff <- c()
iv.eff <- c()
noise <- c()
eta.true <- c()
for (i in 1:n) {
  # error_i(s) ~ N(0,1)
  err.i <- rnorm(n=n_pos,mean = 0, sd=1)
  eta.i <- lamda_zeta1[i] * phi1 + lamda_zeta2[i] * phi2
  c.eff.i <- D[i]*beta.true[,2]
  iv.eff.i <- as.vector(alpha_mat %*% Z_matrix[i, ])
  noise.i <- eta.i + err.i
  Yi <- beta.true[,1] + c.eff.i + iv.eff.i + noise.i
  Y <- rbind(Y, Yi)  # 50*7505
  c.eff <- rbind(c.eff, c.eff.i)
  iv.eff <- rbind(iv.eff, iv.eff.i)
  noise <- rbind(noise, noise.i)
  eta.true <- rbind(eta.true, eta.i)
}


#### 2SLS ####
start.time <- Sys.time()
SE_list_iter <- c()
beta.2sls <- data.frame(matrix(ncol = 2, nrow = n_pos))
for (j in 1:n_pos) {
  tsls <- ivreg(Y[,j] ~ D | 
                  Z_matrix[,1]+ Z_matrix[,2]+ Z_matrix[,3]+ Z_matrix[,4]+ Z_matrix[,5]+ 
                  Z_matrix[,6]+ Z_matrix[,7]+ Z_matrix[,8]+ Z_matrix[,9]+ Z_matrix[,10])
  SE_list_iter <- c(SE_list_iter, (tsls$coefficients[2]-beta.true[j,2])^2)
  beta.2sls[j,1] <- tsls$coefficients[1]
  beta.2sls[j,2] <- tsls$coefficients[2]
}
end.time <- Sys.time()
end.time - start.time
mean(SE_list_iter)


#### sisVIVE ####
beta.sisv <- data.frame(matrix(ncol = 2, nrow = n_pos))
start.time <- Sys.time()
for (j in 1:n_pos) {
  sisv1 <- sisVIVE(Y[,j], D=D, Z=Z_matrix)
  lambda <- cv.sisVIVE(Y[,j], D=D, Z=Z_matrix, K = 10)$lambda
  sisv2 <- predict(sisv1, lambda, type = "coefficients")
  beta.sisv[j,2] <- sisv2$beta
  # recover beta_0
  alpha.hat <- sisv2$alpha
  beta.hat <- sisv2$beta
  D.bar <- mean(D)
  Z.bar <- colMeans(Z_matrix)
  normZ <- sqrt(colSums(Z_matrix^2))
  Y.bar <- mean(Y[,j])
  beta.sisv[j,1] <- Y.bar-sum((Z.bar * normZ / sqrt(n - 1))* alpha.hat)-D.bar*beta.hat
} 
end.time <- Sys.time()
end.time - start.time

#### BPST ####
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values
D_const.est <- cbind(rep(1,n),D.est)
## old
#data(V1); data(Tr1); # rectangular domain
#data(brain_boundary); # brain imaging boundary of slide 48;
#V=V1; Tr=Tr1;
bb = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
VT = TriMesh(bb, 3)
V <- VT$V
Tr <- VT$Tr
d=2; r=1; #np=3; rho=0.5; 
nfold=10; #alpha0=0.05;
iter=1
lambda=10^(seq(-6,6,by=1))
cv=cv.FDAimage(Y,D_const.est,loc,V,Tr,d,r,lambda,nfold,iter)
lamc=cv$lamc
# Image-on-scalar regression;
mfit0=fit.FDAimage(Y,D_const.est,loc,V,Tr,d,r,lamc)
beta.bpst=mfit0$beta

#### BPST-oracle ####
start.time <- Sys.time()
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values
D_const.est <- cbind(rep(1,n),D.est)
bb = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
VT = TriMesh(bb, 3)
V <- VT$V
Tr <- VT$Tr
d=2; r=1; #np=3; rho=0.5; 
nfold=10; #alpha0=0.05;
iter=1
lambda=10^(seq(-6,6,by=1))
## add invalid effect here:
D_const.est2 <- cbind(D_const.est, Z_matrix[,c(1,2)])
cv=cv.FDAimage(Y,D_const.est2,loc,V,Tr,d,r,lambda,nfold,iter)
lamc=cv$lamc
# Image-on-scalar regression;
mfit0=fit.FDAimage(Y,D_const.est2,loc,V,Tr,d,r,lamc)
beta.bpor=mfit0$beta
end.time <- Sys.time()
end.time - start.time





#### Plot ####
rect.plot <- function(betai,zlim){
  z1 <-u1; z2 <- v1;
  n1 <- length(z1); n2 <- length(z2);
  betai.mtx <- matrix(betai,ncol=n1,nrow=n2)
  image.plot(z2, z1, betai.mtx, zlim = zlim, legend.shrink = 0.8, legend.width = 1.2)
}

setwd("/Users/yitingwang/Desktop/ISR_IV/code/sim_plot")

# true beta_0
pdf("b0_true.pdf", width = 7, height = 7)
rect.plot(beta.true[,1],c(-5,5))
dev.off()
# 2SLS - beta_0 
pdf("b0_2sls.pdf", width = 7, height = 7)
rect.plot(beta.2sls[,1],c(-5,5))
dev.off()
# sisv - beta_0
pdf("b0_sisv.pdf", width = 7, height = 7)
rect.plot(beta.sisv[,1],c(-5,5))
dev.off()
# BPST - beta_0 
pdf("b0_bpst.pdf", width = 7, height = 7)
rect.plot(beta.bpst[,1],c(-5,5))
dev.off()
# BPST_ora -beta_1
pdf("b0_bpor.pdf", width = 7, height = 7)
rect.plot(beta.bpor[,1],c(-5,5))
dev.off()

# true beta_1
pdf("b1_true.pdf", width = 7, height = 7)
rect.plot(beta.true[,2],c(-0.5,9))
dev.off()
# 2SLS - beta_1
pdf("b1_2sls.pdf", width = 7, height = 7)
rect.plot(beta.2sls[,2],c(-0.5,9))
dev.off()
# sisVIVE - beta_1
pdf("b1_sisv.pdf", width = 7, height = 7)
rect.plot(beta.sisv[,2],c(-0.5,9))
dev.off()
# BPST - beta_1
pdf("b1_bpst.pdf", width = 7, height = 7)
rect.plot(beta.bpst[,2],c(-0.5,9))
dev.off()
# BPST_ora -beta_1
pdf("b1_bpor.pdf", width = 7, height = 7)
rect.plot(beta.bpor[,2],c(-0.5,9))
dev.off()









