# Test tunning with BIC and evaluate by MISE
# - Stop by first increased loss
# - weight: unsquared
# - sample size: 50

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
alpha_3 <- 1.5 * ((xx - 2)^2 + 2 * (yy - 0.2)^2)
alpha_4 <-  ((xx - 1)^2 + 2.5 * (yy - 1.2)^2)
alpha_mat[,1] <- alpha_1
alpha_mat[,2] <- alpha_2
#alpha_mat[,3] <- alpha_3
#alpha_mat[,4] <- alpha_3

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

#### sisVIVE ####
beta.sisv <- data.frame(matrix(ncol = 2, nrow = n_pos))
alpha.sisv <- data.frame(matrix(ncol = L, nrow = n_pos))
start.time <- Sys.time()
for (j in 1:n_pos) {
  sisv1 <- sisVIVE(Y[,j], D=D, Z=Z_matrix)
  lambda <- cv.sisVIVE(Y[,j], D=D, Z=Z_matrix, K = 10)$lambda
  sisv2 <- predict(sisv1, lambda, type = "coefficients")
  # sisv3 <- predict(sisv1, lambda, type = "instruments")
  # print(sisv3$instruments)
  beta.sisv[j,2] <- sisv2$beta
  alpha.sisv[j,] <- as.vector(sisv2$alpha)
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

#### Step 0: Initialize ####
# (1) triangulation BQ2
bb = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
VT = TriMesh(bb, 3)
V <- VT$V
Tr <- VT$Tr
d=2; r=1;
Ball <- basis(V, Tr, d, r, loc)
# B: 7505  108
# Q2 108  15
# K: 108 108
# BQ2: 7505 15
Q2 <- Ball$Q2 # Q2
B <- Ball$B   # basis function: dim- n by nT*{(d+1)(d+2)/2}
K <- Ball$K    # thin-plate energy function (Penalty)
BQ2 <- as.matrix(B%*%Q2)
# (2) alpha.sisv : 7505*10
start.time <- Sys.time()
nq <- ncol(Q2)
theta_alpha_init <- data.frame(matrix(ncol = L, nrow = nq))
theta_beta_init <- data.frame(matrix(ncol = 2, nrow = nq))
XtX <- crossprod(BQ2)
for (i in 1:10) {
  alpha_tilde <- alpha.sisv[,i]
  Xty <- crossprod(BQ2, alpha_tilde)   
  theta_alpha_init[,i] <- as.vector(solve(XtX, Xty))
}
Xty <- crossprod(BQ2, beta.sisv[,1])   
theta_beta_init[,1] <- as.vector(solve(XtX, Xty))
Xty <- crossprod(BQ2, beta.sisv[,2])   
theta_beta_init[,2] <- as.vector(solve(XtX, Xty))
end.time <- Sys.time()
end.time - start.time # 0.0236

#### Proposed + Stop + weighted group lasso ####
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values
theta_init <- cbind(theta_beta_init,theta_alpha_init)

PGD_stop <- function(Y, Z_matrix, D.est, theta_init, lambda_1, lambda_2, a=0.0001, 
                     num_iterations=30, stop="loss"){
  ## weight in adaptive lasso calculation
  init_l2norm <- sqrt(colSums(theta_init[, 3:12]^2))
  
  weight.vec <- 1/(init_l2norm + abs(rnorm(L, sd=0.001)))
  
  theta_tilde <- theta_init
  residual.matrix <- data.frame(matrix(ncol = n_pos, nrow = n))
  loss.vec <- rep(NA,num_iterations)
  for (k in 1:num_iterations) {
    theta_tilde.last <- theta_tilde
    theta_tilde_vec <- unlist(as.vector(theta_tilde))
    residual.matrix_last <- residual.matrix
    #### Step 1: Gradient ####
    # (1) Gradient 1
    grad1 <- rep(0, nq * (2+L))
    est_BQ2beta0 <- matrix(rep(BQ2 %*% theta_tilde[,1], each=n), nrow=n, byrow=TRUE)
    est_BQ2beta1 <- tcrossprod(D.est,BQ2 %*% theta_tilde[,2])
    est_ZBQ2alpha <- tcrossprod(Z_matrix,(BQ2 %*% as.matrix(theta_tilde[,3:12])))
    residual.matrix <- Y - est_BQ2beta0 - est_BQ2beta1-est_ZBQ2alpha
    for (i in 1:n) {
      residual <- residual.matrix[i,]
      innergrad <- cbind(BQ2, (D.est[i])*BQ2,
                         do.call(cbind,lapply(1:L, function(ell) Z_matrix[i, ell]*BQ2)))
      grad1 <- grad1 + crossprod(residual,innergrad)
    }
    # (2) Gradient 2
    P <- as.matrix(crossprod(Q2,K)%*%Q2) ## 15*15
    Dlam <- lambda_1*P
    grad2 <- c(Dlam %*% theta_tilde[,1], Dlam %*% theta_tilde[,2],rep(0,nq*L))
    theta_tilde_vec <- theta_tilde_vec + a*(grad1/(n*n_pos)-grad2)
    theta_tilde <- matrix(theta_tilde_vec, nrow = nq, ncol = L+2)
    norm_grp_alpha <- rep(NA, L)
    #### Step 2: Proximal ####
    for (ell in 1:L) {
      grp_alpha <- theta_tilde[,(ell+2)]
      norm_grp_alpha[ell] <- sqrt(sum(grp_alpha^2))
      shrinkage_factor <- max(1 - a * lambda_2 * weight.vec[ell]/norm_grp_alpha[ell], 0)
      theta_tilde[,(ell+2)] <- shrinkage_factor * grp_alpha
    }
    est_BQ2beta0 <- matrix(rep(BQ2 %*% theta_tilde[,1], each=n), nrow=n, byrow=TRUE)
    est_BQ2beta1 <- tcrossprod(D.est,BQ2 %*% theta_tilde[,2])
    est_ZBQ2alpha <- tcrossprod(Z_matrix,(BQ2 %*% theta_tilde[,3:(L+2)]))
    residual.matrix <- Y - est_BQ2beta0 - est_BQ2beta1-est_ZBQ2alpha
    term1 <- mean((residual.matrix)^2)
    term2 <- lambda_1*(crossprod(theta_tilde[,1],P)%*%theta_tilde[,1] +
                         crossprod(theta_tilde[,2],P)%*%theta_tilde[,2])
    term3 <- lambda_2*(weight.vec%*%norm_grp_alpha)
    print(k)
    print(paste("residual^2: ", round(mean((residual.matrix)^2),3)))
    loss.vec[k] <- term1 + term2 +term3
    diff <- sqrt(sum((theta_tilde-theta_tilde.last)^2))
    last_l2 <- sqrt(sum(theta_tilde.last))
    loss_diff.1 <- (loss.vec[k-1]-loss.vec[k])/loss.vec[k-1]
    loss_diff.2 <- (loss.vec[k-2]-loss.vec[k-1])/loss.vec[k-2]
    loss_diff.3 <- (loss.vec[k-3]-loss.vec[k-2])/loss.vec[k-3]
    print(paste("# of 0s in alpha:", sum(theta_tilde[,3:12] == 0)))
    print(paste("estimator_l2_diff =",round(diff, 5)))
    print(paste("last_estimator_l2 =",round(last_l2, 5)))
    if(k>3){
      print(paste0("loss_diff:", round(loss_diff.3,5),sep = " ",round(loss_diff.2,5),sep = " ",round(loss_diff.1,5)))
      print(paste0("loss(-3,-2,-1):", round(loss.vec[k-2],5),sep = " ",round(loss.vec[k-1],5),sep = " ",round(loss.vec[k],5)))
    }
    ### Stop 1: L2 of Difference between theta_tilde/L2 of last theta_tilde
    if (stop == "estimator"){
      if (diff < 0.005*last_l2){
        break
      }
    }
    ### Stop 2: Loss of difference/Last loss
    if (k>2 && stop == "loss"){
      loss_diff.0 <- (loss.vec[1]-loss.vec[2])/loss.vec[2]
      if (loss_diff.1/loss_diff.0 < 0.1){
        break
      }
    }
    ### Stop 3: reach local minimum
    if (k>1 && stop =="loss.min"){
      if(loss.vec[k]>loss.vec[k-1]){
        break
      }
    }
  }
  mse <- mean((residual.matrix)^2)
  list(theta_tilde=theta_tilde, residual=residual.matrix, loss=loss.vec, mse=mse)
}


## Tunning by BIC
PDG.bic <- function(PDG_stop, lambda_1, a, num_iter, nlambda2, type.bic){
  bic.vec <- c()
  if (type.bic==1){
    for (nl in 1:length(nlambda2)) {
      lamc <- nlambda2[nl]
      print(paste0("-----lambda = ",lamc, "-----"))
      m.test <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init, 
                         lambda_1=lambda_1, lambda_2=lamc, a=a, num_iterations = num_iter, stop ="loss.min")
      log.mse <- log(m.test$mse)
      theta_alpha <- m.test$theta_tilde[,3:12]
      n <- dim(m.test$residual)[1]
      n_pos <- dim(m.test$residual)[2]
      #n_basis <- dim(theta_alpha)[1]
      non0.alpha <- sum(theta_alpha == 0)
      #card.I <- sum(theta_alpha == 0)/n_basis
      p.alpha <- non0.alpha*log(n*n_pos)/(n*n_pos)
      bic.vec[nl] <- log.mse + p.alpha
    }
  }
  if (type.bic==2){
    for (nl in 1:length(nlambda2)) {
      lamc <- nlambda2[nl]
      print(paste0("-----lambda = ",lamc, "-----"))
      m.test <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init, 
                         lambda_1=lambda_1, lambda_2=lamc, a=a, num_iterations = num_iter, stop ="estimator")
      log.mse <- log(m.test$mse)
      theta_alpha <- m.test$theta_tilde[,3:12]
      n <- dim(m.test$residual)[1]
      n_pos <- dim(m.test$residual)[2]
      n_basis <- dim(theta_alpha)[1]
      p <- dim(theta_alpha)[2]
      #card.I <- sum(theta_alpha == 0)/n_basis
      non0.alpha <- sum(theta_alpha == 0)
      p.alpha <- non0.alpha*log(p*n_basis)*log(n*n_pos)/(2*n*n_pos)
      bic.vec[nl] <- log.mse + p.alpha
    }
  }
  if (type.bic==3){
    ### \eta unsolved
    eta <- 0.5
    for (nl in 1:length(nlambda2)) {
      lamc <- nlambda2[nl]
      print(paste0("-----lambda = ",lamc, "-----"))
      m.test <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init, 
                         lambda_1=lambda_1, lambda_2=lamc, a=a, num_iterations = num_iter, stop ="estimator")
      log.mse <- log(m.test$mse)
      theta_alpha <- m.test$theta_tilde[,3:12]
      n <- dim(m.test$residual)[1]
      n_pos <- dim(m.test$residual)[2]
      n_basis <- dim(theta_alpha)[1]
      p <- dim(theta_alpha)[2]
      M <- sum(theta_alpha == 0)/n_basis
      KMF <- p*n_basis
      p.alpha <- M/n*(log(n)+2*eta*log(KMF))
      bic.vec[nl] <- log.mse + p.alpha
      print(bic.vec[nl])
    }
  }
  j <- which.min(bic.vec)
  min.lamb <- nlambda2[j]
  list(opt.lambda=min.lamb, bic.vec = bic.vec)
}


#### lambda2 test ####
m.001 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                  lambda_1=10/(n_pos*n), lambda_2=0.001, a=0.2, num_iterations=150, 
                  stop ="loss.min")

# 0.2 - 81 steps - Loss: 16.43365
###  # of 0s in alpha 60:  (initial: 30)
loss.vec <- m.001$loss
plot(x=1:150,loss.vec, type = "l")


nfold = 5
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values
D_const.est <- cbind(rep(1,n),D.est)
D_const.est2 <- cbind(D_const.est, Z_matrix[,c(1,2)])
iter=1
lambda=10^(seq(-6,6,by=1))
Y = as.matrix(Y)
cv=cv.FDAimage(Y,D_const.est2,loc,V,Tr,d,r,lambda,nfold,iter)
lamc=cv$lamc
mfit0=fit.FDAimage.ho.full(Y,D_const.est2[, 1:2],loc,V,Tr,d,r,lamc)
theta_init_2 = matrix(0, ncol = 2 + ncol(Z_matrix), nrow = nq)
theta_init_2[, 1:4] = mfit0$theta.mtx
beta.oracle = mfit0$beta[[1]][, 1:2]

m.01 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                 lambda_1=10/(n_pos*n), lambda_2=0.01, a=0.2, num_iterations=150, 
                 stop ="loss.min")
# 0.2 - 78 steps - Loss: 16.49044
# "# of 0s in alpha: 120" (initial: 90)
loss.vec <- m.01$loss
plot(x=1:150,loss.vec, type = "l")


m.0.1 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                  lambda_1=10/(n_pos*n), lambda_2=0.1, a=0.2, num_iterations=150, 
                  stop ="loss.min")
# 0.2 -  steps - Loss: 16.80303
# "# of 0s in alpha: 120" (initial: 105)


m.1 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                lambda_1=10/(n_pos*n), lambda_2=1, a=0.1, num_iterations=150, 
                stop ="loss.min")
# 0.2 -  59 steps - Loss: 19.70179
# "# of 0s in alpha: 120" (initial: 120)
# 0.1 -  137 steps - Loss: 19.055
# "# of 0s in alpha: 120" (initial: 120)
loss.vec <- m.1$loss[1:10]
plot(x=1:10,loss.vec, type = "l")


m.10 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                 lambda_1=10/(n_pos*n), lambda_2=10, a=0.1, num_iterations=150, 
                 stop ="loss.min")
# 0.2 -  42 steps -  Loss: 39.23
# of 0s in alpha: 135" (initial 120)
# 0.1 - 132 steps - Loss: 32.279
# of 0s in alpha: 135" (initial 120)
loss.vec <- m.10$loss[30:42]
plot(x=30:42,loss.vec, type = "l", ylim = c(35,45))

m.100 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                  lambda_1=10/(n_pos*n), lambda_2=100, a=0.01, num_iterations=150, 
                  stop ="loss.min")
# 0.2 - 18 steps - Loss: 225.11714
# of 0s in alpha: 135" (init: 120)
# 0.1 - 73 steps - Loss: 124.4931
# of 0s in alpha: 150" (init: 120)
# 0.01 - 150 steps - Loss: 94.89121
# of 0s in alpha: 135 (init: 120) --- not convinced
loss.vec <- m.100$loss[1:18]
plot(x=1:18,loss.vec, type = "l",ylim = c(150,285))


m.1000 <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, theta_init=theta_init,  
                   lambda_1=10/(n_pos*n), lambda_2=1000, a=0.01, num_iterations=500, 
                   stop ="loss.min")
# 0.2 - 4 steps - Loss: 2157.87497
# 0.1 - 65 steps - Loss: 1042.75247
# 0.01 - 150 steps - Loss: 127.55319
# 0.001 - 150 steps - loss: 716.76694
# "# of 0s in alpha: 150" (init: 135)
loss.vec <- m.1000$loss
plot(x=1:length(loss.vec),loss.vec, type = "l")


#### BIC calc--separation ####
t1.bic <- function(m.test){
  log.mse <- log(m.test$mse)
  theta_alpha <- m.test$theta_tilde[,3:12]
  n <- dim(m.test$residual)[1]
  n_pos <- dim(m.test$residual)[2]
  #n_basis <- dim(theta_alpha)[1]
  non0.alpha <- sum(theta_alpha == 0)
  #card.I <- sum(theta_alpha == 0)/n_basis
  p.alpha <- non0.alpha*log(n*n_pos)/(n*n_pos)
  bic<- log.mse + p.alpha
  return(bic)
}
t2.bic <- function(m.test){
  log.mse <- log(m.test$mse)
  theta_alpha <- m.test$theta_tilde[,3:12]
  n <- dim(m.test$residual)[1]
  n_pos <- dim(m.test$residual)[2]
  n_basis <- dim(theta_alpha)[1]
  p <- dim(theta_alpha)[2]
  #card.I <- sum(theta_alpha == 0)/n_basis
  non0.alpha <- sum(theta_alpha == 0)
  p.alpha <- non0.alpha*log(p*n_basis)*log(n*n_pos)/(2*n*n_pos)
  bic <- log.mse + p.alpha
  return(bic)
}

t3.bic <- function(m.test, eta = 0.1){
  log.mse <- log(m.test$mse)
  theta_alpha <- m.test$theta_tilde[,3:12]
  n <- dim(m.test$residual)[1]
  n_pos <- dim(m.test$residual)[2]
  n_basis <- dim(theta_alpha)[1]
  p <- dim(theta_alpha)[2]
  card.I <- sum(theta_alpha != 0)/n_basis
  p.alpha <- card.I*(log(n) + 2*eta*log(p*n_basis))/n
  bic <- log.mse + p.alpha
  return(bic)
}


t1.bic(m.001) 
t1.bic(m.01) 
t1.bic(m.0.1) 
t1.bic(m.1)  
t1.bic(m.10) 
t1.bic(m.100) 
t1.bic(m.1000) 

t2.bic(m.001) 
t2.bic(m.01) 
t2.bic(m.0.1) 
t2.bic(m.1)  
t2.bic(m.10) 
t2.bic(m.100) 
t2.bic(m.1000) 

t3.bic(m.001) 
t3.bic(m.01) 
t3.bic(m.0.1) 
t3.bic(m.1)  
t3.bic(m.10) 
t3.bic(m.100) 
t3.bic(m.1000) 


### MISE
comp.mise<- function(mtest){
  theta_prop <- mtest$theta_tilde
  beta.prop <- data.frame(matrix(ncol = 2, nrow = n_pos))
  beta.prop[,1] <- BQ2 %*% theta_prop[,1]
  beta.prop[,2] <- BQ2 %*% theta_prop[,2]
  mise <- apply(beta.true - beta.prop, 2, function(x) mean(x^2))
  return(mise)
}

beta.init <- data.frame(matrix(ncol = 2, nrow = n_pos))
beta.init[,1] <- BQ2 %*% theta_init[,1]
beta.init[,2] <- BQ2 %*% theta_init[,2]
apply(beta.true - beta.init, 2, function(x) mean(x^2))
comp.mise(m.001) 
comp.mise(m.01) 
comp.mise(m.0.1) 
comp.mise(m.1)  
comp.mise(m.10) 
comp.mise(m.100) 
comp.mise(m.1000)

