# Test tunning with BIC and evaluate by MISE
# - Stop by first increased loss
# - weight: unsquared
# - sample size: 50


#### Note: BIC does not work for invalid = 4, eta=(0.3 0.5)

library(FDAimage)
library(MASS)
library(matrixcalc)
library(ivreg)
library(sisVIVE)
#library(fields)

#### Data Generation ####
### Setting 1: 
#' (1) endogeneity: sigma_star = 0.3
#' (2) number of invalid IVs: s=2
#' (3) Z_i, Z_j are independent
#' (4) overall strength: 126
#' (5) relative strength: Equal
#' (6) sample size n = 50
#' (7) lambda_1 = 0.1; lambda2 = 0.02

args <- commandArgs(trailingOnly = TRUE)
seed<-as.numeric(args[1])
set.seed(seed)

# settings
n <- 200
sigma_star <- 0.01
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
alpha_mat[,3] <- alpha_3
alpha_mat[,4] <- alpha_4

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
  L = ncol(Z_matrix)
  X = cbind(rep(1,n),D.est,Z_matrix)
  ## weight in adaptive lasso calculation
  init_l2norm <- sqrt(colSums(theta_init[, -(1:2)]^2))
  
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
    Yhat <- tcrossprod(X,(BQ2 %*% as.matrix(theta_tilde)))
    residual.matrix <- Y - Yhat
    grad1 <- rowSums(sapply(1:n, function(iter) 
      as.matrix(kronecker(X[iter, ], crossprod(BQ2, residual.matrix[iter, ])))))
    # (2) Gradient 2
    P <- as.matrix(crossprod(Q2,K)%*%Q2) ## 15*15
    Dlam <- lambda_1*P
    grad2 <- c(Dlam %*% theta_tilde[,1], Dlam %*% theta_tilde[,2],rep(0,nq*L))
    theta_tilde_vec <- theta_tilde_vec + 2*a*(grad1/(n*n_pos)-grad2)
    cat("gradient norm: ", sqrt(sum((grad1/(n*n_pos)-grad2)^2)), "\n")
    theta_tilde <- matrix(theta_tilde_vec, nrow = nq, ncol = L+2)
    norm_grp_alpha <- rep(NA, L)
    #### Step 2: Proximal ####
    for (ell in 1:L) {
      grp_alpha <- theta_tilde[,(ell+2)]
      norm_grp_alpha[ell] <- sqrt(sum(grp_alpha^2))
      shrinkage_factor <- max(1 - a * lambda_2 * weight.vec[ell]/norm_grp_alpha[ell], 0)
      theta_tilde[,(ell+2)] <- shrinkage_factor * grp_alpha
    }
    residual.matrix <- Y - tcrossprod(X,(BQ2 %*% as.matrix(theta_tilde)))
    term1 <- mean((residual.matrix)^2)
    term2 <- lambda_1*(crossprod(theta_tilde[,1],P)%*%theta_tilde[,1] +
                         crossprod(theta_tilde[,2],P)%*%theta_tilde[,2])
    term3 <- lambda_2*(weight.vec%*%norm_grp_alpha)
    print(k)
    print(paste("residual^2: ", round(mean((residual.matrix)^2),3)))
    loss.vec[k] <- term1 + term2 + term3
    
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

### MISE
comp.mise<- function(mtest){
  theta_prop <- mtest$theta_tilde
  beta.prop <- data.frame(matrix(ncol = 2, nrow = n_pos))
  beta.prop[,1] <- BQ2 %*% theta_prop[,1]
  beta.prop[,2] <- BQ2 %*% theta_prop[,2]
  mise <- apply(beta.true - beta.prop, 2, function(x) mean(x^2))
  return(mise)
}

#### lambda2 selection ####
lambda_2_all <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
bic_all <- rep(NA, length(lambda_2_all))
mise_all <- matrix(NA, length(lambda_2_all),2)
invalid_set <- list()

for (iter in 1:length(lambda_2_all)){
  m.fit <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, 
                    theta_init=theta_init,  
                    lambda_1=10/(n_pos*n), lambda_2=lambda_2_all[iter], 
                    a=5, num_iterations=500, 
                    stop ="loss.min")
  theta_norm <- apply(m.fit$theta_tilde[, -c(1:2)], 2, 
                      function(x) sum(x)) 
  invalid_set[[iter]] <- which(theta_norm != 0)
  mise_all[iter,] <- comp.mise(m.fit)
  bic_all[iter] <- t3.bic(m.fit)
  print(iter)
}

ind.min <- which.min(bic_all)
lambda_2<- lambda_2_all[ind.min]

#### variance calculation ####
m.fit <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, 
                  theta_init=theta_init,  
                  lambda_1=10/(n_pos*n), lambda_2=lambda_2, 
                  a=5, num_iterations=500, 
                  stop ="loss.min")

### select invalid effects
theta_norm <- apply(m.fit$theta_tilde[, -c(1:2)], 2, 
                    function(x) sum(x)) 
Z_orac <- Z_matrix[, theta_norm != 0]
a <- dim(Z_orac)[2]

### estimate D.hat
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values
D_const.est <- cbind(rep(1,n),D.est)
X <- cbind(D_const.est, Z_orac)

### \Gamma^{-1}*(n*N)
WW <- kronecker(crossprod(X),crossprod(BQ2))
P <- as.matrix(crossprod(Q2,K)%*%Q2) ## 15*15
math.D <- kronecker(diag( c(1,1,rep(0,a)) ),P)
lamc <- 10/(n_pos*n) 
Dlam <- math.D*lamc
WD.inv <- chol2inv(chol(WW+Dlam)) 

VV <- crossprod(BQ2) ## 15*15
VV.inv <- chol2inv(chol(VV))
R <- m.fit$residual # 50*7505
M <- lapply(1:n, function(iter) WD.inv%*%kronecker(X[iter,],VV)%*%VV.inv) # list 50: 90*15
# 50 individuals with 90 by 15 matrix
eta.theta <- t(R%*%BQ2) # 15 by 50
A <- lapply(1:n, function(iter) tcrossprod(M[[iter]]%*%eta.theta)) # list 50: 90*90
A <- Reduce("+",A)
A <- as.matrix(A)

## multiply \mathbb{B}
Geta <- c()
J <- ncol(BQ2)
nx <- 2+a
for(iter in 1:nx){
  a.begin <- (iter-1)*J+1
  a.end <- iter*J
  A.temp <- A[a.begin:a.end,a.begin:a.end]
  Geta <- cbind(Geta,apply(BQ2,1,function(x) x%*%A.temp%*%x))
}
Sigma <- Geta/n

est.var.beta <- Sigma[1:(n_pos*2)]

save.image(paste0("/home/agy4yy/ISR_IV/output/test_",seed,".RData"))

