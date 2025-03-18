#library(fields)
library(FDAimage)
library(dplyr)

dir_out <- "/home/agy4yy/ISR_IV/output/"
load(paste0(dir_out,"n50_inv2_snr1_seed1.RData")) ### can be found in "Data" folder

### start
P <- as.matrix(crossprod(Q2,K)%*%Q2) ## 15*15
theta_prop <- m.fit$theta_tilde
beta_hat <- BQ2 %*% theta_prop[,2]
theta_norm <- apply(m.fit$theta_tilde[, -c(1:2)], 2, 
                    function(x) sum(x)) 
Z_sub <- Z_matrix[, theta_norm != 0]

#### Calculate sigma_hat(s_j) ####
### ATTENTION: change lamc as an input
comp_sigmahat <- function(D, Z_matrix, Z_orac, m){
  a <- dim(Z_orac)[2]
  lamc <- 10
  ### estimate D.hat
  d.fit <- lm(D~Z_matrix)
  D.est <- d.fit$fitted.values
  
  Z_tilde <- cbind(rep(1,n), Z_orac)
  bb <- crossprod(BQ2)
  ddbb <- kronecker(crossprod(D.est),bb)
  lam1D <- lamc*P
  A.sub <- ddbb + lam1D
  B.sub <- kronecker(crossprod(D.est, Z_tilde),bb)
  C.sub <- kronecker(crossprod(Z_tilde, D.est),bb)
  zzbb <- kronecker(crossprod(Z_tilde),bb)
  lam2D <- lamc * kronecker(diag( c(1,rep(0,a)) ),P)
  D.sub <- zzbb + lam2D
  J <- B.sub %*% solve(D.sub)
  a_bdc <- solve(A.sub - J%*%C.sub)
  
  X.true  = cbind(rep(1,n),D,Z_sub)
  R <- m$Y - tcrossprod(X.true, m$beta)
  G_e <- crossprod(R)/n
  
  bgb <- crossprod(BQ2, G_e%*%BQ2)
  m.ind <- kronecker(crossprod(D.est), bgb) - 
    2*J %*% (kronecker(crossprod(Z_tilde, D.est),bgb)) +
    J %*% (kronecker(crossprod(Z_tilde),bgb)) %*%t(J)
  var.theta <- a_bdc %*% m.ind %*%a_bdc
  var.beta <- c()
  for (i in 1:n_pos) {
    v <- BQ2[i,]%*% tcrossprod(var.theta, t(BQ2[i,]))
    var.beta <- c(var.beta, v)
  }
  sd.beta <- sqrt(var.beta)
  return(sd.beta)
}

#### compute sigma_hat ####
theta_norm <- apply(m.fit$theta_tilde[, -c(1:2)], 2, 
                    function(x) sum(x)) 
Z_orac <- Z_matrix[, theta_norm != 0]
a <- dim(Z_orac)[2]
# Attention: lamc may be wrong
#lamc <- 10/(n_pos*n) 
#lamc <- 10

### estimate D.hat
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values

Z_tilde <- cbind(rep(1,n), Z_orac)
bb <- crossprod(BQ2)
ddbb <- kronecker(crossprod(D.est),bb)
lam1D <- lamc*P
A.sub <- ddbb + lam1D
B.sub <- kronecker(crossprod(D.est, Z_tilde),bb)
C.sub <- kronecker(crossprod(Z_tilde, D.est),bb)
zzbb <- kronecker(crossprod(Z_tilde),bb)
lam2D <- lamc * kronecker(diag( c(1,rep(0,a)) ),P)
D.sub <- zzbb + lam2D
J <- B.sub %*% solve(D.sub)
a_bdc <- solve(A.sub - J%*%C.sub)

theta_tilde <- m.fit$theta_tilde
theta_alpha_all <- theta_tilde[, -c(1:2)]
theta_alpha <- theta_alpha_all[, theta_norm != 0]
D.true_const <- cbind(rep(1,n),D)
R <- Y - tcrossprod(D.true_const, (BQ2 %*% as.matrix(theta_tilde[, c(1:2)]))) - 
  tcrossprod(Z_orac, (BQ2 %*%  as.matrix(theta_alpha)))
G_e <- crossprod(R)/n

bgb <- crossprod(BQ2, G_e%*%BQ2)
m.ind <- kronecker(crossprod(D.est), bgb) - 
  2*J %*% (kronecker(crossprod(Z_tilde, D.est),bgb)) +
  J %*% (kronecker(crossprod(Z_tilde),bgb)) %*%t(J)
var.theta <- a_bdc %*% m.ind %*%a_bdc
var.beta <- c()
for (i in 1:n_pos) {
  v <- BQ2[i,]%*% tcrossprod(var.theta, t(BQ2[i,]))
  var.beta <- c(var.beta, v)
}
sd.beta <- sqrt(var.beta)

#### no bootstrap
Za <- matrix(qnorm((1-0.05/2)),nrow=1)
Za <- matrix(rep(Za,each=n_pos),nrow=n_pos)
beta_alpha <- BQ2 %*% as.matrix(theta_prop)
cc.l.nb <- beta_alpha[,2]-sd.beta*Za
cc.u.nb <- beta_alpha[,2]+sd.beta*Za

test_raw <- between(beta.true[,2], cc.l.nb, cc.u.nb)
cov_raw <- sum(test_raw)/length(test_raw) ### 0.8958
cov_raw

### Data generation process ####
X.true  = cbind(rep(1,n),D,Z_matrix)
R.hat <- Y - tcrossprod(X.true,(BQ2 %*% as.matrix(theta_prop)))# 200*7505
VV <- crossprod(BQ2)
VV.inv <- chol2inv(chol(VV))
P.band <- BQ2%*%VV.inv%*%t(BQ2)
eta.hat <- t(tcrossprod(P.band,R.hat)) # 50*7505
error.hat <- R.hat-eta.hat 
xi.hat <- D - D.est ## 200
#nboot <- 200
nboot <- 300

#### alpha=0.1 ####
alpha0=0.1

a.grid <- c(10^(seq(-10, -4, by=2)), (2:9)*1e-4, seq(0.001,alpha0,0.004))
Zalpha <- matrix(qnorm((1-a.grid/2)),nrow=1)
inside_CI <- rep(0, times = length(a.grid))
beta_b_hat.mat <- array(data = NA, dim = c(nboot, n_pos))
sigma_b_hat.mat <- array(data = NA, dim = c(nboot, n_pos))
beta_alpha <- BQ2 %*% as.matrix(theta_prop)
set.seed(NULL)

#### revise begin 
ind_CI_bool <- array(data = NA, dim = c(nboot, n_pos, length(a.grid)))
#### revise end

for (i in 1:nboot) {
  nu <- sample(c(-1,1),n,replace=TRUE)
  nu_mat <- matrix(sample(c(-1,1),n*n_pos,replace=TRUE), nrow=n, ncol=n_pos)
  D_star <- D.est + nu*xi.hat
  X <- cbind(rep(1,n), D_star, Z_matrix)
  Y_star <- tcrossprod(X,beta_alpha) + nu*eta.hat + nu_mat*error.hat
  
  d.fit_b <- lm(D_star~Z_matrix)
  D.est_b <- d.fit_b$fitted.values
  X_sub <- cbind(rep(1,n), D.est_b, Z_sub)
  
  fit.b <- fit.FDAimage(Y_star,X_sub,loc,V,Tr,d,r,10)
  beta_b_hat <- fit.b$beta[,2]
  mu_b_hat <- fit.b$beta[,1]
  beta_b_hat.mat[i,] <- beta_b_hat
  print(i)
  print(mean(abs(beta_hat-beta_b_hat)))
  sigma_b_hat <- comp_sigmahat(D_star, Z_matrix, Z_orac=Z_sub, m=fit.b)
  sigma_b_hat.mat[i,] <- sigma_b_hat
  print(mean(sigma_b_hat))
  ### step 4: boot SCC
  CI_array <- sapply(Zalpha, function(z) {
    cbind(beta_b_hat-z*sigma_b_hat, beta_b_hat+z*sigma_b_hat)
  }, simplify = "array")
  ### step 5: prep
  dim(CI_array) <- c(n_pos, 2, length(Zalpha))
  ind_CI_bool[i,,] <- apply(CI_array, 3, function(ci) {
    between(beta_hat, ci[,1], ci[,2])
  })
}
### step 5
tau_hat <- apply(ind_CI_bool, c(2,3), function(x){
  sum(x)/nboot
})

### step 6
tmp1 <- tau_hat-1+alpha0
tmp1[tmp1<0] <- 1e6
alpha_adjust_s <- apply(tmp1, 1, function(x){
  a.grid[max(which(x==min(x)))]
})

alpha.adjust <- min(alpha_adjust_s)
alpha.adjust ## 0.009
Za <- matrix(qnorm((1-alpha.adjust/2)),nrow=1)
Za <- matrix(rep(Za,each=n_pos),nrow=n_pos)

cc.l <- beta_alpha[,2]-sd.beta*Za
cc.u <- beta_alpha[,2]+sd.beta*Za

test <- between(beta.true[,2], cc.l, cc.u)
test.ind <- ifelse(test, 1, 0)
cov_boot <- sum(test)/length(test)  # 0.9834777
cov_boot
