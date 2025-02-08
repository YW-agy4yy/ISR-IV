library(fields)
seed <- 1
dir_out <- "/home/agy4yy/ISR_IV/output/"
load(paste0(dir_out,"test_",seed,".RData"))

### choose the lambda_2 selected by BIC 
m.fit <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, 
                  theta_init=theta_init,  
                  lambda_1=10/(n_pos*n), lambda_2=lambda_2, 
                  a=5, num_iterations=500, 
                  stop ="loss.min")
theta_prop <- m.fit$theta_tilde
beta_hat <- BQ2 %*% theta_prop[,2]
theta_norm <- apply(m.fit$theta_tilde[, -c(1:2)], 2, 
                    function(x) sum(x)) 
Z_sub <- Z_matrix[, theta_norm != 0]

#### Calculate sigma_hat(s_j) ####
comp_sigmahat <- function(D, Z_matrix, Z_orac, mfit){
  a <- dim(Z_orac)[2]
  lamc <- 10/(n_pos*n) 
  
  ### estimate D.hat
  d.fit <- lm(D~Z_matrix)
  D.est <- d.fit$fitted.values
  
  Z_tilde <- cbind(rep(1,200), Z_orac)
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
  return(sd.beta)
}


### Data generation process ####
R.hat <- m.fit$residual # 200*7505
xi.hat <- D - D.est ## 200
nboot <- 100

alpha0 <- 0.05
a.grid <- c(10^(seq(-10, -4, by=2)), (2:9)*1e-4, seq(0.001,alpha0,0.004))
Zalpha <- matrix(qnorm((1-a.grid/2)),nrow=1)

set.seed(NULL)
inside_CI <- rep(0, times = 25)
for (i in 1:nboot) {
  nu <- sample(c(-1,1),n,replace=TRUE)
  D_star <- D.est + nu*xi.hat
  X <- cbind(rep(1,n), D_star, Z_matrix)
  Y_star <- tcrossprod(X,(BQ2 %*% as.matrix(theta_prop))) + nu*R.hat
  X_sub <- cbind(rep(1,n), D_star, Z_sub)
  fit.b <- fit.FDAimage(Y_star,X_sub,loc,V,Tr,d,r,lamc)
  
  beta_b_hat <- fit.b$beta[,2]
  print(i)
  print(mean(beta_hat-beta_b_hat))
  sigma_b_hat <- comp_sigmahat(D, Z_matrix, Z_orac=Z_sub, mfit=fit.b)
  print(mean(sigma_b_hat))
  CI_array <- sapply(Zalpha, function(z) {
    cbind(beta_b_hat-z*sigma_b_hat, beta_b_hat+z*sigma_b_hat)
  }, simplify = "array")
  dim(CI_array) <- c(7505, 2, length(Zalpha))
  ins_CI_bool <- apply(CI_array, 3, function(ci) {
    all(beta_hat >= ci[,1] & beta_hat <= ci[,2])
  })
  inside_CI <- inside_CI + as.integer(ins_CI_bool)
}

tau_hat <- inside_CI/nboot
tmp1 <- tau_hat-1+alpha0
tmp1[tmp1<0] <- 1e6
alpha.adjust <- a.grid[max(which.min(tmp1))]

