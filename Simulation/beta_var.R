### choose the lambda_2 selected by BIC 
m.fit <- PGD_stop(Y=Y, Z_matrix=Z_matrix, D.est=D.est, 
                  theta_init=theta_init,  
                  lambda_1=10/(n_pos*n), lambda_2=0.01, 
                  a=5, num_iterations=500, 
                  stop ="loss.min")

### select invalid effects
theta_norm <- apply(m.fit$theta_tilde[, -c(1:2)], 2, 
                    function(x) norm(x, type = "F")^2) 
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
