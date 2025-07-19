# Region-wise IV Regression (Estimation)

library(FDAimage) 
library(MASS)
library(matrixcalc)
library(ivreg)
library(sisVIVE)
library(parallel)
library(fields)

#### Data Generation ####
### Setting 1: 
#' (1) endogeneity: sigma_star = 0.1
#' (2) number of invalid IVs: s=2
#' (3) Z_i, Z_j are independent
#' (4) overall strength: 126
#' (5) relative strength: Equal
#' (6) sample size n = 50
#' (7) lambda1 = 0.1; lambda2 = 0.02

#args <- commandArgs(trailingOnly = TRUE)

# settings
n <- 50
card_alpha <- 10
err_sd <- 0.5
seed <- 12345
set.seed(seed)


print(n)
print(card_alpha)
print(err_sd)
print(seed)


sigma_star <- 0.3
lambda1 <- 0.05; lambda2 <- 0.01
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
#beta0 <- -5*xx^2 + 7*yy^2 #-3*xx^3 + 5*yy^3
beta0 <- -3*xx^3 + 5*yy^3
### revise here from 8* to 4*
beta1 <- 3*(xx+0.2)^2-(yy+0.2)^2  # 7505*1
beta.true <- cbind(beta0, beta1)
alpha_mat <- matrix(0, ncol = L, nrow = n_pos)  # 7505*10
alpha_1 <- 2.5 * ((xx - 0.8)^2 + (yy - 0.9)^2)+1.9  
alpha_2 <- 3 * ((xx - 1)^2 + (yy - 1)^2)+2
alpha_3 <- 1.5 * ((xx - 2)^2 + 2 * (yy - 0.2)^2)
alpha_4 <-  ((xx - 1)^2 + 2.5 * (yy - 1.2)^2)+0.8
alpha_5 <- 3 * ((xx - 1)^2 + (yy - 0.5)^2)+0.5  
alpha_6 <- 4 * ((xx - 1.5)^2 + (yy - 0.5)^2)  
alpha_7 <- 2 * ((xx + 0.5)^2 + (yy+0.3)^2) 
alpha_8 <- 1.2 * ((xx - 1)^2 + (yy - 0.5)^2)
alpha_9 <- 2 * ((xx - 1.5)^2 + (yy - 1.3)^2)  
alpha_10 <- 1.8 * ((xx + 0.5)^2 + (yy+0.5)^2) 


###### Choose sparsity square+circle ######
# (1) triangulation BQ2
bb = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
VT = TriMesh(bb, 5)
V <- VT$V
Tr <- VT$Tr
d=2; r=1;
#d=-1; r=1;
Ball <- basis(V, Tr, d, r, loc)
# B: 7505  108
# Q2 108  15
# K: 108 108
# BQ2: 7505 15
# H: 105 108
Q2 <- Ball$Q2 # Q2
B <- Ball$B   # basis function: dim- n by nT*{(d+1)(d+2)/2}
K <- Ball$K    # thin-plate energy function (Penalty)
#BQ2 <- as.matrix(B%*%Q2)
H <- as.matrix(Ball$H)

#### square+circle
square_size_1 <- 0.2
square_size_2 <- 0.3
square_size_3 <- 0.16
center1 <- c(0.22, 0.8)
center2 <- c(0.8, 0.26)
center3 <- c(0.15, 0.3)
xx <- loc[, 1]
yy <- loc[, 2]

inside_square <- function(x, y, cx, cy, size) {
  x >= (cx - size/2) & x <= (cx + size/2) &
    y >= (cy - size/2) & y <= (cy + size/2)
}
in_square1 <- inside_square(xx, yy, center1[1], center1[2], square_size_1)
in_square2 <- inside_square(xx, yy, center2[1], center2[2], square_size_2)
in_square3 <- inside_square(xx, yy, center3[1], center3[2], square_size_3)
center1 <- c(0.4, 0.2)
center2 <- c(0.8, 0.77)
center3 <- c(0.5, 0.5)
radius_1 <- 0.08
radius_2 <- 0.15
radius_3 <- 0.15
dist1_sq <- (xx-center1[1])^2+(yy-center1[2])^2
dist2_sq <- (xx-center2[1])^2+(yy-center2[2])^2
dist3_sq <- 4*(xx-center3[1])^2+ (yy-center3[2])^2
in_circle1 <- dist1_sq<= radius_1^2
in_circle2 <- dist2_sq <= radius_2^2
in_circle3 <- dist3_sq <= radius_3^2
kr1 <- in_square1
kr2 <- in_square1|in_square2
kr3 <- in_square2|in_square3
kr4 <- in_square1|in_circle1
kr5 <- in_square2|in_circle1
kr6 <- in_square3
kr7 <- in_square3|in_circle1
kr8 <- in_circle2
kr9 <- in_circle2|in_circle3
kr10 <- in_circle3
alpha_1[!kr1] <- 0
alpha_2[!kr2] <- 0
alpha_3[!kr3] <- 0
alpha_4[!kr4] <- 0
alpha_5[!kr5] <- 0
alpha_6[!kr6] <- 0
alpha_7[!kr7] <- 0
alpha_8[!kr8] <- 0
alpha_9[!kr9] <- 0
alpha_10[!kr10] <- 0

##### alpha testing #####
rect.plot <- function(betai,zlim,est.name){
  z1 <-u1; z2 <- v1;
  n1 <- length(z1); n2 <- length(z2);
  betai.mtx <- matrix(betai,ncol=n1,nrow=n2)
  #col <- colorRampPalette(c("red", "white", "blue"),space = "Lab")
  image.plot(z2, z1, betai.mtx, zlim = zlim, legend.shrink = 0.8, cex=10,
             legend.width = 1.2, main = est.name, xlab = "", ylab = "", lwd = 1)
}


#rect.plot(group_labels, zlim = c(1,50), "group")
rect.plot(alpha_1, zlim = c(0.5,5), "signal1")


for (a in 1:card_alpha) {
  alpha_mat[,a] <- get(paste0("alpha_",a))
}

# coefficient: gamma.true (11)
# (5) relative strength = equal 3
### revise
gamma <- rep(1.5,L)
gamma0 <- 0.5
### revise end
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
### revise
Z_matrix <- mvrnorm(n = n, mu = mu_Z, Sigma = diag(0.3, L)) #50x10
Z_matrix <- 0.5*Z_matrix
### revise end
Z_mat_const <- cbind(rep(1,L),Z_matrix) #50x11

# Exposure: D (50)
D <- as.vector(Z_mat_const%*%gamma.true)+xi
# overall strength
# t(D)%*%D/2/sigma_star^2

# components in eta
lamda_zeta1 <- sqrt(lambda1)*zeta_1  # 50*1
lamda_zeta2 <- sqrt(lambda2)*zeta_2
phi1 <-1.42* sin(2 * pi * xx) # 0.56 * sin(2 * pi * xx)  # 7505*1
phi2 <- 1.40 * cos(2 * pi * yy) # 0.61 * cos(2 * pi * yy)

# individual level data
Y <- c()
c.eff <- c()
iv.eff <- c()
noise <- c()
eta.true <- c()
err <- array(data = NA, dim = c(n, n_pos))
for (i in 1:n) {
  # error_i(s) ~ N(0,1)
  # revise sd=1 to 0.5
  err.i <- rnorm(n=n_pos,mean = 0, sd=err_sd)
  err[i,] <- err.i
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



#### SNR computing ####
VAR_1 <- rep(NA, n_pos)
VAR_2 <- rep(NA, n_pos)
var_mu <- var(beta.true[,1])
for (j in 1:n_pos) {
  VAR_1[j] <- var(D*beta.true[j,2]+Z_matrix[,1:card_alpha]%*%alpha_mat[j,1:card_alpha]) 
  VAR_2[j] <- var(lamda_zeta1*phi1[j] + lamda_zeta2*phi2[j] + err[,j])
}
snr <- (sum(VAR_1)+var_mu)/(sum(VAR_2))
snr

VAR_1 <- rep(NA, n_pos)
VAR_2 <- rep(NA, n_pos)
var_mu <- var(beta.true[,1])
for (j in 1:n_pos) {
  VAR_1[j] <- var(D*beta.true[j,2]+Z_matrix[,1:card_alpha]%*%alpha_mat[j,1:card_alpha]) 
  VAR_2[j] <- var(lamda_zeta1*phi1[j] + lamda_zeta2*phi2[j] + err[,j])
}
snr <- (sum(VAR_1)+var_mu)/(sum(VAR_2))
snr

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
SubTriMesh <- function(V,T) {
  edge_key <- function(i,j) paste(sort(c(i, j)), collapse = "-")
  midpoint_cache <- list()
  new_vertices <- list()
  new_vertex_index <- nrow(V)
  get_midpoint <- function(i, j) {
    key <- edge_key(i, j)
    if (!key %in% names(midpoint_cache)) {
      new_vertex_index <<- new_vertex_index + 1
      midpoint <- (V[i,] + V[j,])/2
      midpoint_cache[[key]] <<- new_vertex_index
      new_vertices[[length(new_vertices)+1]] <<- midpoint
    }
    midpoint_cache[[key]]
  }
  
  Tnew <- matrix(0, 4 * nrow(T), 3)
  for (t in seq_len(nrow(T))) {
    v1 <- T[t,1]
    v2 <- T[t,2]
    v3 <- T[t,3]
    m12 <- get_midpoint(v1, v2)
    m23 <- get_midpoint(v2, v3)
    m31 <- get_midpoint(v3, v1)
    #### 4 sub-tri
    base <- 4*(t - 1)
    Tnew[base+1,] <- c(v1, m12, m31)
    Tnew[base+2,] <- c(v2, m23, m12)
    Tnew[base+3,] <- c(v3, m31, m23)
    Tnew[base+4,] <- c(m12, m23, m31)
  }
  
  if (length(new_vertices) > 0) {
    Vnew <- rbind(V, do.call(rbind, new_vertices))
  } else {
    Vnew <- V
  }
  return(list(V=Vnew, Tr=Tnew))
}

# (1) triangulation BQ2
bb = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
VT = TriMesh(bb, 5)
V <- VT$V
Tr <- VT$Tr
#### BPST (d=2) ####
d=2; r=1; # PCST
#d=-1; r=1;
submesh <- SubTriMesh(VT$V, VT$Tr)
Vsub <- submesh$V
Tsub <- submesh$Tr
#valid_tri <- areas >= 1e-12
#Tsub <- Tsub[valid_tri, , drop = FALSE]
Ballsub <- basis(Vsub, Tsub, d, r, loc)
TriPlot(Vsub, Tsub)
Q2sub <- Ballsub$Q2 # Q2
Bsub <- Ballsub$B   # basis function: dim - n by nT*{(d+1)(d+2)/2}
#Ksub <- Ballsub$K    # thin-plate energy function (penalty)
BQ2 <- as.matrix(B%*%Q2)
Hsub <- as.matrix(Ballsub$H)

const.ball.sub <- basis(Vsub, Tsub, -1, r, loc)
const.bi.sub <- as.matrix(const.ball.sub$Bi)
#const.bi.vec <- rowSums(const.bi)
group_labels.sub <- apply(const.bi.sub, 1, function(row) which(row==1)[1])
nTsub <- dim(Tsub)[1] #200
nBgsub <- dim(Bsub)[2]/nTsub #6


# (2) alpha.sisv : 7505*10
nB <- dim(B)[2]
nBsub <- dim(Bsub)[2]
theta_alpha_init <- data.frame(matrix(ncol = L, nrow = nBsub))
theta_beta_init <- data.frame(matrix(ncol = 2, nrow = nB))
XtX <- crossprod(B)
XtXsub <- crossprod(Bsub)
for (i in 1:10) {
  alpha_tilde <- alpha.sisv[,i]
  Xty <- crossprod(Bsub, alpha_tilde)   
  theta_alpha_init[,i] <- as.vector(solve(XtXsub, Xty))
}
Xty <- crossprod(B, beta.sisv[,1])   
theta_beta_init[,1] <- as.vector(solve(XtX, Xty))
Xty <- crossprod(B, beta.sisv[,2])   
theta_beta_init[,2] <- as.vector(solve(XtX, Xty))

###### grouping ######
const.ball <- basis(V, Tr, -1, r, loc)
const.bi <- as.matrix(const.ball$Bi)
#const.bi.vec <- rowSums(const.bi)
group_labels <- apply(const.bi, 1, function(row) which(row==1)[1])
table(group_labels)
nT <- dim(Tr)[1] #50
nBg <- dim(B)[2]/nT #6

groupcomb <- cbind(1:n_pos, group_labels, group_labels.sub)

##### Stage 1 #####
d.fit <- lm(D~Z_matrix)
D.est <- d.fit$fitted.values
#theta_init <- cbind(theta_beta_init,theta_alpha_init)
L = ncol(Z_matrix)
X = cbind(rep(1,n),D.est,Z_matrix)


#### Oracle ####
lambda_1 <- 10
## construct U1,U2
q1 <- ncol(BQ2) # 26
q2_each <- nBgsub  # 6
d <- 2*q1 + ncol(Bsub)
U1 <- kronecker(matrix(1, n, 1), BQ2)  # (nN x q1)
D_vec <- rep(X[, 2] , each = n_pos)  # (nN x 1)
U2 <- U1 * D_vec
## U3
Z <- X[, 3:(L+2)]  
Bsub_dense <- as.matrix(Bsub) 

# COnstruct B_{2,A_l}
kr_list <- mget(paste0("kr", 1:10)) 
cols_to_keep_list <- lapply(1:L, function(l) {
  kr <- kr_list[[l]]
  mathcalA <- tapply(kr, groupcomb[,3], function(x) mean(x) >= 0.1)
  keep_triangles <- which(!is.na(mathcalA) & mathcalA)
  unlist(lapply(keep_triangles, function(tri) {
    start_col <- (tri-1)*q2_each + 1
    end_col <- tri*q2_each
    start_col:end_col
  }))
})
U3_list <- lapply(1:L, function(l) {
  Zl <- Z[, l]  # IV values
  Zl_rep <- rep(Zl, each = n_pos)
  B_l <- Bsub_dense[, cols_to_keep_list[[l]], drop = FALSE]
  Bk_dense <- as.matrix(kronecker(matrix(1, n, 1), B_l))  
  Zl_rep * Bk_dense
})
U3 <- do.call(cbind, U3_list)
sapply(U3_list, ncol)
U_mat <- cbind(U1, U2, U3)
UTU <- crossprod(U_mat)  ### 1354*1354

P <- crossprod(Q2, K%*%Q2)  # q1*q1
d3 <- ncol(U3)

DLam <- bdiag(lambda_1*P, lambda_1*P,Matrix(0, d3, d3))

lhs=UTU+DLam
Y_vec <- as.vector(t(Y))
rhs <- crossprod(U_mat, Y_vec)
theta_tilde <- solve(lhs, rhs) 

theta_0 <- theta_tilde[1:q1]                        
theta_beta<- theta_tilde[(q1+1):(2*q1)] 
theta_alpha <- theta_tilde[(2*q1+1):length(theta_tilde)] 
beta.orac <- BQ2%*%as.matrix(theta_beta)
mise_orac <- mean((beta.true[,2]-beta.orac)^2)


#### Linesearch ####
lambda_2_all <- c(0.001, 0.01, 0.1, 1, 10, 100)
bic_all <- rep(NA, length(lambda_2_all))
mise_all <- matrix(NA, length(lambda_2_all),2)
###### bic-1 ######
bic_iter <- 1
a <- 0.00008
lambda_2 <- lambda_2_all[bic_iter]
lambda_1 <- 10
theta_allg_beta <- theta_beta_init
theta_allg_alpha <- theta_alpha_init
start.time <- Sys.time()
num_iterations <- 300
loss.vec <- rep(NA,num_iterations)
a.vec <- rep(NA,num_iterations)
grad.sum <- 0
step_shink <- 0.5
for (k in 1:num_iterations) {
  term1_nT <- rep(NA,nT)
  term3_nT <- rep(NA,nT)
  ###### update (
  theta_allg_beta.old <- theta_allg_beta
  theta_allg_alpha.old <- theta_allg_alpha
  a.old <- a
  ###### update )
  for (grp_ind in 1:nT) {
    loc_list <- which(group_labels==grp_ind)
    n_pos_g <- length(loc_list)
    bs_list <- ((grp_ind-1)*nBg+1):(grp_ind*nBg)
    bs_list_sub <- ((grp_ind-1)*nBgsub*4+1):(grp_ind*nBgsub*4)
    theta_g_beta <- theta_allg_beta[bs_list,]
    theta_g_alpha <- theta_allg_alpha[bs_list_sub,]
    theta_g_vec.old <- c(as.matrix(theta_g_beta), as.matrix(theta_g_alpha))
    #theta_g_vec <- unlist(as.vector(theta_g))
    B_g <- B[loc_list, bs_list]
    B_g_sub <- Bsub[loc_list, bs_list_sub]
    K_g <- K[bs_list, bs_list]
    init_l2norm_g <- sqrt(colSums(theta_alpha_init[bs_list_sub, ]^2))
    weight.vec <- 1/(init_l2norm_g + abs(rnorm(L, sd=0.0001)))
    ##### Gradient ####
    ###### (1) grad1 
    #Yhat <- tcrossprod(X,(B_g %*% as.matrix(theta_g)))
    Yhat <- tcrossprod(cbind(rep(1,n),D.est),(B_g %*% as.matrix(theta_g_beta))) +
      tcrossprod(Z_matrix,(B_g_sub %*% as.matrix(theta_g_alpha)))
    residual.matrix <- Y[,loc_list] - Yhat
    grad_list <- lapply(1:n, function(i) {
      r_i <- residual.matrix[i, ]
      Z_i <- Z_matrix[i, ]
      Dest_i <- D.est[i]
      g1 <- crossprod(B_g, r_i)            
      g2 <- Dest_i * g1                    
      g_sub <- do.call(c, lapply(1:L, function(j) {
        Z_i[j] * crossprod(B_g_sub, r_i)   
      }))
      c(g1, g2, g_sub)                     
    })
    grad_list_flat <- lapply(grad_list, function(block) {
      as.numeric(do.call(c, lapply(block, as.numeric)))  
    })
    grad1 <- Reduce("+", grad_list_flat)
    #### (2) grad2
    grad2 <- c(lambda_1*as.vector(K_g %*% theta_g_beta[,1]), lambda_1*as.vector(K_g %*% theta_g_beta[,2]), rep(0,4*nBgsub*L))
    #theta_g_vec <- theta_g_vec + 2*a*(grad1/(n*n_pos_g)-grad2)
    theta_g_vec <- c(as.matrix(theta_g_beta), as.matrix(theta_g_alpha))
    giter_grad_old <- -2*(grad1/(n*n_pos)-grad2)
    theta_g_vec <- theta_g_vec - a*giter_grad_old
    theta_g_beta <- matrix(theta_g_vec[1:(2*nBg)], nrow = nBg, ncol = 2)
    theta_g_alpha <- matrix(theta_g_vec[(2*nBg + 1):length(theta_g_vec)], nrow = (4*nBgsub), ncol = L)
    ##### Step 2: Proximal #####
    norm_grp_alpha <- rep(NA, L)
    new.norm_grp_alpha <- rep(NA, L)
    for (ell in 1:L) {
      grp_alpha <- theta_g_alpha[,ell]
      norm_grp_alpha[ell] <- sqrt(sum(grp_alpha^2))
      shrinkage_factor <- max(1 - a * lambda_2 * weight.vec[ell]/norm_grp_alpha[ell], 0)
      theta_g_alpha[,ell] <- shrinkage_factor * grp_alpha
      ###### update (
      new.norm_grp_alpha[ell] <- sqrt(sum(theta_g_alpha[,ell]^2))
      ####### update )
      
    }
    theta_allg_beta[bs_list,] <- theta_g_beta
    theta_allg_alpha[bs_list_sub,] <- theta_g_alpha
    
    term3_nT[grp_ind] <- weight.vec%*%new.norm_grp_alpha
    
    ###### update ( 
    grad.sum <- grad.sum + giter_grad_old%*% (theta_g_vec-theta_g_vec.old)
    ####### update )
  }
  
  #### update (
  Yhat <- tcrossprod(cbind(rep(1,n),D.est),(B %*% as.matrix(theta_allg_beta))) +
    tcrossprod(Z_matrix,(Bsub %*% as.matrix(theta_allg_alpha))) 
  LHS <- sum((Y-Yhat)^2)/(n*n_pos)+lambda_1*(crossprod(theta_allg_beta[,1],K)%*%theta_allg_beta[,1] 
                                               + crossprod(theta_allg_beta[,2],K)%*%theta_allg_beta[,2])
  Yhat.old <- tcrossprod(cbind(rep(1,n),D.est),(B %*% as.matrix(theta_allg_beta.old)))+tcrossprod(Z_matrix,(Bsub %*% as.matrix(theta_allg_alpha.old))) 
  RHS1 <- sum((Y-Yhat.old)^2)/(n*n_pos)+lambda_1*(crossprod(theta_allg_beta.old[,1],K)%*%theta_allg_beta.old[,1] 
                                                  +crossprod(theta_allg_beta.old[,2],K)%*%theta_allg_beta.old[,2])
  RHS2 <- grad.sum
  RHS3 <- 1/(2*a)*(sum((theta_allg_beta.old-theta_allg_beta)^2)+sum((theta_allg_alpha.old-theta_allg_alpha)^2))
  RHS <- RHS1+ RHS2+ RHS3
  #print(LHS)
  #print(RHS)
  if (as.numeric(LHS) < as.numeric(RHS)){
    a <- step_shink*a
  }
  a.vec[k] <- a
  #### update )
  
  proj_mat <- crossprod(H,chol2inv(tcrossprod(H)))%*%H
  theta_allg_beta <- as.matrix(theta_allg_beta) - proj_mat%*%as.matrix(theta_allg_beta)
  
  Yhat <- tcrossprod(cbind(rep(1,n),D.est),(B %*% as.matrix(theta_allg_beta))) +
    tcrossprod(Z_matrix,(Bsub %*% as.matrix(theta_allg_alpha)))
  bic.residual.matrix <- Y-Yhat
  term1 <- sum(bic.residual.matrix^2)/(n*n_pos)
  print(k)
  print(a)
  print(term1)
  
  term2 <- lambda_1*(crossprod(theta_allg_beta[,1],K)%*%theta_allg_beta[,1] + crossprod(theta_allg_beta[,2],K)%*%theta_allg_beta[,2])
  term3 <- lambda_2*(sum(term3_nT))
  loss.vec[k] <- term1 + term2 + term3
}
end.time <- Sys.time()
end.time - start.time

alpha.prop <- Bsub%*%as.matrix(theta_allg_alpha)
beta.prop <- B%*%as.matrix(theta_allg_beta)

eta_bic <- 0.5
A_hat <- sum(theta_allg_alpha != 0)/nBg
bic_all[bic_iter] <- log(term1) + A_hat/(n*nT)*(log(n) + 2*eta_bic*log(L*nBg*nT))
mise_all[bic_iter] <- mean((beta.true[,2]-beta.prop[,2])^2)




bic_lam <- which.min(bic_all)
mise_prop <- mise_all[which.min(bic_all)]
mise_orac <- mean((beta.true[,2]-beta.orac)^2)
mise.list <- list(bic_lam, mise_prop, mise_orac)  
