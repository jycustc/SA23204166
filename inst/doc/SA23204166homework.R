## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
my.sample <- function(x, size, prob = NULL) {
  n=length(x)
  if (is.null(prob)) {
    p=rep(1,n)/n
  } else {
    p=prob
  }
  cp <- cumsum(p)
  U = runif(size) 
  r <- x[findInterval(U,cp)+1]
  return(r)
}

## -----------------------------------------------------------------------------
my.sample(1:3,10,prob = c(.2, .3, .5))
table(my.sample(1:3,10000,prob = c(.2, .3, .5)))

## -----------------------------------------------------------------------------
my.sample(1:3,10)
table(my.sample(1:3,10000))

## -----------------------------------------------------------------------------
my.sample(c(1,2,2,3,3,3),10)
table(my.sample(c(1,2,2,3,3,3),10000))

## -----------------------------------------------------------------------------
my.sample(c(1,2,2,3,3,3),10,prob=c(.3,.1,.05,.3,.15,.1))
table(my.sample(c(1,2,2,3,3,3),10000,prob=c(.3,.1,.05,.3,.15,.1)))

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
F.inv <- function(x) {
  return(ifelse(x<=0.5,log(2*x),-log(2-2*x)))
}
x <- sapply(u,F.inv)
hist(x, prob = TRUE, main = expression(f(x)==frac(1,2)*e^-abs(x))) 
y <- seq(min(x), max(x), .01)
lines(y, exp(-abs(y))/2)

## -----------------------------------------------------------------------------
my.beta <- function(n,a,b) {
  j<-k<-0;y <- numeric(n) 
  while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1)
    if (x^(a-1)*(1-x)^(b-1) > u) {
      k <- k + 1
      y[k] <- x 
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
y<-my.beta(1000,3,2)
hist(y, prob = TRUE, main = expression(f(x)==12*x^2*(1-x))) 
z <- seq(min(y), max(y), .01)
lines(z, 12*z^2*(1-z))

## -----------------------------------------------------------------------------
my.fe <- function(n) {
  y <- numeric(n) 
  k<-0
  while (k < n) {
    u1 <- runif(1,min=-1)
    u2 <- runif(1,min=-1)
    u3 <- runif(1,min=-1)
    if(abs(u3)>=abs(u2)&&abs(u3)>=abs(u1)){y[k] <- u2}else{y[k] <- u3}
    k <- k+1
  }
  return(y)
}

## -----------------------------------------------------------------------------
y <- my.fe(1000)
hist(y, prob = TRUE, main = expression(f(x)==frac(3,4)*(1-x^2))) 
z <- seq(min(y), max(y), .01)
lines(z, 3*(1-z^2)/4)

## -----------------------------------------------------------------------------
var_MC <- function(rho,n = 10^6,K = 100) {
  set.seed(0)
  d <- 1
  l <- d*rho
  pihat <- numeric(K)
  for (i in 1:K) {
    X <- runif(n,0,d/2)
    Y <- runif(n,0,pi/2)
    pihat[i] <- 2*l/d/mean(l/2*sin(Y)>X)
  }
  return(var(pihat))
}

## -----------------------------------------------------------------------------
var_MC(0.1)

## -----------------------------------------------------------------------------
var_MC(0.5)

## -----------------------------------------------------------------------------
var_MC(1)

## -----------------------------------------------------------------------------
var_MC <- function(n = 10^6,K=100) {
  set.seed(12345)
  theta_hat <- numeric(K)
  for (i in 1:K) {
    x <- runif(n)
    theta_hat[i] <- mean(exp(x))
  }
  return(var(theta_hat))
}
var_AVA <- function(n = 10^6,K=100) {
  set.seed(12345)
  theta_hat <- numeric(K)
  for (i in 1:K) {
    x <- runif(n)
    theta_hat[i] <- mean(exp(x)/2+exp(1-x)/2)
  }
  return(var(theta_hat))
}
var_AVA()/var_MC()

## -----------------------------------------------------------------------------
MC <- function(n) {
  set.seed(12345)
  u <- rexp(n)+1
  f <- exp(-u^2/2)*u^2/sqrt(2*pi)
  g <- exp(1-u)
  return(mean(f/g))
}
MC(10000)

## -----------------------------------------------------------------------------
g <- function(x) x^2 * exp(-x^2/2)/sqrt(2 * pi)
integrate(g, lower = 1, upper = Inf)

## -----------------------------------------------------------------------------
set.seed(0)
M <- 10000 #number of replicates
k <- 5 #number of strata
T2 <- numeric(k)
g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1) 
}
sam <- function(M,j,k) {
  u <- runif(M/k,(j-1)/k,j/k) #f3, inverse transform method 
  x <- - log(1 - u * (1 - exp(-1)))
  return(x)
}
for (j in 1:k){
  u<-sam(M/k,j,k)
  fg <- g(u) / (exp(-u) / (1 - exp(-1))) 
  T2[j] <- mean(fg)
}
estimates <- mean(T2)
estimates
var(T2)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 20
rn <- sqrt(n)
t0 <- qt(c(0.025, 0.975), df = n - 1)
CI <- replicate(10000, expr = {
  x <- rchisq(n, df = 2)
  ci <- mean(x) + t0 * sd(x)/rn
})
LCL <- CI[1, ]
UCL <- CI[2, ]
mean(LCL < 2 & UCL > 2)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 20
t0 <- qchisq(c(0.025, 0.975), df = n - 1)
CI <- replicate(10000, expr = {
  x <- rchisq(n, df = 2)
  ci <- (n-1) * var(x) / t0
})
UCL <- CI[1, ]
LCL <- CI[2, ]
mean(LCL < 4 & UCL > 4)

## -----------------------------------------------------------------------------
set.seed(12345)  
alpha <- 0.05  


simulate_t_test <- function(population, sample_size, n_simulations) {
  results <- numeric(n_simulations)
  for (i in 1:n_simulations) {
    sample_data <- switch(population,
      "chi2" = rchisq(sample_size, df = 1),
      "uniform" = runif(sample_size, min = 0, max = 2),
      "exponential" = rexp(sample_size, rate = 1)
    )
    t_test <- t.test(sample_data, mu = 1)
    results[i] <- t_test$p.value <= alpha
  }
  return(mean(results))
}


sample_size <- 30  # 样本大小
n_simulations <- 1000  # 模拟次数

# 模拟并计算Type I错误率
results_chi2 <- simulate_t_test("chi2", sample_size, n_simulations)
results_uniform <- simulate_t_test("uniform", sample_size, n_simulations)
results_exponential <- simulate_t_test("exponential", sample_size, n_simulations)

# 输出结果
cat("Type I Error Rate for Chi-squared(1):", results_chi2*100, "%\n")
cat("Type I Error Rate for Uniform(0,2):", results_uniform*100, "%\n")
cat("Type I Error Rate for Exponential(1):", results_exponential*100, "%\n")


## -----------------------------------------------------------------------------
set.seed(12345) 
m <- 1000
alpha <- 0.1
n_h0 <- 0.95 * m 
n_h1 <- 0.05 * m  

fwer_bonf<-numeric(1000)
fwer_bh<-numeric(1000)
tpr_bonf<-numeric(1000)
tpr_bh<-numeric(1000)
fdr_bonf<-numeric(1000)
fdr_bh<-numeric(1000)

for(i in c(1:1000)){
  p_values_h0 <- runif(n_h0)
  p_values_h1 <- rbeta(n_h1, 0.01, 1)
  p_values <- c(p_values_h0, p_values_h1)
  adjusted_p_bonf <- p.adjust(p_values, method = "bonferroni")
  adjusted_p_bh <- p.adjust(p_values, method = "BH")
  reject_bonf <- adjusted_p_bonf <= alpha
  reject_bh <- adjusted_p_bh <= alpha
  fwer_bonf[i] <- sum(reject_bonf) / m
  fwer_bh[i] <- sum(reject_bh) / m
  tpr_bonf[i] <- sum(reject_bonf[1:n_h1]) / n_h1
  tpr_bh[i] <- sum(reject_bh[1:n_h1]) / n_h1
  fdr_bonf[i] <- sum(reject_bonf[1:n_h1]) / sum(reject_bonf)
  fdr_bh[i] <- sum(reject_bh[1:n_h1]) / sum(reject_bh)
}

cat("Bonferroni校正后的FWER:", mean(fwer_bonf), "\n")
cat("B-H较正后的FWER:", mean(fwer_bh), "\n")
cat("Bonferroni校正后的TPR:", mean(tpr_bonf), "\n")
cat("B-H较正后的TPR:", mean(tpr_bh), "\n")
cat("Bonferroni校正后的FDR:", mean(fdr_bonf), "\n")
cat("B-H较正后的FDR:", mean(fdr_bh), "\n")


## -----------------------------------------------------------------------------
Bootstrap_Bias <- function(n,true_lambda=2,B=1000,m=1000) {
  res<-numeric(m)
  for(i in c(1:m)){
    x <- rexp(n, rate = true_lambda)
    obj <- boot(data=x,statistic=lam,R=B)
    res[i]<-mean(obj$t)-true_lambda
    }
  return(mean(res))
}

Bootstrap_SE <- function(n,true_lambda=2,B=1000,m=1000) {
  res<-numeric(m)
  for(i in c(1:m)){
    x <- rexp(n, rate = true_lambda)
    obj <- boot(data=x,statistic=lam,R=B)
    res[i]<-sd(obj$t)
    }
  return(mean(res))
}

library(boot)
set.seed(12345) 
lam <- function(x,i){1 / mean(x[i])}
true_lambda <- 2
sample_sizes <- c(5, 10, 20)
B <- 1000
m <- 1000
results <- data.frame(Sample_Size = numeric(0), 
                      Bootstrap_Bias = numeric(0), 
                      Bootstrap_SE = numeric(0))

for (n in sample_sizes) {
    results <- rbind(results, data.frame(Sample_Size = n, 
                                        Bootstrap_Bias = Bootstrap_Bias(n), 
                                        Bootstrap_SE = Bootstrap_SE(n)))
}
theoretical_bias <- true_lambda / (sample_sizes - 1)
theoretical_se <- true_lambda * sample_sizes / ((sample_sizes - 1) * sqrt(sample_sizes - 2))

print("Simulation Results:")
print(results)
print("Theoretical Results:")
theoretical_results <- data.frame(Sample_Size = sample_sizes, 
                                  Theoretical_Bias = theoretical_bias, 
                                  Theoretical_SE = theoretical_se)
print(theoretical_results)

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
library(bootstrap)
B <- 1000
n <- nrow(law)
b.cor <- function(x,i) cor(x[i,1],x[i,2])
bt <- boot(data=law,statistic=b.cor,R=B)
alpha <- 0.05
cv <- qt(1 - alpha / 2, df = B - 1)
bt_c<-bt$t0
bt_sd<-sd(bt$t)
lb <- bt_c - cv * (bt_sd / sqrt(B))
ub <- bt_c + cv * (bt_sd / sqrt(B))

cat("Bootstrap t Confidence Interval: [", lb, ", ", ub, "]\n")

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
X<-aircondit[1]
m<-function(x,i){return(mean(as.matrix(x[i, ])))}
b_x<-boot(X, statistic = m, R = 2000) 
b_x
boot.ci(b_x, type = c("norm", "perc", "basic", "bca"))
hist(b_x$t, prob = TRUE) 

## -----------------------------------------------------------------------------
set.seed(12345)
library(bootstrap)
X<-as.matrix(scor)
n<-nrow(X)
t_j<-numeric(n)
lam <- eigen(cov(X))$values
t_h <- max(lam/sum(lam))
for(i in 1:n) {
  Y <- X[-i, ]
  m <- cov(Y)
  lam <- eigen(m)$values
  t_j[i] <- max(lam/sum(lam))
}
b_j <- (n - 1) * (mean(t_j) - t_h)
s_j <- sqrt((n - 1)/n * sum((t_j - mean(t_j))^2)) 
list(estimate = t_h, bias = b_j, standard_error = s_j)

## -----------------------------------------------------------------------------
library(DAAG) 
attach(ironslag)
n <- length(magnetic)
N <- choose(n, 2)
e1 <- e2 <- e3 <- e4 <- e5 <- numeric(N) 
ij <- 1
for(i in 1:(n - 1)){
  for (j in (i + 1):n){ 
    k<-c(i,j)
    y <- magnetic[-k]
    x <- chemical[-k]
    J1<-lm(y~x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
    e1[ij] <- sum((magnetic[k] - yhat1)^2) 
    J2<-lm(y~x+I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
    e2[ij] <- sum((magnetic[k] - yhat2)^2)
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    e3[ij] <- sum((magnetic[k] - yhat3)^2)
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4 <- exp(logyhat4)
    e4[ij] <- sum((magnetic[k] - yhat4)^2)
    c2 <- x^2
    c3 <- x^3
    J5<-lm(y~x+c2+c3)
    yhat5 <- J5$coef[1] + J5$coef[2] * chemical[k] +J5$coef[3] * chemical[k]^2 +J5$coef[4] * chemical[k]^3
    e5[ij] <- sum((magnetic[k] - yhat5)^2)
    ij<-ij+1
  }
}
c(sum(e1), sum(e2), sum(e3), sum(e4), sum(e5))/N

## -----------------------------------------------------------------------------
set.seed(12345)
cm_test <- function(x, y, R = 999) {
  n <- length(x)
  m <- length(y)
  z<-c(x,y)
  N<-n+m
  Fn <- numeric(N)
  Gm <- numeric(N)
  for(i in 1:N){
   Fn[i] <- mean(as.integer(z[i] <= x))
   Gm[i] <- mean(as.integer(z[i] <= y))
  }
  cvm0 <- ((n * m)/N) * sum((Fn - Gm)^2)
  cvm <- replicate(R, expr = { 
    k <- sample(1:N)
    Z <- z[k]
    X <- Z[1:n]
    Y <- Z[(n+1):N]
    for(i in 1:N){
      Fn[i] <- mean(as.integer(Z[i] <= X)) 
      Gm[i] <- mean(as.integer(Z[i] <= Y))
    }
    ((n*m)/N) * sum((Fn - Gm)^2)
  })
  cvm1<-c(cvm,cvm0)
  return(list(statistic = cvm0, p.value = mean(cvm1 >=cvm0)))
}
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"])) 
y <- sort(as.vector(weight[feed == "linseed"])) 
cm_test(x, y)
detach(chickwts)

## -----------------------------------------------------------------------------
maxoutliers <- function(x, y) {
  X<-x-mean(x)
  Y<-y-mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X)) 
  return(max(c(outx, outy)))
}
maxout <- function(x, y, R = 199) { 
  z<-c(x,y)
  n <- length(x)
  N <- length(z)
  stats <- replicate(R, expr = {
    k <- sample(1:N)
    k1 <- k[1:n]
    k2 <- k[(n + 1):N]
    maxoutliers(z[k1], z[k2])
    })
  stat <- maxoutliers(x, y)
  stats1 <- c(stats, stat)
  tab <- table(stats1)/(R + 1)
  return(list(estimate=stat,p=mean(stats1>=stat),freq=tab,cdf=cumsum(tab)))
}

## -----------------------------------------------------------------------------
set.seed(12345)
n1 <- 20
n2 <- 40
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1) 
y <- rnorm(n2, mu2, sigma2) 
maxout(x, y)

## -----------------------------------------------------------------------------
set.seed(12345)
sigma1 <- 1
sigma2 <- 2
x <- rnorm(n1, mu1, sigma1) 
y <- rnorm(n2, mu2, sigma2) 
maxout(x, y)

## -----------------------------------------------------------------------------
al <- function(N,b1,b2,b3,f0) {
  set.seed(1)
  x1 <- rpois(N,lambda=1)
  x2 <- rexp(N,rate = 1)
  x3 <- rbinom(N,size=1,prob=0.5)
  g <- function(alpha){
    t <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+t)
    mean(p) - f0 
    }
  solution <- uniroot(g,c(-20,10))
  round(unlist(solution),5)[1]
}

## -----------------------------------------------------------------------------
al(N=10^6,b1=0,b2=1,b3=-1,f0=0.1)
al(N=10^6,b1=0,b2=1,b3=-1,f0=0.01)
al(N=10^6,b1=0,b2=1,b3=-1,f0=0.001)
al(N=10^6,b1=0,b2=1,b3=-1,f0=0.0001)

## -----------------------------------------------------------------------------
lf<-c(1:10)
res<-numeric(10)
for(i in 1:10){
  res[i]<-al(N=10^6,b1=0,b2=1,b3=-1,f0=exp(-lf[i]))
}
plot(lf,res,xlab=expression(-logf[0]),ylab="a")

## -----------------------------------------------------------------------------
rw<-function(N, x0, sigma) {
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    x_t<-x[i-1]
    y<-rnorm(1, x_t,sigma)
    if(u[i]<=exp(abs(x_t)-abs(y))){
      x[i]<-y
    }else {
      x[i]<-x[i-1]
      k<-k+1
    }
  }
  return(list(x=x,ac=1-k/N))
}

## -----------------------------------------------------------------------------
set.seed(1)
rw1<-rw(N=3000,1,0.2)
rw2<-rw(N=3000,1,1)
rw3<-rw(N=3000,1,5)
rw4<-rw(N=3000,1,25)
cat("sigma=0.2,1,5,25接受率分别为:", (c(rw1$ac, rw2$ac, rw3$ac, rw4$ac)),"\n")

## -----------------------------------------------------------------------------
#par(mfrow = c(2, 2))
plot(rw1$x, type = "l")
plot(rw2$x, type = "l")
plot(rw3$x, type = "l")
plot(rw4$x, type = "l")
#par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
#par(mfrow = c(2, 2))
p<-ppoints(200)
y<-qexp(p, 1)
z<-c(-rev(y), y)
ft<-0.5*exp(-abs(z))
hist(rw1$x, breaks = "Scott", freq = FALSE)
lines(z, ft)
hist(rw2$x, breaks = "Scott", freq = FALSE)
lines(z, ft)
hist(rw3$x, breaks = "Scott", freq = FALSE)
lines(z, ft)
hist(rw4$x, breaks = "Scott", freq = FALSE)
lines(z, ft)
#par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
#par(mfrow = c(2, 2))
q1 <- quantile(rw1$x, p)
qqplot(z, q1)
q2 <- quantile(rw2$x, p)
qqplot(z, q2)
q3 <- quantile(rw3$x, p)
qqplot(z, q3)
q4 <- quantile(rw4$x, p)
qqplot(z, q4)
#par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
N<-3000
burn<-500
X<-matrix(0, N, 2)
r<-0.9
s1<-sqrt(1-r^2)
s2<-sqrt(1-r^2)
X[1,]<-c(0, 0)
set.seed(1)
for(i in 2:N){
  x2<-X[i-1,2]
  m1<-r*x2
  X[i,1]<-rnorm(1,m1,s1)
  x1<-X[i, 1]
  m2<-r*x1
  X[i,2]<-rnorm(1,m2,s2)
}
res<-X[c((burn+1):N),]
Xt<-res[,1]
Yt<-res[,2]
plot(Xt,Yt)

## -----------------------------------------------------------------------------
mean(Xt)
mean(Yt)
var(Xt)
var(Yt)
cov(Xt,Yt)

## -----------------------------------------------------------------------------
re<-lm(Yt~Xt)
summary(re)

## -----------------------------------------------------------------------------
plot(re$fit, re$res)
abline(h = 0)
qqnorm(re$res)
qqline(re$res)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j]) 
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi)
B <- n * var(psi.means)
psi.w <- apply(psi, 1, "var") 
W <- mean(psi.w)
v.hat <- W*(n-1)/n + (B/n) 
r.hat <- v.hat / W 
return(r.hat)
}
set.seed(1)
GR<-function(s,m,x0) {
  x<-numeric(m)
  x[1]<-x0
  u<-runif(m)
  for(i in 2:m){
    x_t<-x[i-1]
    y<-rchisq(1,df=x_t)
    num<-(y/s^2)*exp(-y^2/(2*s^2))*dchisq(x_t,df=y) 
    den<-(x_t/s^2)*exp(-x_t^2/(2*s^2))*dchisq(y,df=x_t) 
    if (u[i]<= num/den){x[i]<-y}else{x[i]<-x_t}
  }
  return(x)
}

## -----------------------------------------------------------------------------
x0<-c(1/4,1,4,16)
X<-matrix(0, nrow = 4, ncol = 3000)
for (i in 1:4) X[i,]<-GR(4,3000,x0[i])
psi<-t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))psi[i,]<-psi[i,]/(1:ncol(psi))
rhat <- Gelman.Rubin(psi)
rhat

## -----------------------------------------------------------------------------
library(coda)
x1<-as.mcmc(X[1,])
x2<-as.mcmc(X[2,])
x3<-as.mcmc(X[3,])
x4<-as.mcmc(X[4,])
Y<-mcmc.list(x1, x2, x3, x4)
print(gelman.diag(Y))
gelman.plot(Y, col = c(1, 1))

## -----------------------------------------------------------------------------
log_likelihood_direct <- function(lambda, u, v) {
  sum(log(exp(-lambda * u) - exp(-lambda * v)))
}
EM_algorithm <- function(u, v, max_iter = 1000, tol = 1e-6) {
  n <- length(u)
  lambda <- 1  
  log_likelihood_prev <- -Inf
  for (iter in 1:max_iter) {
    Q_values <- log(exp(-lambda * u) - exp(-lambda * v))
    lambda_new <- optimize(
      function(l) sum(Q_values),
      interval = c(0, 10),  
      maximum = TRUE
    )$maximum

    log_likelihood <- sum(Q_values)
    if (abs(log_likelihood - log_likelihood_prev) < tol) {
      break
    }
    
    lambda <- lambda_new
    log_likelihood_prev <- log_likelihood
  }
  return(lambda)
}
u_values <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v_values <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
lambda_direct <- optimize(
  function(l) -log_likelihood_direct(l, u_values, v_values),
  interval = c(0, 10), maximum = TRUE)$maximum
lambda_EM <- EM_algorithm(u_values, v_values)
cat("直接最大化观测数据的似然函数得到的lambda的MLE:", lambda_direct, "\n")
cat("EM算法得到的lambda的MLE:", lambda_EM, "\n")

## -----------------------------------------------------------------------------
s_g <- function(X) {
miX <- min(X)
X <- X - miX
maX <- max(X)
X <- X / max(X)
m <- nrow(X)
n <- ncol(X)
iter <- n^3
x <- c(rep(0, m), 1)
X1 <- -cbind(t(X), rep(-1, n)) 
Y1 <- rep(0, n)
X3 <- t(as.matrix(c(rep(1, m), 0))) 
Y3 <- 1
sx <- simplex(a=x, A1=X1, b1=Y1, A3=X3, b3=Y3,maxi=TRUE, n.iter=iter)
x <- c(rep(0, n), 1) 
X1 <- cbind(X, rep(-1, m)) 
Y1 <- rep(0, m)
X3 <- t(as.matrix(c(rep(1, n), 0))) 
Y3 <- 1
sy <- simplex(a=x, A1=X1, b1=Y1, A3=X3, b3=Y3,maxi=FALSE, n.iter=iter)
res <- list("A" = X * maX + miX, "x" = sx$soln[1:m],
             "y" = sy$soln[1:n],
             "v" = sx$soln[m+1] * maX + miX)
return(res)
}

## -----------------------------------------------------------------------------
A <- matrix(c(0,-2,-2,3,0,0,4,0,0, 
              2,0,0,0,-3,-3,4,0,0, 
              2,0,0,3,0,0,0,-4,-4, 
              -3,0,-3,0,4,0,0,5,0, 
              0,3,0,-4,0,-4,0,5,0, 
              0,3,0,0,4,0,-5,0,-5, 
              -4,-4,0,0,0,5,0,0,6, 
              0,0,4,-5,-5,0,0,0,6, 
              0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) 
B <- A + 2
s <- s_g(B) 
s$v
round(cbind(s$x, s$y), 7)
round(s$x , 7)

## -----------------------------------------------------------------------------
round(s$x*61, 7)

## -----------------------------------------------------------------------------
dim(as.vector(c(1:3)))

## -----------------------------------------------------------------------------
df <- data.frame(
  numeric_col = c(1, 2, 3),
  character_col = c("a", "b", "c"),
  logical_col = c(TRUE, FALSE, TRUE)
)
mat <- as.matrix(df)
mat

## -----------------------------------------------------------------------------
empty_df <- data.frame()

## -----------------------------------------------------------------------------
empty_df <- data.frame(matrix(nrow = 0, ncol = 0))

## -----------------------------------------------------------------------------
scale01 <- function(x) {
         rng <- range(x, na.rm = TRUE)
         (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
data_frame <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9)
)
# 对数据框的每一列应用
scaled <- apply(data_frame, MARGIN = 2, scale01)
scaled
# 仅对数据框的数值列应用
data_frame <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  character_col = c("a", "b", "c"),
  logical_col = c(TRUE, FALSE, TRUE)
)
numeric_columns <- sapply(data_frame, is.numeric)
scaled_numeric <- data_frame
scaled_numeric[, numeric_columns] <- apply(data_frame[, numeric_columns], MARGIN = 2, scale01)
scaled_numeric

## -----------------------------------------------------------------------------
data_frame <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9)
)
numeric_std_dev <- vapply(data_frame, sd, numeric(1))
numeric_std_dev

## -----------------------------------------------------------------------------
data_frame <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  character_col = c("a", "b", "c"),
  logical_col = c(TRUE, FALSE, TRUE)
)
numeric_columns <-  vapply(data_frame, is.numeric, numeric(1))
mixed_numeric_std_dev <- vapply(data_frame[, numeric_columns], sd, numeric(1))
mixed_numeric_std_dev


## -----------------------------------------------------------------------------
library(microbenchmark)
gibbs_sampler_r <- function(n_iterations, a, b, n) {
  samples <- matrix(0, nrow = iterations, ncol = 2)
  x <- 0
  y <- 0.5
  for (i in 1:n_iterations) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x + a, n - x + b)
    samples[i, ] <- c(x, y)
  }
  return(samples)
}
library(Rcpp)

## -----------------------------------------------------------------------------
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix gibbs_sampler_Rcpp(int n_iter, double a, double b, int n){
  NumericMatrix samples(n_iter, 2);
  double x = 0;
  double y = R::runif(0, 1);
  for (int iter = 0; iter < n_iter; iter++) {
    x = R::rbinom(n, y);
    y = R::rbeta(x + a, n - x + b);
    samples(iter, 0) = x;
    samples(iter, 1) = y;
  }
  return samples;
}
')

## -----------------------------------------------------------------------------
library(microbenchmark)
iterations <- 1000
a <- 2
b <- 3
n <- 10
mb_result <- microbenchmark(
  gibbs_sampler_r(iterations, a, b, n),
  gibbs_sampler_Rcpp(iterations, a, b, n),
  times = 10 
)
mb_result

