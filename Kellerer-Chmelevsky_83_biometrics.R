# Kellerer and Chmelevsky


seednum <- (runif(1) * 1000) %>% ceiling
set.seed(seednum)

N1<-10; N2<- 10

odds_true <- 1 ; odds_null <- 1

name_dist <- "exponential"


lambda1<-1; lambda2<-1

Nsim <- 100
Nmonte <- 1000

p_asy <-c()
p_ex <-c()
pb <- txtProgressBar(min = 1, max = Nsim, style = 3)

for(i in 1:Nsim){
  x1test<-rexp(n = N1,rate=lambda1)
  x2test<-rexp(n = N2,rate=lambda2)

  d1 <- rep(Inf, N1)
  d2 <- rep(Inf, N2)
  
  x1delta<- x1test < d1; x2delta<- x2test < d2
  x1test <- pmin(x1test, d1); x2test <- pmin(x2test, d2)
  
  res_lr <- survdiff(Surv(c(x1test,x2test), c(x1delta,x2delta)) ~ c(rep(0,N1),rep(1,N2)))
  z_cor <- (abs( res_lr$obs[1]-res_lr$exp[1])-0.5 )/sqrt(res_lr$var[1,1])
  
  p_monte <- 0
  for(j in 1:Nmonte){
    res_lr_t <- survdiff(Surv(c(x1test,x2test), c(x1delta,x2delta)) ~ as.numeric(is.element(seq(1,N1+N2), sample(N1+N2, N2))))
    z_cor_t <- (abs( res_lr_t$obs[1]-res_lr_t$exp[1])-0.5 )/sqrt(res_lr_t$var[1,1])
    if(z_cor>z_cor_t) p_monte<-p_monte+1/Nmonte
  }
  p_asy <- c(p_asy, pnorm(z_cor, lower.tail = FALSE) )
  p_ex <- c(p_ex, p_monte )
  
  setTxtProgressBar(pb, i) 
}



p_asy > qnorm(0.95)
p_ex  > qnorm(0.95)
