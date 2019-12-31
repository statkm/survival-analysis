# Heinze, et al. 2003 


seednum <- (runif(1) * 1000) %>% ceiling
set.seed(seednum)

N1<-10; N2<- 100
odds_true <- 1; odds_null <- 1
lambda<- .04 ; 
lambda1<-lambda*odds_true; lambda2<-lambda

x1test<-rexp(n = N1,rate = lambda1)
x2test<-rexp(n = N2,rate = lambda2)

dd1 <- rexp(n = N1,rate = 0.0000001)
dd2 <- rexp(n = N1,rate = 0.04)



x1delta<- x1test < dd1; x2delta<- x2test < dd2
x1test <- pmin(x1test, dd1); x2test <- pmin(x2test, dd2)

df_surv <- rbind(data.frame(time=x1test, cens=x1delta, gp=1), data.frame(time=x2test, cens=x2delta, gp=2) )



logrank_p <- function(N_monte, df_surv){
  teststat <- numeric(N_monte)
  n_row <- dim(df_surv)[1]
  for(i in 1:N_monte){
    perm_ind <- sample(1:n_row, replace = FALSE)
    
    df_surv_p <- df_surv
    df_surv_p$time_p <- df_surv_p[perm_ind, ]$time
    df_surv_p$cens_p <- df_surv_p[perm_ind, ]$cens
    
    df_surv_p$x_p <- quantile(df_surv_p[which(df_surv_p$cens==1),]$time, runif(n=n_row), type=1)*(1-df_surv_p$cens_p) + df_surv_p$time_p*df_surv_p$cens_p
    df_surv_p$c_p <- 1
    if(dim(df_surv_p[which(df_surv_p$cens_p==0 & df_surv_p$x_p > max(df_surv_p[which(df_surv_p$cens==1),]$time) ),])[1]>0) df_surv_p[which(df_surv_p$cens_p==0 & df_surv_p$x_p > max(df_surv_p[which(df_surv_p$cens==1),]$time) ),]$c_p<-1
    df_surv_p$y_p <- 0
    
    n_row_g <- dim(df_surv_p[which(df_surv_p$gp==1), ])[1]
    df_surv_p[which(df_surv_p$gp==1), ]$y_p <- quantile(df_surv_p[which(df_surv_p$cens==0 & df_surv_p$gp==1 ),]$time, runif(n=n_row_g), type=1)
    n_row_g <- dim(df_surv_p[which(df_surv_p$gp==2), ])[1]
    df_surv_p[which(df_surv_p$gp==2), ]$y_p <- quantile(df_surv_p[which(df_surv_p$cens==0 & df_surv_p$gp==2 ),]$time, runif(n=n_row_g), type=1)
    
    df_surv_p$time2_p <- pmin(df_surv_p$x_p, df_surv_p$y_p)
    df_surv_p$cens2_p <- pmax(df_surv_p$x_p<df_surv_p$y_p, (df_surv_p$x_p==df_surv_p$y_p & df_surv_p$c_p==1))
    
    result_logrank <- survdiff(Surv(time2_p, cens2_p) ~ gp, data=df_surv_p, rho = 0)
    teststat[i] <- ( result_logrank$obs[2]-result_logrank$exp[2]) /sqrt(result_logrank$var[1,1])
  }
  return(teststat)
}


res_lr_p <- logrank_p(1000, df_surv)

result_logrank <- survdiff(Surv(time, cens) ~ gp, data=df_surv, rho = 0)
lr_test <- ( result_logrank$obs[2]-result_logrank$exp[2]) /sqrt(result_logrank$var[1,1])

hist(res_lr_p)
abline(v=lr_test,lwd=5,col="red")

mean(res_lr_p>=lr_test)
pnorm(lr_test,lower.tail = FALSE)


# event distribution
result   <- survfit(Surv(time,cens)~gp, df_surv, type="kaplan-meier")
plot(result, col=1:4, xlab="time", ylab="(%)", main="time to event")

result   <- survfit(Surv(time,cens)~1, df_surv, type="kaplan-meier")
plot(result, col=1:4, xlab="time", ylab="(%)", main="time to event")

# follow-up distribution
df_surv$cens_f <- 1 - df_surv$cens
result_f <- survfit(Surv(time,cens_f)~gp, df_surv, type="kaplan-meier")
plot(result_f, col=1:4, xlab="time", ylab="(%)", main="follow-up time")
