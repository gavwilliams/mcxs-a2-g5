# QUESTION 6
# INSTALLING REQUIRED PACKAGES
library(readrba)
library(readabs)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reticulate)
# CREATING Y (K x (TxN)) VARIABLE =
HOURS               = as.data.frame(read_abs(series_id ="A2304428W"))
HOURS               = HOURS %>%slice(-c(1:104))
HOURS               = ts(HOURS[c(6)], start=c(1984,4), frequency = 4, names=c("Hours Worked"))

HHSAVINGS           = as.data.frame(read_abs(series_id = "A2323382F")) 
HHSAVINGS           = HHSAVINGS %>%slice(-c(1:104))
HHSAVINGS           = ts(HHSAVINGS[c(6)], start=c(1984,4), frequency = 4, names=c("HHSAVINGS"))

RGDP                = as.data.frame(read_abs(series_id = "A2304414J"))
RGDP                = RGDP %>%slice(-c(1:101))
RGDP                = ts(RGDP$value, start=c(1984,4), frequency = 4, names=c("GDP"))

CPI                 = as.data.frame(read_abs(series_id = "A3604506F"))
CPI                 = CPI %>%slice(-c(1:145))
CPI                 = ts(CPI[c(6)], start=c(1984,4), frequency = 4, names=c("CPI"))

VACANCIES           = as.data.frame(read_abs(series_id = "A590698F"))
VACANCIES           = VACANCIES %>%slice(-c(1:22))
VACANCIES           = ts(VACANCIES[c(6)], start=c(1984,4), frequency = 4, names=c("VACANCIES"))

EMPLOYMENT          = as.data.frame(read_abs(series_id = "A84423043C"))
EMPLOYMENT          = EMPLOYMENT %>%slice(-c(1:82))
EMPLOYMENT_M        = ts(EMPLOYMENT$value, start=c(1984,12), frequency = 12, names=c("EMPLOYMENT"))
EMPLOYMENT_M_Q      = aggregate.ts(EMPLOYMENT_M, nfrequency = 4)

WAGES               = as.data.frame(read.csv("TotalComp.csv",header=TRUE))
WAGES               = WAGES %>%slice(-c(1:101))
wAGES               = as.date(WAGES$date, format(%M/%Y))
WAGES               = ts(WAGES$value, start=c(1984,4), frequency = 4, names=c("WAGES"))

RETAIL_SALES        = as.data.frame(read.csv("RetSal.csv",header=TRUE))
RETAIL_SALES        = RETAIL_SALES %>%slice(-c(1:179))
RETAIL_SALES_M      = ts(RETAIL_SALES$value, start=c(1984,12), frequency = 12, names=c("RETAIL_SALES_Q"))
RETAIL_SALES_Q      = aggregate.ts(RETAIL_SALES_M, nfrequency = 4)

BANKBILL_90D        = as.data.frame(read.csv("90_Day_Bank_Bills.csv",header=TRUE))
BANKBILL_90D        = BANKBILL_90D %>% slice(-c(1:179))
BANKBILL_90D_M      = ts(BANKBILL_90D$value, start=c(1984,12), frequency = 12, names=c("BANKBILL_90D"))
BANKBILL_90D_M_Q    = aggregate.ts(BANKBILL_90D_M, nfrequency = 4)

TOT                 = as.data.frame(read_abs(series_id = "A2304400V"))
TOT                 = TOT %>% slice(-c(1:101))
TOT                 = ts(TOT$value, start=c(1984,4), frequency = 4, names=c("TOT"))

Y1                  = inner_join(HOURS[c(4,6)], HHSAVINGS[c(4,6)], by = 'date')
Y1.1                = inner_join(HOURS[c(4,6)], HHSAVINGS[c(4,6)], by = 'date')

Y2                  = inner_join(RGDP[c(4,6)], CPI[c(4,6)], by = 'date')
Y2.1                = inner_join(RGDP[c(4,6)], CPI[c(4,6)], by = 'date')

Y3 = inner_join(Y1, Y1.1, by = 'date')
Y4 = inner_join(Y2, Y2.1, by = 'date')

Y5 = inner_join(Y3, Y4, by = 'date')

Y  = inner_join(Y5, Y1, by = 'date')
Y  = Y %>% slice(-c(1:8))
Y  = ts(Y[c(2:11)], start=c(1986,4), frequency = 4, names=c("VAR1","VAR2","VAR3","VAR4","VAR5","VAR6","VAR7","VAR8","VAR9","VAR10"))
Y   = Y[,c(1:10)] # GETS RID OF TIME

nrow(Y)
Y_VEC  = ts(Y[(4+1):nrow(Y),], start=c(1988,4), frequency=4)
X_VEC  = matrix(1,nrow(Y_VEC),1)

p = 4
for (i in 1:p){
  X    = cbind(X, Y[(p+1):nrow(Y)-i,])
}

GermanGNP           = ts(as.data.frame(read.csv("GermanGNP.csv")), start=c(1975,1), frequency=4, names="GermanGNP")
GerGNP              = matrix(GermanGNP)
colnames(GerGNP)    = "German GNP"
plot.ts(GermanGNP, lwd=3, col="purple", main="")


# INPUTTING DATA 
x                   = matrix (1, K, 1)
X                   = matrix(1, T, K)
Y                   = matrix(1, T, N)

# SPECIFYING HYPER-PARAMETERS OF PRIOR DISTRIBUTION   
S = c(100,1000)
a0 = 2                      # TO BE CHANGED?
e0 = 2                      # TO BE CHANGED?
intial_kappa_a       = a0
intial_kappa_e       = e0

posterior1           = matrix(NA, sum(S), 1)
posterior2           = matrix(NA, sum(S), 1)
posterior3           = matrix(NA, sum(S), 2)
colnames(posterior1) = c("A")
colnames(posterior2) = c("E")
colnames(posterior3) = c("Ka","Ke")

# SPECIFYING HYPER-PARAMETERS
hyper             = c(1,2,3,3,1,3,1)
names(hyper)      = c("S_BAR(Ka)","V_BAR(KA)","ALPHA_BAR(Ke)","BETA_BAR(Ke)","A_uBAR","S_BAR", "V_BAR")

T = 170             #=nrow(Y)
K = 1+(4*N)         #=ncol(X)
N = 10

# STARTING VALUES 
aux_ka          = 2
aux_ke          = 2
aux_A           = matrix(0, K, N)
aux_E           = diag(1, N)
aux_E           = rep(NA,N)
for (n in 1:N){
  aux_E[n]  = var(ar(x=Y[,n], aic=FALSE, order.max=8, method="ols")$resid[9:T])
}

# POSTERIOR DRAWS
posterior_A     = matrix(NA, S, K)
posterior_E     = matrix(NA, S, K)
posterior_ka    = matrix(NA, S)
posterior_ke    = matrix(NA, S)

for (s in 1:S) {
  
  # SAMPLE Ke
  alpha_bar_ke    = hyper[3] - ((hyper[7]*N))/2
  beta_bar_ke     = hyper[4] + 0.5*sum(diag(solve(aux_E)*hyper[6]))
  
  # SAMPLE Ka
  s_bar_ka        = as.numeric(sum(diag((solve(aux_E)%*%t(aux_A-hyper[5])%*%(aux_A-hyper[5])) + hyper[1])))
  v_bar_ka        = hyper[2] + N*K
  
  # PARAMETERS OF MVNIW POSTERIOR
  V_bar_inv       = as.numeric(crossprod(X) + solve(aux_ka*hyper[7]))
  V_bar           = as.numeric(solve(V_bar_inv))
  A_bar           = as.numeric(V_bar %*% (t(X)%*%Y) + solve(aux_ka*V_bar) %*% hyper[5])
  s_bar           = as.numeric(crossprod(Y) + aux_ke*hyper[6] + t(aux_A)*solve((aux_ka*hyper[7]))*aux_A + t(A_bar)%*%solve(V_bar)%*%A_bar )
  
  # SAMPLE -- aux_ka
  aux_ka          = s_bar_ka/rchisq(1, v_bar_ka)
  posterior_ka    = aux_ka
  
  # SAMPLE -- aux_ke
  aux_ke          = rgamma(1, shape = alpha_bar_ke, scale = 1/beta_bar_ke)
  posterior_ke    = aux_ke
  
  
  # SAMPLE -- aux_A & aux_E
  aux_E           = rWishart(S, df = V_bar, Sigma = s_bar)
  aux_E           = apply(aux_E, 3, solve)
  aux_E           = array(aux_E, c(N,N,S))
  
  aux_A           = array(rnorm(prod(c(dim(A_bar),S))), c(dim(A_bar), S)
                          L           = t(chol(V_bar))
                          
                          for (s in 1:S){
                            aux_A[,,s]    = A_bar + L %*%aux_A[,,s]%*%chol(aux_E[,,s])
                          }
                          
                          round(apply(aux_A,1:2,mean),3)
                          
                          posterior_ka[s] = aux_ka
                          posterior_ke[s] = aux_ke
                          posterior_A[s,] = aux_A
                          posterior_E[s,] = aux_E
}

# OUTPUT  
return(
  list(
    A           = posterior_A
    E           = posterior_E
    Ka          = posterior_ka
    Ke          = posterior_ke
  )
)
}


# QUESTION 9 
# SAMPLE -- PREDICTIVE DENSITY 
h               = 100
Y.h             = array(NA,c(h,N.S))

for (s in 1:S){
  x.Ti          = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti          = x.Ti[4:1,]                  #MAY NEED TO CHANGE DIMENSIONS
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1,
                          mean  = x.T%*%posterior_A[,,s],
                          sigma = posterios_E[,,s])
    x.Ti        = rbind(Y.h[i,,s], x.Ti[1:3,])
  }
}

# 1 PERIOD AHEAD -- JOINT PREDICTIVE DENSITY
limits.VAR1     = range(Y.h[1,1,])
bands           = 100

