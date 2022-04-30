# QUESTION 6
# INSTALLING REQUIRED PACKAGES
library(readrba)
library(readabs)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reticulate)
# DATA COLLECTION AND TRANSFORMING
############################################################

HOURS               = read_abs(series_id ="A2304428W")
HOURS               = HOURS %>%slice(-c(1:101))
HOURS               = HOURS[c(4,6)]

HHSAVINGS           = read_abs(series_id = "A2323382F")
HHSAVINGS           = HHSAVINGS %>%slice(-c(1:101))
HHSAVINGS           = HHSAVINGS[c(4,6)]

RGDP                = read_abs(series_id = "A2304414J")
RGDP                = RGDP %>%slice(-c(1:101))
RGDP                = RGDP[c(4,6)]

CPI                 = read_abs(series_id = "A3604506F")
CPI                 = CPI %>%slice(-c(1:145))
CPI                 = CPI[c(4,6)]

VACANCIES           = read_abs(series_id = "A590698F")
VACANCIES           = VACANCIES %>%slice(-c(1:22))
VACANCIES           = VACANCIES[c(4,6)]

VALUE_ADDED         = read_abs(series_id = "A3606058X")
VALUE_ADDED         = VALUE_ADDED %>%slice(-c(1:101))
VALUE_ADDED         = VALUE_ADDED[c(4,6)]

GDP_PHW             = read_abs(series_id = "A2304424L")
GDP_PHW             = GDP_PHW %>%slice(-c(1:101))
GDP_PHW             = GDP_PHW[c(4,6)]

NET_SAVINGS         = read_abs(series_id = "A2304424L")
NET_SAVINGS         = NET_SAVINGS %>%slice(-c(1:101))
NET_SAVINGS         = NET_SAVINGS[c(4,6)]

DOMESTIC_DEMAND     = read_abs(series_id = "A2304198A")
DOMESTIC_DEMAND     = DOMESTIC_DEMAND %>%slice(-c(1:101))
DOMESTIC_DEMAND     = DOMESTIC_DEMAND[c(4,6)]

DISP_INC_PC         = read_abs(series_id = "A2304416L")
DISP_INC_PC         = DISP_INC_PC %>% slice(-c(1:101))
DISP_INC_PC         = DISP_INC_PC[c(4,6)]

TOT                 = read_abs(series_id = "A2304400V")
TOT                 = TOT %>% slice(-c(1:101))
TOT                 = TOT[c(4,6)]

# setup
############################################################
N                   = 10
p                   = 4
K                   = 1+p*N
S                   = 50000
h                   = #TBD
  
TT                  = nrow(Y)
T                   = TT - 1

# Create Y and X
############################################################
y_VEC1              = merge(HOURS, HHSAVINGS, by = 'date')
y_VEC2              = merge(RGDP, CPI, by = 'date')
y_VEC3              = merge(GDP_PHW, NET_SAVINGS, by ='date')
y_VEC4              = merge(DOMESTIC_DEMAND, DISP_INC_PC, by ='date')
y_VEC5              = merge(TOT, VALUE_ADDED, by ='date')

y1                  = merge(y_VEC1, y_VEC2, by = 'date')
y2                  = merge(y_VEC3, y_VEC4, by = 'date')
y3                  = merge(y1, y2, by ='date')

y                   = merge(y3, y_VEC5, by = 'date') 
y                   = y %>% slice(-c(1:39))
y                   = ts(y[,c(2:11)], start = c(1994,3), frequency = 4, names=c("VAR1","VAR2","VAR3","VAR4","VAR5","VAR6","VAR7","VAR8","VAR9","VAR10"))

Y                   = ts(y[(p+1):nrow(y),],start = c(1994,3), frequency = 4)
X                   = matrix(1,nrow(Y),1)
X                   = cbind(X,
                            y[5:nrow(y)-1,],
                            y[5:nrow(y)-2,],
                            y[5:nrow(y)-3,],
                            y[5:nrow(y)-4,]
                            )
          
##GRAPHING CODE
#GermanGNP           = ts(as.data.frame(read.csv("GermanGNP.csv")), start=c(1975,1), frequency=4, names="GermanGNP")
#GerGNP              = matrix(GermanGNP)
#colnames(GerGNP)    = "German GNP"
#plot.ts(GermanGNP, lwd=3, col="purple", main="")

## MLE
A_HAT                 = (1/(t(X)%*%X))%*%t(X)%*%Y
SIGMA_HAT             = t(Y-X%*%A_HAT)%*%(Y-X%*%A_HAT)
round(A_HAT,3)
round(SIGMA_HAT,3)
round(cov2cor(SIGMA_HAT),3)


# PRIOR DISTRIBUTION COPIED FROM L10
############################################################
KAPPA_P_A   = 0.02^2
KAPPA_P_E   = 100
A_PRIOR     = matrix(0,K,N)
V_PRIOR     = diag(c(KAPPA_P_E,KAPPA_P_A*((1:p)^(-2))%x%rep(1,N)))
V_PRIOR_INV = diag(1/c(KAPPA_P_E,KAPPA_P_A*((1:p)^(-2))%x%rep(1,N)))
S_PRIOR     = diag(diag(SIGMA_HAT))
NU_PRIOR    = N+1



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


# STARTING VALUES 
aux_ka          = 2
aux_ke          = 2
aux_A           = matrix(0, K, N)
aux_E           = rep(1,N)
aux_E           = diag(aux_E)


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
  V_bar_inv       = crossprod(X) + 1/(aux_ka*hyper[7])
  V_bar_inv_chol  = chol(V_bar_inv)
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

