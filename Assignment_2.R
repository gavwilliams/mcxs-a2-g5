# QUESTION 6
# INSTALLING REQUIRED PACKAGES
library(readrba)
library(readabs)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)

# CREATING Y (K x (TxN)) VARIABLE =

HOURS               = read_abs(series_id ="A2304428W")
HHSAVINGS           = read_abs(series_id = "A2323382F") 
GDP                 = read_abs(series_id = "A2304370T")
CPI                 = read_abs(series_id = "A3604506F")
VACANCIES           = read_abs(series_id = "A590698F") 
EMPLOYMENT          = read_abs(series_id = "A84932381A") 

HOURS               = ts(HOURS[c(6)], start=c(1959,3), frequency = 4, names=c("Hours Worked"))

y                   = full_join(HOURS, HHSAVINGS[c(6)], GDP[c(6)], CPI[c(6)],VACANCIES[c(6)], EMPLOYMENT[c(6)])

HHSAVINGS           = read_abs(series_id = "A2323382F") 
HHSAVINGS           = ts(HHSAVINGS[c(6)], start=c(1959,3), frequency = 4, names=c("HHSAVINGS"))

GDP                 = read_abs(series_id = "A2304370T")
GDP                 = ts(GDP[c(6)], start=c(1959,3), frequency = 4, names=c("GDP"))

CPI                 = read_abs(series_id = "A3604506F")
CPI                 = ts(CPI[c(6)], start=c(1959,3), frequency = 4, names=c("CPI"))

VACANCIES           = read_abs(series_id = "A590698F") 
VACANCIES           = ts(VACANCIES[c(6)], start=c(1959,3), frequency = 4, names=c("CPI"))


EMPLOYMENT          = read_rba(series_id = "A84932381A") 
EMPLOYMENT          = ts(EMPLOYMENT[c(6)], start=c(1959,3), frequency = 4, names=c("CPI"))





VAR6                = read_rba(series_id = "GLFSURSA") 
VAR7                = read_rba(series_id = "GLFSURSA") 
VAR8                = read_rba(series_id = "GLFSURSA") 
VAR9                = read_rba(series_id = "GLFSURSA") 
VAR10               = read_rba(series_id = "GLFSURSA") 

GermanGNP           = ts(as.data.frame(read.csv("GermanGNP.csv")), start=c(1975,1), frequency=4, names="GermanGNP")
GerGNP              = matrix(GermanGNP)
colnames(GerGNP)    = "German GNP"
plot.ts(GermanGNP, lwd=3, col="purple", main="")


# INPUTTING DATA 
X                   = as.matrix(NA, c(1:100))
Y                   = as.matrix(rw[1:100])

# SPECIFYING HYPER-PARAMETERS OF PRIOR DISTRIBUTION   
priors = c(1,2,3,4,5,6,7)
names(priors) = c("S_BAR(Ka)","V_BAR(KA)","ALPHA_BAR(Ke)","BETA_BAR(Ke)","A_BAR","S_BAR", "V_BAR")

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
hyper             = c(1,2,3,4,5,6,7)
names(hyper)      = c("S_BAR(Ka)","V_BAR(KA)","ALPHA_BAR(Ke)","BETA_BAR(Ke)","A_uBAR","S_BAR", "V_BAR")

T = 170
K = 10

# STARTING VALUES 
aux_A           = rep(0, K)
aux_E           = matrix(0, K, K)
aux_ka          = 2
aux_ke          = 2

# POSTERIOR DRAWS
posterior_A     = matrix(NA, S, K)
posterior_E     = matrix(NA, S, K)
posterior_ka    = matrix(NA, S)
posterior_ke    = matrix(NA, S)

for (s in 1:S) {
  
  # SAMPLE Ke
  alpha_bar_ke    = hyper[3] - ((hyper[7]*K))/2
  beta_bar_ke     = hyper[4] + as.numeric(0.5*trace(solve(aux_E)%*%hyper[6]))
  
  # SAMPLE Ka
  s_bar_ka        = as.numeric(trace(solve(aux_E)%*%t(aux_A-hyper[5])%*%(aux_A-hyper[5])) + hyper[1])
  v_bar_ka        = hyper[2] + N*K
  
  # PARAMETERS OF MVNIW POSTERIOR
  V_bar_inv       = crossprod(X) + solve(aux_ka*hyper[7])
  V_bar           = as.numeric(solve(V_bar_inv))
  A_bar           = as.numeric(V_bar %*% (t(X)%*%Y) + solve(aux_ka*V_bar) %*% hyper[5])
  s_bar           = as.numeric(crossprod(Y) + aux_ke*hyper[6] + t(aux_A)%*%(aux_ka*hyper[7])%*%aux_A + t(A_bar)%*%solve(V_bar)%*%A_bar )
  
  # SAMPLE -- aux_ka
  aux_ka          = s_bar_ka/rchisq(1, v_bar_ka)
  posterior_ka    = aux_ka
  
  # SAMPLE -- aux_ke
  aux_ke          = rgamma(1, shape = alpha_bar_ke, scale = 1/beta_bar_ke)
  
  
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

