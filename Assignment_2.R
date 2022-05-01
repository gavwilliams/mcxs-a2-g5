# INSTALLING REQUIRED PACKAGES
library(readrba)
library(readabs)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reticulate)
library(mvtnorm)
search(gsl)
set.seed(123456)


# EXERCISE 1 DATA COLLECTION AND TRANSFORMING
############################################################
HOURS               = read_abs(series_id ="A2304428W")
HOURS               = HOURS %>%slice(-c(1:101))
HOURS               = HOURS[c(4,6)]

HHSAVINGS           = read_abs(series_id = "A2323382F")
HHSAVINGS           = HHSAVINGS %>%slice(-c(1:101))
HHSAVINGS           = HHSAVINGS[c(4,6)]

RGNI                = read_abs(series_id = "A2304412C")
RGNI                = RGNI %>%slice(-c(1:101))
RGNI                = RGNI[c(4,6)]
RGNI$value          = log(RGNI$value)

CPI                 = read_abs(series_id = "A3604506F")
CPI                 = CPI %>%slice(-c(1:145))
CPI                 = CPI[c(4,6)]

GDPPERCAP           = read_abs(series_id = "A2304372W")
GDPPERCAP           = GDPPERCAP %>%slice(-c(1:101))
GDPPERCAP           = GDPPERCAP[c(4,6)]

GDP                 = read_abs(series_id = "A2304402X")
GDP                 = GDP %>%slice(-c(1:101))
GDP                 = GDP[c(4,6)]
GDP$value           = log(GDP$value)


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
S                   = 11000
  
TT                  = nrow(Y)
T                   = TT - 1

# Create Y and X
############################################################
y_VEC1              = merge(HOURS, HHSAVINGS, by = 'date')
y_VEC2              = merge(RGNI, CPI, by = 'date')
y_VEC3              = merge(GDP, NET_SAVINGS, by ='date')
y_VEC4              = merge(DOMESTIC_DEMAND, DISP_INC_PC, by ='date')
y_VEC5              = merge(TOT, GDPPERCAP, by ='date')

y1                  = merge(y_VEC1, y_VEC2, by = 'date')
y2                  = merge(y_VEC3, y_VEC4, by = 'date')
y3                  = merge(y1, y2, by ='date')

y                   = merge(y3, y_VEC5, by = 'date') 
y                   = y %>% slice(-c(1:8))
y                   = ts(y[,c(2:11)], start = c(1986,4), frequency = 4, names=c("Hours","Household_Savings","Real_Gross_National_Income","CPI","Gross_Domestic_Product","Net_Savings","Domestic_Demand","Diposable_Income","Terms_of_Trade","GDP_Per_Capita"))

Y                   = ts(y[(p+1):nrow(y),],start = c(1986,4), frequency = 4)
X                   = matrix(1,nrow(Y),1)
X                   = cbind(X,
                            y[5:nrow(y)-1,],
                            y[5:nrow(y)-2,],
                            y[5:nrow(y)-3,],
                            y[5:nrow(y)-4,]
)

#Plot macro data and comment on trends
plot(y)

##GRAPHING CODE
#GermanGNP           = ts(as.data.frame(read.csv("GermanGNP.csv")), start=c(1975,1), frequency=4, names="GermanGNP")
#GerGNP              = matrix(GermanGNP)
#colnames(GerGNP)    = "German GNP"
#plot.ts(GermanGNP, lwd=3, col="purple", main="")

## MLE
A_HAT                 = (1/((t(X)%*%X)))%*%t(X)%*%Y
SIGMA_HAT             = t(Y-X%*%A_HAT)%*%(Y-X%*%A_HAT)/T
round(A_HAT,3)
round(SIGMA_HAT,3)
round(cov2cor(SIGMA_HAT),3)

# PRIOR DISTRIBUTIONS SPECIFICATION
############################################################
KAPPA_P_A             = 100
KAPPA_P_E             = 0.02^2
A_MEAN_PRI            = matrix(0, nrow(A_HAT),ncol(A_HAT))
A_MEAN_PRI[2:11,]     = diag(10)
V_PRIOR               = diag(c(KAPPA_P_E,KAPPA_P_A*((1:p)^(-2))%x%t(rep(1,N))))   # COL SPECIFIC VAR
V_PRIOR_INV           = diag(1/c(KAPPA_P_E,KAPPA_P_A*((1:p)^(-2))%x%t(rep(1,N))))

S_PRIOR               = diag(diag(SIGMA_HAT)) 
NU_PRIOR              = N+1

A_PRIOR               = matrix(0, nrow(A_HAT),ncol(A_HAT))
SIGMA_PRIOR           = diag(10)
SIGMA_PRIOR_INV       = solve(SIGMA_PRIOR)

# SPECIFYING HYPER-PARAMETERS
hyper                 = list(1,1,1,1,A_MEAN_PRI,S_PRIOR,V_PRIOR, NU_PRIOR)
names(hyper)          = c("S_PRIOR_Ka","V_PRIOR_KA","ALPHA_PRIOR_Ke","BETA_PRIOR_Ke","A_MEAN_PRI","S_PRIOR", "V_PRIOR", "NU_PRIOR")


# POSTERIOR DRAWS
posterior_A           = array(NA, c(S, K, N))
posterior_E           = array(NA, c(N, N, S))
posterior_ka          = matrix(NA, S)
posterior_ke          = matrix(NA, S)

GIBBS_SAMPLER = function(S, Y, X , hyper){
  
  for (s in 1:S){
    # HYPER-PARAMETER ESTIMATION -- K_E
    alpha_bar_ke        = hyper$ALPHA_PRIOR_Ke + ((hyper$V_PRIOR_KA*N))/2
    beta_bar_ke         = solve(hyper$BETA_PRIOR_Ke) + 0.5*sum(diag(solve(SIGMA_PRIOR)%*%hyper$S_PRIOR))
    
    # SAMPLING -- K_E
    KAPPA_P_E           = rgamma(1, shape = alpha_bar_ke, scale = beta_bar_ke)
    
    help(rgamma)
    # HYPER-PARAMETER ESTIMATION -- K_A
    s_bar_ka            = sum(diag((solve(SIGMA_PRIOR)%*%t(A_PRIOR-hyper$A_MEAN_PRI))%*%(A_PRIOR-hyper$A_MEAN_PRI))) + hyper$S_PRIOR_Ka
    v_bar_ka            = hyper$V_PRIOR_KA + N*K
    
    # SAMPLING -- K_A
    KAPPA_P_A           = s_bar_ka/rchisq(1, v_bar_ka)
    
    # PARAMETERS OF MVNIW POSTERIOR
    V_bar_inv           = crossprod(X) + solve(KAPPA_P_A*hyper$V_PRIOR)
    V_bar               = solve(V_bar_inv)
    A_bar               = V_bar%*%(t(X)%*%Y + solve(KAPPA_P_A*V_bar)%*%hyper$A_MEAN_PRI)
    
    S_bar               = crossprod(Y) + KAPPA_P_E*hyper$S_PRIOR + t(hyper$A_MEAN_PRI)%*%solve((KAPPA_P_A*hyper$V_PRIOR))%*%hyper$A_MEAN_PRI - t(A_bar)%*%V_bar_inv%*%A_bar
    S_bar_inv           = solve(S_bar)
    #S_bar              = 0.5*(S_bar + t(S_bar))
    #S_bar_chol         = chol(S_bar)
    #S_bar_inv          = backsolve(S_bar_chol, forwardsolve(t(S_bar_chol), diag(N)))
    nu_bar              = T + NU_PRIOR
    
    # SAMPLE -- aux_A & aux_E
    L                   = t(solve(chol(V_bar_inv)))
    SIGMA_POST          = solve(rWishart(1, df = nu_bar, Sigma =S_bar_inv)[,,1])
    draw.norm           = array(rnorm(prod(N*K)),c(K,N))
    A_POST              = A_bar + L%*%draw.norm%*% chol(SIGMA_POST)
    
    round(apply(A_POST,1:2, mean),3)
    
    posterior_ka[s]     = KAPPA_P_A
    posterior_ke[s]     = KAPPA_P_E
    posterior_A[s,,]    = A_POST
    posterior_E[,,s]    = SIGMA_POST
  }
  
  # OUTPUT  
  return(
    list(
      A           = posterior_A,
      E           = posterior_E,
      Ka          = posterior_ka,
      Ke          = posterior_ke
    )
  )
}


SAMPLER_OUTPUT = GIBBS_SAMPLER(S, Y , X, hyper)
plot.ts(SAMPLER_OUTPUT$A[1:S], main = "trace plots", xlab = "", ylab = "A")
plot.ts(SAMPLER_OUTPUT$E[1:S], main = "trace plots", xlab = "", ylab = "E")
plot.ts(SAMPLER_OUTPUT$Ka[1:S], main = "trace plots", xlab = "", ylab = "Ka")
plot.ts(SAMPLER_OUTPUT$Ke[1:S], main = "trace plots", xlab = "", ylab = "Ke")

round(SAMPLER_OUTPUT$A[1,,], mean)
mean(SAMPLER_OUTPUT$A[1,,],1:10)


# report posterior means and sd of parameters
A.E         = apply(A.posterior,1:2,mean)
A.sd        = apply(A.posterior,1:2,sd)
Sigma.E     = apply(Sigma.posterior,1:2,mean)
Sigma.sd    = apply(Sigma.posterior,1:2,sd)


# QUESTION 9 
# SAMPLE -- PREDICTIVE DENSITY 
h               = 4
Y.h             = array(NA,c(h,10,S))
Y.h.m           = array(NA,c(h,10))


FORECAST_FUNCTION = function(h,SAMPLER_OUTPUT){
  for (s in 1:S){
    if (p==1){
      x.Ti          = matrix(Y[(nrow(Y)-p+1):nrow(Y),],nrow=1)
      x.Ti.m        = x.Ti
    } else {
      x.Ti          = Y[(nrow(Y)-p+1):nrow(Y),]
      x.Ti          = x.Ti[p:1,]
      x.Ti.m        = x.Ti
    }
    for (i in 1:h){
      x.T           = c(1,as.vector(t(x.Ti)))
      x.T.m         = c(1,as.vector(t(x.Ti.m)))
      Y.f           = rmvnorm(1, mean = x.T%*%SAMPLER_OUTPUT$A[s,,], sigma=SAMPLER_OUTPUT$E[,,s])
      Y.f.m         = x.T.m%*%A_POST
      if (p==1){
        x.Ti        = Y.f
        x.Ti.m      = Y.f.m
      } else {
        x.Ti        = rbind(Y.f,x.Ti[1:(p-1),])
        x.Ti.m      = rbind(Y.f.m,x.Ti.m[1:(p-1),])
      }
      Y.h[i,,s]     = Y.f[1:2]
      Y.h.m[i,]     = Y.f.m[1:2]
    }
  }
  
  # OUTPUT  
  return(
    list(
      Y_FORECAST    = Y.h,
      Y_M_FORECAST  = Y.h.m
    )
  )
}

FORECAST = FORECAST_FUNCTION(h, SAMPLER_OUTPUT)

# QUESTION 10 
# SAMPLE -- PREDICTIVE DENSITY FUNCTION
