# INSTALLING REQUIRED PACKAGES
library(readrba)
library(readabs)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reticulate)
library(mvtnorm)
library(HDInterval)
search(gsl)
set.seed(123456)


# EXERCISE 1 DATA COLLECTION AND TRANSFORMING
############################################################
##VAR1
HOURS               = read_abs(series_id ="A2304428W")
HOURS               = HOURS %>%slice(-c(1:101))
HOURS               = HOURS[c(4,6)]

##VAR2
HHSAVINGS           = read_abs(series_id = "A2323382F")
HHSAVINGS           = HHSAVINGS %>%slice(-c(1:101))
HHSAVINGS           = HHSAVINGS[c(4,6)]

##VAR3
EMPLOYMENT          = read_abs(series_id = "A2454521V")
EMPLOYMENT          = EMPLOYMENT %>%slice(-c(1:109))
EMPLOYMENT          = EMPLOYMENT[c(4,6)]

##VAR4
CPI_CHANGE_Q        = read_abs(series_id = "A2325850V")
CPI_CHANGE_Q        = CPI_CHANGE_Q %>%slice(-c(1:153))
CPI_CHANGE_Q        = CPI_CHANGE_Q[c(4,6)]

##VAR5
REAL_LABOUR_COST    = read_abs(series_id = "A2433071F")
REAL_LABOUR_COST    = REAL_LABOUR_COST %>%slice(-c(1:5))
REAL_LABOUR_COST    = REAL_LABOUR_COST[c(4,6)]

##VAR6
GDP_PCT_CHANGE      = read_abs(series_id = "A2304370T")
GDP_PCT_CHANGE      = GDP_PCT_CHANGE %>%slice(-c(1:101))
GDP_PCT_CHANGE      = GDP_PCT_CHANGE[c(4,6)]

##VAR7
PUBLIC_CAPITAL      = read_abs(series_id = "A2454459T")
PUBLIC_CAPITAL      = PUBLIC_CAPITAL %>%slice(-c(1:109))
PUBLIC_CAPITAL      = PUBLIC_CAPITAL[c(4,6)]

##VAR8
INVENTORIES         = read_abs(series_id = "A3538852F")
INVENTORIES         = INVENTORIES %>%slice(-c(1:5))
INVENTORIES         = INVENTORIES[c(4,6)]

##VAR9
DISP_INC_PC         = read_abs(series_id = "A2304416L")
DISP_INC_PC         = DISP_INC_PC %>% slice(-c(1:101))
DISP_INC_PC         = DISP_INC_PC[c(4,6)]

##VAR10
TOT                 = read_abs(series_id = "A2304400V")
TOT                 = TOT %>% slice(-c(1:101))
TOT                 = TOT[c(4,6)]


# setup
############################################################
set.seed(2711197)
N                   = 10
p                   = 4
K                   = 1+p*N
<<<<<<< Updated upstream
S                   = 1100
=======
S                   = 11000
>>>>>>> Stashed changes
  
TT                  = nrow(Y)
T                   = TT - 1

# Create Y and X
############################################################
y_VEC1              = merge(HOURS, HHSAVINGS, by = 'date')
y_VEC2              = merge(EMPLOYMENT, CPI_CHANGE_Q, by = 'date')
y_VEC3              = merge(REAL_LABOUR_COST, GDP_PCT_CHANGE, by ='date')
y_VEC4              = merge(PUBLIC_CAPITAL, INVENTORIES, by ='date')
y_VEC5              = merge(DISP_INC_PC, TOT, by ='date')

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
plot(Y)


## MLE
<<<<<<< Updated upstream
A_HAT                 = solve(t(X)%*%X)%*%t(X)%*%Y
=======
A_HAT                 = solve((t(X)%*%X))%*%t(X)%*%Y
>>>>>>> Stashed changes
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

s.ols       = rep(NA,N)
for (n in 1:N){
  s.ols[n]  = var(ar(x=Y[,n], aic=FALSE, order.max=8, method="ols")$resid[9:T])
}

#S_PRIOR               = diag(diag(SIGMA_HAT)) 
S_PRIOR               = diag(s.ols)
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
    
    # HYPER-PARAMETER ESTIMATION -- K_A
    s_bar_ka            = sum(diag((solve(SIGMA_PRIOR)%*%t(A_PRIOR-hyper$A_MEAN_PRI))%*%(A_PRIOR-hyper$A_MEAN_PRI))) + hyper$S_PRIOR_Ka
    v_bar_ka            = hyper$V_PRIOR_KA + N*K
    
    # SAMPLING -- K_A
    KAPPA_P_A           = s_bar_ka/rchisq(1, v_bar_ka)
    
    # PARAMETERS OF MVNIW POSTERIOR
    V_bar_inv           = crossprod(X) + diag(1/diag(hyper$V_PRIOR))*(1/KAPPA_P_A)
    V_bar               = solve(V_bar_inv)
    A_bar               = V_bar%*%(t(X)%*%Y + (diag(1/diag(V_bar))*(1/KAPPA_P_A))%*%hyper$A_MEAN_PRI)
    
    S_bar               = KAPPA_P_E*hyper$S_PRIOR +crossprod(Y)+ t(hyper$A_MEAN_PRI)%*%((diag(1/diag(hyper$V_PRIOR)))*KAPPA_P_A)%*%hyper$A_MEAN_PRI - t(A_bar)%*%V_bar_inv%*%A_bar
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

plot.ts(SAMPLER_OUTPUT$A[1000:S], main = "trace plots", xlab = "", ylab = "A")
plot.ts(SAMPLER_OUTPUT$E[1000:S], main = "trace plots", xlab = "", ylab = "E")
plot.ts(SAMPLER_OUTPUT$Ka[1000:S], main = "trace plots", xlab = "", ylab = "Ka")
plot.ts(SAMPLER_OUTPUT$Ke[1000:S], main = "trace plots", xlab = "", ylab = "Ke")

round(SAMPLER_OUTPUT$A[1,,], mean)
mean(SAMPLER_OUTPUT$A[1,,],1:10)


# report posterior means and sd of parameters
A.E         = apply(A_POST,1:2,mean)
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
h               = 4
Y.h             = array(NA,c(h,2,S))
Y.h.m           = array(NA,c(h,2))
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
# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

mcxs1.rgb   = col2rgb(mcxs1)
mcxs1.shade1= rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=50, maxColorValue=255)
mcxs2.rgb   = col2rgb(mcxs2)
mcxs2.shade1= rgb(mcxs2.rgb[1],mcxs2.rgb[2],mcxs2.rgb[3], alpha=50, maxColorValue=255)
mcxs3.rgb   = col2rgb(mcxs3)
mcxs3.shade1= rgb(mcxs3.rgb[1],mcxs3.rgb[2],mcxs3.rgb[3], alpha=50, maxColorValue=255)


plot.ts(FORECAST$Y_FORECAST[4,1,1000:S], main = "trace plots", xlab = "", ylab = "FORECAST - HOURS WORKED")
hist(FORECAST$Y_FORECAST[4,1,1000:S], main = "histograms", xlab = "",ylab = "FORECAST - HOURS WORKED")

Y.h         = FORECAST$Y_FORECAST
Y.h.m       = FORECAST$Y_M_FORECAST


hour_point_forecast         = apply(Y.h[,1,],1,mean)
hour_interval_forecast      = apply(Y.h[,1,],1,hdi,credMass=0.90)
hour_range                  = range(y[,2],hour_interval_forecast)

pdf(file="HOURS_Forecast.pdf", width=15,height=6)

plot(1:(length(y[,1])+h),c(y[,1],hour_point_forecast), type="l", ylim=rgdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=mcxs1)
axis(1,c(3,43,83,123,163,203,243,nrow(y),nrow(y)+h),c("1960","1970","1980","1990","2000","2010","2020","",""), col=mcxs1)
axis(2,c(hour_range[1],mean(hour_range),hour_range[2]),c("","HOURS",""), col=mcxs1)
abline(v=246, col=mcxs1)
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[21:1]),
        c(y[90,1],hour_interval_forecast[1,],hour_interval_forecast[2,20:1],y[90,1]),
        col=mcxs1.shade1, border=mcxs1.shade1)






Y_recent        = ts(rbind(Y[92:117,1:2],matrix(NA,h,2)),start=c(2009,3),frequency=4)
interval_fs1    = apply(Yh1,1,hdi,credMass=0.68)
lims_hours      = range(Y_recent[1:26,1],Yhm1[,1])

plot(1:nrow(Y_recent),Y_recent[,1], type="l", ylim=lims_hours, axes=FALSE, xlab="", ylab="", main=expression(HOURS[t]), lwd=2, col=mcxs1)
axis(1,c(seq(from=1, to=30, by=4),30),c("","2010","","2012","","2014","","2016",""), col=mcxs1)
axis(2,c(lims_hours[1],0,lims_hours[2]),round(c(lims_hours[1],0,lims_hours[2]),2), col=mcxs1)
polygon(c(26:30,30:26), c(Y_recent[26,1],interval_fs1[1,,1],interval.fs2[2,h:1,1],Y_recent[26,1]), col=mcxs1.shade1,border=mcxs1.shade1)
abline(h=0, col=mcxs1)
lines(26:30, c(Y_recent[26,1],Yhm1[,1]),lwd=2,col=mcxs3)
lines(26:30, c(Y_recent[26,1],Yhm2[,1]),lwd=2,col=mcxs2)
lines(26:30, c(Y.recent[26,1],Yhm4[,1]),lwd=2,col=mcxs1)



dim(FORECAST$Y_FORECAST)
