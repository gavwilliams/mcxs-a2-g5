# EXERCISE 1
# INSTALLING REQUIRED PACKAGES
library(readrba)
library(readabs)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reticulate)

rm(list=ls()) #remove current memory of the environment


# DOWNLOADING, TRANSFORMING AND PLOTTING TEN VARIABLES
HOURS               = as.data.frame(read_abs(series_id ="A2304428W"))
HOURS               = HOURS %>%slice(-c(1:152))
HOURS               = ts(HOURS[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("Hours Worked"))
plot(HOURS)

HOURLYPAY           = as.data.frame(read_abs(series_id ="A2713849C"))
HOURLYPAY           = ts(HOURLYPAY[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("Hourly Pay"))
plot(HOURLYPAY) 

HHSAVINGS           = as.data.frame(read_abs(series_id = "A2323382F"))
HHSAVINGS           = HHSAVINGS %>%slice(-c(1:152))
HHSAVINGS           = ts(HHSAVINGS[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("HHSAVINGS"))
HHSAVINGS
plot(HHSAVINGS)

RGNI                = as.data.frame(read_abs(series_id = "A2304412C"))
RGNI                = RGNI %>%slice(-c(1:152))
RGNI                = ts(RGNI[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("Real Gross National Income"))
plot(RGNI)

GDP                 = as.data.frame(read_abs(series_id = "A2304370T"))
GDP                 = GDP %>%slice(-c(1:152))
GDP                 = ts(GDP[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("Gross Domestic Product: Chain VOlume"))
plot(GDP)


GDPPERCAP           = as.data.frame(read_abs(series_id = "A2304372W"))
GDPPERCAP           = GDPPERCAP %>%slice(-c(1:152))
GDPPERCAP           = ts(GDPPERCAP[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("Gross Domestic Product Per Capita: Chain VOlume"))
plot(GDPPERCAP)

CPI                 = as.data.frame(read_abs(series_id = "A3604506F"))
CPI                 = CPI %>%slice(-c(1:197))
CPI                 = ts(CPI[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("CPI"))
plot(CPI)

VACANCIES           = as.data.frame(read_abs(series_id = "A590698F"))
VACANCIES           = VACANCIES %>%slice(-c(1:74))
VACANCIES           = ts(VACANCIES[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("VACANCIES"))
plot(VACANCIES) # THERE ARE SOME N/A VALUES - NEED TO FIND BETTER SOURCE

EMPLOYMENT          = as.data.frame(read_abs(series_id = "A84932381A"))
EMPLOYMENT          = EMPLOYMENT %>%slice(-c(1:51,150))
EMPLOYMENT          = ts(EMPLOYMENT[c(6)], start=c(1997,3), end=c(2021,4), frequency = 4, names=c("EMPLOYMENT"))
plot(EMPLOYMENT)


TOT                 = as.data.frame(read_abs(series_id = "A2304400V"))
TOT                 = TOT %>% slice(-c(1:152))
TOT                 = ts(TOT[c(6)], start=c(1997,3), end=c(2021,4),frequency = 4, names=c("Terms of Trade"))
plot(TOT)





# CREATING Y (K x (TxN)) VARIABLE

MACRODATA           = as.data.frame(cbind(HOURS, GDP, GDPPERCAP, RGNI, HHSAVINGS, HOURLYPAY, CPI, EMPLOYMENT, VACANCIES, TOT),10)

TT = nrow(MACRODATA)
T = TT - 1
Y = as.matrix(MACRODATA[,2:TT])
X_rgdp = cbind(rep(1, T), rgdp[1:T])


nrow(Y)
Y_VEC  = ts(Y[(4+1):nrow(Y),], start=c(1988,4), frequency=4)
X_VEC  = matrix(1,nrow(Y_VEC),1)

p = 4
for (i in 1:p){
  X    = cbind(X, Y[(p+1):nrow(Y)-i,])
}



# INPUTTING DATA 
x                   = matrix (1, K, 1)
X                   = matrix(1, T, K)
Y                   = matrix(1, T, N)