cat("\f")
rm(list=ls())
#pakotnes:
# readxl, forecast, tseries, astsa, urca, HiddenMarkov, depmixS4, ggplo2
library(readxl)

setwd("C:\\Users\\GS00183S\\Downloads\\")
Dati<-read_excel("Example for Mark2.xlsx")
head(Dati)
tail(Dati)
plot(Dati)
is.ts(Dati)
tas