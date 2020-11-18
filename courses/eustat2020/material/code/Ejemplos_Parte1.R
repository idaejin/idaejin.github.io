## ---- message=FALSE, warning=FALSE,error=FALSE--------------------------------
rm(list=ls())
library(lubridate) # Librería para fechas
library(tidyverse)
library(dplyr)
library(deSolve)


## -----------------------------------------------------------------------------
library(coronavirus)
?coronavirus
data(coronavirus)
head(coronavirus)


## -----------------------------------------------------------------------------
range(coronavirus$date)


## ---- echo=TRUE, message=FALSE, warning=FALSE---------------------------------
# Get top confirmed cases by state
coronavirus %>%
  filter(type == "confirmed") %>%
  group_by(country) %>%
  summarise(total = sum(cases)) %>%
  arrange(-total) %>%
  head(20)

# Get the number of recovered cases in Australia by state
coronavirus %>%
  filter(type == "recovered", country == "Australia") %>%
  group_by(province) %>%
  summarise(total = sum(cases)) %>%
  arrange(-total)


## ---- message = FALSE---------------------------------------------------------
`%>%` <- magrittr::`%>%`

datos <- coronavirus %>%
  dplyr::filter(country == "Belgium") %>%
  dplyr::group_by(date, type) %>%
  dplyr::summarise(total = sum(cases, na.rm = TRUE)) %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = total
  ) %>%
  dplyr::arrange(date) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(active = confirmed - death - recovered) %>%
  dplyr::mutate(
    confirmed_cum = cumsum(confirmed),
    death_cum = cumsum(death),
    recovered_cum = cumsum(recovered),
    active_cum = cumsum(active)
  )


## -----------------------------------------------------------------------------
sir_start_date <- "2020-02-04"
  sir_end_date <- "2020-03-30"


## ----fig.width=12,fig.height=15-----------------------------------------------
par(mfrow=c(2,1))
with(datos,{plot(date, confirmed, t = "h")})
abline(v = c(ymd(sir_start_date),ymd(sir_end_date)),col= 2)
with(datos,{plot(date, active_cum, t = "l")})
abline(v = c(ymd(sir_start_date),ymd(sir_end_date)),col= 2)


## -----------------------------------------------------------------------------
Infected <- subset(datos, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$active_cum
Day <- 1:length(Infected)


## -----------------------------------------------------------------------------
N <- 11515793 # Población en Belgica


## -----------------------------------------------------------------------------
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


## -----------------------------------------------------------------------------
init <- c(S = N - Infected[1],I = Infected[1],R= 0)

RSS <- function(parameters){
  names(parameters) <- c("beta","gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[,3]
  sum((Infected-fit)^2)
}

params_init <- c(0.5,0.5)
opt <- optim(params_init, RSS, method = "L-BFGS-B", lower = c(0,0), upper = c(1,1))
opt$message
opt$par <- setNames(opt$par, c("beta","gamma"))
opt$par

fit <- data.frame(ode(y = init, times = Day, func = SIR, parms = opt$par))
head(fit)


Day <- 1:(length(Infected))
Dates <- seq(as.Date(sir_start_date),as.Date(sir_end_date),by = 1)

plot(Dates, Infected, type ="p",cex=.65)
lines(Dates,fit$I,col=2)


## ----fig.width=10, fig.height=12----------------------------------------------
h = 30 # Horizonte de predicción a h días
Infected_long <- subset(datos, date >= ymd(sir_start_date) & date <= ymd(sir_end_date) + h)$active_cum
Day_long <- 1:length(Infected_long)
Dates_long <- seq(as.Date(sir_start_date),as.Date(sir_end_date)+h,by=1)

fit_long <- data.frame(ode(y = init, times = Day_long, func = SIR, parms = opt$par))
head(fit_long)

par(mfrow=c(2,1))
plot(Dates_long,Infected_long,xlab = "",cex=.5,pch=19)
lines(Dates_long,fit_long$I,col = 4,lwd=2)
lines(Dates,fit$I,col=3, lty = 2,lwd=2)
abline(v = ymd(sir_end_date), col= 2)
legend("topleft", c("Modelo SIR (ajuste)", "Modelo SIR (predicción)"), col = c("green","blue"),lty = 1:2,lwd=3, cex = .85)
plot(Dates_long,Infected_long,xlab = "",cex=.5,pch=19,log = "y")
abline(v = ymd(sir_end_date), col= 2)


## ---- echo = FALSE, eval = TRUE, results='hide', message=FALSE, warning=FALSE----
knitr::purl("Ejemplos_Parte1.Rmd", "Ejemplos_Parte1.R")

