#################################################################################################
# Modelos mecanicistas y estadísticos para brotes epidémicos: el caso de la COVID-19 en Euskadi 
# Nazioarteko Mintegia / Seminario Internacional Eustat 2020
# Vitoria-Gasteiz. 19-20 Noviembre 2020
#################################################################################################
# Tutorial en R, parte 1
#################################################################################################
# Dae-Jin Lee (dlee@bcamath.org)
# 
# Material disponible en:
# web: https://idaejin.github.io/courses/eustat2020/material/index.html

# R
# web: https://www.r-project.org

# Rstudio
# web: https://www.rstudio.com


# Instalar librerías de R ----------------------------------------------------------------------
# install.packages(c("deSolve","bblme","tidyverse","ggplot2","dplyr","lubridate","shiny","shinySIR","outbreaks","coronavirus"))

# Recursos web sobre COVID-19

# R epidemics consortium (RECON)
# Enlace: https://www.repidemicsconsortium.org/

# R Views is the RStudio blog devoted to the R Community and the R Language.
# Enlace https://rviews.rstudio.com/tags/covid-19/

library("deSolve")
?ode

## Paso 1: escribir la ODE en R
sir_equations <- function(time, variables, parameters){
  with(as.list(c(variables, parameters)), { 
    dS <- -beta * S * I/N
    dI <-  beta * S * I/N - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))})
}

## Paso 2: definir unos valores para los parámetros
# N población total
#  beta: tasa de contagio
# gamma: tasa de recuperación
N <- 1000
parameters_values <- c(
   beta = 0.4,
  gamma = .05
)


## Paso 3: definir unos valores iniciales para las variables
initial_values <- c(
  S = N-1,  # Nº de susceptibles en t = 0
  I =   1,  # Nº de infecciosos en t = 0
  R =   0   # Nº de recuperados (en inmunes) en t = 0
)


## Paso 4: crear una secuencia de puntos
time_values <- 1:100 # Días

parameters_values
initial_values
time_values


## Paso 5: resolver la ecuación diferencial con ode()
sir_values_1 <- ode(
      y = initial_values,
  times = time_values,
   func = sir_equations,
  parms = parameters_values 
)

head(sir_values_1)


## Veamos el resultado gráficamente
matplot(sir_values_1[,1],sir_values_1[,2:4], col = c("blue","red","green"),
        type="l",lty=1,lwd=3, xlab = "Tiempo (días)", ylab = "Nº de personas")
legend("right", c("Susceptibles", "Infecciosos", "Recuperados"),
       col = c("blue", "red", "green"), lty = 1, lwd = 3, bty = "n")
abline(h = N, lty = 2)


## R0 (Número reproductivo básico)
R0 <- setNames(parameters_values["beta"] / parameters_values["gamma"],"R0")
R0 


## Creando una función propia SIR en R
mySIR <- function(beta, gamma, S0, I0, R0, times) {
 
  # the differential equations:
  sir_equations <- function(time, variables, parameters){
    with(as.list(c(variables, parameters)), { 
      dS <- -beta * S * I / N
      dI <- beta * S * I / N - gamma * I
      dR <- gamma * I
      return(list(c(dS, dI, dR))) })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}

mySIR(beta = 0.4, gamma=0.05, S0=N-1, I0=1, R0=0, times = seq(0,100))


## Ajuste de un modelo SIR con diferentes valores de beta y gamma = 0.05
betas <- c(0.2,0.4,0.75,1.25)
gammas <- rep(0.05,4)
mod <- list()

par(mfrow=c(2,2))
for(i in 1:4){
  mod[[i]]<-mySIR(beta = betas[i], gamma = gammas[i], S0 = 999, I0 = 1, R0 = 0, times = seq(0, 100, l = 200))
  m = betas[i]
  with(mod[[i]],{
    plot(time,S, col = "blue", type = "l", xlab = "Días", ylab = "Nº de personas", lwd = 3, 
         main = bquote(~ beta == .(betas[i]) ~ "," ~gamma == .(gammas[i])))
    lines(time,I, col = "green", lwd = 3)
    lines(time,R, col = "red", lwd = 3)}
  )
  legend("right", c("Susceptibles", "Infecciosos", "Recuperados"),
         col = c("blue", "red", "green"), lty = 1, lwd = 3, bty = "n",cex=.65)
  abline(h = N, lty = 2)
}
dev.off()

## Librería shiny
# install.packages("shinySIR")

library(shinySIR)
run_shiny(model = "SIR")


# Modelo SIR y ajuste de datos
## Ejemplo:
library(outbreaks)
data(influenza_england_1978_school)
head(influenza_england_1978_school)


## Gráficamente
datos <- influenza_england_1978_school
plot(datos[,c(1,2)], xlab = "Día", 
     ylab = "Nº de alumnos en cama", type = "h", col = "lightgrey")
points(datos[,c(1,2)],pch=19)


 dias <- seq(1,14)
casos <- datos$in_bed

N <- 763
predictions <- mySIR(beta = 4, gamma = 0.5, 
                     S0 = N-1, I0 = 1, R0 = 0, times = dias)

head(predictions)


## Gráficamente:
matplot(dias,predictions[,2:4],t="l",lwd=3, ylab ="")


## Gráficamente:
plot(datos[,c(1,2)], xlab = "Día", ylab = "Nº de alumnos en cama", type = "h", 
     col = "lightgrey",ylim = c(0,600))
points(datos[,c(1,2)],pch=19)
lines(datos$date,predictions$I, col = "red", lwd=2)


## Suma de Cuadraro de los Residuos:
SCR <- function(beta, gamma, data = datos, N = 763) {
  I0 <- datos$casos[1] # Valores iniciales de I
  times <- datos$dias
  predictions <- mySIR(beta = beta, gamma = gamma,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0, # variables' intial values
                       times = times)                # time points
  out <- sum((predictions$I[-1] - datos$casos[-1])^2)
  return(out)
}


## Renombro los datos:
datos <- influenza_england_1978_school
datos$dias <- seq(1,14) # vector de días (numérico)
colnames(datos)[2] <- "casos" # renombro la columna in_bed -> casos

head(datos)


## Cálculo SCR:
res <- SCR(beta = 4, gamma = 0.5, data = datos, N = 763)
res


## Secuencia de valores posibles de beta y gamma:
 beta_val <- seq(from = 1, to = 3, l = 100)
gamma_val <- seq(from = 0.4, to = 0.5, l = 100)

## suma de cuadrados para los valores de beta_val con gamma = 0.5
ss_val_beta <- sapply(beta_val, SCR, gamma = 0.5)
min(ss_val_beta) # Valor mínimo de la SCR 


## Valor de beta_hat tal que la suma de cuadrados es mínima
beta_hat <- beta_val[which.min(ss_val_beta)]
beta_hat


## Gráficamente:
plot(beta_val, ss_val_beta, 
     xlab = expression(paste(beta)), 
     ylab = "SCR", lwd=3, type="l")
abline(v = beta_hat,h = min(ss_val_beta), lty=2, lwd = 3, col = "grey")
points(beta_hat,min(ss_val_beta),col="red",pch=19)


## Calculamos el valor de gamma que minimiza la scr:
ss_val_gamma <- sapply(gamma_val, function(x) SCR(beta_hat, x))
min(ss_val_gamma)
gamma_hat <- gamma_val[which.min(ss_val_gamma)]
gamma_hat 


## Gráficamente:
plot(gamma_val, ss_val_gamma, 
     xlab = expression(paste(gamma)), 
     ylab = "SCR", lwd=3, type="l")
abline(v = gamma_hat,h = min(ss_val_gamma), lty=2, lwd=3, col = "grey")
points(gamma_hat,min(ss_val_gamma),col="red",pch=19)


## Vamos a buscar de manera simultánea los valores de beta y gamma
        n <- 50
 beta_val <- seq(from = 0.8, to = 3, l = n)
gamma_val <- seq(from = 0.25, to = 0.8, l = n)
param_val <- expand.grid(beta=beta_val,gamma=gamma_val)
ss_val <- with(param_val, Map(SCR, beta, gamma))
ss_val <- matrix(unlist(ss_val), n)

min(ss_val)

param_val[which.min(ss_val),]


## Perspective plot
persp(beta_val, gamma_val, ss_val, 
      theta = 40, phi = 20, 
      xlab = expression(beta), ylab = expression(gamma),
      zlab = "SCR",
      main = paste0("SCR = ",round(min(ss_val),2)))


## Image plot
image(beta_val, gamma_val, ss_val, xlab = expression(beta),col= rev(hcl.colors(12, "Teal")),
      ylab = expression(gamma),main = paste0("SCR = ",round(min(ss_val),2)))
contour(beta_val, gamma_val, ss_val,add=TRUE,col="black",lwd=2,nlevels = 10)
points(beta_hat, gamma_hat, pch = 3, lwd = 4, col="red")
points(param_val[which.min(ss_val),],pch=3,col="blue")

## Modelo SIR: datos gripe en el colegio inglés
 N <- 763              # Población total
I0 <- datos$casos[1]  # Infecciosos (valor inicial)
S0 <-  N-I0           # Susceptibles (valor inicial)      

time_points <- seq(1,14,l=100)  # 100 valores entre 1 y 14 días

# Ajuste del modelo SIR
fit <- mySIR(beta = beta_hat, gamma = gamma_hat, 
             S0 = S0, I0 = I0, R0 = 0, times = time_points)


## Gráficamente:
# vector de fechas (para el gráfico)
dates_points <- seq(min(datos$date),max(datos$date),l=100)

plot(datos$date,datos$casos,col="black",pch=19, 
     xlab = "Días", ylab = "Nº de casos")
lines(dates_points,fit$I,col="red", lwd = 2)


## Gráficamente (representamos la dinámica del SIR)
plot(dates_points,fit$S,t="l",col=c("blue"),lwd=3,
     xlab = "Días", ylab = "Nº Individuos")
lines(dates_points,fit$I,col="red",lwd=3)
lines(dates_points,fit$R,col="green",lwd=3)
points(datos$date,datos$casos,pch=19,col="black")


## Modelo SIR: con optim()
SCR2<- function(params){
  out <- SCR(beta = params[1], gamma = params[2], data=datos, N=763)
  return(out)
}

params <- c(1, 0.5)
SCR2(params)


## Valores iniciales y la función optim()
init_params <- c(1,0.45)
ss_opt <- optim(init_params, SCR2)
ss_opt


## Extraemos los parámetros de optim
 beta_hat2 <- ss_opt$par[1]
gamma_hat2 <- ss_opt$par[2]

# Ajuste del modelo SIR
fit2 <- mySIR(beta = beta_hat2, gamma = gamma_hat2, 
             S0 = S0, I0 = I0, R0 = 0, times = time_points)

# Comparamos los ajustes con el método de búsqueda y con optim()
plot(dates_points,fit$I,t="l",col="red",lwd=3, xlab = "Días", ylab = "Nº Individuos")
lines(dates_points,fit2$I,col="blue",lty=2,lwd=3)
points(datos$date,datos$casos,pch=19,col="black")
legend("topright", c("búsqueda parámetros", "optim()"), 
       lwd=3, lty=1:3, col=c("red","blue"))


## Modelo SIR: estimación por Máxima Verosimilitud
library(bbmle)

## Definimos una función de la log-verosimilud 
mLL <- function(beta, gamma, sigma, dias, casos, N = 763) {
   beta <- exp(beta)    # para asegurar parámetros >0 
  gamma <- exp(gamma) 
  sigma <- exp(sigma)
  I0 <- casos[1] 
  observations <- casos[-1] # ajustamos los datos restantes
  predictions <- mySIR(beta = beta, gamma = gamma,
                       S0 = N - I0, I0 = I0, R0 = 0, times = dias)
  predictions <- predictions$I[-1] # eliminamos el primer valor ajustado
# - log-verosimilitud de una Distribución Normal:
  -sum(dnorm(x = observations, mean = predictions, sd = sigma, log = TRUE))
}

## Lista de valores iniciales de los parámetros
starting_param_val <- list(beta = 1, gamma = 0.30, sigma = 1)

# Estimamos con mle2
estimates <- mle2(minuslogl = mLL, start = lapply(starting_param_val, log),
                  method = "Nelder-Mead", data = c(datos, N = 763))
estimates


## Resumen
summary(estimates)


## Transformamos los coeficientes a la escala original
exp(coef(estimates))


## Calculamos con confint() los IC95%
exp(confint(estimates))


## Obtenemos la matriz de var-cov
vcov(estimates)


## menos log-verosimilitud
-logLik(estimates)


## Akaike Information Criteria (AIC) y AIC condicional
AIC(estimates)
AICc(estimates)


## Bayesian Information Criteria (BIC)
BIC(estimates)


## Vamos a ajustar el modelo SIR con los parámetros estimados
 beta_hat_mle <- exp(coef(estimates))[1]
gamma_hat_mle <- exp(coef(estimates))[2]
sigma_hat_lme <- exp(coef(estimates))[3]

# Ajustamos el modelo SIR con los parámetros obtenidos por MV
fit_mle <- mySIR(beta = beta_hat_mle,
                 gamma = gamma_hat_mle, 
                 S0 = S0, I0 = I0, R0 = 0, times = time_points)


## Cálculo de los límites superior e inferior
low <- qnorm(p = 1-0.975,mean=fit_mle$I,sd = sigma_hat_lme)
 up <- qnorm(p = 0.975,mean=fit_mle$I,sd = sigma_hat_lme)


## Gráfico del ajuste e IC95%
plot(datos$dias,datos$casos,col="red",pch=19, ylim = c(0,max(up)),
     xlab = "Días", ylab = "Nº de casos")
polygon(c(time_points, rev(time_points)), c(up, rev(low)),
        border = NA, col = adjustcolor("red", alpha.f = 0.1))
lines(time_points,fit_mle$I,col="red")


## Modelo SIR: estimación por MV caso Poisson
# Definimos una función de la log-verosimilud 
mLL_pois <- function(beta, gamma, dias, casos, N = 763) {
   beta <- exp(beta)    # para asegurar parametros >0 
  gamma <- exp(gamma)   #  parametros en escala log

  I0 <- casos[1] # 
  observations <- casos[-1] # ajustamos los datos restantes
  predictions <- mySIR(beta = beta, gamma = gamma,
                       S0 = N - I0, I0 = I0, R0 = 0, times = dias)
  predictions <- predictions$I[-1] # eliminamos el primer valor ajustado
  if (any(predictions < 0)) return(NA) # En caso de valores negativos =NA
# - log-verosimilitud de una Distribución Poisson:
  -sum(dpois(x = observations, lambda = predictions, log = TRUE))
}


## Lista de valores iniciales de los parámetros
starting_param_val <- list(beta = 0.004, gamma = 0.5)

## Estimamos con mle2
estimates_pois <- mle2(minuslogl = mLL_pois, 
                       start = lapply(starting_param_val, log),
                       method = "Nelder-Mead", data = c(datos, N = 763))
estimates_pois


## Resumen y coeficientes
summary(estimates_pois)
exp(coef(estimates_pois)) # Distr. Poisson


## coeficientes en el caso Normal
exp(coef(estimates))      # Distr. Normal

## IC95%
exp(confint(estimates))      # Distr. Normal
exp(confint(estimates_pois)) # Distr. Poisson


## Coeficientes beta y gamma del modelo Poisson
 beta_hat_pois <- exp(coef(estimates_pois))[1]
gamma_hat_pois <- exp(coef(estimates_pois))[2]

# Ajustamos el modelo SIR con los parámetros obtenidos por MV
fit_pois <- mySIR(beta = beta_hat_pois,
                  gamma = gamma_hat_pois, 
                  S0 = S0, I0 = I0, R0 = 0, times = time_points)
head(fit_pois)


## Cálculo de los límites superior e inferior
low_pois <- qpois(p = 1-0.975,lambda=fit_pois$I)
 up_pois <- qpois(p = 0.975,lambda=fit_pois$I)


## Estimación e IC95%
plot(datos$dias,datos$casos,col="red",pch=19, ylim = c(0,max(up_pois)),
     xlab = "Días", ylab = "Nº de casos")
polygon(c(time_points, rev(time_points)), c(up_pois, rev(low_pois)),
        border = NA, col = adjustcolor("green", alpha.f = 0.3))
lines(time_points,fit_pois$I,col="darkgreen")


## Modelo SIR: estimación por MV (caso Bimomial Negativa)
mLL_nb <- function(beta, gamma, theta, dias, casos, N = 763) {
  beta <- exp(beta)    # para asegurar parametros >0 
gamma <- exp(gamma)   #  parametros en escala log
theta <- exp(theta)
 size <- 1/(theta)    # dispersion parameter
  I0 <- casos[1] # 
  observations <- casos[-1] # ajustamos los datos restantes
  predictions <- mySIR(beta = beta, gamma = gamma,
                       S0 = N - I0, I0 = I0, R0 = 0, times = dias)
  predictions <- predictions$I[-1] # eliminamos el primer valor ajustado
  if (any(predictions < 0)) return(NA) # En caso de valores negativos =NA
  mu <- predictions
  prob <- size/(mu+size)
  # - log-verosimilitud de una Distribución Binomial Negativa
  -sum(dnbinom(x = observations, size = size, mu = mu, log = TRUE))
}

# Lista de valores iniciales de los parámetros
starting_param_val <- list(beta = 4, gamma = 0.5, theta = 1)

# Estimamos con mle2
estimates_nb <- mle2(minuslogl = mLL_nb, 
                     start = lapply(starting_param_val, log),
                     method = "Nelder-Mead", data = c(datos, N = 763))
estimates_nb


## Resumen
summary(estimates_nb)


## Coeficientes en la escala original
exp(coef(estimates_nb))


## Intervalos de confianza 95%
exp(confint(estimates_nb))


## Ajuste del Modelo SIR (Binomial Negativa)
 beta_hat_nb <- exp(coef(estimates_nb))[1]
gamma_hat_nb <- exp(coef(estimates_nb))[2]
theta_hat_nb <- exp(coef(estimates_nb))[3]

# Parámetro de dispersión 
1/theta_hat_nb

# Ajustamos el modelo SIR con los parámetros obtenidos por MV
fit_nb <- mySIR(beta = beta_hat_nb,
                gamma = gamma_hat_nb, 
                S0 = S0, I0 = I0, R0 = 0, times = time_points)


## ---- echo = TRUE------------------------------------------------------------------------------------------
# Calculamos los IC al 95%
 size = 1/theta_hat_nb
 mu = fit_nb$I
 prob = size/(mu+size)
 
 low_nb <- qnbinom(p = 1-0.975,size = size,prob = prob)
  up_nb <- qnbinom(p = 0.975,size = size,prob = prob)


## ---- echo = TRUE , fig.align='center', out.width="65%"----------------------------------------------------
plot(datos$dias,datos$casos,col="red",pch=19, ylim = c(0,max(up_nb)),
     xlab = "Días", ylab = "Nº de casos")
polygon(c(time_points, rev(time_points)), c(up_nb, rev(low_nb)),
        border = NA, col = adjustcolor("blue", alpha.f = 0.3))
lines(time_points,(fit_nb$I),col="blue")


## Vamos a comparar la menos log-verosimilitud y los AIC
-logLik(estimates)

-logLik(estimates_pois)

-logLik(estimates_nb)

AICtab(estimates,estimates_pois, estimates_nb,logLik = TRUE)


## Nº reproductivo básico
R0 <- beta_hat/gamma_hat
R0 
infectious_period <- 1/gamma_hat
infectious_period 

# visualización con shinySIR
run_shiny(model = "SIR", ics = c(S=762,I=1,R=0), tstart = 1, tmax=14,
          timestep = 0.01,
          parm0 = c(R0 = R0, Ip = infectious_period),
          parm_min = c(R0 = 0, Ip = 1),parm_max = c(R0 = 8, Ip = 10),
          parm_names = c("R0", "Periodo de Infección"),
          sigfigs = 3, slider_steps = c(0.001,0.001))

## Número reproductivo básico Rt
library(EpiEstim)

# Data on the 2009 H1N1 influenza pandemic in a school in Pennsylvania.
data("Flu2009")

names(Flu2009)

# method "non_parametric_si"
res <- estimate_R(incid = Flu2009$incidence, method = "non_parametric_si",
                  config = make_config(list(si_distr = Flu2009$si_distr)))
plot(res)

# method "parametric_si", con supuestos sobre el intervalo serial
res <- estimate_R(Flu2009$incidence, method = "parametric_si",
                  config = make_config(list(mean_si = 2.6, std_si = 1.5)))
plot(res)

# estimate the reproduction number (method "uncertain_si")
res <- estimate_R(Flu2009$incidence, method = "uncertain_si",
                  config = make_config(list(
                  mean_si = 2.6, std_mean_si = 1,
                  min_mean_si = 1, max_mean_si = 4.2,
                  std_si = 1.5, std_std_si = 0.5,
                  min_std_si = 0.5, max_std_si = 2.5,
                  n1 = 100, n2 = 100)))
plot(res)

