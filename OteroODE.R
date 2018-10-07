################################################################################
# 
#   Otero et al (2006) model
#   Kolmogorov Forward Equations (ODEs)
#   Sean Wu & Fausto Bustos
#   October 2018
# 
################################################################################

################################################################################
# load packages
################################################################################

rm(list=ls());gc()
library(compiler)
library(deSolve)
enableJIT(3)


################################################################################
# model functions
################################################################################

# temperature dependent constants
theta_temp <- function(){
  list(
    # universal gas constant
    R = 1.9858775,
    # rates at 298K
    RdK_e = 0.24,
    RdK_l = 0.2088,
    RdK_p = 0.384,
    RdK_a1 = 0.216,
    RdK_a2 = 0.372,
    # thermodynamic parameters
    Ha_e = 10798,
    Ha_l = 26018,
    Ha_p = 14931,
    Ha_a1 = 15725,
    Ha_a2 = 15725,
    Hh_e = 100000,
    Hh_l = 55990,
    Hh_p = -472379,
    Hh_a1 = 1756481,
    Hh_a2 = 1756481,
    # half-maximum
    T05_e = 14184,
    T05_l = 304.6,
    T05_p = 148,
    T05_a1 = 447.2,
    T05_a2 = 447.2
  )
}

# other biological parameters
theta_const <- function(){
  list(
    # eggs/oviposition
    beta = 63,
    muE = 0.01,
    ef = 0.83, # emergence fraction
    muA =0.09,
    a0 = 1.5,
    BS = 24
  )
}

# t: temp in celsius
schoolfield <- function(t,R,RdK,Ha,Hh,T05){
  RdK * ((((t+273.15)/298) * exp((Ha/R)*((1/298) - (1/(t+273.15))))) / (1 + exp((Hh/R) * ((1/T05) - (1/(t+273.15))))))
}

# density-dependent egg hatching inhibition
gamma <- function(L,BS,a0){
  ifelse(L/BS < a0,0,0.63)
}

# temperature dependent larval and pupal death
muL_t <- function(t){
  0.01 + (0.9725*exp(-(t-278)/2.7035))
}

muP_t <- function(t){
  0.01 + (0.9725*exp(-(t-278)/2.7035))
}

# plug in a time (in days) t and get a temperature in celsius
temp <- function(t,a=18,b=6.7,c=9.2){
  a + b*cos(((2*pi*t)/365.24) + c)
}

# alpha parameter
alpha <- function(a0,BS){
  a0/BS
}

# vector of dx/dt for all states x
mod_dx <- function(time,state,theta){
  with(as.list(c(state,theta)),{
    
    # calc temp-dependent things
    t <- temp(time)
    elr <- schoolfield(t,R,RdK_e,Ha_e,Hh_e,T05_e)
    lpr <- schoolfield(t,R,RdK_l,Ha_l,Hh_l,T05_l)
    par <- schoolfield(t,R,RdK_p,Ha_p,Hh_p,T05_p)
    ovr1 <- schoolfield(t,R,RdK_a1,Ha_a1,Hh_a1,T05_a1)
    ovr2 <- schoolfield(t,R,RdK_a2,Ha_a2,Hh_a2,T05_a2)
    
    # mortality & aquatic
    muL <- muL_t(t)
    muP <- muP_t(t)
    gammaL <- gamma(L,BS,a0)
    
    # dx/dt
    dE <- beta*((ovr1*A1) + (ovr2*A2)) - (muE*E) - (elr*(1-gammaL)*E)
    dL <- (elr*(1-gammaL)*E) - (muL*L) - (alpha*(L^2)) - (lpr*L)
    dP <- (lpr*L) - (muP*P) - (par*P)
    dA1 <- (par*ef*P/2) - (muA*A1) - (ovr1*A1)
    dA2 <- (ovr1*A1) - (muA*A2)
    
    return(list(c(dE,dL,dP,dA1,dA2)))
  })
}

# equilibrium pop sizes
mod_eq <- function(t,L0,theta){
  with(theta,{
    
    elr <- schoolfield(t,R,RdK_e,Ha_e,Hh_e,T05_e)
    lpr <- schoolfield(t,R,RdK_l,Ha_l,Hh_l,T05_l)
    par <- schoolfield(t,R,RdK_p,Ha_p,Hh_p,T05_p)
    ovr1 <- schoolfield(t,R,RdK_a1,Ha_a1,Hh_a1,T05_a1)
    ovr2 <- schoolfield(t,R,RdK_a2,Ha_a2,Hh_a2,T05_a2)
    
    
    E0 <- L0 * ( (beta*(ovr2 + muA)*ovr1*par*ef*lpr) / )
  })
}


theta <- c(theta_const(),theta_temp())
theta$alpha <- alpha(a0 = theta$a0,BS = theta$BS)
