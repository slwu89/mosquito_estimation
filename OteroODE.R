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
    alpha0 = 1.5,
    BS = 25,
    a0 = 0.75 # the critical value resulting from the product of the critical density times the estimated average volume of the breeding sites
  )
}

# t: temp in celsius
schoolfield <- function(t,R,RdK,Ha,Hh,T05){
  tk <- t + 273.15 # temp in kelvin
  RdK * (((tk/298) * exp((Ha/R)*((1/298) - (1/tk)))) / (1 + exp((Hh/R) * ((1/T05) - (1/tk)))))
}

# density-dependent egg hatching inhibition
gamma <- function(L,BS,a0){
  ifelse(L/BS < a0,0,0.63)
}

# smooth functions to replace the step-function
gamma2 <- function(L,BS,a0,a=0.63){
  (a*(L/BS)) / (a0 + (L/BS))
}

gamma3 <- function(L,BS,a0,a=0.63){
  (a*((L/BS)^2)) / ((a0^2) + ((L/BS)^2))
} 

# temperature dependent larval and pupal death
muL_t <- function(t){
  tk <- t + 273.15 # temp in kelvin
  0.01 + (0.9725*exp(-(tk-278)/2.7035))
}

muP_t <- function(t){
  tk <- t + 273.15 # temp in kelvin
  0.01 + (0.9725*exp(-(tk-278)/2.7035))
}

# plug in a time (in days) t and get a temperature in celsius
temp <- function(t,a=18,b=6.7,c=9.2){
  a + b*cos(((2*pi*t)/365.24) + c)
}

# alpha parameter
alpha <- function(alpha0,BS){
  alpha0/BS
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
    
    # cat("gammaL: ",gammaL,"\n")
    
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
    
    muL <- muL_t(t)
    muP <- muP_t(t)
    gammaL <- gamma(L0,BS,a0)
    
    # equilibria
    E0 <- L0 * ((beta*(ovr2 + muA)*ovr1*par*ef*lpr) / (2*muA*(muE + (elr*(1-gammaL)))*((muP*muA) + (muP*ovr1) + (par*muA) + (ovr1*par))))
    P0 <- L0 * (lpr / (muP + par))
    A10 <- L0 * ((par*ef*lpr) / (2*((muP*muA) + (muP*ovr1) + (par*muA) + (ovr1*par))))
    A20 <- L0 * ((ovr1*par*ef*lpr) / (2*muA*(muP*muA) + (muP*ovr1) + (par*muA) + (ovr1*par)))
    
    # check positivity constraints
    pos <- (elr * (1 - gammaL) * (E0 / L0)) - (lpr + muL)
    if(pos < 0){
      stop("positivity constraints not fulfilled")
    }
    
    return(
      list(
        E0=E0,
        P0=P0,
        A10=A10,
        A20=A20
      )
    )
  })
}

theta <- c(theta_const(),theta_temp())
theta$alpha <- alpha(alpha0 = theta$alpha0,BS = theta$BS)

L0 <- 1000
state_eq <- mod_eq(t = temp(1),L0 = L0,theta = theta)

mod_y <- with(state_eq,{
  c(E=E0,L=L0,P=P0,A1=A10,A2=A20)
})

mod_out <- lsoda(y = mod_y,times = 1:(365.25*2),func = mod_dx,parms = theta,
                 verbose = F,rtol = 1e-2)
plot(mod_out)