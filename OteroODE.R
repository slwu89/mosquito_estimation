
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

# t: temp in celsius
schoolfield <- function(t,R,RdK,Ha,Hh,T05){
  RdK * ((((t+273.15)/298) * exp((Ha/R)*((1/298) - (1/(t+273.15))))) / (1 + exp((Hh/R) * ((1/T05) - (1/(t+273.15))))))
}
