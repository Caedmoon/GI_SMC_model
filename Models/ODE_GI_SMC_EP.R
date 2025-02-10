#Model notes
#Based on Corrias and Buist SMC model - biophysical model

#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)
      # Rate equations
      #Influx of Ca2+ through voltage gated, DHP sensitive L-type Ca2+ 
      I_CaL <- G_CaL * d * f * f_Ca * (y["V_m"] - E_Ca)
      #Ca2+ dependent inactivation variable
      f_Ca_inf <- 1 - (1/(1 + exp(-((y["Ca_i"] - h_Ca)/s_Ca))))
      #low voltage-activated, fast-inactivating DHP-insensitive component of the inward current
      I_LVA <- G_LVA * d_LVA * f_LVA * (y["V_m"] - E_Ca)
      #Delayed rectifier potassium channels
      I_Kr <- G_Kr * x_r1 * x_r2 * (y["V_m"] - E_k)
      #A type potassium channels
      I_KA <- G_KA * x_A1 * X_A2 - (y["V_m"] - E_k)
      #Large conductance Calcium-Activated Potassium channels
      I_BK <- G_BK * P_0 * (y["V_m"] - E_k)
      #Open probability which has been quantitatively described for the GI tract
      P_0 <- 1 / (1 + exp(y["V_m"]/K_BK - h_BK * log(Ca_i/Ca_set)))
      #TTX resistant sodium channel
      I_Na <- G_Na * d_Na * f_Na * (y["V_m"] - E_Na)
      #Non-Selective Cationic Channels
      I_NSCC <- G_NSCC * m_nscc * r_lig * h_Ca * (y["V_m"] - E_NSCC)
      #effect of intracellular Ca2+ has on I_NSCC
      h_Ca <- 1/(1 + (Ca_i/K_CaNSCC) * n_Ca)
      #background potassium conductance
      I_Kb <- G_Kb * (y["V_m"] - E_K)
      #Titak rate of Ca2+ uptake by the SR, mito and extrusion via PMCA
      I_Ca_Ext <- 0.317 * Ca_i^1.34
      
      
      # Rate equations
      I_ion <- I_CaL + I_LVA + I_Kr + I_Ka + I_BK + I_Kb + I_Na + I_NSCC
      
      if(t < t_ICC_peak){
        I_stim <- G_couple * delta_V_ICC
      }
      if(t_ICC_peak < t < t_ICC_plateau){
        I_stim <- G_couple * delta_V_ICC * 1/(1 + exp(t-8000/1000))
      }
      if(t_ICC_plateau < t < t_ICC){
        I_stim_t_ICC_plateau
        I_stim <- I_stim_t_ICC_plateau * 1.3/(1 + exp(t-8000/150))
      }
      #Activation
      
        
      #ODE
      #Cell membrane
      d[1] <- 1/C_m * (I_ion + I_stim)
      #CaI change in time
      d[2] <- (I_CaL + I_CaT) / (2 * Faraday * Vc) - I_Ca_Ext
      #gating variable
      d[3] <- (g_inf - y["g"]) / tau_g
      # Return the rates of change
      return(list(d))
    })
  }
  
  # Initial values for A, B, and C
  A_0 = 10
  B_0 = 0
  C_0 = 0
  initial <- c(A = A_0, B = B_0, C = C_0)
  
  #Fixed parameters
  fixed <- list(
    
  )
  
  # Time sequence for the simulation
  times <- seq(0, 100, by = 1)
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  
  return(output_df)
}

#Parameters to fit -----
parms <- list(
  k1 = 0.15,
  k2 = 0.02
)



#initial estimate -----
Initial_out <- Model(parms = parms)


test <- ggplot() +
  geom_line(data = Initial_out, aes(x = time, y = A, colour = "A")) +
  geom_line(data = Initial_out, aes(x = time, y = B, colour = "B")) +
  geom_line(data = Initial_out, aes(x = time, y = C, colour = "C")) +
  geom_point(data = obs_A, aes(x = time, y = A, colour = "A")) + 
  geom_point(data = obs_B, aes(x = time, y = B, colour = "B")) + 
  geom_point(data = obs_C, aes(x = time, y = C, colour = "C")) + 
  scale_colour_manual(values = c("A" = "blue", "B" = "red", "C" = "green")) +
  xlab("time") +
  ylab("amount") +
  ggtitle("Plot of chemicals")

print(test)
