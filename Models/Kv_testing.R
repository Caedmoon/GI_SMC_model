#Model notes------
#Based on Poh et al (2012)
#Testing summaritive performance then investigate the individual ion current predictions using Fig 1 2 3 etc - 
#Create additional scripts to solve for each voltage supplied 
#maybe add all together and validate individual bits? 

#Fixes, set values to the x of the observed data



#addd in calcium buffering as it gives Ca2+ ree ? may be important?
#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
Kv_IV <- read.csv("GI SMC EP/data/Kv_IV.csv")
# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)
      
      # Rate equations
      # Nernst Potentials----
      E_K <- ((R * Temp) / (2 * Faraday)) * log(K_o / K_i)
      #Single slow wave-----
      if (times >= 0 & times < t_peak_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (times/f2)
      } 
      if (times >= t_peak_ICC & times < t_plateau_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (1 + exp(-f1/(2 * t_slope))) * (1/(1 + exp((times - f2 - 0.5 * f1)/t_slope)))
      } 
      # Voltage-dependent Potassium Current (IKv)----
      #to replicate patch clamp simulation
      if (times < holding_period){
        Vm_eq <- Hv
      } 
      if (times >= holding_period){
        Vm_eq <- Vm
      } 
      
      I_Kv <- G_Kv * y["x_Kv"] * y["y_Kv"] * (Vm_eq - E_K)
      
      x_Kv_inf <- x_Kv_inf(Vm_eq)
      y_Kv_inf <- y_Kv_inf(Vm_eq)
      #ODE-----
      # Voltage-dependent Potassium Channel ODEs
      d[1] <- (x_Kv_inf - y["x_Kv"]) / tau_x_Kv
      d[2] <- (y_Kv_inf - y["y_Kv"]) / tau_y_Kv
      return(list(d))
    })
  }
  
  #Fixed parameters
  fixed <- list()
  R <- 8.314  # J/(molK) - Ideal gas constant
  Faraday <- 96.48534  # C/mmol - Faraday’s constant
  Temp <- 310  # K - Temperature
  Cm <- 50  # pF - Cell membrane Capacitance
  Vcell <- 3.5e-12  # l - Cell volume
  Ca_o <- 2  # mM - Extracellular Calcium concentration
  K_o <- 5.4  # mM - Extracellular potassium concentration
  Na_o <- 140  # mM - Extracellular sodium concentration
  Ca_i_total <- 0.004914  # mM - Initial total intracellular Calcium concentration
  Ca_i_free <- 1.26e-4  # mM - Initial free intracellular Calcium concentration
  K_i <- 150  # mM - Intracellular potassium concentration
  Na_i <- 10.5  # mM - Intracellular sodium concentration
  CaQ10 <- 2.1  # - Q10 for Calcium channels
  KQ10 <- 3.1  # - Q10 for potassium channels
  NaQ10 <- 2.45  # - Q10 for sodium channels
  Gcouple <- 2.6  # nS - Coupling conductance between ICC and SMC
  V_rest_ICC <- -57  # mV - Resting membrane potential of ICC
  V_peak_ICC <- -23.5  # mV - Peak membrane potential of ICC
  V_amp_ICC <- 33.5  # mV - Amplitude of ICC membrane potential
  t_period <- 10000  # ms - Period of single ICC slow wave
  t_peak_ICC <- 300  # ms - Duration of ICC slow wave upstroke
  t_plateau_ICC <- 9700  # ms - Duration of ICC slow wave plateau phase
  t_slope <- 600  # ms - Slope factor in ICC membrane equation
  f1 <- 12000  # ms - ICC conditioning factor 1
  f2 <- 300  # ms - ICC conditioning factor 2
  CRT_total <- 0.034  # mM - Total Calreticulin concentration
  CRT_n <- 1  # - Hill coefficient for Calreticulin
  CRT_KD <- 0.0009  # mM - Dissociation constant for Calreticulin
  CaM_total <- 0.012  # mM - Total Calmodulin concentration
  CaM_n <- 4  # - Hill coefficient for Calmodulin
  CaM_KD <- 0.0001  # mM - Dissociation constant for Calmodulin
  G_CaL <- 1.44  # nS - Maximum conductance of CaL
  G_CaT <- 0.0425  # nS - Maximum conductance of CaT
  G_Kv <- 1.0217  # nS - Maximum conductance of Kv
  tau_x_Kv <- 4.7803  # ms - Time constant for Kvx
  tau_y_Kv <- 763.7564  # ms - Time constant for Kvy
  G_BK <- 80  # nS - Maximum conductance of BK
  G_Na <- 25.1  # nS - Maximum conductance of Na
  P_NCX <- 39.8437 # pA/pF - Maximum NCX
  K_mCa <- 1.38  # mM - Ca_i half saturation constant of NCX
  K_mNa <- 87.5  # mM - Na_i half saturation constant of NCX
  k_sat <- 0.1 # - Saturation factor for NCX
  gamma <- 0.35  # - Voltage dependence parameter of NCX
  P_NaK <- 0.1852  # pA/pF - Maximum NaK
  K_mK <- 1  # mM - K_o half saturation constant of NaK
  K_mNa <- 40  # mM - Na_i half saturation constant of NaK
  G_NSNa <- 0.022488  # nS - Maximum conductance of non-selective Na current
  G_NSK <- 0.017512  # nS - Maximum conductance of non-selective K current
  tau_d_CaT <- 1.9058 #ms
  sigma_LCaL <- 0.01 #FOR ICaL
  k_on <- 40633 #BK on rate
  k_c_off <- 11 #BK closed off rate
  k_o_off <- 1.1 #BK open off rate
  # Initial values for A, B, and C
  
  x_Kv_inf <- function(Vm) {
    return(1/ (1 + exp(-((Vm + 43)/17.36))))
  }
  
  y_Kv_inf <- function(Vm) {
    return(1/ (1 + exp(((Vm + 44.9)/12.0096))))
  }
  
  
  initial <- c(
    x_Kv = x_Kv_inf(parms[["Hv"]]),
    y_Kv = y_Kv_inf(parms[["Hv"]])
    )
  
  # Time sequence for the simulation
  times <- seq(750, 1000, by = 1) #miliseconds 
  
  #events
  #make it so when time is 820, voltage changes to the Vm

  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$E_K <- ((R * Temp) / (2 * Faraday)) * log(K_o / K_i)
  output_df$I_K <- G_Kv * output_df$x_Kv * output_df$y_Kv * 
    (ifelse(output_df$time < parms[["holding_period"]], parms[["Hv"]], parms[["Vm"]]) - output_df$E_K)
  # output_df$I_K_norm <- output_df$I_K / output_df$
  output_df$I_K_shifted <- output_df$I_K - output_df$I_K[1]
  return(output_df)
}





#Parameters to fit -----


mv_tests <- Kv_IV$x  # Voltage steps
peak_store <- numeric(length(mv_tests))  # Pre-allocate storage
sim_v <- list()
i <- 1
for(i in seq_along(mv_tests)) {  # Correct looping over mv_tests
  parms <- list(Vm = mv_tests[i], Hv = -70, holding_period = 820)  # Set voltage parameter
  sim_temp <- Model(parms = parms)  # Run simulation
  sim_temp$Vm <- parms[["Vm"]]
  sim_v[[i]] <- as.data.frame(sim_temp)
  peak_store[i] <- max(sim_temp$I_K_shifted[sim_temp$time > parms[["holding_period"]]])
  
}
peak_min <- min(peak_store, na.rm = TRUE)
peak_max <- max(peak_store, na.rm = TRUE)

x <- 1
for(x in 1:length(sim_v)) {  # Correct looping over mv_tests
   sim_v[[x]]$I_K_norm <- sim_v[[x]]$I_K_shifted / max(peak_store)
}

# Combine all simulations into one dataframe while keeping track of Vm
sim_v_df <- bind_rows(sim_v)  # Convert list to dataframe

# Plot all conditions with color mapped to Vm
S6.plot <- ggplot(sim_v_df, aes(x = time, y = I_K_norm, color = factor(Vm), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "Time vs Normalized Current",
    x = "Time (msec)",
    y = "Normalized Current",
    color = "Voltage (mV)"  # Legend label
  ) +
  xlim(breaks = c(750,1000)) +
  scale_y_reverse(breaks = seq(0,1, by = 0.1))+
  theme_minimal()  

# Display the plot
print(S6.plot)



#IV plotting
IV.df <- data.frame(mV = mv_tests, I = peak_store)

#IV.df <- IV.df[IV.df$mV >= -70,]

#IV.df <- IV.df[IV.df$mV <= 20,]



IV.df$I_norm <- (IV.df$I - min(IV.df$I)) / (max(IV.df$I) - min(IV.df$I))


#IV.df$I_norm[IV.df$I_norm < 0] <- 0

IV.plot <- ggplot() +
  geom_point(data = Kv_IV, aes(x = x, y = y), color = "red", size = 3) +  # Scatter points for experimental data
  geom_line(data = IV.df, aes(x = mV, y = I_norm), color = "blue", linewidth = 1) +  # Model I-V curve
  labs(
    title = "Normalized I-V Curve for K",
    x = "Membrane Voltage (mV)",
    y = "Normalized Peak I-K"
  ) +
  scale_x_continuous(breaks = c(-70, -70, -30, -10, 10,20)) +
  theme_minimal()  # Clean theme

# Display the plot
print(IV.plot)