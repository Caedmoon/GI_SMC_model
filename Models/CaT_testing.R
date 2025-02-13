#Model notes------
#Based on Poh et al (2012)

#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
T_type <- read.csv("data/T_type_IV.csv")
T_type$x <- round(T_type$x,0)

#normalisation
norm_stat_val <- function(b){
  sum_b <- sum(b, na.rm = TRUE)
  if(sum_b == 0) return(b)
  b / sum_b
}

# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)
      
      # Rate equations
      # Nernst Potentials----
      E_Ca <- ((R * Temp) / (2 * Faraday)) * log(Ca_o / Ca_i_total)
      #Single slow wave-----
      if (times >= 0 & times < t_peak_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (times/f2)
      } 
      if (times >= t_peak_ICC & times < t_plateau_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (1 + exp(-f1/(2 * t_slope))) * (1/(1 + exp((times - f2 - 0.5 * f1)/t_slope)))
      } 
      
      if (times < clamp_start | times > clamp_end ){
        Vm_eq <- Hv
      } else{
        Vm_eq <- Vm
      } 
      # Temp-type Calcium Current (ICaT)----
      d_CaT_inf <- d_CaT_inf(Vm_eq)
      f_CaT_inf <- f_CaT_inf(Vm_eq)
      
      tau_f_CaT <- 0.38117 * (8.6 + 14.7 * exp(-((Vm_eq +50)^2/900))) / T_correction_Ca
      #ODE-----
      # Temp-type Calcium Channel ODEs-----
      d[1] <- (d_CaT_inf - y["d_CaT"]) / tau_d_CaT
      d[2] <- (f_CaT_inf - y["f_CaT"]) / tau_f_CaT
      return(list(d))
    })
  }
  
  #Fixed parameters
  fixed <- list()
  R <- 8.314  # J/(molK) - Ideal gas constant
  Faraday <- 96.48534  # C/mmol - Faraday’s constant
  Temp <- 310  # K - Temperature
  Texp <- 297.0 #experimental temperatures
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
  T_correction_Ca <- 1.0 * CaQ10^((Temp - Texp) / 10.0) # experimental correction
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
  tau_d_CaT <- 1.9058 / T_correction_Ca #for ICaL
  sigma_LCaL <- 0# 0.01 #FOR ICaL set to 0 to replicate EGTA
  k_on <- 40633 #BK on rate
  k_c_off <- 11 #BK closed off rate
  k_o_off <- 1.1 #BK open off rate
  SMC_resting <- -60 #
  # Initial values for A, B, and C
  
  d_CaT_inf <- function(Vm) {
    return(1 / (1 + exp(-((Vm + 60.5) / 5.3))))
  }
  
  f_CaT_inf <- function(Vm) {
    return(1 / (1 + exp(((Vm + 75.5) / 4.0))))
  }
  
  initial <- c(
    d_CaT = 0.0202,
    f_CaT = 1.0
  )
  
  # Time sequence for the simulation
  times <- seq(0, 1000, by = 1)
  
  #history
  
  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$E_Ca <- ((R * Temp) / (2 * Faraday)) * log(Ca_o / Ca_i_total)
  output_df$I_CaT <- G_CaT * output_df$d_CaT * output_df$f_CaT * 
    (ifelse(output_df$time < parms[["clamp_start"]] | output_df$time > parms[["clamp_end"]], parms[["Hv"]], parms[["Vm"]]) - output_df$E_Ca)
  output_df$Vm <- ifelse(
    output_df$time < parms[["clamp_start"]] | output_df$time > parms[["clamp_end"]], 
    parms[["Hv"]],  # Use Hv before clamp_start or after clamp_end
    parms[["Vm"]]   # Use Vm between clamp_start and clamp_end
  )
  
  return(output_df)
}

#Parameters to fit -----
mv_tests <- seq(-90, 35, by = 5)#T_type$x#seq(-90, 35, by = 5)
sim_v <- list()
peak_store <- numeric(length(mv_tests))
i <- 1
for(i in seq_along(mv_tests)) {  # Correct looping over mv_tests
  parms <- list(Vm = mv_tests[i], Hv = -100, clamp_start = 500, clamp_end = 900)  # Set voltage parameter
  sim_temp <- Model(parms = parms)  # Run simulation
  sim_temp$Vm_identity <- parms[["Vm"]]
  sim_v[[i]] <- as.data.frame(sim_temp)
  sim_temp_peak <- subset(sim_temp, Vm == mv_tests[i])
  peak_store[i] <- min(sim_temp_peak$I_CaT)
}

x <- 1
for(x in 1:length(sim_v)) {  # Correct looping over mv_tests
  sim_v[[x]]$I_CaT_norm <- sim_v[[x]]$I_CaT / min(peak_store)
  #sim_v[[x]]$I_CaT_norm <- norm_stat_val(sim_v[[x]]$I_CaT)
}


sim_v_df <- bind_rows(sim_v)  # Convert list to dataframe

# Plot all conditions with color mapped to Vm
S6.plot <- ggplot(sim_v_df, aes(x = time, y = I_CaT_norm, color = factor(Vm_identity), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "Time vs Normalized Current",
    x = "Time (msec)",
    y = "Normalized Current",
    color = "Voltage (mV)"  # Legend label
  ) +
  scale_y_reverse()+
  xlim(480,580) +
  theme_minimal()  

# Display the plot
print(S6.plot)



# Create dataframe with results
IV.df <- data.frame(mV = mv_tests, I = peak_store)
IV.df$I_norm <- IV.df$I / min(IV.df$I)

IV.plot <- ggplot() +
  geom_point(data = T_type, aes(x = x, y = y), color = "red", size = 3) +  # Scatter points for experimental data
  geom_line(data = IV.df, aes(x = mV, y = I_norm), color = "blue", linewidth = 1) +  # Model I-V curve
  labs(
    title = "Normalized I-V Curve for I_CaT",
    x = "Membrane Voltage (mV)",
    y = "Normalized Peak I_CaT"
  ) +
  scale_y_reverse() +
  scale_x_continuous() +
  theme_minimal()  # Clean theme

# Display the plot
print(IV.plot)