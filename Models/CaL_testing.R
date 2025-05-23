#Model notes------
#Based on Poh et al (2012)

#Kv - add in the holding voltage stuff following the patch clamp protocol
#This requires patch cycling will need to include the cycling of it - 


#has to be problem within equations as
#Breaking down ICAL equation, 0 Vm - constant E_Ca  will always be larger than 10 - E_Ca
#difference in has to be the effect of voltage activation/inactivation of channels
#As Calcium effects are presumed to be 0
#Maybe get someone to check?
#move on and test other channels. 



#addd in calcium buffering as it gives Ca2+ ree ? may be important?
#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
L_type <- read.csv("data/L_type_IV.csv")
L_type$x <- round(L_type$x,0)

#normalisation(not for peak)----
norm_stat_val <- function(b){
  sum_b <- sum(b)
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
      
      # L-type Calcium Current (ICaL) - need to add markov model to solve for C and O as well as for IBK ----
      # Voltage-dependent opening/closing rates
      a_CaL <- T_correction_CaL * 0.7310 * exp(Vm_eq / 30) #When Vm is > 0, A increases
      b_CaL <- T_correction_CaL * 0.2149 * exp(-Vm_eq / 40) #when Vm is > 0, B decreases
      a0_LCaL <- 4 * a_CaL  # S-15
      a1_LCaL <- 3 * a_CaL  # S-16
      a2_LCaL <- 2 * a_CaL  # S-17
      a3_LCaL <- a_CaL  # S-18
      
      b0_LCaL <- b_CaL  # S-19
      b1_LCaL <- 2 * b_CaL  # S-20
      b2_LCaL <- 3 * b_CaL  # S-21
      b3_LCaL <- 4 * b_CaL  # S-22
      # Fast and slow inactivation rates (S-23, S-24)
      phi_f_LCaL <- 0.4742 * exp(Vm_eq / 10)   # Fast inactivation rate (S-23) - Increase
      phi_s_LCaL <- 0.005956 * exp(-Vm_eq / 40)  # Slow inactivation rate (S-24) - Decreases
      
      # Additional inactivation rate equations (S-25, S-26, S-27, S-28)
      xi_f_LCaL <- 0.01407 * exp(-Vm_eq / 300)  # Fast inactivation recovery rate (S-25)
      xi_s_LCaL <- 0.01213 * exp(Vm_eq / 500)   # Slow inactivation recovery rate (S-26)
      psi_f_LCaL <- 0.02197 * exp(Vm_eq / 500)  # Fast inactivation development rate (S-27)
      psi_s_LCaL <- 0.00232 * exp(-Vm_eq / 280) # Slow inactivation development rate (S-28)
      
      # Transition probabilities (S-29, S-30)
      omega_f_LCaL <- (b3_LCaL * xi_f_LCaL * phi_f_LCaL) / (a3_LCaL * psi_f_LCaL)  # Fast inactivation transition probability (S-29)
      omega_s_LCaL <- (b3_LCaL * xi_s_LCaL * phi_s_LCaL) / (a3_LCaL * psi_s_LCaL)  # Slow inactivation transition probability (S-30)
      
      # Transition between fast and slow inactivation (S-31, S-32)
      omega_sf_LCaL <- xi_s_LCaL * psi_f_LCaL / xi_f_LCaL  # Transition from fast to slow inactivation (S-31)
      omega_fs_LCaL <- psi_s_LCaL  # Transition from slow to fast inactivation (S-32)
      
      # Calcium-dependent inactivation rate
      
      if(CalciumDependence == 1){
        theta_LCaL <- 4/(1 + (1 / Ca_i_free)) 
      }else{
        theta_LCaL <- 0 #- this is set to 0 to replicate EGTA present
      }
      
      #ODE-----
      #L-type Calcium current markov map
      # Closed states without Calcium binding
      d[1] <- b0_LCaL * y["C1_I_LCaL"] + sigma_LCaL * y["C0Ca_I_LCaL"] - (theta_LCaL + a0_LCaL) * y["C0_I_LCaL"]
      d[2] <- a0_LCaL * y["C0_I_LCaL"] + b1_LCaL * y["C2_I_LCaL"] + sigma_LCaL * y["C1Ca_I_LCaL"] - (b0_LCaL + a1_LCaL + theta_LCaL) * y["C1_I_LCaL"]
      d[3] <- a1_LCaL * y["C1_I_LCaL"] + b2_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["C2Ca_I_LCaL"] - (b1_LCaL + a2_LCaL + theta_LCaL ) * y["C2_I_LCaL"]
      d[4] <- a2_LCaL * y["C2_I_LCaL"] + b3_LCaL * y["O_I_LCaL"] + sigma_LCaL * y["C3Ca_I_LCaL"] + omega_s_LCaL * y["IVS_I_LCaL"] + omega_f_LCaL * y["IVF_I_LCaL"] - (b2_LCaL + a3_LCaL + phi_s_LCaL + phi_f_LCaL + theta_LCaL) * y["C3_I_LCaL"]
      
      # Closed states with Calcium binding
      d[5] <- b0_LCaL * y["C1Ca_I_LCaL"] + theta_LCaL * y["C0_I_LCaL"] - (a0_LCaL + sigma_LCaL) * y["C0Ca_I_LCaL"]
      d[6] <- a0_LCaL * y["C0Ca_I_LCaL"] + b1_LCaL * y["C2Ca_I_LCaL"] + theta_LCaL * y["C1_I_LCaL"] - (b0_LCaL + a1_LCaL + sigma_LCaL) * y["C1Ca_I_LCaL"]
      d[7] <- a1_LCaL * y["C1Ca_I_LCaL"] + b2_LCaL * y["C3Ca_I_LCaL"] + theta_LCaL * y["C2_I_LCaL"] - (b1_LCaL + a2_LCaL + sigma_LCaL) * y["C2Ca_I_LCaL"]
      d[8] <- a2_LCaL * y["C2Ca_I_LCaL"] + omega_f_LCaL * y["IVF_Ca_I_LCaL"]  + theta_LCaL * y["C3_I_LCaL"] + omega_s_LCaL * y["IVS_Ca_I_LCaL"] + b3_LCaL * y["ICa_I_LCaL"] - (b2_LCaL + a3_LCaL + phi_f_LCaL + phi_s_LCaL + sigma_LCaL) * y["C3Ca_I_LCaL"]
      
      # Open state
      d[9] <- a3_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["ICa_I_LCaL"]  + xi_s_LCaL * y["IVS_I_LCaL"] + xi_f_LCaL * y["IVF_I_LCaL"] - (b3_LCaL + theta_LCaL + psi_s_LCaL + psi_f_LCaL) * y["O_I_LCaL"]
      
      # Inactivation states
      d[10] <- psi_f_LCaL * y["O_I_LCaL"] + phi_f_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["IVF_Ca_I_LCaL"] + omega_sf_LCaL * y["IVS_I_LCaL"] - (xi_f_LCaL + theta_LCaL + omega_fs_LCaL + omega_f_LCaL) * y["IVF_I_LCaL"] # Voltage dependent - fast 
      d[11] <- phi_s_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["IVS_Ca_I_LCaL"] + omega_fs_LCaL * y["IVF_I_LCaL"] + psi_s_LCaL * y["O_I_LCaL"] - (omega_s_LCaL + xi_s_LCaL + omega_sf_LCaL + theta_LCaL) * y["IVS_I_LCaL"] # Voltage-dependent - slow
      d[12] <- a3_LCaL * y["C3Ca_I_LCaL"] + theta_LCaL * y["O_I_LCaL"] + xi_f_LCaL * y["IVF_Ca_I_LCaL"] + xi_s_LCaL * y["IVS_Ca_I_LCaL"] - (b3_LCaL + sigma_LCaL + psi_f_LCaL + psi_s_LCaL) * y["ICa_I_LCaL"] # Calcium dependent
      d[13] <- phi_f_LCaL * y["C3Ca_I_LCaL"] + psi_f_LCaL * y["ICa_I_LCaL"] + theta_LCaL * y["IVF_I_LCaL"] + omega_sf_LCaL * y["IVS_Ca_I_LCaL"] - (sigma_LCaL + xi_f_LCaL + omega_f_LCaL + omega_fs_LCaL) * y["IVF_Ca_I_LCaL"] # fast voltage-dependent Calcium-dependent
      d[14] <- phi_s_LCaL * y["C3Ca_I_LCaL"] + psi_s_LCaL * y["ICa_I_LCaL"] + theta_LCaL * y["IVS_I_LCaL"] + omega_fs_LCaL * y["IVF_Ca_I_LCaL"] - (sigma_LCaL + omega_s_LCaL + xi_s_LCaL + omega_sf_LCaL) * y["IVS_Ca_I_LCaL"] # slow voltage dependent Calcium dependent
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
  NaQ10 <- 2.45  # - Q10 for sodium channels
  T_correction_CaL <- 1.0 * CaQ10^((Temp - 310) / 10.0) # experimental correction
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
  tau_d_CaT <- 1.9058 #for ICaL
  if(parms[["CalciumDependence"]] == 1){
    sigma_LCaL <- 0.01 #FOR ICaL set to 0 to replicate EGTA 
  }else{
    sigma_LCaL <- 0 #- this is set to 0 to replicate EGTA present
  }
  k_on <- 40633 #BK on rate
  k_c_off <- 11 #BK closed off rate
  k_o_off <- 1.1 #BK open off rate
  SMC_resting <- -60 #
  # Initial values for A, B, and C
  initial <- c(
    C0_I_LCaL = 0.563168,
    C1_I_LCaL = 0.167247,
    C2_I_LCaL = 0.0186255,
    C3_I_LCaL = 0.000921887,
    C0Ca_I_LCaL = 0.0202722,
    C1Ca_I_LCaL = 0.003,
    C2Ca_I_LCaL = 0.00602034,
    C3Ca_I_LCaL = 0.000670459,
    O_I_LCaL = 1.71111e-05,
    IVF_I_LCaL = 1.85565e-05,
    IVS_I_LCaL = 4.68845e-06,
    ICa_I_LCaL = 3.3185e-05,
    IVF_Ca_I_LCaL = 6.15944e-07,
    IVS_Ca_I_LCaL = 0
  )
  
  # Time sequence for the simulation
  times <- seq(0, 1000, by = 1)
  
  #history
  
  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$P_total <- rowSums(output_df[2:15])
  output_df$E_Ca <- ((R * Temp) / (2 * Faraday)) * log(Ca_o / Ca_i_total)
  output_df$I_CaL <- G_CaL * output_df$O_I_LCaL * 
    (ifelse(output_df$time < parms[["clamp_start"]] | output_df$time > parms[["clamp_end"]], parms[["Hv"]], parms[["Vm"]]) - output_df$E_Ca)
  
  output_df$Vm <- ifelse(
    output_df$time < parms[["clamp_start"]] | output_df$time > parms[["clamp_end"]], 
    parms[["Hv"]],  # Use Hv before clamp_start or after clamp_end
    parms[["Vm"]]   # Use Vm between clamp_start and clamp_end
  )
  
  return(output_df)
}

#Parameters to fit -----
mv_tests <- L_type$x
sim_v <- list()
peak_store <- numeric(length(mv_tests))
i <- 1
for(i in seq_along(mv_tests)) {  # Correct looping over mv_tests
  parms <- list(Vm = mv_tests[i], Hv = -100, clamp_start = 500, clamp_end = 1012, CalciumDependence = 0)  # Set voltage parameter
  sim_temp <- Model(parms = parms)  # Run simulation
  sim_temp$Vm_identity <- parms[["Vm"]]
  sim_v[[i]] <- as.data.frame(sim_temp)
  sim_temp_peak <- subset(sim_temp, Vm == mv_tests[i])
  peak_store[i] <- min(sim_temp_peak$I_CaL)
}

x <- 1
for(x in 1:length(sim_v)) {  # Correct looping over mv_tests
  sim_v[[x]]$I_CaL_norm <- sim_v[[x]]$I_CaL / min(peak_store)
}


sim_v_df <- bind_rows(sim_v)  # Convert list to dataframe

# Plot all conditions with color mapped to Vm
S6.plot <- ggplot(sim_v_df, aes(x = time, y = I_CaL_norm, color = factor(Vm_identity), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "Time vs Normalized Current",
    x = "Time (msec)",
    y = "Normalized Current",
    color = "Voltage (mV)"  # Legend label
  ) +
  scale_y_reverse()+
  xlim(480,700) +
  theme_minimal()  

# Display the plot
print(S6.plot)



# Create dataframe with results
IV.df <- data.frame(mV = mv_tests, I = peak_store)

IV.df$I_norm <- IV.df$I / min(IV.df$I)

IV.plot <- ggplot() +
  geom_point(data = L_type, aes(x = x, y = y), color = "red", size = 3) +  # Scatter points for experimental data
  geom_line(data = IV.df, aes(x = mV, y = I_norm), color = "blue", linewidth = 1) +  # Model I-V curve
  labs(
    title = "Normalized I-V Curve for I_CaL",
    x = "Membrane Voltage (mV)",
    y = "Normalized Peak I_CaL"
  ) +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(-90, -70, -70, -30, -10, 10, 30)) +
  theme_minimal()  # Clean theme

# Display the plot
print(IV.plot)

